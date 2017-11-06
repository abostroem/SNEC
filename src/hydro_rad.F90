subroutine hydro_rad

  use blmod
  use parameters
  use physical_constants
  implicit none

  integer :: i,k
  integer :: keytemp,keyerr
  real*8 :: dtv

  ! Variables for inversion of Jacobian
  external :: dgbsv
  integer :: info
  integer, parameter :: kl = 1
  integer, parameter :: ku = 1
  integer, parameter :: ldab=2*kl+ku+1
  real*8 :: ab(ldab,imax-1), b(imax-1)
  integer :: ipiv(imax-1)
  
  real*8 :: delta_max
  integer :: location_max
  
  real*8 :: Aarray(imax-1), Barray(imax-1), Carray(imax-1), Darray(imax-1)
  real*8 :: p_temp(imax), eps_temp(imax), lum_temp(imax), temp_temp(imax)
  real*8 :: lambda_temp(imax)

  real*8, parameter :: EPSTOL = 1.0d-7
  integer, parameter :: ITMAX = 300


!------------------------------------------------------------------------------

  ! follow dt convention of MB93:
  dtv = 0.5d0 * (dtime + dtime_p)

  ! copy over data into _p arrays
  rho_p(:) = rho(:)
  vel_p(:) = vel(:)
  r_p(:)   = r(:)
  cr_p(:)  = cr(:)
  eps_p(:) = eps(:)
  p_p(:) = p(:)
  temp_p(:) = temp(:)
  kappa_p(:) = kappa(:)


!---------------------------- update velocities -------------------------------
  if(do_piston .and. time.ge.piston_tstart .and. time.le.piston_tend) then
        vel(1) = piston_vel
        vel(2) = piston_vel
  endif

  do i=2,imax
    vel(i) = vel_p(i) &

     ! gravity
     - dtv * ggrav*mass(i) / r(i)**2 *gravity_switch   &

     ! pressure
     - dtv * 4.0d0*pi*r(i)**2 * (p(i) - p(i-1)) / delta_cmass(i-1)   &

     ! artificial viscosity
     - dtv * 4.0d0*pi * (cr(i)**2 * Q(i) - cr(i-1)**2 * Q(i-1))/delta_cmass(i-1)
  enddo

  if(do_piston.and.time.ge.piston_tend) then
    vel(1) = 0.0d0
  else if(do_bomb) then
    vel(1) = 0.0d0
  endif

!----------------------- update the radial coordinates-------------------------
  do i=1,imax
   r(i) = r_p(i) + dtime * vel(i)
   
   if(i.gt.1) then
       if (r(i).lt.r(i-1)) then
           write(*,*) 'radius of a gridpoint', i, 'is less than preceding'
           stop
       end if
   end if
  enddo

!------------------------- update the zone densities --------------------------
  do i=1,imax-1
     rho(i) = delta_mass(i) / (4.0d0*pi * (r(i+1)**3 - r(i)**3)/3.0d0)
  enddo
  rho(imax) = 0.0d0 !passive boundary condition


!------------------------- update zone center radius --------------------------
  do i=1,imax-1
    cr(i) = ( ( r(i)**3 + r(i+1)**3 ) / 2.0d0 )**(1.0d0/3.0d0)
  enddo
  cr(imax) = r(imax) + (r(imax) - cr(imax-1))
  !passive boundary condition, used in the expression for the velocity update,
  !but multiplied by the artificial viscosity, which is zero at the last point


  ! update the artificial viscosity
  call artificial_viscosity


!----------- update the temperature, pressure and internal energy -------------

  !calculate heating term due to Ni
  if(time.ge.time_Ni) then
      time_Ni = time_Ni + Ni_period
      call nickel_heating
  endif

  !calculate heating term due to bomb
  if(do_bomb .and. time.ge.bomb_tstart .and. time.le.bomb_tend) then
      call bomb_pattern
  else
      bomb_heating(:) = 0.0d0
  endif

  !Initial guess for the quantities at the next time step
  p_temp(1:imax) = p(1:imax)
  eps_temp(1:imax) = eps(1:imax)
  temp_temp(1:imax) = temp(1:imax)

  do k=1, ITMAX

    keytemp = 1
    call eos(rho(1:imax-1),temp_temp(1:imax-1),ye(1:imax-1), &
         abar(1:imax-1),p_temp(1:imax-1),eps_temp(1:imax-1), &
         cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), entropy(1:imax-1), &
         p_rad(1:imax-1),keyerr,keytemp,eoskey)


    call luminosity(r(:), temp_temp(:), kappa_p(:), &
                    lambda_temp(:), inv_kappa(:), lum_temp(:))

    !calculate the coefficients of the equation:
    ! A(i)*\delta T(i+1) + B(i)*\delta T(i) + C(i)*\delta T(i-1) = D(i)
    call matrix_arrays(temp_temp(:), lambda_temp(:), inv_kappa(:), &
        eps_temp(:), p_temp(:), lum_temp(:), &
        Aarray(:), Barray(:), Carray(:), Darray(:))

             
    !assemble the matrix in the form used by lapack
    ab(2,2:imax-1) = Aarray(1:imax-2)
    ab(3,1:imax-1) = Barray(1:imax-1)
    ab(4,1:imax-2) = Carray(2:imax-1)
 
    b(1:imax-1) = Darray(1:imax-1)

    !invert the matrix
    !if the inversion fails, the whole matrix is written in 'failed_matrix.dat'
    info = 0
    call dgbsv(imax-1,kl,ku,1,ab,ldab,ipiv,b,imax-1,info)

    if(info.ne.0) then
       open(unit=666, &
           file=trim(adjustl(trim(adjustl(outdir))//"/failed_matrix.dat")), &
           status="unknown",form='formatted',position="append")
       do i=1,imax-1
           write(666,*) ab(2,i), ab(3,i), ab(4,i), Darray(i)
       enddo
       close(666)
       stop "problem in the matrix inversion (see Data/failed_matrix.dat)"
    endif

    !check if the iteration procedure converged
    delta_max = 0.0d0
    do i=1,imax-1
        if(abs(b(i)/temp_temp(i)).gt.delta_max) then
            delta_max = abs(b(i)/temp_temp(i))
            location_max = i
        endif
    enddo
    if(delta_max.le.EPSTOL) goto 101

    !add the increment to the temperature
    do i=1, imax-1
        temp_temp(i) = temp_temp(i) + b(i)
        if(temp_temp(i).lt.0.0d0) then
            goto 100
        end if
    end do

  enddo

  100 continue

  write(6,*) "EOS problem", delta_max, location_max
  scratch_step = .true.

  101 continue

        
  eps(1:imax-1) = eps_temp(1:imax-1)
  p(1:imax-1)   = p_temp(1:imax-1)
  temp(1:imax-1)  = temp_temp(1:imax-1)


  !passive boundary conditions, do not participate in the evolution
  temp(imax) = 0.0d0
  eps(imax) = 0.0d0

  !active boundary condition, used in the velocity update
  p(imax) = 0.0d0


  call opacity(rho(:),temp_temp(:),kappa(:),kappa_table(:),dkappadt(:))
    
  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))


end subroutine hydro_rad
