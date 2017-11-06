subroutine hydro

  use blmod
  use parameters
  use physical_constants
  implicit none

  integer :: i,ie
  integer :: keytemp,keyerr
  real*8 :: dtv

  ! Vaiables for implementing Newton-Raphson method
  real*8 :: function_star
  real*8 :: function_deriv
  real*8 :: delta_temp(imax)

  real*8 :: delta_max
  integer :: location_max

  real*8 :: p_temp(imax), eps_temp(imax), lum_temp(imax), temp_temp(imax)

  real*8, parameter :: EPSTOL = 1.0d-7
  integer, parameter :: ITMAX = 100

!------------------------------------------------------------------------------

  ! follow dt convention of MB93:
  dtv = 0.5d0 * (dtime + dtime_p)

  ! copy over data into _p arrays
  rho_p(:) = rho(:)
  vel_p(:) = vel(:)
  r_p(:)   = r(:)
  cr_p(:)  = cr(:)
  eps_p(:) = eps(:)

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

  !initial guess for the temperature T_n+1 = T_n
  temp_temp(1:imax) = temp(1:imax)
  delta_temp(1:imax) = temp_temp(1:imax)
  ie=0
  delta_max = 1

  do while(delta_max > EPSTOL)

    keytemp = 1
    call eos(rho(1:imax-1),temp_temp(1:imax-1),ye(1:imax-1), & 
             abar(1:imax-1),p_temp(1:imax-1),eps_temp(1:imax-1), &
             cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), & 
             entropy(1:imax-1),p_rad(1:imax-1),keyerr,keytemp,eoskey)
    
    delta_max = 0

    do i=1,imax-1
        
        function_star = eps_temp(i) - eps(i) + 0.5d0*(p_temp(i) + p(i)) & 
            * (1.0d0/rho(i)-1.0d0/rho_p(i)) - bomb_heating(i)*dtime &
            - Ni_heating(i)*dtime &
            - Qterm(i)

        function_deriv = dedt(i) + 0.5d0 * dpdt(i)*(1.0d0/rho(i)-1.0d0/rho_p(i))
        delta_temp(i) = -function_star/function_deriv
        temp_temp(i) = temp_temp(i) + delta_temp(i)

        temp(i) = temp_temp(i)
        
        if(abs(delta_temp(i)/temp(i)).gt.delta_max) then
            delta_max = abs(delta_temp(i)/temp(i))
            location_max = i
        end if

    enddo                 

    ie=ie+1

    if(ie.gt.ITMAX) then
        write(6,*) "EOS problem", delta_max, location_max
        scratch_step = .true.
        exit
    endif

  enddo

  keytemp = 1
  call eos(rho(1:imax-1),temp(1:imax-1),ye(1:imax-1), &
           abar(1:imax-1),p(1:imax-1),eps(1:imax-1), &
           cs2(1:imax-1), dpdt(1:imax-1), dedt(1:imax-1), & 
           entropy(1:imax-1),p_rad(1:imax-1),keyerr,keytemp,eoskey)

  !passive boundary conditions, does not participate in the evolution
  temp(imax) = 0.0d0
  eps(imax) = 0.0d0

  !active boundary condition, used in the velocity update
  p(imax) = 0.0d0

  call opacity(rho(:),temp_temp(:),kappa(:),kappa_table(:),dkappadt(:))

  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))


end subroutine hydro
