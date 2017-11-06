subroutine nickel_heating

  use blmod, only: ye, comp, Ni_heating, Ni_total_luminosity, Ni_number, &
                   r, rho, time, Ni_deposit_function, Ni_energy_rate, delta_mass
  use parameters
  use physical_constants
  implicit none

!local:

  real*8 :: r_Ni !boundary of the region with non-negligible mass fraction of Ni
  real*8 :: delta_r
  real*8 :: r_x, r_max
  real*8 :: r_j, comp_Ni_j
  real*8 :: delta_tau_j

  real*8 :: th, delta_th
  real*8 :: th_min(imax)

  real*8 :: I_prime, I_prime_av

  integer :: i, i_Ni

  integer :: index, ibuffer

  ! parameter constants
  real*8, parameter :: nimin = 1.0d-5
  real*8, parameter :: th_max = pi
  !numbers of points per radial/angular integration
  integer, parameter :: npoints_radial_integration  = 150
  integer, parameter :: npoints_angular_integration = 150

!------------------------------------------------------------------------------
!based on the work of Swartz et al., ApJ 446:766 (1995)

  
  if(Ni_switch.eq.0) then 
     !****** heating by Ni is not taken into account *******

     Ni_heating(:) = 0.0d0
     Ni_deposit_function(:) = 0.0d0
     Ni_total_luminosity = 0.0d0

  else 
     !*** solve for the local heating due to the radioactive decay of Ni ****

     i_Ni = 0
     do i=imax, 1, -1
        if(comp(i,Ni_number).gt.nimin) then
           r_Ni = r(i)
           i_Ni = i
           exit
        endif
     enddo
     
     ! find the limits of integration with respect to the polar angle
     if(i_Ni.eq.0) then
        write(*,*) 'mass fraction of Ni is lower than NIMIN in every grid point'
        write(*,*) 'try reducing NIMIN'
        stop
     else if(i_Ni.eq.1) then
        th_min(i_Ni) = pi*0.5d0
        do i = i_Ni+1, imax
           th_min(i) = acos(-sqrt(r(i)*r(i)-r_Ni*r_Ni)/r(i))
        end do
     else if(i_Ni.eq.imax) then
        th_min(1:i_Ni-1) = 0.0d0
        th_min(i_Ni) = pi*0.5d0
     else
        th_min(1:i_Ni-1) = 0.0d0
        th_min(i_Ni) = pi*0.5d0
        do i = i_Ni+1, imax
           th_min(i) = acos(-sqrt(r(i)*r(i)-r_Ni*r_Ni)/r(i))
        end do
     end if


     ! find the deposition function at each grid point
     do i=1, imax
        th = th_max
        I_prime_av = 0
        delta_th = (th_max-th_min(i))/npoints_angular_integration

        do while(th.gt.(th_min(i)+1.d-14))
           r_max = -r(i)*cos(th) + sqrt((r(i)*cos(th))**2-(r(i)**2-r_Ni**2))
           delta_r = r_max/npoints_radial_integration
           I_prime = 0
           r_x = r_max
           
           do while(r_x.gt.0)
              r_j = sqrt(r(i)*r(i) + r_x*r_x + 2.0d0*r(i)*r_x*cos(th))
              
              if(r_j.le.r(1)) then !inside the excised region
                 delta_tau_j = 0.0d0
                 comp_Ni_j = 0.0d0
              else if(r_j.ge.r(imax-1)) then
                 delta_tau_j = delta_r * ye(imax-1) * 0.06d0 * rho(imax-1)
                 comp_Ni_j = comp(imax-1,Ni_number)
              else
                 call map_find_index(imax,r,r_j,ibuffer,index)
                 delta_tau_j = delta_r * ye(index) * 0.06d0 * rho(index)
                 comp_Ni_j = comp(index,Ni_number)
              end if

              I_prime = (I_prime-comp_Ni_j)*exp(-delta_tau_j) + comp_Ni_j
              r_x = r_x - delta_r
           end do

           I_prime_av = I_prime_av + I_prime*sin(th)*delta_th*0.5d0
           th = th - delta_th
        end do
        
        Ni_deposit_function(i) = I_prime_av
     end do

     ! rate of energy release per gram of radioactive material
     Ni_energy_rate = &
          3.24d10*exp(-time*overtau_Ni) + 7.29d9*exp(-time*overtau_Co)

     ! local rate of gamma-ray energy deposition at a given grid point
     Ni_heating(1:imax) = Ni_energy_rate*Ni_deposit_function(1:imax)

     ! total energy per second, deposited to the model by gamma-rays
     Ni_total_luminosity = &
          Ni_energy_rate*sum(Ni_deposit_function(1:imax)*delta_mass(1:imax))

  end if

end subroutine nickel_heating
