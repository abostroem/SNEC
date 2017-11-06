subroutine simple_saha(ncomps_consider,k,tempx,rhox,ne)

  use blmod, only: comp, comp_details, ncomps, zav, ion_fractions
  use parameters
  use physical_constants
  use eosmodule, only: xxip, stat_weight_p1_ratio
  implicit none

!input:
  integer :: k ! grid point we are working on
  integer :: ncomps_consider
  real*8  :: tempx, rhox

!output:
  real*8 :: ne

!local:
  integer, parameter :: maxstates = 30
  integer, parameter :: itmax = 500
  real*8,  parameter :: zavtol = 1.0d-15
  real*8,  parameter :: zavtolinv = 1.0d15

  real*8 :: y(maxstates)
  real*8 :: f(maxstates)
  real*8 :: potential(maxstates)

  real*8 :: stat_weight, ionpot, zion
  real*8 :: coef
  real*8 :: n_k, zav_star, delta_zav
  real*8 :: function_x, function_deriv
  real*8 :: ktinv
  real*8 :: sum0, sum1, sum2, sum3
  real*8 :: productx, denom
  real*8 :: zavstar_x_nk_inv,zavtol_inv

  integer :: i,j,ie,izion
  integer :: min_state
  integer :: max_state
  integer :: number_states

!------------------------------------------------------------------------------
! based on the work of Zaghloul et al., J. Phys. D: Appl. Phys. 33:977 (2000)
! here the state number 1 is neutral and the total number of states is Z+1

  !abort if the temperature is negative
  !(this should never happen)
  if(tempx.lt.0.0d0) then
     write(*,*) 'temperature is negative', k
     stop
  end if

  ne = 0.0d0

  ktinv = 1.0d0/(tempx*kboltz)
  ! saha_coeff defined in module physical_constants
  coef = saha_coeff * tempx**(1.5d0)

  do j=1, ncomps

     zion = comp_details(j,2)
     izion = int(zion)
     number_states = izion + 1
     n_k = rhox*comp(k,j)/(comp_details(j,1)*mproton)
     
     if(j.le.ncomps_consider.and.comp_details(j,2).gt.0.0d0) then

        !the initial guess for zav is the value from the previous timestep
        zav_star = zav(j,k)
       
        !look for the lowest and highest ionization states to solve
        !(it doesn't make sense to solve for the states with the fraction
        !of atoms potentially less than the chosen tolerance level)
        min_state = 1
        max_state = number_states

        zavstar_x_nk_inv = 1.0d0/(zav_star*n_k)
        do i = 1, number_states - 1

            f(i) = stat_weight_p1_ratio(j,i) &
                                * coef * exp(-xxip(izion,i)*ktinv + 10.0d0*rhox)

            if( (f(i)*zavstar_x_nk_inv) .gt. zavtolinv ) then
                y(i) = 0.0d0
                min_state = i + 1
            end if

            if( (f(i)*zavstar_x_nk_inv) .lt. zavtol ) then
                max_state = i
                y(i+1:number_states) = 0.0d0
                exit
            end if
        end do

        !find the fractions of atoms at different ionization states
        if(max_state.eq.1) then                     !all atoms are neutral
            y(1) = 1.0d0
            zav_star = 1.0d-6 ! avoiding zero, not to get it in denominator
        else if(min_state.eq.number_states) then    !full ionization
            y(number_states) = 1.0d0
            zav_star = zion
        else if(min_state.eq.max_state) then        !only one state is possible
            y(min_state) = 1.0d0
            zav_star = dble(min_state) - 1.0d0
        else                                        !use Newton-Raphson method

            delta_zav = zav_star
            ie = 0

            do while(abs(delta_zav/zav_star) > zavtol)
            
                ie = ie+1
            
                productx = 1.0d0

                sum0 = 0.0d0
                sum1 = 0.0d0
                sum2 = 0.0d0
                sum3 = 0.0d0

                do i = min_state, max_state - 1
                    productx = productx*f(i)/(zav_star*n_k)

                    sum0 = sum0 + productx
                    sum1 = sum1 + productx * i
                    sum2 = sum2 + productx * (i - dble(min_state) + 1.0d0)
                    sum3 = sum3 + productx * i * (i - dble(min_state) + 1.0d0)

                end do

                denom = 1.0d0/(dble(min_state) - 1.0d0 + sum1)

                function_x = 1.0d0 - zav_star * (1.0d0 + sum0) * denom

                function_deriv = (sum2 - (1+sum0) * (1+sum3*denom)) * denom

                delta_zav = -function_x/function_deriv
                zav_star = zav_star + delta_zav

                if(zav_star.gt.zion) then
                    zav_star = 0.5d0*zion
                end if

                if(ie.gt.ITMAX) then
                    write(6,*) "convergence problem in saha solver", k
                    exit
                endif

            end do

            y(min_state) = zav_star * denom
            do i=min_state+1, max_state
                y(i) = y(i-1)*f(i-1)/(zav_star*n_k)
            end do

        end if
        
        
    else !if the element is not considered, it is assumed to be fully ionized

        y(1:number_states-1) = 0.0d0
        y(number_states) = 1.0d0
        zav_star = zion
        
    end if

    ! save the ionization fractions
    do i = 1, number_states
        ion_fractions(j,i,k) = y(i)
    end do

    ne = ne + zav_star*n_k

    zav(j,k) = zav_star

  end do


end subroutine simple_saha
