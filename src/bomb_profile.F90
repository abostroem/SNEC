subroutine bomb_pattern
  
  use blmod, only: bomb_heating, bomb_total_energy, time, mass, delta_mass, &
                    bomb_spread
  use parameters
  use physical_constants
  implicit none

  real*8 :: exponent_array(bomb_spread)
  real*8 :: coef_A
  real*8 :: coef_B
  real*8 :: coef_C
  real*8 :: coef_D
  real*8 :: current_luminosity

  integer :: bomb_end_point
  integer :: i

  real*8, parameter :: ratio_time = 100.0d0
  real*8, parameter :: ratio_mass = 100.0d0


!------------------------------------------------------------------------------

!energy of the bomb is injected into the model exponentially in time:
!          energy per unit time = coef_D * EXP(- coef_C * time)
!where               bomb_tstart < time < bomb_tend 

!and exponentially in mass coordinate:
!      energy per unit time per unit mass = coef_B * EXP(- coef_A * mass(i))
!where     mass(bomb_start_point) < mass(i) < mass(bomb_end_point)

!ratio_time gives the ratio between the bomb luminosity at time = bomb_tstart
!and the bomb luminosity at time = bomb_tend

!ratio_mass gives the ratio between the bomb heating at 
!mass = mass(bomb_start_point) and the bomb heating at mass(bomb_end_point)

  bomb_heating(:) = 0.0d0

  coef_C = log(ratio_time)/(bomb_tend - bomb_tstart)

  coef_D = coef_C*bomb_total_energy &
                      /(exp(-coef_C*bomb_tstart)-exp(-coef_C*bomb_tend))

  current_luminosity = coef_D*exp(-coef_C*time)

  if (bomb_spread.gt.1) then

      bomb_end_point = bomb_start_point + bomb_spread - 1

      coef_A = log(ratio_mass)/(mass(bomb_end_point) - mass(bomb_start_point))

      do i=bomb_start_point, bomb_end_point
        exponent_array(i) = delta_mass(i) * exp( - coef_A * mass(i) )
      end do

      coef_B = current_luminosity/sum(exponent_array(1:bomb_spread))

      do i=bomb_start_point, bomb_end_point
        bomb_heating(i) = coef_B * exp( - coef_A * mass(i) )
      end do

  else

      bomb_heating(bomb_start_point) = current_luminosity &
                                            /delta_mass(bomb_start_point)

  end if

end subroutine bomb_pattern
