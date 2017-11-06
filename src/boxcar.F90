subroutine boxcar
! This routine smoothes the abundance profiles
! to account for "mixing" (i.e. smoothing of compositional
! gradients) during the explosion

  use blmod, only: comp, ncomps, delta_mass
  use parameters
  use physical_constants
  implicit none

  real*8, parameter  :: boxcar_nominal_mass = 0.4*msun
  integer, parameter :: number_iterations = 4

  !actual boxcar width is different from nominal, as the mass is discrete
  real*8 :: boxcar_actual_mass

  real*8 :: mass_element !mass of an element within the actual boxcar mass
  integer :: i, l, l_max, k, n
  integer :: el

!------------------------------------------------------------------------------

  do n=1, number_iterations !repeat until the desired smoothness

      do i=1, imax

          do l=i, imax
              if( sum(delta_mass(i:l)).gt.boxcar_nominal_mass ) then
                l_max = l
                exit
              endif
              if( l.eq.imax ) then
                l_max = imax
                exit
              endif
          enddo

          boxcar_actual_mass = sum(delta_mass(i:l_max))

          do el=1, ncomps
              mass_element = sum(delta_mass(i:l_max)*comp(i:l_max,el))
              do k = i, l_max
                  comp(k,el) = mass_element/boxcar_actual_mass
              enddo
          enddo

          !stop the averaging, if the boxcar reached the outer boundary
          if( boxcar_actual_mass .lt. boxcar_nominal_mass ) exit

      enddo

  enddo

end subroutine boxcar
