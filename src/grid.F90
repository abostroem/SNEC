subroutine grid
    
  use blmod, only: mass, cmass, delta_mass, delta_cmass
  use parameters
  use physical_constants
  implicit none

  real*8 :: dmass
  real*8 :: grid_pattern(imax)

  integer :: i
  integer :: number_lines_GridPattern

!------------------------------------------------------------------------------
    
  if(gridding.eq.'uniform_in_mass') then
      
      dmass = (mass(imax) - mass(1))/(imax-1)
      do i=2,imax-1
         mass(i) = mass(i-1) + dmass
      enddo

  else if(gridding.eq.'from_file_by_mass') then

      open(666,file=trim("tables/GridPattern.dat"),status='unknown', &
            form='formatted',action='read')
      number_lines_GridPattern = 0
      do
        read(666,*,end=15)
        number_lines_GridPattern = number_lines_GridPattern + 1
      end do
      15 close(666)

      if(number_lines_GridPattern.ne.imax) then
        write(*,*) '******* Number of lines in the file GridPattern.dat'
        write(*,*) '******* does not coincide with the number of grid points.'
        write(*,*) '******* Please, adjust one of the two.'
        stop
      end if

      open(666,file=trim("tables/GridPattern.dat"),status='unknown', &
           form='formatted',action='read')
      do i=1,imax
         read(666,*) grid_pattern(i)
      enddo
      close(666)
      
      do i = 2, imax
          mass(i) = mass(i-1) &
              + (grid_pattern(i)-grid_pattern(i-1))*(mass(imax) - mass(1))
      enddo
      
  else
      
      write(*,*) "the chosen type of gridding is not implemented"
      stop
      
  end if

  do i=1,imax-1
     cmass(i) = mass(i) + 0.5d0*(mass(i+1)-mass(i))
  enddo
  cmass(imax) = mass(imax) + 0.5d0*(mass(imax)-mass(imax-1))

  do i = 1, imax-1
      delta_mass(i) = mass(i+1) - mass(i)
      delta_cmass(i) = cmass(i+1) - cmass(i)
  enddo
  delta_mass(imax) = delta_mass(imax-1)
  delta_cmass(imax) = delta_cmass(imax-1)
    
end subroutine grid
