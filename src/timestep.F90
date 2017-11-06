subroutine timestep
  ! this routine sets the timestep

  use blmod, only: cs2, r, vel, dtime, dtime_p, nt, delta_time
  use parameters
  use physical_constants
  implicit none

  integer :: i

  real*8 :: sound
  real*8 :: dttrans
  real*8, parameter :: dtfac = 0.95d0

!------------------------------------------------------------------------------

  dttrans = dtmax

  if(nt.le.1) then
    dttrans = 1.0d-8
    dtime = 1.0d-8
  endif

  !copy old timestep
  dtime_p = dtime


  dtime = dtmax
  do i=1,imax-1
    sound = sqrt(cs2(i))
    delta_time(i) = (r(i+1) - r(i)) / &
         max(abs(vel(i)+sound),abs(vel(i)-sound))
    dtime = min(dtime, (r(i+1) - r(i)) / &
         max(abs(vel(i)+sound),abs(vel(i)-sound)))
  enddo

  dtime = dtfac * dtime

  dtime = min(1.025d0*dtime_p,dtime)

  if(dttrans.lt.0.75d0*dtime) then
    dtime = dttrans
  endif

  dtime = min(max(dtime,dtmin),dtmax)


end subroutine timestep
