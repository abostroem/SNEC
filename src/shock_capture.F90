subroutine shock_capture
  ! this routine finds the location of th eshock wave
  
  use blmod, only: shockpos_stop, shockpos, shockpos_prev, rho, mass, vel, &
            tau, r, bomb_total_energy, time, time_prev, radshock, &
            radshock_prev, breakoutflag
  use parameters
  use physical_constants
  implicit none

  integer :: i

  real*8 :: veldiff(imax)
  real*8 :: veldiffmin
  real*8 :: rhoshock, mshock, taushock
  real*8, save :: velshock
  real*8 :: velshock_analyt

  character(len=1024) :: filename

!------------------------------------------------------------------------------

  !shockpos tracks the position of maximum negative difference in velocity
  veldiffmin = 1.0d0
  do i=shockpos_prev,imax-1
    veldiff(i) = (vel(i+1) - vel(i))
    if(veldiff(i).lt.veldiffmin) then
        veldiffmin = veldiff(i)
        shockpos = i
    endif
  enddo

  if(veldiffmin.eq.1.0d0) then
    shockpos = imax
    shockpos_stop = 1
    write(*,*) 'finished tracing the shock position'
  endif

  rhoshock = rho(shockpos)
  mshock = mass(shockpos)
  radshock = r(shockpos)
  taushock = tau(shockpos)
  !equation (19) from Matzner & McKee, ApJ 510:379-403
  velshock_analyt = 0.794*(bomb_total_energy/mshock)**(0.5d0) &
      *(mshock/(rhoshock*radshock*radshock*radshock))**0.19d0


  !velshock is not the velocity of matter at the shock, but the shock itself
  if(shockpos.ne.shockpos_prev) then

    velshock = (radshock - radshock_prev)/(time-time_prev)
    time_prev = time
    radshock_prev = radshock

    if(shockpos.ne.imax) then
        filename = trim(adjustl(outdir))//"/velshock_index.dat"
        open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
        write(666,"(I5.4, 5E15.6)") shockpos,time,radshock,velshock, &
                                            velshock_analyt,taushock
        close(666)
    endif

  endif

  !fix the moment of shock breakout
  if(breakoutflag.eq.0 .and. taushock.lt.(clite/velshock) &
                                            .and. shockpos.gt.(imax/2)) then
    
    write(*,*) "Time of breakout is", time

    open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
        status="unknown",form='formatted',position="append")
    write(666,*) 'Time of breakout = ', time, 'seconds'
    write(666,*) 'Gridpoint of breakout = ', shockpos
    close(666)
    
    breakoutflag = 1
    
  end if

  shockpos_prev = shockpos


end subroutine shock_capture
