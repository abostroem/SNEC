module outinfomod
  implicit none
  integer, parameter :: out_info_string_every = 200
  integer :: outinfo_count

end module outinfomod

subroutine outinfo
     
  use blmod, only: nt, time, dtime, rho, ye, cr, shockpos
  use outinfomod
  use parameters
  implicit none
  
  real*8  :: shockr

!------------------------------------------------------------------------------
  
  shockr = cr(shockpos)

  if( mod(outinfo_count,out_info_string_every).eq.0 ) then
    outinfo_count = 0
    write(6,"(A8,A15,A15,A15,A15)") "nt","time","dtime","max(rho)","shockpos"
  endif
  write(*,"(i8,1P10E15.6)") nt,time,dtime,maxval(rho),shockr

  outinfo_count = outinfo_count + 1
  
end subroutine outinfo
