subroutine eos(rhox, t, y, abar, &
                px, ex, cs2x, dpdtx, dedtx, sx, pradx,&
                keyerr, keytemp, eoskeyx)

  use parameters
  use physical_constants
  implicit none

! external vars
  real*8 :: rhox(imax-1), t(imax-1), y(imax-1), abar(imax-1)
  real*8 :: px(imax-1), ex(imax-1), cs2x(imax-1), dpdtx(imax-1), dedtx(imax-1),&
         sx(imax-1), pradx(imax-1)
  integer :: keyerr,keytemp,eoskeyx

! internal vars
  integer :: k(imax-1)
  integer :: i

!-----------------------------------------------------------------------------

!****************************** ideal EOS ************************************
  if(eoskeyx.eq.1) then

    call ideal_eos(rhox,t,px,ex,cs2x,dpdtx,dedtx,imax-1)

    !WARNING: this EOS does not calculate entropy and radiation pressure
    sx = 0.0d0
    pradx = 0.0d0


!***************************** PACZYNSKI EOS *********************************
  else if(eoskeyx.eq.2) then

    do i=1, imax-1
        k(i) = i
    end do

    call paczynski_eos(rhox,t,y,abar,px,ex,cs2x,dpdtx,dedtx, &
                        pradx,k,imax-1)

    !WARNING: this EOS does not calculate entropy
    sx = 0.0d0

 else
    stop "eos choice not implemented, check parameter eoskey"
    
 endif
 
end subroutine eos
