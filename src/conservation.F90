subroutine conservation_compute_energies
! this routine computes the various energies to check
! for energy conservation
  
  use blmod, only: nt, delta_mass, cmass, cr, gravity_switch, eps, vel, &
                    total_initial_energy, time, tdump_scalar
  use parameters
  use physical_constants
  implicit none

  integer :: i
  logical :: force

  real*8 :: egrav, eint, ekin

!------------------------------------------------------------------------------
  
  egrav = 0.0d0
  eint = 0.0d0
  ekin = 0.0d0

  do i=1,imax
     egrav = egrav - ggrav*delta_mass(i)*cmass(i)/cr(i) *gravity_switch
     eint = eint + eps(i)*delta_mass(i)
  enddo
  
  do i=1,imax-1
     ekin = ekin + 0.5d0*(0.50d0*(vel(i+1)+vel(i)))**2 * delta_mass(i)
  enddo
  ekin = ekin + 0.5d0*(vel(imax))**2 * delta_mass(imax)
  
  if(time.eq.0.0d0) then
      total_initial_energy = egrav+eint+ekin
  endif

#if 0
  if(mod(nt,1000).eq.0.or.force) then
     write(6,"(A18,A18,A18,A18)") "egrav","eint","ekin","etot"
     write(6,"(1P10E18.9)") egrav,eint,ekin,egrav+eint+ekin
  endif
#endif

  if(time.eq.0.0d0.or.time.gt.tdump_scalar) then
     open(666,file=trim(adjustl(outdir))//"/conservation.dat",&
          status='unknown',position='append')
     write(666,"(1P10E18.9)") time,egrav,eint,ekin,egrav+eint+ekin, & 
         egrav+eint+ekin-total_initial_energy
     close(666)
  endif

end subroutine conservation_compute_energies

