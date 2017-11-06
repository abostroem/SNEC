subroutine output_all(modeflag)

  use blmod
  use parameters
  use physical_constants
  implicit none

  character(len=1024) :: filename
  character(len=256) :: basename
  integer :: modeflag

!------------------------------------------------------------------------------

  if(modeflag.eq.0) then

  ! meant for checkpoints; not used at the moment

  else if(modeflag.eq.1) then
         
    filename = trim(adjustl(outdir))//"/vel.xg"
    call output_single_mass(vel,filename)

    filename = trim(adjustl(outdir))//"/rho.xg"
    call output_single_mass(rho,filename)

    filename = trim(adjustl(outdir))//"/ye.xg"
    call output_single_mass(ye,filename)

    filename = trim(adjustl(outdir))//"/press.xg"
    call output_single_mass(p,filename)

    filename = trim(adjustl(outdir))//"/cs2.xg"
    call output_single_mass(cs2,filename)

    filename = trim(adjustl(outdir))//"/Q.xg"
    call output_single_mass(Q,filename)

    filename = trim(adjustl(outdir))//"/eps.xg"
    call output_single_mass(eps,filename)

    filename = trim(adjustl(outdir))//"/mass.xg"
    call output_single_radius(mass,filename)

    filename = trim(adjustl(outdir))//"/temp.xg"
    call output_single_mass(temp,filename)

    filename = trim(adjustl(outdir))//"/lum.xg"
    call output_single_mass(lum,filename)

    filename = trim(adjustl(outdir))//"/tau.xg"
    call output_single_mass(tau,filename)

    filename = trim(adjustl(outdir))//"/delta_time.xg"
    call output_single_mass(delta_time,filename)

    filename = trim(adjustl(outdir))//"/radius.xg"
    call output_single_mass(r,filename)

    filename = trim(adjustl(outdir))//"/kappa.xg"
    call output_single_mass(kappa,filename)

    filename = trim(adjustl(outdir))//"/kappa_table.xg"
    call output_single_mass(kappa_table,filename)

    filename = trim(adjustl(outdir))//"/logR_op.xg"
    call output_single_mass(logR_op,filename)

    filename = trim(adjustl(outdir))//"/logT.xg"
    call output_single_mass(logT,filename)

    filename = trim(adjustl(outdir))//"/p_rad.xg"
    call output_single_mass(p_rad,filename)

    filename = trim(adjustl(outdir))//"/Ni_deposit_function.xg"
    call output_single_mass(Ni_deposit_function,filename)

    filename = trim(adjustl(outdir))//"/He_1.xg"
    call output_single_mass(ion_fractions(He_number,1,:),filename)

    filename = trim(adjustl(outdir))//"/He_2.xg"
    call output_single_mass(ion_fractions(He_number,2,:),filename)

    filename = trim(adjustl(outdir))//"/He_3.xg"
    call output_single_mass(ion_fractions(He_number,3,:),filename)

    filename = trim(adjustl(outdir))//"/H_1.xg"
    call output_single_mass(ion_fractions(H_number,1,:),filename)

    filename = trim(adjustl(outdir))//"/H_2.xg"
    call output_single_mass(ion_fractions(H_number,2,:),filename)

    filename = trim(adjustl(outdir))//"/free_electron_frac.xg"
    call output_single_mass(free_electron_frac,filename)

    filename = trim(adjustl(outdir))//"/E_shell.xg"
    call output_single_mass(E_shell,filename)

    filename = trim(adjustl(outdir))//"/time_diff.xg"
    call output_single_mass(time_diff,filename)

    filename = trim(adjustl(outdir))//"/time_exp.xg"
    call output_single_mass(time_exp,filename)

    filename = trim(adjustl(outdir))//"/photosphere_tracer.xg"
    call output_single_mass(photosphere_tracer,filename)

  else if(modeflag.eq.2) then
     
     filename = trim(adjustl(outdir))//"/T_eff.dat"
     call output_scalar(T_eff,filename)
     
     filename = trim(adjustl(outdir))//"/Ni_total_luminosity.dat"
     call output_scalar(Ni_total_luminosity,filename)
     
     filename = trim(adjustl(outdir))//"/lum_observed.dat"
     call output_scalar(lum_observed,filename)
         
     filename = trim(adjustl(outdir))//"/index_photo.dat"
     call output_integer(index_photo,filename)
     
     filename = trim(adjustl(outdir))//"/lum_photo.dat"
     call output_scalar(lum_photo,filename)
 
     filename = trim(adjustl(outdir))//"/mass_photo.dat"
     call output_scalar(mass_photo,filename)
 
     filename = trim(adjustl(outdir))//"/vel_photo.dat"
     call output_scalar(vel_photo,filename)
 
     filename = trim(adjustl(outdir))//"/rad_photo.dat"
     call output_scalar(rad_photo,filename)
     
     filename = trim(adjustl(outdir))//"/opacity_corrupted.dat"
     call output_integer(opacity_corrupted,filename)
     
     filename = trim(adjustl(outdir))//"/index_lumshell.dat"
     call output_integer(index_lumshell,filename)
 
     filename = trim(adjustl(outdir))//"/mass_lumshell.dat"
     call output_scalar(mass_lumshell,filename)

  endif


end subroutine output_all

! *******************************************************************

subroutine output_single_mass(var,filename)
  
  use blmod, only: mass,time
  use parameters

  implicit none
  real*8 var(*)
  character(len=100) filename
  integer nt
  integer i



  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  write(666,*) '"Time = ',time

  do i=1,imax
     write(666,"(1P20E29.20E3)") mass(i),var(i)
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)



end subroutine output_single_mass

! *******************************************************************

subroutine output_single_radius(var,filename)
  
  use blmod, only: r,time
  use parameters

  implicit none
  real*8 var(*)
  character(*) filename
  integer nt
  integer i


  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  write(666,*) '"Time = ',time

  do i=1,imax
     write(666,"(1P20E19.10E3)") r(i),var(i)
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)



end subroutine output_single_radius
    
! *******************************************************************

subroutine output_single_mass_integer(var,filename)

  use blmod, only: mass,time
  use parameters

  implicit none
  integer var(*)
  character(*) filename
  integer nt
  integer i


  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  write(666,*) '"Time = ',time

  do i=1,imax
     write(666,"(E19.10E3, I5.4)") mass(i), var(i)
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)


end subroutine output_single_mass_integer
! ******************************************************************

subroutine output_central(var,filename)
  
  use blmod, only: r,mass,time
  use parameters

  implicit none
  real*8 var(*)
  character(*) filename
  integer nt
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")

  write(666,"(1P20E19.10E3)") time,var(1)
  
  close(666)


end subroutine output_central
    
! ******************************************************************

subroutine output_outer(var,filename)
  
  use blmod, only: r,mass,time
  use parameters

  implicit none
  real*8 var(*)
  character(*) filename
  integer nt
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                        form='formatted',position="append")

  write(666,"(1P20E19.10E3)") time,var(imax)

  close(666)

end subroutine output_outer


! *******************************************************************
subroutine output_scalar(var,filename)

  use blmod, only: time

  implicit none
  real*8 var
  character(len=100) filename
  integer nt
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")

  write(666,"(1P20E19.10E3)") time,var
  
  close(666)

end subroutine output_scalar
    
! *******************************************************************
subroutine output_integer(var,filename)

  use blmod, only: time

  implicit none
  integer var
  character(len=100) filename
  integer nt
  integer i

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
                                            form='formatted',position="append")
  
  write(666,"(E19.10E3, I5.4)") time,var

  close(666)

end subroutine output_integer


!******************************************************************************
!output of the variable versus grid point number for a given moment of time
!used to output the initial values of some variables
subroutine output_screenshot(var,filename,imaximum)

  implicit none
  real*8 var(*)
  character(len=100) filename
  integer nt
  integer i
  integer imaximum

  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
       form='formatted',position="append")
  
  do i=1, imaximum
      write(666,"(I5.4, E25.16E3)") i,var(i)
  enddo

  close(666)

end subroutine output_screenshot


! *******************************************************************
    subroutine generate_filename(varname,outdir,time,nt,suffix,fname)

      implicit none


      real*8 time
      integer nt
      character(*) varname
      character(len=256) outdir
      character*(*) suffix
      character*(*) fname
      character*(400) aa
      character(len=100) outtime
      character(len=20) outnt
      integer i,ii

      aa=" "
      fname=" "
!      write(aa,"(a32,'_',f10.8,'_nt',i6,'.dat')") varname,time,nt
      write(outnt,"(i10.10)") nt
!      write(aa,"(a32,'_nt',i6)") varname,nt

      fname = trim(adjustl(outdir))//"/"//trim(adjustl(varname))//"_nt_"//outnt
      write(outtime,"(f14.7)") time
!      write(*,*) aa

!      ii=0
!      do i=1,80
!         if(aa(i:i).ne.' ') then
!            ii=ii+1
!            fname(ii:ii)=aa(i:i)
!         endif
!      enddo

      fname = trim(adjustl(fname))//"_time_"//trim(adjustl(outtime))

    end subroutine generate_filename




