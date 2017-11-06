subroutine analysis

  use blmod
  use parameters
  use physical_constants
  implicit none

  character(len=1024) :: filename

  integer :: i

  real*8 :: T_eff_for_BC
  real*8 :: bol_corr_used(11)
  real*8, parameter :: T_eff_min = 5000.0d0

!------------------------------------------------------------------------------
!---------- Calculating optical depth and tracing the photosphere -------------

  !kappa_table (without the opacity floor) is used to trace the photosphere
  call optical_depth(rho(:), r(:), kappa_table(:), tau(:))

  !find the grid point, where the photosphere is located
  do i=imax-1, 1, -1
     if(tau(i).gt.0.66d0) then
        index_photo = i + 1
        exit
     endif
  enddo

  !fix the moment, when the photosphere reaches the inner boundary, if it does
  if(tau(1).lt.0.66d0) then
     index_photo = 1
     if(photosphere_fell_on_the_center.eq.0) then
        photosphere_fell_on_the_center = 1
        open(unit=666, &
             file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
             status="unknown",form='formatted',position="append")
        write(666,*) 'Photosphere reached the center at ', time, 'seconds'
        close(666)
     end if
  end if

  !find the values of some variables at the photosphere
  if(photosphere_fell_on_the_center.eq.0) then
     call map_map(lum_photo,  0.66d0, lum(imax:1:-1),  tau(imax:1:-1),imax)
     call map_map(mass_photo, 0.66d0, mass(imax:1:-1), tau(imax:1:-1),imax)
     call map_map(vel_photo,  0.66d0, vel(imax:1:-1),  tau(imax:1:-1),imax)
     call map_map(rad_photo,  0.66d0, r(imax:1:-1),    tau(imax:1:-1),imax)
  else
     lum_photo   = lum(1)
     mass_photo  = mass(1)
     vel_photo   = vel(1)
     rad_photo   = r(1)
  end if

  !photosphere tracer is equal to 1 at the photosphere, and 0 everywhere else
  !used for visualization
  photosphere_tracer(:) = 0.0d0
  photosphere_tracer(index_photo) = 1.0d0

  !check if the photosphere moves through the regions with wrong opacity
  !log10(T) = 3.75 is the lower boundary of the OPAL opacity tables
  if(metallicity(index_photo).gt.envelope_metallicity .and. &
       log10(temp(index_photo)).lt.3.75d0 .and. &
       kappa_table(index_photo).gt.kappa(index_photo) &
       .and. index_photo.gt.1 ) then
     opacity_corrupted = 1
  else
     opacity_corrupted = 0
  end if

!------------------------ Tracing the luminosity shell ------------------------
  index_lumshell = imax
  do i=1, imax - 1
     if(tau(imax-i).gt.(clite/vel(imax-i))) then
        index_lumshell = imax - i + 1
        exit
     end if
  end do
  if (tau(2).le.(clite/vel(2))) then
     index_lumshell = 1
  end if
  mass_lumshell = mass(index_lumshell)
  
  !characteristic diffusion and expansion times for different shells
  do i=1, imax-1
     time_diff(i) = kappa(i)*rho(i)*(r(imax)-r(i))**2.0/clite
     time_exp(i) = (r(imax)-r(i))/(vel(imax)-vel(i))
  end do

  !internal energies of shells from the given radius out to the surface
  do i=1, imax
     E_shell(i) = sum(eps(i:imax)*delta_mass(i:imax))
  end do

!------------------- Calculate the observed luminosity ------------------------

  !observed luminosity is the sum of lum_photosphere and Ni contribution
  if(photosphere_fell_on_the_center.eq.0) then
     lum_observed = lum_photo + sum(Ni_energy_rate* &
          Ni_deposit_function(index_photo:imax)*delta_mass(index_photo:imax))
  else
     lum_observed = lum(1) + &
          sum(Ni_energy_rate*Ni_deposit_function(1:imax)*delta_mass(1:imax))
  end if

  !write down the time when the contribution of the Ni above the
  !photosphere to the luminosity is greater than 5%
  if(shockpos_stop.eq.1 .and. &
       abs((lum_observed - lum_photo)/lum_photo).gt.0.05 .and. &
       Ni_contributes_five_percents.eq.0) then
     
     Ni_contributes_five_percents = 1
     open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
          status="unknown",form='formatted',position="append")
     write(666,*) 'Ni contribution to the luminosity is 5% at ', time, 'seconds'
     close(666)

  end if

!-------- Calculation of color magnitudes using bolometric corrections --------

  T_eff = (lum_photo/(4.0d0*pi*sigma_SB*rad_photo**2))**0.25d0

  !see Eq.(3) of Swartz et al., ApJ 374:266 (1991) and explanation there
  T_eff_for_BC = MAX(T_eff,T_eff_min)

  !here in cases, when the effective temperature goes beyond the boundaries
  !of the table BolCorr.dat, the linear extrapolation is used
  do i=1, 11
     call map_map(bol_corr_used(i),T_eff_for_BC,bol_corr(1:nlines_bol_corr,i), &
          temp_bol_corr,nlines_bol_corr)
     magnitudes(i) = sun_mag - bol_corr_used(i) - 2.5d0*log10(lum_photo/sun_lum)
  end do

  !write down the magnitudes in a file
  if(time.eq.0.0d0.or.time.gt.tdump_scalar) then
     open(666,file=trim(adjustl(outdir))//"/magnitudes.dat",&
          status='unknown',position='append')
     write(666,"(13E18.9)") time, T_eff_for_BC, magnitudes(1:11)
     close(666)
  endif

end subroutine analysis
