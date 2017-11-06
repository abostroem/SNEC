subroutine problem

  use blmod, only: ntstart, tstart, ncomps, mass, r, Rstar, &
            opacity_floor, metallicity, envelope_metallicity, do_piston, &
            do_bomb, total_initial_energy, bomb_total_energy, bomb_spread, &
            rho, kappa, kappa_table, dkappadt, tau, temp, p, delta_mass, &
            lambda, inv_kappa, lum, eos_gamma1
  use parameters
  use eosmodule
  use physical_constants
  implicit none

  real*8 :: buffer(imax)

  integer :: i
  character(len=256) :: filename

!for OPAL interpolation routine
  real*4 :: opact,dopact,dopacr,dopactd
  common/e/ opact,dopact,dopacr,dopactd

!------------------------------------------------------------------------------

  ntstart       = 0
  tstart        = 0.0d0

!****************************** EOS setup *************************************

  if(eoskey.eq.1) then
      write(6,*) "Using the ideal EOS."
     write(6,*) "WARNING: this EOS does not return the radiation pressure."
      eos_gamma1 = 1.4d0 !for Sedov
  else if(eoskey.eq.2) then
      write(6,*) "Using the Paczynski EOS!"
  else
     stop "Choice of EOS not available, check the parameter eoskey."
  endif

!********************* Allocate and initialize variables **********************

  call get_ncomps_from_profile(composition_profile_name,ncomps)

  call allocate_vars

  call initialize_vars

!***************************** Grid setup *************************************

  call get_inner_outer_mass_from_profile(profile_name,mass(1),mass(imax))

  if(mass_excision) then
      mass(1) = mass_excised*msun
  endif

  call get_inner_outer_radius_from_profile(profile_name,mass(1),r(1),r(imax))

  Rstar=r(imax)

  open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
      status="unknown",form='formatted',position="append")
  write(666,*) 'Mass of the model = ', mass(imax)/msun, 'solar masses'
  write(666,*) 'Initial radius = ', Rstar/rsun, 'solar radii'
  close(666)

  ! set up the grid
  call grid

!*************************** Read the profile *********************************

  write(*,*) "Profile file: ",trim(profile_name)
  call read_profile(profile_name)

  !set up radius coordinates based on mass and density
  call integrate_radius_initial

!************************ Set up the opacity floor ****************************

  do i=1, imax
    opacity_floor(i) = (envelope_metallicity*of_core - of_env      &
                + metallicity(i)*(of_env - of_core))/(envelope_metallicity - 1)
  end do

  filename = trim(adjustl(outdir))//"/opacity_floor.dat"
  call output_screenshot(opacity_floor,filename,imax)


!******************** Set up the energy of the thermal bomb *******************

  call conservation_compute_energies

  if(initial_data.eq."Piston_Explosion") then

    do_piston = .true.

  else if(initial_data.eq."Thermal_Bomb") then

    do_bomb = .true.

    if(bomb_mode.eq.1) then
       bomb_total_energy = final_energy - total_initial_energy
    else if (bomb_mode.eq.2) then
       bomb_total_energy = final_energy
    else
       write(6,"(A10,I3,A12)") "bomb_mode ",bomb_mode," not defined!"
       stop
    endif

    write(6,"(A60)") "***************************************************************************"
    write(6,"(A31,I3)") "Operating in Thermal Bomb Mode", bomb_mode

    do i=1, imax
        if(mass(i).ge.(mass(bomb_start_point)+bomb_mass_spread*msun)) then
            bomb_spread = i - bomb_start_point
            open(unit=666, &
                file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
                status="unknown",form='formatted',position="append")
            write(666,*) 'Bomb energy is spread over = ', bomb_spread, 'points'
            write(6,*) 'Bomb energy is spread over = ', bomb_spread, 'points'
            close(666)
            exit
        end if
    end do


    open(unit=666,file=trim(adjustl(trim(adjustl(outdir))//"/info.dat")), &
        status="unknown",form='formatted',position="append")
    write(666,*) 'Total energy of the model = ', total_initial_energy, ' ergs'
    write(666,*) 'Total energy of the bomb = ', bomb_total_energy, ' ergs'
    write(6,*) 'Total energy of the model = ', total_initial_energy, ' ergs'
    write(6,*) 'Total energy of the bomb = ', bomb_total_energy, ' ergs'
    close(666)

  else
    stop "Wrong type of explosion, check the parameter 'initial_data'"
  endif
  write(6,"(A60)") "***************************************************************************"


!****************** initialize some vairables *********************************

  call compose_opacity_tables_OPAL

  call opacity(rho(:),temp(:),kappa(:),kappa_table(:),dkappadt(:))

  call optical_depth(rho(:), r(:), kappa_table(:), tau(:))

  call luminosity(r(:),temp(:),kappa(:),lambda(:),inv_kappa(:),lum(:))

  call read_BolCorr

!*********** output the initial values of some variables for analysis *********

  filename = trim(adjustl(outdir))//"/rho_initial.dat"
  call output_screenshot(rho,filename,imax)

  filename = trim(adjustl(outdir))//"/rad_initial.dat"
  call output_screenshot(r,filename,imax)

  filename = trim(adjustl(outdir))//"/mass_initial.dat"
  call output_screenshot(mass,filename,imax)

  filename = trim(adjustl(outdir))//"/delta_mass_initial.dat"
  call output_screenshot(delta_mass,filename,imax)

  filename = trim(adjustl(outdir))//"/press_initial.dat"
  call output_screenshot(p,filename,imax)

  !density as a function of distance from the surface inwards
  filename = trim(adjustl(outdir))//"/density_profile.dat"
  open(unit=666,file=trim(adjustl(filename)),status="unknown", &
          form='formatted',position="append")
  do i=1, imax-1
      write(666,"(1P20E19.10E3)") r(imax)-r(imax-i), &
          rho(imax-i),mass(imax)-mass(imax-i)
  enddo
  close(666)

end subroutine problem



