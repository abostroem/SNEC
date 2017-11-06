subroutine allocate_vars

  use blmod
  use parameters
  implicit none

  allocate(mass(imax))
  allocate(cmass(imax))
  allocate(rho(imax))
  allocate(rho_p(imax))
  allocate(eps(imax))
  allocate(eps_p(imax))
  allocate(temp(imax))
  allocate(temp_p(imax))
  allocate(r(imax))
  allocate(r_p(imax))
  allocate(cr(imax))
  allocate(cr_p(imax))
  allocate(vel(imax))
  allocate(vel_p(imax))
  allocate(p(imax))
  allocate(p_p(imax))
  allocate(Q(imax))
  allocate(Qterm(imax))
  allocate(ye(imax))
  allocate(abar(imax))
  allocate(cs2(imax))
  allocate(entropy(imax))
  allocate(dedt(imax))
  allocate(dpdt(imax))
  allocate(tau(imax))
  allocate(p_rad(imax))
  allocate(lum(imax))
  allocate(lambda(imax))
  allocate(kappa_table(imax))
  allocate(kappa(imax))
  allocate(kappa_p(imax))
  allocate(dkappadt(imax))
  allocate(inv_kappa(imax))
  allocate(free_electron_frac(imax))
  allocate(delta_time(imax))
  allocate(delta_mass(imax))
  allocate(delta_cmass(imax))
  allocate(metallicity(imax))
  allocate(bomb_heating(imax))

  if(ncomps.gt.0) then
     allocate(comp(imax,ncomps))
     allocate(comp_details(ncomps,2))
     allocate(zav(ncomps,imax))
     allocate(ion_fractions(ncomps,30,imax))
     comp(:,:) = 0.0d0
     comp_details(:,:) = 0.0d0
     zav(:,:) = 1.0d0
     ion_fractions(:,:,:) = 0.0d0
  endif

  allocate(photosphere_tracer(imax))

!variables, used in analysis
  allocate(E_shell(imax))
  allocate(time_diff(imax))
  allocate(time_exp(imax))

!quantities related to the heating by radioactive Ni
  allocate(Ni_deposit_function(imax))
  allocate(Ni_heating(imax))

!used for the calculations of the opacity (see opacity.F90)
  allocate(xxc(imax))
  allocate(xxo(imax))
  allocate(logR_op(imax))
  allocate(logT(imax))
  allocate(opacity_floor(imax))

end subroutine allocate_vars

!****** initialization *************
subroutine initialize_vars

    use blmod
    use parameters
    use outinfomod, only: outinfo_count
    implicit none

    mass(:)     = 0.0d0
    cmass(:)    = 0.0d0
    rho(:)      = 0.0d0
    rho_p(:)    = 0.0d0
    eps(:)      = 0.0d0
    eps_p(:)    = 0.0d0
    temp(:)     = 0.0d0
    temp_p(:)   = 0.0d0
    r(:)        = 0.0d0
    r_p(:)      = 0.0d0
    cr(:)       = 0.0d0
    cr_p(:)     = 0.0d0
    vel(:)      = 0.0d0
    vel_p(:)    = 0.0d0
    p(:)        = 0.0d0
    p_p(:)      = 0.0d0
    Q(:)        = 0.0d0
    Qterm(:)    = 0.0d0
    ye(:)       = 0.0d0
    abar(:)     = 0.0d0
    cs2(:)      = 0.0d0
    entropy(:)  = 0.0d0
    dedt(:)     = 0.0d0
    dpdt(:)     = 0.0d0
    tau(:)      = 0.0d0
    p_rad(:)    = 0.0d0
    lum(:)      = 0.0d0
    lambda(:)   = 0.0d0
    kappa_table(:) = 0.2d0
    kappa(:)       = 0.2d0
    kappa_p(:)     = 0.2d0
    dkappadt(:)    = 0.0d0
    inv_kappa(:)   = 0.0d0
    free_electron_frac(:) = 1.0d0
    delta_time(:)  = 0.0d0
    delta_mass(:)  = 0.0d0
    delta_cmass(:) = 0.0d0
    metallicity(:) = 0.0d0
    bomb_heating(:) = 0.0d0

    photosphere_tracer(:) = 0.0d0

    E_shell(:)   = 0.0d0
    time_diff(:) = 0.0d0
    time_exp(:)  = 0.0d0

    Ni_deposit_function(:) = 0.0d0
    Ni_heating(:) = 0.0d0

    xxc(:)     = 0.0d0
    xxo(:)     = 0.0d0
    logR_op(:) = 0.0d0
    logT(:)    = 0.0d0
    opacity_floor(:) = 0.0d0

    time        = 0.0d0
    dtime       = 0.0d0
    dtime_p     = 0.0d0
    nt          = 0
    time_Ni      = 0.0d0

    scratch_step = .false.

    photosphere_fell_on_the_center = 0
    Ni_contributes_five_percents = 0

    !initial shock position
    shockpos = 1
    shockpos_prev = 1
    shockpos_stop = 0
    breakoutflag = 0

    !reset counter for screen output
    outinfo_count = 0

end subroutine initialize_vars


        
