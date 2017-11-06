module blmod

 implicit none

! some infrastructure stuff
 logical :: wipe_outdir = .true.    !used in input_parser.F90
 logical :: switch1,switch2,switch3 !for output control
 
! physics controls
 logical :: scratch_step
 integer :: gravity_switch !gravity may be turned off (like for Sedov test)

!some profile parameters
 integer :: ncomps              !number of isotopes
 real*8  :: Rstar                   !radius of the model
 real*8  :: total_initial_energy    !initial energy of the model

!explosion parameters
 logical :: do_piston
 logical :: do_bomb
 real*8  :: bomb_total_energy
 integer :: bomb_spread
 real*8, allocatable :: bomb_heating(:)
 
 ! Time centering term for radiative transfer
 real*8, parameter :: theta = 0.5d0

!envelope metallicity is the nominal metallicity of the OPAL Type II tables
!changing the tables, don't forget to change this number
 real*8, parameter :: envelope_metallicity = 0.02d0
 
 !evolution vars
 real*8 :: time, dtime, dtime_p

 real*8 :: tstart
 real*8 :: tdump_check
 real*8 :: tdump,tdump_scalar

 integer :: nt !timestep
 integer :: ntstart

!_p variables contain data from previous timestep
!variables defined at cell inner boundaries
real*8, allocatable :: mass(:)                         !mass coordinate
real*8, allocatable :: delta_cmass(:)
real*8, allocatable :: r(:), r_p(:)                    !radius
real*8, allocatable :: lambda(:)                       !flux-limiter
real*8, allocatable :: lum(:)                          !luminosity
real*8, allocatable :: vel(:), vel_p(:)                !velocity
real*8, allocatable :: inv_kappa(:)                    !inverse opacity

!variables defined at cell centers
real*8, allocatable :: cmass(:)
real*8, allocatable :: delta_mass(:)
real*8, allocatable :: cr(:),   cr_p(:)
real*8, allocatable :: p(:),    p_p(:)                 !pressure
real*8, allocatable :: eps(:),  eps_p(:)               !specific internal energy
real*8, allocatable :: temp(:), temp_p(:)              !temperature
real*8, allocatable :: rho(:),  rho_p(:)               !density
real*8, allocatable :: kappa_table(:)                  !tabular opacity
real*8, allocatable :: kappa(:),kappa_p(:)             !opacity (floor applied)
real*8, allocatable :: dkappadt(:)
real*8, allocatable :: ye(:)                           !electron fraction
real*8, allocatable :: abar(:)
real*8, allocatable :: cs2(:)                          !speed of sound squared
real*8, allocatable :: comp(:,:)                       !mass fractions of elem-s
real*8, allocatable :: comp_details(:,:)               !mass and atomic numbers
real*8, allocatable :: metallicity(:)                  !metallicity

real*8, allocatable :: Q(:)                            !artificial viscosity
!term in the energy equation with artificial viscosity
real*8, allocatable :: Qterm(:)

!derivatives of energy/pressure with respect to temperature
real*8, allocatable :: dedt(:)
real*8, allocatable :: dpdt(:)

real*8, allocatable :: entropy(:)                      !entropy
real*8, allocatable :: tau(:)                          !optical depth
real*8, allocatable :: p_rad(:)                        !radiation pressure

 real*8, allocatable :: delta_time(:)

 real*8, allocatable :: zav(:,:) !average charge of heavy particles
 real*8, allocatable :: ion_fractions(:,:,:) !ionization fractions 
                                             !(ion:state:gridpoint)
 real*8, allocatable :: free_electron_frac(:)

 !quantities at the photosphere
 integer :: index_photo
 real*8 :: lum_photo, mass_photo, vel_photo, rad_photo
 integer :: photosphere_fell_on_the_center
 integer :: Ni_contributes_five_percents
 real*8, allocatable :: photosphere_tracer(:)

 integer :: index_lumshell !index of the luminosity shell
 real*8 :: mass_lumshell   !mass coordinate of the luminosity shell

 real*8 :: lum_observed    !observed luminosity

!numbers of some important elements in the composition profile
 integer :: H_number, He_number, C_number, O_number, Ni_number

!variables, used in analysis
 real*8, allocatable :: E_shell(:)    !internal energy of the shells
 real*8, allocatable :: time_diff(:)  !characteristic diffusion time
 real*8, allocatable :: time_exp(:)   !characteristic expansion time

!variables, needed for the shock capturing (shock_capture.F90)
 integer :: shockpos
 integer :: shockpos_prev
 integer :: shockpos_stop
 integer :: breakoutflag
 real*8 :: radshock
 real*8 :: radshock_prev, time_prev

!quantities needed for the color magnitude calculations (analysis.F90)
 real*8, allocatable :: bol_corr(:,:)
 real*8, allocatable :: temp_bol_corr(:)
 integer :: nlines_bol_corr
 real*8 :: magnitudes(11)
 real*8 :: T_eff

!quantities related to the heating by radioactive Ni
 real*8, allocatable :: Ni_deposit_function(:)
 real*8, allocatable :: Ni_heating(:)
 real*8 :: Ni_total_luminosity
 real*8 :: time_Ni
 real*8 :: Ni_energy_rate

!quantities used for the calculations of the opacity (see opacity.F90)
 real*8, allocatable :: xxc(:), xxo(:)
 real*8, allocatable :: logR_op(:)
 real*8, allocatable :: logT(:)
 real*8, allocatable :: opacity_floor(:) !(set in problem.F90)
 integer :: opacity_corrupted

 integer :: rpoints, tpoints
 real*8, allocatable :: opacity_tables(:,:,:) !log kappa (gridpoint:logT:LogR)
 real*8, allocatable :: logT_array(:) !values of log T in the tables
 real*8, allocatable :: logR_array(:) !values of log R in the tables
 integer, allocatable :: op_rows(:,:)
 integer, allocatable :: op_cols(:,:)

!other
 real*8 :: eos_gamma1

end module blmod

!############################# PARAMETERS MODULE ##############################
       
module parameters

   implicit none

!-------------------- launch --------------------------------

  character(len=256) :: outdir

!-------------------- profile -------------------------------

  character(len=256) :: profile_name
  character(len=256) :: composition_profile_name

!------------------- explosion ------------------------------

  character(len=256) :: initial_data

  !piston stuff
  real*8  :: piston_vel
  real*8  :: piston_tstart
  real*8  :: piston_tend

  !thermal bomb stuff
  real*8  :: final_energy
  real*8  :: bomb_tstart
  real*8  :: bomb_tend
  real*8  :: bomb_mass_spread
  integer :: bomb_start_point

  ! bomb_mode 1 -- (default) "final_energy" is asymptotic energy,
  !           2 -- select "final_energy" as bomb input energy
  integer :: bomb_mode

!-------------------- grid ------------------------------------

  integer :: imax !number of zones
  character(len=256) :: gridding

  logical :: mass_excision
  real*8  :: mass_excised

!-------------------- evolution -------------------------------

  logical :: radiation

  integer :: eoskey ! eos to use

  integer :: Ni_switch
  integer :: Ni_by_hand
  real*8 :: Ni_mass
  real*8 :: Ni_boundary_mass
  real*8 :: Ni_period

  integer :: saha_ncomps

  logical :: boxcar_smoothing

  real*8 :: of_env
  real*8 :: of_core

!--------------------- timing -----------------------------------

  integer :: ntmax

  real*8 :: tend

  real*8 :: dtout, dtout_scalar, dtout_check
  integer :: ntout, ntout_scalar, ntout_check

  integer :: ntinfo

  real*8 :: dtmin
  real*8 :: dtmax

!---------------------- test -------------------------------------

  logical :: sedov = .false.

end module parameters

!############################# EOS MODULE #####################################

module eosmodule

 implicit none

 ! array for ionization potential
 real*8 :: xxip(30,30)
 real*8,allocatable :: stat_weight_p1_ratio(:,:)

 contains
   subroutine init_ionpot
  !#############################################################################
  !
  ! These data are taken from the Timmes EOS with Saha ionization:
  !
  ! http://cococubed.asu.edu/code_pages/eos_ionize.shtml
  !
  ! For more information on the Timmes EOS see:
  !    timmes & arnett, apj supp. 125, 277, 1999
  !    timmes & swesty, apj supp. 126, 501, 2000
  !
  !#############################################################################

     use blmod, only: comp_details,ncomps
     use parameters, only: saha_ncomps
     implicit none
     integer :: k,j,i
     real*8, parameter :: ev2erg  = 1.602d-12
     real*8, parameter :: off_table = 1.0d30
     real*8  :: zion
     real*8  :: stat_weight ! function!
     integer :: izion

     xxip(:,:) = 0.0d0

#if 1
     ! set up ratio of statistical weights needed
     ! in simple_saha.F90
     allocate(stat_weight_p1_ratio(ncomps,30))
     stat_weight_p1_ratio(:,:) = 0.0d0
     do j=1,saha_ncomps
        zion = comp_details(j,2)
        izion = int(zion)
        do i=1,izion
           stat_weight_p1_ratio(j,i) = stat_weight(zion,i+1)/stat_weight(zion,i)
        enddo
     enddo
#endif

     !hydrogen
     xxip(1,1)    = 13.59844d0
     xxip(1,2:30) = off_table

     !helium
     xxip(2,1)    = 24.58741d0
     xxip(2,2)    = 54.41778d0
     xxip(2,3:30) = off_table

     !lithium
     xxip(3,1)    = 5.39172d0
     xxip(3,2)    = 75.64018d0
     xxip(3,3)    = 122.45429d0
     xxip(3,4:30) = off_table

     !berrylium
     xxip(4,1)    = 9.3227d0
     xxip(4,2)    = 18.21116d0
     xxip(4,3)    = 153.89661d0
     xxip(4,4)    = 217.71865d0
     xxip(4,5:30) = off_table

     !boron
     xxip(5,1)    = 8.29803d0
     xxip(5,2)    = 25.15484d0
     xxip(5,3)    = 37.93064d0
     xxip(5,4)    = 259.37521d0
     xxip(5,5)    = 340.22580d0
     xxip(5,6:30) = off_table

     !carbon
     xxip(6,1)    = 11.26030d0
     xxip(6,2)    = 24.38332d0
     xxip(6,3)    = 47.8878d0
     xxip(6,4)    = 64.4939d0
     xxip(6,5)    = 392.087d0
     xxip(6,6)    = 489.99334d0
     xxip(6,7:30) = off_table

     !nitrogen
     xxip(7,1)    = 14.53414d0
     xxip(7,2)    = 29.6013d0
     xxip(7,3)    = 47.44924d0
     xxip(7,4)    = 77.4735d0
     xxip(7,5)    = 97.8902d0
     xxip(7,6)    = 552.0718d0
     xxip(7,7)    = 667.046d0
     xxip(7,8:30) = off_table

     !oxygen
     xxip(8,1)    = 13.61806d0
     xxip(8,2)    = 35.11730d0
     xxip(8,3)    = 54.9355d0
     xxip(8,4)    = 77.41353d0
     xxip(8,5)    = 113.8990d0
     xxip(8,6)    = 138.1197d0
     xxip(8,7)    = 739.29d0
     xxip(8,8)    = 871.4101d0
     xxip(8,9:30) = off_table

     !fluorine
     xxip(9,1)      = 17.42282d0
     xxip(9,2)      = 34.97082d0
     xxip(9,3)      = 62.7084d0
     xxip(9,4)      = 87.1398d0
     xxip(9,5)      = 114.2428d0
     xxip(9,6)      = 157.1651d0
     xxip(9,7)      = 185.18d0
     xxip(9,8)      = 953.9112d0
     xxip(9,9)      = 1103.1176d0
     xxip(9,10:30)  = off_table

     !neon
     xxip(10,1)      = 21.5646d0
     xxip(10,2)      = 40.96328d0
     xxip(10,3)      = 63.45d0
     xxip(10,4)      = 97.12d0
     xxip(10,5)      = 126.21d0
     xxip(10,6)      = 157.93d0
     xxip(10,7)      = 207.2759d0
     xxip(10,8)      = 239.0989d0
     xxip(10,9)      = 1195.8286d0
     xxip(10,10)     = 1362.1995d0
     xxip(10,11:30)  = off_table

      !sodium
      xxip(11,1)      = 5.13908d0
      xxip(11,2)      = 47.2864d0
      xxip(11,3)      = 71.6200d0
      xxip(11,4)      = 98.91d0
      xxip(11,5)      = 138.40d0
      xxip(11,6)      = 172.18d0
      xxip(11,7)      = 208.50d0
      xxip(11,8)      = 264.25d0
      xxip(11,9)      = 299.864d0
      xxip(11,10)     = 1465.121d0
      xxip(11,11)     = 1648.702d0
      xxip(11,12:30)  = off_table

      !magnesium
      xxip(12,1)      = 7.64624d0
      xxip(12,2)      = 15.03528d0
      xxip(12,3)      = 80.1437d0
      xxip(12,4)      = 109.2655d0
      xxip(12,5)      = 141.27d0
      xxip(12,6)      = 186.76d0
      xxip(12,7)      = 225.02d0
      xxip(12,8)      = 265.96d0
      xxip(12,9)      = 328.06d0
      xxip(12,10)     = 367.50d0
      xxip(12,11)     = 1761.805d0
      xxip(12,12)     = 1962.6650d0
      xxip(12,13:30)  = off_table

      !aluminum
      xxip(13,1)      = 5.98577d0
      xxip(13,2)      = 18.82856d0
      xxip(13,3)      = 28.44765d0
      xxip(13,4)      = 119.992d0
      xxip(13,5)      = 153.825d0
      xxip(13,6)      = 190.49d0
      xxip(13,7)      = 241.76d0
      xxip(13,8)      = 284.66d0
      xxip(13,9)      = 330.13d0
      xxip(13,10)     = 398.75d0
      xxip(13,11)     = 442.00d0
      xxip(13,12)     = 2085.98d0
      xxip(13,13)     = 2304.1410d0
      xxip(13,14:30)  = off_table

      !silicon
      xxip(14,1)      = 8.15169d0
      xxip(14,2)      = 16.34585d0
      xxip(14,3)      = 33.49302d0
      xxip(14,4)      = 45.14181d0
      xxip(14,5)      = 166.767d0
      xxip(14,6)      = 205.27d0
      xxip(14,7)      = 246.5d0
      xxip(14,8)      = 303.54d0
      xxip(14,9)      = 351.12d0
      xxip(14,10)     = 401.37d0
      xxip(14,11)     = 476.36d0
      xxip(14,12)     = 523.42d0
      xxip(14,13)     = 2437.63d0
      xxip(14,14)     = 2673.182d0
      xxip(14,15:30)  = off_table

      !phosphorous
      xxip(15,1)      = 10.48669d0
      xxip(15,2)      = 19.7694d0
      xxip(15,3)      = 30.2027d0
      xxip(15,4)      = 51.4439d0
      xxip(15,5)      = 65.0251d0
      xxip(15,6)      = 220.421d0
      xxip(15,7)      = 263.57d0
      xxip(15,8)      = 309.60d0
      xxip(15,9)      = 372.13d0
      xxip(15,10)     = 424.4d0
      xxip(15,11)     = 479.46d0
      xxip(15,12)     = 560.8d0
      xxip(15,13)     = 611.74d0
      xxip(15,14)     = 2816.91d0
      xxip(15,15)     = 3069.842d0
      xxip(15,16:30)  = off_table

      !sulfur
      xxip(16,1)      = 10.36001d0
      xxip(16,2)      = 23.3379d0
      xxip(16,3)      = 34.79d0
      xxip(16,4)      = 47.222d0
      xxip(16,5)      = 72.5945d0
      xxip(16,6)      = 88.0530d0
      xxip(16,7)      = 280.948d0
      xxip(16,8)      = 328.7d0
      xxip(16,9)      = 379.55d0
      xxip(16,10)     = 447.5d0
      xxip(16,11)     = 504.8d0
      xxip(16,12)     = 564.44d0
      xxip(16,13)     = 652.2d0
      xxip(16,14)     = 707.01d0
      xxip(16,15)     = 3223.78d0
      xxip(16,16)     = 3494.1892d0
      xxip(16,17:30)  = off_table

      !chlorine
      xxip(17,1)      = 12.96764d0
      xxip(17,2)      = 23.814d0
      xxip(17,3)      = 39.61d0
      xxip(17,4)      = 53.4652d0
      xxip(17,5)      = 67.8d0
      xxip(17,6)      = 97.03d0
      xxip(17,7)      = 114.1958d0
      xxip(17,8)      = 348.28d0
      xxip(17,9)      = 400.06d0
      xxip(17,10)     = 455.63d0
      xxip(17,11)     = 529.28d0
      xxip(17,12)     = 591.99d0
      xxip(17,13)     = 656.71d0
      xxip(17,14)     = 749.76d0
      xxip(17,15)     = 809.40d0
      xxip(17,16)     = 3658.521d0
      xxip(17,17)     = 3946.2960d0
      xxip(17,18:30)  = off_table

      !argon
      xxip(18,1)      = 15.75962d0
      xxip(18,2)      = 27.62967d0
      xxip(18,3)      = 40.74d0
      xxip(18,4)      = 59.81d0
      xxip(18,5)      = 75.02d0
      xxip(18,6)      = 91.009d0
      xxip(18,7)      = 124.323d0
      xxip(18,8)      = 143.460d0
      xxip(18,9)      = 422.45d0
      xxip(18,10)     = 478.69d0
      xxip(18,11)     = 538.96d0
      xxip(18,12)     = 618.26d0
      xxip(18,13)     = 686.10d0
      xxip(18,14)     = 755.74d0
      xxip(18,15)     = 854.77d0
      xxip(18,16)     = 918.03d0
      xxip(18,17)     = 4120.8857d0
      xxip(18,18)     = 4426.2296d0
      xxip(18,19:30)  = off_table

      !pottasium
      xxip(19,1)      = 4.34066d0
      xxip(19,2)      = 31.63d0
      xxip(19,3)      = 45.806d0
      xxip(19,4)      = 60.91d0
      xxip(19,5)      = 82.66d0
      xxip(19,6)      = 99.4d0
      xxip(19,7)      = 117.56d0
      xxip(19,8)      = 154.88d0
      xxip(19,9)      = 175.8174d0
      xxip(19,10)     = 503.8d0
      xxip(19,11)     = 564.7d0
      xxip(19,12)     = 629.4d0
      xxip(19,13)     = 714.6d0
      xxip(19,14)     = 786.6d0
      xxip(19,15)     = 861.1d0
      xxip(19,16)     = 968.0d0
      xxip(19,17)     = 1033.4d0
      xxip(19,18)     = 4610.8d0
      xxip(19,19)     = 4934.046d0
      xxip(19,20:30)  = off_table

      !calcium
      xxip(20,1)      = 6.11316d0
      xxip(20,2)      = 11.87172d0
      xxip(20,3)      = 50.9131d0
      xxip(20,4)      = 67.27d0
      xxip(20,5)      = 84.50d0
      xxip(20,6)      = 108.78d0
      xxip(20,7)      = 127.2d0
      xxip(20,8)      = 147.24d0
      xxip(20,9)      = 188.54d0
      xxip(20,10)     = 211.275d0
      xxip(20,11)     = 591.9d0
      xxip(20,12)     = 657.2d0
      xxip(20,13)     = 726.6d0
      xxip(20,14)     = 817.6d0
      xxip(20,15)     = 894.5d0
      xxip(20,16)     = 974.0d0
      xxip(20,17)     = 1087.0d0
      xxip(20,18)     = 1157.8d0
      xxip(20,19)     = 5128.8d0
      xxip(20,20)     = 5469.864d0
      xxip(20,21:30)  = off_table

      !scandium
      xxip(21,1)      = 6.5615d0
      xxip(21,2)      = 12.79967d0
      xxip(21,3)      = 24.75666d0
      xxip(21,4)      = 73.4894d0
      xxip(21,5)      = 91.65d0
      xxip(21,6)      = 110.68d0
      xxip(21,7)      = 138.0d0
      xxip(21,8)      = 158.1d0
      xxip(21,9)      = 180.03d0
      xxip(21,10)     = 225.18d0
      xxip(21,11)     = 249.798d0
      xxip(21,12)     = 687.36d0
      xxip(21,13)     = 756.7d0
      xxip(21,14)     = 830.8d0
      xxip(21,15)     = 927.5d0
      xxip(21,16)     = 1009.0d0
      xxip(21,17)     = 1094.0d0
      xxip(21,18)     = 1213.0d0
      xxip(21,19)     = 1287.97d0
      xxip(21,20)     = 5674.8d0
      xxip(21,21)     = 6033.712d0
      xxip(21,22:30)  = off_table

      !titanium
      xxip(22,1)      = 6.8281d0
      xxip(22,2)      = 13.5755d0
      xxip(22,3)      = 27.4917d0
      xxip(22,4)      = 43.2672d0
      xxip(22,5)      = 99.30d0
      xxip(22,6)      = 119.53d0
      xxip(22,7)      = 140.8d0
      xxip(22,8)      = 170.4d0
      xxip(22,9)      = 192.1d0
      xxip(22,10)     = 215.92d0
      xxip(22,11)     = 265.07d0
      xxip(22,12)     = 291.500d0
      xxip(22,13)     = 787.84d0
      xxip(22,14)     = 863.1d0
      xxip(22,15)     = 941.9d0
      xxip(22,16)     = 1044.0d0
      xxip(22,17)     = 1131.0d0
      xxip(22,18)     = 1221.0d0
      xxip(22,19)     = 1346.0d0
      xxip(22,20)     = 1425.4d0
      xxip(22,21)     = 6249.0d0
      xxip(22,22)     = 6625.82d0
      xxip(22,24:30)  = off_table

      !vanadium
      xxip(23,1)      = 6.7463d0
      xxip(23,2)      = 14.66d0
      xxip(23,3)      = 29.311d0
      xxip(23,4)      = 46.709d0
      xxip(23,5)      = 65.2817d0
      xxip(23,6)      = 128.13d0
      xxip(23,7)      = 150.6d0
      xxip(23,8)      = 173.4d0
      xxip(23,9)      = 205.8d0
      xxip(23,10)     = 230.5d0
      xxip(23,11)     = 255.7d0
      xxip(23,12)     = 308.1d0
      xxip(23,13)     = 336.277d0
      xxip(23,14)     = 896.0d0
      xxip(23,15)     = 976.0d0
      xxip(23,16)     = 1060.0d0
      xxip(23,17)     = 1168.0d0
      xxip(23,18)     = 1260.0d0
      xxip(23,19)     = 1355.0d0
      xxip(23,20)     = 1486.0d0
      xxip(23,21)     = 1569.6d0
      xxip(23,22)     = 6851.3d0
      xxip(23,23)     = 7246.12d0
      xxip(23,24:30)  = off_table

      !chromium
      xxip(24,1)      = 6.7665d0
      xxip(24,2)      = 16.4857d0
      xxip(24,3)      = 30.96d0
      xxip(24,4)      = 49.16d0
      xxip(24,5)      = 69.46d0
      xxip(24,6)      = 90.6349d0
      xxip(24,7)      = 160.18d0
      xxip(24,8)      = 184.7d0
      xxip(24,9)      = 209.3d0
      xxip(24,10)     = 244.4d0
      xxip(24,11)     = 270.8d0
      xxip(24,12)     = 298.0d0
      xxip(24,13)     = 354.8d0
      xxip(24,14)     = 384.168d0
      xxip(24,15)     = 1010.6d0
      xxip(24,16)     = 1097.0d0
      xxip(24,17)     = 1185.0d0
      xxip(24,18)     = 1299.0d0
      xxip(24,19)     = 1396.0d0
      xxip(24,20)     = 1496.0d0
      xxip(24,21)     = 1634.0d0
      xxip(24,22)     = 1721.4d0
      xxip(24,23)     = 7481.7d0
      xxip(24,24)     = 7894.81d0
      xxip(24,25:30)  = off_table

      !manganese
      xxip(25,1)      = 7.43402d0
      xxip(25,2)      = 15.63999d0
      xxip(25,3)      = 33.668d0
      xxip(25,4)      = 51.2d0
      xxip(25,5)      = 72.4d0
      xxip(25,6)      = 95.6d0
      xxip(25,7)      = 119.203d0
      xxip(25,8)      = 194.5d0
      xxip(25,9)      = 221.8d0
      xxip(25,10)     = 248.3d0
      xxip(25,11)     = 286.0d0
      xxip(25,12)     = 314.4d0
      xxip(25,13)     = 343.6d0
      xxip(25,14)     = 403.0d0
      xxip(25,15)     = 435.163d0
      xxip(25,16)     = 1134.7d0
      xxip(25,17)     = 1224.0d0
      xxip(25,18)     = 1317.0d0
      xxip(25,19)     = 1437.0d0
      xxip(25,20)     = 1539.0d0
      xxip(25,21)     = 1644.0d0
      xxip(25,22)     = 1788.0d0
      xxip(25,23)     = 1879.9d0
      xxip(25,24)     = 8140.6d0
      xxip(25,25)     = 8571.94d0
      xxip(25,26:30)  = off_table

      !iron
      xxip(26,1)      = 7.9024d0
      xxip(26,2)      = 16.1878d0
      xxip(26,3)      = 30.652d0
      xxip(26,4)      = 54.8d0
      xxip(26,5)      = 75.0d0
      xxip(26,6)      = 99.1d0
      xxip(26,7)      = 124.98d0
      xxip(26,8)      = 151.06d0
      xxip(26,9)      = 233.6d0
      xxip(26,10)     = 262.1d0
      xxip(26,11)     = 290.2d0
      xxip(26,12)     = 330.8d0
      xxip(26,13)     = 361.0d0
      xxip(26,14)     = 392.2d0
      xxip(26,15)     = 457.0d0
      xxip(26,16)     = 489.256d0
      xxip(26,17)     = 1266.0d0
      xxip(26,18)     = 1358.0d0
      xxip(26,19)     = 1456.0d0
      xxip(26,20)     = 1582.0d0
      xxip(26,21)     = 1689.0d0
      xxip(26,22)     = 1799.0d0
      xxip(26,23)     = 1950.0d0
      xxip(26,24)     = 2023.0d0
      xxip(26,25)     = 8828.0d0
      xxip(26,26)     = 9277.69d0
      xxip(26,27:30)  = off_table

      !cobalt
      xxip(27,1)      = 7.8810d0
      xxip(27,2)      = 17.083d0
      xxip(27,3)      = 33.50d0
      xxip(27,4)      = 51.3d0
      xxip(27,5)      = 79.5d0
      xxip(27,6)      = 102.0d0
      xxip(27,7)      = 128.9d0
      xxip(27,8)      = 157.8d0
      xxip(27,9)      = 186.13d0
      xxip(27,10)     = 275.4d0
      xxip(27,11)     = 305.0d0
      xxip(27,12)     = 336.0d0
      xxip(27,13)     = 379.0d0
      xxip(27,14)     = 411.0d0
      xxip(27,15)     = 444.0d0
      xxip(27,16)     = 511.96d0
      xxip(27,17)     = 546.58d0
      xxip(27,18)     = 1397.2d0
      xxip(27,19)     = 1504.6d0
      xxip(27,20)     = 1603.0d0
      xxip(27,21)     = 1735.0d0
      xxip(27,22)     = 1846.0d0
      xxip(27,23)     = 1962.0d0
      xxip(27,24)     = 2119.0d0
      xxip(27,25)     = 2219.0d0
      xxip(27,26)     = 9544.1d0
      xxip(27,27)     = 10012.12d0
      xxip(27,28:30)  = off_table

      !nickel
      xxip(28,1)      = 7.6398d0
      xxip(28,2)      = 18.16884d0
      xxip(28,3)      = 35.19d0
      xxip(28,4)      = 54.9d0
      xxip(28,5)      = 76.06d0
      xxip(28,6)      = 108.0d0
      xxip(28,7)      = 133.0d0
      xxip(28,8)      = 162.0d0
      xxip(28,9)      = 193.0d0
      xxip(28,10)     = 224.6d0
      xxip(28,11)     = 321.0d0
      xxip(28,12)     = 352.0d0
      xxip(28,13)     = 384.0d0
      xxip(28,14)     = 430.0d0
      xxip(28,15)     = 464.0d0
      xxip(28,16)     = 499.0d0
      xxip(28,17)     = 571.08d0
      xxip(28,18)     = 607.06d0
      xxip(28,19)     = 1541.0d0
      xxip(28,20)     = 1648.0d0
      xxip(28,21)     = 1756.0d0
      xxip(28,22)     = 1894.0d0
      xxip(28,23)     = 2011.0d0
      xxip(28,24)     = 2131.0d0
      xxip(28,25)     = 2295.0d0
      xxip(28,26)     = 2399.2d0
      xxip(28,27)     = 10288.8d0
      xxip(28,28)     = 10775.40d0
      xxip(28,29:30)  = off_table

      !copper
      xxip(29,1)     = 7.72638d0
      xxip(29,2)     = 20.29240d0
      xxip(29,3)     = 36.841d0
      xxip(29,4)     = 57.38d0
      xxip(29,5)     = 79.8d0
      xxip(29,6)     = 103.0d0
      xxip(29,7)     = 139.0d0
      xxip(29,8)     = 166.0d0
      xxip(29,9)     = 199.0d0
      xxip(29,10)    = 232.0d0
      xxip(29,11)    = 265.3d0
      xxip(29,12)    = 369.0d0
      xxip(29,13)    = 401.0d0
      xxip(29,14)    = 435.0d0
      xxip(29,15)    = 484.0d0
      xxip(29,16)    = 520.0d0
      xxip(29,17)    = 557.0d0
      xxip(29,18)    = 633.0d0
      xxip(29,19)    = 670.588d0
      xxip(29,20)    = 1697.0d0
      xxip(29,21)    = 1804.0d0
      xxip(29,22)    = 1916.0d0
      xxip(29,23)    = 2060.0d0
      xxip(29,24)    = 2182.0d0
      xxip(29,25)    = 2308.0d0
      xxip(29,26)    = 2478.0d0
      xxip(29,27)    = 2587.5d0
      xxip(29,28)    = 11062.38d0
      xxip(29,29)    = 11567.617d0
      xxip(29,30)    = off_table

      !zinc
      !the last nine are not in the cited reference
      xxip(30,1)     = 9.3942d0
      xxip(30,2)     = 17.96440d0
      xxip(30,3)     = 39.723d0
      xxip(30,4)     = 59.4d0
      xxip(30,5)     = 82.6d0
      xxip(30,6)     = 108.0d0
      xxip(30,7)     = 134.0d0
      xxip(30,8)     = 174.0d0
      xxip(30,9)     = 203.0d0
      xxip(30,10)    = 238.0d0
      xxip(30,11)    = 274.0d0
      xxip(30,12)    = 310.8d0
      xxip(30,13)    = 419.7d0
      xxip(30,14)    = 454.0d0
      xxip(30,15)    = 490.0d0
      xxip(30,16)    = 542.0d0
      xxip(30,17)    = 579.0d0
      xxip(30,18)    = 619.0d0
      xxip(30,19)    = 698.0d0
      xxip(30,20)    = 738.0d0
      xxip(30,21)    = 1856.0d0
      xxip(30,22)    = 1920.0d0
      xxip(30,23)    = 2075.0d0
      xxip(30,24)    = 2190.0d0
      xxip(30,25)    = 2350.0d0
      xxip(30,26)    = 2500.0d0
      xxip(30,27)    = 2600.0d0
      xxip(30,28)    = 11080.0d0
      xxip(30,29)    = 11600.0d0
      xxip(30,30)    = 12000.0d0

      do k=1,30
         do j=1,30
            if (xxip(j,k) .ne. off_table) xxip(j,k) = xxip(j,k) * ev2erg
         enddo
      enddo


   end subroutine init_ionpot

end module eosmodule


!######################## PHYSICAL CONSTANTS MODULE ###########################

module physical_constants

   implicit none
   
   real*8, parameter :: msun = 1.98892d33              !solar mass
   real*8, parameter :: rsun = 6.96d10                 !solar radius
   real*8, parameter :: clite = 2.99792458d10          !speed of light
   real*8, parameter :: ggrav = 6.6742d-8              !gravitational constant
   real*8, parameter :: kboltz = 1.380662d-16
   real*8, parameter :: mev_to_erg = 1.6022d-6
   real*8, parameter :: pi = 3.14159265358979d0 
   real*8, parameter :: emev = 1.60219d-6
   real*8, parameter :: avo_real = 6.0221415d23
   real*8, parameter :: h_cgs = 6.626058d-27
   real*8, parameter :: a_rad = 7.5657d-15
   real*8, parameter :: mproton = 1.67262178d-24
   real*8, parameter :: melectron = 9.1093897d-28
   real*8, parameter :: kdnr = 9.91d12
   real*8, parameter :: kdr = 1.231d15
   real*8, parameter :: sigma_SB = 5.6704d-5
   real*8, parameter :: tau_Ni = 760320.0d0 !8.8 days in seconds
   real*8, parameter :: tau_Co = 9616320.0d0    !111.3 days in seconds
   real*8, parameter :: overtau_Ni = 1.0d0/760320.0d0
   real*8, parameter :: overtau_Co = 1.0d0/9616320.0d0
   !relative solar mass fraction of C (G&N'93)
   real*8, parameter :: C_frac_sol = 0.173285d0 
   !relative solar mass fraction of O (G&N'93)
   real*8, parameter :: O_frac_sol = 0.482273d0 
   !absolute bolometric magnitude of sun
   real*8, parameter :: sun_mag = 4.75d0        
   !bolometric luminosity of sun in erg/s
   real*8, parameter :: sun_lum = 3.846d33

   ! constants used in Saha solver
   real*8, parameter :: saha_coeff = 2.0d0*(2.0d0*pi*melectron*kboltz/(h_cgs**2))**(1.5d0)

   ! helpers
   real*8, parameter :: fivethirds = 5.0d0/3.0d0
   real*8, parameter :: fourthirds = 4.0d0/3.0d0
   real*8, parameter :: overthree  = 1.0d0/3.0d0

end module physical_constants
