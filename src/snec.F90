program snec

  use blmod, only: dtime, dtime_p, time, nt, ntstart, tstart,   &
     tdump, tdump_scalar, rho, tdump_check
  use parameters
  use outinfomod, only: outinfo_count
  implicit none

  logical :: OutputFlag = .false.
  logical :: OutputFlagScalar = .false.
  logical :: OutputFlagCheck = .false.

!------------------------------------------------------------------------------

  write(*,*)
  
  write(*,*) "***********************************"
  write(*,*) "* Supernova Explosion Code (SNEC) *"
  write(*,*) "***********************************"

  write(*,*)
! *****************************************************
! INITIALIZATION
! *****************************************************        

  call input_parser

  call problem

  call artificial_viscosity

! output before first timestep
  call output_all(0)
  call output_all(1)
  call output_all(2)
 
  call timestep

  tdump_check = tstart+dtout_check
  tdump_scalar = tstart+dtout_scalar
  tdump = tstart+dtout
  time = tstart
  nt = ntstart
  

! *****************************************************
! MAIN LOOP
! *****************************************************


  IntegrationLoop: do
     
     dtime_p = dtime
     ! determine dt
     call timestep
     

     if(ntinfo.gt.0) then
        if(mod(nt,ntinfo).eq.0) then
           ! print useful info to stdout
           call outinfo
        endif
     endif
     
     if((time+dtime).gt.tend) dtime = tend-time

     ! actual integration step
     call blstep
     
     ! increment timestep
     nt = nt + 1
 
     ! various output related things
     if (ntout.gt.0) then
        if ( mod(nt,ntout) .eq. 0) OutputFlag = .true.
     endif
   
     if (ntout_scalar.gt.0) then
        if ( mod(nt,ntout_scalar) .eq. 0 ) OutputFlagScalar = .true.
     endif

     if (ntout_check.gt.0) then
        if ( mod(nt,ntout_check) .eq. 0 ) OutputFlagCheck = .true.
     endif
     
     if ( time.ge.tdump) then
        tdump=tdump+dtout
        OutputFlag = .true.
     endif

     if ( time.ge.tdump_scalar) then
        tdump_scalar=tdump_scalar+dtout_scalar
        OutputFlagScalar = .true.
     endif

     if ( time.ge.tdump_check) then
        tdump_check=tdump_check+dtout_check
        OutputFlagCheck = .true.
     endif
     
     ! increment time
     time = time+dtime
     
     if (OutputFlag) then
        call output_all(0)
        call output_all(1)
        OutputFlag = .false.
     endif
     
     if (OutputFlagScalar) then
        call output_all(2)
        OutputFlagScalar = .false.
     endif
     
     if (OutputFlagCheck) then
        OutputFlagCheck = .false.
     endif

     if (time.eq.tend) then
        write(*,*) "Done! :-) tend reached"
        call output_all(0)
        call output_all(1)
        call output_all(2)
        exit
     endif

     if (nt.ge.ntmax) then
        write(*,*) "Done! :-) ntmax reached"
        call output_all(0)
        call output_all(1)
        call output_all(2)
        exit
     endif
   
  enddo IntegrationLoop

end program snec
