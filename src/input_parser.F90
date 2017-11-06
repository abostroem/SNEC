subroutine input_parser
!
! This routine parses the "parameters" file and sets various flags
! that are kept in the modules blmod and parameters.
!
  use blmod, only: gravity_switch, wipe_outdir
  use parameters
  implicit none

  character(len=128) :: cpstring
  character(len=500) :: rmstring
  logical :: opt
  logical :: outdirthere

!------------------------------------------------------------------------------

  opt = .false.

!****************************** LAUNCH ****************************************

  call get_string_parameter('outdir',outdir,opt)

!************************* STELLAR PROFILE ************************************


  call get_string_parameter('profile_name',profile_name,opt)
  call get_string_parameter('comp_profile_name',composition_profile_name,opt)


!***************************** EXPLOSION **************************************

  call get_string_parameter('initial_data',initial_data,opt)

  if(initial_data.eq."Piston_Explosion") then
     call get_double_parameter('piston_vel',piston_vel,opt)
     call get_double_parameter('piston_tstart',piston_tstart,opt)
     call get_double_parameter('piston_tend',piston_tend,opt)
  endif

  if(initial_data.eq."Thermal_Bomb") then
     call get_double_parameter('final_energy',final_energy,opt)
     call get_double_parameter('bomb_tstart',bomb_tstart,opt)
     call get_double_parameter('bomb_tend',bomb_tend,opt)
     call get_double_parameter('bomb_mass_spread',bomb_mass_spread,opt)
     call get_integer_parameter('bomb_start_point',bomb_start_point,opt)
     call get_integer_parameter('bomb_mode',bomb_mode,.true.)
     if(bomb_mode.eq.-666) then
        bomb_mode = 1 ! default
     endif
  endif

!******************************** GRID ****************************************

  call get_integer_parameter('imax',imax,opt)
  call get_string_parameter('gridding',gridding,opt)
  call get_logical_parameter('mass_excision',mass_excision,opt)
  if(mass_excision) then
     call get_double_parameter('mass_excised',mass_excised,opt)
  end if

!****************************** EVOLUTION *************************************

  call get_logical_parameter('radiation',radiation,opt)
  call get_integer_parameter('eoskey',eoskey,opt)
  call get_integer_parameter('Ni_switch',Ni_switch,opt)
  call get_double_parameter('Ni_mass',Ni_mass,opt)
  call get_double_parameter('Ni_boundary_mass',Ni_boundary_mass,opt)
  call get_double_parameter('Ni_period',Ni_period,opt)
  call get_integer_parameter('Ni_by_hand',Ni_by_hand,.true.)
  if (Ni_by_hand.eq.-666) Ni_by_hand = 1 ! set to default value

  call get_integer_parameter('saha_ncomps',saha_ncomps,opt)
  
  call get_logical_parameter('boxcar_smoothing',boxcar_smoothing,opt)

  call get_double_parameter('opacity_floor_envelope',of_env,opt)
  call get_double_parameter('opacity_floor_core',of_core,opt)

!********************** WHEN TO DO THINGS *************************************

  call get_integer_parameter('ntmax',ntmax,opt)
  call get_double_parameter('tend',tend,opt)
  call get_double_parameter('dtout',dtout,opt)
  call get_double_parameter('dtout_scalar',dtout_scalar,opt)
  call get_double_parameter('dtout_check',dtout_check,opt)
  call get_integer_parameter('ntout',ntout,opt)
  call get_integer_parameter('ntout_scalar',ntout_scalar,opt)
  call get_integer_parameter('ntout_check',ntout_check,opt)
  call get_integer_parameter('ntinfo',ntinfo,opt)
  call get_double_parameter('dtmin',dtmin,opt)
  call get_double_parameter('dtmax',dtmax,opt)
  
!********************************** TEST **************************************

  call get_logical_parameter('sedov',sedov,opt)
  if(sedov) then
     gravity_switch = 0
  else
     gravity_switch = 1
  end if

!******************************************************************************

! check if output directory exists
#if __INTEL_COMPILER
  inquire(directory=trim(adjustl(outdir)),exist=outdirthere)
#else
  inquire(file=trim(adjustl(outdir)),exist=outdirthere)
#endif
  if(.not.outdirthere) then
     write(6,*) "*** Output directory does not exist."
     write(6,*) "Please create the output directory: ", trim(adjustl(outdir))
     stop
  endif

! wipe output dir if requested:
  if(wipe_outdir) then
     write(*,*) "Removing output directory contents: ", trim(outdir)
     write(rmstring,*) "rm -rf ", trim(outdir), '/*'
     call system(rmstring)
  endif
  
! copy parameter file
  cpstring="cp parameters "//trim(adjustl(outdir))
  call system(cpstring)


!######## subroutines to get parse the parameters file ########################

contains
    subroutine get_string_parameter(parname,par,opt)

        implicit none
        logical opt
        character*(*) parname
        character*(*) par
        character*(200) line_string
        integer i,j,l,ll
        character*(200) temp_string

        open(unit=27,file='parameters',status='unknown')

        10 continue
        read(27,'(a)',end=19,err=10) line_string
        ! separator is an equal sign '=', # is comment
        i = index(line_string,'=')
        j = index(line_string,'#')

        if (i.eq.0.or.j.eq.1) goto 10
        !   if the whole line is a comment or there is no
        !   equal sign, then go on to the next line    

        if(j.gt.0.and.j.lt.i) goto 10
        !   if there is an equal sign, but it is in a comment
        !   then go on to the next line

        ! is this the right parameter? If not, cycle
        temp_string=trim(adjustl(line_string(1:i-1)))
        l=len(parname)
        if(parname.ne.temp_string(1:l)) goto 10

        !  If there is a comment in the line, exclude it!
        l = len(line_string)
        if (j.gt.0) l = j - 1

        par = line_string(i+1:l)
        ! now remove potential crap!
        do ll=1,len(par)
            if(par(ll:ll).eq.'\t') par(ll:ll) = ' '
            if(par(ll:ll).eq.'"') par(ll:ll) = ' '
            if(par(ll:ll).eq."'") par(ll:ll) = ' '
        enddo
        ! adjust left...
        par = adjustl(par)
        ! get rid of trailing blanks
        j = index(par," ")
        par = par(1:j-1)


        ! now look for " or ' and remove them
        j=index(par,'"')
        if(j.ne.0) stop "No quotes in my strings, please!"

        j=index(par,"'")
        if(j.ne.0) stop "No quotes in my strings, please!"

        close(27)
        return

        19 continue
        if(.not.opt) then
            write(6,*) "Fatal problem in input parser:"
            write(6,*) "Parameter ",parname
            write(6,*) "could not be read!"
            write(6,*) 
            call flush(6)
            stop
        else
            par = "NOTTHERE"
            close(27)
        endif

    end subroutine get_string_parameter

    subroutine get_double_parameter(parname,par,opt)

        implicit none
        logical opt
        character(*) parname
        character*256 line_string
        real*8 par

        call get_string_parameter(parname,line_string,opt)

        if(index(line_string,'.').eq.0) then
            write(6,*) "Uh. Bad double parameter ",trim(parname)
            write(6,*) "Please check input file!"
            call flush(6)
            stop
        endif

        read(line_string,"(e20.15)") par

    end subroutine get_double_parameter

    subroutine get_integer_parameter(parname,par,opt)

        implicit none
        logical opt
        character(*) parname
        character*256 line_string
        integer par

        call get_string_parameter(parname,line_string,opt)

        if((opt).and.trim(adjustl(line_string)) &
             .eq."NOTTHERE") then
           par = -666
        else
           read(line_string,"(i10)") par
        endif

    end subroutine get_integer_parameter

    subroutine get_logical_parameter(parname,par,opt)

        implicit none
        logical opt
        character*(*) parname
        character*(50) value_string
        integer temp_par
        logical par


        call get_string_parameter(parname,value_string,opt)

        if(opt) then
            ! don't try to set the parameter if it is not
            ! in the input file
            if(value_string .eq. "NOTTHERE") then
                write(6,*) "*** Parameter ",trim(adjustl(parname)), &
                    "not found in input file. Using default value."
                return
            endif
        endif
        read(value_string,"(i10)") temp_par

        if(temp_par.ne.0) then
            par = .true.
        else
            par = .false.
        endif

    end subroutine get_logical_parameter

end subroutine input_parser
