subroutine read_profile_compositions(prof_name)

  use blmod, only: comp, comp_details, cmass, H_number, He_number, C_number, &
                    O_number, Ni_number, mass, ye, abar, metallicity, ncomps
  use parameters
  use physical_constants
  implicit none

!input:
  character(*) :: prof_name

!local:
  real*8,allocatable :: pmass(:),pradius(:),pcomp(:,:)
  real*8 :: buffer
  real*8 :: Ni_mass_fraction
  real*8 :: sum_initial

  integer :: profile_zones
  integer :: i,l
  integer :: ibuffer
  integer :: Ni_boundary_index

  character(len=1024) :: filename


!------------------------------------------------------------------------------

  open(666,file=trim(prof_name),status='unknown',form='formatted',action='read')

  read(666,*) profile_zones,ibuffer
  profile_zones = profile_zones
  write(*,*) "We have ",profile_zones, "composition profile zones."

!------------ read the composition profile and map it on the grid -------------

  allocate(pmass(profile_zones))
  allocate(pradius(profile_zones))
  allocate(pcomp(profile_zones,ncomps))

  !read in A's
  read(666,*) comp_details(1:ncomps,1)

  !read in Z's
  read(666,*) comp_details(1:ncomps,2)

  do i=1,profile_zones
     read(666,*) pmass(i),pradius(i),pcomp(i,1:ncomps)
  enddo
  close(666)

  !mass fractions live at the cell centers
  do i=1,imax-1
      do l=1,ncomps
        call map_map(comp(i,l), cmass(i), pcomp(:,l), pmass, profile_zones)
     enddo
  enddo
  do l=1,ncomps
      comp(imax,l) = comp(imax-1,l)
  enddo

  deallocate(pmass)
  deallocate(pradius)
  deallocate(pcomp)

!------------------------------------------------------------------------------

  !extract the numbers of some important elements in the composition profile
  H_number = 0
  He_number = 0
  C_number = 0
  O_number = 0
  Ni_number = 0

  do l=1,ncomps 
  if(comp_details(l,2).eq.1.0d0 .and. comp_details(l,1).eq.1.0d0) then
      H_number = l  !hydrogen
  else if(comp_details(l,2).eq.2.0d0 .and. comp_details(l,1).eq.4.0d0) then
      He_number = l !helium
  else if(comp_details(l,2).eq.6.0d0 .and. comp_details(l,1).eq.12.0d0) then
      C_number = l  !carbon
  else if(comp_details(l,2).eq.8.0d0 .and. comp_details(l,1).eq.16.0d0) then
      O_number = l  !oxygen
  else if(comp_details(l,2).eq.28.0d0 .and. comp_details(l,1).eq.56.0d0) then
      Ni_number = l !radioactive nickel
  end if
  enddo

  if(Ni_switch.eq.1 .and. Ni_number.eq.0) then
      write(*,*) 'radioactive Ni is absent in the composition profile'
      write(*,*) 'please, add a column for it (see documentation of the code)'
      write(*,*) 'or put Ni_switch = 0'
      stop
  endif

  !seed Ni by hand as a step function before renormalization and boxcar
  if (Ni_by_hand.ne.0) then
    Ni_mass_fraction = Ni_mass/(Ni_boundary_mass-mass(1)/msun)

    call map_find_index(imax,mass,Ni_boundary_mass*msun,Ni_boundary_index,ibuffer)

    if(Ni_number.ne.0) then
        do i=1, imax
            if(i.lt.Ni_boundary_index) then
                comp(i,Ni_number) = Ni_mass_fraction
            else
                comp(i,Ni_number) = 0.0d0
            end if
        end do
    end if
  end if

  !normalize composition sum to 1
  do i=1,imax
      if(Ni_number.eq.0) then
          comp(i,1:ncomps) = comp(i,1:ncomps)/sum(comp(i,1:ncomps))
      else !the fraction of Ni shouldn't change
          sum_initial = sum(comp(i,1:ncomps))
          do l=1, ncomps
              if(l.ne.Ni_number) then
                  comp(i,l) = comp(i,l)*(1-comp(i,Ni_number)) &
                      /(sum_initial-comp(i,Ni_number))
              endif
          enddo
      endif
  enddo

  !apply boxcar smoothing
  if(boxcar_smoothing) then
      call boxcar
  endif


  !set Y_e and abar
  do i=1,imax
      ye(i) = sum(comp(i,1:ncomps)* &
          (comp_details(1:ncomps,2)/comp_details(1:ncomps,1)))
      abar(i) = 1.0d0/sum(comp(i,1:ncomps)/comp_details(1:ncomps,1))
  enddo


  !output the initial fractions of some elements
  do l=1,ncomps
  if(l.eq.H_number) then !hydrogen
      filename = trim(adjustl(outdir))//"/H_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.He_number) then !helium
      filename = trim(adjustl(outdir))//"/He_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.C_number) then !carbon
      filename = trim(adjustl(outdir))//"/C_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.O_number) then !oxygen
      filename = trim(adjustl(outdir))//"/O_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  else if(l.eq.Ni_number) then !radioactive Ni
      filename = trim(adjustl(outdir))//"/Ni_init_frac.dat"
      call output_screenshot(comp(:,l),filename,imax)
  end if
  enddo


  !calculate metallicity as 1 - He_fraction - H_fraction, output initial value
  do i=1,imax
      if(H_number.ne.0 .and. He_number.ne.0) then
          metallicity(i) = 1 - comp(i,H_number) - comp(i,He_number)
      else if(H_number.ne.0 .and. He_number.eq.0) then
          metallicity(i) = 1 - comp(i,H_number)
      else if(H_number.eq.0 .and. He_number.ne.0) then
          metallicity(i) = 1 - comp(i,He_number)
      else
          metallicity(i) = 1
      end if
  enddo

  filename = trim(adjustl(outdir))//"/metallicity_init.dat"
  call output_screenshot(metallicity,filename,imax)


end subroutine read_profile_compositions
