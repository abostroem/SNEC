!##############################################################################
subroutine ideal_eos(rhox,tempx,px,ex,cs2x,dpdtx,dedtx,np)

  use blmod, only: eos_gamma1
  use physical_constants
  implicit none

!Input:
  integer :: np
  real*8 :: rhox(np),tempx(np)

!Output:
  real*8 :: px(np),ex(np),cs2x(np),dpdtx(np),dedtx(np)

!Local:
  integer :: i

!------------------------------------------------------------------------------

  do i=1, np

    px(i) = (kboltz*avo_real)*rhox(i)*tempx(i)

    ex(i) = (kboltz*avo_real)/(eos_gamma1-1.0d0)*tempx(i)

    cs2x(i) = eos_gamma1*(kboltz*avo_real)*tempx(i)

    dpdtx(i) = (kboltz*avo_real)*rhox(i)

    dedtx(i) = (kboltz*avo_real)/(eos_gamma1-1.0d0)

  end do

end subroutine ideal_eos


!##############################################################################
subroutine paczynski_eos(rhox,tempx,yex,abarx,px,ex,cs2x,dpdtx,dedtx, &
    pradx,k,np)

  use blmod, only: free_electron_frac, ion_fractions, comp_details, comp, ncomps
  use parameters
  use physical_constants
  use eosmodule, only: xxip
  implicit none

!input:
  integer :: np
  integer :: k(np)
  real*8 :: rhox(np),tempx(np),yex(np),abarx(np)

!output:
  real*8 :: px(np),ex(np),cs2x(np),dpdtx(np),dedtx(np),pradx(np)

!local:
  real*8 :: N     !Number of atoms and ions per unit mass
  real*8 :: ybar  !Number of free electrons / Total number of atoms and ions
  real*8 :: pion, pe, ped, pend, pednr, pedr   !pressure components
  real*8 :: ne                                 !free electron concentration

  real*8 :: invrhox

  real*8 :: e_ioncorr  !Ionization energy
  real*8 :: ffactor !Factor to switch between degenerate and nondegenerate cases

  real*8 :: sum1, sum2, sum3, sum_pot, sum_pot_frac
  real*8 :: ymax
  real*8 :: zion
  real*8 :: stat_weight
  real*8 :: y_r, chi_r, chi_T, chi_rho, c_V
  real*8 :: nu_j
  real*8 :: Gamma1

  integer :: i,j
  integer :: l,lmax


!------------------------------------------------------------------------------
!variables ybar, N, y_r, chi_r, chi_T, chi_rho, c_V, nu_j, Gamma1 are named as 
!the analogous quantities in paragraph 9.18 of "Cox & Giuli's principles of
!stellar structure" by Weiss et al. (2004). However, the quantities presented
!here are obtained without assuming that only one element at a time changes
!ionization state.
!The degenerate electron gas is treated as 
!in Paczynski B., ApJ, 267:315 (1983).

  do i = 1, np

    call simple_saha(saha_ncomps,k(i),tempx(i),rhox(i),ne)


    N = 1.0d0 / (mproton * abarx(i))

    invrhox = 1.0d0/rhox(i)
    ybar = ne / N * invrhox

    sum1 = 0.0d0
    sum2 = 0.0d0
    sum3 = 0.0d0
    e_ioncorr = 0.0d0

    do j=1, ncomps

      zion = comp_details(j,2)

      if(zion.lt.(1.0d0-1.0d-6)) cycle !don't consider neutrons

      sum_pot = 0.0d0
      sum_pot_frac = 0.0d0

      ymax = 0.0d0
      do l=1, int(zion) + 1
        if(ion_fractions(j,l,k(i)).gt.ymax) then
            ymax = ion_fractions(j,l,k(i))
            lmax = l
        end if
        if(l.gt.1) then
            sum_pot = sum_pot + xxip(int(zion),l-1)
            sum_pot_frac = sum_pot_frac + ion_fractions(j,l,k(i)) * sum_pot
        end if
      enddo

      !Here we look for the two most populated ionization states of
      !each element. y_r is the fraction of atoms in the highest of
      !these two states. chi_r is the ionization energy between these
      !two states
      if(lmax.eq.1) then
        y_r = ion_fractions(j,lmax+1,k(i))
        chi_r = xxip(int(zion),lmax)
      else if(lmax.eq.(int(zion) + 1)) then
        y_r = ion_fractions(j,lmax,k(i))
        chi_r = xxip(int(zion),lmax-1)
      else
        if(ion_fractions(j,lmax+1,k(i)).gt.ion_fractions(j,lmax-1,k(i))) then
            y_r = ion_fractions(j,lmax + 1,k(i))
            chi_r = xxip(int(zion),lmax)
        else
            y_r = ion_fractions(j,lmax,k(i))
            chi_r = xxip(int(zion),lmax-1)
        end if
      end if

      !number abundance of the j-th element
      nu_j = abarx(i) * comp(k(i),j) / comp_details(j,1)

      sum1 = sum1 + nu_j * y_r * (1 - y_r)

      sum2 = sum2 + nu_j * chi_r * y_r * (1 - y_r)

      sum3 = sum3 + nu_j * chi_r * chi_r * y_r * (1 - y_r)

      e_ioncorr = e_ioncorr + N * nu_j * sum_pot_frac

    enddo

    !Ion pressure
    pion = N * rhox(i) * kboltz * tempx(i)

    !Electron pressure
    pend = ybar * pion
    pednr = kdnr * (yex(i)*rhox(i))**(fivethirds)
    pedr = kdr * (yex(i)*rhox(i))**(fourthirds)
    ped = (1.0d0/(pednr*pednr)+1.0d0/(pedr*pedr))**(-0.5d0)
    pe = (pend*pend + ped*ped)**(0.5d0)

    !Radiation pressure
    pradx(i) = a_rad * tempx(i)**4 * overthree

    ffactor = (fivethirds)*ped*ped/(pednr*pednr) &
              +(fourthirds)*ped*ped/(pedr*pedr)

    !Total pressure
    px(i) = pion + pe + pradx(i)


    !Total specific internal energy
    ex(i) = 1.5d0 * pion * invrhox + 3.0d0 * pradx(i) * invrhox  &
              + pe / (rhox(i) * (ffactor - 1.0d0)) + e_ioncorr

    dpdtx(i) = pion / tempx(i) + 4.0d0 * pradx(i) / tempx(i)     &
          + pend*pend / (pe * tempx(i)) * (1.0d0                 &
              + ( 1.5d0*sum1 + sum2/(kboltz*tempx(i)) )/(ybar + sum1) )

    dedtx(i) = 1.5d0 * N * kboltz  +  4.0d0 * a_rad * tempx(i)**3 * invrhox    &
          + pend*pend / (pe * tempx(i) * rhox(i) * (ffactor - 1.0d0)) * (1.0d0 &
              + ( 1.5d0*sum1 + sum2/(kboltz*tempx(i)) )/(ybar + sum1) )        &
          + N * ( sum2 * ( 1.5d0*ybar - sum2/(kboltz*tempx(i)) )/(ybar + sum1) &
              + sum3/(kboltz*tempx(i)) ) / tempx(i)

    chi_T = tempx(i) * dpdtx(i) / px(i)

    chi_rho = rhox(i) / px(i) * ( N * kboltz * tempx(i)   &
       + (ybar * pend*pend /(ybar + sum1) + ffactor * ped*ped )/(pe * rhox(i)) )

    c_V = dedtx(i)

    Gamma1 = chi_T*chi_T * px(i) / (c_V * rhox(i) * tempx(i)) + chi_rho

    cs2x(i) = Gamma1 * px(i) * invrhox

    free_electron_frac(k(i)) = ne * mproton / (yex(i) * rhox(i))

  end do

end subroutine paczynski_eos
