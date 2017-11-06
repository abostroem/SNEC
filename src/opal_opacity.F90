      subroutine instruct
!
! The vast majority of the code in this file is from the OPAL opacity
! project: http://opalopacity.llnl.gov/
!
!---------------------- WARNING from the SNEC team: ---------------------------
!In the version of the OPAL interpolator, adopted
!for the need of SNEC, the set of routines for the initial smoothing of the
!data is absent and the parameter ismdata in readco is set to 1.
!------------------------------------------------------------------------------
!
!-----VERSION of November 20, 1995-----------------------------------------
!-----------------------------------------------------------------------
!
!     This subroutine contains instructions for using the subroutine
!OPAC( z, xh, xxc, xxo, t6, r ).
!     The purpose of the subroutine OPAC is to perform 4 or 5 variable
!interpolation on  log10(kappa).  The opacity tables to be interpolated
!are known to have somewhat random numerical errors of a few percent.
!Consequently, adjusting the data prior to performing the interpolation
!is justified at this level.  The code first reads the original(unsmoothed)
!tabular data, this data is then passed through a smoothing filter; using
!a set of routines developed by Mike Seaton (see M.J. Seaton,MNRAS 265,
!L25,1993). It is the adjusted data that is actually used in carrying out the
!interpolations in OPAC.  Steps are taken, as described below to insure
!that the interpolated data is also smooth. The initial adjustment step
!helps improve the smoothness of the OPAC output,particularly at the
!smallest values of R. The medium to large R output is only slightly
!effected by this step. It takes only a few seconds to carryout the initial
!data smoothing step.  This step can be skipped by setting the parameter
!ismdata=1 in subroutine readco.
!
!     The interpolation variables are :
!
!     xh     The hydrogen mass fraction, X
!     xxc    The enhanced carbon mass fraction, delta Xc.  
!            The total C mass fraction, Xc, is the sum of the initial 
!            amount included in the metal mass fraction, Z, and 
!            delta Xc
!     xxo    The enhanced oxygen mass fraction, delta Xo
!     t6     The temperature in millions of degrees Kelvin, T6
!     r      =rho(gm/cc)/T6**3, R
!
!Additional input to OPAC is:
!
!      z      The metallicity, Z
!
!     An interpolation between overlapping quadratics is used to obtain
!smoothed results.  A 4x4 grid in logT6 and logR is used to interpolate
!in four different 3x3 sub-grids. Linear interpolation between quadratic
!fits in these different  sub-grids gives smoothed results in both log T6
!and Log R directions. This procedure produces results that are similar
!to bicubic spline interpolation, but require storage of only local 
!information.  
!    Each of the individual tables in a file Gx**x**z covers 70 temperatures
!in the range logT=3.75[T6=0.0056341325]( referred to as temperature 1) to
!logT=8.7[T6=501.187] and 19 values of log R in the range -8 (referred to as 1)
!to +1; at half-integer steps.  (NOTE: earlier tables were explicitly in
!terms of T6. For convenience the present tables tabulate log Kappa vs logT. The
!interpolation however still uses T6 for the temperature variable)
!For specialized problems, if storage space is a problem, a reduced set of
!data can be input .  This requires a recompilation with altered parameter
!values.  In order to limit the range of R, set the parameter nrb= index of
!lowest value of log R to use(count from log R=-8).  Then set the parameter
!nre to the index of the largest value of log R to use.  (NOTE: nre-nrb must
!be at least 5). To ignore the first ntb-1 values of T6 (starting from
!temperature 1) set the parameter ntb to desired value.  No other parameters
!should be modified.
!
!----------------------- WARNING from the SNEC team: --------------------------
!The version of the OPAL interpolator, adopted for
!the need of SNEC, does not allow changing the parameter ntb, it should stay
!equal to 1 everywhere in the routine.
!------------------------------------------------------------------------------
!
!     A five variable interpolation is done when X is greater than 0.  In this
!case it is assumed that variable values of X are also needed.  We have provided
!sets of tables for X=0, 0.03, 0.1, 0.35 and 0.7, for each of the metallicities
!0.0,0.001, 0.004,0.01,0.02,0.03,05 and .1.  The five sets of tables associated with
!a particular value of Z (eg. Gx0z01, Gx03x01,Gx10z01, Gx35z01,
!Gx70z01) should be placed in files named codataa, codatab, codatac, codatad,
!and codatae, respectively.  To create storage for this data the routines must
!be recompiled with the parameter mx=5.  Again if storage is a problem the T6
!and log R ranges can be restricted.  This version of the interpolation routine
!does not interpolate in Z.  Values of Z not in the table can be obtained by
!interpolating the existing tables to produce similar tables for the Z of interest.
!Interpolation in xh,xxo,xxc,t6, and r can be obtained as just described
!    A 4 variable interpolation in xxc, xxo, t6 and r is performed in the special
!case  X=0.  The set of 60 data tables in (xxc,xxo) for a given Z that have
!been provided for X=0 should be placed in a file called 'codataa'.  This file
!will be read from unit 2 in the subroutine readco.  In this special case the
!set of routines provided should be compiled with the parameter mx=1 (there are
!4 occurances).  (NOTE: The version of the code described above, intended
!for X> 0, also  handles this case, but takes more storage space since mx=5).
!     If you want to work with a single value of X which is not zero, i.e.,
!X=0.03, 0.1, 0.35, or 0.70, then compile the code with mx=1 but change the statement
!
!     data (xa(i), i=1,5 )/0.0, 0.03, 0.1, 0.35, 0.7/
!
!If for example you want to use only the table for X=.70, then
!
!     data (xa(i), i=1,5 )/0.70, 0.03, 0.1, 0.35, 0.0/
!
!You also need to place the tables for X=0.7 into a file named
!codataa.
!***CAUTION***
!    As a result of the mixing procedure used to calculate the data a few
!X=0.0, low T-small R, table values fell outside the range of T and R accessible
!from the X=0.35 data directly calculated for this purpose.  These T-R 
!locations are filled in with 9.99 (or for diagnostic purposes in some cases
!larger values.  At the same locations the derivatives are set to 99.9.  When
!T-R falls in this region a message is issued by the interpolation code to 
!inform the user of this situation.  Presumable very few users will have
!applications that take them into this region.
!
!     Your routine that performs the call to OPAC should include the
! statement:
!
!    common/e/ opact,dopact,dopacr,dopactd
!
!    These variables have the following meanings:
!   
!    OPACT        Is the Log of the Rosseland mean opacity: Log(kappa)  
!    DOPACT      Is Dlog(kappa)/Dlog(T6)   !at constant R
!    DOPACTD     Is Dlog(kappa)/Dlog(T6)   !at constant density
!    DOPACR      Is Dlog(kappa)/Dlog(R),   !at constant T6
!
      dum=0.0
      return
      end
!
!***********************************************************************
      subroutine opac(z,xh,xxc,xxo,t6,r)
!..... The purpose of this subroutine is to interpolate log kappa
!          and obtain smooth derivatives.
!      in C/O abundance and in T6, R,i.e. (xc,xo,T6,R)
!      xxc=Xc=carbon mass fraction
!      xxo=Xo=oxygen mass abundance
!      t6=T6=temperature in millions of degrees kelvin
!      r=R=density(g/cm**3)/T6**3
!..... to use opac insert common/e/ in the calling routine.
!      This common contains interpolated values for kappa and its
!      first derivities.
!
      save
      integer w
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb &
          ,ntabs=60,ntm=134,ntb=1,nt=ntm+1-ntb)
      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc), &
          opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo) &
          ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx) &
          ,nc,no
      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr), index(101), &
          t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr) &
          ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) &
          ,t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs), &
          zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
      common/d/dkap
      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
!..... OPACT- opacity obtained from a quadraric interpolation at
!      fixed log T6 at three values of log R; followed by quadratic
!      interpolation along log T6. Results smoothed by mixing
!      overlapping quadratics.
!..... DOPACT- is Dlog(k)/Dlog(T6) smoothed by mixing quadratics.
!..... DOPACR- is  Dlog(k)/Dlog(R) smoothed by mixing quadratics.
      common/e/ opact,dopact,dopacr,dopactd
      common/recoin/ itimeco,mxzero
!
!..... INDEX refers to the index,i, of xc(i) or xo(i); xc(i) and xo(i)
!      are the abundance grid points.
      iop=1   ! provides smoothed interpolations; iop=0 gives no smoothing
      if(nr .lt. 6) go to 65
      if((xh .gt. 1.e-6) .and. (mx .lt.4)) go to 64
!
!..... set-up C/O axis points
      xxco=xxc+xxo
      if(z+xh+xxco-1.e-6 .gt. 1 ) go to 61
      zzz=z+0.001
      xxh=xh
       xxci=xxc
      xxoi=xxo
      xxi=xh
      t6i=t6
      ri=r
      
!
!..... convert xxc and xxo to logarithmic shifted by Z
      cxx=log10(zzz+xxc)
      oxx=log10(zzz+xxo)
      xxx=log10(0.005+xh)
      slt=log10(t6)
      slr=log10(r)
 
!..... set X indices
        ilo=2
        ihi=mx
    8   if(ihi-ilo .gt. 1) then
          imd=(ihi+ilo)/2
            if(xh .le. xa(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 8
        endif
        i=ihi
        mf=i-2
        mg=i-1
        mh=i
        mi=i+1
        mf2=mi
        istep1=1
        if (mx .gt. 1) then
        istep1=mx-1
        if((xh .le. xa(2)+1.e-7) .or. (xh .ge. xa(istep1)-1.e-7)) mf2=mh
        endif
!
        if ((mx .eq. 1) .or. (xh .lt. 1.e-6)) then
          mf=1
          mg=1
          mh=1
          mi=2
          mf2=1
        endif

      if (itime(1) .ne. 12345678) then 
          alr(1)=-8.+(nrb-1)*0.5 
          do i=2,nr  
            alr(i)=alr(i-1)+0.5  
          enddo

        alt(1) = -3.3
        do i=2, 5
            alt(i) = alt(i-1)+0.05
        end do
        alt(6) = alt(5)+0.01
        do i=7, 11
            alt(i) = alt(i-1)+0.02
        end do
        do i=12, 60
            alt(i) = alt(i-1)+0.01
        end do
        do i=61, 110
            alt(i) = alt(i-1)+0.05
        end do
        do i=111, 131
            alt(i) = alt(i-1)+0.1
        end do
        do i=132,134
            alt(i) = alt(i-1)+0.2
        end do

          do i=1,nt
          t6list(i)=10.**alt(i)
          enddo
      endif   
      
        ilo=2
        ihi=nr
   12     if(ihi-ilo .gt. 1) then
          imd=(ihi+ilo)/2
            if(slr .le. alr(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 12
          endif
        i=ihi
        l1=i-2
        l2=i-1
        l3=i
        l4=l3+1
        
!
        ilo=2
        ihi=nt
   11     if(ihi-ilo .gt. 1) then
          imd=(ihi+ilo)/2
            if(t6 .le. t6list(imd)+1.e-7) then
              ihi=imd
            else
              ilo=imd
            endif
          go to 11
          endif
        i=ihi
        k1=i-2
        k2=i-1
        k3=i
        k4=k3+1
        k3s=k3+ntb-1
        l3s=l3+nrb-1
        

!-----set-up indices to skip when low T&R data missing for X=0.
      kmin=0
      k1in=k1
      iadvance=0
      mfin=mf
      if ((mfin .eq. 1) .and. (co(1,1,1,k1,l1) .gt. 9.)) then! data missing
      
          do i=1,6
            if (co(1,1,1,i,l1) .gt. 9.)  then
              if (xh .lt. .1) then
               kmin=i+1
              else

                if (iadvance .eq. 0) then  ! sfift X index to avoid X=0.
                iadvance=iadvance+1
                mf=mf+1 
                mg=mg+1
                mh=mh+1
                mi=mi+1
                mf2=mf2+1
                endif
              endif
            endif
          enddo
          
          if ((iadvance .eq. 0) .and. (k1 .le. kmin) .and. &
             (slt .le. alt(kmin))) then
          k1=kmin
              if ((co(1,1,1,kmin,l1+1) .lt. 9.) .and. &
                ((slr+.01) .gt. alr(l1+1))) then
              l1=l1+1
              kmin=0
              k1=k1in
                  do i=1,6
                  if (co(1,1,1,i,l1) .gt. 9.) kmin=i+1
                  enddo
              if ((kmin .ne. 0) .and. (k1in .lt. kmin)) k1=kmin
              endif
          endif
      
          if ((slt+.001) .lt. alt(k1)) then
          write (*,'("OPAL data not available for X=", f7.5," logT6=", f7.3, &
                    " logR=",f7.3)') xh,slt,slr
          opact=30.
          dopact=99.
          dopacr=99.
          dopactd=99.
          return  
          endif
          
      l2=l1+1
      l3=l2+1
      l4=l3+1
      l3s=l3+nrb-1
      k2=k1+1
      k3=k2+1
      k4=k3+1
      k3s=k3+ntb-1
      
      endif
!-----end of check for missing data

      do 123 m=mf,mf2
       if(mx .ge. 4) then
!.....  C and O  fractions determined by the ray through the origin that
!       also passes through the point (Xc,Xo). Specific interpolation 
!       values determined by tabulated X values;i.e. xa(m).  Inter-
!       polation along the ray gives log (kappa(Xc,Xo)).  (Advantage
!       of method: keeps indices within table boundaries)
!      Subtract Z to prevent out-of-range C+O values for small X
           if(1.-xh-z.gt.1.e-6)then
              cmod=(1.-xa(m)-z)/(1.-xh-z)
           else
              cmod=0.
           endif
       xxc=cmod*xxci
       xxo=cmod*xxoi
       cxx=log10(zzz+xxc)
       oxx=log10(zzz+xxo)
       endif
!      ythism=z+xa(m)+xxc+xxo

         do i=1,mc
         xhe=1.-xa(m)-z
         nc=i
         no=i
         xc(i)=xcs(i)
         xo(i)=xos(i)
!          if(xcs(i) .ge. xhe-1.e-6) then
           if(xcs(i) .gt. xhe) then
           xc(i)=xhe
           xo(i)=xhe
           go to 3
           endif
         enddo
    3    continue
    
!
      if(itime(m) .ne. 12345678) then
      itime(m)=12345678
      mxzero=0
        do i=1,mx
        xx(i)=log10(0.005+xa(i))
        if(xa(i) .eq. 0.0) mxzero=i
        enddo
!..... this is the first time throught this m. Calculate the decadic
!      log of the perimeter points shifted by Z+0.001(to avoid divergence 
!      at origin); m refers to xa(m); the hydrogen table value.
!
!      note that the nc-th elements are sometimes needed!
        do i=1,nc
          oxf(m,i)=log10(zzz+xo(i))
          cxf(m,i)=log10(zzz+xc(i))
          xcdf(m,i)=-xo(i)+xo(no)
          xodf(m,i)=-xc(i)+xc(nc)
          cxdf(m,i)=log10(zzz+xcdf(m,i))
          oxdf(m,i)=log10(zzz+xodf(m,i))
        enddo
!
!      note that the nc-th elements are sometimes needed!
        do i=1,nc
          ox(i)=oxf(m,i)
          cx(i)=cxf(m,i)
          xcd(i)=xcdf(m,i)
          xod(i)=xodf(m,i)
          cxd(i)=cxdf(m,i)
          oxd(i)=oxdf(m,i)
        enddo
!
!.....read the data files
        call readco
      
      endif
!
        do i=1,nc
          ox(i)=oxf(m,i)
          cx(i)=cxf(m,i)
          xcd(i)=xcdf(m,i)
          xod(i)=xodf(m,i)
          cxd(i)=cxdf(m,i)
          oxd(i)=oxdf(m,i)
        enddo
!
!..... Determine log R and log T6 grid points to use in the
!      interpolation.
      if((slt .lt. alt(1)).or.(slt .gt. alt(nt))) go to 62
      if((slr .lt. alr (1)).or.(slr .gt. alr(nr))) go to 62
!
      if (m .eq. mf) then  !  calculate table indices
!
      if((mf2 .ne. mxzero) .and. (k3s .gt. ntm)) go to 62
        do i=14,18
          if((l3s .gt. i) .and. (k3s .gt. nta(i+1))) go to 62
        enddo
      ip=3
      iq=3
      ntlimit=nta(l3s)
      if((k3s .eq. ntlimit) .or. (iop .eq. 0)) then
        ip=2
        iq=2
      endif
      if(t6 .le. t6list(2)+1.e-7 .or. iop .eq. 0) ip=2

      if((l3 .eq. nr) .or. (iop .eq. 0)) then ! right edge of full table
        iq=2
        ip=2
      endif
      if(slr .le. alr(2)+1.e-7 .or. iop .eq. 0) iq=2
      endif
!
      xodp=max(-xxc+xc(nc),0.)
      xcdp=max(-xxo+xo(no),0.)
      is=0
!
      call cointerp(xxc,xxo)
  123 continue

        if((zz(mg,1) .ne. zz(mf,1)) .or. (zz(mh,1) .ne. zz(mf,1))) then
          write(*,'("Z does not match Z in codata* files you are"  &
         ," using")')
          stop
        endif
      if(z .ne. zz(mf,1)) go to 66
      xxc=xxci   ! restores input value; necessary if stop replaced 
!                  with return
      xxo=xxoi   ! restores input value
      is=0
      iw=1
      do 45 ir=l1,l1+iq
       do 46 it=k1,k1+ip
         if((mx .eq. 1) .or. (mf2 .eq. 1)) then
           opk(it,ir)=opl(mf,it,ir)
           go to 46
         endif
         opk(it,ir)=quad(is,iw,xxx,opl(mf,it,ir),opl(mg,it,ir) &
            ,opl(mh,it,ir),xx(mf),xx(mg),xx(mh))
         is=1
   46 continue
   45 continue

      if (mi .eq. mf2) then  ! interpolate between quadratics
      is=0
      iw=1
       dixr=(xx(mh)-xxx)*dfsx(mh)
      do 47 ir=l1,l1+iq
        do it=k1,k1+ip
        opk2(it,ir)=quad(is,iw,xxx,opl(mg,it,ir),opl(mh,it,ir) &
           ,opl(mi,it,ir),xx(mg),xx(mh),xx(mi))
        opk(it,ir)=opk(it,ir)*dixr+opk2(it,ir)*(1.-dixr)
        is=1
        enddo
   47 continue
!     interpolate X between two overlapping quadratics
      endif
      is=0
!
!..... completed H,C,O interpolation. Now interpolate T6 and log R on a
!      4x4 grid. (log(T6(i)),i=i1,i1+3),log(R(j)),j=j1,j1+3)).Procedure
!      mixes overlapping quadratics to obtain smoothed derivatives.
!
        
      call t6rinterp(slr,slt)
      return
!
   61 write(*,'(" Mass fractions exceed unity")')
      xxc=xxci   ! restores input value; required if stop replaced 
!                  with a return
      xxo=xxoi   ! restores imput value
      stop
   62 write(*,'(" T6/LogR outside of table range")')
      xxc=xxci   ! restores input value; required if stop replaced
!                  with a return 
      xxo=xxoi   ! restores input value
      stop
   64 write(*,'(" X not equal to zero: To run this case it &
          is necessary"/ "to recompile with parameter (mx=5)")')
      stop
   65 write(*,'("Too few R values; NRE+1-NRB < 6")')
      stop
   66 write(*,'(" Z does not match Z in codata* files you are", &
          " using")')
      stop
      end

!************************************************************************
      subroutine cointerp(xxc,xxo)
!     The purpose of this subroutine is to interpolate in C and O abund-
!     ances.
      save
      integer w
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb &
          ,ntabs=60,ntm=134,ntb=1,nt=ntm+1-ntb)
      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc), &
          opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo) &
          ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx) &
          ,nc,no
      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr), index(101), &
          t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr) &
          ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) &
          ,t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
       is=0
      if(xxco .lt. 1.e-6) then
        do ir=l1,l1+iq
          do it=k1,k1+ip
            opl(m,it,ir)=co(m,1,1,it,ir)
          enddo
        enddo
            is=1
            go to 123
      endif
!     include boundaries that could later cause division by 0!
      if(xxc .gt. xcd(3)-1.e-6) then
!__________
      oxdp=log10(zzz+xodp)
!     handle possibility that xodp=0
      fac=max(min((oxx-ox(1))/max(oxdp-ox(1),1.e-6),1.),0.)
      do 40 ir=l1,l1+iq
      do 41 it=k1,k1+ip
!
!                    interpolation in region c1
!
!     include boundaries that could later cause division by 0!
      if(xxc .gt. xcd(2)-1.e-6) then
      iw=1
      a(1,m)=quad(is,iw,cxx,co(m,nc-2,1,it,ir),co(m,nc-1,1,it,ir), &
          diag(m,1,it,ir),cx(nc-2),cx(nc-1),cx(nc))
      iw=iw+1
      a(2,m)=quad(is,iw,cxx,diag(m,1,it,ir),diag(m,2,it,ir), &
          diag(m,3,it,ir),cxd(1),cxd(2),cxd(3))
         do w=1,2
           b(w)=a(w,m)
         enddo
!     handle possibility that xodp=0
           opl(m,it,ir)=b(1)+(b(2)-b(1))*fac
      is=1
      go to 41
      endif
!                    interpolation in region c2
!
      iw=1
      a(1,m)=quad(is,iw,cxx,co(m,nc-2,1,it,ir),co(m,nc-1,1,it,ir), &
          diag(m,1,it,ir),cx(nc-2),cx(nc-1),cx(nc))
      iw=iw+1
      a(2,m)=quad(is,iw,cxx,co(m,n(m,2)-2,2,it,ir),co(m,n(m,2)-1,2,it, &
          ir),diag(m,2,it,ir),cx(n(m,2)-2),cx(n(m,2)-1),cxd(2))
      iw=iw+1
      a(3,m)=quad(is,iw,cxx,diag(m,1,it,ir),diag(m,2,it,ir) &
         ,diag(m,3,it,ir),cxd(1),cxd(2),cxd(3))
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(1),ox(2),oxdp)
      is=1
   41 continue
   40 continue
      if(is .eq. 1) go to 123
!__________
      endif
!
!                    interpolation in region c3 to c6
      is=0
!
      if(nc .ge. 5) then
!__________
      do 44 i=4,nc-1
!     do not go beyond middle (where c3-c6 overlaps o3-o6), and
        if((xxc .gt. xcd(i)-1.e-6) .and. (xxo .gt. xo(i-1)-1.e-6) .and. &
             (xcd(i-1) .gt. xc(i-1))) then
      do 42 ir=l1,l1+iq
      do 43 it=k1,k1+ip
        oxdp=log10(zzz+xodp)
        iw=1
        m1=i-1
        m2=i-2
        a(1,m)=quad(is,iw,cxx,co(m,n(m,m2)-2,m2,it,ir),co(m,n(m,m2)-1,m2 &
           ,it,ir),diag(m,m2,it,ir),cx(n(m,m2)-2),cx(n(m,m2)-1),cxd(m2))
        iw=iw+1
        a(2,m)=quad(is,iw,cxx,co(m,n(m,m1)-2,m1,it,ir),co(m,n(m,m1)-1,m1 &
           ,it,ir),diag(m,m1,it,ir),cx(n(m,m1)-2),cx(n(m,m1)-1),cxd(m1))
        iw=iw+1
        a(3,m)=quad(is,iw,cxx,diag(m,m2,it,ir),diag(m,m1,it,ir), &
           diag(m,i,it,ir),cxd(m2),cxd(m1),cxd(i))
         do w=1,3
           b(w)=a(w,m)
         enddo
        iw=iw+1
        opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(i-2),ox(i-1),oxdp)
        is=1
   43 continue
   42 continue
      if (is .eq. 1) go to 123
        endif
   44 continue
!__________
      endif
!
      if (is .eq. 1) go to 123
!
!     include boundaries that could later cause division by 0!
      if(xxo .gt. xod(3)-1.e-6) then
!__________
      cxdp=log10(zzz+xcdp)
!     handle possibility that xcdp=0
      fac=max(min((cxx-cx(1))/max(cxdp-cx(1),1.e-6),1.),0.)
      do 140 ir=l1,l1+iq
      do 141 it=k1,k1+ip
!
!                    interpolation in region  o1
!
!     include boundaries that could later cause division by 0!
      if(xxo .gt. xod(2)-1.e-6) then
      iw=1
      a(1,m)=quad(is,iw,oxx,co(m,1,no-2,it,ir),co(m,1,no-1,it,ir), &
          diago(m,no-1,it,ir),ox(no-2),ox(no-1),ox(no))
      iw=iw+1
      a(2,m)=quad(is,iw,oxx,diago(m,no-1,it,ir),diago(m,no-2,it,ir), &
         diago(m,no-3,it,ir),oxd(1),oxd(2),oxd(3))
        do w=1,2
          b(w)=a(w,m)
        enddo
!     handle possibility that xcdp=0
      opl(m,it,ir)=b(1)+(b(2)-b(1))*fac
      is=1
      go to 141
      endif
!                    interpolation in region  o2
!
      iw=1
      a(1,m)=quad(is,iw,oxx,co(m,1,no-2,it,ir),co(m,1,no-1,it,ir), &
          diago(m,no-1,it,ir),ox(no-2),ox(no-1),ox(no))
      iw=iw+1
      a(2,m)=quad(is,iw,oxx,co(m,2,n(m,2)-2,it,ir),co(m,2,n(m,2)-1,it, &
          ir),diago(m,no-2,it,ir),ox(n(m,2)-2),ox(n(m,2)-1),oxd(2))
      iw=iw+1
      a(3,m)=quad(is,iw,oxx,diago(m,no-1,it,ir),diago(m,no-2,it,ir), &
         diago(m,nc-3,it,ir),oxd(1),oxd(2),oxd(3))
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(1),cx(2),cxdp)
      is=1
  141 continue
  140 continue
      if(is .eq. 1) go to 123
!__________
      endif
!
!                    interpolation in region  o3 to o6
      is=0
      if(no .ge. 5) then
!__________
      do 144 i=4,no-1
!     do not go beyond middle (where o3-o6 overlaps c3-c6), and
      if((xxo .gt. xod(i)-1.e-6) .and. (xxc .gt. xc(i-1)-1.e-6) .and. &
          (xod(i-1) .gt. xo(i-1)-1.e-6)) then
      do 142 ir=l1,l1+iq
      do 143 it=k1,k1+ip
      cxdp=log10(zzz+xcdp)
      iw=1
      m2=i-2
      m1=i-1
      a(1,m)=quad(is,iw,oxx,co(m,m2,n(m,m2)-2,it,ir),co(m,m2,n(m,m2)- &
         1,it,ir),diago(m,no-m2,it,ir),ox(n(m,m2)-2),ox(n(m,m2)-1),oxd(m2))
      iw=iw+1
      a(2,m)=quad(is,iw,oxx,co(m,m1,n(m,m1)-2,it,ir),co(m,m1,n(m,m1)-1, &
          it,ir),diago(m,no-m1,it,ir),ox(n(m,m1)-2),ox(n(m,m1)-1),oxd(m1))
      iw=iw+1
      a(3,m)=quad(is,iw,oxx,diago(m,no-m2,it,ir),diago(m,no-m1,it,ir), &
         diago(m,no-i,it,ir),oxd(m2),oxd(m1),oxd(i))
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(m2),cx(m1),cxdp)
      is=1
  143 continue
  142 continue
      if (is .eq. 1) go to 123
      endif
  144 continue
!__________
      endif
!
      if (is .eq. 1) go to 123
!
!.....find index of C grid.
   52 ie=100*xxc+1
      iei=index(ie)+1
!     must also allow index = nc, to avoid extrapolation
      if (iei .gt. nc) iei=nc
!
        if(iei .gt. 3) then
          i1=iei-2
          i2=iei-1
          i3=iei
        else
          i1=1
          i2=2
          i3=3
        endif
!
!.....find index of O grid
      ie=100*xxo+1
      iej=index(ie)+1
!     must also allow index = no, to avoid extrapolation
      if (iej .gt. no) iej=no
!
        if(iej .gt. 3) then
          j1=iej-2
          j2=iej-1
          j3=iej
        else
          j1=1
          j2=2
          j3=3
        endif
!
!     lower-O part of grid: interpolate C before O
      if(j3.lt.no .and. i3.le.n(m,j3) .and. &
            (xxc.lt.xcd(j3)+1.e-6 .or. xxc.ge.xxo))then
      do 20 ir=l1,l1+iq
      do 21 it=k1,k1+ip
      iw=0
        do jx=j1,j1+2
          iw=iw+1
!     if i3=n(m,jx), then must replace cx(i3) with cxd(jx)
          a(iw,m)=quad(is,iw,cxx,co(m,i1,jx,it,ir),co(m,i2,jx,it,ir), &
             co(m,i3,jx,it,ir),cx(i1),cx(i2),min(cx(i3),cxd(jx)))
        enddo
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,oxx,b(1),b(2),b(3),ox(j1),ox(j2),ox(j3))
      is=1
   21 continue
   20 continue
!     else: high-O part of grid: must interpolate O before C
      else
       do ir=l1,l1+iq
       do it=k1,k1+ip
        iw=0
        do ix=i1,i1+2
          iw=iw+1
          if(j3.lt.n(m,ix))then
            a(iw,m)=quad(is,iw,oxx,co(m,ix,j1,it,ir),co(m,ix,j2,it,ir), &
               co(m,ix,j3,it,ir),ox(j1),ox(j2),ox(j3))
          else
            a(iw,m)=quad(is,iw,oxx,co(m,ix,j1,it,ir),co(m,ix,j2,it,ir), &
               diago(m,no-ix,it,ir),ox(j1),ox(j2),oxd(ix))
          endif
        enddo
        do w=1,3
          b(w)=a(w,m)
        enddo
      iw=iw+1
      opl(m,it,ir)=quad(is,iw,cxx,b(1),b(2),b(3),cx(i1),cx(i2),cx(i3))
      is=1
       enddo
       enddo
      endif
  123 continue
      return
      end

!***********************************************************************
      subroutine t6rinterp(slr,slt)
!     The purpose of this subroutine is to interpolate in logT6 and logR
      save
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb &
          ,ntabs=60,ntm=134,ntb=1,nt=ntm+1-ntb)
      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc), &
          opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo) &
          ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx) &
          ,nc,no
      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr), index(101), &
          t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr) &
          ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) &
          ,t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/d/dkap
      common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq,xodp,xcdp,xxco,cxx,oxx
      common/e/ opact,dopact,dopacr,dopactd
!
      is=0
      iu=0
      do kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quad(is,iw,slr,opk(kx,l1),opk(kx,l2),opk(kx,l3), &
           alr(l1),alr(l2),alr(l3))
          if(iq.eq.3) then
            iw=2
            q(iu)=quad(is,iw,slr,opk(kx,l2),opk(kx,l3),opk(kx,l4), &
               alr(l2),alr(l3),alr(l4))
          endif
        is=1
      enddo
!
      is=0
      iw=1
!..... k and Dlog(k)/dlog(T6) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      opact=quad(is,iw,slt,h(1),h(2),h(3),alt(k1),alt(k2),alt(k3))
      dopact=dkap
      dkap1=dkap
        if (iq.eq.3) then
!.....    k and Dlog(k)/Dlog(T6) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          opactq=quad(is,iw,slt,q(1),q(2),q(3),alt(k1),alt(k2),alt(k3))
          dkapq1=dkap
        endif
        if(ip.eq.3) then
!.....    k and Dlog(k)/Dlog(T6) in lower-left 3x3.
          opact2=quad(is,iw,slt,h(2),h(3),h(4),alt(k2),alt(k3),alt(k4))
          dkap2=dkap
!.....    k and Dlog(k)/Dlog(T6) smoothed in left 3x4
          dix=(alt(k3)-slt)*dfs(k3)
          dopact=dkap1*dix+dkap2*(1.-dix)
          opact=opact*dix+opact2*(1.-dix)
        if(iq.eq.3) then
 
!.....    k and Dlog(k)/Dlog(T6) in upper-right 3x3.
          opactq2=quad(is,iw,slt,q(2),q(3),q(4),alt(k2),alt(k3),alt(k4))
          dkapq2=dkap
          dopactq=dkapq1*dix+dkapq2*(1.-dix)
          opactq=opactq*dix+opactq2*(1.-dix)
         endif
        endif
!
      iu=0
      do lx=l1,l1+iq
        iw=1
        iu=iu+1
        h(iu)=quad(is,iw,slt,opk(k1,lx),opk(k2,lx),opk(k3,lx), &
           alt(k1),alt(k2),alt(k3))
          if(ip .eq. 3) then
            iw=2
            q(iu)=quad(is,iw,slt,opk(k2,lx),opk(k3,lx),opk(k4,lx), &
               alt(k2),alt(k3),alt(k4))
          endif
        is=1
      enddo
!
      is=0
      iw=1
!..... k and Dlog(k)/Dlog(R) in lower-left 3x3
      opacr=quad(is,iw,slr,h(1),h(2),h(3),alr(l1),alr(l2),alr(l3))
      dopacr=dkap
        if(ip .eq. 3) then
          opacrq=quad(is,iw,slr,q(1),q(2),q(3),alr(l1),alr(l2),alr(l3))
!.....    k and Dlog(k)/Dlog(R) in upper-left 3x3.
          dopacrq=dkap
        endif
        if(iq .eq. 3) then
!.....    k and Dlog(k)/Dlog(R) in lower-right 3x3.
          opact2=quad(is,iw,slr,h(2),h(3),h(4),alr(l2),alr(l3),alr(l4))
          dix2=(alr(l3)-slr)*dfsr(l3)
          dopacr=dopacr*dix2+dkap*(1.-dix2)
!.....        k and Dlog(k)/Dlog(T6) smoothed in both log(T6) and log(R)
              dopact=dopact*dix2+dopactq*(1.-dix2)
              opact=opact*dix2+opactq*(1.-dix2)
         endif
        if(ip .eq. 3) then
         if(iq .eq. 3) then
!.....    k and Dlog(k)/Dlog(R) in upper-right 3x3.
          opacrq=quad(is,iw,slr,q(2),q(3),q(4),alr(l2),alr(l3),alr(l4))
!.....        Dlog(k)/Dlog(R) smoothed in both log(T6) and Log(R).
              dopacrq=dopacrq*dix2+dkap*(1.-dix2)
            endif
              dopacr=dopacr*dix+dopacrq*(1.-dix)
        endif
      dopactd=dopact-3.*dopacr
        if (opact .gt. 1.e+15) then
          write(*,'("Interpolation indices out of range", &
                   ";please report conditions.")') 
          write(*,*) slt, slr, ip, iq, 10.0**cxx-zzz, 10.0**oxx-zzz, xxh
          write(*,*) xc(:)
          write(*,*) xo(:)
          stop
        endif
      if (opact .gt. 9) then
      opact=30.
      dopact=99.
      dopactr=99.
      dopactd=99.
      endif
      return
      end

!************************************************************************
      subroutine readco
!..... The purpose of this subroutine is to read the data tables
      save
      parameter (ismdata=1)   ! modified
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb &
          ,ntabs=60,ntm=134,ntb=1,nt=ntm+1-ntb)
      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo) &
          ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx) &
          ,nc,no
      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr), index(101), &
          t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr) &
          ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) &
          ,t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs), &
          zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
      common/alink/ NTEMP,NSM,nrlow,nrhigh,RLE,t6arr(200),coff(200,nr)  
      COMMON/CST/NRL,RLS,nset,tmax  ! modified
      common/e/ opact,dopact,dopacr,dopacrd
      character*1 dumarra(250)
      common/recoin/ itimeco,mxzero
      
!
        if (itimeco .ne. 12345678) then
        do i=1,mx
          do j=1,mc 
            do k=1,mo
              do l=1,nt
                do mq=1,nr
                  co(i,j,k,l,mq)=1.e+35
                enddo
              enddo
            enddo
          enddo
        enddo
        do i=1,mx
          do j=1,mc 
            do l=1,nt
              do mq=1,nr
                diag(i,j,l,mq)=1.e+35
                diago(i,j,l,mq)=1.e+35
              enddo
            enddo
          enddo
        enddo
        itimeco=12345678
        endif
      do 20 j=1,nc-1
       do 21 i=1,nc
         if(xcd(j) .ge. xc(i)) then
           n(m,j)=i+1
           if(xcd(j) .lt. xc(i)+1.e-6) n(m,j)=i
         endif
   21  continue
   20 continue
      n(m,nc)=0
!
      close (2)
!..... read X=0.0 tables
      if(m .eq. 1) open(2, FILE='tables/codataa')
!..... read X=0.03 tables
      if(m .eq. 2) open(2, FILE='tables/codatab')
!..... read X=0.10 tables
      if(m .eq. 3) open(2, FILE='tables/codatac')
!..... read X=0.35 tables
      if(m .eq. 4) open(2, FILE='tables/codatad')
!.....read X=0.70 tables
      if(m .eq. 5) open(2, FILE='tables/codatae')
!
!      read header
      read(2,'(a)') (dumarra(i),i=1,240)
!
      int=0
      do 1 j=1,no-1
      do 2 i=1,n(m,j)
        int=int+1
!
        read(2,'(f10.5)') dum
        read (2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)') &
           itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),xoa(m,int)
        xca(m,int)=min(xca(m,int),1.-x(m,int)-zz(m,int)-xoa(m,int))

        read(2,'(f10.5)') dum,dum,dum
        read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
        read(2,'(f10.5)') dum
          do k=1,ntm
            read(2,'(f4.2,19f7.3)') altin,(cof(k,l), l=1,nrm) 
            
            do ll=1,nrm   ! modified
              coff(k,ll)=cof(k,ll)
            enddo
          enddo
          if (isett6 .ne. 1234567) then
          do k=1,ntm  
            t6arr(k)=t6list(k)
          enddo  
          endif  
          isett6=1234567


      ll=1
      do 110 kk=nrb,nre
      alr(ll)=alrf(kk)
        do k=1,nt
            co(m,i,j,k,ll)=coff(k+ntb-1,kk)
        enddo
  110 ll=ll+1
    2 continue
    1 continue
      if(x(m,1) .ne. xa(m)) then
      write(*,'(" X in the codata? file does not match xa(m)")')
      stop
      endif
!
      do i=1, nc-1
       do k=1,nt
        do l=1,nr
          diag(m,i,k,l)=co(m,n(m,i),i,k,l)
        enddo
       enddo
      enddo
!
      do 6 j=1,no-1
        int=int+1
        read(2,'(f10.5)') dum
        read (2,'(7x,i3,26x,f6.4,3x,f6.4,3x,f6.4,5x,f6.4,5x,f6.4)') &
           itab(m,int),x(m,int),y(m,int),zz(m,int),xca(m,int),xoa(m,int)
        read(2,'(f10.5)') dum,dum,dum
        read(2,'(4x,f6.1,18f7.1)') (alrf(kk),kk=1,nrm)
        read(2,'(f10.5)') dum
         do k=1,ntm
           read (2,'(f4.2,19f7.3)') dum, (cof(k,l), l=1,nrm)
!     set up to smooth final "diago" opacity tables
           do l=1,nrm
              coff(k,l)=cof(k,l)
           enddo
        enddo

        ll=1
        do kk=nrb,nre
            do k=1,nt
                diago(m,j,k,ll)=cof(k+ntb-1,kk)
            enddo
        ll=ll+1
        enddo

    6 continue

      do i=2,nt
        dfs(i)=1./(alt(i)-alt(i-1))
      enddo
      do i=2,nr
       dfsr(i)=1./(alr(i)-alr(i-1))
      enddo
      istep=-1
      if (mx .gt. 1 ) then
        istep=1
        do i=2,mx,istep
          dfsx(i)=1./(xx(i)-xx(i-1))
        enddo
      endif
      return
      end
      
!
!************************************************************************
      function quad(ic,i,x,y1,y2,y3,x1,x2,x3)
!..... this function performs a quadratic interpolation.
      save
      common/d/dkap
      common/coquad/ xx(3),yy(3),xx12(30),xx13(30),xx23(30),xx1sq(30) &
          ,xx1pxx2(30)
      xx(1)=x1
      xx(2)=x2
      xx(3)=x3
      yy(1)=y1
      yy(2)=y2
      yy(3)=y3
        if(ic .eq. 0) then
          xx12(i)=1./(xx(1)-xx(2))
          xx13(i)=1./(xx(1)-xx(3))
          xx23(i)=1./(xx(2)-xx(3))
          xx1sq(i)=xx(1)*xx(1)
          xx1pxx2(i)=xx(1)+xx(2)
        endif
      c3=(yy(1)-yy(2))*xx12(i)
      c3=c3-(yy(2)-yy(3))*xx23(i)
      c3=c3*xx13(i)
      c2=(yy(1)-yy(2))*xx12(i)-(xx1pxx2(i))*c3
      c1=yy(1)-xx(1)*c2-xx1sq(i)*c3
      dkap=c2+(x+x)*c3
      quad=c1+x*(c2+x*c3)
      return
      end
!
!************************************************************************
      block data
      parameter (mx=5,mc=8,mo=8,nrm=19,nrb=1,nre=19,nr=nre+1-nrb &
          ,ntabs=60,ntm=134,ntb=1,nt=ntm+1-ntb)
      common/aa/ q(4), h(4), xcd(mc),xod(mc), xc(mc),xo(mo) &
          ,xcs(mc),xos(mo), cxd(mc),oxd(mo),cx(mc),ox(mo),zzz,xxh,xx(mx) &
          ,nc,no
      common/a/ co(mx,mc,mo,nt,nr), diag(mx,mc,nt,nr), index(101), &
          t6list(nt),alr(nr),n(mx,mc),alt(nt),diago(mx,mo,nt,nr),opk(nt,nr) &
          ,dfs(nt),dfsr(nr),a(3,mx),b(3),m,mf,xa(8),alrf(nrm),cof(ntm,nrm) &
          ,t6listf(ntm),opk2(nt,nr),dfsx(mx)
      common/b/ itab(mx,ntabs),nta(nrm),x(mx,ntabs),y(mx,ntabs), &
          zz(mx,ntabs),xca(mx,ntabs),xoa(mx,ntabs)
      common/aaa/ oxf(mx,mc),cxf(mx,mc),xcdf(mx,mc),xodf(mx,mc), &
          opl(mx,nt,nr),itime(mx),cxdf(mx,mc),oxdf(mx,mc)
      common/recoin/ itimeco,mxzero
      data itime/mx*0/,itimeco/0/
      data ( index(i),i=1,101)/1,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4, &
          4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6, &
          6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7, &
          7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
      data (xcs(i),i= 1,mc)/ 0.0,0.01,0.03,0.1,0.2,0.4,0.6,1.0/
      data (xos(i),i= 1,mc)/0.0,0.01,0.03,0.1,0.2,0.4,0.6,1.0/
      data (xa(i),i=1,5)/0.0,0.03,0.1,0.35,0.7/
!      data (nta(i),i=1,nrm)/14*70,69,64,60,58,57/
      data (nta(i),i=1,nrm)/14*134,133,128,124,122,121/
      end
