!#############################################################################
!
! This routine is a part of the Timmes EOS with Saha ionization taken from
!
! http://cococubed.asu.edu/code_pages/eos_ionize.shtml
!
! It is used and distributed as part of SNEC by permission of Frank Timmes
!
! For more information on the Timmes EOS see:
!    timmes & arnett, apj supp. 125, 277, 1999
!    timmes & swesty, apj supp. 126, 501, 2000
!
!#############################################################################

      real*8 function stat_weight(zion,stage)
      include 'implno.dek'
      include 'const.dek'

!returns the statistical weight of the ionization state 

!this routine includes most of the weights from hydrogen to nickel.
!i am unsure of scandium, vanadium, manganese, and cobalt, so i've
!set all those to unity. 


!input is the element charge zion (a convenient index) and the ionization stage.


!declare the pass
      integer          stage
      real*8 zion
      
!locals
      integer          j,k,kmax,kmaxp1,ifirst,izion 
      parameter        (kmax = 30, kmaxp1 = kmax + 1)
      real*8 wgt(kmax,kmaxp1),off_table
      parameter        (off_table = 1.0e30)


!first time flag
      data    ifirst/0/
   
!hydrogen
      data (wgt(1,k), k=1,kmaxp1)/&
           2.0, 1.0,      &
           29*off_table/

!helium
      data (wgt(2,k), k=1,kmaxp1)/&
           1.0, 2.0, 1.0,      &
           28*off_table/

!lithium
      data (wgt(3,k), k=1,kmaxp1)/&
           2.0, 1.0, 2.0, 1.0,&
           27*off_table/

!berrylium
      data (wgt(4,k), k=1,kmaxp1)/&
           1.0, 2.0, 1.0, 2.0, 1.0,&
           26*off_table/

!boron
      data (wgt(5,k), k=1,kmaxp1)/&
           6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           25*off_table/

!carbon
      data (wgt(6,k), k=1,kmaxp1)/&
           9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           24*off_table/

!nitrogen
      data (wgt(7,k), k=1,kmaxp1)/&
           4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           23*off_table/

!oxygen
      data (wgt(8,k), k=1,kmaxp1)/&
           9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           22*off_table/

!fluorine
      data (wgt(9,k), k=1,kmaxp1)/&
           6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           21*off_table/


!neon
      data (wgt(10,k), k=1,kmaxp1)/&
           1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           20*off_table/


!sodium
      data (wgt(11,k), k=1,kmaxp1)/&
           2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           19*off_table/


!magnesium
      data (wgt(12,k), k=1,kmaxp1)/&
           1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, &
           1.0,&
           18*off_table/


!aluminum
      data (wgt(13,k), k=1,kmaxp1)/&
           6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, &
           2.0, 1.0,&
           17*off_table/


!silicon 
      data (wgt(14,k), k=1,kmaxp1)/&
           9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, &
           1.0, 2.0, 1.0,&
           16*off_table/


!phosphorous
      data (wgt(15,k), k=1,kmaxp1)/&
           4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, &
           2.0, 1.0, 2.0, 1.0,&
           15*off_table/


!sulfer
      data (wgt(16,k), k=1,kmaxp1)/&
           9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, &
           1.0, 2.0, 1.0, 2.0, 1.0,&
           14*off_table/


!chlorine
      data (wgt(17,k), k=1,kmaxp1)/&
           6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, &
           6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           13*off_table/


!argon
      data (wgt(18,k), k=1,kmaxp1)/&
           1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, 4.0, &
           9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           12*off_table/


!pottasium
      data (wgt(19,k), k=1,kmaxp1)/&
           2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, 9.0, &
           4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           11*off_table/


!calcium
      data (wgt(20,k), k=1,kmaxp1)/&
           1.0, 2.0, 1.0, 6.0, 9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 6.0, &
           9.0, 4.0, 9.0, 6.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           10*off_table/


!scandium
      data (wgt(21,k), k=1,kmaxp1)/&
           10.0, 15.0, 10.0, 19*1.0,&
           9*off_table/


!titanium
      data (wgt(22,k), k=1,kmaxp1)/&
           21.0, 28.0, 21.0, 10.0, 1.0, 6.0, 9.0, 4.0, 9.0, 2.0, 1.0,&
           2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, &
           8*off_table/


!vanadium
      data (wgt(23,k), k=1,kmaxp1)/&
           28.0, 25.0, 28.0, 21*1.0,&
           7*off_table/


!chromium
      data (wgt(24,k), k=1,kmaxp1)/&
           7.0, 6.0, 25.0, 28.0, 21.0, 10.0, 1.0, 6.0, 9.0, 2.0, 1.0,&
           2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, &
           2.0, 1.0,     & 
           6*off_table/


!manganese
      data (wgt(25,k), k=1,kmaxp1)/&
           6.0, 7.0, 6.0, 23*1.0,&
           5*off_table/


!iron
      data (wgt(26,k), k=1,kmaxp1)/&
           25.0, 30.0, 25.0, 6.0, 25.0, 28.0, 21.0, 10.0, 1.0,&
           6.0, 9.0, 2.0, 9.0, 2.0, 1.0,        &
           2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, &
           4*off_table/


!cobalt
      data (wgt(27,k), k=1,kmaxp1)/&
           28*1.0,&
           3*off_table/


!nickel
      data (wgt(28,k), k=1,kmaxp1)/&
           21.0, 10.0, 21.0, 10.0, 10.0, 10.0, 25.0, 28.0, 6.0, 6.0, &
           6.0, 6.0, 9.0, 6.0, 9.0, 6.0, 1.0, &
           2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0,&
           2*off_table/       



!copper
      data (wgt(29,k), k=1,kmaxp1)/&
           30*1.0,&
           1*off_table/            


!zinc
      data (wgt(30,k), k=1,kmaxp1)/&
           31*1.0/      



!this makes the output simple enough
      izion = int(zion)

      if (izion .gt. 30 .or. stage .gt. izion+1) &
           stop 'bad input in stat_weight'

      stat_weight = wgt(izion,stage)

      return
      end

