      INTEGER :: Jchan , Jday , Jmonth , Jyear,is
      REAL :: S
!
! Local variables
!
      REAL :: basedir , zyear
      INTEGER :: jrd_day , jrd_month , jrd_year
      CHARACTER(200) :: longfn
      CHARACTER(4) :: ssn_st
!
!-----------------------------------------------------------------------
!
!     Sunspot data can be obtained from
!     http://www.oma.be/KSB-ORB/SIDC/
!     now
!     http://sidc.oma.be/index.php3
!
!-----------------------------------------------------------------------
!
jyear=2003
jmonth=12
jday=7
      jchan=1
      longfn = 'I:\AutoChemV6\Data\CrossSections\dayssn.dat'
      OPEN (Jchan,FILE=TRIM(longfn),STATUS='old')
      DO WHILE ( .TRUE. )
         READ (Jchan,'(i4,i2,i2,2x,f8.3,a4)',END=100) jrd_year , jrd_month , jrd_day , zyear ,&
                              & ssn_st
         IF ( (Jyear==jrd_year) .AND. (Jmonth==jrd_month) .AND.         &
            & (Jday==jrd_day) ) then
          print *,':'//trim(ssn_st)//':'
							  read(ssn_st,*)s
		  print *,jrd_year , jrd_month , jrd_day , zyear ,&
                              & s
			EXIT
			endif
      END DO
 100  CONTINUE
      CLOSE (Jchan)
!
!-----------------------------------------------------------------------
!     
      end
