      PROGRAM SR
!
!-----------------------------------------------------------------------
!
!     A program to tabulate the O2 Schumann-Runge cross sections as a
!     function of Temperature using my wavelength intervals based on: 
!     MINSCHWANER K, SALAWITCH RJ, MCELROY MB, ABSORPTION OF 
!     SOLAR-RADIATION BY O2 - IMPLICATIONS FOR O3 AND LIFETIMES OF N2O, 
!     CFCL3, AND CF2CL2, JOURNAL OF GEOPHYSICAL RESEARCH-ATMOSPHERES 
!     98 (D6): 10543-10561 JUN 20 1993 
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER :: NWAVE = 203 , NSR = 16000
!
! Local variables
!
      REAL , DIMENSION(NSR) :: a , b , c , dat , sigma , wn_sr
      REAL :: adev , ave , aveq , curt , d , sdev , skew , t , var
      CHARACTER(400) :: fn
      INTEGER :: j , jchan , jhi , jjt , jjw , jlo , jsum , jt , jw , n
      REAL :: REAL
      REAL , DIMENSION(NWAVE) :: wave_nm_c , wave_nm_hi , wave_nm_lo ,  &
     &                           wave_wn_c , wave_wn_hi , wave_wn_lo
!
!-----------------------------------------------------------------------
!
      jchan = 1
!
!-----------------------------------------------------------------------
!
!     Read in wavelength intervals in microns
      fn = 'intervals-cm.dat'
!     print *,trim(fn)
      OPEN (jchan,FILE=TRIM(fn),STATUS='old')
!
!     Skip header
      DO jw = 1 , 2
         READ (jchan,*)
      ENDDO
!
!     read data
      DO jw = 1 , NWAVE
         READ (jchan,*) wave_wn_c(jw) , wave_wn_lo(jw) , wave_wn_hi(jw)
      ENDDO
!
      CLOSE (jchan)
!
!-----------------------------------------------------------------------
!
!     Read in wavelength intervals in microns
      fn = 'intervals-nm.dat'
!     print *,trim(fn)
      OPEN (jchan,FILE=TRIM(fn),STATUS='old')
!
!     Skip header
      DO jw = 1 , 2
         READ (jchan,*)
      ENDDO
!
!     read data
      DO jw = 1 , NWAVE
         READ (jchan,*) wave_nm_c(jw) , wave_nm_lo(jw) , wave_nm_hi(jw)
      ENDDO
!
      CLOSE (jchan)
!
!-----------------------------------------------------------------------
!
!     Temperature loop
      jjt = 0
      DO jt = 130 , 500
!
!-----------------------------------------------------------------------
!
         jjt = jjt + 1
         t = REAL(jt)
         d = ((t-100.)/10.)**2.
!
!        Read in wavelength intervals in microns
         IF ( t<=190. ) THEN
            fn = 'fitcoef_cold.txt'
         ELSEIF ( t>190 .AND. t<=280. ) THEN
            fn = 'fitcoef_mid.txt'
         ELSEIF ( t>280 .AND. t<=500. ) THEN
            fn = 'fitcoef_hot.txt'
         ENDIF
!        print  *,trim(fn)
         OPEN (jchan,FILE=TRIM(fn),STATUS='old')
!
!        Skip header
         DO jw = 1 , 3
            READ (jchan,*)
         ENDDO
!
!        read data
         DO jw = 1 , NSR
            READ (jchan,*) wn_sr(jw) , a(jw) , b(jw) , c(jw)
!           print *,wn_sr(jw),a(jw),b(jw),c(jw)
         ENDDO
! 
         CLOSE (jchan)
!
!        print  *,t,d 
!
!-----------------------------------------------------------------------
!
!        Wavelength loop
         jjw = 0
         DO jw = 1 , NWAVE
!
!-----------------------------------------------------------------------
!
            CALL POS(wn_sr,NSR,wave_wn_hi(jw),jlo)
            CALL POS(wn_sr,NSR,wave_wn_lo(jw),jhi)
            IF ( jlo/=jhi ) THEN
               jjw = jjw + 1
               jlo = MAX(jlo,1)
               jhi = MIN(jhi,NSR)
!              print      *,'jlo,jhi:',jlo,jhi
!              print      *,'cm-1:',wave_wn_hi(jw),wave_wn_lo(jw)
!              print      *,'nm  :',wave_nm_lo(jw),wave_nm_hi(jw)
               j = 1
               DO jsum = jlo , jhi
                  dat(j) = 1E-20*(a(jsum)*d**2.+b(jsum)*d+c(jsum))
                  j = j + 1
!              print '(a,1x,1p,g15.6,g15.2)','wn_sr(jw),sigma(jw):',wn_sr(jsum),sigma(jsum
               ENDDO
               n = j - 1
               aveq = SUM(dat(1:n))/REAL(n)
!              CALL MOMENT(NSR,dat,n,ave,adev,sdev,var,skew,curt)
!              sigma(jsum) = ave
               sigma(jsum) = aveq
!              print '(a,1x,1p,g13.6,3g10.2)','wn_sr(jw),sigma(jw):',wave_nm_c(jw),sigma(jsum)
!              print '(a,1x,1p,g13.6,3g10.2)','wn_sr(jw),sigma(jw):',wave_nm_c(jw),sigma(jsum)
!               PRINT '(i3,1x,1p,6g13.6)' , jw , wave_nm_c(jw) ,         &
!     &               wave_nm_lo(jw) , wave_nm_hi(jw) , t , sigma(jsum)
               PRINT '(i3,1x,1p,2g13.6)' , jw , t , sigma(jsum)
            ENDIF
!
!-----------------------------------------------------------------------
!
         ENDDO
!        End of wavelength loop
!
!-----------------------------------------------------------------------
!
      ENDDO
!     End of temperature loop
!
!-----------------------------------------------------------------------
!
      PRINT * , '# wavelengths, # temperatures:' , jjw , jjt
!
!-----------------------------------------------------------------------
!
      END PROGRAM SR
!
!-----------------------------------------------------------------------
!
      SUBROUTINE POS(Xx,N,X,J)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: J,N
      REAL :: X
      REAL , DIMENSION(N) :: Xx
!
! Local variables
!
      INTEGER :: jl , jm , ju
!
!-----------------------------------------------------------------------
!
!     Argrument list.
!
!     Name      Type              Description.
!     XX        Array of real     Monotonic array of length N.
!     (Unchanged on exit).
!
!     N         Integer           Length of array XX.
!     (Unchanged on exit).
!
!     X         Real              Value whose position in XX is
!     required.
!     (Unchanged on exit).
!
!     J         Integer           Index of X in array XX.
!     (Contains answer on exit).
!
!-----------------------------------------------------------------------
!
!     Given an array XX of length N, and given a value X, POS returns a
!     value J sucha that X lies between XX(J) and XX(J+1). XX must be
!     monotonic, either increasing or decreasing. J=0 or J=N is returned
!     to indicate that X is out of range.
!
!     The table entry J is found by bisection.
!
!     Based on Numerical Recipes, The art of scientific computing,
!     section 3.4, by Press, Flannery, Teukolsky & Vetterling,
!     Cambridge University Press, 1987.
!
!     Modified by: David Lary
!     ----------
!
!     Date Started : 7/2/1990
!
!     Last modified: 27/9/1991
!
!-----------------------------------------------------------------------
!
!     Initialize upper & lower limits.
      jl = 0
      ju = N + 1
!
!----------------------------------------------------------------------
!
      DO WHILE ( .TRUE. )
!
!----------------------------------------------------------------------
!
         IF ( .NOT..TRUE. ) THEN
            RETURN
!
         ELSEIF ( ju-jl>1 ) THEN
!
            jm = (ju+jl)/2
            IF ( Xx(N)>Xx(1) .EQV. X>Xx(jm) ) THEN
               jl = jm
            ELSE
               ju = jm
            ENDIF
!
!           Repeat untill the test condition is satisfied.
            CYCLE
         ENDIF
!
!        Set the output.
         J = jl
         EXIT
!
!----------------------------------------------------------------------
!
      ENDDO
!
!----------------------------------------------------------------------
!
      END SUBROUTINE POS
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE MOMENT(NMAX,Dat,N,Ave,Adev,Sdev,Var,Skew,Curt)
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL :: Adev , Ave , Curt , Sdev , Skew , Var
      INTEGER :: N
      INTEGER :: NMAX
      REAL , DIMENSION(NMAX) :: Dat
!
! Local variables
!
      REAL :: ep , p , s
      INTEGER :: j
!
!-------------------------------------------------------------------------------
!
      IF ( N==1 ) THEN
!          write(6,*)'Moment n=1'
         Ave = Dat(1)
         Adev = 0.
         Sdev = 0.
         Var = 0.
         Skew = 0.
         Curt = 0.
      ELSE
         s = 0.
         DO j = 1 , N
            s = s + Dat(j)
         ENDDO
         Ave = s/N
         Adev = 0.
         Var = 0.
         Skew = 0.
         Curt = 0.
         ep = 0.
         DO j = 1 , N
            s = Dat(j) - Ave
            ep = ep + s
            Adev = Adev + ABS(s)
            p = s*s
            Var = Var + p
            p = p*s
            Skew = Skew + p
            p = p*s
            Curt = Curt + p
         ENDDO
         Adev = Adev/N
         Var = (Var-ep**2/N)/(N-1)
         Sdev = SQRT(Var)
         IF ( Var/=0. ) THEN
            Skew = Skew/(N*Sdev**3)
            Curt = Curt/(N*Var**2) - 3.
         ELSE
            Skew = 0.
            Curt = 0.
!           WRITE (6,*) 'Sdev:' , Sdev
!           WRITE (6,*) 'Var :' , Var
!           WRITE (6,*)
!           &                'no skew or kurtosis when zero variance in moment'
         ENDIF
      ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE MOMENT
