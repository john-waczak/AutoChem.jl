      PROGRAM SR
!
!-----------------------------------------------------------------------
!
!     A program to read a tabulation of the O2 Schumann-Runge cross 
!     sections as a function of Temperature using my wavelength intervals 
!     based on: MINSCHWANER K, SALAWITCH RJ, MCELROY MB, ABSORPTION OF 
!     SOLAR-RADIATION BY O2 - IMPLICATIONS FOR O3 AND LIFETIMES OF N2O, 
!     CFCL3, AND CF2CL2, JOURNAL OF GEOPHYSICAL RESEARCH-ATMOSPHERES 
!     98 (D6): 10543-10561 JUN 20 1993 
!
!-----------------------------------------------------------------------
!
      USE DIRECTORIES
      USE KINETIC
      USE KINETIC_PARAM
      USE PHYSICAL
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER :: NWAVESR = 16 , NTSR = 371
      INTEGER , PARAMETER :: JPWAVE = 203
!
! Local variables
!
      REAL , DIMENSION(JPWAVE) :: ao2sr
      REAL , DIMENSION(NWAVESR,NTSR) :: sigma_sr
      INTEGER , DIMENSION(NWAVESR) :: iwave_sr
      REAL , DIMENSION(NTSR) :: t_sr , part
      CHARACTER(400) :: fn
      INTEGER :: j , jchan , jhi , jjt , jjw , jlo , jsum , jt , jw , n , ifl
      real :: twant,Value
      COMMON /KEEP_SR/ sigma_sr,iwave_sr,t_sr
!
      DATA ifl/0/
!
!-----------------------------------------------------------------------
!
      IF ( ifl==0 ) THEN
         ifl = 1
         jchan = 99
         PRINT * , ' ACSSR Minschwaner ( July 2005 )'
         fn = TRIM(directory%basedir)                                   &
             &//'/Data/CrossSections/Minschwaner/table-o2-sr-cross.dat'
		  print *,trim(fn)
		  OPEN (jchan,FILE=TRIM(fn),STATUS='old')
!
!        read data
		  DO jt = 1 , NTSR
			  DO jw = 1 , NWAVESR
				 READ (jchan,*) iwave_sr(jw) , t_sr(jt) , sigma_sr(jw,jt)
			  ENDDO
		  ENDDO
!
		  CLOSE (jchan)
!
		  if(verb%lverb)then
			  DO jt = 1 , NTSR
				 DO jw = 1 , NWAVESR
				   PRINT '(i3,1x,1p,6g13.6)' , iwave_sr(jw) , t_sr(jt) , sigma_sr(jw,jt)
				 ENDDO
			  ENDDO
		  endif
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      twant=273.15
!
!     For a given temperature find the absorption crosssection.
      ao2sr=0.
!
!     Lyman-alpha, First interval
      IF ( Tc>=0.0 ) THEN
		ao2sr(1) = 2.115E-18*tc**(-0.1145)
          if(verb%lverb)print *,jw,iwave_sr(jw),value
      endif
!
!     For each wavelength interval that we have O2 Schumann-Runge cross 
!     section data for interpolate to required tempearture.
      do jw=1,NWAVESR
          part=sigma_sr(jw,:)
          call FIND1D(twant,Value,t_sr,part,NTSR)
          ao2sr(iwave_sr(jw))=value
          if(verb%lverb)print *,jw,iwave_sr(jw),value
      enddo
!
!-----------------------------------------------------------------------
!
      END PROGRAM SR
!
!-----------------------------------------------------------------------
!
      SUBROUTINE FIND1D(Xwant,Value,X,Array,NX)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: NX
      REAL :: Value , Xwant
      REAL , DIMENSION(NX) :: Array , X
      INTENT (IN) Array , NX
      INTENT (OUT) Value
!
! Local variables
!
      REAL :: c
      INTEGER :: jx , jxp1 , m
!
!-----------------------------------------------------------------------
!
      CALL POS(X,NX,Xwant,jx)
      IF ( jx==0 ) THEN
         jx = 1
      ELSEIF ( jx>=NX ) THEN
         jx = NX - 1
      ENDIF
      jxp1 = jx + 1
!
      m = (Array(jx+1)-Array(jx))/(X(jx+1)-X(jx))
      c = Array(jx) - (m*X(jx))
      Value = Xwant*m + c
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FIND1D
!
!----------------------------------------------------------------------
!
      SUBROUTINE POS(Xx,N,X,J)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: J
      INTEGER :: N
      REAL :: X
      REAL , DIMENSION(N) :: Xx
      INTENT (IN) N , X , Xx
      INTENT (OUT) J
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
