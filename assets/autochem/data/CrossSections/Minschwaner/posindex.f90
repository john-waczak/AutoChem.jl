      SUBROUTINE POSINDEX(Xx,Index,N,X,J)
      IMPLICIT NONE
!
! Dummy arguments
!
      INTEGER :: J
      INTEGER :: N
      REAL :: X
      INTEGER , DIMENSION(N) :: Index
      REAL , DIMENSION(N) :: Xx
      INTENT (IN) Index , N , X , Xx
      INTENT (OUT) J
!
! Local variables
!
      INTEGER :: jl , jm , ju
!
!-----------------------------------------------------------------------
!
!     written by:   David Lary
!
!     started:      7/1/1993
!
!     last updated: 22/1/2004
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
!     based on Numerical Recipes, The art of scientific computing,
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
            IF ( Xx(Index(N))>Xx(Index(1)) .EQV. X>Xx(Index(jm)) ) THEN
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
      END SUBROUTINE POSINDEX
