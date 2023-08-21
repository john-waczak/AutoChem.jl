!-------------------------------------------------------------------------------
!
      SUBROUTINE MOMENT(Nmax,Dat,N,Ave,Adev,Sdev,Var,Skew,Curt)
      IMPLICIT NONE
!
! Dummy arguments
!
      REAL :: Adev , Ave , Curt , Sdev , Skew , Var
      INTEGER :: N , Nmax
      REAL , DIMENSION(Nmax) :: Dat
!
! Local variables
!
      REAL :: ep , p , s
      INTEGER :: j
!
!-------------------------------------------------------------------------------
!
      IF ( N==1 ) THEN
!	   write(6,*)'Moment n=1'
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
!            WRITE (6,*) 'Sdev:' , Sdev
!            WRITE (6,*) 'Var :' , Var
!            WRITE (6,*)
!     &                'no skew or kurtosis when zero variance in moment'
         ENDIF
      ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE MOMENT
