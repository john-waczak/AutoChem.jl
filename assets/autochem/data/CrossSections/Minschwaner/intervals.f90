      PROGRAM INTERVALS
!
!-----------------------------------------------------------------------
!
!     A program to write out my wavelength intervals, start, stop and 
!     centre in nm, wave numbers from an input in microns.
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
! PARAMETER definitions
!
      INTEGER , PARAMETER :: NWAVE = 203
!
! Local variables
!
      CHARACTER(400) :: fn
      INTEGER :: jchan , jw
      REAL , DIMENSION(NWAVE) :: wave_micron_c , wave_micron_hi ,       &
     &                           wave_micron_lo , wave_nm_c ,           &
     &                           wave_nm_hi , wave_nm_lo , wave_wn_c ,  &
     &                           wave_wn_hi , wave_wn_lo
!
!-----------------------------------------------------------------------
!
      jchan = 1
!
!-----------------------------------------------------------------------
!
!     Read in wavelength intervals in microns
      fn = 'intervals-microns.dat'
      PRINT * , TRIM(fn)
      OPEN (jchan,FILE=TRIM(fn),STATUS='old')
!
!     Skip header
      DO jw = 1 , 2
         READ (jchan,*)
      ENDDO
!
!     read data
      DO jw = 1 , NWAVE
         READ (jchan,*) wave_micron_c(jw) , wave_micron_lo(jw) ,        &
     &                  wave_micron_hi(jw)
      ENDDO
!
      CLOSE (jchan)
!
!-----------------------------------------------------------------------
!
      wave_nm_c = 1E3*wave_micron_c
      wave_nm_lo = 1E3*wave_micron_lo
      wave_nm_hi = 1E3*wave_micron_hi
!
!     Write out wavelength intervals in nms
      fn = 'intervals-nm.dat'
      PRINT * , TRIM(fn)
      OPEN (jchan,FILE=TRIM(fn),STATUS='unknown')
!
!     Write header
      WRITE (jchan,*) 'wavelength bands (nm) '
      WRITE (jchan,*) '       Center WL       Start WL       Stop WL'
!
!     write data
      DO jw = 1 , NWAVE
         WRITE (jchan,FMT='(1p,3g15.6)') wave_nm_c(jw) , wave_nm_lo(jw) &
     &          , wave_nm_hi(jw)
      ENDDO
!
      CLOSE (jchan)
!
!-----------------------------------------------------------------------
!
      wave_wn_c = 1./(1E-4*wave_micron_c)
      wave_wn_lo = 1./(1E-4*wave_micron_lo)
      wave_wn_hi = 1./(1E-4*wave_micron_hi)
!
!     Write out wavelength intervals in nms
      fn = 'intervals-cm.dat'
      PRINT * , TRIM(fn)
      OPEN (jchan,FILE=TRIM(fn),STATUS='unknown')
!
!     Write header
      WRITE (jchan,*) 'wavelength bands (cm-1) '
      WRITE (jchan,*) '       Center WL       Start WL       Stop WL'
!
!     write data
      DO jw = 1 , NWAVE
         WRITE (jchan,FMT='(1p,3g15.6)') wave_wn_c(jw) , wave_wn_lo(jw) &
     &          , wave_wn_hi(jw)
      ENDDO
!
      CLOSE (jchan)
!
!-----------------------------------------------------------------------
!
      END PROGRAM INTERVALS
