!
!=======================================================================
!
MODULE EXTREMES
!
!  This module stores extrement values of the paths
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : N_SCAT_M
!
      IMPLICIT NONE
!
      INTEGER                ::  IREF
!
      REAL (WP)              ::  FMIN(0:N_SCAT_M),FMAX(0:N_SCAT_M)
!
END MODULE EXTREMES
!
!=======================================================================
!
MODULE PATH_INFO
!
!  This module stores various path information
!
      USE DIMENSION_CODE,     ONLY : N_SCAT_M
!
      IMPLICIT NONE
!
      INTEGER                ::  NPATH(0:N_SCAT_M),NPATH2(0:N_SCAT_M)
      INTEGER                ::  NTHOF
      INTEGER                ::  NPMA(0:N_SCAT_M),NPMI(0:N_SCAT_M)
!
END MODULE PATH_INFO
!
!=======================================================================
!
MODULE PRINT_PATH_INFO
!
!  This module stores various path printing information
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NPATH_M,N_SCAT_M
!
      IMPLICIT NONE
!
      INTEGER                ::  JON(NPATH_M),JPON(NPATH_M,N_SCAT_M)
!
      REAL (WP)              ::  FMN(NPATH_M),PATH(NPATH_M),DMN(NPATH_M)
!
END MODULE PRINT_PATH_INFO
!
!=======================================================================
!
MODULE TESTPA
!
!  This module stores various path information
!
      USE DIMENSION_CODE,     ONLY : N_SCAT_M
!
      IMPLICIT NONE
!
      INTEGER                ::  IT(0:N_SCAT_M),IN(0:N_SCAT_M)
      INTEGER                ::  IJ
!
END MODULE TESTPA
!
!=======================================================================
!
MODULE TESTPB
!
!  This module stores various path information
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)              ::  TH01,PHI01
!
      COMPLEX (WP)           ::  RHO01
!
END MODULE TESTPB


