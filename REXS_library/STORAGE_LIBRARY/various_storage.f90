!
!=======================================================================
!
MODULE MATRIX_COEF
!
!  This modules contains matrix coefficients
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : LINIMAX,LINFMAX
!
      IMPLICIT NONE
!
      COMPLEX (WP)              ::  COEF(5),ANG(5,3,3,LINIMAX,LINFMAX)
!
END MODULE MATRIX_COEF
!
!=======================================================================
!
MODULE EL_PH_CORRECTION
!
!  This module contains information about correction terms
!
!
      IMPLICIT NONE
!
      INTEGER                   ::  I_CORR
!
END MODULE EL_PH_CORRECTION

