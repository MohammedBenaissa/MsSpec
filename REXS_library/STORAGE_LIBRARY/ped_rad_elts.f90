!
!=======================================================================
!
!
!  These modules contain the radial matrix elements for:
!
!                              * PED
!                              * XAS
!                              * REXS
!
!
!=======================================================================
!
MODULE RADINT_R
!
!  This module contains the REGULAR radial integrals
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M,NSPIN2_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOR(NE_M,10,NSPIN2_M,5)
!
END MODULE RADINT_R
!
!=======================================================================
!
MODULE RADINT_I
!
!  This module contains the IRREGULAR radial integrals
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M,NSPIN2_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOI(NE_M,10,NSPIN2_M,5)
!
END MODULE RADINT_I
!
!=======================================================================
!
MODULE RAD_E1E2
!
!  This module contains the off-diagonal E1-E2 radial integrals
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOI_E1_E2(NE_M,2,3)
!
END MODULE RAD_E1E2
