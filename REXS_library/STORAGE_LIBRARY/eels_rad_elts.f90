!
!=======================================================================
!
!
!  These modules contain the radial matrix elements for:
!
!                              * EELS
!                              * (E,2E)
!
!
!=======================================================================
!
MODULE RADINT_L1
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M,NL_M,NATCLU_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOL_00_RD(NE_M,0:NL_M,0:NL_M,0:NL_M,0:2*NL_M)
!
END MODULE RADINT_L1
!
!=======================================================================
!
MODULE RADINT_L2
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M,NL_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOL_00_RE(NE_M,0:NL_M,0:NL_M,0:NL_M,0:2*NL_M)
!
END MODULE RADINT_L2
!
!=======================================================================
!
MODULE RADINT_L3
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M,NATCLU_M,NL_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOL_0J_RD(NE_M,NATCLU_M,0:NL_M,0:NL_M,0:NL_M,0:2*NL_M,0:2*NL_M)
!
END MODULE RADINT_L3
!
!=======================================================================
!
MODULE RADINT_L4
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M,NATCLU_M,NL_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOL_0J_RE(NE_M,NATCLU_M,0:NL_M,0:NL_M,0:NL_M,0:2*NL_M,0:2*NL_M)
!
END MODULE RADINT_L4
!
!=======================================================================
!
MODULE RADINT_L5
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M,NL_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOL_ID(NE_M,0:NL_M,0:NL_M,0:NL_M,0:2*NL_M)
!
END MODULE RADINT_L5
!
!=======================================================================
!
MODULE RADINT_L6
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NE_M,NL_M,NATCLU_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOL_IE(NE_M,0:NL_M,0:NL_M,0:NL_M,0:2*NL_M)
!
END MODULE RADINT_L6
