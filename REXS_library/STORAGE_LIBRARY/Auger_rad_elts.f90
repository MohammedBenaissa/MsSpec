!
!=======================================================================
!
MODULE RADINT_A
!
!  This module contains the REGULAR radial integrals
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY :  NT_M,NSPIN2_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  RHOA(0:NT_M,0:40,2,NSPIN2_M,5)
!
END MODULE RADINT_A

