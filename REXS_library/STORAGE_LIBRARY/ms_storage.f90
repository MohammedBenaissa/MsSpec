!
!=======================================================================
!
MODULE ROTATION
!
!  This module stores extrement values of the paths
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NL_M
!
      IMPLICIT NONE
!
      REAL (WP)              ::  RLM01(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
!
END MODULE ROTATION
!
!=======================================================================
!
MODULE SCAT_MAT
!
!  This module stores the scattering matrix
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NSPIN2_M,NLAMBDA_M,N_SCAT_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)           ::  F21(NSPIN2_M,NLAMBDA_M,NLAMBDA_M,N_SCAT_M)
!
END MODULE SCAT_MAT
!
!=======================================================================
!
MODULE TEMP_STORAGE_1
!
!  This module stores extrement values of the paths
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : N_SCAT_M
!
      IMPLICIT NONE
!
      REAL (WP)              ::  DW(0:N_SCAT_M)
!
      COMPLEX (WP)           ::  CEX(0:N_SCAT_M),CEXDW(0:N_SCAT_M)
!
END MODULE TEMP_STORAGE_1
