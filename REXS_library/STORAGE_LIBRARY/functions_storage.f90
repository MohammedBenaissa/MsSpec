!
!=======================================================================
!
MODULE LOGAMMA
!
!  This module stores the logarithm of the Gamma function
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,      ONLY : L_MAX,N_GAUNT
!
      IMPLICIT NONE
!
      REAL (WP)             ::  GLG(0:N_GAUNT)
      REAL (WP)             ::  GLD(0:N_GAUNT,2)
!
END MODULE LOGAMMA
!
!=======================================================================
!
MODULE STORE_COEF
!
!  This modules contains stored coefficients
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NL_M,N_GAUNT
!
      IMPLICIT NONE
!
      REAL (WP)             ::  CF(0:2*NL_M-2,0:2*NL_M-2,0:2*NL_M-2)
      REAL (WP)             ::  EXPF(0:2*NL_M-2,0:2*NL_M-2)
      REAL (WP)             ::  EXPR(0:2*NL_M-2,0:2*NL_M-2)
      REAL (WP)             ::  EXPF2(0:2*NL_M-2,0:2*NL_M-2)
      REAL (WP)             ::  FSQ(0:2*NL_M-2)
!
END MODULE STORE_COEF
!
!=======================================================================
!
MODULE CLEBSCH_GORDAN
!
!  This module stores the Clebsch-Gordan coefficients
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : LI_M,NCG_M
!
      IMPLICIT NONE
!
      REAL (WP)             ::  CG(-LI_M:LI_M+1,2,2)
      REAL (WP)             ::  CG_S(2,2,2)
      REAL (WP)             ::  CGA(0:NCG_M,-NCG_M:NCG_M,0:NCG_M,-NCG_M:NCG_M,0:2*NCG_M)
!
END MODULE CLEBSCH_GORDAN
!
!=======================================================================
!
MODULE GAUNT_C
!
!  This module stores the Gaunt coefficients
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : LINMAX,N_GAUNT
!
      IMPLICIT NONE
!
      REAL (WP)             ::  GNT(0:N_GAUNT,LINMAX,LINMAX)
!
END MODULE GAUNT_C
