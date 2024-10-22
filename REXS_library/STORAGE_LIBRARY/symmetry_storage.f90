!
!=======================================================================
!
MODULE ROT_CUB
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NL_M
!
      IMPLICIT NONE
!
      REAL   (WP)       ::  R_PIS2(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
!
END MODULE ROT_CUB
!
!=======================================================================
!
MODULE SYM_OP
!
!
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
!      INTEGER           ::  IZ(64)
!
!      REAL   (WP)       ::  ZL(64)
!
!      COMPLEX (WP)      ::  ZM1(64),ZM2(64)
!
END MODULE SYM_OP
!
!=======================================================================
!
MODULE TAU_PROT
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER           ::  I_Z_P(NATP_M)
      INTEGER           ::  ISTEP_L(NATP_M)
      INTEGER           ::  ISTEP_M(NATP_M)
      INTEGER           ::  I_LM(NATP_M)
      INTEGER           ::  I_REL_MP(NATP_M)
!
      REAL   (WP)       ::  Z_L_P(NATP_M)
!
      COMPLEX (WP)      ::  Z_M_P(NATP_M)
!
END MODULE TAU_PROT
!
!=======================================================================
!
MODULE TAU_SYM_OP
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NAT_EQ_M,NATP_M
!
      IMPLICIT NONE
!
      INTEGER           ::  I_Z(NAT_EQ_M,NATP_M)
      INTEGER           ::  ISYM(NAT_EQ_M,NATP_M)
      INTEGER           ::  NSYM_G(48)
      INTEGER           ::  NSIZE_GR
!
      REAL   (WP)       ::  Z_L(NAT_EQ_M,NATP_M)
      REAL   (WP)       ::  S_INV(64,3,3)
!
      COMPLEX (WP)      ::  Z_M1(NAT_EQ_M,NATP_M)
      COMPLEX (WP)      ::  Z_M2(NAT_EQ_M,NATP_M)
!
END MODULE TAU_SYM_OP
