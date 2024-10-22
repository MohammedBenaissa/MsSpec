MODULE SYMMETRIZATION
!
!  This module provides the tools to symmetrized the cluster
!
      USE ACCURACY_REAL
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE SYM_CLUS(CUNIT,COORD,CHEM_OLD,NZ_AT,NTYP,IPHA)
!
!  This subroutine reorganizes the cluster into equivalence classes
!    taking advantage of its symmetry properties. The corresponding
!    symmetry operations are stored for further use.
!
!
!   Author   :  M. Gavaza      : 15 Feb 2000
!   Modified :  D. Sébilleau   : 05 Nov 2003
!
!
!
!                                       Last modified (DS) : 18 Jun 2021
!
!
      USE DIMENSION_CODE,         ONLY : NATCLU_M,LI_M
!
      USE REAL_NUMBERS,           ONLY : ZERO,ONE,HALF
      USE COMPLEX_NUMBERS,        ONLY : ONEC,IC
      USE PI_ETC,                 ONLY : PI
      USE SQUARE_ROOTS,           ONLY : SQR2,SQR3
!
      USE ATOMS
      USE CLUSTER_LIMITS
      USE CLUSTER_COORD
      USE OUTFILES
      USE OUTUNITS
!
      USE CALC_TYPE,              ONLY : I_TEST
      USE ROT_CUB
      !USE SYM_OP
      USE TAU_PROT
      USE TAU_SYM_OP
      USE TESTS
!
      USE SORT1
      USE SORT2
!
!!
      USE WIGNER_ROTATIONS
!!
      IMPLICIT NONE
!
      CHARACTER (LEN =  2), INTENT(INOUT) ::  CHEM_OLD(NATCLU_M)
!
      CHARACTER (LEN =  3)             ::  NAME_G(0:32)
      CHARACTER (LEN =  4)             ::  ANG_ROT(4)
      CHARACTER (LEN =  5)             ::  NAME_S(64)
      CHARACTER (LEN =  8)             ::  AT_ADD(2)
!
      INTEGER, INTENT(INOUT)           ::  NTYP(NATCLU_M)
      INTEGER, INTENT(IN)              ::  IPHA
!
      INTEGER                          ::  JAT,JAT1,JAT2
      INTEGER                          ::  I_ROT,JROT,I_SB
      INTEGER                          ::  NAB,NABOR,NABOR1,NABOR2
      INTEGER                          ::  JTYP_C_M,NEQAT_M
      INTEGER                          ::  NASHUFF,NABORS
      INTEGER                          ::  JNUM,JTYP,NA1P1,N1,N2
      INTEGER                          ::  NZT,NROTCLUS
      INTEGER                          ::  NSIZE_C
      INTEGER                          ::  JNUM1,JNUM2
      INTEGER                          ::  I_THREE,I_INCRG,I_INVG
      INTEGER                          ::  JSYM
      INTEGER                          ::  NEQATS,I_CALC_ROT
      INTEGER                          ::  JGROUP,JPROT,JPLAN
      INTEGER                          ::  NSIZE_RE
      INTEGER                          ::  JSYM1,JSYM2,JSYM3
      INTEGER                          ::  JS1,JS2,JS3,JS4
      INTEGER                          ::  I_EQ,I_NEW,I_OLD,J_OLD
      INTEGER                          ::  NB_GR,NB_GR1,JG,IWRONG
      INTEGER                          ::  JR,J_RE,JCHANGE
      INTEGER                          ::  JS,JS_MAX,JS_RE,JS_M
      INTEGER                          ::  JS_GR1,JS_GR2,JS_NE1,JS_NE2
      INTEGER                          ::  ICHECK_S
      INTEGER                          ::  IOLD,JG_MIN,I_MINUS
      INTEGER                          ::  JGROUP_M,I_END,IROT,NAT
      INTEGER                          ::  J_AT,J_AT_P,NZ_NEW,I_IN_OK
      INTEGER                          ::  I_BACK,NAT_NEW
      INTEGER                          ::  JDIST,JCLASS,NMAX
      INTEGER                          ::  NAP,JSYM_0,NA
      INTEGER                          ::  NSYM_P,NSYM_PT,JAT_P
      INTEGER                          ::  ISL,ISM
      INTEGER                          ::  NZ_OLD,NDIST,NB_G
      INTEGER                          ::  I,J,K,N,JOP
      INTEGER                          ::  NSYM_M(64),CONT_G(48,32),SIZE_G(32)
      INTEGER                          ::  INV_P(NATP_M,64),INV_PT(NATP_M,64)
      INTEGER                          ::  NAT_SYM(4)
      INTEGER                          ::  I_SET(NATP_M),GR(NATP_M)
      INTEGER                          ::  NEQAT(NATCLU_M),NZ_AT(NATCLU_M)
      INTEGER                          ::  NBRZ(NATCLU_M,NATCLU_M)
      INTEGER                          ::  IUSED(NATCLU_M),NGAIN_B(NATP_M)
      INTEGER                          ::  NQAT(NATCLU_M),N_NEW_OLD(NATCLU_M)
      INTEGER                          ::  I_33(32),IS_WRONG(64)
      INTEGER                          ::  OLD_G(32),NEW_G(32)
! Symmetry ops ----------------------
     INTEGER           ::  IZ(64)
!
      REAL   (WP)       ::  ZL(64)

      COMPLEX (WP)      ::  ZM1(64),ZM2(64)
!-----------------------------------


!
      REAL (WP), INTENT(IN)            ::  CUNIT
      REAL (WP), INTENT(INOUT)         ::  COORD(3,NATCLU_M)
!
      REAL (WP), PARAMETER             ::  SQ2O2 = SQR2 * HALF
      REAL (WP), PARAMETER             ::  SQ3O2 = SQR3 * HALF
!
      REAL (WP)                        ::  SSMALL,RAB,RABT
      REAL (WP)                        ::  AD,S1,S2
      REAL (WP)                        ::  D1,D2
      REAL (WP)                        ::  X_NEW,Y_NEW,Z_NEW
      REAL (WP)                        ::  Z_DIF,Z_TYP
      REAL (WP)                        ::  X,Y,Z,X_A,Y_A,Z_A
      REAL (WP)                        ::  GAIN_B,GAIN_G
      REAL (WP)                        ::  X1,Y1,Z1
      REAL (WP)                        ::  COORD1(3,NATCLU_M)
      REAL (WP)                        ::  SYM_AT1(3,NATCLU_M)
      REAL (WP)                        ::  DNBR(NATCLU_M,NATCLU_M)
      REAL (WP)                        ::  DIST0(NATCLU_M),DIST1(NATCLU_M)
      REAL (WP)                        ::  S_M(3,3,64)
      REAL (WP)                        ::  S_MUL(3,3)
      REAL (WP)                        ::  PIS2
!
      REAL (WP)                        ::  SQRT,FLOAT
!
      COMPLEX (WP), PARAMETER          ::  ONEC_ = - ONEC
      COMPLEX (WP), PARAMETER          ::  IC_   = - IC
      COMPLEX (WP), PARAMETER          ::  JC1   = (HALF,SQ3O2)
      COMPLEX (WP), PARAMETER          ::  JC1_  = - JC1
      COMPLEX (WP), PARAMETER          ::  JC2   = JC1 * JC1
      COMPLEX (WP), PARAMETER          ::  JC4   = JC2 * JC2
      COMPLEX (WP), PARAMETER          ::  JC5   = JC4 * JC1
      COMPLEX (WP), PARAMETER          ::  JC2_  = ONEC_ * JC2
      COMPLEX (WP), PARAMETER          ::  JC4_  = ONEC_ * JC4
      COMPLEX (WP), PARAMETER          ::  JC5_  = ONEC_ * JC5
!
!-------------------------------
! negative values
      COMPLEX (WP), PARAMETER          ::  ONE_    = - ONE
      COMPLEX (WP), PARAMETER          ::  HALF_   = - HALF
      COMPLEX (WP), PARAMETER          ::  SQ3O2_  = - SQ3O2
!
!----------------------------------
      LOGICAL                          ::  EQUIV,MATCH
!
!  Matrices of the 64 symmetry operations in the form S_M(I,J)
!     I is chosen as the line index and J as the column one
!
!  Matrices for the rotations of the group Oh
!
!     The greek indices for the cubic rotations have been written as :
!
!                     alpha ----> l
!                     beta  ----> m
!                     gamma ----> n
!                     delta ----> o
!
      DATA ((S_M(I,J,1), J=1,3), I=1,3)                                    & ! E
                           / ONE,ZERO,ZERO,ZERO, ONE,ZERO,ZERO,ZERO, ONE/    !
!
      DATA ((S_M(I,J,2), J=1,3), I=1,3)                                    & ! C2X
                           / ONE,ZERO,ZERO,ZERO,ONE_,ZERO,ZERO,ZERO,ONE_/    !
!
      DATA ((S_M(I,J,3), J=1,3), I=1,3)                                    & ! C2Y
                           /ONE_,ZERO,ZERO,ZERO, ONE,ZERO,ZERO,ZERO,ONE_/    !
!
      DATA ((S_M(I,J,4), J=1,3), I=1,3)                                    & ! C2Z
                           /ONE_,ZERO,ZERO,ZERO,ONE_,ZERO,ZERO,ZERO, ONE/    !
!
      DATA ((S_M(I,J,5), J=1,3), I=1,3)                                    & ! C3l
                           /ZERO,ZERO, ONE, ONE,ZERO,ZERO,ZERO, ONE,ZERO/    !
!
      DATA ((S_M(I,J,6), J=1,3), I=1,3)                                    & ! C3m
                           /ZERO,ONE_,ZERO,ZERO,ZERO, ONE,ONE_,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,7), J=1,3), I=1,3)                                    & ! C3n
                           /ZERO,ZERO,ONE_, ONE,ZERO,ZERO,ZERO,ONE_,ZERO/    !
!
      DATA ((S_M(I,J,8), J=1,3), I=1,3)                                    & ! C3o
                           /ZERO,ONE_,ZERO,ZERO,ZERO,ONE_, ONE,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,9), J=1,3), I=1,3)                                    & ! C3l2
                           /ZERO, ONE,ZERO,ZERO,ZERO, ONE, ONE,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,10), J=1,3), I=1,3)                                   & ! C3m2
                           /ZERO,ZERO,ONE_,ONE_,ZERO,ZERO,ZERO, ONE,ZERO/    !
!
      DATA ((S_M(I,J,11), J=1,3), I=1,3)                                   & ! C3n2
                           /ZERO, ONE,ZERO,ZERO,ZERO,ONE_,ONE_,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,12), J=1,3), I=1,3)                                   & ! C3o2
                           /ZERO,ZERO, ONE,ONE_,ZERO,ZERO,ZERO,ONE_,ZERO/    !
!
      DATA ((S_M(I,J,13), J=1,3), I=1,3)                                   & ! C4X
                           / ONE,ZERO,ZERO,ZERO,ZERO,ONE_,ZERO, ONE,ZERO/    !
!
      DATA ((S_M(I,J,14), J=1,3), I=1,3)                                   & ! C4Y
                           /ZERO,ZERO, ONE,ZERO, ONE,ZERO,ONE_,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,15), J=1,3), I=1,3)                                   & ! C4Z
                           /ZERO,ONE_,ZERO, ONE,ZERO,ZERO,ZERO,ZERO, ONE/    !
!
      DATA ((S_M(I,J,16), J=1,3), I=1,3)                                   & ! C4X3
                           / ONE,ZERO,ZERO,ZERO,ZERO, ONE,ZERO,ONE_,ZERO/    !
!
      DATA ((S_M(I,J,17), J=1,3), I=1,3)                                   & !  C4Y3
                           /ZERO,ZERO,ONE_,ZERO, ONE,ZERO, ONE,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,18), J=1,3), I=1,3)                                   & ! C4Z3
                           /ZERO, ONE,ZERO,ONE_,ZERO,ZERO,ZERO,ZERO, ONE/    !
!
      DATA ((S_M(I,J,19), J=1,3), I=1,3)                                   & ! C2a
                           /ZERO, ONE,ZERO, ONE,ZERO,ZERO,ZERO,ZERO,ONE_/    !
!
      DATA ((S_M(I,J,20), J=1,3), I=1,3)                                   & ! C2b
                           /ZERO,ONE_,ZERO,ONE_,ZERO,ZERO,ZERO,ZERO,ONE_/    !
!
      DATA ((S_M(I,J,21), J=1,3), I=1,3)                                   & ! C2c
                           /ZERO,ZERO, ONE,ZERO,ONE_,ZERO, ONE,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,22), J=1,3), I=1,3)                                   & ! C2d
                           /ZERO,ZERO,ONE_,ZERO,ONE_,ZERO,ONE_,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,23), J=1,3), I=1,3)                                   & ! C2e
                           /ONE_,ZERO,ZERO,ZERO,ZERO,ONE_,ZERO,ONE_,ZERO/    !
!
      DATA ((S_M(I,J,24), J=1,3), I=1,3)                                   & ! C2f
                           /ONE_,ZERO,ZERO,ZERO,ZERO, ONE,ZERO, ONE,ZERO/    !
!
!  Matrices for the rotations of the group D6h
!
      DATA ((S_M(I,J,25), J=1,3), I=1,3)                                   & ! C3Z
                 /HALF_,SQ3O2_,ZERO, SQ3O2,HALF_,ZERO,ZERO,ZERO, ONE/        !
!
      DATA ((S_M(I,J,26), J=1,3), I=1,3)                                   & ! C3Z2
                 /HALF_, SQ3O2,ZERO,SQ3O2_,HALF_,ZERO,ZERO,ZERO, ONE/        !
!
      DATA ((S_M(I,J,27), J=1,3), I=1,3)                                   & ! C6Z
                 / HALF,SQ3O2_,ZERO, SQ3O2, HALF,ZERO,ZERO,ZERO, ONE/        !
!
      DATA ((S_M(I,J,28), J=1,3), I=1,3)                                   & ! C6Z5
                 / HALF, SQ3O2,ZERO,SQ3O2_, HALF,ZERO,ZERO,ZERO, ONE/        !
!
      DATA ((S_M(I,J,29), J=1,3), I=1,3)                                   & ! C2A
                 /HALF_,SQ3O2_,ZERO,SQ3O2_, HALF,ZERO,ZERO,ZERO,ONE_/        !
!
      DATA ((S_M(I,J,30), J=1,3), I=1,3)                                   & ! C2B
                 /HALF_, SQ3O2,ZERO, SQ3O2, HALF,ZERO,ZERO,ZERO,ONE_/        !
!
      DATA ((S_M(I,J,31), J=1,3), I=1,3)                                   & ! C2C
                 / HALF,SQ3O2_,ZERO,SQ3O2_,HALF_,ZERO,ZERO,ZERO,ONE_/        !
!
      DATA ((S_M(I,J,32), J=1,3), I=1,3)                                   & ! C2D
                 / HALF, SQ3O2,ZERO, SQ3O2,HALF_,ZERO,ZERO,ZERO,ONE_/        !
!
!  Matrices for the roto-inversions of the group Oh
!
      DATA ((S_M(I,J,33), J=1,3), I=1,3)                                   & ! I
                           /ONE_,ZERO,ZERO,ZERO,ONE_,ZERO,ZERO,ZERO,ONE_/    !
!
      DATA ((S_M(I,J,34), J=1,3), I=1,3)                                   & ! IC2X
                           /ONE_,ZERO,ZERO,ZERO, ONE,ZERO,ZERO,ZERO, ONE/    !
!
      DATA ((S_M(I,J,35), J=1,3), I=1,3)                                   & ! IC2Y
                           / ONE,ZERO,ZERO,ZERO,ONE_,ZERO,ZERO,ZERO, ONE/    !
!
      DATA ((S_M(I,J,36), J=1,3), I=1,3)                                   & ! IC2Z
                           / ONE,ZERO,ZERO,ZERO, ONE,ZERO,ZERO,ZERO,ONE_/    !
!
      DATA ((S_M(I,J,37), J=1,3), I=1,3)                                   & ! IC3l
                           /ZERO,ZERO,ONE_,ONE_,ZERO,ZERO,ZERO,ONE_,ZERO/    !
!
      DATA ((S_M(I,J,38), J=1,3), I=1,3)                                   & ! IC3m
                           /ZERO, ONE,ZERO,ZERO,ZERO,ONE_, ONE,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,39), J=1,3), I=1,3)                                   & ! IC3n
                           /ZERO,ZERO, ONE,ONE_,ZERO,ZERO,ZERO, ONE,ZERO/    !
!
      DATA ((S_M(I,J,40), J=1,3), I=1,3)                                   & ! IC3o
                           /ZERO, ONE,ZERO,ZERO,ZERO, ONE,ONE_,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,41), J=1,3), I=1,3)                                   & ! IC3l2
                           /ZERO,ONE_,ZERO,ZERO,ZERO,ONE_,ONE_,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,42), J=1,3), I=1,3)                                   & ! IC3m2
                           /ZERO,ZERO, ONE, ONE,ZERO,ZERO,ZERO,ONE_,ZERO/    !
!
      DATA ((S_M(I,J,43), J=1,3), I=1,3)                                   & ! IC3n2
                           /ZERO,ONE_,ZERO,ZERO,ZERO, ONE, ONE,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,44), J=1,3), I=1,3)                                   & ! IC3o2
                           /ZERO,ZERO,ONE_, ONE,ZERO,ZERO,ZERO, ONE,ZERO/    !
!
      DATA ((S_M(I,J,45), J=1,3), I=1,3)                                   & ! IC4X
                           /ONE_,ZERO,ZERO,ZERO,ZERO, ONE,ZERO,ONE_,ZERO/    !
!
      DATA ((S_M(I,J,46), J=1,3), I=1,3)                                   & ! IC4Y
                           /ZERO,ZERO,ONE_,ZERO,ONE_,ZERO, ONE,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,47), J=1,3), I=1,3)                                   & ! IC4Z
                           /ZERO, ONE,ZERO,ONE_,ZERO,ZERO,ZERO,ZERO,ONE_/    !
!
      DATA ((S_M(I,J,48), J=1,3), I=1,3)                                   & ! IC4X3
                           /ONE_,ZERO,ZERO,ZERO,ZERO,ONE_,ZERO, ONE,ZERO/    !
!
      DATA ((S_M(I,J,49), J=1,3), I=1,3)                                   & ! IC4Y3
                           /ZERO,ZERO, ONE,ZERO,ONE_,ZERO,ONE_,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,50), J=1,3), I=1,3)                                   & ! IC4Z3
                           /ZERO,ONE_,ZERO, ONE,ZERO,ZERO,ZERO,ZERO,ONE_/    !
!
      DATA ((S_M(I,J,51), J=1,3), I=1,3)                                   & ! IC2a
                           /ZERO,ONE_,ZERO,ONE_,ZERO,ZERO,ZERO,ZERO, ONE/    !
!
      DATA ((S_M(I,J,52), J=1,3), I=1,3)                                   & ! IC2b
                           /ZERO, ONE,ZERO, ONE,ZERO,ZERO,ZERO,ZERO, ONE/    !
!
      DATA ((S_M(I,J,53), J=1,3), I=1,3)                                   & ! IC2c
                           /ZERO,ZERO,ONE_,ZERO, ONE,ZERO,ONE_,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,54), J=1,3), I=1,3)                                   & ! IC2d
                           /ZERO,ZERO, ONE,ZERO, ONE,ZERO, ONE,ZERO,ZERO/    !
!
      DATA ((S_M(I,J,55), J=1,3), I=1,3)                                   & ! IC2e
                           / ONE,ZERO,ZERO,ZERO,ZERO, ONE,ZERO, ONE,ZERO/    !
!
      DATA ((S_M(I,J,56), J=1,3), I=1,3)                                   & ! IC2f
                           / ONE,ZERO,ZERO,ZERO,ZERO,ONE_,ZERO,ONE_,ZERO/    !
!
!  Matrices for the roto-inversions of the group D6h
!
      DATA ((S_M(I,J,57), J=1,3), I=1,3)                                   & ! IC3Z
                 / HALF, SQ3O2,ZERO,SQ3O2_, HALF,ZERO,ZERO,ZERO,ONE_/        !
      DATA ((S_M(I,J,58), J=1,3), I=1,3)                                   & ! IC3Z2
                 / HALF,SQ3O2_,ZERO, SQ3O2, HALF,ZERO,ZERO,ZERO,ONE_/        !
      DATA ((S_M(I,J,59), J=1,3), I=1,3)                                   & ! IC6Z
                 /HALF_, SQ3O2,ZERO,SQ3O2_,HALF_,ZERO,ZERO,ZERO,ONE_/        !
      DATA ((S_M(I,J,60), J=1,3), I=1,3)                                   & ! IC6Z5
                 /HALF_,SQ3O2_,ZERO, SQ3O2,HALF_,ZERO,ZERO,ZERO,ONE_/        !
      DATA ((S_M(I,J,61), J=1,3), I=1,3)                                   & ! IC2A
                 / HALF, SQ3O2,ZERO, SQ3O2,HALF_,ZERO,ZERO,ZERO, ONE/        !
      DATA ((S_M(I,J,62), J=1,3), I=1,3)                                   & ! IC2B
                 / HALF,SQ3O2_,ZERO,SQ3O2_,HALF_,ZERO,ZERO,ZERO, ONE/        !
      DATA ((S_M(I,J,63), J=1,3), I=1,3)                                   & ! IC2C
                 /HALF_, SQ3O2,ZERO, SQ3O2, HALF,ZERO,ZERO,ZERO, ONE/        !
      DATA ((S_M(I,J,64), J=1,3), I=1,3)                                   & ! IC2D
                 /HALF_,SQ3O2_,ZERO,SQ3O2_, HALF,ZERO,ZERO,ZERO, ONE/        !
!
!  For a given symmetry operation S, IZ can have three values :
!
!         IZ =  1  ---->  delta_{L, L'} in S_{L L'}
!         IZ =  0  ---->  delta_{l, l'} only in S_{L L'}
!         IZ = -1  ---->  delta_{L,-L'} in S_{L L'}
!
      DATA IZ / 1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 1, 0,       &
                0, 1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1,       &
                1,-1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,       &
                0, 1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1 /
!
!  For a given symmetry operation S, S_{L L'} is proportional to
!                      ZL**L and ZM1**M (and ZM2**M' if IZ=0)
!
      DATA ZL / ONE,ONE_,ONE_, ONE, ONE, ONE, ONE, ONE, ONE, ONE,     &
                ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE,ONE_,ONE_,     &
                ONE, ONE, ONE, ONE, ONE, ONE, ONE, ONE,ONE_,ONE_,     &
               ONE_,ONE_,ONE_, ONE, ONE,ONE_,ONE_,ONE_,ONE_,ONE_,     &
               ONE_,ONE_,ONE_,ONE_,ONE_,ONE_,ONE_,ONE_,ONE_,ONE_,     &
                ONE, ONE,ONE_,ONE_,ONE_,ONE_,ONE_,ONE_,ONE_,ONE_,     &
                ONE, ONE, ONE, ONE /
!
      DATA ZM1 /ONEC,ONEC,ONEC_,ONEC_,ONEC,IC_,ONEC_,IC,IC_,          &
               ONEC_,IC,ONEC,IC,ONEC,IC_,IC_,ONEC_,IC,                &
               IC_,IC,ONEC,ONEC_,IC,IC_,JC4,JC2,JC5,JC1,JC5_,         &
               JC1_,JC4_,JC2_,ONEC,ONEC,ONEC_,ONEC_,ONEC,IC_,         &
               ONEC_,IC,IC_,ONEC_,IC,ONEC,IC,ONEC,IC_,IC_,ONEC_,      &
               IC,IC_,IC,ONEC,ONEC_,IC,IC_,JC4,JC2,JC5,JC1,JC5_,      &
               JC1_,JC4_,JC2_/
      DATA ZM2 /ONEC,ONEC,ONEC,ONEC,IC_,ONEC,IC,ONEC_,ONEC_,          &
                IC_,ONEC,IC,IC_,ONEC,ONEC,IC,ONEC_,ONEC,              &
                ONEC,ONEC,ONEC_,ONEC,IC,IC_,ONEC,ONEC,ONEC,ONEC,      &
                ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,ONEC,IC_,ONEC,     &
                IC,ONEC_,ONEC_,IC_,ONEC,IC,IC_,ONEC,ONEC,ONEC,ONEC,   &
                ONEC,ONEC,ONEC,ONEC_,ONEC,IC,IC_,ONEC,ONEC,IC,ONEC_,  &
                ONEC,ONEC,ONEC,ONEC/
!
!  Name of the crystallographic point-groups
!
      DATA NAME_G /'---',                                             &
                   ' C1',' Ci',' C2','C1h','C2h',' D2','C2v','D2h',   &
                   ' C4',' S4','C4h',' D4','C4v','D2d','D4h',' C3',   &
                   'C3i',' D3','C3v','D3d',' C6','C3h','C6h',' D6',   &
                   'C6v','D3h','D6h','  T',' Th',' Td','  O',' Oh'/
!
!  Content of the crystallographic point-groups : CONT_G(JROT,JGROUP)
!
!     In some cases, two contents are given. They correspond to a rotation of
!         the x and y axes about 0z by pi/4 in the cube (D2,C2v,D2h,D2d)
!         or pi/6 in the hexagon (D3,C3v,D3d,D3h). In this case, x and y
!          are respectively transformed into a and b in the first case,
!           and D and A in the second. The cube is invariant by any pi/2
!             rotation and the hexagon by any pi/3 rotation about 0z.
!
!
      DATA (CONT_G(I,1), I=1,48) /1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,    & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! C1
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,2), I=1,48) /1,33,0,0,0,0,0,0,0,0,0,0,0,0,0,0,   & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! Ci
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,3), I=1,48) /1,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,    & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! C2
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,4), I=1,48) /1,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0,   & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! C1h
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,5), I=1,48) /1,4,33,36,0,0,0,0,0,0,0,0,0,0,0,0,  & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! C2h
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,6), I=1,48) /1,2,3,4,0,0,0,0,0,0,0,0,0,0,0,0,    & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! D2
                           0,0,0,0,0,0,0,0,0,0,0,0,1,4,19,20/           !
!
      DATA (CONT_G(I,7), I=1,48) /1,4,34,35,0,0,0,0,0,0,0,0,0,0,0,0,  & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! C2v
                           0,0,0,0,0,0,0,0,0,0,0,0,1,4,51,52/           !
!
      DATA (CONT_G(I,8), I=1,48) /1,2,3,4,33,34,35,36,0,0,0,0,0,0,0,0,& !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! D2h
                           0,0,0,0,0,0,0,0,1,4,19,20,33,36,51,52/       !
!
      DATA (CONT_G(I,9), I=1,48) /1,4,15,18,0,0,0,0,0,0,0,0,0,0,0,0,  & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! C4
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,10), I=1,48)/1,4,47,50,0,0,0,0,0,0,0,0,0,0,0,0,  & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! S4
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,11), I=1,48)/1,4,15,18,33,36,47,50,0,0,0,0,0,0,  & !
                           0,0,                                       & ! C4h
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,12), I=1,48)/1,2,3,4,15,18,19,20,0,0,0,0,0,0,0,0,& !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! D4
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,13), I=1,48)/1,4,15,18,34,35,51,52,0,0,0,0,0,0,  & !
                           0,0,                                       & ! C4v
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,14), I=1,48)/1,2,3,4,47,50,51,52,0,0,0,0,0,0,0,0,& !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! D2d
                           0,0,0,0,0,0,0,0,1,4,19,20,34,35,47,50/       !
!
      DATA (CONT_G(I,15), I=1,48)/1,2,3,4,15,18,19,20,33,34,35,36,47, & !
                           50,51,52,                                  & ! D4h
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,16), I=1,48)/1,25,26,0,0,0,0,0,0,0,0,0,0,0,0,0,  & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! C3
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,17), I=1,48)/1,25,26,33,57,58,0,0,0,0,0,0,0,0,   & !
                           0,0,                                       & ! C3i
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,18), I=1,48)/1,3,25,26,31,32,0,0,0,0,0,0,0,0,0,0,& !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! D3
                           0,0,0,0,0,0,0,0,0,0,1,2,25,26,29,30/         !
!
      DATA (CONT_G(I,19), I=1,48)/1,25,26,34,61,62,0,0,0,0,0,0,0,0,   & !
                           0,0,                                       & ! C3v
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,1,25,26,35,63,64/        !
!
      DATA (CONT_G(I,20), I=1,48)/1,3,25,26,31,32,33,35,57,58,63,64,  & !
                           0,0,0,0,                                   & ! D3d
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,1,2,25,26,29,30,33,34,57,58,61,62/   !
!
      DATA (CONT_G(I,21), I=1,48)/1,4,25,26,27,28,0,0,0,0,0,0,0,0,0,0,& !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! C6
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,22), I=1,48)/1,25,26,36,59,60,0,0,0,0,0,0,0,0,   & !
                           0,0,                                       & ! C3h
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,23), I=1,48)/1,4,25,26,27,28,33,36,57,58,59,60,  & !
                           0,0,0,0,                                   & ! C6h
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,24), I=1,48)/1,2,3,4,25,26,27,28,29,30,31,32,0,  & !
                           0,0,0,                                     & ! D6
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,25), I=1,48)/1,4,25,26,27,28,34,35,61,62,63,64,  & !
                           0,0,0,0,                                   & ! C6v
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,26), I=1,48)/1,3,25,26,31,32,34,36,59,60,61,62,  & !
                           0,0,0,0,                                   & ! D3h
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & !
                           0,0,0,0,1,2,25,26,29,30,35,36,59,60,63,64/   !
!
      DATA (CONT_G(I,27), I=1,48)/1,2,3,4,25,26,27,28,29,30,31,32,33, & !
                           34,35,36,                                  & ! D6h
                           57,58,59,60,61,62,63,64,0,0,0,0,0,0,0,0,   & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,28), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,0,0,0,0, & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,           & ! T
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,29), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,33,34,   & !
                           35,36,                                     & ! Th
                           37,38,39,40,41,42,43,44,0,0,0,0,0,0,0,0,   & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,30), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,45,46,   & !
                           47,48,                                     & ! Td
                           49,50,51,52,53,54,55,56,0,0,0,0,0,0,0,0,   & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,31), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,13,14,   & !
                           15,16,                                     & ! O
                           17,18,19,20,21,22,23,24,0,0,0,0,0,0,0,0,   & !
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/             !
!
      DATA (CONT_G(I,32), I=1,48)/1,2,3,4,5,6,7,8,9,10,11,12,13,14,   & !
                      15,16,                                          & ! Oh
                      17,18,19,20,21,22,23,24,33,34,35,36,37,38,39,40,& !
                      41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56/  !
!
!  Size of the point-groups
!
      DATA SIZE_G /1,2,2,2,4,4,4,8,4,4,8,8,8,8,16,3,                  &
                   6,6,6,12,6,6,12,12,12,12,24,12,24,24,24,48/
!
!  Groups containing the space inversion I
!
      DATA I_33   /0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,                   &
                   1,0,0,1,0,0,1,0,0,0,1,0,1,0,0,1/
!
!  Name of the symmetry operations
!
!  Note : the three-fold rotation axes alpha, beta, gamma and delta
!         of the cubic groups have been coded respectively
!         l, m, n and o
!
      DATA NAME_S /'    E','  C2X','  C2Y','  C2Z','  C3l','  C3m',   &
                   '  C3n','  C3o',' C3l2',' C3m2',' C3n2',' C3o2',   &
                   '  C4X','  C4Y','  C4Z',' C4X3',' C4Y3',' C4Z3',   &
                   '  C2a','  C2b','  C2c','  C2d','  C2e','  C2f',   &
                   '  C3Z',' C3Z2','  C6Z',' C6Z5','  C2A','  C2B',   &
                   '  C2C','  C2D','    I',' IC2X',' IC2Y',' IC2Z',   &
                   ' IC3l',' IC3m',' IC3n',' IC3o','IC3l2','IC3m2',   &
                   'IC3n2','IC3o2',' IC4X',' IC4Y',' IC4Z','IC4X3',   &
                   'IC4Y3','IC4Z3',' IC2a',' IC2b',' IC2c',' IC2d',   &
                   ' IC2e',' IC2f',' IC3Z','IC3Z2',' IC6Z','IC6Z5',   &
                   ' IC2A',' IC2B',' IC2C',' IC2D'/
!
      DATA ANG_ROT  /'PI/4','PI/2','PI/6','PI/3'/
!
      DATA AT_ADD   /'        ','<--- NEW'/
!
      DATA SSMALL   / 0.001E0_WP /
!
      PRINT 444                                                     !
!
      DO JS = 1, 64                                                 !
        IS_WRONG(JS) = 0                                            !
      END DO                                                        !
      IROT    = 0                                                   !
      NAT_NEW = NATCLU                                              !
!
      IF(I_GR == 2) THEN                                            !
        I_SB = 1                                                    !
      ELSE                                                          !
        I_SB = 0                                                    !
      END IF                                                        !
!
 111  DO JAT = 1, NAT_NEW                                           !
        NEQAT(JAT) = 0                                              !
      END DO                                                        !
!
!    Calculates the distance between the atoms and all
!    their neighbors
!
      DIST0(1) = ZERO                                               !
      DO JAT1 = 1, NAT_NEW                                          !
        INEW_AT(JAT1) = 0                                           !
        NABOR         = 0                                           !
        DO JAT2 = 1, NAT_NEW                                        !
          IF(JAT1 == JAT2) GO TO 10                                 !
          NABOR            = NABOR + 1                              !
          NBRZ(NABOR,JAT1) = NZ_AT(JAT2)                            !
          RAB = SQRT( (COORD(1,JAT1)-COORD(1,JAT2)) *             & !
                      (COORD(1,JAT1)-COORD(1,JAT2)) +             & !
                      (COORD(2,JAT1)-COORD(2,JAT2)) *             & !
                      (COORD(2,JAT1)-COORD(2,JAT2)) +             & !
                      (COORD(3,JAT1)-COORD(3,JAT2)) *             & !
                      (COORD(3,JAT1)-COORD(3,JAT2)) )               !
          IF((JAT2 > JAT1) .AND. (RAB < SSMALL)) GO TO 895          !
          IF(JAT1 == 1) DIST0(JAT2) = RAB                           !
          DNBR(NABOR,JAT1) = RAB                                    !
 10       CONTINUE                                                  !
        END DO                                                      !
      END DO                                                        !
!
!   Generates the set of class distances to the absorber
!
      CALL SORT_DEC(NAT_NEW,DIST0,NDIST,DIST1)                      !
!
!     Determines the prototypical and equivalent atoms:
!     two atoms are equivalent if they have the same type,
!     they are at the same distance from the emitter
!     and they have the same geometrical environment.
!
!     This part of the routine has been adapted from
!     R. Gunnella and C. R. Natoli.
!
      JTYP_C_M = 0                                                  !
      NEQAT_M  = 0                                                  !
      JTYP     = 2                                                  !
      NASHUFF  = 2                                                  !
      NABORS   = NAT_NEW - 1                                        !
!
      NATYP(1)     = 1                                              !
      NQAT(1)      = 1                                              !
      NCORR(1,1)   = 1                                              !
      NEQAT(1)     = 1                                              !
      NCHTYP(1)    = 1                                              !
      I_Z(1,1)     = 1                                              !
      Z_L(1,1)     = ONE                                            !
      Z_M1(1,1)    = ONEC                                           !
      Z_M2(1,1)    = ONEC                                           !
      SYM_AT(1,1)  = COORD(1,1)                                     !
      SYM_AT(2,1)  = COORD(2,1)                                     !
      SYM_AT(3,1)  = COORD(3,1)                                     !
      NZ_AT(1)     = NZ_AT(1)                                       !
      CHEM(1)      = CHEM_OLD(1)                                    !
      INEW_AT(1)   = 0                                              !
      N_NEW_OLD(1) = 1                                              !
!
      NATYP(2)     = 1                                              !
      NQAT(2)      = 2                                              !
      NCORR(1,2)   = 2                                              !
      NEQAT(2)     = 2                                              !
      NCHTYP(2)    = 2                                              !
      I_Z(1,2)     = 1                                              !
      Z_L(1,2)     = ONE                                            !
      Z_M1(1,2)    = ONEC                                           !
      Z_M2(1,2)    = ONEC                                           !
      SYM_AT(1,2)  = COORD(1,2)                                     !
      SYM_AT(2,2)  = COORD(2,2)                                     !
      SYM_AT(3,2)  = COORD(3,2)                                     !
      NZ_AT(2)     = NZ_AT(2)                                       !
      CHEM(2)      = CHEM_OLD(2)                                    !
      INEW_AT(2)   = 0                                              !
      N_NEW_OLD(2) = 2                                              !
!
      JNUM = 0                                                      !
      DO JAT1 = 2, NAT_NEW                                          !
        IF(JNUM > NEQAT_M) THEN                                     !
          NEQAT_M  = JNUM                                           !
          JTYP_C_M = JTYP                                           !
        END IF                                                      !
        JNUM  = 1                                                   !
        NA1P1 = JAT1 + 1                                            !
        IF(NA1P1 > NAT_NEW) GO TO 32                                !
        DO JAT2 = NA1P1 ,NAT_NEW                                    !
          IF(NZ_AT(JAT1) /= NZ_AT(JAT2)) GO TO 30                   !
          N1 = JAT1 - 1                                             !
          N2 = JAT2 - 1                                             !
          IF( ABS(DNBR(N1,1) - DNBR(N2,1) ) > SSMALL) GO TO 30      !
          DO NAB = 1, NABORS                                        !
            IUSED(NAB) = 0                                          !
          END DO                                                    !
          DO NABOR1 = 1, NABORS                                     !
            NZT   = NBRZ(NABOR1,JAT1)                               !
            RABT  = DNBR(NABOR1,JAT1)                               !
            EQUIV = .FALSE.                                         !
            DO NABOR2 = 1, NABORS                                   !
              IF(IUSED(NABOR2) == 1) GO TO 22                       !
              IF(NBRZ(NABOR2,JAT2) /= NZT) GO TO 22                 !
              IF( ABS(DNBR(NABOR2,JAT2) - RABT) > SSMALL ) GO TO 22 !
              EQUIV         = .TRUE.                                !
              IUSED(NABOR2) = 1                                     !
              GO TO 23                                              !
 22           CONTINUE                                              !
            END DO                                                  !
 23         IF(.NOT.EQUIV) GO TO 30                                 !
          END DO                                                    !
          IF(NEQAT(JAT1) == 0) THEN                                 !
            JTYP               = JTYP + 1                           !
            NASHUFF            = NASHUFF + 1                        !
            NCORR(JNUM,JTYP)   = NASHUFF                            !
            I_Z(JNUM,JTYP)     = 1                                  !
            Z_L(JNUM,JTYP)     = ONE                                !
            Z_M1(JNUM,JTYP)    = ONEC                               !
            Z_M2(JNUM,JTYP)    = ONEC                               !
            SYM_AT(1,NASHUFF)  = COORD(1,JAT1)                      !
            SYM_AT(2,NASHUFF)  = COORD(2,JAT1)                      !
            SYM_AT(3,NASHUFF)  = COORD(3,JAT1)                      !
            N_NEW_OLD(NASHUFF) = JAT1                               !
            NEQAT(JAT1)        = JTYP                               !
            NCHTYP(JTYP)       = NTYP(JAT1)                         !
            NZ_AT(NASHUFF)     = NZ_AT(JAT1)                        !
            CHEM(NASHUFF)      = CHEM_OLD(JAT1)                     !
            NQAT(NASHUFF)      = JTYP                               !
            IF(JAT1 > NATCLU) INEW_AT(NASHUFF) = 1                  !
          END IF                                                    !
          IF(NEQAT(JAT2) == 0) THEN                                 !
            JNUM               = JNUM + 1                           !
            NATYP(JTYP)        = JNUM                               !
            NASHUFF            = NASHUFF + 1                        !
            NCORR(JNUM,JTYP)   = NASHUFF                            !
            SYM_AT(1,NASHUFF)  = COORD(1,JAT2)                      !
            SYM_AT(2,NASHUFF)  = COORD(2,JAT2)                      !
            SYM_AT(3,NASHUFF)  = COORD(3,JAT2)                      !
            N_NEW_OLD(NASHUFF) = JAT2                               !
            NEQAT(JAT2)        = JTYP                               !
            NCHTYP(JTYP)       = NTYP(JAT2)                         !
            NZ_AT(NASHUFF)     = NZ_AT(JAT2)                        !
            CHEM(NASHUFF)      = CHEM_OLD(JAT2)                     !
            NQAT(NASHUFF)      = JTYP                               !
            IF(JAT2 > NATCLU) INEW_AT(NASHUFF) = 1                  !
         END IF                                                     !
30       CONTINUE                                                   !
        END DO                                                      !
32      IF(NEQAT(JAT1) == 0) THEN                                   !
          JTYP               = JTYP + 1                             !
          NATYP(JTYP)        = JNUM                                 !
          NASHUFF            = NASHUFF + 1                          !
          NCORR(JNUM,JTYP)   = NASHUFF                              !
          I_Z(JNUM,JTYP)     = 1                                    !
          Z_L(JNUM,JTYP)     = ONE                                  !
          Z_M1(JNUM,JTYP)    = ONEC                                 !
          Z_M2(JNUM,JTYP)    = ONEC                                 !
          SYM_AT(1,NASHUFF)  = COORD(1,JAT1)                        !
          SYM_AT(2,NASHUFF)  = COORD(2,JAT1)                        !
          SYM_AT(3,NASHUFF)  = COORD(3,JAT1)                        !
          N_NEW_OLD(NASHUFF) = JAT1                                 !
          NEQAT(JAT1)        = JTYP                                 !
          NCHTYP(JTYP)       = NTYP(JAT1)                           !
          NZ_AT(NASHUFF)     = NZ_AT(JAT1)                          !
          CHEM(NASHUFF)      = CHEM_OLD(JAT1)                       !
          NQAT(NASHUFF)      = JTYP                                 !
          IF(JAT1 > NATCLU) INEW_AT(NASHUFF) = 1                    !
        END IF                                                      !
        CONTINUE                                                    !
      END DO                                                        !
      N_PROT = JTYP                                                 !
!
!     Stop if the maximal number of prototypical and equivalent
!                atoms are not correctly dimensionned
!
      IF(N_PROT > NATP_M) THEN                                      !
        WRITE(IUO1,897) N_PROT                                      !
        STOP                                                        !
      END IF                                                        !
      IF(NEQAT_M > NAT_EQ_M) THEN                                   !
        WRITE(IUO1,898) NEQAT_M                                     !
        STOP                                                        !
      END IF                                                        !
!
      DO JAT = 1, NAT_NEW                                           !
        COORD1(1,JAT) = SYM_AT(1,JAT) - SYM_AT(1,1)                 !
        COORD1(2,JAT) = SYM_AT(2,JAT) - SYM_AT(2,1)                 !
        COORD1(3,JAT) = SYM_AT(3,JAT) - SYM_AT(3,1)                 !
      END DO                                                        !
!
!     Test of all symmetry operations for the largest
!              symmetry class
!
      NROTCLUS  = 0                                                 !
 556  NSYM_M(1) = 1                                                 !
      NSIZE_C   = 1                                                 !
      DO JROT = 2, 64                                               !
        DO JNUM1 = 1, NEQAT_M                                       !
          JAT1  = NCORR(JNUM1,JTYP_C_M)                             !
          MATCH = .FALSE.                                           !
          SYM_AT1(1,JAT1) = S_M(1,1,JROT) * COORD1(1,JAT1) +      & !
                            S_M(1,2,JROT) * COORD1(2,JAT1) +      & !
                            S_M(1,3,JROT) * COORD1(3,JAT1)          !
          SYM_AT1(2,JAT1) = S_M(2,1,JROT) * COORD1(1,JAT1) +      & !
                            S_M(2,2,JROT) * COORD1(2,JAT1) +      & !
                            S_M(2,3,JROT) * COORD1(3,JAT1)          !
          SYM_AT1(3,JAT1) = S_M(3,1,JROT) * COORD1(1,JAT1) +      & !
                            S_M(3,2,JROT) * COORD1(2,JAT1) +      & !
                            S_M(3,3,JROT) * COORD1(3,JAT1)          !
          DO JNUM2 = 1, NEQAT_M                                     !
            JAT2 = NCORR(JNUM2,JTYP_C_M)                            !
            AD   = ABS( COORD1(1,JAT2) - SYM_AT1(1,JAT1) ) +      & !
                   ABS( COORD1(2,JAT2) - SYM_AT1(2,JAT1) ) +      & !
                   ABS( COORD1(3,JAT2) - SYM_AT1(3,JAT1) )          !
            IF(AD < SSMALL) THEN                                    !
              MATCH = .TRUE.                                        !
            END IF                                                  !
          END DO                                                    !
          IF(.NOT.MATCH) GO TO 333                                  !
        END DO                                                      !
        NSIZE_C         = NSIZE_C + 1                               !
        NSYM_M(NSIZE_C) = JROT                                      !
 333    CONTINUE                                                    !
      END DO                                                        !
!
!      Test on all classes of the symmetry operations that work
!            for the largest class
!
      NSYM_G(1) = 1                                                 !
      NSIZE_GR  = 1                                                 !
      I_THREE   = 0                                                 !
      I_INCRG   = 0                                                 !
      I_INVG    = 0                                                 !
      DO JSYM = 2, NSIZE_C                                          !
        JROT = NSYM_M(JSYM)                                         !
        DO JTYP = 1, N_PROT                                         !
          NEQATS = NATYP(JTYP)                                      !
          DO JNUM1 = 1, NEQATS                                      !
            JAT1 = NCORR(JNUM1,JTYP)                                !
            MATCH = .FALSE.                                         !
            SYM_AT1(1,JAT1) = S_M(1,1,JROT) * COORD1(1,JAT1) +    & !
                              S_M(1,2,JROT) * COORD1(2,JAT1) +    & !
                              S_M(1,3,JROT) * COORD1(3,JAT1)        !
            SYM_AT1(2,JAT1) = S_M(2,1,JROT) * COORD1(1,JAT1) +    & !
                              S_M(2,2,JROT) * COORD1(2,JAT1) +    & !
                              S_M(2,3,JROT) * COORD1(3,JAT1)        !
            SYM_AT1(3,JAT1) = S_M(3,1,JROT) * COORD1(1,JAT1) +    & !
                              S_M(3,2,JROT) * COORD1(2,JAT1) +    & !
                              S_M(3,3,JROT) * COORD1(3,JAT1)        !
            DO JNUM2 = 1, NEQATS
              JAT2 = NCORR(JNUM2,JTYP)
              AD = ABS( COORD1(1,JAT2) - SYM_AT1(1,JAT1) ) +      & !
                   ABS( COORD1(2,JAT2) - SYM_AT1(2,JAT1) ) +      & !
                   ABS( COORD1(3,JAT2) - SYM_AT1(3,JAT1) )          !
              IF(AD < SSMALL) THEN                                  !
                MATCH = .TRUE.                                      !
              END IF                                                !
            END DO                                                  !
            IF(.NOT.MATCH) GO TO 335                                !
          END DO                                                    !
        END DO                                                      !
        IF(MATCH) THEN                                              !
          NSIZE_GR         = NSIZE_GR + 1                           !
          NSYM_G(NSIZE_GR) = JROT                                   !
          IF(JROT == 25) I_THREE = 1                                !
          IF(JROT == 33) I_INVG  = 1                                !
        END IF                                                      !
 335    CONTINUE                                                    !
      END DO                                                        !
      IF(NSIZE_GR > 48) GO TO 998                                   !
!
!      Set up of the parameter used to check if a larger
!         group containing the inversion I can be constructed
!         if the inversion is physically acceptable as an approximation
!
      IF((I_GR == 2) .AND. (I_INVG == 0)) THEN                      !
        IF(I_INV == 1) THEN                                         !
          I_INCRG = 1                                               !
        END IF                                                      !
      END IF                                                        !
!
      IF(NROTCLUS <= 2) THEN                                        !
        IF(NSIZE_GR > 1) THEN                                       !
          WRITE(IUO1,699) NSIZE_GR                                  !
        ELSE                                                        !
          WRITE(IUO1,698) NSIZE_GR                                  !
        END IF                                                      !
        IF(NSIZE_GR == 1) THEN                                      !
          WRITE(IUO1,705) NAME_S(NSYM_G(1))                         !
        ELSE IF(NSIZE_GR == 2) THEN                                 !
          WRITE(IUO1,704) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2))       !
        ELSE IF(NSIZE_GR == 3) THEN                                 !
          WRITE(IUO1,703) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2)),    & !
                          NAME_S(NSYM_G(3))                         !
        ELSE IF(NSIZE_GR == 4) THEN                                 !
          WRITE(IUO1,702) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2)),    & !
                          NAME_S(NSYM_G(3)),NAME_S(NSYM_G(4))       !
        ELSE IF(NSIZE_GR == 6) THEN                                 !
          WRITE(IUO1,701) (NAME_S(NSYM_G(JROT)), JROT=1,6)          !
        ELSE IF(NSIZE_GR >= 8) THEN                                 !
          WRITE(IUO1,700) (NAME_S(NSYM_G(JROT)), JROT=1,NSIZE_GR)   !
        END IF                                                      !
        IF(NROTCLUS > 0) THEN                                       !
          IF(I_THREE == 0) THEN                                     !
            WRITE(IUO1,706) ANG_ROT(NROTCLUS)                       !
          ELSE                                                      !
            WRITE(IUO1,706) ANG_ROT(NROTCLUS+2)                     !
          END IF                                                    !
        END IF                                                      !
      END IF                                                        !
!
!      Finding the cluster's symmetry group
!
 555  JGROUP = 0                                                    !
      IF(NSIZE_GR == 1) THEN                                        !
        JGROUP = 1                                                  !
      ELSE IF(NSIZE_GR == 2) THEN                                   !
        IF(NSYM_G(2) == 33) THEN                                    !
          JGROUP = 2                                                !
        ELSE IF(NSYM_G(2) == 4) THEN                                !
          JGROUP = 3                                                !
        ELSE IF(NSYM_G(2) == 36) THEN                               !
          JGROUP = 4                                                !
        END IF                                                      !
      ELSE IF(NSIZE_GR == 3) THEN                                   !
        JGROUP = 16                                                 !
      ELSE IF(NSIZE_GR == 4) THEN                                   !
        IF(NSYM_G(3) == 33) THEN                                    !
          JGROUP = 5                                                !
        ELSE IF((NSYM_G(3) == 3) .OR. (NSYM_G(3) == 19)) THEN       !
          JGROUP = 6                                                !
        ELSE IF((NSYM_G(3) == 34) .OR. (NSYM_G(3) == 51)) THEN      !
          JGROUP = 7                                                !
        ELSE IF(NSYM_G(3) == 15) THEN                               !
          JGROUP = 9                                                !
        ELSE IF(NSYM_G(3) == 47) THEN                               !
          JGROUP = 10                                               !
        END IF                                                      !
      ELSE IF(NSIZE_GR == 6) THEN                                   !
        IF(NSYM_G(4) == 26) THEN                                    !
          JGROUP = 18                                               !
        ELSE IF((NSYM_G(4) == 34) .OR. (NSYM_G(4) == 35)) THEN      !
          JGROUP = 19                                               !
        ELSE IF(NSYM_G(4) == 33) THEN                               !
          JGROUP = 17                                               !
        ELSE IF(NSYM_G(4) == 26) THEN                               !
          JGROUP = 21                                               !
        ELSE IF(NSYM_G(4) == 36) THEN                               !
          JGROUP = 22                                               !
        END IF                                                      !
      ELSE IF(NSIZE_GR == 8) THEN                                   !
        IF(NSYM_G(4) == 33) THEN                                    !
          IF(NSYM_G(8) == 50) THEN                                  !
            JGROUP = 11                                             !
          ELSE                                                      !
            JGROUP = 8                                              !
          END IF                                                    !
        ELSE                                                        !
          IF(NSYM_G(5) == 15) THEN                                  !
            JGROUP = 12                                             !
          ELSE                                                      !
            IF(NSYM_G(3) == 15) THEN                                !
              JGROUP = 13                                           !
            ELSE IF((NSYM_G(3) == 3).OR.(NSYM_G(3) == 19)) THEN     !
              JGROUP = 14                                           !
            END IF                                                  !
          END IF                                                    !
        END IF                                                      !
      ELSE IF(NSIZE_GR == 12) THEN                                  !
        IF(NSYM_G(5) == 5) THEN                                     !
          JGROUP = 28                                               !
        ELSE IF(NSYM_G(7) == 33) THEN                               !
          IF(NSYM_G(12) == 60) THEN                                 !
            JGROUP = 23                                             !
          ELSE                                                      !
            JGROUP = 20                                             !
          END IF                                                    !
        ELSE                                                        !
          IF(NSYM_G(3) == 3) THEN                                   !
            JGROUP = 24                                             !
          ELSE IF(NSYM_G(9) == 59) THEN                             !
            JGROUP = 26                                             !
          ELSE IF(NSYM_G(5) == 27) THEN                             !
            JGROUP = 25                                             !
          END IF                                                    !
        END IF                                                      !
      ELSE IF(NSIZE_GR == 16) THEN                                  !
        JGROUP = 15                                                 !
      ELSE IF(NSIZE_GR == 24) THEN                                  !
        IF(NSYM_G(17) == 57) THEN                                   !
          JGROUP = 27                                               !
        ELSE IF(NSYM_G(17) == 17) THEN                              !
          JGROUP = 31                                               !
        ELSE IF(NSYM_G(17) == 49) THEN                              !
          JGROUP = 30                                               !
        ELSE IF(NSYM_G(17) == 37) THEN                              !
          JGROUP = 29                                               !
        END IF                                                      !
      ELSE IF(NSIZE_GR == 48) THEN                                  !
        JGROUP = 32                                                 !
      END IF                                                        !
!
      IF(JGROUP > 0) THEN                                           !
        WRITE(IUO1,886)                                             !
        WRITE(IUO1,880) NAME_G(JGROUP)                              !
      ELSE                                                          !
!
!  If no group is found, the cluster is rotated by pi/4 ,or pi/6 for
!    hexagonal structures (i. e. when C3z=25 is present), around Oz
!    at most twice to account for possible misnaming of the Ox axis.
!    Then, the search for the point-group is restarted.
!
        IF(NROTCLUS < 2) THEN                                       !
          NROTCLUS = NROTCLUS + 1                                   !
          DO JPROT = 1 ,N_PROT                                      !
            NEQATS = NATYP(JPROT)                                   !
            DO JNUM = 1, NEQATS                                     !
              JAT=NCORR(JNUM,JPROT)                                 !
              X = COORD1(1,JAT)                                     !
              Y = COORD1(2,JAT)                                     !
              IF(I_THREE == 0) THEN                                 !
                COORD1(1,JAT) = SQ2O2 * (X - Y)                     !
                COORD1(2,JAT) = SQ2O2 * (X + Y)                     !
              ELSE                                                  !
                COORD1(1,JAT) = SQ3O2 * X - HALF   * Y              !
                COORD1(2,JAT) = HALF  * X +  SQ3O2 * Y              !
              END IF                                                !
              COORD1(3,JAT) = COORD1(3,JAT)                         !
              SYM_AT(1,JAT) = COORD1(1,JAT) + SYM_AT(1,1)           !
              SYM_AT(2,JAT) = COORD1(2,JAT) + SYM_AT(2,1)           !
              SYM_AT(3,JAT) = COORD1(3,JAT) + SYM_AT(3,1)           !
            END DO                                                  !
          END DO                                                    !
          GO TO 556                                                 !
        ELSE                                                        !
          IF((JGROUP == 0) .AND. (NROTCLUS == 2))  WRITE(IUO1,881)  !
        END IF                                                      !
      END IF                                                        !
      IF(JGROUP >= 28) WRITE(IUO1,697)                              !
      IF((JGROUP > 0) .AND. (I_INCRG == 1)) WRITE(IUO1,722)         !
!
!  Recovery of the original cluster when no group has been found
!
      IF((NROTCLUS == 2) .AND. (JGROUP == 0)) THEN                  !
        DO JPROT = 1, N_PROT                                        !
          NEQATS = NATYP(JPROT)                                     !
          DO JNUM = 1, NEQATS                                       !
            JAT = NCORR(JNUM,JPROT)                                 !
            X   = COORD1(1,JAT)                                     !
            Y   = COORD1(2,JAT)                                     !
            IF(I_THREE == 0) THEN                                   !
              COORD1(1,JAT) =   Y                                   !
              COORD1(2,JAT) = - X                                   !
            ELSE                                                    !
              COORD1(1,JAT) =   HALF  * X + SQ3O2 * Y               !
              COORD1(2,JAT) = - SQ3O2 * X + HALF  * Y               !
            END IF                                                  !
            COORD1(3,JAT) = COORD1(3,JAT)                           !
            SYM_AT(1,JAT) = COORD1(1,JAT) + SYM_AT(1,1)             !
            SYM_AT(2,JAT) = COORD1(2,JAT) + SYM_AT(2,1)             !
            SYM_AT(3,JAT) = COORD1(3,JAT) + SYM_AT(3,1)             !
          END DO                                                    !
        END DO                                                      !
        NROTCLUS = NROTCLUS + 1                                     !
        GO TO 556                                                   !
      END IF                                                        !
!
!  If still no group is found, or if the group can be augmented by
!   the inversion (as an approximation), check of the other symmetries that
!   should be present to account for group properties (i. e. those
!   obtained by multiplication of those that are effectively present)
!
      IF((JGROUP == 0) .OR. (I_INCRG == 1)) THEN                    !
        IF(I_INCRG == 1) THEN                                       !
          NSIZE_GR         = NSIZE_GR + 1                           !
          NSYM_G(NSIZE_GR) = 33                                     !
        END IF                                                      !
        NSIZE_RE = NSIZE_GR                                         !
 553    I_NEW    = 0                                                !
        DO JSYM1 = 2, NSIZE_RE                                      !
          JS1 = NSYM_G(JSYM1)                                       !
          DO JSYM2 = JSYM1,NSIZE_RE                                 !
            JS2 = NSYM_G(JSYM2)                                     !
            DO I = 1 ,3                                             !
              DO J = 1, 3                                           !
                S_MUL(I,J) = ZERO                                   !
                DO K = 1, 3                                         !
                  S_MUL(I,J) = S_MUL(I,J) +                       & !
                               S_M(I,K,JS1) * S_M(K,J,JS2)          !
                END DO                                              !
              END DO                                                !
            END DO                                                  !
            I_EQ = 0                                                !
            DO JSYM3 = 1, NSIZE_RE                                  !
              JS3   = NSYM_G(JSYM3)                                 !
              I_OLD = 0                                             !
              DO I = 1, 3                                           !
                DO J = 1, 3                                         !
                  S1 = S_MUL(I,J)                                   !
                  S2 = S_M(I,J,JS3)                                 !
                  IF(ABS(S1 - S2) < SSMALL) I_OLD = I_OLD + 1       !
                END DO                                              !
              END DO                                                !
              IF(I_OLD == 9) I_EQ = 1                               !
            END DO                                                  !
            IF(I_EQ == 0) THEN                                      !
              I_NEW = I_NEW + 1                                     !
              J_RE  = NSIZE_GR + I_NEW                              !
              DO JS4 = 2 ,64                                        !
                I_OLD = 0                                           !
                DO I = 1 ,3                                         !
                  DO J = 1, 3                                       !
                    S1 = S_MUL(I,J)                                 !
                    S2 = S_M(I,J,JS4)                               !
                    IF(ABS(S1 - S2) < SSMALL) I_OLD = I_OLD + 1     !
                  END DO                                            !
                END DO                                              !
                IF(I_OLD == 9) THEN                                 !
                  NSYM_G(J_RE) = JS4                                !
                  NSIZE_RE     = NSIZE_RE + 1                       !
                END IF                                              !
              END DO                                                !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
        IF(I_NEW > 0) GO TO 553                                     !
!
        IF(NSIZE_RE > NSIZE_GR) THEN                                !
          WRITE(IUO1,696) NSIZE_RE                                  !
          IF(NSIZE_RE == 2) THEN                                    !
            WRITE(IUO1,704) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2))     !
          ELSE IF(NSIZE_RE == 3) THEN                               !
            WRITE(IUO1,703) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2)),  & !
                            NAME_S(NSYM_G(3))                       !
          ELSE IF(NSIZE_RE == 4) THEN                               !
            WRITE(IUO1,702) NAME_S(NSYM_G(1)),NAME_S(NSYM_G(2)),  & !
                            NAME_S(NSYM_G(3)),NAME_S(NSYM_G(4))     !
          ELSE IF(NSIZE_RE == 6) THEN                               !
            WRITE(IUO1,701) (NAME_S(NSYM_G(JROT)), JROT=1,6)        !
          ELSE IF(NSIZE_RE >= 8) THEN                               !
              WRITE(IUO1,700) (NAME_S(NSYM_G(JROT)), JROT=1,NSIZE_RE)
          END IF                                                    !
        END IF                                                      !
      END IF                                                        !
!
!  If no group has been found, it means that the cluster does not
!   retain the whole symmetry of the crystal. When I_GR = 1 or 2,
!   the largest group consistent with the symmetry operations
!   found is sought. Then a new cluster is built to support
!             all the symmetries of this group.
!
      I_BACK = 0                                                    !
      IF(((JGROUP == 0) .OR. (I_INCRG == 1)) .AND. (I_GR >= 1)) THEN!
!
!  Search for the different point-groups containing the NSIZE_RE
!   symmetries found. If the cluster can not support the
!   inversion I (I_INV=0), a test is made on the content
!   of the point-groups to suppress those that contain I.
!
        WRITE(IUO1,707)                                             !
        WRITE(IUO1,708)                                             !
!
!  Input cluster
!
        NB_GR = 0                                                   !
        DO JG = 1, 32                                               !
          IF((I_INV == 0) .AND. (I_33(JG) == 1)) GO TO 223          !
          JS_MAX = SIZE_G(JG)                                       !
          IF(JS_MAX < NSIZE_RE) GO TO 223                           !
!
          ICHECK_S = 0                                              !
          DO JSYM = 1, NSIZE_RE                                     !
            JS_RE = NSYM_G(JSYM)                                    !
!
            DO JS_GR1 = 1, JS_MAX                                   !
              JS_NE1 = CONT_G(JS_GR1,JG)                            !
              IWRONG = 0                                            !
              IF(IROT > 0) THEN                                     !
                DO JS = 1, IROT                                     !
                  IF(JS_NE1 == IS_WRONG(JS)) THEN                   !
                    IWRONG = 1                                      !
                  END IF                                            !
                END DO                                              !
              END IF                                                !
              IF(JS_NE1 == JS_RE) THEN                              !
                ICHECK_S = ICHECK_S + 1                             !
              END IF                                                !
            END DO                                                  !
!
          END DO                                                    !
!
          IF(IWRONG == 1) GO TO 223                                 !
          IF(ICHECK_S == NSIZE_RE) THEN                             !
            NB_GR        = NB_GR + 1                                !
            OLD_G(NB_GR) = JG                                       !
          END IF                                                    !
 223      CONTINUE                                                  !
!
        END DO                                                      !
        NB_GR1 = NB_GR                                              !
        WRITE(IUO1,709) (NAME_G(OLD_G(J)), J=1,NB_GR1)              !
        WRITE(IUO1,710)                                             !
!
!  Input cluster rotated
!
        DO JG = 1, 32                                               !
          IF((I_INV == 0) .AND. (I_33(JG) == 1)) GO TO 225          !
          IOLD = 0                                                  !
          DO J_OLD = 1 ,NB_GR1                                      !
            IF(OLD_G(J_OLD) == JG) THEN                             !
               IOLD = IOLD + 1                                      !
            END IF                                                  !
          END DO                                                    !
          IF(IOLD /= 0) GO TO 225                                   !
          JS_MAX = SIZE_G(JG)                                       !
          IF(JS_MAX < NSIZE_RE) GO TO 225                           !
!
          ICHECK_S = 0                                              !
          DO JSYM = 1, NSIZE_RE                                     !
            JS_RE = NSYM_G(JSYM)                                    !
!
            DO JS_GR1 = 1, JS_MAX                                   !
              JS_GR2 = 49 - JS_GR1                                  !
              JS_NE2 = CONT_G(JS_GR2,JG)                            !
              IWRONG = 0                                            !
              IF(IROT > 0) THEN                                     !
                DO JS = 1, IROT                                     !
                  IF(JS_NE2 == IS_WRONG(JS)) THEN                   !
                    IWRONG = 1                                      !
                  END IF                                            !
                END DO                                              !
              END IF                                                !
              IF(JS_NE2 == JS_RE) THEN                              !
                ICHECK_S = ICHECK_S + 1                             !
              END IF                                                !
            END DO                                                  !
!
          END DO                                                    !
!
          IF(IWRONG == 1) GO TO 225                                 !
          IF(ICHECK_S == NSIZE_RE) THEN                             !
            NB_GR        = NB_GR + 1                                !
            OLD_G(NB_GR) = JG                                       !
          END IF                                                    !
 225      CONTINUE                                                  !
!
        END DO                                                      !
        IF(NB_GR1 < NB_GR) THEN                                     !
          WRITE(IUO1,713)                                           !
          WRITE(IUO1,711) (NAME_G(OLD_G(J)), J=NB_GR1+1,NB_GR)      !
          WRITE(IUO1,712)                                           !
        END IF                                                      !
        IF(I_INV == 0) THEN                                         !
          IF(I_GR == 2) THEN                                        !
            WRITE(IUO1,721)                                         !
          ELSE                                                      !
            WRITE(IUO1,726)                                         !
          END IF                                                    !
        END IF                                                      !
        CALL SORT_INC_I(NB_GR,OLD_G,NB_G,NEW_G)                     !
!
!  Construction of the augmented cluster consistent with the group JGROUP_M.
!    Note that this group must be consistent with the original cluster.
!    This is checked through 3 criteria :
!
!         1. No atom must be generated above the surface plane
!
!         2. If an atom generated by a symmetry of JGROUP_M coincides
!             with one of the original cluster, it must have the same
!             atomic number
!
!         3. Every new atom generated by a symmetry of JGROUP_M must
!             be outside the original cluster except those meeting
!             the previous criterion
!
!  When one of this criteria is not satisfied, it means that the group
!   JGROUP_M can not accomodate the original cluster. Hence, a smaller
!   group must be looked for
!
!  An extra criterion is used when I_GR = 1 :
!
!         4. No surface atom can be transformed into a bulk atom
!              and likewise no bulk atom into a surface atom
!
!
        JG_MIN  = NEW_G(1)                                          !
        I_MINUS  = 0                                                !
 222    JGROUP_M = NEW_G(NB_G-I_MINUS)                              !
        IF(JGROUP_M < JG_MIN) THEN                                  !
          IF(I_GR >= 1) THEN                                        !
            I_GR = 0                                                !
            WRITE(IUO1,723)                                         !
            GO TO 111                                               !
          ELSE                                                      !
            WRITE(IUO1,719)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
        JS_M = SIZE_G(JGROUP_M)                                     !
        WRITE(IUO1,714) NAME_G(JGROUP_M)                            !
        I_END = 0                                                   !
        IROT  = 0                                                   !
        NAT   = NAT_NEW                                             !
!
        DO JS = 2, JS_M                                             !
          JROT = CONT_G(JS,JGROUP_M)                                !
!
          DO J_AT = 2, NAT                                          !
            IF(I_END == 1) GO TO 224                                !
            X_NEW = S_M(1,1,JROT) * (COORD(1,J_AT) - COORD(1,1)) +& !
                    S_M(1,2,JROT) * (COORD(2,J_AT) - COORD(2,1)) +& !
                    S_M(1,3,JROT) * (COORD(3,J_AT) - COORD(3,1)) +& !
                                     COORD(1,1)                     !
            Y_NEW = S_M(2,1,JROT) * (COORD(1,J_AT) - COORD(1,1)) +& !
                    S_M(2,2,JROT) * (COORD(2,J_AT) - COORD(2,1)) +& !
                    S_M(2,3,JROT) * (COORD(3,J_AT) - COORD(3,1)) +& !
                                     COORD(2,1)                     !
            Z_NEW = S_M(3,1,JROT) * (COORD(1,J_AT) - COORD(1,1)) +& !
                    S_M(3,2,JROT) * (COORD(2,J_AT) - COORD(2,1)) +& !
                    S_M(3,3,JROT) * (COORD(3,J_AT) - COORD(3,1)) +& !
                                     COORD(3,1)                     !
!
!  Check for criterion 1
!
            IF(Z_NEW > VALZ_MAX) THEN                               !
              WRITE(IUO1,715) (COORD(J,J_AT), J=1,3),X_NEW,Y_NEW, & !
                                                     Z_NEW          !
              WRITE(IUO1,716) NAME_S(JROT),VALZ_MAX                 !
              IROT           = IROT + 1                             !
              IS_WRONG(IROT) = JROT                                 !
            END IF                                                  !
            IF(IROT > 0) THEN                                       !
              I_END = 1                                             !
              GO TO 224                                             !
            END IF                                                  !
            NZ_NEW  = NZ_AT(J_AT)                                   !
            I_OLD   = 0                                             !
            I_IN_OK = 0                                             !
!
!  Check for criterion 2
!
            DO J_AT_P = 2, NAT                                      !
              D2 = (X_NEW - COORD(1,J_AT_P)) *                    & !
                   (X_NEW-COORD(1,J_AT_P))   +                    & !
                   (Y_NEW - COORD(2,J_AT_P)) *                    & !
                   (Y_NEW-COORD(2,J_AT_P))   +                    & !
                   (Z_NEW - COORD(3,J_AT_P)) *                    & !
                   (Z_NEW-COORD(3,J_AT_P))                          !
              NZ_OLD = NZ_AT(J_AT_P)                                !
              IF(D2 < SSMALL) THEN                                  !
                I_OLD = I_OLD + 1                                   !
                IF(NZ_NEW /= NZ_OLD) THEN                           !
                  IROT           = IROT + 1                         !
                  IS_WRONG(IROT) = JROT                             !
                  WRITE(IUO1,715) (COORD(J,J_AT), J=1,3),         & !
                                   X_NEW,Y_NEW,Z_NEW                !
                  WRITE(IUO1,717)  CHEM_OLD(J_AT),J_AT,           & !
                                   NAME_S(JROT),CHEM_OLD(J_AT_P), & !
                                   J_AT_P                           !
                ELSE                                                !
                  I_IN_OK = 1                                       !
                END IF                                              !
              END IF                                                !
            END DO                                                  !
            IF(I_IN_OK == 1) GO TO 226                              !
!
!  Check for criterion 3
!
            DO JPLAN = 1, NPLAN                                     !
              IF(ABS(Z_NEW-VAL(JPLAN)) < SSMALL) THEN               !
                IF(X_NEW > X_MIN(JPLAN)) THEN                       !
                  IF(X_NEW < X_MAX(JPLAN)) THEN                     !
                    IF(Y_NEW > Y_MIN(JPLAN)) THEN                   !
                      IF(Y_NEW < Y_MAX(JPLAN)) THEN                 !
                        IROT           = IROT + 1                   !
                        IS_WRONG(IROT) = JROT                       !
                        WRITE(IUO1,715) (COORD(J,J_AT), J=1,3),   & !
                                         X_NEW,Y_NEW,Z_NEW          !
                        WRITE(IUO1,720) NAME_S(JROT)                !
                      END IF                                        !
                    END IF                                          !
                  END IF                                            !
                END IF                                              !
              END IF                                                !
            END DO                                                  !
!
!  Check for criterion 4
!
           IF(I_SB == 0) THEN                                       !
            Z_DIF = ABS(Z_NEW - COORD(3,J_AT))                      !
            Z_TYP = ABS(Z_NEW - VALZ_MAX)                           !
            IF(Z_DIF > SSMALL) THEN                                 !
              WRITE(IUO1,715) (COORD(J,J_AT), J=1,3),X_NEW,Y_NEW, & !
                                                     Z_NEW          !
              IF(Z_TYP < SSMALL) THEN                               !
                WRITE(IUO1,725) NAME_S(JROT),I_GR                   !
              ELSE                                                  !
                WRITE(IUO1,724) NAME_S(JROT),I_GR                   !
              END IF                                                !
              IROT           = IROT + 1                             !
              IS_WRONG(IROT) = JROT                                 !
            END IF                                                  !
            IF(IROT > 0) THEN                                       !
              I_END = 1                                             !
              GO TO 224                                             !
            END IF                                                  !
           END IF                                                   !
!
 226       IF(I_OLD == 0) THEN                                      !
              NAT = NAT + 1                                         !
              IF(NAT > NATCLU_M) THEN                               !
                WRITE(IUO1,718) NAT                                 !
                STOP                                                !
              END IF                                                !
              COORD(1,NAT)  = X_NEW                                 !
              COORD(2,NAT)  = Y_NEW                                 !
              COORD(3,NAT)  = Z_NEW                                 !
              VALZ(NAT)     = Z_NEW                                 !
              CHEM_OLD(NAT) = CHEM_OLD(J_AT)                        !
              NZ_AT(NAT)    = NZ_AT(J_AT)                           !
              NTYP(NAT)      = NTYP(J_AT)                           !
            END IF                                                  !
 224        CONTINUE                                                !
          END DO                                                    !
!
        END DO                                                      !
!
        I_BACK = 1                                                  !
        IF(IROT > 0) THEN                                           !
          I_MINUS = I_MINUS + 1                                     !
          GO TO 222                                                 !
        END IF                                                      !
      END IF                                                        !
      IF(I_BACK == 1) THEN                                          !
        NAT_NEW = NAT                                               !
        GO TO 111                                                   !
      END IF                                                        !
!
!  Writes the classes of atoms by increasing distance
!
      WRITE(IUO1,888)                                               !
      DO JDIST = NDIST, 1, -1                                       !
        DO JTYP = 1, N_PROT                                         !
          NMAX = NATYP(JTYP)                                        !
          DO JCLASS = 1, NMAX                                       !
            N  = NCORR(JCLASS,JTYP)                                 !
            D1 = SQRT( COORD1(1,N) * COORD1(1,N) +                & !
                       COORD1(2,N) * COORD1(2,N) +                & !
                       COORD1(3,N) * COORD1(3,N) )                  !
            IF(ABS(D1-DIST1(JDIST)) < SSMALL) THEN                  !
              WRITE(IUO1,557) N,SYM_AT(1,N),SYM_AT(2,N),          & !
                              SYM_AT(3,N),NQAT(N),JCLASS,CHEM(N), & !
                              DIST1(JDIST),AT_ADD(INEW_AT(N)+1)     !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
!
!  Writes the augmented symmetrized cluster into the output
!          file OUTFILE4 for further use if necessary
!
      OPEN(UNIT=IUO4, FILE=OUTFILE4, STATUS='UNKNOWN')              !
      WRITE(IUO4,778) IPHA                                          !
      DO JTYP = 1, N_PROT                                           !
        NMAX = NATYP(JTYP)                                          !
        DO JCLASS = 1, NMAX                                         !
          N = NCORR(JCLASS,JTYP)                                    !
          X = SYM_AT(1,N) / CUNIT                                   !
          Y = SYM_AT(2,N) / CUNIT                                   !
          Z = SYM_AT(3,N) / CUNIT                                   !
          IF(IPHA == 0) THEN                                        !
            WRITE(IUO4,125) N,CHEM(N),NZ_AT(N),X,Y,Z,JTYP           !
          ELSE IF(IPHA == 1) THEN                                   !
            WRITE(IUO4,779) N,CHEM(N),NZ_AT(N),X,Y,Z,JTYP           !
          ELSE IF(IPHA == 2) THEN                                   !
            WRITE(IUO4,*) NZ_AT(N),X,Y,Z,JTYP                       !
          END IF                                                    !
        END DO                                                      !
      END DO                                                        !
      CLOSE(IUO4)                                                   !
!
!  Construction of the symmetry operations list.
!    Associates the appropriate symmetry operation
!         for all equivalent atoms
!
      IF(IPRINT >= 2) WRITE(IUO1,887)                               !
 438  CONTINUE                                                      !
      I_CALC_ROT = 0                                                !
      DO JTYP = 1, N_PROT                                           !
        ISYM(1,JTYP) = 1                                            !
        NEQATS       = NATYP(JTYP)                                  !
        IF(NEQATS == 1) THEN                                        !
          IF(JTYP > N_PROT) GO TO 555                               !
          GO TO 338                                                 !
        END IF                                                      !
        NAP = NCORR(1,JTYP)                                         !
        DO JNUM = 1, NEQATS                                         !
          IF(JNUM >= 2) THEN                                        !
            JSYM_0 = 2                                              !
          ELSE                                                      !
            JSYM_0 = 1                                              !
          END IF                                                    !
          NA = NCORR(JNUM,JTYP)                                     !
          X  = COORD1(1,NAP)                                        !
          Y  = COORD1(2,NAP)                                        !
          Z  = COORD1(3,NAP)                                        !
          DO JSYM = JSYM_0, NSIZE_GR                                !
            JROT           = NSYM_G(JSYM)                           !
            SYM_AT1(1,NAP) = S_M(1,1,JROT) * X +                  & !
                             S_M(1,2,JROT) * Y +                  & !
                             S_M(1,3,JROT) * Z                      !
            SYM_AT1(2,NAP) = S_M(2,1,JROT) * X +                  & !
                             S_M(2,2,JROT) * Y +                  & !
                             S_M(2,3,JROT) * Z                      !
            SYM_AT1(3,NAP) = S_M(3,1,JROT) * X +                  & !
                             S_M(3,2,JROT) * Y +                  & !
                             S_M(3,3,JROT) * Z                      !
            AD             = ABS(COORD1(1,NA) - SYM_AT1(1,NAP)) + & !
                             ABS(COORD1(2,NA) - SYM_AT1(2,NAP)) + & !
                             ABS(COORD1(3,NA) - SYM_AT1(3,NAP))     !
            IF(AD < SSMALL)THEN                                     !
              ISYM(JNUM,JTYP) = JROT                                !
              I_Z(JNUM,JTYP)  = IZ(JROT)                            !
              Z_L(JNUM,JTYP)  = ZL(JROT)                            !
              Z_M1(JNUM,JTYP) = ZM1(JROT)                           !
              Z_M2(JNUM,JTYP) = ZM2(JROT)                           !
              IF(IZ(JROT) == 0) THEN                                !
                I_CALC_ROT = I_CALC_ROT + 1                         !
              END IF                                                !
              GO TO 404                                             !
            END IF                                                  !
          END DO                                                    !
          IF(ISYM(JNUM,JTYP) == 0) THEN                             !
            N_PROT          = N_PROT + 1                            !
            NCHTYP(N_PROT)  = NCHTYP(JTYP)                          !
            NCORR(1,N_PROT) = NA                                    !
            NATYP(N_PROT)   = 1                                     !
            I_Z(1,N_PROT)   = 1                                     !
            Z_L(1,N_PROT)   = ONE                                   !
            Z_M1(1,N_PROT) = ONEC                                   !
            Z_M2(1,N_PROT) = ONEC                                   !
            DO JCHANGE = JNUM, NEQATS-1                             !
              NCORR(JCHANGE,JTYP) = NCORR(JCHANGE+1,JTYP)           !
            END DO                                                  !
            NATYP(JTYP) = NATYP(JTYP) - 1                           !
            GO TO 438                                               !
          END IF                                                    !
 404      CONTINUE                                                  !
          IF((IPRINT >= 2) .AND. (NSIZE_GR > 1)) THEN               !
            JR = ISYM(JNUM,JTYP)                                    !
            WRITE(IUO1,849) JTYP,JNUM,NCORR(JNUM,JTYP),           & !
                            NAME_S(JR),JR,Z_L(JNUM,JTYP),         & !
                            Z_M1(JNUM,JTYP),Z_M2(JNUM,JTYP),      & !
                            I_Z(JNUM,JTYP)                          !
          END IF                                                    !
        END DO                                                      !
        WRITE(IUO1,*) '     '                                       !
 338    CONTINUE                                                    !
      END DO                                                        !
!
      GAIN_G = FLOAT(NAT_NEW) / FLOAT(N_PROT)                       !
      WRITE(IUO1,854) GAIN_G                                        !
!
!     Test of the symmetry operations leaving the prototypical
!       atoms invariant. Associates the apropriate symmetry
!          relation and/or selection rule for each atom.
!
!            NAT_SYM(J) is the number of prototypical atoms in the
!                       various symmetry sets :
!
!                             J = 1 : atom 0
!                             J = 2 : z axis
!                             J = 3 : x0y plane
!                             J = 4 : other atoms
!
      IF(IPRINT >= 2) WRITE(IUO1,889)                               !
      NAT_SYM(1) = 1                                                !
      NAT_SYM(2) = 0                                                !
      NAT_SYM(3) = 0                                                !
      NAT_SYM(4) = 0                                                !
      I_SET(1)   = 1                                                !
!
!   Loop on the prototypical atoms
!
      DO JTYP = 1, N_PROT                                           !
!
        ISTEP_L(JTYP)  = 1                                          !
        ISTEP_M(JTYP)  = 1                                          !
        I_REL_MP(JTYP) = 0                                          !
        I_LM(JTYP)     = 0                                          !
        I_Z_P(JTYP)    = 0                                          !
        Z_L_P(JTYP)    = ONE                                        !
        Z_M_P(JTYP)    = ONEC                                       !
!
        NSYM_P               = 1                                    !
        NSYM_PT              = 1                                    !
        INV_P(JTYP,NSYM_P)   = 1                                    !
        INV_PT(JTYP,NSYM_PT) = 1                                    !
!
        JAT_P = NCORR(1,JTYP)                                       !
        X     = COORD1(1,JAT_P)                                     !
        Y     = COORD1(2,JAT_P)                                     !
        Z     = COORD1(3,JAT_P)                                     !
        X_A   = ABS(X)                                              !
        Y_A   = ABS(Y)                                              !
        Z_A   = ABS(Z)                                              !
        IF(JTYP > 1) THEN                                           !
          IF((X_A+Y_A < SSMALL) .AND. (Z_A >= SSMALL)) THEN         !
            NAT_SYM(2)  = NAT_SYM(2) + 1                            !
            I_SET(JTYP) = 2                                         !
          ELSE IF((Z_A < SSMALL) .AND. (X_A+Y_A >= SSMALL)) THEN    !
            NAT_SYM(3)  = NAT_SYM(3) + 1                            !
            I_SET(JTYP) = 3                                         !
          ELSE IF(((X_A+Y_A) >= SSMALL) .AND. (Z_A >= SSMALL)) THEN !
            NAT_SYM(4)  = NAT_SYM(4) + 1                            !
            I_SET(JTYP) = 4                                         !
          END IF                                                    !
        END IF                                                      !
!
!     Loop on the symmetries keeping the cluster unchanged
!
        DO JSYM = 2 ,NSIZE_GR                                       !
!
          JROT = NSYM_G(JSYM)                                       !
          X1   = S_M(1,1,JROT) * X +                              & !
                 S_M(1,2,JROT) * Y +                              & !
                 S_M(1,3,JROT) * Z                                  !
          Y1   = S_M(2,1,JROT) * X +                              & !
                 S_M(2,2,JROT) * Y +                              & !
                 S_M(2,3,JROT) * Z                                  !
          Z1   = S_M(3,1,JROT) * X +                              & !
                 S_M(3,2,JROT) * Y +                              & !
                 S_M(3,3,JROT) * Z                                  !
          AD   = ABS(X - X1) + ABS(Y - Y1) + ABS(Z - Z1)            !
!
!       Case of an atom invariant by the symmetry JROT
!
          IF(AD < SSMALL) THEN                                      !
            NSYM_PT              = NSYM_PT + 1                      !
            INV_PT(JTYP,NSYM_PT) = JROT                             !
            IF(IZ(JROT) /= 0) THEN                                  !
              I_Z_P(JTYP)        = IZ(JROT)                         !
              NSYM_P             = NSYM_P + 1                       !
              INV_P(JTYP,NSYM_P) = JROT                             !
              ISL                = ISTEP_L(JTYP)                    !
              ISM                = ISTEP_M(JTYP)                    !
!
!       Case of an atom off the z axis
!
              IF((ABS(X) >= SSMALL) .OR. (ABS(Y) >= SSMALL)) THEN   !
!
!         Symmetry = IC2z
!
                IF(JROT == 36) THEN                                 !
                  ISTEP_M(JTYP) = MAX(ISM,2)                        !
                  I_LM(JTYP)    = 1                                 !
!
!         Symmetry = C2u or IC2u_perp
!
                ELSE                                                !
                  I_REL_MP(JTYP) = 1                                !
                  Z_L_P(JTYP)    = ZL(JROT)                         !
                  Z_M_P(JTYP)    = ZM1(JROT)                        !
                END IF                                              !
!
!       Case of an atom on the z axis but different from the absorber
!
              ELSE                                                  !
                IF(ABS(Z) >= SSMALL) THEN                           !
!
!         Symmetry = C2z
!
                  IF(JROT == 4) THEN                                !
                    ISTEP_M(JTYP) = MAX(ISM,2)                      !
!
!         Symmetry = C4z
!
                  ELSE IF(JROT == 15) THEN                          !
                    ISTEP_M(JTYP) = MAX(ISM,4)                      !
!
!         Symmetry = C4z3
!
                  ELSE IF(JROT == 18) THEN                          !
                    ISTEP_M(JTYP) = MAX(ISM,4)                      !
!
!         Symmetry = C3z
!
                  ELSE IF(JROT == 25) THEN                          !
                    ISTEP_M(JTYP) = MAX(ISM,3)                      !
!
!         Symmetry = C3z2
!
                  ELSE IF(JROT == 26) THEN                          !
                    ISTEP_M(JTYP) = MAX(ISM,3)                      !
!
!         Symmetry = C6z
!
                  ELSE IF(JROT == 27) THEN                          !
                    ISTEP_M(JTYP) = MAX(ISM,6)                      !
!
!         Symmetry = C6z5
!
                  ELSE IF(JROT == 28) THEN                          !
                    ISTEP_M(JTYP) = MAX(ISM,6)                      !
!
!         Symmetry = IC2u
!
                  ELSE IF(JROT > 33) THEN                           !
                    I_REL_MP(JTYP) = 1                              !
                    I_Z_P(JTYP)    = IZ(JROT)                       !
                    Z_L_P(JTYP)    = ZL(JROT)                       !
                    Z_M_P(JTYP)    = ZM1(JROT)                      !
                  END IF                                            !
!
!       Case of  atom 0 (the absorber)
!
                ELSE                                                !
!
!         Symmetry = C2z or IC2z
!
                  IF((JROT == 4) .OR. (JROT == 36)) THEN            !
                    ISTEP_M(JTYP) = MAX(ISM,2)                      !
                    IF(JROT == 36) THEN                             !
                      I_LM(JTYP) = 1                                !
                    END IF                                          !
!
!         Symmetry = C4z or IC4z
!
                  ELSE IF((JROT == 15) .OR. (JROT == 47)) THEN      !
                    ISTEP_M(JTYP) = MAX(ISM,4)                      !
                    IF(JROT == 47) THEN                             !
                      I_LM(JTYP) = 1                                !
                    END IF                                          !
!
!         Symmetry = C4z3 or IC4z3
!
                  ELSE IF((JROT == 18) .OR. (JROT == 50)) THEN      !
                    ISTEP_M(JTYP) = MAX(ISM,4)                      !
                    IF(JROT == 50) THEN                             !
                      I_LM(JTYP) = 1                                !
                    END IF
!
!         Symmetry = C3z or IC3z
!
                  ELSE IF((JROT == 25) .OR. (JROT == 57)) THEN      !
                    ISTEP_M(JTYP) = MAX(ISM,3)                      !
                    IF(JROT == 57) THEN                             !
                      ISTEP_L(JTYP) = MAX(ISL,2)                    !
                    END IF                                          !
!
!         Symmetry = C3z2 or IC3z2
!
                  ELSE IF((JROT == 26) .OR. (JROT == 58)) THEN      !
                    ISTEP_M(JTYP) = MAX(ISM,3)                      !
                    IF(JROT == 58) THEN                             !
                      ISTEP_L(JTYP) = MAX(ISL,2)                    !
                    END IF                                          !
!
!         Symmetry = C6z or IC6z
!
                  ELSE IF((JROT == 27) .OR. (JROT == 59)) THEN      !
                    ISTEP_M(JTYP) = MAX(ISM,6)                      !
                    IF(JROT == 59) THEN                             !
                      I_LM(JTYP) = 1                                !
                    END IF                                          !
!
!         Symmetry = C6z5 or IC6z5
!
                  ELSE IF((JROT == 28) .OR. (JROT == 60)) THEN      !
                    ISTEP_M(JTYP) = MAX(ISM,6)                      !
                    IF(JROT == 60) THEN                             !
                      I_LM(JTYP) = 1                                !
                    END IF                                          !
!
!         Symmetry = I
!
                  ELSE IF(JROT == 33) THEN                          !
                    ISTEP_L(JTYP) = MAX(ISL,2)                      !
!
!         Symmetry = C2u or IC2u_perp with u within (x0y)
!
                  ELSE                                              !
                    IF(IZ(JROT) == -1) THEN                         !
                      I_REL_MP(JTYP) = 1                            !
                      Z_L_P(JTYP)    = ZL(JROT)                     !
                      Z_M_P(JTYP)    = ZM1(JROT)                    !
                    END IF                                          !
                  END IF                                            !
                END IF                                              !
              END IF                                                !
            END IF                                                  !
          END IF                                                    !
        END DO                                                      !
!
!   Finding the symmetry group (if any) associated to each prototypical atom
!
        JGROUP = 0                                                  !
        IF(NSYM_PT == 1) THEN                                       !
          JGROUP = 1                                                !
        ELSE IF(NSYM_PT == 2) THEN                                  !
          IF(INV_PT(JTYP,2) == 33) THEN                             !
            JGROUP = 2                                              !
          ELSE IF(INV_PT(JTYP,2) == 4) THEN                         !
            JGROUP = 3                                              !
          ELSE IF(INV_PT(JTYP,2) == 36) THEN                        !
            JGROUP = 4                                              !
          END IF                                                    !
        ELSE IF(NSYM_PT == 3) THEN                                  !
          JGROUP = 16                                               !
        ELSE IF(NSYM_PT == 4) THEN                                  !
          IF(INV_PT(JTYP,3) == 33) THEN                             !
            JGROUP = 5                                              !
          ELSE IF( (INV_PT(JTYP,3) == 3) .OR.                     & !
                   (INV_PT(JTYP,3) == 19) ) THEN                    !
            JGROUP = 6                                              !
          ELSE IF( (INV_PT(JTYP,3) == 34) .OR.                    & !
                   (INV_PT(JTYP,3) == 51) ) THEN                    !
            JGROUP = 7                                              !
          ELSE IF(INV_PT(JTYP,3) == 15) THEN                        !
            JGROUP = 9                                              !
          ELSE IF(INV_PT(JTYP,3) == 47) THEN                        !
            JGROUP = 10                                             !
          END IF                                                    !
        ELSE IF(NSYM_PT == 6) THEN                                  !
          IF(INV_PT(JTYP,4) == 26) THEN                             !
            JGROUP = 18                                             !
          ELSE IF( (INV_PT(JTYP,4) == 34) .OR.                    & !
                   (INV_PT(JTYP,4) == 35) ) THEN                    !
            JGROUP = 19                                             !
          ELSE IF(INV_PT(JTYP,4) == 33) THEN                        !
            JGROUP = 17                                             !
          ELSE IF(INV_PT(JTYP,4) == 26) THEN                        !
            JGROUP = 21                                             !
          ELSE IF(INV_PT(JTYP,4) == 36) THEN                        !
            JGROUP = 22                                             !
          END IF                                                    !
        ELSE IF(NSYM_PT == 8) THEN                                  !
          IF(INV_PT(JTYP,4) == 33) THEN                             !
            IF(INV_PT(JTYP,8) == 50) THEN                           !
              JGROUP = 11                                           !
            ELSE                                                    !
              JGROUP = 8                                            !
            END IF                                                  !
          ELSE                                                      !
            IF(INV_PT(JTYP,5) == 15) THEN                           !
              JGROUP = 12                                           !
            ELSE                                                    !
              IF(INV_PT(JTYP,3) == 15) THEN                         !
                JGROUP = 13                                         !
              ELSE IF( (INV_PT(JTYP,3) == 3) .OR.                 & !
                       (INV_PT(JTYP,3) == 19) ) THEN                !
                JGROUP = 14                                         !
              END IF                                                !
            END IF                                                  !
           END IF                                                   !
        ELSE IF(NSYM_PT == 12) THEN                                 !
          IF(INV_PT(JTYP,5) == 5) THEN                              !
            JGROUP = 28                                             !
          ELSE IF(INV_PT(JTYP,7) == 33) THEN                        !
            IF(INV_PT(JTYP,12) == 60) THEN                          !
              JGROUP = 23                                           !
            ELSE                                                    !
              JGROUP = 20                                           !
            END IF                                                  !
          ELSE                                                      !
            IF(INV_PT(JTYP,3) == 3) THEN                            !
              JGROUP = 24                                           !
            ELSE IF(INV_PT(JTYP,9) == 59) THEN                      !
              JGROUP = 26                                           !
            ELSE IF(INV_PT(JTYP,5) == 27) THEN                      !
              JGROUP = 25                                           !
            END IF                                                  !
          END IF                                                    !
        ELSE IF(NSYM_PT == 16) THEN                                 !
          JGROUP = 15                                               !
        ELSE IF(NSYM_PT == 24) THEN                                 !
          IF(INV_PT(JTYP,17) == 57) THEN                            !
            JGROUP = 27                                             !
          ELSE IF(INV_PT(JTYP,17) == 17) THEN                       !
            JGROUP = 31                                             !
          ELSE IF(INV_PT(JTYP,17) == 49) THEN                       !
            JGROUP = 30                                             !
          ELSE IF(INV_PT(JTYP,17) == 37) THEN                       !
            JGROUP = 29                                             !
          END IF                                                    !
        ELSE IF(NSYM_PT == 48) THEN                                 !
          JGROUP = 32                                               !
        END IF                                                      !
        GR(JTYP) = JGROUP                                           !
!
        IF((IPRINT >= 2) .AND. (NSIZE_GR > 1)) THEN                 !
          WRITE(IUO1,851) JTYP,JAT_P,SYM_AT(1,JAT_P),             & !
                          SYM_AT(2,JAT_P),SYM_AT(3,JAT_P),        & !
                          I_SET(JTYP),NSYM_PT,NSYM_P,             & !
                          NAME_G(GR(JTYP)),                       & !
                          (NAME_S(INV_P(JTYP,JS)),JS=1,NSYM_P)      !
        END IF                                                      !
      END DO                                                        !
      WRITE(IUO1,852)                                               !
      GAIN_B = ZERO                                                 !
      DO JTYP = 1, N_PROT                                           !
        NGAIN_B(JTYP) = ISTEP_L(JTYP) * ISTEP_M(JTYP) *           & !
                        (I_REL_MP(JTYP) + 1)                        !
        GAIN_B        = GAIN_B + NGAIN_B(JTYP)                      !
        WRITE(IUO1,853) JTYP,I_Z_P(JTYP),INT(Z_L_P(JTYP)),        & !
                        Z_M_P(JTYP),ISTEP_L(JTYP),ISTEP_M(JTYP),  & !
                        I_LM(JTYP),I_REL_MP(JTYP),NGAIN_B(JTYP)     !
      END DO
      GAIN_B = GAIN_B / FLOAT(N_PROT)                               !
      WRITE(IUO1,855) GAIN_B                                        !
!
!  Calculation and storage of r^{l}_{m m'}(pi/2) for the specific
!            cubic group symmetry operations
!
      IF(I_CALC_ROT > 0) THEN                                       !
        PIS2 = PI * HALF                                            !
        CALL DJMN2(PIS2,R_PIS2,LI_M+1,0)                            ! adding the switch 0
        !CALL DJMN2(PIS2,R_PIS2,LI_M+1,2)
      END IF                                                        !
!
!  Construction of the inverse matrices used for point-groups
!
      DO I = 1, 3                                                   !
        DO J = 1, 3                                                 !
          DO JOP = 1, NSIZE_GR                                      !
            JSYM             = NSYM_G(JOP)                          !
            S_INV(JSYM,J,I) = S_M(I,J,JSYM)                         !
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
      GO TO 999                                                     !
!
 895  WRITE(IUO1,896) JAT1,JAT2                                     !
      STOP                                                          !
 998  WRITE(IUO1,997) NSIZE_GR                                      !
      STOP                                                          !
 999  PRINT 445                                                     !
!
!  Input/output formats
!
 125  FORMAT(2X,I4,5X,A2,5X,I2,3F10.4,12X,I4)
 444  FORMAT(////,5X,'++++++++++++++++++++  SYMMETRIZATION OF THE ',     &
             'CLUSTER   +++++++++++++++++++',/)
 445  FORMAT(//,5X,'+++++++++++++++++++++++++++++++++++++++++++++++',    &
                '+++++++++++++++++++++++++',//)
 557  FORMAT(10X,I3,3X,'(',F7.3,',',F7.3,',',F7.3,')',3X,I3,3X,I3,       &
             3X,A2,3X,F7.3,2X,A8)
 696  FORMAT(///,6X,I2,' SYMMETRY OPERATIONS AT LEAST SHOULD BE ',       &
             'PRESENT IN THE CLUSTER :',/)
 697  FORMAT(///,16X,' THE THREE-FOLD AXES HAVE BEEN CODED AS : ',//,    &
             30X,'alpha  ----->  l',/,                                   &
             30X,'beta   ----->  m',/,                                   &
             30X,'gamma  ----->  n',/,                                   &
             30X,'delta  ----->  o')
 698  FORMAT(///,17X,I1,' SYMMETRY OPERATION FOUND IN THE CLUSTER :',/)
 699  FORMAT(///,16X,I2,' SYMMETRY OPERATIONS FOUND IN THE CLUSTER :',/)
 700  FORMAT(12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,         &
             12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,         &
             12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,         &
             12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,         &
             12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,         &
             12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,         &
             12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,/,         &
             12X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5)
 701  FORMAT(19X,A5,2X,A5,2X,A5,2X,A5,2X,A5,2X,A5)
 702  FORMAT(26X,A5,2X,A5,2X,A5,2X,A5)
 703  FORMAT(30X,A5,2X,A5,2X,A5)
 704  FORMAT(33X,A5,2X,A5)
 705  FORMAT(37X,A5)
 706  FORMAT(/,19X,'AFTER ROTATION OF THIS CLUSTER BY ',A4,/)
 707  FORMAT(//,'    ----------',9X,'CONSTRUCTION OF A NEW CLUSTER :',   &
             8X,'----------')
 708  FORMAT(//,10X,'THE DIFFERENT GROUPS THAT COULD SUPPORT ',          &
             'THE CLUSTER ARE :',/)
 709  FORMAT(28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,               &
             28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,               &
             28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X))
 710  FORMAT(/,28X,'FOR THE INPUT CLUSTER ')
 711  FORMAT(/,28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,             &
               28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X),/,             &
               28X,4(A3,3X),/,28X,4(A3,3X),/,28X,4(A3,3X))
 712  FORMAT(/,24X,'FOR THE ROTATED INPUT CLUSTER')
 713  FORMAT(//,36X,'AND :')
 714  FORMAT(//,25X,'---> TRYING THE GROUP : ',A3)
 715  FORMAT(/,13X,'GROUP IMPOSSIBLE : ATOM (',                          &
             F7.3,',',F7.3,',',F7.3,')',/,19X,' TRANSFORMED INTO (',     &
             F7.3,',',F7.3,',',F7.3,')')
 716  FORMAT(20X,'THE NEW ATOM GENERATED BY ',A5,' IS ABOVE',/,20X,      &
             'THE SURFACE PLANE LOCATED AT Z = ',F7.3)
 717  FORMAT(27X,'NEW ATOM OF TYPE ',A2,' (No ',I4,') GENERATED BY ',A5, &
             /,27X,'AT THE POSITION OF AN ATOM OF TYPE ',A2,             &
             ' (No ',I4,')',/,27X,'IN THE ORIGINAL CLUSTER')
 718  FORMAT(///,5X,'<<<<<<<<<<  NATCLU_M TOO SSMALL IN THE .inc FILE ',  &
             'TO INCREASE  >>>>>>>>>>',/,5X,'<<<<<<<<<<        THE ',    &
             'CLUSTER. SHOULD BE AT LEAST ',I4,'       >>>>>>>>>>')
 719  FORMAT(///,3X,'<<<<<<<<<<  ERROR : NO GROUP WAS FOUND TO ',        &
             'ACCOMODATE THE CLUSTER  >>>>>>>>>>')
 720  FORMAT(20X,'THE NEW ATOM GENERATED BY ',A5,' IS INSIDE',/,25X,     &
             'THE ORIGINAL CLUSTER, BUT NOT',/,29X,                      &
             'AT AN ATOMIC POSITION')
 721  FORMAT(//,11X,' THE INVERSION I IS NOT CONSISTENT WITH THIS ',     &
            'CLUSTER',/,17X,'AS THE ABSORBER IS CLOSER TO THE SURFACE',  &
             /,21X,'THAN TO THE BOTTOM OF THE CLUSTER')
 722  FORMAT(//,16X,'THE CLUSTER CAN ACCOMODATE THE INVERSION I',//,     &
             13X,'---> BUILDING A LARGER CLUSTER INVARIANT BY I')
 723  FORMAT(//,4X,'--------   IMPOSSIBLE TO AUGMENT THE CLUSTER TO ',   &
             'SUPPORT I    --------',//,15X,                             &
             '---> RESTARTING WITH THE ORIGINAL CLUSTER ...')
 724  FORMAT(20X,'THE NEW ATOM GENERATED BY ',A5,' IS OF BULK TYPE',     &
             /,25X,'WHILE THE ORIGINAL ONE IS OF SURFACE TYPE',/,29X,    &
             '---> IMPOSSIBLE WITH I_GR = ',I1)
 725  FORMAT(20X,'THE NEW ATOM GENERATED BY ',A5,' IS OF SURFACE TYPE',  &
             /,25X,'WHILE THE ORIGINAL ONE IS OF BULK TYPE',/,29X,       &
             '---> IMPOSSIBLE WITH I_GR = ',I1)
 726  FORMAT(//,11X,' THE INVERSION I IS NOT CONSISTENT WITH THIS ',     &
            'CALCULATION',/,10X,'AS SURFACE AND BULK ATOMS HAVE TO ',    &
             'BE DISCRIMINATED (I_GR=1)')
 778  FORMAT(30X,I1)
 779  FORMAT(2X,I4,5X,A2,5X,I2,3F10.4,I5)
 849  FORMAT(1X,I4,4X,I2,3X,I4,2X,A5,' = ',I2,2X,F6.3,2X,                &
             '(',F6.3,',',F6.3,')',2X,'(',F6.3,',',F6.3,')',4X,I2)
 851  FORMAT(1X,I4,3X,I4,2X,'(',F7.3,',',F7.3,',',F7.3,')',2X,I1,        &
                          3X,I2,1X,I2,2X,A3,2X,                          &
                          4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X,  &
                          4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X,  &
                          4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X,  &
                          4(1X,A5),/,57X,4(1X,A5),/,57X,4(1X,A5),/,57X)
 852  FORMAT(///,10X,'SELECTION RULE AND/OR RELATION ON THE MATRIX',     &
                ' ELEMENTS OF TAU :',//,                                 &
                4X,'CLASS',2X,'I_Z',2X,'Z_L',8X,'Z_M',9X,'ISTEP_L',2X,   &
                'ISTEP_M',2X,'I_LM',2X,'I_REL_MP',2X,'GAIN',/)
 853  FORMAT(4X,I3,4X,I2,3X,I2,3X,'(',F6.3,',',F6.3,')',                 &
             6X,I1,8X,I1,7X,I1,7X,I1,6X,I2)
 854  FORMAT(//,16X,'-----> EXPECTED GAIN FOR THE GLOBAL LEVEL : ',      &
             F5.2,//)
 855  FORMAT(///,18X,'-----> EXPECTED GAIN FOR THE BASIS LEVEL : ',      &
             F5.2,//)
 880  FORMAT(33X,A3)
 881  FORMAT(/,'    ----------  THESE OPERATIONS DON''T FORM A',         &
               ' SYMMETRY GROUP  ----------', /                          &
               '    ----------  i. e. THE CLUSTER IS NOT CORRECTLY',     &
               ' TRUNCATED  ----------',/,                               &
               '    ----------       FROM THE POINT OF VIEW OF ',        &
               'SYMMETRY       ----------')
 886  FORMAT(///,25X,'CLUSTER SYMMETRY GROUP :',/,33X,A4)
 887  FORMAT(///,10X,'SYMMETRY OPERATIONS ASSOCIATED WITH THE ',         &
             'EQUIVALENT ATOMS : ',//,2X,'CLASS',2X,'ATOM',2X,           &
             ' No ',3X,'SYMMETRY',3X,'Z**L',                             &
             8X,'Z**M',13X,'Z**M''',7X,'DELTA',/)
 888  FORMAT(///,20X,'CONTENTS OF THE SYMMETRIZED CLUSTER',/,            &
                 18X,'BY INCREASING DISTANCE TO THE ABSORBER :',//,      &
                 11X,'No',13X,'(X,Y,Z)',11X,'CLASS',1X,'ATOM',           &
                 1X,'ChSp',4X,'DIST',//)
 889  FORMAT(///,10X,'SYMMETRY OPERATIONS LEAVING THE PROTOTYPICAL ',    &
                     'ATOMS INVARIANT : ',//,                            &
             27X,'G = 1  ----->  atom 0 only',/,                         &
             27X,'G = 2  ----->  atom along 0z',/,                       &
             27X,'G = 3  ----->  atom within x0y',/,                     &
             27X,'G = 4  ----->  other atom',//,                         &
             19X,'ST : total number of symmetries leaving',/,            &
             29X,'the prototypical atom invariant',//,                   &
             19X,'SE : number of symmetries taken into account',/,       &
             29X,'(Euler angle BETA = 0 or pi)',                         &
                     //,2X,'CLASS',2X,' No ',11X,'(X,Y,Z)',              &
                     10X,'G',3X,'ST',1X,'SE',1X,'GROUP',7X,              &
                     'SE SYMMETRIES',/)
 896  FORMAT(///,'<<<<<<<<<<  ERROR IN THE COORDINATES OF THE  ATOMS',   &
             '  >>>>>>>>>>',/,'<<<<<<<<<<  ATOMS ',I4,' AND ',I4,        &
             ' ARE IDENTICAL  >>>>>>>>>>')
 897  FORMAT(///,11X,'<<<<<<<<<<  NATP_M TOO SSMALL IN THE INCLUDE ',     &
             'FILE  >>>>>>>>>>',/,11X,'<<<<<<<<<<        SHOULD BE',     &
             ' AT LEAST ',I4,'         >>>>>>>>>>')
 898  FORMAT(///,12X,'<<<<<<<<<<  NAT_EQ_M TOO SSMALL IN THE INCLUDE',    &
             ' FILE  >>>>>>>>>>',/,12X,'<<<<<<<<<< SHOULD BE AT LEAST ', &
             I4,'  >>>>>>>>>>')
 997  FORMAT(///,'<<<<<<<<<<  ',I2,' SYMMETRIES HAVE BEEN FOUND. THIS ', &
             '>>>>>>>>>>',/,'<<<<<<<<<<  EXCEEDS THE SIZE OF THE ',      &
             'LARGEST POINT-GROUP  >>>>>>>>>>')
!
      END SUBROUTINE SYM_CLUS
!
END MODULE SYMMETRIZATION
