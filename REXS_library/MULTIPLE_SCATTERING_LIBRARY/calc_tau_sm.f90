!
!======================================================================
!
MODULE CALC_TAU0J_SM
!
!  This module provides the subroutine to compute
!    Tau^{0 j} or Tau^{j 0} in matrix series expansion approach
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE CALC_TAU_SM(TAU_TYPE,IATL,JE,TAU)
!
!  This subroutine computes the scattering path operator Tau^{0 0},
!    Tau^{0 J} or Tau^{j 0}
!
!  Input variables:
!
!     TAU_TYPE : scattering path operator type:
!                   'ABS' --> absorption case       Tau^{0 0}
!                   'INC' --> incoming beam case    Tau^{j 0}
!                   'OUT' --> outgoing beam case    Tau^{0 J}
!     IATL     : index of the atom for which BL is computed
!     JE       : energy point for which BL is computed
!     VAL      : plane indices array
!
!  Output variable:
!
!     TAU     : scattering path operator
!                   Note: Tau^{j 0}, Tau^{0 j} and Tau^{0 0} are stored
!                         in the same way to avoid having 3 TAU arrays
!
!  This is the A X = B version (partial inversion)
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 16 Jun 2021
!
!
      USE COMPLEX_NUMBERS
!
      USE CLUSTER_COORD
      USE CURRENT_BEAM,             ONLY : N_SCAT
      USE CURRENT_T_MATRIX
      USE INIT_L
!
      USE ANGULAR_MOMENTUM,         ONLY : GAUNT
      USE SPHERICAL_BESSEL
      USE SPHERICAL_HARMONICS
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)      ::  TAU_TYPE
!
      INTEGER, INTENT(IN)       ::  IATL,JE
!
      INTEGER                   ::  IBESS
      INTEGER                   ::  JLIN,JTYP,JNUM,JATL
      INTEGER                   ::  NBTYPJ,LMJ
      INTEGER                   ::  LJ,ILJ,MJ,INDJ
      INTEGER                   ::  KLIN,KTYP,KNUM,KATL
      INTEGER                   ::  NBTYPK,LMK
      INTEGER                   ::  LK,ILK,MK,INDK
      INTEGER                   ::  L_MIN,L_MAX
      INTEGER                   ::  L,M
      INTEGER                   ::  JS
      INTEGER                   ::  M_SIZE
!
      INTEGER                   ::  ABS
!
      REAL (WP)                 ::  PI4
      REAL (WP)                 ::  GNT(0:N_GAUNT)
      REAL (WP)                 ::  XJ,YJ,ZJ
      REAL (WP)                 ::  XKJ,YKJ,ZKJ,RKJ
      REAL (WP)                 ::  ZDKJ,KRKJ,ATTKJ
!
      REAL (WP)                 ::  SQRT,EXP
!
      COMPLEX (WP), INTENT(OUT) ::  TAU(LINMAX,LINFMAX,NATCLU_M)
!
      COMPLEX (WP)              ::  GOT(LINMAX*NATCLU_M,LINMAX*NATCLU_M)  ! GoT matrix
      COMPLEX (WP)              ::  MUL(LINMAX*NATCLU_M,LINMAX*NATCLU_M)  ! matrix GoT is multiplied by
      COMPLEX (WP)              ::  GMU(LINMAX*NATCLU_M,LINMAX*NATCLU_M)  ! result of matrix multiplication
      COMPLEX (WP)              ::  SSE(LINMAX*NATCLU_M,LINMAX*NATCLU_M)  ! sum of series expansion
      COMPLEX (WP)              ::  HL1(0:NLTWO)                          ! Hankel polynomials
      COMPLEX (WP)              ::  YLM(0:NLTWO,-NLTWO:NLTWO)             ! spherical harmonics
      COMPLEX (WP)              ::  TLJ,TLK,EXPKJ
      COMPLEX (WP)              ::  SUM_L
!
      DATA PI4 / 12.5663706143591729538505735331180E0_WP /          ! 4 pi
!
      IBESS = 3                                                     ! Hankel function of first kind
!
!  Construction of the GoT matrix: Elements are stored using
!    (i)  an angular momentum linear index LINJ representing (J,LJ)
!    (ii) a total linear index JLIN representing (J,LINJ)
!
!  Initialization of series expansion sum matrix SSE to I + GoT (single scattering result)
!  Initialization of MUL matrix to GoT
!
      JLIN = 0                                                      !
      DO JTYP = 1, N_PROT                                           !
        NBTYPJ = NATYP(JTYP)                                        !
        LMJ    = LMAX(JTYP,JE)                                      !
        DO JNUM = 1, NBTYPJ                                         !
          JATL = NCORR(JNUM,JTYP)                                   !
          XJ   = SYM_AT(1,JATL)                                     !
          YJ   = SYM_AT(2,JATL)                                     !
          ZJ   = SYM_AT(3,JATL)                                     !
!
          DO LJ = 0, LMJ                                            !
            ILJ = LJ * LJ + LJ + 1                                  !
            TLJ = TL(LJ,1,JTYP,JE)                                  !
            DO MJ = -LJ, LJ                                         !
              INDJ = ILJ + MJ                                       !
              JLIN = JLIN + 1                                       !
!
              KLIN = 0                                              !
              DO KTYP = 1, N_PROT                                   !
                NBTYPK = NATYP(KTYP)                                !
                LMK    = LMAX(KTYP,JE)                              !
                DO KNUM = 1, NBTYPK                                 !
                  KATL = NCORR(KNUM,KTYP)                           !
!
                  IF(KATL /= JATL) THEN                             !
                    XKJ = SYM_AT(1,KATL) - XJ                       !
                    YKJ = SYM_AT(2,KATL) - YJ                       !
                    ZKJ = SYM_AT(3,KATL) - ZJ                       !
                    RKJ = SQRT(XKJ * XKJ + YKJ * YKJ + ZKJ * ZKJ)   !
!
                    KRKJ  = VK(JE) * RKJ                            !
                    ATTKJ = EXP(- AIMAG(VK(JE)) * RKJ)              !
                    EXPKJ = (XKJ + IC * YKJ) / RKJ                  !
                    ZDKJ  = ZKJ / RKJ                               !
!
                    CALL SPH_HAR(2*NL_M,ZDKJ,EXPKJ,YLM,LMJ+LMK)     !
                    CALL BESPHE(LMJ+LMK,IBESS,KRKJ,HL1)             !
!
                  END IF                                            !
!
                  DO LK = 0, LMK                                    !
                    ILK   = LK * LK + LK + 1                        !
                    L_MIN = ABS(LK - LJ)                            !
                    L_MAX = LK + LJ                                 !
                    TLK   = TL(LK,1,KTYP,JE)                        !
                    DO MK = - LK ,LK
                      INDK          = ILK + MK                      !
                      KLIN          = KLIN + 1                      !
!
                      SSE(KLIN,JLIN) = ZEROC                        !
                      MUL(KLIN,JLIN) = ZEROC                        ! initialization
                      GOT(KLIN,JLIN) = ZEROC                        ! of the matrices
                      GMU(KLIN,JLIN) = ZEROC                        !
!
                      SUM_L         = ZEROC                         !
                      IF(KATL /= JATL) THEN                         !
                        CALL GAUNT(LK,MK,LJ,MJ,GNT)                 !
!
                        DO L = L_MIN , L_MAX, 2                     !
                          M = MJ - MK                               ! calculation of
                          IF(ABS(M) <= L) THEN                      ! the matrix elements
                            SUM_L = SUM_L + (IC**L) * HL1(L) *    & ! of the free
                                            YLM(L,M) * GNT(L)       ! electron propagator
                          END IF                                    !
                        END DO                                      ! G_{LK LJ}^{K J}
                        SUM_L = SUM_L * ATTKJ * PI4 * IC            !
                      ELSE                                          !
                        SUM_L = ZEROC                               !
                      END IF                                        !
!
                      IF(KLIN == JLIN) THEN                         !
                         SSE(KLIN,JLIN) = ONEC + TLK * SUM_L        ! initialization of
                      ELSE                                          ! SSE to I + GoT
                         SSE(KLIN,JLIN) = TLK * SUM_L               ! (single scattering)
                      END IF                                        !
!
                      MUL(KLIN,JLIN)=TLK*SUM_L                      ! GoT + initialization
                      GOT(KLIN,JLIN)=TLK*SUM_L                      ! of to MUL to GoT
!
                    END DO                                          !
                  END DO                                            !
!
                END DO                                              !
              END DO                                                !
!
            END DO                                                  !
          END DO                                                    !
!
        END DO                                                      !
      END DO                                                        !
!                                                                   ! single scattering already
      IF(N_SCAT > 1) THEN                                           ! contained in SSE
!
!  Size of the square matrices to be multiplied
!
        M_SIZE = JLIN                                               !
!
!   Loop over the scattering orders starting with 2 (order 1 = GoT)
!
        DO JS = 2, N_SCAT                                           !
!
!  Matrix multiplication of MUL by GoT: FORTRAN intrinsic MATMUL version
!
          GMU = MATMUL(MUL,GOT)                                     !
!
!  Storage of series expansion in SSE and MUL set to GMU
!
          DO JLIN = 1, M_SIZE                                       !
            DO KLIN = 1 ,M_SIZE                                     !
              SSE(KLIN,JLIN) = SSE(KLIN,JLIN) + GMU(KLIN,JLIN)      !
              MUL(KLIN,JLIN) = GMU(KLIN,JLIN)                       !
              GMU(KLIN,JLIN) = ZEROC                                !
            END DO                                                  !
          END DO                                                    !
!
        END DO                                                      !
!
      END IF                                                        !
!
!   Storage of the Tau matrix (first block column of SSE)
!
      JLIN = 0                                                      !
      DO JTYP = 1, N_PROT                                           !
        NBTYPJ = NATYP(JTYP)                                        !
        LMJ = LMAX(JTYP,JE)                                         !
        DO JNUM = 1, NBTYPJ                                         !
          JATL = NCORR(JNUM,JTYP)                                   !
!
          DO LJ = 0, LMJ                                            !
            ILJ = LJ * LJ + LJ+  1                                  !
            TLJ = TL(LJ,1,JTYP,JE)                                  !
            DO MJ = - LJ, LJ                                        !
              INDJ = ILJ + MJ                                       !
              JLIN = JLIN + 1                                       !
!
              KLIN = 0                                              !
              DO KTYP = 1, N_PROT                                   !
                NBTYPK = NATYP(KTYP)                                !
                LMK    = LMAX(KTYP,JE)                              !
                DO KNUM = 1, NBTYPK                                 !
                  KATL = NCORR(KNUM,KTYP)                           !
!
                  DO LK = 0, LMK                                    !
                    ILK = LK * LK + LK + 1                          !
                    DO MK = - LK, LK                                !
                      INDK = ILK + MK                               !
                      KLIN = KLIN + 1                               !
                      IF((JATL == 1) .AND. (LJ == LF2)) THEN        !
                        TAU(INDK,INDJ,KATL) = SSE(KLIN,JLIN) * TLJ  !
                      END IF                                        !
                    END DO                                          !
                  END DO                                            !
!
                END DO                                              !
              END DO                                                !
!
            END DO                                                  !
          END DO                                                    !
!
        END DO                                                      !
      END DO                                                        !
!
      END SUBROUTINE CALC_TAU_SM
!
END MODULE CALC_TAU0J_SM

