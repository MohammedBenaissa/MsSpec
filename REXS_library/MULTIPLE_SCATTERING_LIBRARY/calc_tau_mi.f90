!
!======================================================================
!
MODULE CALC_TAU0J_MI
!
!  This module provides the subroutine to compute
!    Tau^{0 j} or Tau^{j 0} in the partial inversion algorithm (A X = B version)
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE CALC_TAU_MI(TAU_TYPE,IATL,JE,TAU)
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
      CHARACTER (LEN =  1)      ::  CH
!
      INTEGER, INTENT(IN)       ::  IATL,JE
!
      INTEGER                   ::  LOGF
      INTEGER                   ::  IPIV(LINMAX*NATCLU_M)
      INTEGER                   ::  IBESS
      INTEGER                   ::  JLIN,JTYP,JNUM,JATL
      INTEGER                   ::  NBTYPJ,LMJ
      INTEGER                   ::  LJ,ILJ,MJ,INDJ
      INTEGER                   ::  KLIN,KTYP,KNUM,KATL
      INTEGER                   ::  NBTYPK,LMK
      INTEGER                   ::  LK,ILK,MK,INDK
      INTEGER                   ::  L_MIN,L_MAX
      INTEGER                   ::  L,M
      INTEGER                   ::  LW2
      INTEGER                   ::  INFO,INFO1
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
      COMPLEX (WP)              ::  HL1(0:NLTWO)
      COMPLEX (WP)              ::  SM(LINMAX*NATCLU_M,LINMAX*NATCLU_M)
      COMPLEX (WP)              ::  IN(LINMAX*NATCLU_M,LINMAX)
      COMPLEX (WP)              ::  YLM(0:NLTWO,-NLTWO:NLTWO)
      COMPLEX (WP)              ::  TLJ,TLK,EXPKJ
      COMPLEX (WP)              ::  SUM_L
!
      DATA PI4 / 12.5663706143591729538505735331180E0_WP /          ! 4 pi
!
      LOGF = 6                                                      ! log file unit
!
      IBESS = 3                                                     ! Hankel function of first kind
      CH    = 'N'                                                   !
!
!  Construction of the multiple scattering matrix MS = (I-GoT).
!    Elements are stored using a linear index LINJ representing
!    (J,LJ)
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
                    CALL SPH_HAR(2*NL_M,ZDKJ,EXPKJ,YLM,LMJ+LMK)     !
                    CALL BESPHE(LMJ+LMK,IBESS,KRKJ,HL1)             !
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
                      SM(KLIN,JLIN) = ZEROC                         !
                      SUM_L         = ZEROC                         !
!
                      IF(KATL /= JATL) THEN                         !
                        CALL GAUNT(LK,MK,LJ,MJ,GNT)                 !
!
                        DO L = L_MIN , L_MAX, 2                     !
                          M = MJ - MK                               !
                          IF(ABS(M) <= L) THEN                      !
                            SUM_L = SUM_L + (IC**L) * HL1(L) *    & !
                                            YLM(L,M) * GNT(L)       !
                          END IF                                    !
                        END DO                                      !
                        SUM_L = SUM_L * ATTKJ * PI4 * IC            !
                      ELSE                                          !
                        SUM_L = ZEROC                               !
                      END IF                                        !
!
                      IF(KLIN == JLIN) THEN                         !
                         SM(KLIN,JLIN) = ONEC - TLK * SUM_L         !
                         IF(JTYP.EQ.1) THEN                         !
                            IN(KLIN,JLIN) = ONEC                    !
                         END IF                                     !
                      ELSE                                          !
                         SM(KLIN,JLIN) = - TLK * SUM_L              !
                         IF(JTYP == 1) THEN                         !
                            IN(KLIN,JLIN) = ZEROC                   !
                         END IF                                     !
                      END IF                                        !
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
!
      LW2 = (LMAX(1,JE) + 1 ) * (LMAX(1,JE) + 1)                    !
!
!   Partial inversion of the multiple scattering matrix MS and
!     multiplication by T : the LAPACK subroutine performing
!
!                          A * X = B
!
!     is used where B is the block column corresponding to
!     the absorber 0 in the identity matrix. X is then TAU^{j 0}.
!
!  1) LU factorization of matrix A = SM (scattering matrix)
!
      CALL ZGETRF(JLIN,JLIN,SM,LINMAX*NATCLU_M,IPIV,INFO1)          !
!
      IF(INFO1 /= 0) THEN                                           !
         WRITE(LOGF,10) INFO1                                       ! error in the LU factorization
      ELSE                                                          !
!
!  2) Solution of A * X = B for A LU-factorized:
!
!       IN : right-hand side of matrix B on entry
!            matrix X on exit i.e. TAU^{j 0}
!
        CALL ZGETRS(CH,JLIN,LW2,SM,LINMAX*NATCLU_M,IPIV,          & !
                    IN,LINMAX*NATCLU_M,INFO)                        !
        IF(INFO /= 0) THEN                                          !
          WRITE(LOGF,20) INFO                                       ! error in the inversion
        END IF                                                      !
!
      END IF                                                        !
!
!   Storage of the Tau matrix
!
      JLIN = 0                                                      !
      DO JTYP = 1, N_PROT                                           !
        NBTYPJ = NATYP(JTYP)                                        !
        LMJ    = LMAX(JTYP,JE)                                      !
        DO JNUM = 1, NBTYPJ                                         !
          JATL = NCORR(JNUM,JTYP)                                   !
!
          DO LJ = 0, LMJ                                            !
            ILJ = LJ * LJ + LJ + 1                                  !
            TLJ = TL(LJ,1,JTYP,JE)                                  !
            DO MJ = - LJ, LJ                                        !
              INDJ = ILJ + MJ                                       !
              JLIN = JLIN + 1                                       !
!
              KLIN = 0                                              !
              DO KTYP=1,N_PROT                                      !
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
                        TAU(INDK,INDJ,KATL) = IN(KLIN,JLIN) * TLJ   !
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
!  Formats:
!
  10  FORMAT(//,5X,'    --->  INFO  = ',I2)
  20  FORMAT(//,5X,'    --->  INFO1 = ',I2)
!
      END SUBROUTINE CALC_TAU_MI
!
END MODULE CALC_TAU0J_MI
