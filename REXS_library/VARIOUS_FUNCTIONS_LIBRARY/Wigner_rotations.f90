!
!=======================================================================
!
MODULE WIGNER_ROTATIONS
!
!  This module contains the Wigner's rotations subroutine
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NL_M
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE DJMN2(RBETA,R,LMAX,ISWITCH)
!
!  This routine calculates Wigner rotation matrices R^{L}_{M1 M2} up to
!             order LMAX, following Messiah's convention.
!  They are stored as R(M2,M1,L) and multiplied (ISWITCH = 1) or divided
!             by EXPF.
!
!
!  Input variables:
!
!         RBETA       :  Euler angle beta
!         LMAX        :  l maximum
!         ISWITCH     :  switch
!                            ISWITCH = 0 : R(M1,M2,L)
!                            ISWITCH = 1 : (2l+1) * R(M1,M2,L) * EXPF(ABS(M2),L)
!                            ISWITCH = 2 : (2l+1) * R(M1,M2,L) / EXPF(ABS(M2),L)
!
!  Output variable:
!
!         FSPH        :  R(M1,M2,L)
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 19 May 2021
!
!
      USE STORE_COEF
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE,HALF
      USE SQUARE_ROOTS,      ONLY : SQR2
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  LMAX,ISWITCH
!
      INTEGER                  ::  EPS0
      INTEGER                  ::  L,M,M1,M2,L1
!
      INTEGER                  ::  MOD
!
      REAL (WP), INTENT(IN)    ::  RBETA
      REAL (WP), INTENT(OUT)   ::  R(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
!
      REAL (WP)                ::  SSMALL
      REAL (WP)                ::  C,S,CC,CMUL,S1,C1
      REAL (WP)                ::  PRODL,COEF,CL,FLL1,FLL2
      REAL (WP)                ::  PRODM,CM,FLLM
      REAL (WP)                ::  R_1,R_2,R_A,R_B
      REAL (WP)                ::  CM2,CF1,CF2,D0
!
      REAL (WP)                ::  COS,SIN,ABS,FLOAT,SQRT
!
      DATA SSMALL / 0.001E0_WP /
!
      C    = COS(RBETA) * HALF                                      !
      S    = SIN(RBETA) * HALF                                      !
      CC   = C + C                                                  !
      CMUL = - ONE                                                  !
!
      IF(ABS(S) < SSMALL) THEN                                      !
!
        IF(C > ZERO) EPS0 =   ONE                                   !
        IF(C < ZERO) EPS0 = - ONE                                   !
        DO L = 0, LMAX                                              !
          DO M1 = - L, L                                            !
            DO M2= - L, L                                           !
              IF(M1 /= M2 * EPS0) THEN                              !
                R(M2,M1,L) = ZERO                                   !
              ELSE                                                  !
                IF(EPS0 == 1) THEN                                  !
                  R(M2,M1,L) = ONE                                  !
                ELSE                                                !
                  IF(MOD(L + M1,2) == 0) THEN                       !
                    R(M2,M1,L) =   ONE                              !
                  ELSE                                              !
                    R(M2,M1,L) = - ONE                              !
                  END IF                                            !
                END IF                                              !
              END IF                                                !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
!
      ELSE                                                          !
!
        S1 = S * SQR2                                               !
        C1 = HALF + C                                               !
!
        R(0,0,0)   = ONE                                            !
        R(-1,-1,1) = C1                                             !
        R(0,-1,1)  = S1                                             !
        R(1,-1,1)  = ONE - C1                                       !
        R(-1,0,1)  = - S1                                           !
        R(0,0,1)   = CC                                             !
        R(1,0,1)   = S1                                             !
        R(-1,1,1)  = ONE - C1                                       !
        R(0,1,1)   = - S1                                           !
        R(1,1,1)   = C1                                             !
!
        PRODL = - S                                                 !
        COEF  = - S / C1                                            !
        CL    = - ONE                                               !
!
        DO L = 2, LMAX                                              !
          CL    = - CL                                              !
          L1    = L - 1                                             !
          FLL1  = CC * FLOAT(L + L1)                                !
          FLL2  = ONE / (FLOAT(L * L1) * CC)                        !
          PRODL = - PRODL * S                                       !
!
!  Case M = 0
!
          R_1 = EXPR(0,L) * PRODL                                   !
!
          R(-L,0,L)  = R_1                                          !
!
          R(L,0,L)   = R_1 * CL                                     !
          R(0,-L,L ) = R_1 * CL                                     !
!
          R(0,L,L)   = R_1                                          !
!
          CM2 = CL                                                  !
!
          DO M2 = - L1, -1                                          !
            CM2 = CM2 * CMUL                                        !
            CF1 = CF(L1,0,-M2) / FLL1                               !
            CF2 = FLL1 / CF(L,0,-M2)                                !
            IF(- M2 < L1) THEN                                      !
              R_A = CF2 * (R(M2,0,L1) - R(M2,0,L-2) * CF1)          !
            ELSE                                                    !
              R_A = CF2 * R(M2,0,L1)                                !
            END IF                                                  !
!
            R(M2,0,L)  = R_A                                        !
            R(-M2,0,L) = R_A * CM2                                  !
            R(0,M2,L)  = R_A * CM2                                  !
            R(0,-M2,L) = R_A                                        !
!
          END DO                                                    !
!
          R(0,0,L) = FLL1 * R(0,0,L1) / CF(L,0,0) -               & !
                     R(0,0,L-2) * CF(L1,0,0) / CF(L,0,0)            !
!
!  Case M > 0
!
          PRODM = ONE                                               !
          CM    = CL                                                !
          FLLM  = ZERO                                              !
!
          DO M = 1, L1                                              !
            CM    = - CM                                            !
            PRODM = PRODM * COEF                                    !
            FLLM  = FLLM + FLL2                                     !
!
            R_1 = EXPR(M,L) * PRODL * PRODM                         !
            R_2 = R_1 / (PRODM * PRODM)                             !
!
            R(-L,M,L)  = R_1                                        !
            R(-L,-M,L) = R_2                                        !
!
            R(L,-M,L)  = R_1 * CM                                   !
            R(M,-L,L)  = R_1 * CM                                   !
            R(L,M,L)   = R_2 * CM                                   !
            R(-M,-L,L) = R_2 * CM                                   !
!
            R(-M,L,L)  = R_1                                        !
            R(M,L,L)   = R_2                                        !
!
            CM2 = CM                                                !
!
            DO M2 = - L1, - M                                       !
              CM2 = - CM2                                           !
              D0  = FLOAT(M2) * FLLM                                !
              CF1 = CF(L1,M,-M2) / FLL1                             !
              CF2 = FLL1 / CF(L,M,-M2)                              !
              IF((M < L1) .AND. (-M2 < L1)) THEN                    !
                R_A = CF2 * ( (ONE - D0) * R(M2,M,L1)  -          & !
                               R(M2,M,L-2)  * CF1 )                 !
                R_B = CF2 * ( (ONE + D0) * R(M2,-M,L1) -          & !
                               R(M2,-M,L-2) * CF1 )                 !
              ELSE                                                  !
                R_A = CF2 * (ONE - D0) * R(M2,M,L1)                 !
                R_B = CF2 * (ONE + D0) * R(M2,-M,L1)                !
              END IF                                                !
!
              R(M2,M,L)   = R_A                                     !
              R(M2,-M,L)  = R_B                                     !
!
              R(-M2,-M,L) = R_A * CM2                               !
              R(M,M2,L)   = R_A * CM2                               !
              R(-M,M2,L)  = R_B * CM2                               !
              R(-M2,M,L)  = R_B * CM2                               !
!
              R(-M,-M2,L) = R_A                                     !
              R(M,-M2,L)  = R_B                                     !
!
            END DO                                                  !
          END DO                                                    !
!
          PRODM = PRODM * COEF                                      !
!
          R_1 = PRODL * PRODM                                       !
          R_2 = PRODL / PRODM                                       !
!
          R(-L,L,L)  = R_1                                          !
          R(L,-L,L)  = R_1                                          !
          R(L,L,L)   = R_2                                          !
          R(-L,-L,L) = R_2                                          !
!
        END DO                                                      !
!
      END IF                                                        !
!
      IF(ISWITCH == 1) THEN                                         !
        DO L = 0, LMAX                                              !
          DO M1 = - L, L                                            !
             DO M2 = - L, L                                         !
               R(M2,M1,L) = SQRT(FLOAT(L + L + 1)) *              & !
                            R(M2,M1,L) * EXPF(ABS(M2),L)            !
             END DO                                                 !
          END DO                                                    !
        END DO                                                      !
      ELSE IF(ISWITCH == 2) THEN                                    !
        DO L = 0, LMAX                                              !
          DO M1 = - L, L                                            !
             DO M2 = - L, L                                         !
               R(M2,M1,L) = SQRT(FLOAT(L + L + 1)) *              & !
                            R(M2,M1,L) / EXPF(ABS(M2),L)            !
             END DO                                                 !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
      END SUBROUTINE DJMN2
!
END MODULE WIGNER_ROTATIONS
