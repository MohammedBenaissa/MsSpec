!
!=======================================================================
!
MODULE SPHERICAL_HARMONICS
!
!  This module contains routines to compute the spherical harmonics
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE HARSPH0(NL,THETA,PHI,YLM,NC)
!
!  This routine computes the complex spherical harmonics using Condon and
!                  Shortley phase convention.
!
!  Input variables:
!
!     NL        : dimension of YLM array
!     THETA     : polar angle     in radians
!     PHI       : azimuthal angle in radians
!     NC        : maximum l value computed
!
!  Output variable:
!
!     YLM       : spherical harmonics array Y^l_m(theta,phi)
!                   stored as YLM(L,M)
!
!
!   Author :  D. Sébilleau
!
!
!                                           Last modified : 26 May 2021
!
!
      USE STORE_COEF,           ONLY : EXPF2,FSQ
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE,TWO,FOUR,HALF
      USE COMPLEX_NUMBERS
      USE PI_ETC,            ONLY : PI
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  NL,NC
!
      INTEGER                  ::  L,M
      INTEGER                  ::  ICOS
!
      REAL (WP), INTENT(IN)    ::  THETA,PHI
      REAL (WP)                ::  SQ4PI_INV,SQR3_INV,SSMALL
      REAL (WP)                ::  X,Y,C2,C1,SL,M_ONE
!
      REAL (WP)                ::  COS,ABS,FLOAT,SQRT
!
      COMPLEX(WP), INTENT(OUT) ::  YLM(0:NL,-NL:NL)
!
      COMPLEX(WP)              ::  COEF,YMM,YMMP,C
!
      DATA SQ4PI_INV  / 0.28209479177387814347403972578E0_WP /      ! sqrt(1 / 4*pi)
      DATA SQR3_INV   / 0.48860251190291992158638462284E0_WP /      ! sqrt(3 / 4*pi)
      DATA SSMALL     / 0.00010E0_WP                         /
!
      COMPLEX(WP)              ::  CMPLX,CONJG,EXP
!
!  Checking for pathological values of cos(theta)
!
      X = COS(THETA)                                                !
!
      ICOS = 2                                                      !
!
      IF(ABS(X) < SSMALL) THEN                                      ! theta = pi/2
        X     = ZERO                                                !
        ICOS  = 0                                                   !
      END IF                                                        !
      IF(ABS(X + ONE) < SSMALL) THEN                                ! theta = pi
        X     = - ONE                                               !
        ICOS  = - 1                                                 !
        SL    = - ONE                                               !
        M_ONE = - ONE                                               !
      END IF                                                        !
      IF(ABS(X - ONE) < SSMALL) THEN                                ! theta = 0
        X     = ONE                                                 !
        ICOS  = 1                                                   !
        SL    = ONE                                                 !
        M_ONE = ONE                                                 !
      END IF                                                        !
!
!  Initialization of YLM
!
      DO L = 0, NC                                                  !
        DO M = -L, L                                                !
          YLM(L,M) = ZEROC                                          !
        END DO                                                      !
      END DO                                                        !
!
!  Particular cases: theta = 0, pi --> X = +/- 1
!
      IF(ABS(ICOS) == 1) THEN                                       !
        DO L = 0, NC                                                !
          SL       = SL * M_ONE                                     !
          YLM(L,0) = SQRT(TWO * FLOAT(L) + ONE) * SQ4PI_INV * SL    !
        END DO                                                      !
        GO TO 10                                                    !
      END IF                                                        !
!
!  YLM for M = 0:
!
      YLM(0,0) = CMPLX(SQ4PI_INV)                                   !
      YLM(1,0) = CMPLX(X * SQR3_INV)                                !
!
      DO L = 2, NC
        Y        = ONE / FLOAT(L)
        YLM(L,0) = X * SQRT(FOUR - Y * Y) * YLM(L-1,0) -          & !
                   (ONE - Y) * SQRT( ONE + TWO /                  & !
                                     (FLOAT(L) - 1.5E0_WP)        & !
                                   ) * YLM(L-2,0)                   !
      END DO                                                        !
!
      C2 = - ONE                                                    !
      IF((THETA >= ZERO) .AND. (THETA <= PI)) THEN                  !
        C = - HALF * SQRT(ONE - X * X) * EXP( IC * PHI)             !
      ELSE                                                          !
        C =   HALF * SQRT(ONE - X * X) * EXP( IC * PHI)             !
      END IF                                                        !
!
!  YLM for |M| > 0
!
      C1   = ONE                                                    !
      COEF = ONEC                                                   !
!
      DO M = 1, NC                                                  !
        C1          = C1 * C2                                       !
        COEF        = COEF * C                                      !
        YMM         = SQ4PI_INV * COEF * FSQ(M)                     !
        YLM(M,M)    = YMM                                           !
        YLM(M,-M)   = C1 * CONJG(YMM)                               !
        YMMP        = X * SQRT(FLOAT(M + M + 3)) * YMM              !
        YLM(M+1,M)  = YMMP                                          !
        YLM(M+1,-M) = C1*  CONJG(YMMP)                              !
        IF(M < NC-1) THEN                                           !
          DO L = M + 2, NC                                          !
            YLM(L,M)  = ( X * FLOAT(L + L - 1) * EXPF2(L-1,M) *   & !
                          YLM(L-1,M) - FLOAT(L+M-1) *             & !
                          EXPF2(L-2,M) * YLM(L-2,M)               & !
                        ) / ( EXPF2(L,M) * FLOAT(L-M) )             !
            YLM(L,-M) = C1 * CONJG(YLM(L,M))                        !
          END DO                                                    !
        END IF                                                      !
      END DO                                                        !
!
  10  RETURN                                                        !
!
      END SUBROUTINE HARSPH0

!
!=======================================================================
!
      SUBROUTINE HARSPH(NL,THETA,PHI,YLM,NC)
!
!  This routine computes the complex spherical harmonics using Condon and
!                  Shortley phase convention.
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 31 May 2021
!
!
!
      USE STORE_COEF,           ONLY : EXPF2,FSQ
!
      USE REAL_NUMBERS,         ONLY : ZERO,ONE,TWO,FOUR,HALF
      USE COMPLEX_NUMBERS
      USE PI_ETC,               ONLY : PI
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  NL,NC
!
      INTEGER                  ::  L,M
      INTEGER                  ::  ICOS
!
      REAL (WP), INTENT(IN)    ::  THETA,PHI
      REAL (WP)                ::  SQ4PI_INV,SQR3_INV,SSMALL
      REAL (WP)                ::  X,Y,C2,C1
      REAL (WP)                ::  SL,M_ONE
!
      REAL (WP)                ::  COS,ABS,FLOAT,SQRT
!
      COMPLEX(WP), INTENT(OUT) ::  YLM(0:NL,-NL:NL)
!
      COMPLEX(WP)              ::  COEF,YMM,YMMP,C
!
      COMPLEX(WP)              ::  CMPLX,CONJG,EXP
!
      DATA SQ4PI_INV  / 0.28209479177387814347403972578E0_WP /      ! sqrt(1 / 4*pi)
      DATA SQR3_INV   / 0.48860251190291992158638462284E0_WP /      ! sqrt(3 / 4*pi)
      DATA SSMALL     / 0.00010E0_WP                         /
!
      X = COS(THETA)                                                !
!
!  Checking for pathological values of cos(theta)
!
      IF(ABS(X) < SSMALL)       THEN                                !
        X     = ZERO                                                !
        ICOS  = 0                                                   !
      END IF                                                        !
      IF(ABS(X + ONE) < SSMALL) THEN                                !
        X     = - ONE                                               !
        ICOS  = - 1                                                 !
        SL    = - ONE                                               !
        M_ONE = - ONE                                               !
      END IF                                                        !
      IF(ABS(X - ONE) < SSMALL) THEN                                !
        X     =   ONE                                               !
        ICOS  =   1                                                 !
        SL    =   ONE                                               !
        M_ONE =   ONE                                               !
      END IF                                                        !
!
!  Initialization of YLM
!
      DO L = 0, NC                                                  !
        DO M = -L, L                                                !
          YLM(L,M) = ZEROC                                          !
        END DO                                                      !
      END DO                                                        !
!
!  Particular cases: theta = 0, pi --> X = +/- 1
!
      IF(ABS(ICOS) == 1) THEN                                       !
        DO L = 0, NC                                                !
          SL       = 1 !SL * M_ONE                                     !
          YLM(L,0) = SQRT( TWO * FLOAT(L) + ONE ) * SQ4PI_INV * SL  !
        END DO                                                      !
        GO TO 10                                                    !
      END IF                                                        !
!
!  YLM for M = 0:
!
      YLM(0,0) = CMPLX(SQ4PI_INV)                                   !
      YLM(1,0) = X * SQR3_INV                                       !
!
      DO L = 2, NC
        Y        = ONE / FLOAT(L)                                   !
        YLM(L,0) = X * SQRT(FOUR - Y * Y) * YLM(L-1,0) -          & !
                  (ONE - Y) * SQRT( ONE + TWO / (FLOAT(L) -       & !
                                    1.5E0_WP) ) * YLM(L-2,0)        !
      END DO                                                        !
!
      C2 = - ONE                                                    !
      IF((THETA >= ZERO) .AND. (THETA <= PI)) THEN                  !
        C = - HALF * SQRT(ONE - X * X) * EXP(IC * PHI)              !
      ELSE                                                          !
        C =   HALF * SQRT(ONE - X * X) * EXP(IC * PHI)              !
      END IF                                                        !
!
      C1   = ONE                                                    !
      COEF = ONEC                                                   !
!
!  YLM for |M| > 0
!
      DO M = 1, NC                                                  !
        C1          = C1 * C2                                       !
        COEF        = COEF * C                                      !
        YMM         = SQ4PI_INV * COEF * FSQ(M)                     !
        YMMP        = X * SQRT( FLOAT(M + M + 3) ) * YMM            !
        YLM(M, M)   = YMM                                           !
        YLM(M,-M)   = C1 * CONJG(YMM)                               !
        YLM(M+1, M) = YMMP                                          !
        YLM(M+1,-M) = C1 * CONJG(YMMP)                              !
        IF(M < NC-1) THEN                                           !
          DO L = M+2, NC                                            !
            YLM(L,M)  = ( X * FLOAT(L + L - 1) * EXPF2(L-1,M) *   & !
                          YLM(L-1,M) - FLOAT(L + M - 1) *         & !
                          EXPF2(L-2,M) * YLM(L-2,M)               & !
                        ) / ( EXPF2(L,M) * FLOAT(L - M) )           !
            YLM(L,-M) = C1 * CONJG(YLM(L,M))                        !
          END DO                                                    !
        END IF                                                      !
      END DO                                                        !
!
  10  RETURN                                                        !
!
      END SUBROUTINE HARSPH
!
!=======================================================================
!
      SUBROUTINE HARSPH2(NL,THETA,YLM,NC)
!
!  This routine computes the complex spherical harmonics using Condon and
!          Shortley phase convention.
!
!
!  This version for m=0 only
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 19 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE,TWO,FOUR
      USE COMPLEX_NUMBERS
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  NL,NC
!
      INTEGER                  ::  L
!
      REAL (WP), INTENT(IN)    ::  THETA
      REAL (WP)                ::  SQ4PI_INV,SQR3_INV,SSMALL
      REAL (WP)                ::  X,Y
!
      REAL (WP)                ::  COS,ABS,FLOAT,SQRT
!
      COMPLEX(WP), INTENT(OUT) ::  YLM(0:NL,-NL:NL)
!
      COMPLEX(WP)              ::  CMPLX
!
      DATA SQ4PI_INV  / 0.28209479177387814347403972578E0_WP /      ! sqrt(1 / 4*pi)
      DATA SQR3_INV   / 0.48860251190291992158638462284E0_WP /      ! sqrt(3 / 4*pi)
      DATA SSMALL     / 0.00010E0_WP                         /
!
      X = COS(THETA)                                                !
!
      IF(ABS(X) < SSMALL)       X =   ZERO                          !
      IF(ABS(X + ONE) < SSMALL) X = - ONE                           !
      IF(ABS(X - ONE) < SSMALL) X =   ONE                           !
!
      YLM(0,0) = CMPLX(SQ4PI_INV)                                   !
      YLM(1,0) = X * SQR3_INV                                       !
!
      DO L = 2, NC                                                  !
        Y        = ONE / FLOAT(L)                                   !
        YLM(L,0) = X * SQRT(FOUR - Y * Y) * YLM(L-1,0) -          & !
                   (ONE - Y) * SQRT( ONE + TWO / (FLOAT(L) -      & !
                                     1.5E0_WP) )*  YLM(L-2,0)       !
      END DO                                                        !
!
      END SUBROUTINE HARSPH2
!
!=======================================================================
!
      SUBROUTINE HARSPH3(NL,THETA,PHI,YLM2,NC)
!
!  This routine computes the complex spherical harmonics using Condon and
!                  Shortley phase convention.
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 26 May 2021
!
!
      USE DIMENSION_CODE,       ONLY : LINMAX
!
      USE STORE_COEF,           ONLY : EXPF2,FSQ
!
      USE REAL_NUMBERS,         ONLY : ZERO,ONE,TWO,FOUR,HALF,SMALL
      USE COMPLEX_NUMBERS
      USE PI_ETC,               ONLY : PI
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  NL,NC
!
      INTEGER                  ::  L,M,IL,IND
!
      REAL (WP), INTENT(IN)    ::  THETA,PHI
      REAL (WP)                ::  SQ4PI_INV,SQR3_INV,SSMALL
      REAL (WP)                ::  X,Y,C2,C1
!
      REAL (WP)                ::  COS,ABS,FLOAT,SQRT
!
      COMPLEX(WP), INTENT(OUT) ::  YLM2(LINMAX)
!
      COMPLEX(WP)              ::  YLM(0:NL,-NL:NL)
      COMPLEX(WP)              ::  COEF,YMM,YMMP,C
!
      COMPLEX(WP)              ::  CMPLX,CONJG,EXP
!
      DATA SQ4PI_INV  / 0.28209479177387814347403972578E0_WP /      ! sqrt(1 / 4*pi)
      DATA SQR3_INV   / 0.48860251190291992158638462284E0_WP /      ! sqrt(3 / 4*pi)
      DATA SSMALL     / 0.00010E0_WP                         /
!
      X = COS(THETA)                                                !
!
      IF(ABS(X) < SSMALL)       X =   ZERO                          !
      IF(ABS(X + ONE) < SSMALL) X = - ONE                           !
      IF(ABS(X - ONE) < SSMALL) X =   ONE                           !
!
      YLM(0,0) = CMPLX(SQ4PI_INV)                                   !
      YLM(1,0) = X * SQR3_INV                                       !
!
      DO L = 2, NC
        Y        = ONE / FLOAT(L)                                   !
        YLM(L,0) = X * SQRT(FOUR - Y * Y) * YLM(L-1,0) -          & !
                  (ONE - Y) * SQRT( ONE + TWO / (FLOAT(L) -       & !
                                    1.5E0_WP) ) * YLM(L-2,0)        !
      END DO                                                        !
!
      C2 = - ONE                                                    !
      IF((THETA >= ZERO) .AND. (THETA <= PI)) THEN                  !
        C = - HALF * SQRT(ONE - X * X) * EXP(IC * PHI)              !
      ELSE                                                          !
        C =   HALF * SQRT(ONE - X * X) * EXP(IC * PHI)              !
      END IF                                                        !
!
      C1   = ONE                                                    !
      COEF = ONEC                                                   !
!
      DO M = 1, NC                                                  !
        C1          = C1 * C2                                       !
        COEF        = COEF * C                                      !
        YMM         = SQ4PI_INV * COEF * FSQ(M)                     !
        YLM(M,M)    = YMM                                           !
        YLM(M,-M)   = C1 * CONJG(YMM)                               !
        YMMP        = X * SQRT( FLOAT(M + M + 3) ) * YMM            !
        YLM(M+1,M)  = YMMP                                          !
        YLM(M+1,-M) = C1 * CONJG(YMMP)                              !
        IF(M < NC-1) THEN                                           !
          DO L = M+2, NC                                            !
            YLM(L,M)  = ( X * (L + L - 1) * EXPF2(L-1,M) *        & !
                          YLM(L-1,M) -                            & !
                          (L + M - 1) * EXPF2(L-2,M) * YLM(L-2,M) & !
                        ) / ( EXPF2(L,M) * (L - M) )                !
            YLM(L,-M) = C1 * CONJG(YLM(L,M))                        !
          END DO                                                    !
        END IF                                                      !
      END DO                                                        !
!
      DO L = 0, NC                                                  !
        IL = L * L + L + 1                                          !
        DO M = - L, L                                               !
          IND       = IL + M                                        !
          YLM2(IND) = YLM(L,M)                                      !
        END DO                                                      !
      END DO                                                        !
!
      END SUBROUTINE HARSPH3
!
!=======================================================================
!
      SUBROUTINE SPH_HAR(NL,X,CF,YLM,NC)
!
!  This routine computes the complex spherical harmonics using Condon and
!                  Shortley phase convention.
!
!  If the angular direction R=(THETAR,PHIR)  is given in cartesian
!      coordinates by (XR,YR,ZR), the arguments of the subroutine are :
!
!                    X  = ZR         = cos(THETAR)
!                    CF = XR + i YR  = sin(THETAR)*exp(i PHIR)
!
!          NL is the dimensioning of the YLM array and NC is
!                   the maximum l value to be computed.
!
!
!  Input parameters:
!
!       * NL       : dimensioning of YLM
!       * X        : cos(THETAR)
!       * CF       : sin(THETAR)*exp(i PHIR)
!       * NC       : l cut-off
!
!
!
!  Output parameters:
!
!       * YLM      : spherical harmonics
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 26 May 2021
!
!
      USE STORE_COEF,           ONLY : EXPF2,FSQ
!
      USE REAL_NUMBERS,         ONLY : ZERO,ONE,TWO,FOUR,HALF
      USE COMPLEX_NUMBERS
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  NL,NC
!
      INTEGER                   ::  L,M
!
      REAL (WP), INTENT(IN)     ::  X
!
      REAL (WP)                 ::  SQ4PI_INV,SQR3_INV
      REAL (WP)                 ::  Y,C2,C1
!
      REAL (WP)                 ::  COS,ABS,FLOAT,SQRT
!
      COMPLEX (WP), INTENT(IN)  ::  CF
      COMPLEX (WP), INTENT(OUT) ::  YLM(0:NL,-NL:NL)
!
      COMPLEX (WP)              ::  COEF,YMM,YMMP,C
!
      COMPLEX (WP)              ::  CMPLX,CONJG
!
      DATA SQ4PI_INV  / 0.28209479177387814347403972578E0_WP /      ! sqrt(1 / 4*pi)
      DATA SQR3_INV   / 0.48860251190291992158638462284E0_WP /      ! sqrt(3 / 4*pi)
!
      YLM(0,0) = CMPLX(SQ4PI_INV)                                   !
      YLM(1,0) = CMPLX(X * SQR3_INV)                                !
!
      DO L = 2, NC
        Y        = ONE / FLOAT(L)                                   !
        YLM(L,0) = X * SQRT(FOUR - Y * Y) * YLM(L-1,0) -          & !
                  (ONE - Y) * SQRT( ONE + TWO / (FLOAT(L) -       & !
                                    1.5E0_WP) ) * YLM(L-2,0)        !
      END DO                                                        !
!
      C2 = - ONE                                                    !
      C  = - HALF * CF                                              !
!
      C1   = ONE                                                    !
      COEF = ONEC                                                   !
!
      DO M = 1, NC                                                  !
        C1          = C1 * C2                                       !
        COEF        = COEF * C                                      !
        YMM         = SQ4PI_INV * COEF * FSQ(M)                     !
        YLM(M,M)    = YMM                                           !
        YLM(M,-M)   = C1 * CONJG(YMM)                               !
        YMMP        = X * SQRT( FLOAT(M + M + 3) ) * YMM            !
        YLM(M+1,M)  = YMMP                                          !
        YLM(M+1,-M) = C1 * CONJG(YMMP)                              !
        IF(M < NC-1) THEN                                           !
          DO L = M+2, NC                                            !
            YLM(L,M)  = ( X * (L + L - 1) * EXPF2(L-1,M) *        & !
                          YLM(L-1,M) -                            & !
                          (L + M - 1) * EXPF2(L-2,M) * YLM(L-2,M) & !
                        ) / ( EXPF2(L,M) * (L - M) )                !
            YLM(L,-M) = C1 * CONJG(YLM(L,M))                        !
          END DO                                                    !
        END IF                                                      !
      END DO                                                        !
!
      END SUBROUTINE SPH_HAR
!
END MODULE SPHERICAL_HARMONICS
