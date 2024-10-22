!
!=======================================================================
!
MODULE ANGULAR_MOMENTUM
!
!  This module provides subroutine to computes angular momentum
!    functions, namely:
!
!  1)  SUBROUTINE GAUNT(L2,M2,L1,M1,GNT)                      <-- Gaunt coefficient
!  2)  SUBROUTINE STORE_GAUNT                                 <-- storage of Gaunt coefficient
!  3)  SUBROUTINE THREE_J(L2,ML2,L1,ML1,W3J)                  <-- 3j symbols
!  4)  FUNCTION SIXJ_IN(J1,J2,L1,L2,L3)                       <-- 6j initial value for J=J1+J2
!  5)  SUBROUTINE N_J(J1,MJ1,J2,MJ2,MJ6,NJ,I_INT,N_IN)        <-- Clebsch-Gordan, 3j and 6j
!  6)  SUBROUTINE CG_6J(J1,J2,JA,JC,J,CG6J)                   <-- 2 Clebsch-Gordan * 6j
!  7)  SUBROUTINE NINE_J(J1,J2,J12,J3,J4,J34,J13,J24,J,NINEJ) <-- 9j symbols
!  8)  SUBROUTINE TENS_PROD(F1,NL1,F2,NL2,F)                  <-- tensor product
!
!
      USE ACCURACY_REAL
      USE LOGAMMA
      USE CALC_LOGAMMA
      USE DIMENSION_CODE,      ONLY : L_MAX
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE GAUNT(L2,M2,L1,M1,GNT)
!
!   This subroutine calculates the Gaunt coefficient G(L2,L3|L1)
!    using a downward recursion scheme due to Schulten and Gordon
!    for the Wigner's 3j symbols. The result is stored as GNT(L3),
!    making use of the selection rule M3 = M1 - M2.
!
!
!   Ref. :  K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
!
!
!  Input variables:
!
!     L2        :  \
!     M2        :   \  arguments of the Gaunt coefficient
!     L1        :   /          G(L2,L3|L1)
!     M1        :  /
!
!
!  Output variable:
!
!     GNT       : Gaunt coefficient G(L2,L3|L1)
!
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,TWO,HALF
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  L2,M2,L1,M1
!
      INTEGER                  ::  L12,K12,LM1,LM2,KM1,KM2
      INTEGER                  ::  M3,JMIN
      INTEGER                  ::  J,J1,J2
!
      INTEGER                  ::  MOD
!
      REAL (WP), INTENT(OUT)   ::  GNT(0:N_GAUNT)
!
      REAL (WP)                ::  F(0:N_GAUNT),G(0:N_GAUNT),A(0:N_GAUNT)
      REAL (WP)                ::  A1(0:N_GAUNT),B(0:N_GAUNT)
      REAL (WP)                ::  PI4
      REAL (WP)                ::  COEF,GND
      REAL (WP)                ::  D1,D2,D3,D4,D5,D6
!
      REAL (WP)                ::  ABS,SQRT,FLOAT,EXP,MAX
!
      DATA PI4 / 12.5663706143591729538505735331180E0_WP /          ! 4 pi
!
      L12 = L1 + L2                                                 !
      K12 = L1 - L2                                                 !
!
      DO J = 1 ,N_GAUNT                                             !
        GNT(J) = ZERO                                               !
      END DO                                                        !
!
      IF((ABS(M1) > L1) .OR. (ABS(M2) > L2)) GO TO 10               !
!
      M3  = M1 - M2                                                 !
      LM1 = L1 + M1                                                 !
      LM2 = L2 + M2                                                 !
      KM1 = L1 - M1                                                 !
      KM2 = L2 - M2                                                 !
!
      IF(MOD(M1,2) == 0) THEN                                       !
        COEF =   SQRT(FLOAT((2 * L1 + 1) * (2 * L2 + 1)) / PI4)     !
      ELSE                                                          !
        COEF = - SQRT(FLOAT((2 * L1 + 1) * (2 * L2 + 1)) / PI4)     !
      END IF                                                        !
!
      F(L12+1)  = ZERO                                              !
      G(L12+1)  = ZERO                                              !
      A(L12+1)  = ZERO                                              !
      A1(L12+1) = ZERO                                              !
!
      D1 = GLD(2*L2+1,1)   - GLD(2*L12+2,1)                         !
      D2 = GLD(2*L1+1,1)   - GLD(LM2+1,1)                           !
      D3 = GLD(L12+M3+1,1) - GLD(KM2+1,1)                           !
      D4 = GLD(L12-M3+1,1) - GLD(LM1+1,1)   - GLD(KM1+1,1)          !
!
      IF(MOD(KM1-KM2,2) == 0) THEN                                  !
        F(L12) =   SQRT( EXP(D1 + D2 + D3 + D4) )                   !
      ELSE                                                          !
        F(L12) = - SQRT( EXP(D1 + D2 + D3 + D4) )                   !
      END IF                                                        !
!
      D5 = HALF * ( GLD(2*L1+1,1) + GLD(2*L2+1,1) - GLD(2*L12+2,1) )!
      D6 = GLD(L12+1,1) - GLD(L1+1,1) - GLD(L2+1,1)                 !
!
      IF(MOD(K12,2) == 0) THEN                                      !
        G(L12) =   EXP(D5 + D6)                                     !
      ELSE                                                          !
        G(L12) = - EXP(D5 + D6)                                     !
      END IF                                                        !
!
      A(L12)  = TWO * SQRT(                                       & !
                FLOAT( L1 * L2 * (2*L12 + 1) *                    & !
                       (L12 * L12 - M3 * M3)                      & !
                     )                                            & !
                          )                                         !
      B(L12)  = - FLOAT( (2 * L12 + 1) * (                        & !
                 (L2 * L2 - L1 * L1 - K12) * M3 + L12 * (L12 + 1) & !
                  * (M2 + M1)                                     & !
                                         )                        & !
                       )                                            !
      A1(L12) = TWO * FLOAT(L12) * SQRT(                          & !
                FLOAT(L1 * L2 * (L12 + L12 + 1))                  & !
                                       )                            !
!
      IF(ABS(M3) <= L12) THEN                                       !
        GNT(L12) = COEF * F(L12) * G(L12) * SQRT(FLOAT(2 * L12 + 1))!
      ELSE                                                          !
        GNT(L12) = ZERO                                             !
      END IF                                                        !
!
      JMIN = MAX(ABS(K12),ABS(M3))                                  !
!
      DO J = L12-1, JMIN, -1                                        !
        J1    = J + 1                                               !
        J2    = J + 2                                               !
        A(J)  = SQRT( FLOAT( (J * J - K12 * K12) ) *              & !
                      FLOAT( (L12 + 1) * (L12 + 1)- J * J ) *     & !
                      FLOAT( J * J - M3 * M3 )                    & !
                    )                                               !
        B(J)  = - FLOAT( (J + J + 1) *                            & !
                         ( L2 * (L2 + 1) * M3 -                   & !
                           L1 * (L1+1) * M3 +                     & !
                           J * J1 * (M2 + M1) )                   & !
                       )                                            !
        A1(J) = FLOAT(J) * SQRT( FLOAT( (J * J - K12 * K12 ) *    & !
                           ( (L12 + 1) * (L12 + 1) - J * J)       & !
                                      )                           & !
                               )                                    !
        F(J)  = - ( FLOAT(J1) * A(J2) * F(J2) +                   & !
                    B(J1) * F(J1)                                 & !
                  ) / ( FLOAT(J2) * A(J1) )                         !
        G(J)  = - ( FLOAT(J1) * A1(J2) * G(J2) ) /                & !
                  (FLOAT(J2) * A1(J1) )                             !
        GND   = COEF * F(J) * G(J) * SQRT( FLOAT(J + J + 1) )       !
!
        IF(ABS(M3) <= J) THEN                                       !
          GNT(J) = GND                                              !
        ELSE                                                        !
          GNT(J) = ZERO                                             !
        END IF                                                      !
!
      END DO                                                        !
!
  10  RETURN                                                        !
!
      END SUBROUTINE GAUNT
!
!=======================================================================
!
      SUBROUTINE STORE_GAUNT
!
!   This subroutine calculates and stores the Gaunt coefficient G(L2,L3|L1)
!    using a downward recursion scheme due to Schulten and Gordon
!    for the Wigner's 3j symbols. The result is stored as GNT(L3,L2,L1)
!    and makes use of the selection rule M3 = M1 - M2.
!
!
!   Ref. :  K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
!
!
!   The Gaunt coefficients are stored in the module GAUNT_C
!
!
!   Authors: H. F. Zhao and D. Sebilleau
!
!                                          Last modified :  4 Jun 2021
!
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,TWO,HALF
!
      USE LOGAMMA,           ONLY : GLD
      USE GAUNT_C
!
      IMPLICIT NONE
!
      INTEGER                  ::  L2,M2,L1,M1
      INTEGER                  ::  IL1,IND1,IL2,IND2
      INTEGER                  ::  L12,K12,LM1,LM2,KM1,KM2
      INTEGER                  ::  L12_1,L12_2,K12_2,L12_21
      INTEGER                  ::  M3,JMIN
      INTEGER                  ::  J,J1,J2,JJ
!
      INTEGER                  ::  MOD
!
      REAL (WP)                ::  F(0:N_GAUNT),G(0:N_GAUNT),A(0:N_GAUNT)
      REAL (WP)                ::  A1(0:N_GAUNT),B(0:N_GAUNT)
      REAL (WP)                ::  PI4
      REAL (WP)                ::  COEF,GND
      REAL (WP)                ::  D1,D2,D3,D4,D5,D6
!
      REAL (WP)                ::  ABS,SQRT,FLOAT,EXP,MAX
!
      DATA PI4 / 12.5663706143591729538505735331180E0_WP /          ! 4 pi
!
      DO L1 = 0, L_MAX                                              !
        IL1 = L1 * L1 + L1 + 1                                      !
        DO M1 = -L1, L1                                             !
          IND1 = IL1 + M1                                           !
          LM1  = L1 + M1                                            !
          KM1  = L1 - M1                                            !
          DO L2 = 0, L_MAX                                          !
            IL2 = L2 * L2 + L2 + 1                                  !
!
            IF(MOD(M1,2) == 0) THEN                                 !
              COEF =   SQRT( FLOAT((L1+L1+1) * (L2+L2+1)) / PI4 )   !
            ELSE                                                    !
              COEF = - SQRT( FLOAT((L1+L1+1) * (L2+L2+1)) / PI4 )   !
            END IF                                                  !
!
            L12    = L1 + L2                                        !
            K12    = L1 - L2                                        !
            L12_1  = L12 + L12 + 1                                  !
            L12_2  = L12 * L12                                      !
            L12_21 = L12 * L12 + L12 + L12 + 1                      !
            K12_2  = K12 * K12                                      !
!
            F(L12+1)  = ZERO                                        !
            G(L12+1)  = ZERO                                        !
            A(L12+1)  = ZERO                                        !
            A1(L12+1) = ZERO                                        !
            A1(L12)   = TWO * SQRT( FLOAT(L1 * L2 * L12_1 * L12_2) )!
            D1        = GLD(L2+L2+1,1) - GLD(L12_1+1,1)             !
            D5        = HALF * ( GLD(L1+L1+1,1) +                 & !
                                 GLD(L2+L2+1,1) - GLD(L12_1+1,1) )  !
            D6        = GLD(L12+1,1) - GLD(L1+1,1) - GLD(L2+1,1)    !
!
            IF(MOD(K12,2) == 0) THEN                                !
              G(L12) =   EXP(D5 + D6)                               !
            ELSE                                                    !
              G(L12) = - EXP(D5 + D6)                               !
            END IF                                                  !
!
            DO M2 = - L2, L2                                        !
              IND2 = IL2 + M2                                       !
!
              M3  = M1 - M2                                         !
              LM2 = L2 + M2                                         !
              KM2 = L2 - M2                                         !
!
              DO J = 1, N_GAUNT                                     !
                GNT(J,IND2,IND1) = ZERO                             !
              END DO                                                !
!
              IF((ABS(M1) > L1) .OR. (ABS(M2) > L2)) GO TO 10       !
!
              D2 = GLD(L1+L1+1,1)  - GLD(LM2+1,1)                   !
              D3 = GLD(L12+M3+1,1) - GLD(KM2+1,1)                   !
              D4 = GLD(L12-M3+1,1) - GLD(LM1+1,1) - GLD(KM1+1,1)    !
!
              IF(MOD(KM1-KM2,2) == 0) THEN                          !
                F(L12) =   SQRT( DEXP(D1 + D2 + D3 + D4) )          !
              ELSE                                                  !
                F(L12) = - SQRT( DEXP(D1 + D2 + D3 + D4) )          !
              END IF                                                !
!
              A(L12) = TWO * SQRT( FLOAT( L1 * L2 * L12_1 *       & !
                                          (L12_2-M3*M3) ) )         !
              B(L12) = - FLOAT( L12_1 * (                         & !
                                          (L2*L2-L1*L1-K12) *     & !
                                          M3 + L12 * (L12 + 1) *  & !
                                          (M2+M1)                 & !
                                        ) )                         !
!
              IF(ABS(M3) <= L12) THEN                               !
                GNT(L12,IND2,IND1) = COEF * F(L12) * G(L12) *     & !
                                     SQRT( FLOAT(L12_1) )           !
              END IF                                                !
!
              JMIN = MAX(ABS(K12),ABS(M3))                          !
!
              DO J = L12-1, JMIN, -1                                !
                J1    = J + 1                                       !
                J2    = J + 2                                       !
                JJ    = J * J                                       !
                A1(J) = SQRT( FLOAT( JJ * (JJ - K12_2) *          & !
                                          (L12_21 - JJ)           & !
                                   ) )                              !
                A(J)  = SQRT( FLOAT( (JJ - K12_2)  *              & !
                                     (L12_21 - JJ) *              & !
                                     (JJ - M3 * M3)               & !
                                   ) )                              !
                B(J)  = - FLOAT( (J + J1) * ( L2 * (L2 + 1) * M3 -& !
                                              L1 * (L1 + 1) * M3 +& !
                                              J * J1 * (M2 + M1)  & !
                                            ) )                     !
                F(J)  = - ( FLOAT(J1) * A(J2) * F(J2) +           & !
                            B(J1) * F(J1) ) / ( FLOAT(J2) * A(J1) ) !
                G(J)  = - ( FLOAT(J1) * A1(J2) * G(J2) ) /        & !
                          ( FLOAT(J2) * A1(J1) )                    !
!
                IF(ABS(M3) <= J) THEN
                  GNT(J,IND2,IND1) = COEF * F(J) * G(J) *         & !
                                     SQRT( FLOAT(J + J1) )          !
                END IF                                              !
              END DO                                                !
!
            END DO                                                  ! end of M2 loop
          END DO                                                    ! end of L2 loop
        END DO                                                      ! end of M1 loop
      END DO                                                        ! end of L1 loop
!
  10  RETURN
!
      END SUBROUTINE STORE_GAUNT
!
!=======================================================================
!
      SUBROUTINE THREE_J(L2,M2,L1,M1,W3J)
!
!  This function calculates the Wigner's 3j coefficient
!
!            (  L1 L2  L )
!            (  M1 M2  M )
!
!    using a downward recursion scheme due to Schulten and Gordon.
!
!  The result is stored as W3J(L3), making use of the selection rule
!
!                        ML3 = ML1 - MJ2.
!
!  Warning: this is the INTEGER version
!
!
!   Ref. :  K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO,TWO,HALF
!
      IMPLICIT NONE
!
      INTEGER               ::  L2,M2,L1,M1
      INTEGER               ::  L12,K12
      INTEGER               ::  M3
      INTEGER               ::  LM1,LM2,KM1,KM2
      INTEGER               ::  J1,J2,JMIN
      INTEGER               ::  J,LMIN
!
      INTEGER               ::  ABS,MOD
!
      REAL (WP)             ::  W3J(0:N_GAUNT)
      REAL (WP)             ::  D1,D2,D3,D4,D5,D6
      REAL (WP)             ::  F(0:N_GAUNT)
      REAL (WP)             ::  A(0:N_GAUNT),B(0:N_GAUNT)
!
      REAL (WP)             ::  SQRT,EXP,FLOAT
!
      L12 = L1+L2                                                   !
      K12 = L1-L2                                                   !
!
!  Initialization
!
      DO J=1,N_GAUNT                                                !
        W3J(J) = ZERO                                               !
      END DO                                                        !
!
      IF((ABS(M1) > L1) .OR. (ABS(M2) > L2)) GO TO 10               !
!
      M3  = M1 - M2                                                 !
!
      LM1 = L1 + M1                                                 !
      LM2 = L2 + M2                                                 !
      KM1 = L1 - M1                                                 !
      KM2 = L2 - M2                                                 !
!
      F(L12+1)  = ZERO                                              !
      A(L12+1)  = ZERO                                              !
!
      D1 = GLD(2*L2+1,1)   - GLD(2*L12+2,1)                         !
      D2 = GLD(2*L1+1,1)   - GLD(LM2+1,1)                           !
      D3 = GLD(L12+M3+1,1) - GLD(KM2+1,1)                           !
      D4 = GLD(L12-M3+1,1) - GLD(LM1+1,1)   - GLD(KM1+1,1)          !
!
      IF(MOD(KM1-KM2,2) == 0) THEN                                  !
        F(L12) =   SQRT( EXP(D1 + D2 + D3 + D4) )                   !
      ELSE                                                          !
        F(L12) = - SQRT( EXP(D1 + D2 + D3 + D4) )                   !
      END IF                                                        !
!
      A(L12)  = TWO * SQRT( FLOAT( L1 * L2 * (1 + 2 * L12) *      & !
                                  (L12 * L12 - M3 * M3) )         & !
                          )                                         !
      B(L12)  = - FLOAT( (2 * L12 + 1) * (                        & !
                         (L2 * L2 - L1 * L1 - K12) * M3 +         & !
                          L12 * (L12 + 1) * (M2 + M1) )           & !
                       )                                            !
!
      IF(ABS(M3) <= L12) THEN                                       !
        W3J(L12) = F(L12)                                           !
      ELSE                                                          !
        W3J(L12) = ZERO                                             !
      END IF                                                        !
!
      LMIN = MAX(ABS(K12),ABS(M3))                                  !
!
      DO J=L12-1,LMIN,-1                                            !
        J1 = J + 1                                                  !
        J2 = J + 2                                                  !
        A(J)  = SQRT( FLOAT( (J * J - K12 * K12) ) *              & !
                      FLOAT( (L12 + 1) * (L12 + 1) - J * J) *     & !
                      FLOAT( J * J - M3 * M3 )                    & !
                     )                                              !
        B(J)  = - FLOAT( (2 * J + 1) *                            & !
                         ( L2 * (L2+1) * M3 - L1 * (L1+1) * M3 +  & !
                          J * J1 * (M2 + M1) )                    & !
                       )                                            !
        F(J)  = -( FLOAT(J1) * A(J2) * F(J2) + B(J1) * F(J1) ) /  & !
                 ( FLOAT(J2) * A(J1) )                              !
!
        IF(ABS(M3) <= J) THEN                                       !
          W3J(J) = F(J)                                             !
        ELSE                                                        !
          W3J(J) = ZERO                                             !
        END IF                                                      !
!
      END DO                                                        !
!
  10  RETURN                                                        !
!
      END SUBROUTINE THREE_J
!
!=======================================================================
!
      FUNCTION SIXJ_IN(J1,J2,L1,L2,L3)
!
!  This function calculates the initial value {J1 J2 L1+L2}
!                                             {L1 L2   L3 }
!
!  A 6j symbol {J1 J2 J3} is non zero only if
!              {J4 J5 J6}
!
!   (J1,J2,J3),(J4,J5,J3),(J2,J4,J6) and (J1,J5,J6) satisfy the triangular inequality :
!
!       (a,b,c) non zero if |a-b| <= c <= (a+b) . This means also that (a+b) and c must
!       have the same nature (integer or half-integer).
!
!   (J1,J2,J3) and (J4,J5,J3) are taken care of by the bounds of J3, JJ_MIN and JJ_MAX,
!       as chosen in the N_J routine. Here we check the two last ones.
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,         ONLY : ZERO,ONE
!
      IMPLICIT NONE
!
      INTEGER               ::  IZERO
      INTEGER               ::  LL1_2,LL2_2,MSIGN
!
      INTEGER               ::  INT,MOD
!
      REAL (WP), INTENT(IN) ::  J1,J2,L1,L2,L3
      REAL (WP)             ::  SIXJ_IN
!
      REAL (WP)             ::  SMALL,SIGNE
      REAL (WP)             ::  D1,D2,D3,D4,D5,D6,D7
!
      REAL (WP)             ::  ABS,SIGN,SQRT,EXP
!
      DATA SMALL / 0.0001E0_WP /
!
      IZERO = 0                                                     !
!
!  Check for unphysical values of L3
!
      IF(ABS(J2-L1) > L3)     IZERO = 1                             !
      IF(J2 + L1 < L3)        IZERO = 1                             !
      IF(IG(J2+L1) /= IG(L3)) IZERO = 1                             !
      IF(ABS(J1 - L2) > L3)   IZERO = 1                             !
      IF(J1 + L2 < L3)        IZERO = 1                             !
      IF(IG(J1+L2) /= IG(L3)) IZERO = 1                             !
!
      IF(IZERO == 1) THEN                                           !
        SIXJ_IN = ZERO                                              !
      ELSE                                                          !
!
!  Storage indices of the angular momenta.
!
        LL1_2 = INT(L1 + L1 + SIGN(SMALL,L1))                       !
        LL2_2 = INT(L2 + L2 + SIGN(SMALL,L2))                       !
!
        MSIGN = INT(J1 + J2 + L1 + L2 + SIGN(SMALL,J1+J2+L1+L2))    !
        IF(MOD(MSIGN,2) == 0) THEN                                  !
          SIGNE =   ONE                                             !
        ELSE                                                        !
          SIGNE = - ONE                                             !
        END IF                                                      !
!
        D1 = GLD(LL1_2+1,1) + GLD(LL2_2+1,1) - GLD(LL1_2+LL2_2+2,1) !
        D2 = GLD(INT(J1+J2+L1+L2)+2,IG(J1+J2+L1+L2)) -            & !
             GLD(INT(J1+J2-L1-L2)+1,IG(J1+J2-L1-L2))                !
        D3 = GLD(INT(J1-J2+L1+L2)+1,IG(J1-J2+L1+L2)) -            & !
             GLD(INT(J1+L2-L3)+1,IG(J1+L2-L3))                      !
        D4 = GLD(INT(J2-J1+L1+L2)+1,IG(J2-J1+L1+L2)) -            & !
             GLD(INT(-J1+L2+L3)+1,IG(-J1+L2+L3))                    !
        D5 = GLD(INT(J1-L2+L3)+1,IG(J1-L2+L3)) -                  & !
             GLD(INT(J1+L2+L3)+2,IG(J1+L2+L3))                      !
        D6 = GLD(INT(J2+L3-L1)+1,IG(J2+L3-L1)) -                  & !
             GLD(INT(J2-L3+L1)+1,IG(J2-L3+L1))                      !
        D7 = GLD(INT(L1+L3-J2)+1,IG(L1+L3-J2)) +                  & !
             GLD(INT(L1+J2+L3)+2,IG(L1+J2+L3))                      !
!
        SIXJ_IN = SIGNE * SQRT(                                   & !
                  EXP(D1 + D2 + D3 + D4 + D5 + D6 - D7)           & !
                              )                                     !
!
      END IF                                                        !
!
      END FUNCTION SIXJ_IN
!
!=======================================================================
!
      FUNCTION IG(J)
!
!  This function is returns the value 1 if J is an integer
!   and 2 if it is a half-integer
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 11 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : SMALL
!
      IMPLICIT NONE
!
      INTEGER               ::  IG
!
      INTEGER               ::  LL
!
      INTEGER               ::  INT
!
      REAL (WP), INTENT(IN) ::  J
!
      REAL (WP)             ::  JJ
!
      JJ = ABS(J+J)                                                 !
!
      LL = INT(JJ + SMALL)                                          !
!
      IF(MOD(LL,2) == 0) THEN                                       !
        IG = 1                                                      !
      ELSE                                                          !
        IG = 2                                                      !
      END IF                                                        !
!
      END FUNCTION IG
!
!=======================================================================
!
      SUBROUTINE N_J(J1,MJ1,J2,MJ2,MJ6,NJ,I_INT,N_IN)
!
!   This subroutine calculates Wigner's 3j and 6j coefficients
!    using a downward recursion scheme due to Schulten and Gordon.
!    The 3j are defined as (J1 J2 J) where in fact L1=MJ1, etc are
!                          (L1 L2 L)
!    azimuthal quantum numbers, and the 6j as {J1 J2 J} where now
!                                             {L1 L2 L}
!    J1, L1, etc are the same kind of orbital quantum numbers.
!    The result is stored as NJ(J).
!
!    The parameter N allows to choose between 3j and 6j calculation, and
!    Clebsch-Gordan. It can take the values :
!
!                    N = 2 ----> Clebsch-Gordan
!                    N = 3 ----> Wigner's 3j
!                    N = 6 ----> Wigner's 6j
!
!    The Clebsch-Gordan coefficients are related to Wigner's 3j through :
!
!       CG(J1,M1,J2,M2|J,MJ) = ( J1 J2  J  )*sqrt(2*J+1)*(-1)**(J1-J2+MJ)
!                              ( M1 M2 -MJ )
!    I_INT is a flag that returns 1 if the index J of the nj symbol
!      is integer and 0 if it is a half integer.
!
!   Note : For 3j, MJ6 is ignored while for 6j, we have :
!
!                J1=J1   MJ1=L1   J2=J2   MJ2=L2   MJ6=L
!
!   Ref. : K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
!
!
!   Author :  D. Sébilleau
!
!                      Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,   ONLY : ZERO,ONE,TWO,HALF,SMALL
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)    ::  N_IN
      INTEGER, INTENT(OUT)   ::  I_INT
!
      INTEGER                ::  N_OU,I_CG
      INTEGER                ::  L1,L2
      INTEGER                ::  IG1,IG2,IGG
      INTEGER                ::  L12,K12,LM1,LM2,KM1,KM2
      INTEGER                ::  M,L12M,K12M,L1_2,L2_2,L12_2
      INTEGER                ::  L,N12,LMIN,LP1,LP2
      INTEGER                ::  LJ12,LL12
      INTEGER                ::  LJJ_MIN,LJJ_MAX
!
      INTEGER                ::  INT
!
      REAL (WP), INTENT(IN)  ::  J1,MJ1,J2,MJ2,MJ6
      REAL (WP), INTENT(OUT) ::  NJ(0:N_GAUNT)
!
      REAL (WP)              ::  J,MJ,JP1,JP2
      REAL (WP)              ::  JL12,JK12,SIG
      REAL (WP)              ::  JJ1,JJ2,JL1,JL2,JL3
      REAL (WP)              ::  JJ12,JJ_MIN,JJ_MAX
      REAL (WP)              ::  F(0:N_GAUNT),A(0:N_GAUNT),B(0:N_GAUNT)
      REAL (WP)              ::  DIF1,DIF2
      REAL (WP)              ::  D1,D2,D3,D4
!
      REAL (WP)              ::  FLOAT,SIGN,SQRT,EXP,MIN,MAX
!
      IF(N_IN == 2) THEN                                            !
        N_OU = 3                                                    !
        I_CG = 1                                                    !
      ELSE                                                          !
        N_OU = N_IN                                                 !
        I_CG = 0                                                    !
      END IF                                                        !
!
      IF(N_OU == 3) THEN                                            !
!
!------------------------------  3j case  ---------------------------------
!
!
!  Test to check if J1 and J2 are integer or semi-integer
!
!     Integer      : IG=1
!     Half-integer : IG=2
!
!  Each angular momentum J is represented by the integer index L
!
        L1   = INT(J1 + SMALL)                                      !
        L2   = INT(J2 + SMALL)                                      !
        DIF1 = J1 - FLOAT(L1)                                       !
        DIF2 = J2 - FLOAT(L2)                                       !
!
!  IGx is a flag telling the code which case of Gamma function to use :
!
!     IGx = 1 : integer case
!     IGx = 2 : half-integer case
!
        IF(ABS(DIF1) < SMALL) THEN                                  !
          IG1 = 1                                                   !
        ELSE                                                        !
          IG1 = 2                                                   !
        END IF                                                      !
        IF(ABS(DIF2) < SMALL) THEN                                  !
          IG2 = 1                                                   !
        ELSE                                                        !
          IG2 = 2                                                   !
        END IF                                                      !
        IF(IG1 == IG2) THEN                                         !
          IGG = 1                                                   !
        ELSE                                                        !
          IGG = 2                                                   !
        END IF                                                      !
!
!  Here, we assume that (J1,J2) are both either integer or half-integer
!  If J is integer, the corresponding index is L = j (for loops or storage)
!    while if J is an half-integer, this index is L=  j - 1/2 = int(j)
!
!  Integer indices are used for loops and for storage while true values
!    are used for the initial values. When J1 and J2 are both half-integers,
!    the values of J are integer and L should be increased by 1
!
        JL12 = J1 + J2                                              !
        JK12 = J1 - J2                                              !
!
        L12 = INT(JL12 + SIGN(SMALL,JL12))                          !
        K12 = INT(JK12 + SIGN(SMALL,JK12))                          !
!
        LM1 = INT(J1+MJ1 + SIGN(SMALL,J1+MJ1))                      !
        LM2 = INT(J2+MJ2 + SIGN(SMALL,J2+MJ2))                      !
        KM1 = INT(J1-MJ1 + SIGN(SMALL,J1-MJ1))                      !
        KM2 = INT(J2-MJ2 + SIGN(SMALL,J2-MJ2))                      !
!
        MJ = - MJ1 - MJ2                                            !
!
        M     = INT(MJ + SIGN(SMALL,MJ))                            !
        L12M  = INT(JL12 + MJ + SIGN(SMALL,JL12+MJ))                !
        K12M  = INT(JL12 - MJ + SIGN(SMALL,JL12-MJ))                !
        L1_2  = INT(J1 + J1 + SIGN(SMALL,J1))                       !
        L2_2  = INT(J2 + J2 + SIGN(SMALL,J2))                       !
        L12_2 = INT(JL12 + JL12 + SIGN(SMALL,JL12))                 !
!
        IF(IG(JL12) == 1) THEN                                      !
          I_INT = 1                                                 !
        ELSE                                                        !
          I_INT = 0                                                 !
        END IF                                                      !
!
!  Initialisation of the 3j symbol NJ(J) = (J1  J2   J)
!                                          (MJ1 MJ2 MJ)
!
        DO L = 0, L12                                               !
          NJ(L) = ZERO                                              !
        END DO                                                      !
!
        IF((ABS(MJ1) > J1) .OR. (ABS(MJ2) > J2)) GO TO 10           !
!
!  Initial values (J1+J2+1) and (J1+J2) for J to be used in the downward
!     recursion scheme. This scheme writes as
!
!     J A(J+1) NJ(J+1) + B(J) NJ(J) + (J+1) A(J) NJ(J-1) = 0
!
        F(L12+1) = ZERO                                             !
        A(L12+1) = ZERO                                             !
        D1       = GLD(L2_2+1,1) - GLD(L12_2+2,1)                   !
        D2       = GLD(L1_2+1,1) - GLD(LM2+1,1)                     !
        D3       = GLD(L12M+1,1) - GLD(KM2+1,1)                     !
        D4       = GLD(K12M+1,1) - GLD(LM1+1,1) - GLD(KM1+1,1)      !
!
        N12 = INT(JK12-MJ + SIGN(SMALL,JK12-MJ))                    !
!
        IF(I_CG == 1) THEN                                          !
          IF(MOD(N12,2) == 0) THEN                                  !
            SIG =   ONE                                             !
          ELSE                                                      !
            SIG = - ONE                                             !
          END IF                                                    !                                        !
        END IF                                                      !
!
        IF(MOD(N12,2) == 0) THEN                                    !
          F(L12 )=   SQRT(EXP(D1 + D2 + D3 + D4))                   !
        ELSE                                                        !
          F(L12) = - SQRT(EXP(D1 + D2 + D3 + D4))                   !
       END IF                                                       !
!
        A(L12) = TWO * SQRT( J1 * J2 * (ONE + TWO * JL12) *       & !
                             (JL12 * JL12 - MJ * MJ) )              !
        B(L12) = - (TWO * JL12 + ONE) * (                         & !
                 (J1 * J1 - J2 * J2 + JK12) * MJ - JL12 *         & !
                 (JL12 + ONE ) * (MJ2 - MJ1)                      & !
                                        )                           !
!
        IF(ABS(M).LE.L12) THEN                                      !
          IF(I_CG == 0) THEN                                        !
            NJ(L12) = F(L12)                                        !
          ELSE                                                      !
            NJ(L12) = F(L12) * SIG * SQRT(JL12 + JL12 + ONE)        !
          END IF                                                    !
        ELSE                                                        !
          NJ(L12) = ZERO                                            !
        END IF                                                      !
!
        LMIN = MAX(ABS(K12),ABS(M))                                 !
!
!  Downward recursion for NJ(J)
!
        DO L = L12 - 1, LMIN, -1                                    !
          LP1 = L + 1                                               !
          LP2 = L + 2                                               !
!
!  Value of the angular momentum J corresponding to the loop index L
!
          IF(IGG == 1) THEN                                         !
            J   = FLOAT(L)                                          !
            JP1 = FLOAT(LP1)                                        !
            JP2 = FLOAT(LP2)                                        !
          ELSE                                                      !
            J   = FLOAT(L)   + HALF                                 !
            JP1 = FLOAT(LP1) + HALF                                 !
            JP2 = FLOAT(LP2) + HALF                                 !
          END IF                                                    !
!
          A(L) =  SQRT( (J * J - JK12 * JK12) *                   & !
                        ((JL12 + ONE) * (JL12 + ONE)- J * J) *    & !
                        (J * J - MJ * MJ)                         & !
                      )                                             !
          B(L) = - (TWO * J + ONE) * ( J1 * (J1 + ONE) * MJ -     & !
                    J2 * (J2 + ONE) * MJ - J * JP1 * (MJ2 - MJ1)  & !
                                     )                              !
          F(L) = - ( JP1 * A(LP2) * F(LP2) + B(LP1) * F(LP1) ) /  & !
                   ( JP2 * A(LP1) )                                 !
!
          IF(ABS(MJ) <= J) THEN                                     !
            IF(I_CG == 0) THEN                                      !
              NJ(L) = F(L)                                          !
            ELSE                                                    !
              NJ(L) = F(L) * SIG * SQRT(J + J + ONE)                !
            END IF                                                  !
          ELSE                                                      !
            NJ(L) = ZERO                                            !
          END IF                                                    !
!
        END DO                                                      !
!
  10    CONTINUE                                                    !
!
      ELSE IF(N_OU == 6) THEN                                       !
!
!------------------------------  6j case  ---------------------------------
!
!  Change of notation for greater readability ---> NJ(JJ)
!
!    True angular momentum value : begins with a J (JJn,JLn)
!    Corresponding integer storage and loop index : begins by L (LJn,LLn)
!
        JJ1 = J1                                                    !
        JJ2 = J2                                                    !
        JL1 = MJ1                                                   !
        JL2 = MJ2                                                   !
        JL3 = MJ6                                                   !
!
        JJ12 = JJ1 - JJ2                                            !
        JL12 = JL1 - JL2                                            !
!
        LJ12 = INT(JJ12 + SIGN(SMALL,JJ12))                         !
        LL12 = INT(JL12 + SIGN(SMALL,JL12))                         !
!
        JJ_MIN  = MAX(ABS(LJ12),ABS(LL12))                          !
        JJ_MAX  = MIN(JJ1 + JJ2,JL1 + JL2)                          !
        LJJ_MIN = INT(JJ_MIN + SIGN(SMALL,JJ_MIN))                  !
        LJJ_MAX = INT(JJ_MAX + SIGN(SMALL,JJ_MAX))                  !
!
!  Initialisation of the 6j symbol NJ(J) = {J1  J2  J }
!                                          {L1  L2  L3}
!
        DO L = 0, LJJ_MAX                                           !
          NJ(L) = ZERO                                              !
        END DO                                                      !
!
!  Initial values (J1+J2+1) and (J1+J2) for J to be used in the downward
!     recursion scheme. This scheme writes as
!
!     J A(J+1) NJ(J+1) + B(J) NJ(J) + (J+1) A(J) NJ(J-1) = 0
!
!  There are two possible initial values as max(|J1-J2|,|L1-L2|) <= J <=
!    min(J1+J2,L1+L2) :
!
!    {J1  J2  L1+L2}  and  {J1  J2  J1+J2}  =  {L1 L2 J1+J2}
!    {L1  L2    L3 }       {L1  L2    L3 }     {J1 J2   L3 }
!
!    They can be calculated from equation (6.3.1) of Edmonds page 97
!
        F(LJJ_MAX+1) = ZERO                                         !
        A(LJJ_MAX+1) = ZERO                                         !
!
        IF(ABS(JJ_MAX-JL1-JL2) < SMALL) THEN                        !
          F(LJJ_MAX) = SIXJ_IN(JJ1,JJ2,JL1,JL2,JL3)                 !
        ELSE                                                        !
          F(LJJ_MAX) = SIXJ_IN(JL1,JL2,JJ1,JJ2,JL3)                 !
        END IF                                                      !
        NJ(LJJ_MAX) = F(LJJ_MAX)                                    !
!
        A(LJJ_MAX) = SQRT( (JJ_MAX * JJ_MAX -                     & !
                           (JJ1 - JJ2) * (JJ1 - JJ2)) *           & !
                           ( (JJ1 + JJ2 + ONE) *                  & !
                             (JJ1 + JJ2 + ONE) - JJ_MAX * JJ_MAX  & !
                           ) *                                    & !
                           ( JJ_MAX * JJ_MAX -                    & !
                             (JL1 - JL2) * (JL1 - JL2)            & !
                           ) *                                    & !
                           ( (JL1 + JL2 + ONE) *                  & !
                             (JL1 + JL2 + ONE) - JJ_MAX * JJ_MAX  & !
                           )                                      & !
                         )                                          !
        B(LJJ_MAX) = (JJ_MAX + JJ_MAX + ONE) *                    & !
                     ( JJ_MAX * (JJ_MAX + ONE) *                  & !
                      ( - JJ_MAX * (JJ_MAX + ONE) +               & !
                        JJ1 * (JJ1 + ONE) +                       & !
                        JJ2 * (JJ2 + ONE)                         & !
                      ) +                                         & !
                      JL1 * (JL1 + ONE) * (                       & !
                         JJ_MAX * (JJ_MAX + ONE) +                & !
                         JJ1 * (JJ1 + ONE) - JJ2 * (JJ2 + ONE)    & !
                                          ) +                     & !
                      JL2 * (JL2 + ONE) * (                       & !
                         JJ_MAX * (JJ_MAX + ONE) -                & !
                         JJ1 * (JJ1 + ONE) + JJ2 * (JJ2 + ONE)    & !
                                          ) -                     & !
                      (JJ_MAX + JJ_MAX) * (JJ_MAX + ONE) *        & !
                      JL3 * (JL3 + ONE)                           & !
                     )                                              !
!
      IF(IG(JJ_MAX) == 1) THEN                                      !
        I_INT = 1                                                   !
      ELSE                                                          !
        I_INT = 0                                                   !
      END IF                                                        !
!
!  Downward recurrence relation
!
        DO L = LJJ_MAX-1, LJJ_MIN, -1                               !
          LP1 = L + 1                                               !
          LP2 = L + 2                                               !
!
!  Value of the angular momentum J corresponding to the loop index L
!
          IF(IG(JJ_MAX) == 1) THEN                                  !
            J   = FLOAT(L)                                          !
            JP1 = FLOAT(LP1)                                        !
            JP2 = FLOAT(LP2)                                        !
          ELSE                                                      !
            J   = FLOAT(L)   + HALF                                 !
            JP1 = FLOAT(LP1) + HALF                                 !
            JP2 = FLOAT(LP2) + HALF                                 !
          END IF                                                    !
!
          A(L) = SQRT( ( J * J - (JJ1 - JJ2) * (JJ1 - JJ2) ) *    & !
                       ( (JJ1 + JJ2 + ONE) * (JJ1 + JJ2 + ONE) -  & !
                         J * J                                    & !
                       ) *                                        & !
                       ( J * J - (JL1 - JL2) * (JL1 - JL2) ) *    & !
                       ( (JL1 + JL2 + ONE) * (JL1 + JL2 + ONE) -  & !
                         J * J                                    & !
                       )                                          & !
                     )                                              !
          B(L) = (J + J + 1) * (                                  & !
                     J * JP1 * ( - J * JP1 + JJ1 * (JJ1 + ONE) +  & !
                                  JJ2 * (JJ2 + ONE)               & !
                               ) +                                & !
                     JL1 * (JL1 + ONE) * ( J * JP1 + JJ1 *        & !
                                           (JJ1 + ONE) -          & !
                                           JJ2 * (JJ2 + ONE)      & !
                                         ) +                      & !
                     JL2 * (JL2 + ONE) * ( J * JP1 - JJ1 *        & !
                                           (JJ1 + ONE) +          & !
                                           JJ2 * (JJ2 + ONE)      & !
                                         ) -                      & !
                     (J + J) * JP1 * JL3 * (JL3 + ONE)            & !
                               )                                    !
!
          F(L)  = -(JP1 * A(LP2) * F(LP2) + B(LP1) * F(LP1)) /    & !
                   (JP2 * A(LP1))                                   !
          NJ(L) =  F(L)                                             !
!
        END DO                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE N_J
!
!=======================================================================
!
      SUBROUTINE CG_6J(J1,J2,JA,JC,J,CG6J)
!
!   This subroutine calculates the product of two Clebsch-Gordan
!                coefficient by a 6j :
!
!      CG6J(JK) = CG(J1,0,JK,0|JA,0)*CG(J2,0,JK,0|JC,0)*{ J2 JC JK }
!                                                       { JA J1 J  }
!
!   Note : we store the CG and the 6j as a function of JK. For this we rewrite
!            the CG making use of the relation (Varshalovich et al p. 245) :
!
!   CG(J1,0,JK,0|JA,0) = (-1)**J1 * sqrt[(2*JA+1)/(2*JK+1)] * CG(J1,0,JA,0|JK,0)
!
!   Because the MJ = 0, the J must be integer
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO,ONE,TWO,SMALL
!
      IMPLICIT NONE
!
      INTEGER               ::  N1,N2,N12
      INTEGER               ::  L
      INTEGER               ::  LK
      INTEGER               ::  I_INT1,I_INT2,I_INT3
      INTEGER               ::  LK_MIN,LK_MAX
!
      INTEGER               ::  INT,MOD
!
      REAL (WP)             ::  J1,J2,JA,JC,J
      REAL (WP)             ::  CG6J(0:N_GAUNT)
      REAL (WP)             ::  CG1(0:N_GAUNT),CG2(0:N_GAUNT)
      REAL (WP)             ::  SIXJ(0:N_GAUNT)
      REAL (WP)             ::  SIG
      REAL (WP)             ::  MJ
      REAL (WP)             ::  JK,JK_MIN,JK_MAX
!
      REAL (WP)             ::  SIGN,SQRT,FLOAT
      REAL (WP)             ::  ABS,MAX,MIN
!
      N1 = 2                                                        !
      N2 = 6                                                        !
!
!  Initialization of the result array to zero
!
      DO L = 0, L_MAX                                               !
        CG6J(L) = ZERO                                              !
      END DO                                                        !
!
!  Calculation of the Clebsch-Gordan coefficients
!
!    Sign coming from relation Varshalovich et al p. 245
!
      N12 = INT( J1 + J2 + SIGN(SMALL,J1+J2) )                      !
      IF(MOD(N12,2) == 0) THEN                                      !
        SIG =  ONE                                                  !
      ELSE                                                          !
        SIG = -ONE                                                  !
      END IF                                                        !
!
      MJ = ZERO                                                     !
!
      CALL N_J(J1,MJ,JA,MJ,MJ,CG1,I_INT1,N1)                        !
      CALL N_J(J2,MJ,JC,MJ,MJ,CG2,I_INT2,N1)                        !
!
!  The angular momentum JK must satisfy both Clebsch-Gordan. Therefore,
!
!            JK in [max(|J1-JA|,|J2-JC|),min(J1+JA,J2+JC)]
!
      JK_MIN = MAX(ABS(J1-JA),ABS(J2-JC))                           !
      JK_MAX = MIN(J1+JA,J2+JC)                                     !
!
      LK_MIN = INT( JK_MIN + SIGN(SMALL,JK_MIN) )                   !
      LK_MAX = INT( JK_MAX + SIGN(SMALL,JK_MAX) )                   !
!
!  Calculation of the 6j coefficients
!
      CALL N_J(J2,JA,JC,J1,J,SIXJ,I_INT3,N2)                        !
!
!  Computation of the coefficients
!
      DO LK = LK_MIN, LK_MAX                                        !
        JK       = FLOAT(LK)                                        !
        CG6J(LK) = SIXJ(LK) * CG1(LK) * CG2(LK) * SIG *           & !
                   SQRT( (TWO * JA + ONE) *                       & !
                         (TWO * JC + ONE) ) /                     & !
                   (TWO * JK + ONE)                                 !
      END DO                                                        !
!
      END SUBROUTINE CG_6J
!
!=======================================================================
!
      SUBROUTINE NINE_J(J1,J2,J12,J3,J4,J34,J13,J24,J,NINEJ)
!
!   This subroutine calculates Wigner's 9j coefficients using
!     the standard expansion in terms of products of 6j
!     coefficients (Varshalovich, Moskalev and Khersonskii p.340,
!     equation (20))
!
!   This formula is used here under the form:
!
!     { J1 J2 J12 }
!     { J3 J4 J34 } = sum_{X} (-1)**(2X) (2X+1) { J1  J    X } { J3 J24 X } { J24 J3 X }
!     { J13 J24 J }                             { J34 J2 J12 } { J2 J34 J4} { J1 J J13 }
!
!   Note real angular mometum values start by a J while the
!     corresponding integer index used for storage starts with
!     a L ---> LJ1 = INT(J1)
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO,ONE,TWO,HALF,SMALL
!
      IMPLICIT NONE
!
      INTEGER               ::  NC
      INTEGER               ::  LJ1J,LJ34J2
      INTEGER               ::  LJJ1_MIN,LJJ1_MAX
      INTEGER               ::  LJ3J24,LJ2J34
      INTEGER               ::  LJJ2_MIN,LJJ2_MAX
      INTEGER               ::  LJ24J3
      INTEGER               ::  LJJ3_MIN,LJJ3_MAX
      INTEGER               ::  LX,LX_MIN,LX_MAX
      INTEGER               ::  I_INT1,I_INT2,I_INT3
!
      INTEGER               ::  INT,ABS,MAX,MIN
!
      REAL (WP)             ::  J1,J2,J12,J3,J4,J34,J13,J24,J
      REAL (WP)             ::  NJ1(0:N_GAUNT),NJ2(0:N_GAUNT),NJ3(0:N_GAUNT)
      REAL (WP)             ::  NINEJ
      REAL (WP)             ::  J1J,J34J2
      REAL (WP)             ::  JJ1_MIN,JJ1_MAX
      REAL (WP)             ::  J3J24,J2J34
      REAL (WP)             ::  JJ2_MIN,JJ2_MAX
      REAL (WP)             ::  J24J3
      REAL (WP)             ::  JJ3_MIN,JJ3_MAX
      REAL (WP)             ::  X
      REAL (WP)             ::  SIG
!
      REAL (WP)             ::  SIGN,FLOAT
!
!   Call of the first 6j NJ1(J) with J varying between
!                 JJ1_MIN and JJ1_MAX
!
      NC = 6                                                        !
!
      CALL N_J(J1,J34,J,J2,J12,NJ1,I_INT1,NC)                       !
!
      J1J   = J1  - J                                               !
      J34J2 = J34 - J2                                              !
!
      LJ1J   = INT( J1J   + SIGN(SMALL,J1J)   )                     !
      LJ34J2 = INT( J34J2 + SIGN(SMALL,J34J2) )                     !
!
      JJ1_MIN  = MAX(ABS(J1J),ABS(J34J2))                           !
      JJ1_MAX  = MIN(J1+J,J34+J2)                                   !
      LJJ1_MIN = INT( JJ1_MIN + SIGN(SMALL,JJ1_MIN))                !
      LJJ1_MAX = INT( JJ1_MAX + SIGN(SMALL,JJ1_MAX))                !
!
!   Call of the second 6j NJ2(J) with J varying between
!                 JJ2_MIN and JJ2_MAX
!
      CALL N_J(J3,J2,J24,J34,J4,NJ2,I_INT2,NC)                      !
!
      J3J24 = J3 - J24                                              !
      J2J34 = J2 - J34                                              !
!
      LJ3J24 = INT( J3J24 + SIGN(SMALL,J3J24))                      !
      LJ2J34 = INT( J2J34 + SIGN(SMALL,J2J34))                      !
!
      JJ2_MIN  = MAX(ABS(J3J24),ABS(J2J34))                         !
      JJ2_MAX  = MIN(J3+J24,J2+J34)                                 !
      LJJ2_MIN = INT(JJ2_MIN + SIGN(SMALL,JJ2_MIN))                 !
      LJJ2_MAX = INT(JJ2_MAX + SIGN(SMALL,JJ2_MAX))                 !
!
!   Call of the third 6j NJ3(J) with J varying between
!                 JJ3_MIN and JJ3_MAX
!
      CALL N_J(J24,J1,J3,J,J13,NJ3,I_INT3,NC)                       !
!
      J24J3 = J24 - J3                                              !
      J1J   = J1  - J                                               !
!
      LJ24J3 = INT( J24J3 + SIGN(SMALL,J24J3))                      !
      LJ1J   = INT( J1J   + SIGN(SMALL,J1J)  )                      !
!
      JJ3_MIN  = MAX(ABS(J24J3),ABS(J1J))                           !
      JJ3_MAX  = MIN(J24+J3,J1+J)                                   !
      LJJ3_MIN = INT( JJ3_MIN + SIGN(SMALL,JJ3_MIN))                !
      LJJ3_MAX = INT( JJ3_MAX + SIGN(SMALL,JJ3_MAX))                !
!
!   General minimum value index LX for the loop index X
!
      LX_MIN = MAX(LJJ1_MIN,LJJ2_MIN,LJJ3_MIN)                      !
      LX_MAX = MIN(LJJ1_MAX,LJJ2_MAX,LJJ3_MAX)                      !
!
!   Calculation of the 9j symbol
!
      NINEJ = ZERO                                                  !
      DO LX = LX_MIN,LX_MAX                                         !
        IF(I_INT1 == 1) THEN                                        !
          X = FLOAT(LX)                                             !
        ELSE                                                        !
          X = FLOAT(LX) + HALF                                      !
        END IF                                                      !
        SIG   = (-ONE)**(TWO*X)                                     !
        NINEJ = NINEJ + SIG * (TWO*X+ONE) * NJ1(LX) *             & !
                                            NJ2(LX) * NJ3(LX)       !
      END DO                                                        !
!
      END SUBROUTINE NINE_J
!
!=======================================================================
!
      SUBROUTINE TENS_PROD(F1,NL1,F2,NL2,F)
!
!  This subroutine computes the tensor product of two
!    angular momentum dependent functions F1 and F2
!    according to :
!
!     F_{L,M} = sum_{M1,M2} CG(L1,M1,L2,M2|L,M) * F1_{L1,M1} * F2_{L2,M2}
!
!    where CG(L1,M1,L2,M2|L,M) is a Clebsch-Gordan coefficient
!
!  Note that it is written for integer values of the angular momenta.
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO
      USE COMPLEX_NUMBERS,     ONLY : ZEROC
!
      IMPLICIT NONE
!
      INTEGER               ::  NL1,NL2
      INTEGER               ::  N
      INTEGER               ::  L1,L2,L
      INTEGER               ::  LMIN,LMAX
      INTEGER               ::  M1,M2,M
      INTEGER               ::  I_INT
!
      INTEGER               ::  ABS
!
      REAL (WP)             ::  J1,MJ1,J2,MJ2,MJ
      REAL (WP)             ::  CG(0:N_GAUNT)
!
      REAL (WP)             ::  FLOAT
!
      COMPLEX (WP)          ::  F_1,F_2
      COMPLEX (WP)          ::  F1(0:L_MAX,-L_MAX:L_MAX)
      COMPLEX (WP)          ::  F2(0:L_MAX,-L_MAX:L_MAX)
      COMPLEX (WP)          ::  F(0:L_MAX,0:L_MAX,0:L_MAX,-L_MAX:L_MAX)
!
      MJ = ZERO                                                     !
      N  = 2                                                        !
!
      DO L1 = 0,NL1                                                 !
        J1 = FLOAT(L1)                                              !
        DO L2 = 0,NL2                                               !
          J2    = FLOAT(L2)                                         !
          LMIN = ABS(L1-L2)                                         !
          LMAX = L1+L2                                              !
!
          DO L = LMIN,LMAX                                          !
            DO M=-L,L                                               !
            F(L1,L2,L,M) = ZEROC                                    !
            END DO                                                  !
            DO M1 = -L1,L1                                          !
              F_1 = F1(L1,M1)                                       !
              MJ1 = FLOAT(M1)                                       !
              DO M2 = -L2,L2                                        !
                F_2 = F2(L2,M2)                                     !
                MJ2 = FLOAT(M2)                                     !
                CALL N_J(J1,MJ1,J2,MJ2,MJ,CG,I_INT,N)               !
                M = M1 + M2                                         !
                F(L1,L2,L,M) = F(L1,L2,L,M) + F_1 * F_2 * CG(L)     !
              END DO                                                !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
!
      END SUBROUTINE TENS_PROD
!
END MODULE ANGULAR_MOMENTUM

