MODULE RENORMALIZATION
!
!  This module computes the renormalization coefficients
!
      USE ACCURACY_REAL
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE COEF_RENORM(N_SCAT)
!
!  This subroutine computes the coefficients for the renormalization
!     of the multiple scattering series. These coefficients are
!     expressed as C_REN(K) where K is the multiple scattering order.
!     REN2 is the value of the mixing (or renormalization) parameter.
!
!     N_SCAT is the scattering order at which the series is truncated,
!     so that K varies from 0 to N_SCAT.
!
!  MODULE CURRENT_RENORM :
!
!     I_REN = 1 : renormalization in terms of G_n matrices (n : N_REN)
!           = 2 : renormalization in terms of the Sigma_n matrices
!           = 3 : renormalization in terms of the Z_n matrices
!           = 4 : Löwdin renormalization in terms of the Pi_1 matrices
!           = 5 : Löwdin renormalization in terms of the L_n  matrices
!
!     N_REN = renormalization order n
!
!     REN   = REN_R + i * REN_I : omega
!
!
!
!  Authors : D. Sébilleau, A. Takatsu, M. Terao-Dunseath, K. Dunseath
!
!
!                                       Last modified (DS): 14 Jun 2021
!
!
      USE DIMENSION_CODE,      ONLY : N_SCAT_M
!
      USE REAL_NUMBERS,        ONLY : ZERO,ONE
      USE COMPLEX_NUMBERS
!
      USE CURRENT_RENORM
!
      USE CURRENT_COEF_RENORM
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)    ::  N_SCAT
!
      INTEGER                ::  J,K,M,N
      INTEGER                ::  IN,N_SCAT2
!
      INTEGER                ::  INT
!
      REAL (WP)              ::  C(0:N_SCAT_M,0:N_SCAT_M)
!
      REAL (WP)              ::  FLOAT
!
      COMPLEX (WP)           ::  REN,REN2,COEF1,COEF2
      COMPLEX (WP)           ::  SUM_L,POWER
      COMPLEX (WP)           ::  X(0:N_SCAT_M,0:N_SCAT_M)
      COMPLEX (WP)           ::  Y1(0:N_SCAT_M,0:N_SCAT_M)
!
      REN = REN_R + IC * REN_I                                      ! omega
!
!  Initialisation of renormalization coefficients
!
      DO J = 0, N_SCAT                                              !
        C_REN(J) = ZEROC                                            !
      END DO                                                        !
!
!  Computing the binomial coefficients C(N,K) = (N) = N! / K! (N-K)!
!                                               (K)
      C(0,0) = ONE                                                  !
      C(1,0) = ONE                                                  !
      C(1,1) = ONE                                                  !
      DO N = 2, N_SCAT                                              !
        C(N,0) = ONE                                                !
        C(N,N) = ONE                                                !
        DO K = 1, N-1                                               !
          C(N,K) = C(N-1,K) +C (N-1,K-1)                            !
        END DO                                                      !
      END DO                                                        !
!
      IF(I_REN <= 3) THEN                                           !
!
!  Computing the modified renormalization parameter REN2 (g_n,s_n,zeta_n)
!
        IF(I_REN == 1) THEN                                         !
!
!.....(g_n,G_n) renormalization
!
          REN2 = REN**N_REN                                         ! g_n = omega^n
!
        ELSE IF(I_REN == 2) THEN                                    !
!
!.....(s_{n},Sigma_n) renormalization
!
          REN2 = ( ONEC - REN**(N_REN+1) ) /                      & !
                 ( FLOAT(N_REN+1) * (ONEC-REN) )                    ! s_n
!
        ELSE IF(I_REN == 3) THEN                                    !
!
!.....(zeta_{n},Z_n) renormalization
!
           REN2 = - (ONEC - REN)**(N_REN+1)                         ! zeta_n
!
        END IF                                                      !
!
        DO K=0,N_SCAT                                               !
          C_REN(K) = ZEROC                                          !
          DO J = K, N_SCAT                                          !
            C_REN(K) = C_REN(K) + C(J,K) * (ONEC - REN2)**(J-K)     !
          END DO                                                    !
          C_REN(K) = C_REN(K) * REN2**(K+1)                         !
        END DO                                                      !
!
        ELSE IF(I_REN == 4) THEN                                    !
!
!     Loewdin (Pi_1) renormalization for n = 1
!
!     Notation: Y1(K,M) : [Y_1^k]_m
!
        N_SCAT2 = (N_SCAT - 1) / 2                                  !
!
        COEF1 = ONEC - REN                                          ! (1 - omega)
!
        Y1      = ZEROC                                             !
        Y1(0,0) = ONEC                                              !
        Y1(0,1) = REN                                               !
!
        DO K = 1, N_SCAT2                                           !
          Y1(K,K) = COEF1**K                                        !
            DO M = K+1, 2*K                                         !
              COEF2   = (REN**(M-K)) * (COEF1**(2*K-M))             !
              Y1(K,M) = COEF2 * ( C(K,M-K) + COEF1 * C(K,M-K-1) )   !
            END DO                                                  !
            Y1(K,2*K+1) = REN**(K+1)                                !
        END DO
!                                                                   !
        C_REN(0) = ONEC                                             !
        C_REN(1) = ONEC                                             !
        DO K = 2, N_SCAT                                            !
           IN = INT(K / 2)                                          !
           C_REN(K) = ZEROC                                         !
           DO M = IN, N_SCAT2                                       !
              C_REN(K) = C_REN(K) + Y1(M,K)                         !
           END DO                                                   !
        END DO                                                      !
!
        ELSE IF(I_REN == 5) THEN
!
!     Loewdin L_n(omega,N_SCAT) renormalization
!
!     Notation: X(K,N) = X_n(omega,k)
!
!
!  Computing the X(N,K) coefficients, with K <= N
!
        POWER = ONEC / REN                                          !
        DO N = 0, N_SCAT                                            !
          POWER = POWER * REN                                       ! omega^n
          IF(N == 0) THEN                                           !
            X(N,0) = ONEC                                           !
          ELSE                                                      !
            X(N,0) = ZEROC                                          !
          END IF                                                    !
          DO K = 1, N_SCAT                                          !
            IF(K > N) THEN                                          !
              X(N,K) = ZEROC                                        !
            ELSE IF(K == N) THEN                                    !
              X(N,K) = POWER * X(N-1,K-1)                           !
            ELSE                                                    !
              X(N,K) = X(N-1,K) * (REN - POWER) + POWER * X(N-1,K-1)!
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
!
!  Calculation of L_n(omega,NDIF)
!
        DO N = 0, N_SCAT                                            !
          SUM_L = ZEROC                                             !
          DO K = N, N_SCAT                                          !
            SUM_L = SUM_L + X(K,N)                                  !
          END DO                                                    !
          C_REN(N) = SUM_L                                          !
        END DO                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE COEF_RENORM
!
END MODULE RENORMALIZATION
