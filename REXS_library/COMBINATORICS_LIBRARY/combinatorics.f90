!
!=======================================================================
!
MODULE CALC_LOGAMMA
!
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,      ONLY : L_MAX,N_GAUNT
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE CALC_LOG_GAMMA
!
!  This subroutine computes the logarithm of the Gamma function
!    and stores the result in a common block:
!
!   GLD contains the logarithm of the Gamma function
!              for integers and semi-integers
!
!   GLD contains the logarithm of the Gamma function
!              for integers only
!
!  Author : D. Sébilleau
!
!                                          Last modified : 21 May 2021
!
      USE LOGAMMA
      USE REAL_NUMBERS,        ONLY : ZERO,TWO,HALF
      USE PI_ETC,              ONLY : PI
!
      IMPLICIT NONE
!
      INTEGER               ::  I,J
!
      REAL (WP)             ::  LOG,SQRT,FLOAT
!
!
!  Storage of the logarithm of the Gamma function GLD(N+1,N_INT)
!  for integer (N_INT=1) and half-integer (N_INT=2) values :
!
!    GLD(N+1,1)   =   Log(N!) for N integer
!    GLD(N+1/2,2) =   Log(N!) for N half-integer
!
      DO I = 0, N_GAUNT                                             !
        GLG(I)   = ZERO                                             !
      END DO                                                        !
      DO I = 0, L_MAX                                               !
        GLD(I,1) = ZERO                                             !
        GLD(I,2) = ZERO                                             !
      END DO                                                        !
!
      GLD(1,2) = LOG(SQRT(PI)/TWO)                                  !
!
      DO I = 2, N_GAUNT                                             !
        J        = I - 1                                            !
        GLG(I)   = GLG(J)   + LOG(FLOAT(J))                         !
      END DO                                                        !
      DO I = 2, L_MAX                                               !
        J        = I - 1                                            !
        GLD(I,1) = GLD(J,1) + LOG(FLOAT(J))                         !
        GLD(I,2) = GLD(J,2) + LOG(FLOAT(J) + HALF)                  !
      END DO                                                        !
!
      END SUBROUTINE CALC_LOG_GAMMA
!
END MODULE CALC_LOGAMMA
!
!=======================================================================
!
MODULE COMBINATORICS
!
!  This module provides various combinatorics functions
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,      ONLY : L_MAX
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE LOG_GAMMA(NMAX,LG)
!
!  This subroutine computes the logarithm of the Gamma function for
!           integer values (i.e. Log(n!))
!
!
!  Input variables :
!
!                       NMAX      :  upper value of n
!
!  Output variables :
!
!                       LG        :  array containing Log(n!)
!
!
!  Author : D. Sébilleau
!
!                                          Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO
!
      IMPLICIT NONE
!
      REAL (WP)             ::  LG(0:L_MAX)
!
      REAL (WP)             ::  FLOAT,LOG
!
      INTEGER               ::  NMAX,I,J
!
!  Initialization to zero
!
      DO I = 0, L_MAX                                               !
        LG(I) = ZERO                                                !
      END DO                                                        !
!
      DO I = 1, NMAX                                                !
        J = I - 1                                                   !
        LG(I) = LG(J) + LOG(FLOAT(J))                               !
      END DO                                                        !
!
      END SUBROUTINE LOG_GAMMA
!
!=======================================================================
!
      FUNCTION BIN(N,K)
!
!  Tnis function computes the binomial coefficient
!
!           C_n^k = ( n ) = n! / ( k! (n-k)! )
!                   ( k )
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 21 May 2021
!
!
      IMPLICIT NONE
!
      INTEGER               ::  N,K
!
      REAL (WP)             ::  LG(0:L_MAX)
      REAL (WP)             ::  BIN
!
      REAL (WP)             ::  EXP
!
!  Computing the logarithm of the Gamma function
!
      CALL LOG_GAMMA(N,LG)                                          !
!
      BIN = EXP(LG(N) - LG(K) - LG(N-K))                            !
!
      END FUNCTION BIN
!
!=======================================================================
!
      FUNCTION NIB(N,K)
!
!  Tnis function computes the coefficient
!
!           (n+k)! / ( k! (n-k)! )
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 21 May 2021
!
!
      IMPLICIT NONE
!
      INTEGER               ::  N,K
!
      REAL (WP)             ::  LG(0:L_MAX)
      REAL (WP)             ::  NIB
!
      REAL (WP)             ::  EXP
!
!  Computing the logarithm of the Gamma function
!
      CALL LOG_GAMMA(N+K,LG)                                        !
!
      NIB = EXP(LG(N+K) - LG(K) - LG(N-K))                          !
!
      END FUNCTION NIB
!
!=======================================================================
!
      SUBROUTINE COMBINATORIAL(NMAX,NUMBER,CN)
!
!  This subroutine computes numbers resulting from combinatorics
!    This version if for integers only.
!
!
!
!  Input variables :
!
!                       NMAX      :  upper value of n
!                       NUMBER    :  type of numbers computed
!                                       ---> 'BINOMIAL  ' : binomial coefficients
!                                       ---> 'POCHHAMMER' : Pochhammer coefficients
!                                       ---> 'STIRLING1S' : signed Stirling numbers of 1st kind
!                                       ---> 'STIRLING1U' : unsigned Stirling numbers of 1st kind
!                                       ---> 'STIRLING2N' : Stirling numbers of 2nd kind
!
!  Output variables :
!
!
!                       CN        :  resulting numbers
!
!
!  Author : D. Sébilleau
!
!                                          Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO,ONE
!
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 10)  ::  NUMBER
!
      INTEGER               ::  NMAX
      INTEGER               ::  I,J,N,K
!
      REAL (WP)             ::  X
      REAL (WP)             ::  CN(0:NMAX,0:NMAX)
      REAL (WP)             ::  LG(0:L_MAX)
!
      REAL (WP)             ::  EXP,FLOAT
!
!  Initialization of the array
!
      DO I = 0, NMAX                                                !
        DO J = 0, NMAX                                              !
          CN(I,J) = ZERO                                            !
        END DO                                                      !
      END DO                                                        !
!
      IF(NUMBER == 'BINOMIAL  ') THEN                               !  ( N )
!                                                                   !  ( K )
        CALL LOG_GAMMA(NMAX,LG)                                     !
!
        CN(0,0) = ONE                                               !
        DO N = 1, NMAX                                              !
          DO K = 1, NMAX-N                                          !
            X       = LG(N) - LG(K) - LG(N-K)                       !
            CN(N,K) = EXP(X)                                        !
          END DO                                                    !
        END DO                                                      !
!
      ELSE IF(NUMBER == 'POCHHAMMER') THEN                          !  (N)_K
!
        CALL LOG_GAMMA(NMAX,LG)                                     !
!
        CN(0,0) = ONE                                               !
        DO N = 1, NMAX                                              !
          DO K = 1, NMAX-N                                          !
            X       = LG(N+K) - LG(N)                               !
            CN(N,K) = EXP(X)                                        !
          END DO                                                    !
        END DO                                                      !
!
      ELSE IF(NUMBER == 'STIRLING1U') THEN                          !  c(N,K)
!
        CN(0,0)    = ONE                                            !
        CN(NMAX,0) = ZERO                                           !
!
        DO N = 1, NMAX-1                                            !
          CN(N,0) = ZERO                                            !
          DO K=1, NMAX-N+1                                          !
            CN(N+1,K) = FLOAT(N) * CN(N,K) + CN(N,K-1)              !
          END DO                                                    !
        END DO                                                      !
!
      ELSE IF(NUMBER == 'STIRLING1S') THEN                          !  s(N,K)
!
        CN(0,0)    = ONE                                            !
        CN(NMAX,0) = ZERO                                           !
!
        DO N = 1, NMAX-1                                            !
          CN(N,0) = ZERO                                            !
          DO K = 1, NMAX-N+1                                        !
            CN(N+1,K) = - FLOAT(N) * CN(N,K) + CN(N,K-1)            !
          END DO                                                    !
        END DO                                                      !
!
      ELSE IF(NUMBER == 'STIRLING2N') THEN                          !  S(N,K)
!
        CN(0,0)    = ONE                                            !
        CN(NMAX,0) = ZERO                                           !
!
        DO N = 1, NMAX-1                                            !
          CN(N,0) = ZERO                                            !
          DO K =1, NMAX-N+1                                         !
            CN(N+1,K) = FLOAT(K) * CN(N,K) + CN(N,K-1)              !
          END DO                                                    !
        END DO                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE COMBINATORIAL
!
!=======================================================================
!
      FUNCTION GAMMLN(XX)
!
!     Logarithm of Gamma function
!
!
!  Adapted from the Fortran 77 version in:
!
!    "Numerical Recipes : The Art of Scientific Computing"
!  by W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vetterling
!               (Cambridge University Press 1992)
!
!
!                                     Last modified (DS) : 21 May 2021
!
!
      USE REAL_NUMBERS,     ONLY : ONE,HALF
!
      IMPLICIT NONE
!
      INTEGER                ::  J
!
      REAL (WP), INTENT(IN)  ::  XX
      REAL (WP)              ::  GAMMLN
!
      REAL (WP)              ::  SER,STP,TMP,X,Y
      REAL (WP)              ::  COF(6)
!
      REAL (WP)              ::  LOG
!
      SAVE COF,STP
!
      DATA COF / 76.18009172947146E0_WP, -86.50532032941677E0_WP, &
                 24.01409824083091E0_WP, -1.231739572450155E0_WP, &
                .1208650973866179E-2_WP,  -.5395239384953E-5_WP   /
      DATA STP /2.5066282746310005E0_WP                           /
!
      X   = XX                                                      !
      Y   = X                                                       !
      TMP = X + 5.5E0_WP                                            !
      TMP = (X + HALF) * LOG(TMP) - TMP                             !
      SER = 1.000000000190015E0_WP                                  !
!
      DO J = 1, 6                                                   !
        Y   = Y + ONE                                               !
        SER = SER + COF(J) / Y                                      !
      END DO                                                        !
!
      GAMMLN = TMP + LOG(STP * SER / X)                             !
!
      END FUNCTION GAMMLN
!
!=======================================================================
!
      SUBROUTINE GAMMA(X,GA)
!
!*****************************************************************************80
!
!  GAMMA evaluates the Gamma function.
!
!  Licensing:
!
!    The original FORTRAN77 version of this routine is copyrighted by
!    Shanjie Zhang and Jianming Jin.  However, they give permission to
!    incorporate this routine into a user program that the copyright
!    is acknowledged.
!
!  Modified:
!
!    08 September 2007
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!    X must not be 0, or any negative integer.
!
!    Output, real ( kind = 8 ) GA, the value of the Gamma function.
!
!
!                                     Last modified (DS) : 20 Aug 2020
!
!
      USE REAL_NUMBERS,     ONLY : ZERO,ONE
      USE PI_ETC,           ONLY : PI
!
      IMPLICIT NONE
!
      REAL (WP), DIMENSION(26) :: G = (/                        &
                                     1.0E+00_WP               , &
                                     0.5772156649015329E+00_WP, &
                                    -0.6558780715202538E+00_WP, &
                                    -0.420026350340952E-01_WP , &
                                     0.1665386113822915E+00_WP, &
                                    -0.421977345555443E-01_WP , &
                                    -0.96219715278770E-02_WP  , &
                                     0.72189432466630E-02_WP  , &
                                    -0.11651675918591E-02_WP  , &
                                    -0.2152416741149E-03_WP   , &
                                     0.1280502823882E-03_WP   , &
                                    -0.201348547807E-04_WP    , &
                                    -0.12504934821E-05_WP     , &
                                     0.11330272320E-05_WP     , &
                                    -0.2056338417E-06_WP      , &
                                     0.61160950E-08_WP        , &
                                     0.50020075E-08_WP        , &
                                    -0.11812746E-08_WP        , &
                                     0.1043427E-09_WP         , &
                                     0.77823E-11_WP           , &
                                    -0.36968E-11_WP           , &
                                     0.51E-12_WP              , &
                                    -0.206E-13_WP             , &
                                    -0.54E-14_WP              , &
                                     0.14E-14_WP              , &
                                     0.1E-15_WP                 &
                                      /)
!
      INTEGER                   ::  K,M,M1
!
      INTEGER                   ::  INT
!
      REAL (WP), INTENT(IN)     ::  X
!
      REAL (WP), INTENT(OUT)    ::  GA
!
      REAL (WP)                 ::  GR,R,Z
!
      REAL (WP)                 ::  ABS,SIN,REAL
!
      IF(X == AINT(X)) THEN                                         !
!
        IF(ZERO < X) THEN                                           !
          GA = ONE                                                  !
          M1 = INT(X) - 1                                           !
          DO K = 2, M1                                              !
            GA = GA * K                                             !
          END DO                                                    !
        ELSE                                                        !
          GA = 1.0E+300_WP                                          !
        END IF                                                      !
!
      ELSE                                                          !
!
        IF(ONE < ABS(X)) THEN                                       !
          Z = ABS(X)                                                !
          M = INT(Z)                                                !
          R = ONE                                                   !
          DO K = 1, M                                               !
            R = R * (Z - REAL(K,KIND = 8))                          !
          END DO                                                    !
          Z = Z - REAL(M,KIND = 8)                                  !
        ELSE                                                        !
          Z = X                                                     !
        END IF                                                      !
!
        GR = G(26)                                                  !
        DO K = 25, 1, -1                                            !
          GR = GR * Z + G(K)                                        !
        END DO                                                      !
!
        GA = ONE / (GR * Z)                                         !

        IF (ONE < ABS(X)) THEN                                      !
          GA = GA * R                                               !
          IF(X < ZERO) THEN                                         !
            GA = - PI / (X* GA * SIN(PI * X))                       !
          END IF                                                    !
        END IF                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE GAMMA
!
END MODULE COMBINATORICS
