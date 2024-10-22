!
!=======================================================================
!
MODULE FACTORIALS
!
!  This module provides factorials and other related numbers
!
!
      USE ACCURACY_REAL
      USE REAL_NUMBERS,             ONLY : ZERO,ONE
!
CONTAINS
!
!=======================================================================
!
      FUNCTION FAC(N)
!
!  This function computes the factorial of n
!
!  Input parameters:
!
!       * N        : integer
!
!  Output variables :
!
!       * FAC      : n!
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 21 May 2021
!
      IMPLICIT NONE
!
      REAL (WP)             ::  FAC
      REAL (WP)             ::  FACT(50)
!
      REAL (WP)             ::  FLOAT
!
      INTEGER               ::  N,K
      INTEGER               ::  LOGF
!
      LOGF=6                                                        !
!
      IF(N > 50) THEN                                               !
        WRITE(LOGF,10)                                              !
        STOP                                                        !
      END IF                                                        !
!
      FACT(1) = ONE                                                 !
!
      DO K = 2, N                                                   !
        FACT(K) = FACT(K-1) * FLOAT(K)                              !
      END DO                                                        !
!
      FAC = FACT(N)                                                 !
!
  10  FORMAT(5X,'<<<<<  DIMENSION ERROR IN FAC FUNCTION  >>>>>',/, &!
             5X,'<<<<<  N SHOULD BE <= 50 OR REDIMENSION >>>>>',//) !
!
      END FUNCTION FAC
!
!=======================================================================
!
      FUNCTION FACTLN(N)
!
!     Logarithm of factorial function
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
      USE REAL_NUMBERS,     ONLY : ZERO,ONE
      USE COMBINATORICS,    ONLY : GAMMLN
!
      IMPLICIT NONE
!
      INTEGER, PARAMETER     ::  NMAX = 500
!
      INTEGER, INTENT(IN)    ::  N
!
      INTEGER                ::  J
!
      REAL (WP)              ::  FACTLN
!
      REAL (WP)              ::  A(NMAX)
!
      REAL (WP)             ::  FLOAT
!
      IF (N < 0) THEN                                               !
        WRITE(6,10) N                                               !
        STOP                                                        !
      END IF                                                        !
!
!  A initialized to -1
!
      DO J = 1,NMAX                                                 !
        A(J)= - ONE                                                 !
      END DO                                                        !
!
      IF(N <= NMAX-1) THEN                                          !
        IF(A(N+1) < ZERO) A(N+1) = GAMMLN(FLOAT(N) + ONE)           !
        FACTLN = A(N+1)                                             !
      ELSE                                                          !
        FACTLN = GAMMLN(FLOAT(N) + ONE)                             !
      END IF                                                        !
!
!  Formats:
!
  10  FORMAT('Negative factorial in FACTLN : N = ',I5)
!
      END FUNCTION FACTLN
!
END MODULE FACTORIALS
