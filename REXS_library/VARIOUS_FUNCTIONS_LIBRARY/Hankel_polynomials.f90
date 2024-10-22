!
!=======================================================================
!
MODULE HANKEL_POLYNOMIALS
!
!  This module provides the Hankel polynomials
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE POLHAN(I_BASIS,NO,NC,RHO,HLM)
!
!  This routine calculates a function HLM(L,M), related to the the Hankel
!      polynomials and their derivatives with respect to z=1/ikr,
!      necessary for the Rehr-Albers expansion of the propagator.
!
!
!  Input variables:
!
!     I_BASIS   :  approximation index
!                    I_BASIS = 1  --> Hankel polynomials set to 1.0
!                    I_BASIS > 1  --> standard calculation
!     NO        : Rehr-Albers index
!     NC        : maximum l value computed
!     RHO       : current value of k * r
!
!  Output variable:
!
!     HLM       :  Hankel polynomial array
!                    stored as HLM(M,L)
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 17 May 2021
!
!
      USE COMPLEX_NUMBERS,        ONLY : ONEC,IC
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)        ::  I_BASIS,NO,NC
!
      INTEGER                    ::  L,M
!
      COMPLEX (WP), INTENT(IN )  ::  RHO
      COMPLEX (WP), INTENT(OUT)  ::  HLM(0:NO_ST_M,0:NL_M-1)
!
      REAL    (WP)               ::  FLOAT
!
      COMPLEX (WP)               ::  Z
!
      IF(I_BASIS >= 1) THEN                                         !
!
        Z = - IC / RHO                                              !
!
!  Case M = 0
!
        HLM(0,0) = ONEC                                             !
        HLM(0,1) = ONEC - Z                                         !
        DO L = 2, NC                                                !
          HLM(0,L) = HLM(0,L-2) - FLOAT(L + L - 1) * Z * HLM(0,L-1) !
        END DO                                                      !
!
!  Case M > 0
!
        IF(NO >= 1) THEN                                            !
          DO M = 1, NO                                              !
            HLM(M,M)   = - Z * HLM(M-1,M-1) * FLOAT(M + M - 1)      !
            HLM(M,M+1) = HLM(M,M) * FLOAT(M + M + 1) *            & !
                         ( ONEC - Z * FLOAT(M + 1) )                !
            DO L = M + 2, NC                                        !
              HLM(M,L) = HLM(M,L-2) - FLOAT(L + L - 1) * Z *      & !
                         (HLM(M,L-1) + HLM(M-1,L-1))                !
            END DO                                                  !
          END DO                                                    !
        END IF                                                      !
!
      ELSE                                                          !
!
        DO M = 0, NO                                                !
          DO L = M, NC                                              !
            HLM(M,L) = ONEC                                         !
          END DO                                                    !
        END DO                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE POLHAN
!
END MODULE HANKEL_POLYNOMIALS
