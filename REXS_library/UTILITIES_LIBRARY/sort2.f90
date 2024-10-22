!
!=======================================================================
!
MODULE SORT2
!
!  This module provides several routines for sorting INTEGER arrays
!
!
!   1) SUBROUTINE SORT_INC_I(NINI,VALIN,NFIN,VALFIN)
!
!   2) SUBROUTINE SORT_DEC_I(NINI,VALIN,NFIN,VALFIN)
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE SORT_INC_I(NINI,VALIN,NFIN,VALFIN)
!
!   Given a set of **integer** numbers VALINI, this routine orders them and
!       suppresses the values appearing more than once. The remaining
!       values are stored in VALFIN.
!
!       VALINI(K+1) > VALINI(K) : decreasing order
!       VALINI(K+1) < VALINI(K) : increasing order
!
!
!                                     Last modified (DS) : 18 Jun 2021
!
!
      USE REAL_NUMBERS,     ONLY : SMALL
!
      IMPLICIT NONE
!
      INTEGER, INTENT (IN)   ::  NINI
      INTEGER, INTENT (OUT)  ::  NFIN
      INTEGER, INTENT(IN)    ::  VALIN(NINI)
      INTEGER, INTENT(OUT)   ::  VALFIN(NINI)
!
      INTEGER                ::  I,J,K
      INTEGER                ::  JFIN
      INTEGER                ::  VALINI(NINI)
      INTEGER                ::  R1
!
      LOGICAL                ::  BUBBLE
!
!...Store input array
!
      DO I = 1,NINI                                                 !
        VALINI(I) = VALIN(I)                                        !
      END DO                                                        !
!
      DO J = 1,NINI-1                                               !
        K      = J                                                  !
        BUBBLE = .TRUE.                                             !
!
  150   IF(K >= 1 .AND. BUBBLE) THEN                                !
          IF(VALINI(K+1) < VALINI(K)) THEN                          !
            R1          = VALINI(K)                                 !
            VALINI(K)   = VALINI(K+1)                               !
            VALINI(K+1) = R1                                        !
          ELSE                                                      !
            BUBBLE = .FALSE.                                        !
          END IF                                                    !
          K = K - 1                                                 !
          GO TO 150                                                 !
        END IF                                                      !
      END DO                                                        !
!
      JFIN      = 1                                                 !
      VALFIN(1) = VALINI(1)                                         !
      DO J = 1,NINI-1                                               !
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)) > SMALL) THEN              !
          JFIN         = JFIN+1                                     !
          VALFIN(JFIN) = VALINI(J+1)                                !
        END IF                                                      !
      END DO                                                        !
      NFIN = JFIN                                                   !
!
      RETURN                                                        !
!
      END SUBROUTINE SORT_INC_I
!
!=======================================================================
!
      SUBROUTINE SORT_DEC_I(NINI,VALIN,NFIN,VALFIN)
!
!   Given a set of **integer** numbers VALINI, this routine orders them and
!       suppresses the values appearing more than once. The remaining
!       values are stored in VALFIN.
!
!       VALINI(K+1) > VALINI(K) : decreasing order
!       VALINI(K+1) < VALINI(K) : increasing order
!
!
!                                     Last modified (DS) : 18 Jun 2021
!
!
      USE REAL_NUMBERS,     ONLY : SMALL
!
      IMPLICIT NONE
!
      INTEGER, INTENT (IN)   ::  NINI
      INTEGER, INTENT (OUT)  ::  NFIN
      INTEGER, INTENT(IN)    ::  VALIN(NINI)
      INTEGER, INTENT(OUT)   ::  VALFIN(NINI)
!
      INTEGER                ::  I,J,K
      INTEGER                ::  JFIN
      INTEGER                ::  VALINI(NINI)
      INTEGER                ::  R1
!
      LOGICAL                ::  BUBBLE
!
!...Store input array
!
      DO I = 1,NINI                                                 !
        VALINI(I) = VALIN(I)                                        !
      END DO                                                        !
!
      DO J = 1,NINI-1                                               !
        K      = J                                                  !
        BUBBLE = .TRUE.                                             !
!
  150   IF(K >= 1 .AND. BUBBLE) THEN                                !
          IF(VALINI(K+1) > VALINI(K)) THEN                          !
            R1          = VALINI(K)                                 !
            VALINI(K)   = VALINI(K+1)                               !
            VALINI(K+1) = R1                                        !
          ELSE                                                      !
            BUBBLE = .FALSE.                                        !
          END IF                                                    !
          K = K - 1                                                 !
          GO TO 150                                                 !
        END IF                                                      !
      END DO                                                        !
!
      JFIN      = 1                                                 !
      VALFIN(1) = VALINI(1)                                         !
      DO J = 1,NINI-1                                               !
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)) > SMALL) THEN              !
          JFIN         = JFIN+1                                     !
          VALFIN(JFIN) = VALINI(J+1)                                !
        END IF                                                      !
      END DO                                                        !
      NFIN = JFIN                                                   !
!
      RETURN                                                        !
!
      END SUBROUTINE SORT_DEC_I
!
END MODULE SORT2

