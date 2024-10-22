!
!=======================================================================
!
MODULE LOCATE_MOD 
!
!  This module provides several routines for locating 
!    within an ordered array
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE LOCATE(XX,N,X,J)
!
!
!            This subroutine is taken from the book :
!          "Numerical Recipes : The Art of Scientific
!           Computing" par W.H. PRESS, B.P. FLANNERY,
!               S.A. TEUKOLSKY et W.T. VETTERLING 
!               (Cambridge University Press 1992)
!
!    It performs a search in an ordered table using a bisection method.
!    Given a monotonic array XX(1:N) and a value X, it returns J such 
!                  that X is between XX(J) and XX(J+1). 
!
      INTEGER, INTENT(IN)   ::  N
      INTEGER, INTENT(OUT)  ::  J
!
      INTEGER               ::  JL,JM,JU
!
      REAL (WP), INTENT(IN) ::  XX(N),X
!
      JL = 0                                                        !
      JU = N + 1                                                    !
!
  10  IF(JU-JL > 1)THEN                                             !
        JM = (JU+JL) / 2                                            !
        IF((XX(N) > XX(1)) .EQV. (X > XX(JM)))THEN                  !
          JL = JM                                                   !
        ELSE
          JU = JM                                                   !
        END IF                                                      !
      GO TO 10                                                      !
      END IF                                                        !
      J = JL                                                        !
!
      END SUBROUTINE LOCATE
!
END MODULE LOCATE_MOD
