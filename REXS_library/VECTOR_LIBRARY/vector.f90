!
!=======================================================================
!
MODULE VECTOR
!
!  This module contains vector operation routines
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      FUNCTION PRSCAL(A1,A2)
!
!  This function computes the dot product of the two vectors A1 and A2
!
!
!                                           Last modified:  7 May 2021
!
!
      IMPLICIT NONE
!
      REAL (WP), INTENT(IN)    ::  A1(3),A2(3)
!
      REAL (WP)                ::  PRSCAL
!
      PRSCAL = A1(1) * A2(1) + A1(2) * A2(2) + A1(3)*  A2(3)        !
!
      END FUNCTION PRSCAL
!
!=======================================================================
!
      SUBROUTINE PRVECT(A1,A2,A3,C)
!
!  This function computes the vector product of the two vectors A1 and A2.
!            The result is A3; C is a scaling factor
!
!
!                                           Last modified:  7 May 2021
!
!
      IMPLICIT NONE
!
      REAL (WP), INTENT(IN)    ::  A1(3),A2(3),C
      REAL (WP), INTENT(OUT)   ::  A3(3)
!
      A3(1) = (A1(2) * A2(3) - A1(3) * A2(2)) / C                   !
      A3(2) = (A1(3) * A2(1) - A1(1) * A2(3)) / C                   !
      A3(3) = (A1(1) * A2(2) - A1(2) * A2(1)) / C                   !
!
      END SUBROUTINE PRVECT
!
!=======================================================================
!
      SUBROUTINE XYZ_TO_R_THETA_PHI(X,Y,Z,R,TH,PH)
!
!  This subroutine computes the r,theta,phi of a vector from
!    the knowledge ot is cartesian coordinates
!
!
!  Input variables :
!
!                       X         :  \
!                       Y         :   >  Cartesian coordinates of vector
!                       Z         :  /
!
!
!  Output variables :
!
!                       R         :  \
!                       TH        :   >  r, theta and phi
!                       PH        :  /
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified :  2 Jun 2021
!
!
      USE REAL_NUMBERS,          ONLY : ZERO,ONE,HALF,SMALL
      USE COMPLEX_NUMBERS,       ONLY : IC
      USE PI_ETC,                ONLY : PI
!
      USE ARC_SIN
!
      IMPLICIT NONE
!
      REAL (WP), INTENT(IN)      ::  X,Y,Z
      REAL (WP), INTENT(OUT)     ::  R,TH,PH
!
      REAL (WP)                  ::  COS_TH
!
      REAL (WP)                  ::  ACOS,SQRT
!
      COMPLEX (WP)               ::  XPIY
!
      R = SQRT( X * X + Y * Y + Z * Z )                             !
!
      IF(R > SMALL) THEN                                            !
        COS_TH = Z / R                                              !
      ELSE                                                          !
        TH = HALF * PI                                              !
        PH = ZERO                                                   !
        GO TO 10                                                    !
      END IF
!
      IF(COS_TH > ONE)  THEN                                        !
         COS_TH =   ONE                                             !
      ELSE IF(COS_TH < - ONE) THEN                                  !
         COS_TH = - ONE                                             !
      END IF                                                        !
      TH = ACOS(COS_TH)                                             !
!
      XPIY = X + IC * Y                                             !
      CALL ARCSIN(XPIY,COS_TH,PH)                                   !
!
  10  RETURN                                                        !
!
      END SUBROUTINE XYZ_TO_R_THETA_PHI
!
END MODULE VECTOR
