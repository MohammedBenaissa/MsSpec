!
!=======================================================================
!
MODULE ARC_SIN
!
!  This module contains the arcsin function
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE ARCSIN(U,CST,RANGLE)
!
!   For a given complex number U, this subroutine calculates its phase
!   Warning : it is implicitely assumed that U = r sin(theta) exp(i*phi)
!   with theta > or = to 0 which is always the case when theta is obtained
!   from the coordinates of a given vector r by the ACOS intrinsic function.
!
!
!  Input variables :
!
!                       U         :  complex number
!                       CST       :  cos(theta)
!
!
!  Output variables :
!
!                       RANGLE    :  phase phi (in radian)
!
!
!
!   When sin(theta) = 0, then phi = 0  if cos(theta) =  1
!                             phi = pi if cos(theta) = -1
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified :  7 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE,TWO,FOUR,HALF,SMALL
      USE COMPLEX_NUMBERS,   ONLY : IC
      USE PI_ETC,            ONLY : PI
!
      IMPLICIT NONE
!
      REAL (WP), INTENT(IN)    ::  CST
      REAL (WP), INTENT(OUT)   ::  RANGLE
!
      REAL (WP)                ::  REAL
!
      COMPLEX(WP), INTENT(IN ) ::  U
!
      COMPLEX(WP)              ::  CANGLE
      COMPLEX(WP)              ::  ABS,LOG
!
      IF(ABS(U) < SMALL) THEN                                       !
        IF(CST > ZERO) THEN                                         !
          RANGLE = ZERO                                             !
        ELSE IF(CST <  ZERO) THEN                                   !
          RANGLE = PI                                               !
        END IF                                                      !
      ELSE                                                          !
        CANGLE = - IC * LOG(U / ABS(U))                             !
        RANGLE = REAL(CANGLE,KIND=WP)                               !
      END IF                                                        !
!
      END SUBROUTINE ARCSIN
!
END MODULE ARC_SIN
