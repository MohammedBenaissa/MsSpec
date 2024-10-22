!
!=======================================================================
!
MODULE EULER_ANGLES
!
!  This module computes the Enler angles
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE EULER(RTHETA1,RPHI1,RTHETA2,RPHI2,                 &
                       RALPHA,RBETA,RGAMMA,IROT)
!
!  This routine calculates the Euler angles RALPHA,RBETA,RGAMMA corresponding
!       to the rotation r1(RTHETA1,RPHI1) ----> r2(RTHETA2,RPHI2)
!
!       IROT = 1 : r ---> z represented by (0,RTHETA,PI-RPHI)
!       IROT = 0 : r ---> z represented by (0,-RTHETA,-RPHI)
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 19 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE
      USE COMPLEX_NUMBERS,   ONLY : IC
      USE PI_ETC,            ONLY : PI
!
      USE ARC_SIN
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  IROT
!
      INTEGER                  ::  EPS
!
      REAL (WP), INTENT(IN)    ::  RTHETA1,RPHI1
      REAL (WP), INTENT(IN)    ::  RTHETA2,RPHI2
      REAL (WP), INTENT(OUT)   ::  RALPHA,RBETA,RGAMMA
!
      REAL (WP)                ::  SMALL
      REAL (WP)                ::  DPHI
      REAL (WP)                ::  A1,A2,A3,A4
      REAL (WP)                ::  U3
!
      REAL (WP)                ::  SIN,COS,ACOS,ABS
!
      COMPLEX(WP)              ::  U1,U2
!
      DATA SMALL  / 0.0001E0_WP /
!
      IF(IROT == 1) THEN                                            !
        EPS =   1                                                   !
      ELSE                                                          !
        EPS = - 1                                                   !
      END IF                                                        !
!
      DPHI = RPHI2 - RPHI1                                          !
!
      A1 = SIN(RTHETA1) * COS(RTHETA2)                              !
      A2 = COS(RTHETA1) * SIN(RTHETA2)                              !
      A3 = COS(RTHETA1) * COS(RTHETA2)                              !
      A4 = SIN(RTHETA1) * SIN(RTHETA2)                              !
!
      U1 = A1 - A2 * COS(DPHI) - IC * SIN(RTHETA2) * SIN(DPHI)      !
      U2 = A1 * COS(DPHI) - A2 + IC * SIN(RTHETA1) * SIN(DPHI)      !
      U3 = A3 + A4 * COS(DPHI)                                      !
!
      IF(U3 > ONE)   U3 =   ONE                                     !
      IF(U3 < - ONE) U3 = - ONE                                     !
      RBETA = ACOS(U3)                                              !
!
      IF(ABS(SIN(RBETA)) > SMALL) THEN                              !
!
        U1 = EPS * U1 / SIN(RBETA)                                  !
        U2 = EPS * U2 / SIN(RBETA)                                  !
!
        CALL ARCSIN(U1,U3,RALPHA)                                   !
        CALL ARCSIN(U2,U3,RGAMMA)                                   !
!
      ELSE
!
        RALPHA = ZERO                                               !
        IF(ABS(U3 - ONE) < SMALL) THEN                              !
          RGAMMA = ZERO                                             !
        ELSE                                                        !
          RGAMMA = PI                                               !
        END IF                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE EULER
!
END MODULE EULER_ANGLES
