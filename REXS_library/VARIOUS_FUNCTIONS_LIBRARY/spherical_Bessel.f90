!
!=======================================================================
!
MODULE SPHERICAL_BESSEL
!
!  This module contains routines to compute the spherical Bessel functions
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,        ONLY : N_BESS,NL_M
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE BESPHE(L_MAX,IBES,X,FL)
!
!  This routine computes the spherical Bessel functions for
!                       a real argument X.
!
!  Input variables:
!
!     L_MAX     : upper dvalue of BESPHE index
!     IBES      : type of spherical Bessel function calculated
!                   IBES = 1 : Bessel function
!                   IBES = 2 : Neumann function
!                   IBES = 3 : Hankel function of the first kind
!                   IBES = 4 : Hankel function of the second kind
!                   IBES = 5 : Modified Bessel function
!                   IBES = 6 : Modified Neumann function
!                   IBES = 7 : Modified Hankel function
!     X         : argument of spherical Bessel function
!
!  Output variable:
!
!     FL        : spherical Bessel function array
!
!
!
!   Author :  D. Sébilleau
!
!                                         Last modified : 16 Jun 2021
!
      USE OUTUNITS
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE,SMALL,TTINY,LARGE
      USE COMPLEX_NUMBERS
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  L_MAX,IBES
!
      INTEGER                  ::  ITEST,IBES1
      INTEGER                  ::  NX,NREC,L
!
      INTEGER                  ::  INT
!
      REAL (WP), INTENT(IN)    ::  X
!
      REAL (WP)                ::  SCALN(0:N_BESS)
      REAL (WP)                ::  ECH,COMP,COMM,DEB,A,B
      REAL (WP)                ::  ECHEL,REN
!
      REAL (WP)                ::  MAX,ABS,SIN,COS,SINH,COSH,EXP,FLOAT
!
      COMPLEX(WP), INTENT(OUT) ::  FL(0:2*NL_M)
!
      COMPLEX(WP)              ::  FLNN(0:N_BESS),GL(0:N_BESS)
      COMPLEX(WP)              ::  C1,C2,CNORM
!
      ECH  = 30.0E0_WP                                              !
      COMP = LARGE                                                  !
      COMM = TTINY                                                  !
!
      NX   = INT(X)                                                 !
      NREC = 5 * MAX(L_MAX,NX)                                      !
!
      IF(NREC > N_BESS) GO TO 16                                    !
!
      ITEST = 0                                                     !
      C1    = ONEC                                                  !
      C2    = IC                                                    !
      DEB   = ONE                                                   !
      IF((IBES == 3) .OR. (IBES == 4)) THEN                         !
        IBES1 = 1                                                   !
        IF(IBES == 4) C2 = - IC                                     !
      ELSE IF(IBES == 7) THEN                                       !
        IBES1 = 5                                                   !
        C2    = - ONEC                                              !
      ELSE                                                          !
        IBES1 = IBES                                                !
      END IF                                                        !
!
!  Initialization of the spherical Bessel function
!
      DO L = 0, 2 * NL_M                                            !
        FL(L) = ZEROC                                               !
      END DO                                                        !
!
!   Case where the argument is zero
!
      IF(ABS(X) < SMALL) THEN                                       !
        IF((IBES == 1) .OR. (IBES == 5)) THEN                       !
          FL(0) = ONEC                                              !
          DO L = 1, L_MAX                                           !
            FL(L) = ZEROC                                           !
          END DO                                                    !
          ITEST =   1                                               !
        ELSE                                                        !
          ITEST = - 1                                               !
        END IF                                                      !
      END IF                                                        !
!
      IF(ITEST) 11,12,13                                            !
!
  11  WRITE(IUO1,14)                                                !
      STOP                                                          !
  16  WRITE(IUO1,17) NREC                                           !
      STOP                                                          !
  15  IBES1 = IBES1 + 1                                             !
!
!   Initial values
!
  12  A = - ONE                                                     !
      B =   ONE                                                     !
!
      IF(IBES1 == 1) THEN                                           !
        FL(0)         = ONEC * SIN(X) / X                           !
        FLNN(NREC)    = ZEROC                                       !
        SCALN(NREC)   = ZERO                                        !
        FLNN(NREC-1)  = ONEC * DEB                                  !
        SCALN(NREC-1) = ZERO                                        !
      ELSE IF(IBES1 == 2) THEN                                      !
        GL(0)         = - ONEC * COS(X) / X                         !
        GL(1)         =   GL(0) / X - SIN(X) / X                    !
      ELSE IF(IBES1 == 5) THEN                                      !
        A             =   ONE                                       !
        B             = - ONE                                       !
        FL(0)         = ONEC * SINH(X) / X                          !
        FLNN(NREC)    = ZEROC                                       !
        SCALN(NREC)   = ZERO                                        !
        FLNN(NREC-1)  = ONEC * DEB                                  !
        SCALN(NREC-1) = ZERO                                        !
      ELSE IF(IBES1 == 6) THEN                                      !
        A             =   ONE                                       !
        B             = - ONE                                       !
        GL(0)         = ONEC * COSH(X) / X                          !
        GL(1)         = (SINH(X) - GL(0)) / X                       !
      END IF                                                        !
!
!   Downward reccurence for the spherical Bessel function
!
      IF((IBES1 == 1) .OR. (IBES1 == 5)) THEN                       !
!
        DO L = NREC - 1, 1, -1                                      !
          ECHEL      = ZERO                                         !
          SCALN(L-1) = SCALN(L)                                     !
          REN        = EXP( SCALN(L) - SCALN(L+1) )                 !
          FLNN(L-1)  = A * ( REN * FLNN(L+1) -                    & !
                             B * FLOAT(L + L + 1) * FLNN(L) / X   & !
                           )                                        !
          IF(CDABS(FLNN(L-1)) > COMP) THEN                          !
            ECHEL = - ECH                                           !
          ELSE IF(CDABS(FLNN(L-1)) < COMM) THEN                     !
            ECHEL =   ECH                                           !
          END IF                                                    !
          IF(ECHEL /= ZERO) SCALN(L-1) = ECHEL + SCALN(L-1)         !
          FLNN(L-1) = FLNN(L-1) * EXP(ECHEL)                        !
        END DO                                                      !
!
        CNORM = FL(0) / FLNN(0)                                     !
!
        DO L = 1, L_MAX                                             !
          FL(L) = CNORM * FLNN(L) * EXP( SCALN(0) - SCALN(L) )      !
        END DO                                                      !
!
      ELSE                                                          !
!
!   Upward recurrence for the spherical Neumann function
!
        DO L = 1, L_MAX                                             !
          IF(IBES == 7) C1 = (- ONEC)**(L+2)                        !
          GL(L+1) = A * GL(L-1) + B * FLOAT(L + L + 1) * GL(L) / X  !
          IF(IBES1 /= IBES) THEN                                    !
!
!   Calculation of the spherical Hankel function
!
            FL(L+1) = C1 * ( FL(L+1) + C2 * GL(L+1) )               !
          ELSE                                                      !
            FL(L+1) = GL(L+1)                                       !
          END IF                                                    !
        END DO                                                      !
        IF(IBES1 == IBES) THEN                                      !
          FL(0) = GL(0)                                             !
          FL(1) = GL(1)                                             !
        ELSE                                                        !
          FL(0) = C1 * ( FL(0) + C2 * GL(0) )                       !
          FL(1) = C1 * ( FL(1) + C2 * GL(1) )                       !
        END IF                                                      !
        IBES1 = IBES                                                !
!
      END IF                                                        !
!
      IF(IBES /= IBES1) GO TO 15                                    !
!
  13  RETURN                                                        !
!
! Formats:
!
  14  FORMAT(/////,3X,'<<<<<<<<<< THE ARGUMENT OF THE BESSEL ',  &
            'FONECCTIONS IS NUL >>>>>>>>>>')
  17  FORMAT(/////,3X,'<<<<<<<<<< THE DIMENSIONNING N_BESS ',    &
            'IS NOT CORRECT FOR SUBROUTINE BESPHE >>>>>>>>>>',//,&
            15X,'<<<<<<<<<< IT SHOULD BE AT LEAST : ',I5,        &
            ' >>>>>>>>>>')
!
      END SUBROUTINE BESPHE
!
END MODULE SPHERICAL_BESSEL
