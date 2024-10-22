!
!=======================================================================
!
MODULE VIBRATIONS
!
!  This module contains functions to compute the effect of vibrations
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!=======================================================================
!
      FUNCTION UJ_SQ(JTYP)
!
!  This routine evaluates the mean square displacements UJ_SQ,
!    first along une direction (x, y or z): UJ2 within the Debye model,
!              using the Debye function formulation
!
!  X1 is the Debye function phi_1
!  UJ_SQ is given in unit of the square of the lattice parameter A0
!  Temperatures are expressed in Kelvin
!
!  The coefficient COEF equals:
!
!           3 hbar^{2} N_A 10^{3} / (4 k_B)
!
!    where N_A is the Avogadro number, k_B is Boltzmann's constant
!    and 10^3 arises from the fact that the atomic mass is
!    expressed in grams
!
!  Then UJ_SQ is obtained as UJ_SQ = (2 + RSJ) UJJ for surface atoms
!                            UJ_SQ = 3 UJJ for bulk atoms
!
!
!  For empty spheres, two possibilities are provided. By construction,
!    they are very light (their mass is taken as 1/1836 of the mass
!    of a H atom) and so they will vibrate a lot (IDCM = 1). When
!    setting IDCM = 2, their mean square displacement is set to a
!    tiny value so that they hardly vibrate (frozen empty spheres)
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified :  7 Jun 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE,THREE,FOUR,TTINY
!
      USE ATOMIC_MASS
      USE DEB_WAL_CLU
      USE CLUSTER
      USE VIBR_TYPE
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  JTYP
!
      INTEGER               ::  N,N_MAX
!
      REAL (WP)             ::  UJ_SQ
!
      REAL (WP)             ::  COEF,RZ2,LITTLE
      REAL (WP)             ::  AT,MJ,C,X1
      REAL (WP)             ::  Z,P1,UJJ
!
      REAL (WP)             ::  FLOAT,EXP
!
      DATA COEF    / 36.38156264E0_WP /                             ! 3 hbar^{2} N_A / (4 k_B)
      DATA RZ2     /  1.6449340668482264364724151666460E0_WP /      ! Pi^2 / 6
      DATA LITTLE  /  0.01E0_WP /                                   ! lowest temperature for
!                                                                   ! calculation  of phi_1
      N_MAX = 20                                                    !
!
!  Computation of the 1D mean square displacement UJ2
!
      AT = TD / T                                                   !
      MJ = XMT(JTYP)                                                !
      C  = COEF / (MJ * TD)                                         !
      X1 = ZERO                                                     !
!
      DO N = 1, N_MAX                                               !
        Z  = FLOAT(N)                                               !
        X1 = X1 + EXP(- Z * A) * (A + ONE / Z) / Z                  !
      END DO                                                        !
!
      P1  = ONE + FOUR * (RZ2 - X1) / (AT * AT)                     !
      UJJ = C * P1 / (A * A)                                        !
!
!  3D mean square displacement UJ_SQ
!
      IF(IMSD == 1) THEN                                            !
        UJ_SQ = ( THREE + FLOAT(I_FREE(JTYP)) * (RSJ - ONE) ) * UJJ !
      ELSE IF(IMSD == 2) THEN                                       !
        UJ_SQ = TTINY                                               !
      END IF                                                        !
!
      END FUNCTION UJ_SQ
!
!=======================================================================
!
      SUBROUTINE AV_T_MATRIX(JTYP,JE,X,TLT)
!
!  This routine recomputes the T-matrix elements taking into account the
!        mean square displacements (averaging over vibrations)
!
!  When the argument X is tiny, no vibrations are taken into account
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 16 Jun 2021
!
!
      USE REAL_NUMBERS,         ONLY : ZERO
      USE COMPLEX_NUMBERS,      ONLY : ZEROC
!
      USE CURRENT_T_MATRIX
!
      USE ANGULAR_MOMENTUM,     ONLY : GAUNT
      USE SPHERICAL_BESSEL
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  JTYP,JE
!
      INTEGER                   ::  LM2,IBESP
      INTEGER                   ::  MG1,MG2
      INTEGER                   ::  L,L1,L2
      INTEGER                   ::  L2MIN,L2MAX
!
      INTEGER                   ::  ABS
!
      REAL (WP), INTENT(IN)     ::  X
!
      REAL (WP)                 ::  PI4,EPS
      REAL (WP)                 ::  GNT(0:N_GAUNT)
      REAL (WP)                 ::  COEF,C
      REAL (WP)                 ::  XL,XL1,XL2
      REAL (WP)                 ::  SL2
!
      REAL (WP)                 ::  EXP,FLOAT,SQRT,REAL
!
      COMPLEX (WP), INTENT(OUT) ::  TLT(0:NT_M,4,NATM,NE_M)
!
      COMPLEX (WP)              ::  SL1
      COMPLEX (WP)              ::  FFL(0:2*NL_M)
!
      DATA PI4    / 12.5663706143591729538505735331180E0_WP /          ! 4 pi
      DATA EPS    / 1.0E-10_WP /
!
      IF(X > EPS) THEN                                              !
!
!  Standard case: vibrations
!
        COEF  = PI4 * EXP(- X)                                      !
        LM2   = 2 * LMAX(JTYP,JE)                                   !
        IBESP = 5                                                   !
        MG1   = 0                                                   !
        MG2   = 0                                                   !
!
        CALL BESPHE(LM2,IBESP,X,FFL)                                !
!
        DO L = 0, LMAX(JTYP,JE)                                     !
          XL  = FLOAT(L + L + 1)                                    !
          SL1 = ZEROC                                               !
!
          DO L1 = 0, LMAX(JTYP,JE)                                  !
            XL1 = FLOAT(L1 + L1 + 1)                                !
            CALL GAUNT(L,MG1,L1,MG2,GNT)                            !
            L2MIN = ABS(L1 - L)                                     !
            L2MAX = L1 + L                                          !
            SL2   = ZERO                                            !
!
            DO L2 = L2MIN, L2MAX, 2                                 !
              XL2 = FLOAT(L2 + L2 + 1)                              !
              C   = SQRT( XL1 * XL2 / (PI4 * XL) )                  !
              SL2 = SL2 + C * GNT(L2) * REAL(FFL(L2),KIND=WP)       !
            END DO                                                  !
!
            SL1 = SL1 + SL2 * TL(L1,1,JTYP,JE)                      !
          END DO                                                    !
!
          TLT(L,1,JTYP,JE) = COEF * SL1                             !
!
        END DO                                                      !
!
      ELSE                                                          !
!
!  Argument X tiny: no vibrations
!
        DO L = 0, LMAX(JTYP,JE)                                     !
!
          TLT(L,1,JTYP,JE) = TL(L,1,JTYP,JE)                        !
!
        END DO                                                      !
!
      END IF                                                        !

!
      END SUBROUTINE AV_T_MATRIX
!
END MODULE VIBRATIONS
