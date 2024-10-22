!
!=======================================================================
!
MODULE MEAN_FREE_PATH
!
!
!
      USE ACCURACY_REAL
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE MFP(E,XMFP,*)
!
!  This routine generates the electron mean free path
!
!  Input variables:
!
!     E        : electron kinetic energy (eV)
!
!
!  Output variable:
!
!     XMFP     : mean free path
!
!
!       IMFP = -1 : XMFP is set to 1.E+30
!       IMFP =  0 : XMFP is the value given in the input data file
!       IMFP =  1 : XMFP computed from Tokutaka et al, Surf. Sci. 149,349 (1985)
!       IMFP =  2 : XMFP computed from the Seah and Dench expression
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 28 May 2021
!
!
      USE REAL_NUMBERS,                ONLY : LARGE
!
      USE CLUSTER
      USE CURRENT_INIT_VAL
      USE CURRENT_MEAN_FREE_PATH
      USE OUTUNITS
      USE TESTS
!
      IMPLICIT NONE
!
      REAL (WP), INTENT(INOUT)  ::  E
      REAL (WP), INTENT(OUT)    ::  XMFP
!
      REAL (WP)                 ::  Q
      REAL (WP)                 ::  CSTE1,CSTE2,CSTE3
      REAL (WP)                 ::  A1,A2,A3
      REAL (WP)                 ::  XLN
!
      REAL (WP)                 ::  FLOAT,LOG,SQRT
!
      E = E + VINT                                                  !
!
      IF(IMFP == -1) THEN                                           !
        XMFP = LARGE                                                !
      ELSE IF(IMFP == 0) THEN                                       !
        XMFP = XMFP0                                                !
      ELSE IF(IMFP == 1) THEN                                       !
        Q     = FLOAT(NZ_A) * RHOT_A / XMT_A                        !
        CSTE1 = LOG(Q / 4.50E0_WP) / (LOG(7.74E0_WP / 4.50E0_WP))   !
        CSTE2 = LOG(Q / 3.32E0_WP) / (LOG(7.74E0_WP / 3.32E0_WP))   !
        CSTE3 = LOG(Q / 3.32E0_WP) / (LOG(4.50E0_WP / 3.32E0_WP))   !
        A1    =   0.7271E0_WP + 0.2595E0_WP * LOG(E)                !
        A2    = - 3.2563E0_WP + 0.9395E0_WP * LOG(E)                !
        A3    = - 2.5716E0_WP + 0.8226E0_WP * LOG(E)                !
!
        IF(E >= 350.E0_WP) GO TO 10                                 !
!
        XLN = CSTE1 * (0.0107E0_WP - 0.0083E0_WP * LOG(E)) + A1     !
        GO TO 20                                                    !
!
  10    IF((NZ_A >= 24) .AND. (NZ_A <= 74)) GO TO 30                !
        XLN = CSTE2 * (1.6551E0_WP - 0.2890E0_WP * LOG(E)) + A2     !
        GO TO 20                                                    !
!
  30    IF(NZ_A >= 42) GO TO 40                                     !
        XLN = CSTE3 * (0.6847E0_WP - 0.1169E0_WP * LOG(E)) + A2     !
        GO TO 20                                                    !
!
  40    XLN = CSTE1 * (0.9704E0_WP - 0.1721E0_WP * LOG(E)) + A3     !
!
  20    XMFP = EXP(XLN)                                             !
      ELSE IF(IMFP == 2) THEN                                       !
        XMFP = 1430.E0_WP / (E**2) + 0.54E0_WP * SQRT(E)            !
      ELSE                                                          !
        RETURN 1                                                    !
      END IF                                                        !
!
      E = E - VINT                                                  !
      IF(IPRINT > 0) WRITE(IUO1,80) E,XMFP                          !
!
!  Fomat:
!
  80  FORMAT(/////,2X,'=========  E = ',F7.2,' eV',5X,'MEAN',  &
             ' FREE PATH = ',F6.3,' ANGSTROEMS  ','=========')
!
      END SUBROUTINE MFP
!
END MODULE MEAN_FREE_PATH
