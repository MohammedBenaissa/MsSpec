!
!=======================================================================
!
MODULE CROSS_SECTIONS
!
!  This module provides the routines to compute the cross-section
!    of several spectroscopies
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE CROSS_SECTION(N_INFILES,J_FILE,IPHA,ZEM,NP)
!
!  This subroutine calls the calculation of the cross-section corresponding
!    to the spectroscopy selected
!
!
!  Input variables :
!
!                       N_INFILES :  number of input data files
!                       J_FILE    :  current input data file index
!                       IPHA      :  phagen_scf.f index
!                       ZEM       :  z position of absorber
!
!
!  Output variables :
!
!                       NP        :
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 28 May 2021
!
!
      USE EXTERNAL_O1
      USE EXTERNAL_O2
!
      USE PHI_O1
      USE PHI_O2
!
      USE THETA_O1
      USE THETA_O2
!
      USE ABSORBER
      USE CLUSTER
      USE CLUSTER_COORD
      USE CLUSTER_LIMITS
      USE CLUS_ELEC
      USE CURRENT_CALC
      USE CURRENT_F_TH
      USE CURRENT_SPIN
      USE EXP_TYPE
      USE TESTS
!
      USE CALC_PED_MS
      USE ELECTRON_CHOICE
      USE INITIALIZE_CALC
      USE SCATTERING_AMPLITUDE
      USE TREATMENT_PED
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  N_INFILES
      INTEGER, INTENT(IN)       ::  J_FILE,IPHA
      INTEGER, INTENT(OUT)      ::  NP
!
      INTEGER                   ::  JEL,J_EL
      INTEGER                   ::  I_EXT_A,NPHI_A,NTHETA_A
!
      REAL (WP), INTENT(IN)     ::  ZEM
!
!   Initialization of algorithm-dependent variables
!
      CALL INIT_ALGO(NATCLU)                                        !
!
!  Setting the Auger electron
!
      I_EXT_A  = 0                                                  !
      NPHI_A   = 1                                                  !
      NTHETA_A = 1                                                  !
!
      IF(SPECTRO == 'AED') THEN                                     !
        I_EXT_A  = I_EXT_O1                                         !
        NPHI_A   = NPHI_O1                                          !
        NTHETA_A = NTHETA_O1                                        !
      ELSEIF(SPECTRO == 'APC') THEN                                 !
        I_EXT_A  = I_EXT_O2                                         !
        NPHI_A   = NPHI_O2                                          !
        NTHETA_A = NTHETA_O2                                        !
      END IF                                                        !

!
!..........  Loop on the electrons involved in the spectroscopy       ..........
!
      DO JEL = 1, NEL                                               !
!
!..........  Check of the dimensioning of the treatment routine  ..........
!
        CALL STOP_TREAT(N_INFILES,NPLAN,NEMET,NE,NTHETA,NTHETA_A, & !
                        NPHI,NPHI_A,ISOM,I_EXT,I_EXT_A,SPECTRO)     !
!
!..........  Choosing the electron's beam and properties  ..........
!
        J_EL = NTAB_EL(JEL)                                         !
!
        CALL CHOOSE_ELEC(J_EL,NAT)                                  !
        CALL CHG_MOD(J_EL,NAT)                                      !
!
!  Setting the angular parameters for the current beam
!                     and checking them
!
        CALL INIT_CHECK_ANGLES(SPECTRO)                             !
!
!  Calculation of the scattering amplitude for electron J_EL
!
        IF(IFTHET > 0) THEN                                         !
          CALL CALC_SF(J_EL)                                        !
        END IF                                                      !
!
        IF(ISPIN == 0) THEN                                         !
!
!  Calling the MS description of the spectroscopy
!
          IF(SPECTRO == 'PED') THEN                                 !
            CALL PED_MS(NPLAN,ZEM,IPHA,N_INFILES,J_FILE,NP)         !
          ELSE IF(SPECTRO == 'LED') THEN                            !
!           CALL LED_MS(NPLAN,VAL,ZEM,IPHA,NAT,COORD,NATYP,        & !
!                          NATCLU,N_INFILES,J_FILE,NP)               !
          CONTINUE                                                  !
          ELSE IF(SPECTRO == 'AED') THEN                            !
!           CALL AED_MS(NPLAN,VAL,ZEM,IPHA,NATYP,N_INFILES,J_FILE,NP)!
          CONTINUE                                                  !
          ELSE IF(SPECTRO == 'XAS') THEN                            !
!           CALL XAS_MS(NPLAN,VAL,ZEM,IPHA,RHOR,N_INFILES,J_FILE,NP) !
          CONTINUE                                                  !
          ELSE IF(SPECTRO == 'APC') THEN                            !
!           IF(J_EL == 1) THEN                                       !
!             CALL PED_MS(NPLAN,VAL,ZEM,IPHA,NAT,COORD,NATYP,      & !
!                            NATCLU,N_INFILES,J_FILE,NP)             !
!           ELSE IF(J_EL == 2) THEN                                  !
!             CALL AED_MS(NPLAN,VAL,ZEM,IPHA,NAT,COORD,NATYP,      & !
!                           NATCLU,N_INFILES,J_FILE,NP,LE_MIN,LE_MAX)!
!           END IF
          CONTINUE                                                  !
          END IF                                                     !
!
        ELSE IF(ISPIN == 1) THEN                                    !
!
!         IF(SPECTRO == 'PED') THEN                                  !
!           CALL PED_MS(NPLAN,VAL,ZEM,IPHA,NAT,COORD,NATYP,        & !
!                          NATCLU,N_INFILES,J_FILE,NP)               !
!         ELSE IF(SPECTRO == 'AED') THEN                             !
!           CALL AED_MS                                              !
!         ELSE IF(SPECTRO == 'XAS') THEN                             !
!           CALL XAS_MS                                              !
!         END IF                                                     !
          CONTINUE                                                  !
!
        END IF                                                      !
!
!..........  End of the loop on the electrons  ..........
!
      END DO                                                        !
!
      END SUBROUTINE CROSS_SECTION
!
END MODULE CROSS_SECTIONS
