!
!=======================================================================
!
MODULE INITIALIZE_CALC
!
!  This module provides set of routine to initialize various
!    variables, functions, ... necessary to perform the calculations.
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE INIT_CORE_STATE(NLI,I_SO,S_O,LI,IRET,*)
!
!  This subroutine initializes cote state parameters
!
!
!  Input variables :
!
!                       NLI       :  core state orbital (1-character string: s,p,d,f,g)
!                       I_SO      :  switch for spin-orbit-resolved calculation
!                       S_O       :  spin-orbit component (3-character string: 1/2,3/2,...)
!
!
!  Output variables :
!
!                       LI        :  core state angular momentum
!                       IRET      :  error return code
!                       *         :  label return for error in input data file
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified :  5 May 2021
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  1)  ::  NLI
      CHARACTER (LEN =  3)  ::  S_O
!
      INTEGER, INTENT(IN)   ::  I_SO
      INTEGER, INTENT(OUT)  ::  LI,IRET
!
      IF((NLI == 's') .OR. (NLI == 'S')) THEN                       !
        LI = 0                                                      !
      ELSE IF((NLI == 'p') .OR. (NLI == 'P')) THEN                  !
        LI = 1                                                      !
      ELSE IF((NLI == 'd') .OR. (NLI == 'D')) THEN                  !
        LI = 2                                                      !
      ELSE IF((NLI == 'f') .OR. (NLI == 'F')) THEN                  !
        LI = 3                                                      !
      ELSE IF((NLI == 'g') .OR. (NLI == 'G')) THEN                  !
        LI = 4                                                      !
      ELSE                                                          !
        IRET = 5                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
      IF(LI > LI_M) THEN                                            !
        IRET = 6                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
      IF(I_SO <= 0) THEN                                            !
        S_O = '   '                                                 !
      ELSE IF(I_SO == 1) THEN                                       !
        IF(S_O == '1/2') THEN                                       !
          IF(LI > 1) IRET = 7                                       !
        ELSE IF(S_O == '3/2') THEN                                  !
          IF((LI < 1) .OR. (LI > 2)) IRET = 7                       !
        ELSE IF(S_O == '5/2') THEN                                  !
          IF((LI < 2) .OR. (LI > 3)) IRET = 7                       !
        ELSE IF(S_O == '7/2') THEN                                  !
          IF((LI < 3) .OR. (LI > 4)) IRET = 7                       !
        ELSE IF(S_O == '9/2') THEN                                  !
          IF(LI /= 4) IRET = 7                                      !
        END IF                                                      !
      ELSE IF(I_SO == 2) THEN                                       !
        S_O = '   '                                                 !
      END IF                                                        !
!
      END SUBROUTINE INIT_CORE_STATE
!
!=======================================================================
!
      SUBROUTINE INIT_MSD
!
!  This subroutine initializes the mean square displacements
!
!
!
!  internal variable (contained in module DEB_WAL_CLU) :
!
!                       I_MSD     :  index characterizing the mean square displacements
!                                       ---> 0 : read-in from input data file
!                                       ---> 1 : computed internally
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 31 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO
!
      USE CLUSTER
      USE DEB_WAL_CLU
      USE VIBRATIONS
!
      IMPLICIT NONE
!
      INTEGER               ::  JTYP
!
      REAL (WP)             ::  UJ2_OLD(NATM)
!
!  Storing read-in values in UJ2_OLD
!
      DO JTYP = 1, NAT                                              !
        UJ2_OLD(JTYP) = UJ2(JTYP)                                   !
      END DO                                                        !
!
      IF (IMSD == 1) THEN                                           !
        UJ2(JTYP) = UJ_SQ(JTYP)                                     !
      END IF                                                        !
!
      END SUBROUTINE INIT_MSD
!
!=======================================================================
!
      SUBROUTINE INIT_SPECTRO(SPECTRO,SELRULE,STEREO,I_SPEC,        &
                              EXCITATION,INTERACT_I,INTERACT_O)
!
!  This subroutine initializes various spectroscopy-dependent parameters
!
!
!  Input variables :
!
!                       SPECTRO   :  type of spectroscopy:
!                                       ---> 'PED' : photoelectron diffraction
!                                            'LED' : low-energy electron diffraction
!                                            'XAS' : X-ray absorption spectroscopy
!                                            'AED' : Auger electron diffraction
!                                            'APC' : Auger photoelectron coincidence spectroscopy
!                                            'E2E' : (e,2e) spectroscopy
!                                            'E3E' : (e,3e) spectroscopy
!                                            'ELS' : electron energy loss spectroscopy
!                                            'RES' : resonant elastic X-ray scattering
!                                            'STM' : scanning tunneling microscopy
!                                            'PLS' : photoelectron energy loss spectroscopy
!                                            'BEM' : ballistic energy electron microscopy
!                                            'EIG' : eigenvalue calculation
!                       SELRULE   :  type of selection rule (a 5-character string)
!                       STEREO    :  switch for stereographic projection: 'YES'/' NO'
!                       I_SPEC    :  index characterizing the spectroscopy-dependent blocks
!                                    to read in the input data file
!
!  Output variables :
!
!                       EXCITATION:  way the electron is created:
!                                       ---> 'NOINTER' : electron already exists
!                                       ---> 'PH-ELEC' : electron excited by photon
!                                       ---> 'COUL-LL' : Coulomb with two electrons localized
!                                       ---> 'COUL-DL' : Coulomb with one electron delocalized
!                                       ---> 'COUL-DD' : Coulomb with two electrons delocalized
!                       INTERACT_I:  excitation of the incoming beam (REXS case only)
!                       INTERACT_O:  excitation of the outgoing beam (REXS case only)
!
!
!  Here, we have:
!                       NCL       :  number of clusters
!                       NEL       :  number of electrons involved in the spectroscopy
!                       NTAB_EL   :  array containing the type of electrons involved
!                                       ---> 1 : incoming electron
!                                       ---> 2 : excited electron (non detected)
!                                       ---> 3 : first outgoing electron (detected)
!                                       ---> 4 : second outgoing electron (detected)
!                                       ---> 5 : third outgoing electron (detected)
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 28 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO
      USE CLUS_ELEC
      USE CURRENT_CALC
      USE CURRENT_FINA_VAL
      USE CURRENT_INIT_VAL
      USE DICHROISM
      USE EXP_TYPE,            ONLY : IMOD
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  1)  ::  CHAN
      CHARACTER (LEN =  2)  ::  MPOLE
      CHARACTER (LEN =  3)  ::  SPECTRO,STEREO
      CHARACTER (LEN =  5)  ::  SELRULE
      CHARACTER (LEN =  7)  ::  EXCITATION(5),INTERACT_I,INTERACT_O
!
      INTEGER, INTENT(IN)   ::  I_SPEC
!
      INTEGER               ::  JEL
!
!   Initialization of the switches giving the kind
!     of electron involved (incoming, excited,
!     1st outgoing, 2nd outgoing, 3rd outgoing)
!     and the possible existence of a tip over the cluster
!
      INC  = 0                                                      !
      EXC  = 0                                                      !
      OUT1 = 0                                                      !
      OUT2 = 0                                                      !
      OUT3 = 0                                                      !
      TIP  = 0                                                      !
      CLU  = 1                                                      !
      EIG  = 0                                                      !
      PLA  = 0                                                      !
!
!   Spectroscopy-dependent variables
!
      DO JEL = 1, 5                                                 !
        EXCITATION(JEL) = 'NOINTER'                                 !
        NTAB_EL(JEL)   = 0                                          !
      END DO                                                        !
!
      IF(SPECTRO == 'PED') THEN                                     !
        EXCITATION(3) = 'PH-ELEC'                                   !
        NCL           = 1                                           !
        NEL           = 1                                           !
        NTAB_EL(1)    = 3                                           !
        OUT1          = 1                                           !
      ELSE IF(SPECTRO == 'LED') THEN                                !
        NCL           = 1                                           !
        NEL           = 2                                           !
        NTAB_EL(1)    = 1                                           !
        NTAB_EL(2)    = 3                                           !
        INC           = 1                                           !
        OUT1          = 1                                           !
      ELSE IF(SPECTRO == 'XAS') THEN                                !
        EXCITATION(3) = 'PH-ELEC'                                   !
        IF(I_DICHR > 1) THEN                                        !
          PRINT 100                                                 !
          STOP                                                      !
        END IF                                                      !
        NCL           = 1                                           !
        NEL           = 1                                           !
        NTAB_EL(1)    = 2                                           !
        EXC           = 1                                           !
      ELSE IF(SPECTRO == 'AED') THEN                                !
        EXCITATION(3) = 'COUL-LL'                                   !
        NCL           = 1                                           !
        NEL           = 1                                           !
        NTAB_EL(1)    = 3                                           !
        OUT1          = 1                                           !
      ELSE IF(SPECTRO == 'APC') THEN                                !
        EXCITATION(3) = 'PH-ELEC'                                   !
        EXCITATION(4) = 'COUL-LL'                                   !
        NCL           = 1                                           !
        NEL           = 2                                           !
        NTAB_EL(1)    = 3                                           !
        NTAB_EL(2)    = 4                                           !
        OUT1          = 1                                           !
        OUT2          = 1                                           !
      ELSE IF(SPECTRO == 'E2E') THEN                                !
        EXCITATION(4) = 'COUL-DL'                                   !
        NCL           = 1                                           !
        NEL           = 3                                           !
        NTAB_EL(1)    = 1                                           !
        NTAB_EL(2)    = 3                                           !
        NTAB_EL(3)    = 4                                           !
        INC           = 1                                           !
        OUT1          = 1                                           !
        OUT2          = 1                                           !
      ELSE IF(SPECTRO == 'E3E') THEN                                !
        EXCITATION(4) = 'COUL-DL'                                   !
        EXCITATION(5) = 'COUL-DL'                                   !
        NCL           = 1                                           !
        NEL           = 4                                           !
        NTAB_EL(1)    = 1                                           !
        NTAB_EL(2)    = 3                                           !
        NTAB_EL(3)    = 4                                           !
        NTAB_EL(4)    = 5                                           !
        INC           = 1                                           !
        OUT1          = 1                                           !
        OUT2          = 1                                           !
        OUT3          = 1                                           !
      ELSE IF(SPECTRO == 'ELS') THEN                                !
        EXCITATION(2) = 'COUL-DL'                                   !
        NCL           = 1                                           !
        NEL           = 3                                           !
        NTAB_EL(1)    = 1                                           !
        NTAB_EL(2)    = 2                                           !
        NTAB_EL(3)    = 3                                           !
        INC           = 1                                           !
        EXC           = 1                                           !
        OUT1          = 1                                           !
      ELSE IF(SPECTRO == 'RES') THEN                                !
        INTERACT_I = 'PH-ELEC'                                      !
        INTERACT_O = 'PH-ELEC'                                      !
        NCL           = 1                                           !
        NEL           = 1                                           !
        NTAB_EL(1)    = 2                                           !
        EXC           = 1                                           !
      ELSE IF(SPECTRO == 'STM') THEN                                !
        NCL           = 2                                           !
        NEL           = 2                                           !
        NTAB_EL(1)    = 1                                           !
        NTAB_EL(2)    = 3                                           !
        INC           = 1                                           !
        OUT1          = 1                                           !
        TIP           = 1                                           !
      ELSE IF(SPECTRO == 'PLS') THEN                                !
        EXCITATION(3) = 'PH-ELEC'                                   !
        IF(I_DICHR > 1) THEN                                        !
          PRINT 101                                                 !
          STOP                                                      !
        END IF                                                      !
        NCL           = 1                                           !
        NEL           = 1                                           !
        NTAB_EL(1)    = 3                                           !
        OUT1          = 1                                           !
        PLA           = 1                                           !
      ELSE IF(SPECTRO == 'EIG') THEN                                !
        NCL           = 1                                           !
        NEL           = 1                                           !
        NTAB_EL(1)    = 2                                           !
        EXC           = 1                                           !
        EIG           = 1                                           !
      END IF                                                        !

!   Check of the validity of the selection rule in the single channel case
!
      IF(SELRULE(3:3) == '>') THEN                                  !
        CHAN  = SELRULE(2:2)                                        !
        MPOLE = SELRULE(4:5)                                        !
!
        IF(MPOLE == 'E1') THEN                                      !
          IF(CHAN /= '1') THEN                                      !
            SELRULE = '+1>E1'                                       !
            PRINT 102                                               !
          END IF                                                    !
        ELSE IF(MPOLE == 'E2') THEN                                 !
          IF((CHAN /= '0').OR.(CHAN /= '2')) THEN                   !
            SELRULE = '+2>E2'                                       !
            PRINT 103                                               !
          END IF                                                    !
        ELSE IF(MPOLE == 'E3') THEN                                 !
          IF((CHAN /= '1').OR.(CHAN /= '3')) THEN                   !
            SELRULE = '+3>E3'                                       !
            PRINT 104                                               !
          END IF                                                    !
        ELSE IF(MPOLE == 'M1') THEN                                 !
          IF(CHAN /= '0') THEN                                      !
            SELRULE = '+0>M1'                                       !
            PRINT 105                                               !
          END IF                                                    !
        ELSE IF(MPOLE == 'M2') THEN                                 !
        ELSE                                                        !
          SELRULE = '+1>M2'                                         !
          PRINT 106                                                 !
        END IF                                                      !
      END IF                                                        !
!
!  Non Auger electron case
!
      IF(I_SPEC <= 2) THEN                                          !
        IF(IPHI == -1) THEN                                         !
          IPHI   = 1                                                !
          I_EXT  = 0                                                !
          STEREO = 'YES'                                            !
          IF(ABS(PHI1-PHI0) < 0.0001E0_WP) THEN                     !
            PHI0 = ZERO                                             !
            PHI1 = 360.0E0_WP                                       !
            NPHI = 361                                              !
          END IF                                                    !
          IF(ABS(THETA1-THETA0) < 0.0001E0_WP) THEN                 !
            THETA0 = ZERO                                           !
            THETA1 = 88.0E0_WP                                      !
            NTHETA = 89                                             !
          END IF
        ELSE IF(IPHI == 2) THEN                                     !
          IPHI   = 1                                                !
          I_EXT  = 1                                                !
        ELSE IF(IPHI == 3) THEN                                     !
          IPHI   = 1                                                !
          I_EXT  = - 1                                              !
        ELSE IF(ITHETA == 2) THEN                                   !
          ITHETA = 1                                                !
          I_EXT  = 1                                                !
        ELSE IF(ITHETA == 3) THEN                                   !
          ITHETA = 1                                                !
          I_EXT  = - 1                                              !
        ELSE IF(IE == 2) THEN                                       !
          IE     = 1                                                !
          I_EXT  = 1                                                !
        ELSE IF(IE == 3) THEN                                       !
          IE     = 1                                                !
          I_EXT  = - 1                                              !
        ELSE IF(IE == 4) THEN                                       !
          IF(SPECTRO == 'PED') THEN                                 !
            IE    = 1                                               !
            I_EXT = 2                                               !
            IMOD  = 0                                               !
          ELSE                                                      !
            IE    = 1                                               !
            I_EXT = 1                                               !
          END IF                                                    !
        END IF                                                      !
      END IF                                                        !
!
!  Formats
!
 100  FORMAT(///,4X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE A SPIN RESOLVED ', &
             'EXAFS EXPERIMENT : DECREASE I_DICHR  >>>>>>>>>>')
 101  FORMAT(///,4X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE A SPIN RESOLVED ', &
             'REXS EXPERIMENT : DECREASE I_DICHR  >>>>>>>>>>')
 102  FORMAT(///,13X,'---> ERROR IN SELRULE: CONTINUING WITH SELRULE',  &
             ' = +1>E1  ')
 103  FORMAT(///,13X,'---> ERROR IN SELRULE: CONTINUING WITH SELRULE',  &
             ' = +2>E2  ')
 104  FORMAT(///,13X,'---> ERROR IN SELRULE: CONTINUING WITH SELRULE',  &
             ' = +3>E3  ')
 105  FORMAT(///,13X,'---> ERROR IN SELRULE: CONTINUING WITH SELRULE',  &
             ' = +0>M1  ')
 106  FORMAT(///,13X,'---> ERROR IN SELRULE: CONTINUING WITH SELRULE',  &
             ' = +1>M2  ')
!
      END SUBROUTINE INIT_SPECTRO
!
!=======================================================================
!
      SUBROUTINE BEAM_INI(SPECTRO,MODE)
!
!  This subroutine initializes the beam energy and angular
!     parameters and stores them in the corresponding
!     commons
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 26 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO
!
      USE AMPLI_IN
      USE ENERGY_IN
      USE INTEGR_IN
      USE PHI_IN
      USE POLAR_IN
      USE THETA_IN
!
      USE AMPLI_EX
      USE ENERGY_EX
!
      USE AMPLI_O1
      USE ENERGY_O1
      USE INTEGR_O1
      USE PHI_O1
      USE POLAR_O1
      USE THETA_O1
!
      USE AMPLI_O2
      USE ENERGY_O2
      USE INTEGR_O2
      USE PHI_O2
      USE POLAR_O2
      USE THETA_O2
!
      USE AMPLI_O3
      USE ENERGY_O3
      USE INTEGR_O3
      USE PHI_O3
      USE POLAR_O3
      USE THETA_O3
!
      USE CURRENT_CALC
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  SPECTRO,MODE
!
      REAL (WP)             ::  E_RANGE
!
!  Initialization to zero of the beam characteristics
!
      I_AMP_IN = 0                                                  !
      I_AMP_EX = 0                                                  !
      I_AMP_O1 = 0                                                  !
      I_AMP_O2 = 0                                                  !
      I_AMP_O3 = 0                                                  !
!
      I_INT_IN = 0                                                  !
      I_INT_O1 = 0                                                  !
      I_INT_O2 = 0                                                  !
      I_INT_O3 = 0                                                  !
!
      IPOL_IN = 0                                                   !
      IPOL_O1 = 0                                                   !
      IPOL_O2 = 0                                                   !
      IPOL_O3 = 0                                                   !
!
      IE_IN = 0                                                     !
      IE_EX = 0                                                     !
      IE_O1 = 0                                                     !
      IE_O2 = 0                                                     !
      IE_O3 = 0                                                     !
!
      ITHETA_IN = 0                                                 !
      ITHETA_O1 = 0                                                 !
      ITHETA_O2 = 0                                                 !
      ITHETA_O3 = 0                                                 !
!
      IPHI_IN = 0                                                   !
      IPHI_O1 = 0                                                   !
      IPHI_O2 = 0                                                   !
      IPHI_O3 = 0                                                   !
!
      NE_IN = 0                                                     !
      NE_EX = 0                                                     !
      NE_O1 = 0                                                     !
      NE_O2 = 0                                                     !
      NE_O3 = 0                                                     !
!
      NPHI_IN = 0                                                   !
      NPHI_O1 = 0                                                   !
      NPHI_O2 = 0                                                   !
      NPHI_O3 = 0                                                   !
!
      NTHETA_IN = 0                                                 !
      NTHETA_O1 = 0                                                 !
      NTHETA_O2 = 0                                                 !
      NTHETA_O3 = 0                                                 !
!
      E0_IN = ZERO                                                  !
      E0_EX = ZERO                                                  !
      E0_O1 = ZERO                                                  !
      E0_O2 = ZERO                                                  !
      E0_O3 = ZERO                                                  !
!
      E1_IN = ZERO                                                  !
      E1_EX = ZERO                                                  !
      E1_O1 = ZERO                                                  !
      E1_O2 = ZERO                                                  !
      E1_O3 = ZERO                                                  !
!
      PHI0_IN = ZERO                                                !
      PHI0_O1 = ZERO                                                !
      PHI0_O2 = ZERO                                                !
      PHI0_O3 = ZERO                                                !
!
      PHI1_IN = ZERO                                                !
      PHI1_O1 = ZERO                                                !
      PHI1_O2 = ZERO                                                !
      PHI1_O3 = ZERO                                                !
!
      THETA0_IN = ZERO                                              !
      THETA0_O1 = ZERO                                              !
      THETA0_O2 = ZERO                                              !
      THETA0_O3 = ZERO                                              !
!
      THETA1_IN = ZERO                                              !
      THETA1_O1 = ZERO                                              !
      THETA1_O2 = ZERO                                              !
      THETA1_O3 = ZERO                                              !
!
!  Case of specific spectroscopies
!
      IF(SPECTRO == 'ELS') THEN                                     !
!
        IF(MODE == 'CEL') THEN                                      !
!
          E1_IN  = E0_IN + E_RANGE                                  !
          E1_O1  = E0_O1 + E_RANGE                                  !
!
        ELSE IF(MODE == 'CIS') THEN                                 !
!
          E1_EX  = E0_EX + E_RANGE                                  !
          E1_O1  = E0_O1 + E_RANGE                                  !
!
        ELSE IF(MODE == 'CFS') THEN                                 !
!
          E1_EX  = E0_EX + E_RANGE                                  !
          E1_O1  = E0_O1 + E_RANGE                                  !
!
        END IF                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE BEAM_INI
!
!=======================================================================
!
      SUBROUTINE INIT_CHECK_ANGLES(SPECTRO)
!
!  This subroutine initializes and checks angle-dependent variables
!
!
!  Input variables :
!
!                       SPECTRO   :  type of spectroscopy:
!                                       ---> 'PED' : photoelectron diffraction
!                                            'LED' : low-energy electron diffraction
!                                            'XAS' : X-ray absorption spectroscopy
!                                            'AED' : Auger electron diffraction
!                                            'APC' : Auger photoelectron coincidence spectroscopy
!                                            'E2E' : (e,2e) spectroscopy
!                                            'E3E' : (e,3e) spectroscopy
!                                            'ELS' : electron energy loss spectroscopy
!                                            'RES' : resonant elastic X-ray scattering
!                                            'STM' : scanning tunneling microscopy
!                                            'PLS' : photoelectron energy loss spectroscopy
!                                            'BEM' : ballistic energy electron microscopy
!                                            'EIG' : eigenvalue calculation
!                       STEREO    :  switch for stereographic projection: 'YES',' NO'
!
!  Output variables :
!
!                       *         :  label return for error in input data file
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 28 May 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO,FOUR,EIGHT
!
!
      USE CURRENT_AVER
      USE CURRENT_CALC
      USE CURRENT_FINA_VAL
      USE CURRENT_INIT_VAL
      USE CURRENT_FIXSCAN
      USE INIT_A
      USE INIT_L
      USE XAS
!
      USE AVER_IN
      USE ENERGY_IN
      USE PHI_IN
      USE THETA_IN
!
      USE ENERGY_EX
!
      USE AVER_O1
      USE ENERGY_O1
      USE PHI_O1
      USE THETA_O1
!
      USE AVER_O2
      USE ENERGY_O2
      USE PHI_O2
      USE THETA_O2
!
      USE AVER_O3
      USE ENERGY_O3
      USE PHI_O3
      USE THETA_O3
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  SPECTRO
!
      REAL (WP)             ::  SMALL
      REAL (WP)             ::  DPH,DTH
      REAL (WP)             ::  DPH_A,DTH_A
!
      REAL (WP)             ::  ABS
!
      DATA SMALL / 0.0001E0_WP /
!
!  When the direction of the analyzer might be experimentally
!    inaccurate, the calculation will be done for nine
!    direction across the one given in the data file
!    with an increment of one degree.
!
      IF(ICHKDIR == 1) THEN                                         !
        IF((ITHETA == 1) .AND. (IPHI == 0)) THEN                    !
          NPHI = 9                                                  !
          PHI0 = PHI0 - FOUR                                        !
          PHI1 = PHI0 + EIGHT                                       !
        ELSE IF((IPHI == 1) .AND. (ITHETA == 0)) THEN               !
          NTHETA = 9                                                !
          THETA0 = THETA0 - FOUR                                    !
          THETA1 = THETA0 + EIGHT                                   !
        END IF                                                      !
      END IF                                                        !
!
!  Initialization of the values for the scanned angle
!                 and the "fixed" one
!
      IF(IPHI == 1) THEN                                            !
        N_FIXED = NTHETA                                            !
        N_SCAN  = NPHI                                              !
        FIX0    = THETA0                                            !
        FIX1    = THETA1                                            !
        SCAN0   = PHI0                                              !
        SCAN1   = PHI1                                              !
        IPH_1   = 0                                                 !
      ELSE IF(ITHETA == 1) THEN                                     !
        N_FIXED = NPHI                                              !
        N_SCAN  = NTHETA                                            !
        FIX0    = PHI0                                              !
        FIX1    = PHI1                                              !
        SCAN0   = THETA0                                            !
        SCAN1   = THETA1                                            !
        IPH_1   = 1                                                 !
      ELSE IF(IE == 1) THEN                                         !
        IF(NTHETA >= NPHI) THEN                                     !
          N_FIXED = NPHI                                            !
          N_SCAN  = NTHETA                                          !
          FIX0    = PHI0                                            !
          FIX1    = PHI1                                            !
          SCAN0   = THETA0                                          !
          SCAN1   = THETA1                                          !
          IPH_1   = 1                                               !
        ELSE                                                        !
          N_FIXED = NTHETA                                          !
          N_SCAN  = NPHI                                            !
          FIX0    = THETA0                                          !
          FIX1    = THETA1                                          !
          SCAN0   = PHI0                                            !
          SCAN1   = PHI1                                            !
          IPH_1   = 0                                               !
        END IF                                                      !
      END IF                                                        !
!
      IF(SPECTRO == 'XAS') THEN                                     !
        NE   = NE_X                                                 !
      END IF                                                        !
!
      END SUBROUTINE INIT_CHECK_ANGLES
!
!=======================================================================
!
      SUBROUTINE ATDATA
!
!  This routine contains the atomic mass and the density of all the
!    elements,and the equivalence between their atomic number and
!    chemical symbol.
!
!  Value Z = 0 added for empty spheres. The values entered in this
!    case are arbitrary and set to the corresponding Z = 1 value
!    divided by 1836 (the ratio of the mass of the proton and electron).
!
!                                           Last modified :  6 May 2021
!
      USE ATOMIC_PROPERTIES
      USE XM_RHO
!
      IMPLICIT NONE
!
      INTEGER               ::  J
!
      DO J = 0, 105                                                 !
        XM_AT(J)  = MASS_AT(J)                                      !
        RHO_AT(J) = ATOM_DE(J)                                      !
      END DO                                                        !
!
      END SUBROUTINE ATDATA
!
!=======================================================================
!
      SUBROUTINE INIT_TEST_CALC(SPECTRO,EXCITATION,INTERACT_P,INTERACT_E)
!
!  This subroutine initializes the variables used for test calculations
!
!
!  Input variables :
!
!                       SPECTRO   :  type of spectroscopy:
!                                       ---> 'PED' : photoelectron diffraction
!                                            'LED' : low-energy electron diffraction
!                                            'XAS' : X-ray absorption spectroscopy
!                                            'AED' : Auger electron diffraction
!                                            'APC' : Auger photoelectron coincidence spectroscopy
!                                            'E2E' : (e,2e) spectroscopy
!                                            'E3E' : (e,3e) spectroscopy
!                                            'ELS' : electron energy loss spectroscopy
!                                            'RES' : resonant elastic X-ray scattering
!                                            'STM' : scanning tunneling microscopy
!                                            'PLS' : photoelectron energy loss spectroscopy
!                                            'BEM' : ballistic energy electron microscopy
!                                            'EIG' : eigenvalue calculation
!                       EXCITATION:  way the electron is created:
!                                       ---> 'NOINTER' : electron already exists
!                                       ---> 'PH-ELEC' : electron excited by photon
!                                       ---> 'COUL-LL' : Coulomb with two electrons localized
!                                       ---> 'COUL-DL' : Coulomb with one electron delocalized
!                                       ---> 'COUL-DD' : Coulomb with two electrons delocalized
!
!  Output variables :
!
!                       INTERACT_P:  excitation of the photon beam
!                       INTERACT_E:  excitation of the electron beam
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 28 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE
!
      USE CURRENT_CALC
      USE CURRENT_INIT_VAL
      USE CLUSTER
      USE CLUSTER_TIP
      USE EXP_TYPE,          ONLY : IMOD
      USE INIT_A
      USE INIT_L
!
      USE ATOMIC_INDEX
      USE SELECT_RULE
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 3)   ::  SPECTRO
      CHARACTER (LEN = 5)   ::  EXCITATION(5)
      CHARACTER (LEN = 7)   ::  INTERACT_P,INTERACT_E
!
      INTEGER               ::  JEL
!
      IF(I_TEST == 1) THEN                                          !
        INTERACT_P = 'NOINTER'                                      !
        INTERACT_E = 'NOINTER'                                      !
        DO JEL = 1, 5                                               !
          IF(EXCITATION(JEL) == 'PH-ELEC') THEN                     !
            INTERACT_P = 'PH-ELEC'                                  !
          END IF                                                    !
          IF(EXCITATION(JEL) == 'COUL-LL') THEN                     !
            INTERACT_E = 'COUL-LL'                                  !
          END IF                                                    !
          IF(EXCITATION(JEL) == 'COUL-DL') THEN                     !
            INTERACT_E = 'COUL-DL'                                  !
          END IF                                                    !
          IF(EXCITATION(JEL) == 'COUL-DD') THEN                     !
            INTERACT_E = 'COUL-DD'                                  !
          END IF                                                    !
        END DO
        IF(INTERACT_P == 'PH-ELEC') THEN                            !
          LI = 0                                                    !
          IF(SPECTRO /= 'RES') THEN                                 !
            SELRULE   = '+1>E1'                                     !
            IPOL      = 1                                           !
          ELSE                                                      !
            SELRULE_I = '+1>E1'                                     !
            SELRULE_O = '+1>E1'                                     !
          END IF                                                    !
        END IF                                                      !
        IF(INTERACT_E == 'COUL-LL') THEN                            !
          LI_C = 0                                                  !
          LI_I = 0                                                  !
        END IF                                                      !
      END IF                                                        !
!
      IF(I_TEST == 2) THEN                                          !
        IF(ABS(IPOL) == 1) THEN                                     !
          THLUM  = - 90.0E0_WP                                      !
          PHILUM = ZERO                                             !
        ELSE IF(ABS(IPOL) == 2) THEN                                !
          THLUM  = ZERO                                             !
          PHILUM = ZERO                                             !
        END IF                                                      !
        IMOD    = 0                                                 !
        VINT    = ZERO                                              !
        VINT_TI = ZERO                                              !
        A       = ONE                                               !
      END IF                                                        !
!
!..........  Atomic case index  ..........
!
       IF(I_TEST == 2) THEN                                         !
         I_AT = 1                                                   !
       END IF                                                        !
!
      END SUBROUTINE INIT_TEST_CALC
!
!=======================================================================
!
      SUBROUTINE EXTERNAL_ANGLES(SPECTRO,IRET)
!
!  This subroutine initializes the switch controlling the external reading
!   of the detector directions (for angular averaging of an undetected
!   electron)
!
!
!  Input variables :
!
!                       SPECTRO   :  type of spectroscopy:
!                                       ---> 'PED' : photoelectron diffraction
!                                            'LED' : low-energy electron diffraction
!                                            'XAS' : X-ray absorption spectroscopy
!                                            'AED' : Auger electron diffraction
!                                            'APC' : Auger photoelectron coincidence spectroscopy
!                                            'E2E' : (e,2e) spectroscopy
!                                            'E3E' : (e,3e) spectroscopy
!                                            'ELS' : electron energy loss spectroscopy
!                                            'RES' : resonant elastic X-ray scattering
!                                            'STM' : scanning tunneling microscopy
!                                            'PLS' : photoelectron energy loss spectroscopy
!                                            'BEM' : ballistic energy electron microscopy
!                                            'EIG' : eigenvalue calculation
!
!  Output variables :
!
!                       IRET      :  error return code
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 27 May 2021
!
!
      USE CURRENT_AVER
      USE CURRENT_CALC
      USE CURRENT_FIXSCAN
      USE CURRENT_EXT_DIR
!
      USE AVER_O1
      USE AVER_O2
      USE EXTERNAL_O1
      USE EXTERNAL_O2
      USE INTEGR_O1
      USE INTEGR_O2
      USE PHI_O1
      USE PHI_O2
      USE THETA_O1
      USE THETA_O2
!
      USE AVERAGING
      USE INFILES
      USE INUNITS
      USE OUTUNITS
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 3)   ::  SPECTRO
!
      INTEGER, INTENT(OUT)  ::  IRET
!
      INTEGER               ::  JS,NSET_O1,NSET_O2
      INTEGER               ::  IDIR_O1,I_SET_O1,N_POINTS_O1
      INTEGER               ::  IDIR_O2,I_SET_O2,N_POINTS_O2
      INTEGER               ::  I_PH_O1,N_FIXED_O1,N_SCAN_O1
      INTEGER               ::  I_PH_O2,N_FIXED_O2,N_SCAN_O2
!
      IF(SPECTRO == 'APC') THEN                                     !
        IF((I_EXT_O1 == -1) .OR. (I_EXT_O2 == -1)) THEN             !
          IF(I_EXT_O1 * I_EXT_O2 == 0) THEN                         !
            WRITE(IUO1,523)                                         !
            I_EXT_O1 = - 1                                          !
            I_EXT_O2 = - 1                                          !
            OPEN(UNIT=IUI11, FILE = INFILE11, STATUS='OLD')         !
            OPEN(UNIT=IUI14, FILE = INFILE14, STATUS='OLD')         !
            READ(IUI11,713) IDIR_O1,NSET_O1                         !
            READ(IUI14,713) IDIR_O2,NSET_O2                         !
            IF(IDIR_O1 == 2) THEN                                   !
              IF(NSET_O1 /= NSET_O2) WRITE(IUO1,524) NSET_O1,NSET_O2!
              STOP                                                  !
            END IF                                                  !
          END IF                                                    !
        END IF                                                      !
        IF(I_INT_O1 == 1) THEN                                      !
          I_EXT_O1 = 2                                              !
        ELSE IF(I_INT_O1 == 2) THEN                                 !
          I_EXT_O2 = 2                                              !
        ELSE IF(I_INT_O1 == 3) THEN                                 !
          I_EXT_O1 = 2                                              !
          I_EXT_O2 = 2                                              !
        END IF                                                      !
      END IF                                                        !
!
      IF(I_EXT_O1 == -1) THEN                                       !
        OPEN(UNIT=IUI11, FILE = INFILE11, STATUS='OLD')             !
        READ(IUI11,701) IDIR_O1,I_SET_O1,N_POINTS_O1                !
        READ(IUI11,702) I_PH_O1,N_FIXED_O1,N_SCAN_O1                !
        DO JS = 1, I_SET_O1                                         !
          READ(IUI11,703) TH_0(JS),PH_0(JS)                         !
        END DO                                                      !
        CLOSE(IUI11)                                                !
        IF(IDIR_O1 /= 2) IRET = 12                                  !
        IF(I_PH_O1 /= IPH_1) IPH_1 = I_PH_O1                        !
        IF((SPECTRO == 'PED') .OR. (SPECTRO == 'APC')) THEN         !
          IF(I_PH_O1 == 0) THEN                                     !
            NTHETA = N_FIXED                                        !
            NPHI   = N_SCAN                                         !
          ELSE                                                      !
            NTHETA = N_SCAN                                         !
            NPHI   = N_FIXED                                        !
          END IF                                                    !
          ICHKDIR = 2                                               !
        END IF                                                      !
      END IF                                                        !
      IF(I_EXT_O1 >= 1) THEN                                        !
        OPEN(UNIT=IUI11, FILE = INFILE11, STATUS='OLD')             !
        READ(IUI11,701) IDIR,I_SET,N_POINTS                         !
        CLOSE(IUI11)                                                !
        IF((IDIR /= 1) .AND. (I_EXT == 2)) IRET = 12                !
        N_FIXED = N_POINTS                                          !
        N_SCAN  = 1                                                 !
        NTHETA  = N_POINTS                                          !
        NPHI    = 1                                                 !
      END IF                                                        !
      IF(I_EXT_O2 >= 1) THEN                                        !
        IF(SPECTRO == 'APC') THEN                                   !
          OPEN(UNIT=IUI14, FILE = INFILE14, STATUS='OLD')           !
          READ(IUI14,701) IDIR_O2,I_SET_O2,N_POINTS_O2              !
          CLOSE(IUI14)                                              !
        ELSE                                                        !
          OPEN(UNIT=IUI11, FILE = INFILE11, STATUS='OLD')           !
          READ(IUI11,701) IDIR_O2,I_SET_O2,N_POINTS_O2              !
          CLOSE(IUI11)                                              !
        END IF                                                      !
        IF((IDIR_O2 /= 1) .AND. (I_EXT_O2 == 2)) IRET = 12          !
        N_FIXED_O2 = N_POINTS_O2                                    !
        N_SCAN_O2  = 1                                              !
        NTHETA_O2  = N_POINTS_O2                                    !
        NPHI_O2    = 1                                              !
      END IF                                                        !
!
      IF(I_EXT_O1 == -1) THEN                                       !
        IF(SPECTRO == 'APC') THEN                                   !
          OPEN(UNIT=IUI14, FILE = INFILE14, STATUS='OLD')           !
          READ(IUI14,701) IDIR_O2,I_SET_O2,N_POINTS_O2              !
          READ(IUI14,702) I_PH_O2,N_FIXED_O2,N_SCAN_O2              !
          CLOSE(IUI14)                                              !
        ELSE                                                        !
          OPEN(UNIT=IUI11, FILE = INFILE11, STATUS='OLD')           !
          READ(IUI11,701) IDIR_O2,I_SET_O2,N_POINTS_O2              !
          READ(IUI11,702) I_PH_O2,N_FIXED_O2,N_SCAN_O2              !
          CLOSE(IUI11)                                              !
        END IF                                                      !
        IF(IDIR_O2 /= 2) IRET = 12                                  !
        IF(I_PH_O2 == 0) THEN                                       !
          NTHETA_O2 = N_FIXED_O2                                    !
          NPHI_O2   = N_SCAN_O2                                     !
        ELSE                                                        !
          NTHETA_O2 = N_SCAN_O2                                     !
          NPHI_O2   = N_FIXED_O2                                    !
        END IF                                                      !
        ICHKDIR_O2 = 2                                              !
      END IF                                                        !
!
!  Formats
!
 523  FORMAT(///,4X,' <<<<<<<<<<  BOTH DETECTOR DIRECTIONS MUST BE ', &
             'EITHER INTERNAL OR EXTERNAL  >>>>>>>>>>',/,8X,          &
                    ' -----> PROCEEDING WITH EXTERNAL DIRECTIONS',/)
 524  FORMAT(///,4X,' <<<<<<<<<<  AVERAGING OVER ',I3,' DOMAINS ',    &
             'FOR PHOTOELECTRON   >>>>>>>>>>',/,4X,                   &
             ' <<<<<<<<<<  AVERAGING OVER ',I3,' DOMAINS ',           &
             'FOR AUGER ELECTRON   >>>>>>>>>>',/,8X,                  &
             ' -----> IMPOSSIBLE : CHECK INPUT FILES !')
!
 701  FORMAT(6X,I1,1X,I3,2X,I4)
 702  FORMAT(6X,I1,1X,I3,3X,I3)
 703  FORMAT(15X,F8.3,3X,F8.3)
 713  FORMAT(6X,I1,1X,I3)
!
      END SUBROUTINE EXTERNAL_ANGLES
!
!=======================================================================
!
      SUBROUTINE INIT_ALGO(NATCLU)
!
!  This subroutine initializes what is needed for the different algorithms
!    that compute the multiple scattering problem.
!
!  Input variables :
!
!                       NATCLU    :  number of atoms in the cluster
!
!
!  Varia           :
!                       ALGO(J)   :  two-character string representing
!                       ALGO_XX      the algorithm used to compute the
!                                    multiple scattering problem for
!                                    the electron J (XX) considered:
!
!                                    ---> 'MI' : matrix inversion
!                                    ---> 'CE' : correlation expansion
!                                    ---> 'RE' : renormalized expansion
!                                    ---> 'SE' : Rehr-Albers series expansion
!                                    ---> 'SM' : series matrix expansion
!
!                                   with:
!
!                                   ---> J = 1 : incoming electron
!                                   ---> J = 2 : excited electron (not detected)
!                                   ---> J = 3 : first outgoing electron
!                                   ---> J = 4 : second outgoing electron
!                                   ---> J = 5 : third outgoing electron
!
!
!
!   Author: Didier Sébilleau
!
!                                        Last modified : 10 Jun 2021
!
!
      USE REAL_NUMBERS,            ONLY : ONE
!
      USE ALGORITHMS
      USE BEAMS_ALGO
      USE CURRENT_LINLBD
      USE PATH_INFO
!
      USE APPROXIMATIONS_IN,       ONLY : N_SCAT_IN,NO_IN
      USE APPROXIMATIONS_EX,       ONLY : N_SCAT_EX,NO_EX
      USE APPROXIMATIONS_O1,       ONLY : N_SCAT_O1,NO_O1
      USE APPROXIMATIONS_O2,       ONLY : N_SCAT_O2,NO_O2
      USE APPROXIMATIONS_O3,       ONLY : N_SCAT_O3,NO_O3
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  NATCLU
!
      INTEGER               ::  I_SE,I_MI,I_CE,I_RE,I_SM
      INTEGER               ::  J,LAMBDA0,N_O,NMX,NU,MU,NMU
      INTEGER               ::  N_SCAT,NO
!
      INTEGER               ::  MAX,ABS
!
      N_SCAT = MAX(N_SCAT_IN,N_SCAT_EX,N_SCAT_O1,N_SCAT_O2,N_SCAT_O3)
      NO     = MAX(NO_IN,NO_EX,NO_O1,NO_O2,NO_O3)                   !
!
      ALGO(1) = ALGO_IN                                             !
      ALGO(2) = ALGO_EX                                             !
      ALGO(3) = ALGO_O1                                             !
      ALGO(4) = ALGO_O2                                             !
      ALGO(5) = ALGO_O3                                             !

      I_SE = 0                                                      !
      I_MI = 0                                                      !
      I_CE = 0                                                      !
      I_RE = 0                                                      !
      I_SM = 0                                                      !

      DO J = 1, 5                                                   !
        IF(ALGO(J) == 'CE') I_CE = I_CE + 1                         ! at least one electron using CE
        IF(ALGO(J) == 'MI') I_MI = I_MI + 1                         ! at least one electron using MI
        IF(ALGO(J) == 'RE') I_RE = I_RE + 1                         ! at least one electron using RE
        IF(ALGO(J) == 'SE') I_SE = I_SE + 1                         ! at least one electron using SE
        IF(ALGO(J) == 'SM') I_SM = I_SM + 1                         ! at least one electron using SM
      END DO                                                        !

      IF(I_SE >= 1) THEN                                            ! series expansion case

        NPATH2(0) = ONE                                             !
        NPATH(0)  = ONE                                             !
        NPMA(0)   = ONE                                             !
        NPMI(0)   = ONE                                             !
!
!..........     Construction of the linear index LAMBDA=(MU,NU)     ..........
!
        LAMBDA0 = 0                                                 !
        DO N_O = 0, NO                                              !
          NMX = N_O / 2                                             !
          DO NU = 0, NMX                                            !
            DO MU = -N_O, N_O                                       !
              NMU = 2  *NU + ABS(MU)                                !
              IF(NMU == N_O) THEN                                   !
                LAMBDA0    = LAMBDA0 + 1                            !
                LBD(MU,NU) = LAMBDA0                                !
              END IF                                                !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
        LBDMAX = LAMBDA0                                            !
!
      END IF                                                        !
!
      IF(I_MI >= 1) THEN                                            ! matrix inversion case
!
        CONTINUE                                                    !
!
      END IF                                                        !
!
      IF(I_CE >= 1) THEN                                            ! correlation expansion case
!
        CALL COEFPQ(NATCLU,N_SCAT)                                  !
!
      END IF                                                        !
!
      IF(I_RE >= 1) THEN                                            ! renormalized expansion case
!
        CONTINUE                                                    !
!
      END IF                                                        !
!
      END SUBROUTINE INIT_ALGO
!
!=======================================================================
!
      SUBROUTINE CALC_LOG_GAMMA(I_SPHER)
!
!  This subroutine computes the logarithm of the Gamma function
!    and stores the result in common blocks. It also computes
!            various other coefficients of use
!
!
!  Input variables :
!
!                       I_SPHER   :  maximum of the I_BASIS values
!
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 31 May 2021
!
!
      USE REAL_NUMBERS,         ONLY : ZERO,ONE,HALF
      USE PI_ETC,               ONLY : SQR_PI
!
      USE INIT_M
!
      USE LOGAMMA
      USE STORE_COEF
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  I_SPHER
!
      INTEGER               ::  I,J,L,M1,M2
      INTEGER               ::  NFAC
!
      REAL (WP)             ::  FACT1L,FACT2L
      REAL (WP)             ::  XDEN
!
      REAL (WP)             ::  LOG,EXP,FLOAT,SQRT
!
      NFAC = N_GAUNT                                                !
!
!  Storage of the logarithm of the Gamma function GLD(N+1,N_INT)
!  for integer (N_INT=1) and semi-integer (N_INT=2) values :
!
!    GLD(N+1,1)   =   Log(N!) for N integer
!    GLD(N+1/2,2) =   Log(N!) for N semi-integer
!
      IF((I_SPHER >= 0) .OR. (I_MULT == 1)) THEN                    !
        GLG(1)   = ZERO                                             !
        GLD(1,1) = ZERO                                             !
        GLD(1,2) = LOG(SQR_PI * HALF)                               !
        DO I = 2, NFAC                                              !
          J        = I - 1                                          !
          GLG(I)   = GLG(J)   + LOG(FLOAT(J))                       !
          GLD(I,1) = GLD(J,1) + LOG(FLOAT(J))                       !
          GLD(I,2) = GLD(J,2) + LOG(FLOAT(J) + HALF)                !
        END DO                                                      !
      ELSE                                                          !
        GLG(1) = ZERO                                               !
        DO I = 2, NFAC                                              !
          J      = I - 1                                            !
          GLG(I) = GLG(J) + LOG(FLOAT(J))                           !
        END DO                                                      !
      END IF                                                        !
!
      EXPF(0,0) = ONE                                               !
      EXPR(0,0) = ONE                                               !
      FACT1L    = ZERO                                              !
      DO L = 1, 2 * L_MAX                                           !
        XDEN   = ONE / SQRT(FLOAT(L+L+1))                           !
        FACT1L = FACT1L + LOG(FLOAT(L))                             !
        FACT2L = LOG(FLOAT(L+1))                                    !
        DO M1 = 0, L                                                !
          EXPF(M1,L)  = EXP(HALF * ( GLG(L+M1+1) - GLG(L-M1+1)) )   !
          EXPR(M1,L)  = EXP(HALF * ( GLG(L+L+1) - GLG(L+M1+1) -   & !
                                     GLG(L-M1+1) ))                 !
          EXPF2(L,M1) = EXPF(M1,L) * XDEN                           !
          IF(M1 > 0) THEN                                           !
            FACT2L = FACT2L + LOG(FLOAT(1+L+M1))                    !
          END IF                                                    !
          IF(L < NL_M) THEN                                         !
            DO M2 = 0, L                                            !
              CF(L,M1,M2) = SQRT( FLOAT( (L*L-M1*M1) *            & !
                                         (L*L-M2*M2) )            & !
                                ) / FLOAT(L)                        !
            END DO                                                  !
          END IF                                                    !
        END DO                                                      !
        FSQ(L) = EXP( HALF * (FACT2L - FACT1L) )                    !
      END DO                                                        !
!
      END SUBROUTINE CALC_LOG_GAMMA
!
!=======================================================================
!
      SUBROUTINE CHOOSE_COUMAT(SPECTRO,N_CALL)
!
!  This routine selects the type of optical matrix elements that will
!     be computed. It also sets up the angular momentum selection rules.
!
!  Input variables :
!
!                       SPECTRO   :  a 3-character string representing the
!                                      type of spectroscopy considered
!                       N_CALL    :  index for the number of the call to
!                                      this subroutine (= 1 for incoming beam
!                                      and = 2 for outgoing beam)
!
!  Important parameter :
!
!                       SELRULE   :  a 5-character string indicating
!                                      the type of coupling approximation
!                                      for the photon-electron excitation
!
!                                    it can be either a binary string,
!                                    (for instance '11000' when 1/0 indicates
!                                    that the corresponding approximation is
!                                    activated/disactivated with the order
!                                    E1 E2 E3 M1 M2. In the previous example,
!                                    only E1 and E2 are activated) or a string
!                                    of the form '+1>E1' which indicates that
!                                    only the LI+1 channel of E1 will be computed
!                                    (single channel case). The case '00000'
!                                    will also correspond to a single channel,
!                                    but with NO radial integral taken into
!                                    account (by contrast to '+1>E1')
!
!                       M_R       :  0 to 3 (the absolute value of max|MF-MI|)
!
!  Types of matrix elements calculated:
!
!              * I_C1 = 1 : single channel
!              * I_E1 = 1 : electric dipole E1
!              * I_E1 = 2 : electric dipole E1 + 3rd order correction term from exponential
!              * I_E2 = 1 : electric quadrupole E2
!              * I_E3 = 1 : electric octupole E3
!              * I_M1 = 1 : magnetic dipole M1
!              * I_M2 = 1 : magnetic quadrupole M2 + polar toroidal dipole T1
!
!
!
!   Author :  D. Sébilleau
!
!                                     Last modified :  6 May 2021
!
!
      USE INIT_L
      USE INIT_L_I
      USE INIT_L_O
      USE M_RULE
      USE ONE_CHANNEL
      USE ONE_CHANNEL_I
      USE ONE_CHANNEL_O
      USE OPTICAL_ELTS
      USE OPTICAL_ELTS_I
      USE OPTICAL_ELTS_O
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  SPECTRO
!
      INTEGER, INTENT(IN)   ::  N_CALL
!
      INTEGER               ::  JF
!
!  Initialization of switches to select the type of matrix elements calculated
!
      I_C1 = 0                                                      ! single channel
      I_E1 = 0                                                      ! electric dipole E1
      I_E2 = 0                                                      ! electric quadrupole E2
      I_E3 = 0                                                      ! electric octupole E3
      I_M1 = 0                                                      ! magnetic dipole M1
      I_M2 = 0                                                      ! magnetic quadrupole M2 + T1
!
      IF(SELRULE(3:3) == '>') THEN                                  !
!
!  Single channel case
!
        I_C1    = 1                                                 !
        CHANNEL = SELRULE(1:2)                                      !
        TYPCHAN = SELRULE(4:5)                                      !
!
      ELSE                                                          !
!
!  Standard approximations: E1, E2, E3, M1, M2+T1
!
        IF(SELRULE(1:1) == '1') I_E1 = 1                            !
        IF(SELRULE(1:1) == '2') I_E1 = 2                            !
        IF(SELRULE(2:2) == '1') I_E2 = 1                            !
        IF(SELRULE(3:3) == '1') I_E3 = 1                            !
        IF(SELRULE(4:4) == '1') I_M1 = 1                            !
        IF(SELRULE(5:5) == '1') I_M2 = 1                            !
        CHANNEL = '  '                                              !
        TYPCHAN = '  '                                              !
      END IF                                                        !
!
!  Alternative single channel case (with no radial integral)
!
      IF(SELRULE == '00000') THEN                                   !
        I_C1    = 1                                                 !
        CHANNEL = '  '                                              !
        TYPCHAN = '  '                                              !
      END IF                                                        !
!
!  Setting up of the angular momentum selection rules
!
      ISTEP_LF = 1                                                  !
      I_0      = 0                                                  ! switch to add LF = LI
!                                                                   ! to selection rules
      IF(I_C1 == 1) THEN                                            !
!
!.....  Single channel case  .....
!
        IF(CHANNEL == '-3') THEN                                    !
          LF1 = LI - 3                                              !
          M_R = 3                                                   !
        ELSE IF(CHANNEL == '-2') THEN                               !
          LF1 = LI - 2                                              !
          M_R = 2                                                   !
        ELSE IF(CHANNEL == '-1') THEN                               !
          LF1 = LI - 1                                              !
          M_R = 1                                                   !
        ELSE IF(CHANNEL == '+1') THEN                               !
          LF1 = LI + 1                                              !
          M_R = 1                                                   !
        ELSE IF(CHANNEL == '+2') THEN                               !
          LF1 = LI + 2                                              !
          M_R = 2                                                   !
        ELSE IF(CHANNEL == '+3') THEN                               !
          LF1 = LI + 3                                              !
          M_R = 3                                                   !
        ELSE                                                        ! case LF = LI
          LF1 = LI                                                  !
          M_R = 0                                                   !
        END IF                                                      !
!
        IF(SELRULE == '00000') THEN                                 ! single channel case
          LF1 = LI                                                  ! with no radial integral
          M_R = 0                                                   ! (matrix elements set to 1)
        END IF                                                      !
!
        LF2 = LF1                                                   !
!
      ELSE                                                          !
!
!..... Multichannel case  .....
!
        IF(SELRULE == '10000') THEN                                 ! E1 case
          LF1      = LI - 1                                         !
          LF2      = LI + 1                                         !
          ISTEP_LF = 2                                              !
          M_R      = 1                                              !
        ELSE IF(SELRULE == '20000') THEN                            ! E1 case
          LF1      = LI - 1                                         ! with 3rd order
          LF2      = LI + 1                                         ! correction
          ISTEP_LF = 2                                              !
          M_R      = 1                                              !
        ELSE IF(SELRULE == '01000') THEN                            ! E2 case
          LF1      = LI - 2                                         !
          LF2      = LI + 2                                         !
          ISTEP_LF = 2                                              !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '00100') THEN                            ! E3 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          ISTEP_LF = 2                                              !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '00010') THEN                            ! M1 case
          LF1      = LI                                             !
          LF2      = LF1                                            !
          M_R      = 1                                              !
        ELSE IF(SELRULE == '00001') THEN                            ! M2 case
          LF1      = LI - 1                                         !
          LF2      = LI + 1                                         !
          ISTEP_LF = 2                                              !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '11000') THEN                            ! E1+E2 case
          LF1      = LI - 2                                         !
          LF2      = LI + 2                                         !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '10100') THEN                            ! E1+E3 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          ISTEP_LF = 2                                              !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '10010') THEN                            ! E1+M1 case
          LF1      = LI - 1                                         !
          LF2      = LI + 1                                         !
          M_R      = 1                                              !
        ELSE IF(SELRULE == '10001') THEN                            ! E1+M2 case
          LF1      = LI - 1                                         !
          LF2      = LI + 1                                         !
          ISTEP_LF = 2                                              !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '01100') THEN                            ! E2+E3 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '01010') THEN                            ! E2+M1 case
          LF1      = LI - 2                                         !
          LF2      = LI + 2                                         !
          ISTEP_LF = 2                                              !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '01001') THEN                            ! E2+M2 case
          LF1      = LI - 2                                         !
          LF2      = LI + 2                                         !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '00110') THEN                            ! E3+M1 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          I_0      = 1                                              ! add selection rule LF = LI
          ISTEP_LF = 2                                              !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '00101') THEN                            ! E3+M2 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          ISTEP_LF = 2                                              !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '11100') THEN                            ! E1+E2+E3 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '11010') THEN                            ! E1+E2+M1 case
          LF1      = LI - 2                                         !
          LF2      = LI + 2                                         !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '11001') THEN                            ! E1+E2+M2 case
          LF1      = LI - 2                                         !
          LF2      = LI + 2                                         !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '01110') THEN                            ! E2+E3+M1 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '01101') THEN                            ! E2+E3+M2 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '00111') THEN                            ! E3+M1+M2 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          I_0      = 1                                              ! add selection rule LF = LI
          ISTEP_LF = 2                                              !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '11110') THEN                            ! E1+E2+E3+M1 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '11101') THEN                            ! E1+E2+E3+M2 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '01111') THEN                            ! E2+E3+M1+M2 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '10111') THEN                            ! E1+E3+M1+M2 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          I_0      = 1                                              ! add selection rule LF = LI
          ISTEP_LF = 2                                              !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '11011') THEN                            ! E1+E2+M1+M2 case
          LF1      = LI - 2                                         !
          LF2      = LI + 2                                         !
          M_R      = 2                                              !
        ELSE IF(SELRULE == '11101') THEN                            ! E1+E2+E3+M1 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '11110') THEN                            ! E1+E2+E3+M2 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        ELSE IF(SELRULE == '11111') THEN                            ! E1+E2+E3+M1+M2 case
          LF1      = LI - 3                                         !
          LF2      = LI + 3                                         !
          M_R      = 3                                              !
        END IF                                                      !
      END IF                                                        !
!
!  Change of initial value LF1 when negative
!
      IF(LF1 < 0) THEN                                              !
        DO JF = 1, 10                                               !
          LF1 = LF1 + ISTEP_LF                                      !
          IF(LF1 >= 0) GO TO 100                                    !
        END DO                                                      !
 100    CONTINUE                                                    !
      END IF                                                        !
!
!  REXS case: making the difference between incoming and outgoing beam
!
      IF(SPECTRO == 'RES') THEN                                     !
!
        IF(N_CALL == 1) THEN                                        !
!
          I_C1_I     = I_C1                                         !
          I_E1_I     = I_E1                                         !
          I_E2_I     = I_E2                                         !
          I_E3_I     = I_E3                                         !
          I_M1_I     = I_M1                                         !
          I_M2_I     = I_M2                                         !
!
          CHANNEL_I  = CHANNEL                                      !
          TYPCHAN_I  = TYPCHAN                                      !
!
          LF1_I      = LF1                                          !
          LF2_I      = LF2                                          !
          ISTEP_LF_I = ISTEP_LF                                     !
          I_0_I      = I_0                                          !
!
        ELSE
          I_C1_O     = I_C1                                         !
          I_E1_O     = I_E1                                         !
          I_E2_O     = I_E2                                         !
          I_E3_O     = I_E3                                         !
          I_M1_O     = I_M1                                         !
          I_M2_O     = I_M2                                         !
!
          CHANNEL_O  = CHANNEL                                      !
          TYPCHAN_O  = TYPCHAN                                      !
!
          LF1_O      = LF1                                          !
          LF2_O      = LF2                                          !
          ISTEP_LF_O = ISTEP_LF                                     !
          I_0_O      = I_0                                          !
        END IF                                                      !
      END IF                                                        !
!
      END SUBROUTINE CHOOSE_COUMAT
!
!=======================================================================
!
      SUBROUTINE INIT_SPIN_ORBIT(I_SPIN,IRET,*)
!
!  This subroutine initializes spin-orbit parameters
!
!
!  Input variables :
!
!                       I_SPIN    :  switch for spin-resolved calculation

!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 10 May 2021
!
!
      USE REAL_NUMBERS,   ONLY : ZERO
!
      USE CLEBSCH_GORDAN
      USE HEADER
      USE INIT_J
      USE INIT_L
      USE INIT_L_I
      USE INIT_L_O
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  I_SPIN
      INTEGER, INTENT(OUT)  ::  IRET
!
      INTEGER               ::  JS,JI,MJI
!
      REAL (WP)             ::  FLOAT,SQRT
!
!  Initialization of the values of ji if spin-orbit is taken
!                  into account.
!
!  Here :   JI is the loop index going from JF1 to JF2 with :
!
!                   JI=1    : ji = li + 1/2
!                   JI=2    : ji = li - 1/2
!
      IF(I_SO <= 0) THEN                                            !
        JF1 = 1                                                     !
        JF2 = 2                                                     !
      ELSE IF(I_SO == 1) THEN                                       !
        IF(S_O == '1/2') THEN                                       !
          IF(LI == 0) THEN                                          !
            JF1 = 1                                                 !
            JF2 = 1                                                 !
          ELSE IF(LI == 1) THEN                                     !
            JF1 = 2                                                 !
            JF2 = 2                                                 !
          END IF                                                    !
        ELSE IF(S_O == '3/2') THEN                                  !
          IF(LI == 1) THEN                                          !
            JF1 = 1                                                 !
            JF2 = 1                                                 !
          ELSE IF(LI == 2) THEN                                     !
            JF1 = 2                                                 !
            JF2 = 2                                                 !
          END IF                                                    !
        ELSE IF(S_O == '5/2') THEN                                  !
          IF(LI == 2) THEN                                          !
            JF1 = 1                                                 !
            JF2 = 1                                                 !
          ELSE IF(LI == 3) THEN                                     !
            JF1 = 2                                                 !
            JF2 = 2                                                 !
          END IF                                                    !
        ELSE IF(S_O == '7/2') THEN                                  !
          IF(LI == 3) THEN                                          !
            JF1 = 1                                                 !
            JF2 = 1                                                 !
          ELSE IF(LI == 4) THEN                                     !
            JF1 = 2                                                 !
            JF2 = 2                                                 !
          END IF                                                    !
        ELSE IF(S_O == '9/2') THEN                                  !
          IF(LI == 4) THEN                                          !
            JF1 = 1                                                 !
            JF2 = 1                                                 !
          ELSE                                                      !
            IRET = 7                                                !
            RETURN 1                                                !
          END IF                                                    !
        ELSE                                                        !
          IRET = 7                                                  !
          RETURN 1                                                  !
        END IF                                                      !
      ELSE IF(I_SO == 2) THEN                                       !
        JF1 = 1                                                     !
        JF2 = 2                                                     !
      ELSE                                                          !
        IRET = 7                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
      IF(NI <= 5) THEN                                              !
         NNL = NI * (NI-1) / 2 + LI + 1                             !
      ELSE IF(NI == 6) THEN                                         !
         NNL = NI * (NI-1) / 2 + LI                                 !
      ELSE IF(NI == 7) THEN                                         !
         NNL = NI * (NI-1) / 2 + LI - 3                             !
      END IF                                                        !
      NNL1 = NNL                                                    !
      NNL2 = NNL                                                    !
!
!  Storage of the Clebsch-Gordan coefficients for the spin-orbit
!  dependent coupling matrix elements in the array CG(MJI,JI,JSPIN).
!
!  Here :           JI=1    : ji = li + 1/2
!                   JI=2    : ji = li - 1/2
!                   MJI     : mji + 1/2
!                   JSPIN=1 : msi = +1/2
!                   JSPIN=2 : msi = -1/2
!
!              so that all indices remain integer
!
      IF((I_SO > 0) .OR. (I_SPIN == 1)) THEN                        !
        DO JS = 1, 2                                                !
          DO JI = 1, 2                                              !
            DO MJI = - LI, LI + 1                                   !
              CG(MJI,JI,JS) = ZERO                                  !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
        DO MJI = - LI, LI + 1                                       !
          CG(MJI,1,1) = SQRT(FLOAT(LI+MJI)   / FLOAT(LI+LI+1))      !
          CG(MJI,1,2) = SQRT(FLOAT(LI-MJI+1) / FLOAT(LI+LI+1))      !
          IF((MJI > -LI) .AND. (MJI < LI+1)) THEN                   !
            CG(MJI,2,1) = - SQRT(FLOAT(LI-MJI+1) / FLOAT(LI+LI+1))  !
            CG(MJI,2,2) =   SQRT(FLOAT(LI+MJI)   / FLOAT(LI+LI+1))  !
          END IF                                                    !
        END DO                                                      !
      END IF                                                        !
!
      END SUBROUTINE INIT_SPIN_ORBIT
!
!=======================================================================
!
      SUBROUTINE INIT_MULTIPLET(SPECTRO)
!
!  This subroutine initializes Auger multiplet parameters
!
!
!  Input variables :
!
!                       SPECTRO   :  type of spectroscopy:
!                                       ---> 'PED' : photoelectron diffraction
!                                            'LED' : low-energy electron diffraction
!                                            'XAS' : X-ray absorption spectroscopy
!                                            'AED' : Auger electron diffraction
!                                            'APC' : Auger photoelectron coincidence spectroscopy
!                                            'E2E' : (e,2e) spectroscopy
!                                            'E3E' : (e,3e) spectroscopy
!                                            'ELS' : electron energy loss spectroscopy
!                                            'RES' : resonant elastic X-ray scattering
!                                            'STM' : scanning tunneling microscopy
!                                            'PLS' : photoelectron energy loss spectroscopy
!                                            'BEM' : ballistic energy electron microscopy
!                                            'EIG' : eigenvalue calculation
!
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 21 May 2021
!
!
      USE REAL_NUMBERS,     ONLY : ZERO,ONE,TWO,HALF,SMALL
!
      USE CLEBSCH_GORDAN
      USE INIT_A
      USE INIT_M
      USE OUTUNITS
      USE TEMP01
!
      USE ANGULAR_MOMENTUM, ONLY : N_J
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 3)   ::  SPECTRO
!
      INTEGER               ::  N,LJ_MAX
      INTEGER               ::  LJ1,LMJ1,LJ2,LMJ2
      INTEGER               ::  LJ12,LJJ_MIN,LJJ_MAX
      INTEGER               ::  LJJ,L_EXP,I_INT
!
      INTEGER               ::  INT
!
      REAL (WP)             ::  J1,J2,MJ1,MJ2,MJ3,JJ
      REAL (WP)             ::  JJ_MIN,JJ_MAX,JJ12
      REAL (WP)             ::  NJ(0:N_GAUNT)
!
      REAL (WP)             ::  FLOAT,SIGN,SQRT
!
!  Storage of the Clebsch-Gordan coefficients for the Auger multiplet
!                    dependent coupling matrix elements
!                   in the array CGA(LJ1,MJ1,LJ2,MJ2,LJ).
!
!  Here :       LJ1  is an integer index related to J1 (LJ1=2*J1)
!               LMJ1  is an integer index related to MJ1 (LMJ1=2*MJ1)
!               LJ2  is an integer index related to J2 (LJ2=2*J2)
!               LMJ2  is an integer index related to MJ2 (LMJ2=2*MJ2)
!               LJ is an integer index related to J :
!                     J = FLOAT(LJ) for J integer
!                     J = FLOAT(LJ) + 0.5 for J half integer
!
!              so that all indices remain integer
!
      IF((SPECTRO == 'AED') .OR. (SPECTRO == 'APC')) THEN           !
!
        N       = 3                                                 !
        MJ3    = ZERO                                               !
        LJ_MAX = 2 * (LI_I + LI_A + 1)                              !
        DO LJ1 = 0, LJ_MAX                                          !
          J1 = FLOAT(LJ1) * HALF                                    !
          DO LMJ1 = -LJ1, LJ1, 2                                    !
            MJ1 = FLOAT(LMJ1) * HALF                                !
            DO LJ2 = 0, LJ_MAX                                      !
              J2 = FLOAT(LJ2) * HALF                                !
              DO LMJ2 = -LJ2, LJ2, 2                                !
                MJ2 = FLOAT(LMJ2) * HALF                            !
                CALL N_J(J1,MJ1,J2,MJ2,MJ3,NJ,I_INT,N)              !
!
                JJ12 = J1 - J2                                      !
!
                LJ12 = INT(JJ12 + SIGN(SMALL,JJ12))                 !
!
                JJ_MIN  = ABS(LJ12)                                 !
                JJ_MAX  = J1 + J2                                   !
                LJJ_MIN = INT(JJ_MIN + SIGN(SMALL,JJ_MIN))          !
                LJJ_MAX = INT(JJ_MAX + SIGN(SMALL,JJ_MAX))          !
!
                DO LJJ = LJJ_MIN, LJJ_MAX, 1                        !
                  IF(I_INT == 1) THEN                               !
                    JJ = FLOAT(LJJ)                                 !
                  ELSE                                              !
                    JJ = FLOAT(LJJ) + HALF                          !
                  END IF                                            !
                  L_EXP = INT(J1 - J2 + MJ1 + MJ2)                  !
                  IF(MOD(L_EXP,2) == 0) THEN                        !
                    CGA(LJ1,LMJ1,LJ2,LMJ2,LJJ) = NJ(LJJ)   *      & !
                                                 SQRT(TWO * JJ + ONE)
                  ELSE
                    CGA(LJ1,LMJ1,LJ2,LMJ2,LJJ) = - NJ(LJJ) *      & !
                                                 SQRT(TWO * JJ + ONE)
                  END IF                                            !
                END DO                                              !
              END DO                                                !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
!
      END IF                                                        !
!
!  Storage of another of the spin Clebsch-Gordan used
!    when the Auger line is multiplet-resolved. It
!    originates from the coupling of SA and SC,
!    the spins of the Auger electron of the original
!    core electron (which is supposed to be the same
!    as that of the photoelectron).
!
!    CG_S(I,J,K) with : I = 1 ---> MSA = -1/2
!                       I = 2 ---> MSA =  1/2
!                       J = 1 ---> MSC = -1/2
!                       J = 2 ---> MSC =  1/2
!                       K = 1 ---> S   =  0
!                       K = 2 ---> S   =  1
!
!                       MS = MSA+MSC
!
      CG_S(1,1,1) =   ZERO                                          !
      CG_S(1,1,2) =   ONE                                           !
      CG_S(1,2,1) = - 0.70710678118654752440084436E0_WP             !
      CG_S(1,2,2) =   0.70710678118654752440084436E0_WP             !
      CG_S(2,1,1) =   0.70710678118654752440084436E0_WP             !
      CG_S(2,1,2) =   0.70710678118654752440084436E0_WP             !
      CG_S(2,2,1) =   ZERO                                          !
      CG_S(2,2,2) =   ONE                                           !
!
!  Initialization of the variables used when only one multiplet
!    is taken into account in the Auger peak
!
      IF((SPECTRO == 'AED').OR.(SPECTRO == 'APC')) THEN             !
!
        MULTIPLET=CHAR(48+IM1)//MULT//CHAR(48+IM2)                  !
!
        IF(MOD(IM1,2) == 0) THEN                                    !
          WRITE(IUO1,100) IM1                                       !
          STOP                                                      !
        END IF                                                      !
!
        S_MUL = (IM1 - 1) / 2                                       !
        J_MUL = IM2                                                 !
!
        IF(MULT == 'S') THEN                                        !
          L_MUL = 0                                                 !
        ELSE IF(MULT == 'P') THEN                                   !
          L_MUL = 1                                                 !
        ELSE IF(MULT == 'D') THEN                                   !
          L_MUL = 2                                                 !
        ELSE IF(MULT == 'F') THEN                                   !
          L_MUL = 3                                                 !
        ELSE IF(MULT == 'G') THEN                                   !
          L_MUL = 4                                                 !
        ELSE IF(MULT == 'H') THEN                                   !
          L_MUL = 5                                                 !
        ELSE IF(MULT == 'I') THEN                                   !
          L_MUL = 6                                                 !
        ELSE IF(MULT == 'K') THEN                                   !
          L_MUL = 7                                                 !
        ELSE IF(MULT == 'L') THEN                                   !
          L_MUL = 8                                                 !
        ELSE IF(MULT == 'M') THEN                                   !
          L_MUL = 9                                                 !
        ELSE                                                        !
          WRITE(IUO1,101) MULTIPLET                                 !
          STOP                                                      !
        END IF                                                      !
!
      END IF                                                        !
!
!  Formats
!
 100  FORMAT(///,4X,' <<<<<<<<<<  WRONG NAME FOR THE MULTIPLET',       &
             '  >>>>>>>>>>',/,4X,' <<<<<<<<<<  ODD NUMBER ',           &
             'EXPECTED INSTEAD OF',I2,'  >>>>>>>>>>')
 101  FORMAT(///,4X,' <<<<<<<<<<  ',A3,' IS NOT IMPLEMENTED IN THIS ', &
             'VERSION  >>>>>>>>>>')
!
      END SUBROUTINE INIT_MULTIPLET
!
!=======================================================================
!
      SUBROUTINE INIT_PARAM_CONS(SPECTRO,ITRTL)
!
!  This subroutine initializes various parameters for consistency's sake
!
!
!  Input variables :
!
!                       SPECTRO   :  type of spectroscopy:
!                                       ---> 'PED' : photoelectron diffraction
!                                            'LED' : low-energy electron diffraction
!                                            'XAS' : X-ray absorption spectroscopy
!                                            'AED' : Auger electron diffraction
!                                            'APC' : Auger photoelectron coincidence spectroscopy
!                                            'E2E' : (e,2e) spectroscopy
!                                            'E3E' : (e,3e) spectroscopy
!                                            'ELS' : electron energy loss spectroscopy
!                                            'RES' : resonant elastic X-ray scattering
!                                            'STM' : scanning tunneling microscopy
!                                            'PLS' : photoelectron energy loss spectroscopy
!                                            'BEM' : ballistic energy electron microscopy
!                                            'EIG' : eigenvalue calculation
!
!  Output variables :
!
!                       ITRTL     :  switch for the number of non zero digits above which
!                                    T-matrix elements will be considered to be null. by
!                                    default, these elements are read with nine decimal digits
!                                        ---> 1    : only tl values >= +/- 0.1 will be kept
!                                             2    : only tl values >= +/- 0.01 will be kept
!                                        ...
!
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 19 May 2021
!
!
      USE APPROXIMATIONS_IN
      USE APPROXIMATIONS_EX
      USE APPROXIMATIONS_O1
      USE APPROXIMATIONS_O2
      USE APPROXIMATIONS_O3
!
      USE CLUS_ELEC
      USE CURRENT_SPIN
      USE DAMPING
!
      USE SPIN_IN
      USE SPIN_EX
      USE SPIN_O1
      USE SPIN_O2
      USE SPIN_O3
!
      USE TMP_04
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  SPECTRO
!
      INTEGER, INTENT(OUT)  ::  ITRTL
!
      IF(ISPIN == 0) THEN                                           !
        NSPIN2 = 1                                                  !
        NSPIN  = 1                                                  !
      ELSE IF(ISPIN == 1) THEN                                      !
        NSPIN2 = 4                                                  !
        NSPIN  = 2                                                  !
      END IF                                                        !
!
      IF(SPECTRO == 'STM') THEN                                     !
        IATT_IN = IATT_TI                                           !
        IATT_EX = IATT_TI                                           !
        IATT_O1 = IATT                                              !
        IATT_O2 = IATT                                              !
        IATT_O3 = IATT                                              !
      ELSE                                                          !
        IATT_IN = IATT                                              !
        IATT_EX = IATT                                              !
        IATT_O1 = IATT                                              !
        IATT_O2 = IATT                                              !
        IATT_O3 = IATT                                              !
      END IF                                                        !
!
      ISPIN_IN  = ISPIN                                             !
      NSPIN_IN  = NSPIN                                             !
      NSPIN2_IN = NSPIN2                                            !
      ISPIN_EX  = ISPIN                                             !
      NSPIN_EX  = NSPIN                                             !
      NSPIN2_EX = NSPIN2                                            !
      ISPIN_O1  = ISPIN                                             !
      NSPIN_O1  = NSPIN                                             !
      NSPIN2_O1 = NSPIN2                                            !
      ISPIN_O2  = ISPIN                                             !
      NSPIN_O2  = NSPIN                                             !
      NSPIN2_O2 = NSPIN2                                            !
      ISPIN_O3  = ISPIN                                             !
      NSPIN_O3  = NSPIN                                             !
      NSPIN2_O3 = NSPIN2                                            !
!
      IF(NEL == 1) THEN                                             !
        IF(INC == 1) THEN                                           !
          ITRTL = ITRTL_IN                                          !
        ELSE IF(EXC == 1) THEN                                      !
          ITRTL = ITRTL_EX                                          !
        ELSE IF(OUT1 == 1) THEN                                     !
          ITRTL = ITRTL_O1                                          !
        END IF                                                      !
      ELSE IF(NEL == 2) THEN                                        !
        ITRTL = MAX(ITRTL_O1,ITRTL_O2)                              !
      ELSE IF(NEL == 3) THEN                                        !
        ITRTL = MAX(ITRTL_O1,ITRTL_O2,ITRTL_EX)                     !
      END IF                                                        !
!
      END SUBROUTINE INIT_PARAM_CONS
!
!=======================================================================
!
      SUBROUTINE SET_INT_CALC(NAT)
!
!  This routine sets up of the variables used for an internal calculation
!    of the inelastic mean free path and/or of the mean square displacements
!
!  Input variables :
!
!                       NAT       :  number of prototypical atoms
!
!
!  Output variables :   stored in module CURRENT_MEAN_FREE_PATH
!
!
!
!   Author :  D. Sébilleau
!
!                                            Last modified : 18 May 2021
!
!
      USE ATOMIC_MASS
      USE ATOMS
      USE CURRENT_MEAN_FREE_PATH
      USE DEB_WAL_CLU
      USE XM_RHO
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  NAT
!
      INTEGER               ::  JTYP
!
      REAL (WP)             ::  RHOT(NATM)
!
      IF((IMSD == 1) .OR. (IMFP == 1)) THEN                         !
        DO JTYP = 1, NAT                                            ! \
          XMT(JTYP)  = XM_AT(NZAT(JTYP))                            !  \ coming from
          RHOT(JTYP) = RHO_AT(NZAT(JTYP))                           !  /  XM_RHO
        END DO                                                      ! /
        XMT_A  = XMT(1)                                             !
        RHOT_A = RHOT(1)                                            ! \  stored in
        NZ_A   = NZAT(1)                                            !  > module
      END IF                                                        ! /  CURRENT_MEAN_FREE_PATH
!
      END SUBROUTINE SET_INT_CALC
!
!=======================================================================
!
      SUBROUTINE COEFPQ(NAT,NGR)
!
!  This subroutine computes the P(n,m) and Q(n) coefficients
!   involved in the correlation expansion formulation
!
!  Reference  : equations (2.15) and (2.16) of
!               H. Zhao, D. Sebilleau and Z. Wu,
!               J. Phys.: Condens. Matter 20, 275241 (2008)
!
!
!  Input variables :
!
!          NAT       :  number of atoms
!          NGR       :  correlation order
!
!
!   Author :  H.-F. Zhao 2007
!
!
!                                        Last modified (DS): 10 Jun 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE
      USE Q_ARRAY
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  NAT,NGR
!
      INTEGER               ::  LOGF
      INTEGER               ::  N,M,I
!
      REAL (WP)             ::  CMN(NGR_M,NGR_M)
      REAL (WP)             ::  P(NGR_M,NGR_M)
!
      LOGF = 6                                                      ! log file
!
      IF(NGR > NAT) THEN                                            !
        WRITE(LOGF,10)                                              !
        STOP                                                        !
      END IF                                                        !
!
      CALL CMNGR(NAT,NGR,CMN)                                       !
!
      DO N = 1, NGR                                                 !
        P(N,N) = ONE                                                !
        Q(N)   = P(N,N)                                             !
        DO M = N+1, NGR                                             !
          P(N,M) = ZERO                                             !
          DO I = N, M-1                                             !
            P(N,M) = P(N,M) - P(N,I) * CMN(I,M)                     !
          END DO                                                    !
          Q(N) = Q(N) + P(N,M)                                      !
!
        END DO                                                      !
!
      END DO                                                        !
        stop
!
!  Format:
!
  10  FORMAT(//,5X,'<<<<<  NGR is larger than NAT  >>>>>',//)
!
      END SUBROUTINE COEFPQ
!
!=======================================================================
!
      SUBROUTINE CMNGR(NAT,NGR,CMN)
!
!    This subroutine calculates C(NAT-N,M-N) where,
!
!                    1 <= M <= NGR <= NAT
!                    1 <= N <= M
!
!               C(NAT-N,M-N) is stored as CMN(N,M)
!
!
!  Input variables :
!
!          NAT       :  number of atoms
!          NGR       :  correlation order
!
!
!  Output variables :
!
!          CMN       :
!
!
!   Author :  H.-F. Zhao 2007
!
!
!                                        Last modified (DS): 10 Jun 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)    ::  NAT,NGR
!
      INTEGER                ::  LOGF
      INTEGER                ::  N,M
!
      REAL (WP), INTENT(OUT) ::  CMN(NGR_M,NGR_M)
!
      REAL (WP)              ::  FLOAT
!
      LOGF = 6                                                      ! log file
!
      IF(NGR > NAT) THEN                                            !
        WRITE(LOGF,10)                                              !
        STOP                                                        !
      END IF                                                        !
!
      DO M = 1,NGR                                                  !
        DO N = 1,NGR                                                !
          CMN(N,M) = ZERO                                           !
        END DO                                                      !
        CMN(M,M) = ONE                                              !
      END DO                                                        !
!
      DO M = 1, NGR                                                 !
        DO N = M-1, 1, -1                                           !
          CMN(N,M) = CMN(N+1,M) * FLOAT(NAT-N) / FLOAT(M-N)         !
        END DO                                                      !
      END DO                                                        !
!
!  Format:
!
  10  FORMAT(//,5X,'<<<<<  NGR is larger than NAT  >>>>>',//)
!
      END SUBROUTINE CMNGR
!
END MODULE INITIALIZE_CALC
