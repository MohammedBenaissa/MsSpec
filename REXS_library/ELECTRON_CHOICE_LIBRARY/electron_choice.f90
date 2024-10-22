!
!=======================================================================
!
MODULE ELECTRON_CHOICE
!
!  This module select the electron considered
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE CHOOSE_ELEC(JEL,NAT)
!
!  This routine selects the modules and variable corresponding
!   to the input electron
!
!  Input variables :
!
!                       JEL       :  electron number --> 1 : IN  incoming
!                                                        2 : EX  excited
!                                                        3 : O1  1st outgoing
!                                                        4 : O2  2nd outgoing
!                                                        5 : O3  3rd outgoing
!                       NAT       :  number of prototypical atoms
!
!  Output variables :
!
!                     all the corresponding variables are stored
!                     in the module T_MATRIX
!
!
!  Author : D. Sébilleau
!
!
!                                         Last modified : 27 May 2021
!
!
      USE CURRENT_T_MATRIX
!
      USE ENERGY_IN
      USE ENERGY_EX
      USE ENERGY_O1
      USE ENERGY_O2
      USE ENERGY_O3
!
      USE TLM_IN
      USE TLM_EX
      USE TLM_O1
      USE TLM_O2
      USE TLM_O3
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  JEL,NAT
!
      INTEGER                   ::  NE,JE,JAT,L
!
      IF(JEL == 1) THEN                                             !
        IPOTC = IPOTC_IN                                            !
        ITL   = ITL_IN                                              !
        NE    = NE_IN                                               !
      ELSE IF(JEL == 2) THEN                                        !
        IPOTC = IPOTC_EX                                            !
        ITL   = ITL_EX                                              !
        NE    = NE_EX                                               !
      ELSE IF(JEL == 3) THEN                                        !
        IPOTC = IPOTC_O1                                            !
        ITL   = ITL_O1                                              !
        NE    = NE_O1                                               !
      ELSE IF(JEL == 4) THEN                                        !
        IPOTC = IPOTC_O2                                            !
        ITL   = ITL_O2                                              !
        NE    = NE_O2                                               !
      ELSE IF(JEL == 5) THEN                                        !
        IPOTC = IPOTC_O3                                            !
        ITL   = ITL_O3                                              !
        NE    = NE_O3                                               !
      END IF                                                        !

      DO JE = 1, NE                                                 !

        IF(JEL == 1) THEN                                           !
          VK(JE)  = KL_IN(JE)                                       !
          VK2(JE) = K2L_IN(JE)                                      !
        ELSE IF(JEL == 2) THEN                                      !
          VK(JE)  = KL_EX(JE)                                       !
          VK2(JE) = K2L_EX(JE)                                      !
        ELSE IF(JEL == 3) THEN                                      !
          VK(JE)  = KL_O1(JE)                                       !
          VK2(JE) = K2L_O1(JE)                                      !
        ELSE IF(JEL == 4) THEN                                      !
          VK(JE)  = KL_O2(JE)                                       !
          VK2(JE) = K2L_O2(JE)                                      !
        ELSE IF(JEL == 5) THEN                                      !
          VK(JE)  = KL_O3(JE)                                       !
          VK2(JE) = K2L_O3(JE)                                      !
        END IF                                                      !

        DO JAT = 1, NAT                                             !

          IF(JEL == 1) THEN                                         !
            LMAX(JAT,JE) = LMAX_IN(JAT,JE)                          !
            DO L = 0, LMAX(JAT,JE)                                  !
              TL(L,1,JAT,JE) = TL_IN(L,1,JAT,JE)                    !
            END DO                                                  !
          ELSE IF(JEL == 2) THEN                                    !
            LMAX(JAT,JE) = LMAX_EX(JAT,JE)                          !
            DO L = 0, LMAX(JAT,JE)                                  !
              TL(L,1,JAT,JE) = TL_EX(L,1,JAT,JE)                    !
            END DO                                                  !
          ELSE IF(JEL == 3) THEN                                    !
            LMAX(JAT,JE) = LMAX_O1(JAT,JE)                          !
            DO L = 0, LMAX(JAT,JE)                                  !
              TL(L,1,JAT,JE) = TL_O1(L,1,JAT,JE)                    !
            END DO                                                  !
          ELSE IF(JEL == 4) THEN                                    !
            LMAX(JAT,JE) = LMAX_O2(JAT,JE)                          !
            DO L = 0, LMAX(JAT,JE)                                  !
              TL(L,1,JAT,JE) = TL_O2(L,1,JAT,JE)                    !
            END DO                                                  !
          ELSE IF(JEL == 5) THEN                                    !
            LMAX(JAT,JE) = LMAX_O3(JAT,JE)                          !
            DO L = 0, LMAX(JAT,JE)                                  !
              TL(L,1,JAT,JE) = TL_O3(L,1,JAT,JE)                    !
            END DO                                                  !
          END IF                                                    !

        END DO                                                      !

      END DO                                                        !

      END SUBROUTINE CHOOSE_ELEC
!
!=======================================================================
!
      SUBROUTINE CHG_MOD(JEL,NAT)
!
!   This subroutine redistributes some module blocks depending
!     on the number of electrons involved in the spectroscopy.
!     It allows in particular to simplify the structure when
!     only one electron is involved in the process.
!
!  Input variables :
!
!                       JEL       :  electron number --> 1 : IN  incoming
!                                                        2 : EX  excited
!                                                        3 : O1  1st outgoing
!                                                        4 : 02  2nd outgoing
!                                                        5 : 03  3rd outgoing
!                       NAT       :  number of prototypical atoms
!
!
!
!  Author : D. Sébilleau
!
!
!                                         Last modified : 16 Jun 2021
!
!
!
      USE AMPLI_IN
      USE AMPLI_EX
      USE AMPLI_O1
      USE AMPLI_O2
      USE AMPLI_O3
!
      USE APPROXIMATIONS_IN
      USE APPROXIMATIONS_EX
      USE APPROXIMATIONS_O1
      USE APPROXIMATIONS_O2
      USE APPROXIMATIONS_O3
!
      USE AVER_IN
      USE AVER_O1
      USE AVER_O2
      USE AVER_O3
!
      USE ENERGY_IN
      USE ENERGY_EX
      USE ENERGY_O1
      USE ENERGY_O2
      USE ENERGY_O3
!
      USE EXTERNAL_IN
      USE EXTERNAL_O1
      USE EXTERNAL_O2
      USE EXTERNAL_O3
!
      USE EXT_DIR_IN
      USE EXT_DIR_O1
      USE EXT_DIR_O2
      USE EXT_DIR_O3
!
      USE F_TH_IN
      USE F_TH_EX
      USE F_TH_O1
      USE F_TH_O2
      USE F_TH_O3
!
      USE INTEGR_IN
      USE INTEGR_O1
      USE INTEGR_O2
      USE INTEGR_O3
!
      USE LINLBD_IN
      USE LINLBD_EX
      USE LINLBD_O1
      USE LINLBD_O2
      USE LINLBD_O3
!
      USE MFP_IN
      USE MFP_EX
      USE MFP_O1
      USE MFP_O2
      USE MFP_O3
!
      USE PHI_IN
      USE PHI_O1
      USE PHI_O2
      USE PHI_O3
!
      USE POLAR_IN
      USE POLAR_O1
      USE POLAR_O2
      USE POLAR_O3
!
      USE RA_IN
      USE RA_EX
      USE RA_O1
      USE RA_O2
      USE RA_O3
!
      USE RENORM_IN
      USE RENORM_EX
      USE RENORM_O1
      USE RENORM_O2
      USE RENORM_O3
!
      USE SPIN_IN
      USE SPIN_EX
      USE SPIN_O1
      USE SPIN_O2
      USE SPIN_O3
!
      USE TEST_IN
      USE TEST_EX
      USE TEST_O1
      USE TEST_O2
      USE TEST_O3
!
!
      USE THETA_IN
      USE THETA_O1
      USE THETA_O1
      USE THETA_O3
!
      USE CURRENT_AVER                    !  <--
      USE CURRENT_BEAM                    !  <--
      USE CURRENT_CALC                    !  <--
      USE CURRENT_EXT_DIR                 !  <--
      USE CURRENT_F_TH                    !  <--
      USE CURRENT_FINA_VAL                !  <--
      USE CURRENT_INIT_VAL                !  <--
      USE CURRENT_LINLBD                  !  <--
      USE CURRENT_MEAN_FREE_PATH          !  <--
      USE CURRENT_REHR_ALBERS             !  <--
      USE CURRENT_RENORM                  !  <--
      USE CURRENT_SPIN                    !  <--
!
      IF(JEL == 1) THEN                                             !
!
        I_AMP   = I_AMP_IN                                          !
        IE      = IE_IN                                             !
        NE      = NE_IN                                             !
        I_INT   = I_INT_IN                                          !
        IPHI    = IPHI_IN                                           !
        NPHI    = NPHI_IN                                           !
        IPOL    = IPOL_IN                                           !
        ITHETA  = ITHETA_IN                                         !
        NTHETA  = NTHETA_IN                                         !
        I_TEST  = I_TEST_IN                                         !
        I_EXT   = I_EXT_IN                                          !
        IDIR    = IDIR_IN                                           !
        NSET    = NSET_IN                                           !
        N_POINTS= N_POINTS_IN                                       !
        N_SCAT    = N_SCAT_IN                                       !
        NO      = NO_IN                                             !
        I_BASIS = I_BASIS_IN                                        !
        IFWD    = IFWD_IN                                           !
        NTHOUT  = NTHOUT_IN                                         !
        IPW     = IPW_IN                                            !
        NCUT    = NCUT_IN                                           !
        PCTINT  = PCTINT_IN                                         !
        IPP     = IPP_IN                                            !
        IATTS   = IATT_IN                                           !
        ILENGTH = ILENGTH_IN                                        !
        RLENGTH = RLENGTH_IN                                        !
        I_NO    = I_NO_IN                                           !
        I_RA    = I_RA_IN                                           !
        ISPIN   = ISPIN_IN                                          !
        IDICHR  = IDICHR_IN                                         !
        NSPIN   = NSPIN_IN                                          !
        NSPIN2  = NSPIN2_IN                                         !
        ISFLIP  = ISFLIP_IN                                         !
        IR_DIA  = IR_DIA_IN                                         !
        NSTEP   = NSTEP_IN                                          !
        IMFP    = IMFP_IN                                           !
        XMFP0   = XMFP0_IN                                          !
        IFTHET  = IFTHET_IN                                         !
        NFTHET  = NFTHET_IN                                         !
        R0      = R0_IN                                             !
        R1      = R1_IN                                             !
        E0      = E0_IN                                             !
        E1      = E1_IN                                             !
        PHI0    = PHI0_IN                                           !
        PHI1    = PHI1_IN                                           !
        THETA0  = THETA0_IN                                         !
        THETA1  = THETA1_IN                                         !
        DO JAT = 1 ,NAT                                             !
          RTHFWD(JAT) = RTHFWD_IN(JAT)                              !
          RTHBWD(JAT) = RTHBWD_IN(JAT)                              !
          IBWD(JAT)   = IBWD_IN(JAT)                                !
          NUMAX(JAT)  = NUMAX_IN(JAT)                               !
          N_RA(JAT)   = N_RA_IN(JAT)                                !
        END DO                                                      !
        I_AVER  = I_AVER_IN                                         !
        NDIR    = NDIR_IN                                           !
        ICHKDIR = ICHKDIR_IN                                        !
        ACCEPT  = ACCEPT_IN                                         !
        I_REN   = I_REN_IN                                          !
        N_REN   = N_REN_IN                                          !
        REN_R   = REN_R_IN                                          !
        REN_I   = REN_I_IN                                          !
!
      ELSE IF(JEL == 2) THEN                                        !
!
        I_AMP   = I_AMP_EX                                          !
        IE      = IE_EX                                             !
        NE      = NE_EX                                             !
        I_INT   = I_INT_EX                                          !
        I_TEST  = I_TEST_EX                                         !
        N_SCAT    = N_SCAT_EX                                       !
        NO      = NO_EX                                             !
        I_BASIS = I_BASIS_EX                                        !
        IFWD    = IFWD_EX                                           !
        NTHOUT  = NTHOUT_EX                                         !
        IPW     = IPW_EX                                            !
        NCUT    = NCUT_EX                                           !
        PCTINT  = PCTINT_EX                                         !
        IPP     = IPP_EX                                            !
        IATTS   = IATT_EX                                           !
        ILENGTH = ILENGTH_EX                                        !
        RLENGTH = RLENGTH_EX                                        !
        I_NO    = I_NO_EX                                           !
        I_RA    = I_RA_EX                                           !
        ISPIN   = ISPIN_EX                                          !
        IDICHR  = IDICHR_EX                                         !
        NSPIN   = NSPIN_EX                                          !
        NSPIN2  = NSPIN2_EX                                         !
        ISFLIP  = ISFLIP_EX                                         !
        IR_DIA  = IR_DIA_EX                                         !
        NSTEP   = NSTEP_EX                                          !
        IMFP    = IMFP_EX                                           !
        XMFP0   = XMFP0_EX                                          !
        IFTHET  = IFTHET_EX                                         !
        NFTHET  = NFTHET_EX                                         !
        R0      = R0_EX                                             !
        R1      = R1_EX                                             !
        E0      = E0_EX                                             !
        E1      = E1_EX                                             !
        DO JAT = 1, NAT                                             !
          RTHFWD(JAT) = RTHFWD_EX(JAT)                              !
          RTHBWD(JAT) = RTHBWD_EX(JAT)                              !
          IBWD(JAT)   = IBWD_EX(JAT)                                !
          NUMAX(JAT)  = NUMAX_EX(JAT)                               !
          N_RA(JAT)   = N_RA_EX(JAT)                                !
        END DO                                                      !
        I_REN   = I_REN_EX                                          !
        N_REN   = N_REN_EX                                          !
        REN_R   = REN_R_EX                                          !
        REN_I   = REN_I_EX                                          !
!
      ELSE IF(JEL == 3) THEN                                        !
!
        I_AMP   = I_AMP_O1                                          !
        IE      = IE_O1                                             !
        NE      = NE_O1                                             !
        I_INT   = I_INT_O1                                          !
        IPHI    = IPHI_O1                                           !
        NPHI    = NPHI_O1                                           !
        IPOL    = IPOL_O1                                           !
        ITHETA  = ITHETA_O1                                         !
        NTHETA  = NTHETA_O1                                         !
        I_TEST  = I_TEST_O1                                         !
        I_EXT   = I_EXT_O1                                          !
        IDIR    = IDIR_O1                                           !
        NSET    = NSET_O1                                           !
        N_POINTS= N_POINTS_O1                                       !
        N_SCAT    = N_SCAT_O1                                       !
        NO      = NO_O1                                             !
        I_BASIS = I_BASIS_O1                                        !
        IFWD    = IFWD_O1                                           !
        NTHOUT  = NTHOUT_O1                                         !
        IPW     = IPW_O1                                            !
        NCUT    = NCUT_O1                                           !
        PCTINT  = PCTINT_O1                                         !
        IPP     = IPP_O1                                            !
        IATTS   = IATT_O1                                           !
        ILENGTH = ILENGTH_O1                                        !
        RLENGTH = RLENGTH_O1                                        !
        I_NO    = I_NO_O1                                           !
        I_RA    = I_RA_O1                                           !
        ISPIN   = ISPIN_O1                                          !
        IDICHR  = IDICHR_O1                                         !
        NSPIN   = NSPIN_O1                                          !
        NSPIN2  = NSPIN2_O1                                         !
        ISFLIP  = ISFLIP_O1                                         !
        IR_DIA  = IR_DIA_O1                                         !
        NSTEP   = NSTEP_O1                                          !
        IMFP    = IMFP_O1                                           !
        XMFP0   = XMFP0_O1                                          !
        IFTHET  = IFTHET_O1                                         !
        NFTHET  = NFTHET_O1                                         !
        R0      = R0_O1                                             !
        R1      = R1_O1                                             !
        E0      = E0_O1                                             !
        E1      = E1_O1                                             !
        PHI0    = PHI0_O1                                           !
        PHI1    = PHI1_O1                                           !
        THETA0  = THETA0_O1                                         !
        THETA1  = THETA1_O1                                         !
        DO JAT = 1, NAT                                             !
          RTHFWD(JAT) = RTHFWD_O1(JAT)                              !
          RTHBWD(JAT) = RTHBWD_O1(JAT)                              !
          IBWD(JAT)   = IBWD_O1(JAT)                                !
          NUMAX(JAT)  = NUMAX_O1(JAT)                               !
          N_RA(JAT)   = N_RA_O1(JAT)                                !
        END DO                                                      !
        I_AVER  = I_AVER_O1                                         !
        NDIR    = NDIR_O1                                           !
        ICHKDIR = ICHKDIR_O1                                        !
        ACCEPT  = ACCEPT_O1                                         !
        I_REN   = I_REN_O1                                          !
        N_REN   = N_REN_O1                                          !
        REN_R   = REN_R_O1                                          !
        REN_I   = REN_I_O1                                          !
!
      ELSE IF(JEL == 4) THEN                                        !
!
        I_AMP   = I_AMP_O2                                          !
        IE      = IE_O2                                             !
        NE      = NE_O2                                             !
        I_INT   = I_INT_O2                                          !
        IPHI    = IPHI_O2                                           !
        NPHI    = NPHI_O2                                           !
        IPOL    = IPOL_O2                                           !
        ITHETA  = ITHETA_O2                                         !
        NTHETA  = NTHETA_O2                                         !
        I_TEST  = I_TEST_O2                                         !
        I_EXT   = I_EXT_O2                                          !
        IDIR    = IDIR_O2                                           !
        NSET    = NSET_O2                                           !
        N_POINTS= N_POINTS_O2                                       !
        N_SCAT    = N_SCAT_O2                                       !
        NO      = NO_O2                                             !
        I_BASIS = I_BASIS_O2                                        !
        IFWD    = IFWD_O2                                           !
        NTHOUT  = NTHOUT_O2                                         !
        IPW     = IPW_O2                                            !
        NCUT    = NCUT_O2                                           !
        PCTINT  = PCTINT_O2                                         !
        IPP     = IPP_O2                                            !
        IATTS   = IATT_O2                                           !
        ILENGTH = ILENGTH_O2                                        !
        RLENGTH = RLENGTH_O2                                        !
        I_NO    = I_NO_O2                                           !
        I_RA    = I_RA_O2                                           !
        ISPIN   = ISPIN_O2                                          !
        IDICHR  = IDICHR_O2                                         !
        NSPIN   = NSPIN_O2                                          !
        NSPIN2  = NSPIN2_O2                                         !
        ISFLIP  = ISFLIP_O2                                         !
        IR_DIA  = IR_DIA_O2                                         !
        NSTEP   = NSTEP_O2                                          !
        IMFP    = IMFP_O2                                           !
        XMFP0   = XMFP0_O2                                          !
        IFTHET  = IFTHET_O2                                         !
        NFTHET  = NFTHET_O2                                         !
        R0      = R0_O2                                             !
        R1      = R1_O2                                             !
        E0      = E0_O2                                             !
        E1      = E1_O2                                             !
        PHI0    = PHI0_O2                                           !
        PHI1    = PHI1_O2                                           !
        THETA0  = THETA0_O2                                         !
        THETA1  = THETA1_O2                                         !
        DO JAT = 1, NAT                                             !
          RTHFWD(JAT) = RTHFWD_O2(JAT)                              !
          RTHBWD(JAT) = RTHBWD_O2(JAT)                              !
          IBWD(JAT)   = IBWD_O2(JAT)                                !
          NUMAX(JAT)  = NUMAX_O2(JAT)                               !
          N_RA(JAT)   = N_RA_O2(JAT)                                !
        END DO                                                      !
        I_AVER  = I_AVER_O2                                         !
        NDIR    = NDIR_O2                                           !
        ICHKDIR = ICHKDIR_O2                                        !
        ACCEPT  = ACCEPT_O2                                         !
        I_REN   = I_REN_O2                                          !
        N_REN   = N_REN_O2                                          !
        REN_R   = REN_R_O2                                          !
        REN_I   = REN_I_O2                                          !
!
      ELSE IF(JEL == 5) THEN                                        !
!
        I_AMP   = I_AMP_O3                                          !
        IE      = IE_O3                                             !
        NE      = NE_O3                                             !
        I_INT   = I_INT_O3                                          !
        IPHI    = IPHI_O3                                           !
        NPHI    = NPHI_O3                                           !
        IPOL    = IPOL_O3                                           !
        ITHETA  = ITHETA_O3                                         !
        NTHETA  = NTHETA_O3                                         !
        I_TEST  = I_TEST_O3                                         !
        I_EXT   = I_EXT_O3                                          !
        IDIR    = IDIR_O3                                           !
        NSET    = NSET_O3                                           !
        N_POINTS= N_POINTS_O3                                       !
        N_SCAT    = N_SCAT_O3                                       !
        NO      = NO_O3                                             !
        I_BASIS = I_BASIS_O3                                        !
        IFWD    = IFWD_O3                                           !
        NTHOUT  = NTHOUT_O3                                         !
        IPW     = IPW_O3                                            !
        NCUT    = NCUT_O3                                           !
        PCTINT  = PCTINT_O3                                         !
        IPP     = IPP_O3                                            !
        IATTS   = IATT_O3                                           !
        ILENGTH = ILENGTH_O3                                        !
        RLENGTH = RLENGTH_O3                                        !
        I_NO    = I_NO_O3                                           !
        I_RA    = I_RA_O3                                           !
        ISPIN   = ISPIN_O3                                          !
        IDICHR  = IDICHR_O3                                         !
        NSPIN   = NSPIN_O3                                          !
        NSPIN2  = NSPIN2_O3                                         !
        ISFLIP  = ISFLIP_O3                                         !
        IR_DIA  = IR_DIA_O3                                         !
        NSTEP   = NSTEP_O3                                          !
        IMFP    = IMFP_O3                                           !
        XMFP0   = XMFP0_O3                                          !
        IFTHET  = IFTHET_O3                                         !
        NFTHET  = NFTHET_O3                                         !
        R0      = R0_O3                                             !
        R1      = R1_O3                                             !
        E0      = E0_O3                                             !
        E1      = E1_O3                                             !
        PHI0    = PHI0_O3                                           !
        PHI1    = PHI1_O3                                           !
        THETA0  = THETA0_O3                                         !
        THETA1  = THETA1_O3                                         !
        DO JAT = 1, NAT                                             !
          RTHFWD(JAT) = RTHFWD_O3(JAT)                              !
          RTHBWD(JAT) = RTHBWD_O3(JAT)                              !
          IBWD(JAT)   = IBWD_O3(JAT)                                !
          NUMAX(JAT)  = NUMAX_O3(JAT)                               !
          N_RA(JAT)    =N_RA_O3(JAT)                                !
        END DO                                                      !
        I_AVER  = I_AVER_O3                                         !
        NDIR    = NDIR_O3                                           !
        ICHKDIR = ICHKDIR_O3                                        !
        ACCEPT  = ACCEPT_O3                                         !
        I_REN   = I_REN_O3                                          !
        N_REN   = N_REN_O3                                          !
        REN_R   = REN_R_O3                                          !
        REN_I   = REN_I_O3                                          !
!
      END IF                                                        !
!
      NEPS = 2 - ABS(IPOL_IN)                                       !
!
      END SUBROUTINE CHG_MOD
!
END MODULE ELECTRON_CHOICE
