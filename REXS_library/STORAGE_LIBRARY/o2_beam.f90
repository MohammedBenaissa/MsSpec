!
!  These modules contain the parameters controlling the 2nd outgoing beam
!
!=======================================================================
!
MODULE AMPLI_O2
!
!  This module contains the switch controlling the storage of the amplitudes
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AMP_O2
!
END MODULE AMPLI_O2
!
!=======================================================================
!
MODULE APPROXIMATIONS_O2
!
!  This module contains the variables for the approximations
!    of the second outgoing beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  N_SCAT_O2,NO_O2,I_BASIS_O2
      INTEGER               ::  IFWD_O2,NTHOUT_O2
      INTEGER               ::  IBWD_O2(NATP_M),IPW_O2,NCUT_O2
      INTEGER               ::  IPP_O2,IATT_O2,ILENGTH_O2
!
      REAL (WP)             ::  RTHFWD_O2(NATP_M),RTHBWD_O2(NATP_M)
      REAL (WP)             ::  PCTINT_O2,RLENGTH_O2
!
END MODULE APPROXIMATIONS_O2
!
!=======================================================================
!
MODULE AVER_O2
!
!  This module contains general parameters for the averaging over
!    the 2nd outgoing beam direction
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AVER_O2,NDIR_O2,ICHKDIR_O2
!
      REAL (WP)             ::  ACCEPT_O2
!
END MODULE AVER_O2
!
!=======================================================================
!
MODULE ENERGY_O2
!
!  This module contains the parameters controlling the energy
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IE_O2,NE_O2
!
      REAL (WP)             ::  E0_O2,E1_O2
!
END MODULE ENERGY_O2
!
!=======================================================================
!
MODULE EXTERNAL_O2
!
!  This module contains the switch controlling the external angle values
!
      IMPLICIT NONE
!
      INTEGER               ::  I_EXT_O2
!
END MODULE EXTERNAL_O2
!
!=======================================================================
!
MODULE EXT_DIR_O2
!
!  This module contains parameters concerning the external directions
!
      IMPLICIT NONE
!
      INTEGER               ::  IDIR_O2,NSET_O2,N_POINTS_O2
!
END MODULE EXT_DIR_O2
!
!=======================================================================
!
MODULE F_TH_O2
!
!  This module contains the scattering amplitude parameters
!    for the 2nd outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IFTHET_O2,NFTHET_O2
!
      REAL (WP)             ::  R0_O2,R1_O2
!
END MODULE F_TH_O2
!
!=======================================================================
!
MODULE INTEGR_O2
!
!  This module contains the switch controlling the integration
!
      IMPLICIT NONE
!
      INTEGER               ::  I_INT_O2
!
END MODULE INTEGR_O2
!
!=======================================================================
!
MODULE LINLBD_O2
!
!  This module contains the Rehr-Albers variables for 2nd outgoing beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M,N_MU_M,N_NU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  LBD_O2(-N_MU_M:N_MU_M,0:N_NU_M)
      INTEGER               ::  LBDMAX_O2
      INTEGER               ::  NUMAX_O2(NATP_M)
!
END MODULE LINLBD_O2
!
!=======================================================================
!
MODULE MFP_O2
!
!  This module contains mfp value for the 2nd outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IMFP_O2
!
      REAL (WP)             ::  XMFP0_O2
!
END MODULE MFP_O2
!
!=======================================================================
!
MODULE PHI_O2
!
!  This module contains the parameters controlling the azimuthal angle
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IPHI_O2,NPHI_O2
!
      REAL (WP)             ::  PHI0_O2,PHI1_O2
!
END MODULE PHI_O2
!
!=======================================================================
!
MODULE POLAR_O2
!
!  This module contains the parameters controlling the polarization
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOL_O2
!
END MODULE POLAR_O2
!
!=======================================================================
!
MODULE RA_O2
!
!  This module contains Rehr-Albers parameters for the 2nd outgoing beam
!
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  I_NO_O2,I_RA_O2,N_RA_O2(NATP_M)
!
END MODULE RA_O2
!
!=======================================================================
!
MODULE RENORM_O2
!
!  This module contains renormalization parameters for the 2nd outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_REN_O2,N_REN_O2
!
      REAL (WP)             ::  REN_R_O2,REN_I_O2
!
END MODULE RENORM_O2
!
!=======================================================================
!
MODULE SPIN_O2
!
!  This module contains the spin parameters for the 2nd outgoing beam
!
      IMPLICIT NONE
!
      INTEGER               ::  ISPIN_O2,NSPIN_O2,NSPIN2_O2
      INTEGER               ::  ISFLIP_O2,IR_DIA_O2,NSTEP_O2
!
END MODULE SPIN_O2
!
!=======================================================================
!
MODULE TEST_O2
!
!  This module contains the parameters controlling the tests
!
      IMPLICIT NONE
!
      INTEGER               ::  I_TEST_O2
!
END MODULE TEST_O2
!
!=======================================================================
!
MODULE THETA_O2
!
!  This module contains the parameters controlling the polar angle
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  ITHETA_O2,NTHETA_O2
!
      REAL (WP)             ::  THETA0_O2,THETA1_O2
!
END MODULE THETA_O2


