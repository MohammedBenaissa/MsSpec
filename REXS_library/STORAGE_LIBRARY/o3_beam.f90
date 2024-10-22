!
!  These modules contain the parameters controlling the 3rd outgoing beam
!
!=======================================================================
!
MODULE AMPLI_O3
!
!  This module contains the switch controlling the storage of the amplitudes
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AMP_O3
!
END MODULE AMPLI_O3
!
!=======================================================================
!
MODULE APPROXIMATIONS_O3
!
!  This module contains the variables for the approximations
!    of the thirs outgoing beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  N_SCAT_O3,NO_O3,I_BASIS_O3
      INTEGER               ::  IFWD_O3,NTHOUT_O3
      INTEGER               ::  IBWD_O3(NATP_M),IPW_O3,NCUT_O3
      INTEGER               ::  IPP_O3,IATT_O3,ILENGTH_O3
!
      REAL (WP)             ::  RTHFWD_O3(NATP_M),RTHBWD_O3(NATP_M)
      REAL (WP)             ::  PCTINT_O3,RLENGTH_O3
!
END MODULE APPROXIMATIONS_O3
!
!=======================================================================
!
MODULE AVER_O3
!
!  This module contains general parameters for the averaging over
!    the 3rd outgoing beam direction
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AVER_O3,NDIR_O3,ICHKDIR_O3
!
      REAL (WP)             ::  ACCEPT_O3
!
END MODULE AVER_O3
!
!=======================================================================
!
MODULE ENERGY_O3
!
!  This module contains the parameters controlling the energy
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IE_O3,NE_O3
!
      REAL (WP)             ::  E0_O3,E1_O3
!
END MODULE ENERGY_O3
!
!=======================================================================
!
MODULE EXTERNAL_O3
!
!  This module contains the switch controlling the external angle values
!
      IMPLICIT NONE
!
      INTEGER               ::  I_EXT_O3
!
END MODULE EXTERNAL_O3
!
!=======================================================================
!
MODULE EXT_DIR_O3
!
!  This module contains parameters concerning the external directions
!
      IMPLICIT NONE
!
      INTEGER               ::  IDIR_O3,NSET_O3,N_POINTS_O3
!
END MODULE EXT_DIR_O3
!
!=======================================================================
!
MODULE F_TH_O3
!
!  This module contains the scattering amplitude parameters
!    for the 3rd outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IFTHET_O3,NFTHET_O3
!
      REAL (WP)             ::  R0_O3,R1_O3
!
END MODULE F_TH_O3
!
!=======================================================================
!
MODULE INTEGR_O3
!
!  This module contains the switch controlling the integration
!
      IMPLICIT NONE
!
      INTEGER               ::  I_INT_O3
!
END MODULE INTEGR_O3
!
!=======================================================================
!
MODULE LINLBD_O3
!
!  This module contains the Rehr-Albers variables for 3rd outgoing beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M,N_MU_M,N_NU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  LBD_O3(-N_MU_M:N_MU_M,0:N_NU_M)
      INTEGER               ::  LBDMAX_O3
      INTEGER               ::  NUMAX_O3(NATP_M)
!
END MODULE LINLBD_O3
!
!=======================================================================
!
MODULE MFP_O3
!
!  This module contains mfp value for the 3rd outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IMFP_O3
!
      REAL (WP)             ::  XMFP0_O3
!
END MODULE MFP_O3
!
!=======================================================================
!
MODULE PHI_O3
!
!  This module contains the parameters controlling the azimuthal angle
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IPHI_O3,NPHI_O3
!
      REAL (WP)             ::  PHI0_O3,PHI1_O3
!
END MODULE PHI_O3
!
!=======================================================================
!
MODULE POLAR_O3
!
!  This module contains the parameters controlling the polarization
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOL_O3
!
END MODULE POLAR_O3
!
!=======================================================================
!
MODULE RA_O3
!
!  This module contains Rehr-Albers parameters for the 3rd outgoing beam
!
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  I_NO_O3,I_RA_O3,N_RA_O3(NATP_M)
!
END MODULE RA_O3
!
!=======================================================================
!
MODULE RENORM_O3
!
!  This module contains renormalization parameters for the 3rd outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_REN_O3,N_REN_O3
!
      REAL (WP)             ::  REN_R_O3,REN_I_O3
!
END MODULE RENORM_O3
!
!=======================================================================
!
MODULE SPIN_O3
!
!  This module contains the spin parameters for the 3rd outgoing beam
!
      IMPLICIT NONE
!
      INTEGER               ::  ISPIN_O3,NSPIN_O3,NSPIN2_O3
      INTEGER               ::  ISFLIP_O3,IR_DIA_O3,NSTEP_O3
!
END MODULE SPIN_O3
!
!=======================================================================
!
MODULE TEST_O3
!
!  This module contains the parameters controlling the tests
!
      IMPLICIT NONE
!
      INTEGER               ::  I_TEST_O3
!
END MODULE TEST_O3
!
!=======================================================================
!
MODULE THETA_O3
!
!  This module contains the parameters controlling the polar angle
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  ITHETA_O3,NTHETA_O3
!
      REAL (WP)             ::  THETA0_O3,THETA1_O3
!
END MODULE THETA_O3


