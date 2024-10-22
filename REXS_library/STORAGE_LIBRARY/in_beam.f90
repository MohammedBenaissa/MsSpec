!
!  These modules contain the parameters controlling the incoming beam
!
!=======================================================================
!
MODULE AMPLI_IN
!
!  This module contains the switch controlling the storage of the amplitudes
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AMP_IN
!
END MODULE AMPLI_IN
!
!=======================================================================
!
MODULE APPROXIMATIONS_IN
!
!  This module contains the variables for the approximations
!    of the incoming beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  N_SCAT_IN,NO_IN,I_BASIS_IN
      INTEGER               ::  IFWD_IN,NTHOUT_IN
      INTEGER               ::  IBWD_IN(NATP_M),IPW_IN,NCUT_IN
      INTEGER               ::  IPP_IN,IATT_IN,ILENGTH_IN
!
      REAL (WP)             ::  RTHFWD_IN(NATP_M),RTHBWD_IN(NATP_M)
      REAL (WP)             ::  PCTINT_IN,RLENGTH_IN
!
END MODULE APPROXIMATIONS_IN
!
!=======================================================================
!
MODULE AVER_IN
!
!  This module contains general parameters for the averaging over
!    the incoming beam direction
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AVER_IN,NDIR_IN,ICHKDIR_IN
!
      REAL (WP)             ::  ACCEPT_IN
!
END MODULE AVER_IN
!
!=======================================================================
!
MODULE ENERGY_IN
!
!  This module contains the parameters controlling the energy
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IE_IN,NE_IN
!
      REAL (WP)             ::  E0_IN,E1_IN
!
END MODULE ENERGY_IN
!
!=======================================================================
!
MODULE EXTERNAL_IN
!
!  This module contains the switch controlling the external angle values
!
      IMPLICIT NONE
!
      INTEGER               ::  I_EXT_IN
!
END MODULE EXTERNAL_IN
!
!=======================================================================
!
MODULE EXT_DIR_IN
!
!  This module contains parameters concerning the external directions
!
      IMPLICIT NONE
!
      INTEGER               ::  IDIR_IN,NSET_IN,N_POINTS_IN
!
END MODULE EXT_DIR_IN
!
!=======================================================================
!
MODULE F_TH_IN
!
!  This module contains the scattering amplitude parameters
!    for the incoming beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IFTHET_IN,NFTHET_IN
!
      REAL (WP)             ::  R0_IN,R1_IN
!
END MODULE F_TH_IN
!
!=======================================================================
!
MODULE INTEGR_IN
!
!  This module contains the switch controlling the integration
!
      IMPLICIT NONE
!
      INTEGER               ::  I_INT_IN
!
END MODULE INTEGR_IN
!
!=======================================================================
!
MODULE LINLBD_IN
!
!  This module contains the Rehr-Albers variables for incoming beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M,N_MU_M,N_NU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  LBD_IN(-N_MU_M:N_MU_M,0:N_NU_M)
      INTEGER               ::  LBDMAX_IN
      INTEGER               ::  NUMAX_IN(NATP_M)
!
END MODULE LINLBD_IN
!
!=======================================================================
!
MODULE MFP_IN
!
!  This module contains mfp value for the incoming beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IMFP_IN
!
      REAL (WP)             ::  XMFP0_IN
!
END MODULE MFP_IN
!
!=======================================================================
!
MODULE PHI_IN
!
!  This module contains the parameters controlling the azimuthal angle
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IPHI_IN,NPHI_IN
!
      REAL (WP)             ::  PHI0_IN,PHI1_IN
!
END MODULE PHI_IN
!
!=======================================================================
!
MODULE POLAR_IN
!
!  This module contains the parameters controlling the polarization
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOL_IN
!
END MODULE POLAR_IN
!
!=======================================================================
!
MODULE RA_IN
!
!  This module contains Rehr-Albers parameters for the incoming beam
!
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  I_NO_IN,I_RA_IN,N_RA_IN(NATP_M)
!
END MODULE RA_IN
!
!=======================================================================
!
MODULE RENORM_IN
!
!  This module contains renormalization parameters for the incoming beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_REN_IN,N_REN_IN
!
      REAL (WP)             ::  REN_R_IN,REN_I_IN
!
END MODULE RENORM_IN
!
!=======================================================================
!
MODULE SPIN_IN
!
!  This module contains the spin parameters for the incoming beam
!
      IMPLICIT NONE
!
      INTEGER               ::  ISPIN_IN,NSPIN_IN,NSPIN2_IN
      INTEGER               ::  ISFLIP_IN,IR_DIA_IN,NSTEP_IN
!
END MODULE SPIN_IN
!
!=======================================================================
!
MODULE TEST_IN
!
!  This module contains the parameters controlling the tests
!
      IMPLICIT NONE
!
      INTEGER               ::  I_TEST_IN
!
END MODULE TEST_IN
!
!=======================================================================
!
MODULE THETA_IN
!
!  This module contains the parameters controlling the polar angle
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  ITHETA_IN,NTHETA_IN
!
      REAL (WP)             ::  THETA0_IN,THETA1_IN
!
END MODULE THETA_IN

