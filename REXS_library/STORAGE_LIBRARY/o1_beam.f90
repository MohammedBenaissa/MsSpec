!
!  These modules contain the parameters controlling the 1st outgoing beam
!
!=======================================================================
!
MODULE AMPLI_O1
!
!  This module contains the switch controlling the storage of the amplitudes
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AMP_O1
!
END MODULE AMPLI_O1
!
!=======================================================================
!
MODULE APPROXIMATIONS_O1
!
!  This module contains the variables for the approximations
!    of the first outgoing beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  N_SCAT_O1,NO_O1,I_BASIS_O1
      INTEGER               ::  IFWD_O1,NTHOUT_O1
      INTEGER               ::  IBWD_O1(NATP_M),IPW_O1,NCUT_O1
      INTEGER               ::  IPP_O1,IATT_O1,ILENGTH_O1
!
      REAL (WP)             ::  RTHFWD_O1(NATP_M),RTHBWD_O1(NATP_M)
      REAL (WP)             ::  PCTINT_O1,RLENGTH_O1
!
END MODULE APPROXIMATIONS_O1
!
!=======================================================================
!
MODULE AVER_O1
!
!  This module contains general parameters for the averaging over
!    the 1st outgoing beam direction
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AVER_O1,NDIR_O1,ICHKDIR_O1
!
      REAL (WP)             ::  ACCEPT_O1
!
END MODULE AVER_O1
!
!=======================================================================
!
MODULE ENERGY_O1
!
!  This module contains the parameters controlling the energy
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IE_O1,NE_O1
!
      REAL (WP)             ::  E0_O1,E1_O1
!
END MODULE ENERGY_O1
!
!=======================================================================
!
MODULE EXTERNAL_O1
!
!  This module contains the switch controlling the external angle values
!
      IMPLICIT NONE
!
      INTEGER               ::  I_EXT_O1
!
END MODULE EXTERNAL_O1
!
!=======================================================================
!
MODULE EXT_DIR_O1
!
!  This module contains parameters concerning the external directions
!
      IMPLICIT NONE
!
      INTEGER               ::  IDIR_O1,NSET_O1,N_POINTS_O1
!
END MODULE EXT_DIR_O1
!
!=======================================================================
!
MODULE F_TH_O1
!
!  This module contains the scattering amplitude parameters
!    for the 1st outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IFTHET_O1,NFTHET_O1
!
      REAL (WP)             ::  R0_O1,R1_O1
!
END MODULE F_TH_O1
!
!=======================================================================
!
MODULE INTEGR_O1
!
!  This module contains the switch controlling the integration
!
      IMPLICIT NONE
!
      INTEGER               ::  I_INT_O1
!
END MODULE INTEGR_O1
!
!=======================================================================
!
MODULE LINLBD_O1
!
!  This module contains the Rehr-Albers variables for 1st outgoing beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M,N_MU_M,N_NU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  LBD_O1(-N_MU_M:N_MU_M,0:N_NU_M)
      INTEGER               ::  LBDMAX_O1
      INTEGER               ::  NUMAX_O1(NATP_M)
!
END MODULE LINLBD_O1
!
!=======================================================================
!
MODULE MFP_O1
!
!  This module contains mfp value for the 1sr outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IMFP_O1
!
      REAL (WP)             ::  XMFP0_O1
!
END MODULE MFP_O1
!
!=======================================================================
!
MODULE PHI_O1
!
!  This module contains the parameters controlling the azimuthal angle
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IPHI_O1,NPHI_O1
!
      REAL (WP)             ::  PHI0_O1,PHI1_O1
!
END MODULE PHI_O1
!
!=======================================================================
!
MODULE POLAR_O1
!
!  This module contains the parameters controlling the polarization
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOL_O1
!
END MODULE POLAR_O1
!
!=======================================================================
!
MODULE RA_O1
!
!  This module contains Rehr-Albers parameters for the 1st outgoing beam
!
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  I_NO_O1,I_RA_O1,N_RA_O1(NATP_M)
!
END MODULE RA_O1
!
!=======================================================================
!
MODULE RENORM_O1
!
!  This module contains renormalization parameters for the 1st outgoing beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_REN_O1,N_REN_O1
!
      REAL (WP)             ::  REN_R_O1,REN_I_O1
!
END MODULE RENORM_O1
!
!=======================================================================
!
MODULE SPIN_O1
!
!  This module contains the spin parameters for the 1st outgoing beam
!
      IMPLICIT NONE
!
      INTEGER               ::  ISPIN_O1,NSPIN_O1,NSPIN2_O1
      INTEGER               ::  ISFLIP_O1,IR_DIA_O1,NSTEP_O1
!
END MODULE SPIN_O1
!
!=======================================================================
!
MODULE TEST_O1
!
!  This module contains the parameters controlling the tests
!
      IMPLICIT NONE
!
      INTEGER               ::  I_TEST_O1
!
END MODULE TEST_O1
!
!=======================================================================
!
MODULE THETA_O1
!
!  This module contains the parameters controlling the polar angle
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  ITHETA_O1,NTHETA_O1
!
      REAL (WP)             ::  THETA0_O1,THETA1_O1
!
END MODULE THETA_O1


