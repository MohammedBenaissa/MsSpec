!
!  These modules contain the parameters controlling the excited beam
!
!=======================================================================
!
MODULE AMPLI_EX
!
!  This module contains the switch controlling the storage of the amplitudes
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AMP_EX
!
END MODULE AMPLI_EX
!
!=======================================================================
!
MODULE APPROXIMATIONS_EX
!
!  This module contains the variables for the approximations
!    of the excited beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  N_SCAT_EX,NO_EX,I_BASIS_EX
      INTEGER               ::  IFWD_EX,NTHOUT_EX
      INTEGER               ::  IBWD_EX(NATP_M),IPW_EX,NCUT_EX
      INTEGER               ::  IPP_EX,IATT_EX,ILENGTH_EX
!
      REAL (WP)             ::  RTHFWD_EX(NATP_M),RTHBWD_EX(NATP_M)
      REAL (WP)             ::  PCTINT_EX,RLENGTH_EX
!
END MODULE APPROXIMATIONS_EX
!
!=======================================================================
!
MODULE ENERGY_EX
!
!  This module contains the parameters controlling the energy
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IE_EX,NE_EX
!
      REAL (WP)             ::  E0_EX,E1_EX
!
END MODULE ENERGY_EX
!
!=======================================================================
!
MODULE F_TH_EX
!
!  This module contains the scattering amplitude parameters
!    for the excited beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IFTHET_EX,NFTHET_EX
!
      REAL (WP)             ::  R0_EX,R1_EX
!
END MODULE F_TH_EX
!
!=======================================================================
!
MODULE LINLBD_EX
!
!  This module contains the Rehr-Albers variables for excited beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M,N_MU_M,N_NU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  LBD_EX(-N_MU_M:N_MU_M,0:N_NU_M)
      INTEGER               ::  LBDMAX_EX
      INTEGER               ::  NUMAX_EX(NATP_M)
!
END MODULE LINLBD_EX
!
!=======================================================================
!
MODULE MFP_EX
!
!  This module contains mfp value for the excited beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IMFP_EX
!
      REAL (WP)             ::  XMFP0_EX
!
END MODULE MFP_EX
!
!=======================================================================
!
MODULE RA_EX
!
!  This module contains Rehr-Albers parameters for the excited beam
!
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  I_NO_EX,I_RA_EX,N_RA_EX(NATP_M)
!
END MODULE RA_EX
!
!=======================================================================
!
MODULE RENORM_EX
!
!  This module contains renormalization parameters for the excited beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_REN_EX,N_REN_EX
!
      REAL (WP)             ::  REN_R_EX,REN_I_EX
!
END MODULE RENORM_EX
!
!=======================================================================
!
MODULE SPIN_EX
!
!  This module contains the spin parameters for the excited beam
!
      IMPLICIT NONE
!
      INTEGER               ::  ISPIN_EX,NSPIN_EX,NSPIN2_EX
      INTEGER               ::  ISFLIP_EX,IR_DIA_EX,NSTEP_EX
!
END MODULE SPIN_EX
!
!=======================================================================
!
MODULE TEST_EX
!
!  This module contains the parameters controlling the tests
!
      IMPLICIT NONE
!
      INTEGER               ::  I_TEST_EX
!
END MODULE TEST_EX

