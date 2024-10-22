!
!=======================================================================
!
MODULE CURRENT_AVER
!
!  This module contains general parameters for the averaging over
!    the current beam direction
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AVER,NDIR,ICHKDIR
!
      REAL (WP)             ::  ACCEPT
!
END MODULE CURRENT_AVER
!
!=======================================================================
!
MODULE CURRENT_BEAM
!
!  This module contains the parameters of the current electron beam
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,      ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER                   ::  N_SCAT,NO,I_BASIS,IFWD,NTHOUT
      INTEGER                   ::  IBWD(NATP_M)
      INTEGER                   ::  IPW,NCUT,IPP,IATTS,ILENGTH
!
      REAL (WP)                 ::  RTHFWD(NATP_M),RTHBWD(NATP_M)
      REAL (WP)                 ::  PCTINT,RLENGTH
!
END MODULE CURRENT_BEAM
!
!=======================================================================
!
MODULE CURRENT_CALC
!
!  This module contains the parameters of the current calculation
!
!
      IMPLICIT NONE
!
      INTEGER               ::  IPHI,IE,ITHETA
      INTEGER               ::  NPHI,NE,NTHETA,NEPS
      INTEGER               ::  IPOL
      INTEGER               ::  I_EXT,I_TEST
      INTEGER               ::  I_AMP,I_INT
!
END MODULE CURRENT_CALC
!
!=======================================================================
!
MODULE CURRENT_COEF_RENORM
!
!  This module contains the renormalization coefficients for the current beam
!
!
      USE DIMENSION_CODE,    ONLY : N_SCAT_M
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      COMPLEX (WP)          ::  C_REN(0:N_SCAT_M)
!
END MODULE CURRENT_COEF_RENORM
!
!=======================================================================
!
MODULE CURRENT_CLUSTER
!
!  This module contains the parameters of the current cluster
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,      ONLY :  NATM,NATP_M,NAT_EQ_M,NATCLU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  NATCLU,N_PROT,NATYP(NATM),NCHTYP(NATP_M)
      INTEGER               ::  NCORR(NAT_EQ_M,NATP_M),INEW_AT(NATCLU_M)
!
      REAL (WP)             ::  COORD(3,NATCLU_M)
!
END MODULE CURRENT_CLUSTER
!
!=======================================================================
!
MODULE CURRENT_EXT_DIR
!
!  This module contains parameters concerning the external directions
!
      IMPLICIT NONE
!
      INTEGER               ::  IDIR,NSET,N_POINTS
!
END MODULE CURRENT_EXT_DIR
!
!=======================================================================
!
MODULE CURRENT_FIXSCAN
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  N_FIXED,N_SCAN,IPH_1
!
      REAL (WP)             ::  FIX0,FIX1,SCAN0,SCAN1
!
END MODULE CURRENT_FIXSCAN
!
!=======================================================================
!
MODULE CURRENT_F_TH
!
!  This module contains the scattering amplitude parameters
!    for the current beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IFTHET,NFTHET
!
      REAL (WP)             ::  R0,R1
!
END MODULE CURRENT_F_TH
!
!=======================================================================
!
MODULE CURRENT_FINA_VAL
!
!  This module contains the final values for the current beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)             ::  PHI1,THETA1,E1
!
END MODULE CURRENT_FINA_VAL
!
!=======================================================================
!
MODULE CURRENT_INIT_VAL
!
!  This module contains the initial values for the current beam
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)             ::  PHI0,THETA0,E0
      REAL (WP)             ::  PHILUM,ELUM,THLUM
!
END MODULE CURRENT_INIT_VAL
!
!=======================================================================
!
MODULE CURRENT_LBDM_STORE
!
!  This module contains the Rehr-Albers storage for the current beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : N_SCAT_M
!
      IMPLICIT NONE
!
      INTEGER               ::  LBDM(2,N_SCAT_M)
!
END MODULE CURRENT_LBDM_STORE
!
!=======================================================================
!
MODULE CURRENT_LINLBD
!
!  This module contains the Rehr-Albers variables for the current beam
!
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M,N_MU_M,N_NU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  LBD(-N_MU_M:N_MU_M,0:N_NU_M)
      INTEGER               ::  LBDMAX
      INTEGER               ::  NUMAX(NATP_M)
!
END MODULE CURRENT_LINLBD
!
!=======================================================================
!
MODULE CURRENT_MEAN_FREE_PATH
!
!  This module contains the mfp variables for the current beam
!
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  IMFP,NZ_A
!
      REAL (WP)             ::  XMT_A,RHOT_A,XMFP0
!
END MODULE CURRENT_MEAN_FREE_PATH
!
!=======================================================================
!
MODULE CURRENT_REHR_ALBERS
!
!  This module contains the Rehr-Albers parameters of the current cluster
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,      ONLY :  NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  I_NO,I_RA
      INTEGER               ::  N_RA(NATP_M)
!
END MODULE CURRENT_REHR_ALBERS
!
!=======================================================================
!
MODULE CURRENT_RENORM
!
!  This module contains general parameters for the renormalization of
!    the current beam direction
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_REN,N_REN
!
      REAL (WP)             ::  REN_R,REN_I
!
END MODULE CURRENT_RENORM
!
!=======================================================================
!
MODULE CURRENT_SPIN
!
!  This module contains the spin parameters of the current beam
!
      IMPLICIT NONE
!
      INTEGER               ::  ISPIN,NSPIN,NSPIN2
      INTEGER               ::  ISFLIP,IR_DIA,NSTEP
!
END MODULE CURRENT_SPIN
!
!=======================================================================
!
MODULE CURRENT_T_MATRIX
!
!  This module contains the parameters of the T-matrix elements
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NATM,NE_M,NT_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTC,ITL,LMAX(NATM,NE_M)
!
      REAL (WP)             ::  VK2(NE_M)
!
      COMPLEX (WP)          ::  VK(NE_M)
      COMPLEX (WP)          ::  TL(0:NT_M,4,NATM,NE_M)
!
END MODULE CURRENT_T_MATRIX
