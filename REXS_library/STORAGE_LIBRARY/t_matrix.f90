!
!=======================================================================
!
MODULE TLM_IN
!
!  This module contains the parameters for the incoming beam
!    T-matrix elements in the spherical wave representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NATM,NE_M,NT_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTC_IN,ITL_IN,LMAX_IN(NATM,NE_M)
!
      REAL (WP)             ::  K2L_IN(NE_M)
!
      COMPLEX (WP)          ::  TL_IN(0:NT_M,4,NATM,NE_M)
      COMPLEX (WP)          ::  KL_IN(NE_M)
!
END MODULE TLM_IN
!
!=======================================================================
!
MODULE TLM_EX
!
!  This module contains the parameters for the excited beam
!    T-matrix elements in the spherical wave representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NATM,NE_M,NT_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTC_EX,ITL_EX,LMAX_EX(NATM,NE_M)
!
      REAL (WP)             ::  K2L_EX(NE_M)
!
      COMPLEX (WP)          ::  TL_EX(0:NT_M,4,NATM,NE_M)
      COMPLEX (WP)          ::  KL_EX(NE_M)
!
END MODULE TLM_EX
!
!=======================================================================
!
MODULE TLM_O1
!
!  This module contains the parameters for the 1st outgoing beam
!    T-matrix elements in the spherical wave representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NATM,NE_M,NT_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTC_O1,ITL_O1,LMAX_O1(NATM,NE_M)
!
      REAL (WP)             ::  K2L_O1(NE_M)
!
      COMPLEX (WP)          ::  TL_O1(0:NT_M,4,NATM,NE_M)
      COMPLEX (WP)          ::  KL_O1(NE_M)
!
END MODULE TLM_O1
!
!=======================================================================
!
MODULE TLM_O2
!
!  This module contains the parameters for the 2nd outgoing beam
!    T-matrix elements in the spherical wave representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NATM,NE_M,NT_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTC_O2,ITL_O2,LMAX_O2(NATM,NE_M)
!
      REAL (WP)             ::  K2L_O2(NE_M)
!
      COMPLEX (WP)          ::  TL_O2(0:NT_M,4,NATM,NE_M)
      COMPLEX (WP)          ::  KL_O2(NE_M)
!
END MODULE TLM_O2
!
!=======================================================================
!
MODULE TLM_O3
!
!  This module contains the parameters for the 3rd outgoing beam
!    T-matrix elements in the spherical wave representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NATM,NE_M,NT_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTC_O3,ITL_O3,LMAX_O3(NATM,NE_M)
!
      REAL (WP)             ::  K2L_O3(NE_M)
!
      COMPLEX (WP)          ::  TL_O3(0:NT_M,4,NATM,NE_M)
      COMPLEX (WP)          ::  KL_O3(NE_M)
!
END MODULE TLM_O3
!
!=======================================================================
!
MODULE R_MESH
!
!  This module provides the r-mesh for the impact parameter
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : NATM,N_MESH_M,NE_M
!
      IMPLICIT NONE
!
      REAL (WP)             ::  B_MESH(NATM,N_MESH_M,NE_M)
!
END MODULE R_MESH
!
!=======================================================================
!
MODULE TBP_IN
!
!  This module contains the parameters for the incoming beam
!    T-matrix elements in the impact parameter representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : N_MESH_M,NATM,NE_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTB_IN,ITB_IN,NBMAX_IN(NATM,NE_M)
!
      REAL (WP)             ::  K2B_IN(NE_M)
!
      COMPLEX (WP)          ::  KB_IN(NE_M)
      COMPLEX (WP)          ::  TB_IN(N_MESH_M,4,NATM,NE_M)
!
END MODULE TBP_IN
!
!=======================================================================
!
MODULE TBP_EX
!
!  This module contains the parameters for the excited beam
!    T-matrix elements in the impact parameter representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : N_MESH_M,NATM,NE_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTB_EX,ITB_EX,NBMAX_EX(NATM,NE_M)
!
      REAL (WP)             ::  K2B_EX(NE_M)
!
      COMPLEX (WP)          ::  KB_EX(NE_M)
      COMPLEX (WP)          ::  TB_EX(N_MESH_M,4,NATM,NE_M)
!
END MODULE TBP_EX
!
!=======================================================================
!
MODULE TBP_O1
!
!  This module contains the parameters for the 1st outgoing beam
!    T-matrix elements in the impact parameter representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : N_MESH_M,NATM,NE_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTB_O1,ITB_O1,NBMAX_O1(NATM,NE_M)
!
      REAL (WP)             ::  K2B_O1(NE_M)
!
      COMPLEX (WP)          ::  KB_O1(NE_M)
      COMPLEX (WP)          ::  TB_O1(N_MESH_M,4,NATM,NE_M)
!
END MODULE TBP_O1
!
!=======================================================================
!
MODULE TBP_O2
!
!  This module contains the parameters for the 2nd outgoing beam
!    T-matrix elements in the impact parameter representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : N_MESH_M,NATM,NE_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTB_O2,ITB_O2,NBMAX_O2(NATM,NE_M)
!
      REAL (WP)             ::  K2B_O2(NE_M)
!
      COMPLEX (WP)          ::  KB_O2(NE_M)
      COMPLEX (WP)          ::  TB_O2(N_MESH_M,4,NATM,NE_M)
!
END MODULE TBP_O2
!
!=======================================================================
!
MODULE TBP_O3
!
!  This module contains the parameters for the 3rd outgoing beam
!    T-matrix elements in the impact parameter representation
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : N_MESH_M,NATM,NE_M
!
      IMPLICIT NONE
!
      INTEGER               ::  IPOTB_O3,ITB_O3,NBMAX_O3(NATM,NE_M)
!
      REAL (WP)             ::  K2B_O3(NE_M)
!
      COMPLEX (WP)          ::  KB_O3(NE_M)
      COMPLEX (WP)          ::  TB_O3(N_MESH_M,4,NATM,NE_M)
!
END MODULE TBP_O3
!
!=======================================================================
!
MODULE TL_TRUNC
!
!  This module contains the parameters for the truncation of the t_l
!
!
      IMPLICIT NONE
!
      INTEGER               ::  ITRTL_IN,ITRTL_EX
      INTEGER               ::  ITRTL_O1,ITRTL_O2,ITRTL_O3
!
END MODULE TL_TRUNC
