!
!=======================================================================
!
MODULE WAVE_FUNCTION
!
!  This module stores the core state wave function
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,     ONLY : N_MESH_M
!
      IMPLICIT NONE
!
      REAL (WP)             ::  R(N_MESH_M)
      REAL (WP)             ::  COREWF(N_MESH_M)
!
END MODULE WAVE_FUNCTION
