!
!  These modules contain the parameters controlling the correlation expansion
!
!=======================================================================
!
MODULE CORR_EXP
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : LINMAX,NGR_M
!
      IMPLICIT NONE
!
      COMPLEX (WP)                ::  A(LINMAX*NGR_M,LINMAX*NGR_M)
!
END MODULE CORR_EXP
!
!=======================================================================
!
MODULE Q_ARRAY
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NGR_M
!
      IMPLICIT NONE
!
      REAL (WP)                   :: Q(NGR_M)
!
END MODULE Q_ARRAY
