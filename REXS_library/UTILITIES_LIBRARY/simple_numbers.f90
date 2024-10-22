!
!=======================================================================
!
MODULE REAL_NUMBERS 
!
!  This module defines frequent real numbers
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP), PARAMETER   ::   ZERO      =   0.00E0_WP
      REAL (WP), PARAMETER   ::   ONE       =   1.00E0_WP
      REAL (WP), PARAMETER   ::   TWO       =   2.00E0_WP
      REAL (WP), PARAMETER   ::   THREE     =   3.00E0_WP
      REAL (WP), PARAMETER   ::   FOUR      =   4.00E0_WP
      REAL (WP), PARAMETER   ::   FIVE      =   5.00E0_WP
!
      REAL (WP), PARAMETER   ::   SIX       =   6.00E0_WP
      REAL (WP), PARAMETER   ::   SEVEN     =   7.00E0_WP
      REAL (WP), PARAMETER   ::   EIGHT     =   8.00E0_WP
      REAL (WP), PARAMETER   ::   NINE      =   9.00E0_WP
      REAL (WP), PARAMETER   ::   TEN       =  10.00E0_WP
!
      REAL (WP), PARAMETER   ::   TWENTY    =  20.00E0_WP
!
      REAL (WP), PARAMETER   ::   HALF      =   0.50E0_WP
      REAL (WP), PARAMETER   ::   THIRD     =   0.33333333333333333333333333333333E0_WP
      REAL (WP), PARAMETER   ::   FOURTH    =   0.25E0_WP
      REAL (WP), PARAMETER   ::   FIFTH     =   0.20E0_WP
!
      REAL (WP), PARAMETER   ::   SIXTH     =   0.16666666666666666666666666666667E0_WP
      REAL (WP), PARAMETER   ::   SEVENTH   =   0.14285714285714285714285714285714E0_WP
      REAL (WP), PARAMETER   ::   EIGHTH    =   0.125E0_WP
      REAL (WP), PARAMETER   ::   NINTH     =   0.11111111111111111111111111111111E0_WP
      REAL (WP), PARAMETER   ::   TENTH     =   0.10E0_WP
!
      REAL (WP), PARAMETER   ::   SMALL     =   1.0E-006_WP
      REAL (WP), PARAMETER   ::   TTINY     =   1.0E-030_WP
      REAL (WP), PARAMETER   ::   LARGE     =   1.0E+030_WP
      REAL (WP), PARAMETER   ::   INF       =   1.0E+300_WP
      REAL (WP), PARAMETER   ::   MIC       =   1.0E-300_WP
      REAL (WP), PARAMETER   ::   EPS       =   1.0E-010_WP
!
END MODULE REAL_NUMBERS
!
!=======================================================================
!
MODULE COMPLEX_NUMBERS 
!
!  This module defines frequent complex numbers
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
! 
   COMPLEX (WP), PARAMETER   ::   ZEROC   =   (0.0E0_WP,0.0E0_WP)
   COMPLEX (WP), PARAMETER   ::   ONEC    =   (1.0E0_WP,0.0E0_WP)
   COMPLEX (WP), PARAMETER   ::   IC      =   (0.0E0_WP,1.0E0_WP)
!
END MODULE COMPLEX_NUMBERS
 
