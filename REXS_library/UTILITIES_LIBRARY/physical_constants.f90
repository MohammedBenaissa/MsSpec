!
!=======================================================================
!
MODULE CONSTANTS_P1
!
!  This module defines standard physical constants
!
!  Note : COULOMB = 1 / (4 pi epsilon_0)
!
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)   ::   BOHR   =5.2917721067E-011_WP       ! Bohr radius               ... a_0
      REAL (WP)   ::   H_BAR  =1.054571800E-034_WP        ! reduced Planck constant   ... J s
      REAL (WP)   ::   M_E    =9.10938356E-031_WP         ! electron mass             ... kg
      REAL (WP)   ::   E      =1.6021766208E-019_WP       ! charge of electron        ... C
      REAL (WP)   ::   EPS_0  =8.854187817E-012_WP        ! vacuum permittivity       ... F / m
      REAL (WP)   ::   COULOMB=8.9875517873681764E+009_WP ! Coulomb constant          ... kg m^3 / (s^4 A^2)
      REAL (WP)   ::   K_B    =1.38064852E-023_WP         ! Boltzmann constant        ... J / K
!
END MODULE CONSTANTS_P1
!
!=======================================================================
!
MODULE CONSTANTS_P2
!
!  This module defines standard physical constants
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)   ::   ALPHA  =7.2973525664E-003_WP       ! fine structure constant   ... dimensionless
      REAL (WP)   ::   HARTREE=4.359744650E-018_WP        ! Hartree energy            ... J
      REAL (WP)   ::   RYDBERG=10973731.568508E0_WP       ! Rydberg constant          ... m^{-1}
!
END MODULE CONSTANTS_P2
!
!=======================================================================
!
MODULE CONSTANTS_P3
!
!  This module defines standard physical constants
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)   ::   R_E   =2.8179403227E-015_WP        ! classical electron radius ... m
      REAL (WP)   ::   M_P   =1.672621898E-027_WP         ! proton mass               ... kg
      REAL (WP)   ::   M_N   =1.674927211E-027_WP         ! neutron mass              ... kg
      REAL (WP)   ::   C     =299792458.0E0_WP            ! speed of light in vacuum  ... m/s
      REAL (WP)   ::   G     =6.67408E-011_WP             ! constant of gravitation   ... m^3 / (kg s^2)
      REAL (WP)   ::   PLANCK=6.626070040E-034_WP         ! Planck constant           ... J s
      REAL (WP)   ::   MU_0  =1.256637061E-006_WP         ! vacuum permeability       ... N / A^2
      REAL (WP)   ::   MU_B  =9.274009994E-024_WP         ! Bohr magneton             ... J / T
      REAL (WP)   ::   MU_N  =5.050783699E-027_WP         ! nuclear magneton          ... J / T
      REAL (WP)   ::   N_A   =6.022140857E+023_WP         ! Avogadro constant         ... mol^{-1}
!
END MODULE CONSTANTS_P3
!
!=======================================================================
!
MODULE G_FACTORS
!
!  This module defines standard physical constants
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)   ::   G_E=-2.00231930436182E0_WP         ! electron g-factor         ... dimensionless
      REAL (WP)   ::   G_P=+5.585694702E0_WP              ! proton g-factor           ... dimensionless
      REAL (WP)   ::   G_N=-3.82608545E0_WP               ! neutron g-factor          ... dimensionless
!
END MODULE G_FACTORS
!
!=======================================================================
!
MODULE ENE_CHANGE
!
!  This module defines energy, etc  change factors
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)   ::   EV    =1.6021766208E-019_WP        ! electron volt             ... J
      REAL (WP)   ::   RYD   =13.605693009E0_WP           ! Rydberg energy            ... eV
      REAL (WP)   ::   HAR   =27.21138602E0_WP            ! Hartree energy            ... eV
      REAL (WP)   ::   BOHR2A=0.52917721067E0_WP          ! Bohr radius               ... Angstroem
      REAL (WP)   ::   ANG   =1.0E-010_WP                 ! Angstroem                 ... m
      REAL (WP)   ::   RY2SI =2.17987232488E-18_WP        ! conversion Ryd --> SI
!
      REAL (WP)   ::   ROOM  =273.0E0_WP                  ! room temperature          ... K
!
END MODULE ENE_CHANGE
