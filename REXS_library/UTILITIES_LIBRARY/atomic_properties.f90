!
!=======================================================================
!
MODULE ATOMIC_PROPERTIES
!
!  This module contains physical properties of the chemical elements, 
!    namely: The atomic mass, the density of all the 
!    elements,the various radii available in the literature, 
!    Debye temperatures and bulk, shear, Young moduli and Poisson ratio,  
!    the work function and the valence (integer and noninteger),
!    the electric resistivity and the lattice constants, 
!    the magnetic type.
!
!  Units and main references:  
!
!     atomic mass         --> g/mol         http://periodictable.com/Properties/A/AtomicMass.html
!     density             --> g/cm^3        http://periodictable.com/Properties/A/Density.html
!     atomic radius       --> Angstroem     http://periodictable.com/Properties/A/AtomicRadius.html
!     ionic radius        --> Angstroem     http://chemistry-reference.com/images/ionic%20radius%20table.jpg
!     covalent radius     --> Angstroem     http://periodictable.com/Properties/A/CovalentRadius.html
!     Wigner-Seitz radius --> Angstroem     "Solid State Physics", Ashcroft-Mermin p.5
!     Debye temperature   --> Kelvin        http://www.knowledgedoor.com/2/elements_handbook/debye_temperature.html
!     bulk modulus        --> GPa           http://periodictable.com/Properties/A/BulkModulus.html
!     shear modulus       --> GPa           http://periodictable.com/Properties/A/ShearModulus.html
!     Young modulus       --> GPa           http://periodictable.com/Properties/A/YoungModulus.html
!     Poisson ratio       --> dimensionless http://periodictable.com/Properties/A/PoissonRatio.v.html
!     work function       --> eV            S. Halas, Materials Science-Poland 24, 951 (2006)
!     valence             --> dimensionless https://ptable.com/#Property/Valence
!     noninteger valence  --> dimensionless "The Physics of Solids. Essentials and Beyond"
!     resistivity         --> m Ohm         http://periodictable.com/Properties/A/Resistivity.an.log.html
!     lattice constants   --> Angstroem     http://periodictable.com/Properties/A/LatticeConstants.html
!     crystal structure   --> dimensionless http://periodictable.com/Properties/A/CrystalStructure.html
!     magnetic type       --> dimensionless http://periodictable.com/Properties/A/MagneticType.html
!
!  Additional references: 
!
!     * "Fundamentals of the Physics of Solids", Vol1, Solyom p.596 (Debye)
!     * A. Ruban et al, J. Mol. Cat. A: Chem 115, 421-429 (1997) (WS)
!     * "The Physics of Solids. Essentials and Beyond", E.N. Economou, Springer, table 4.2 p. 89 (WS)
!                                                                                table 4.1 p. 85 (NIBV)
!     * https://i.stack.imgur.com/1mEVV.png (Work function)
!     * http://periodictable.com/Properties/A/Valence.al.html (Valence)
!
!
!  Value Z = 0 added for empty spheres (ES). The values entered in this 
!    case are arbitrary and set to the corresponding Z = 1 value 
!    divided by 1836 (the ratio of the mass of the proton and electron).
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified :  6 Aug 2020
!
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
!  1) Chemical symbol
!
      CHARACTER (LEN = 2), DIMENSION(0:105), PARAMETER  ::  CHEM_SY = (/        & !
                'ES',' H','He','Li','Be',' B',' C',' N',' O',                   & !
                ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',                   & !
                'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',                   & !
                'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',                   & !
                'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',                   & !
                'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',                   & !
                'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',                   & !
                'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',                   & !
                'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',                   & !
                'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',                   & !
                'Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf',                   & !
                'Es','Fm','Md','No','Lr','Rf','Db'                              & !
                                                                      /)          !
!
!  2) Atomic mass --> g/mol 
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  MASS_AT = (/                   & !
         0.00055E0_WP,1.0079E0_WP,4.0026E0_WP,6.941E0_WP,9.0122E0_WP,10.81E0_WP,& ! ES |  H | He | Li | Be |  B |
        12.011E0_WP, 14.0067E0_WP,15.9994E0_WP,18.998403E0_WP,20.179E0_WP ,     & !  C |  N |  O |  F | Ne |
        22.98977E0_WP,24.305E0_WP,26.98154E0_WP,28.0855E0_WP,30.97376E0_WP,     & ! Na | Mg | Al | Si |  P |
        32.06E0_WP,35.453E0_WP,39.948E0_WP,39.0983E0_WP,40.08E0_WP,             & !  S | Cl | Ar |  K | Ca |
        44.9559E0_WP,47.88E0_WP,50.9415E0_WP,51.996E0_WP,54.9380E0_WP,          & ! Sc | Ti |  V | Cr | Mn |
        55.847E0_WP,58.9332E0_WP,58.69E0_WP,63.546E0_WP,65.38E0_WP,             & ! Fe | Co | Ni | Cu | Zn |
        69.72E0_WP,72.59E0_WP,74.9216E0_WP,78.96E0_WP,79.904E0_WP,              & ! Ga | Ge | As | Se | Br |
        83.80E0_WP,85.4678E0_WP,87.62E0_WP,88.9059E0_WP,91.22E0_WP,             & ! Kr | Rb | Sr |  Y | Zr |
        92.9064E0_WP,95.94E0_WP,98.E0_WP, 101.07E0_WP,102.9055E0_WP,            & ! Nb | Mo | Tc | Ru | Rh |
       106.42E0_WP,107.8682E0_WP,112.41E0_WP,114.82E0_WP,118.69E0_WP,           & ! Pd | Ag | Cd | In | Sn |
       121.75E0_WP,127.60E0_WP,126.9045E0_WP,131.29E0_WP,132.9054E0_WP,         & ! Sb | Te |  I | Xe | Cs |
       137.33E0_WP,138.9055E0_WP,140.12E0_WP,140.9077E0_WP,144.24E0_WP,         & ! Ba | La | Ce | Pr | Nd |
       145.E0_WP,150.36E0_WP,151.96E0_WP,157.25E0_WP,158.9254E0_WP,             & ! Pm | Sm | Eu | Gd | Tb |
       162.50E0_WP,164.9304E0_WP,167.26E0_WP,168.9342E0_WP,173.04E0_WP,         & ! Dy | Ho | Er | Tm | Yb |
       174.967E0_WP,178.49E0_WP,180.9479E0_WP,183.85E0_WP,186.207E0_WP,         & ! Lu | Hf | Ta |  W | Re |
       190.2E0_WP,192.22E0_WP,195.08E0_WP,196.9665E0_WP,200.59E0_WP,            & ! Os | Ir | Pt | Au | Hg |
       204.383E0_WP,207.2E0_WP,208.9804E0_WP,209.E0_WP,210.E0_WP,               & ! Tl | Pb | Bi | Po | At |
       222.E0_WP,223.E0_WP,226.0254E0_WP,227.0278E0_WP,232.0381E0_WP,           & ! Rn | Fr | Ra | Ac | Th |
       231.0359E0_WP,238.0289E0_WP,237.0482E0_WP,244.054E0_WP,243.061E0_WP,     & ! Pa |  U | Np | Pu | Am |
       247.070E0_WP,247.070E0_WP,251.080E0_WP,254.E0_WP,257.095E0_WP,           & ! Cm | Bk | Cf | Es | Fm |
       258.1E0_WP,259.101E0_WP,262.E0_WP,261.E0_WP,262.E0_WP                    & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
!  3) Atomic density --> g/cm^3
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  ATOM_DE = (/                   & !
            0.00005E0_WP,0.0899E0_WP,0.122E0_WP,0.535E0_WP,1.848E0_WP,2.46E0_WP,& ! ES |  H | He | Li | Be |  B |
            2.26E0_WP,1.251E0_WP,1.429E0_WP,1.696E0_WP,0.9E0_WP,                & !  C |  N |  O |  F | Ne |
            0.968E0_WP,1.738E0_WP,2.7E0_WP,2.33E0_WP,1.823E0_WP,                & ! Na | Mg | Al | Si |  P |
            1.96E0_WP,3.214E0_WP,1.784E0_WP,0.856E0_WP,1.55E0_WP,               & !  S | Cl | Ar |  K | Ca |
            2.985E0_WP,4.507E0_WP,6.11E0_WP,7.19E0_WP,7.47E0_WP,                & ! Sc | Ti |  V | Cr | Mn |
            7.874E0_WP,8.9E0_WP,8.9E0_WP,8.908E0_WP,7.14E0_WP,                  & ! Fe | Co | Ni | Cu | Zn |
            5.904E0_WP,5.323E0_WP,5.727E0_WP,4.819E0_WP,3.12E0_WP,              & ! Ga | Ge | As | Se | Br |
            3.75E0_WP,1.532E0_WP,2.63E0_WP,4.472E0_WP,6.511E0_WP,               & ! Kr | Rb | Sr |  Y | Zr |
            8.57E0_WP,10.28E0_WP,11.5E0_WP,12.37E0_WP,12.45E0_WP,               & ! Nb | Mo | Tc | Ru | Rh |
            12.023E0_WP,10.49E0_WP,8.65E0_WP,7.31E0_WP,7.31E0_WP,               & ! Pd | Ag | Cd | In | Sn |
            6.697E0_WP,6.24E0_WP,4.94E0_WP,5.9E0_WP,1.879E0_WP,                 & ! Sb | Te |  I | Xe | Cs |
            3.51E0_WP,6.146E0_WP,6.689E0_WP,6.64E0_WP,7.01E0_WP,                & ! Ba | La | Ce | Pr | Nd |
            7.264E0_WP,7.353E0_WP,5.244E0_WP,7.901E0_WP,8.219E0_WP,             & ! Pm | Sm | Eu | Gd | Tb |
            8.551E0_WP,8.795E0_WP,9.066E0_WP,9.32E0_WP,6.57E0_WP,               & ! Dy | Ho | Er | Tm | Yb |
            9.841E0_WP,13.31E0_WP,16.65E0_WP,19.25E0_WP,21.02E0_WP,             & ! Lu | Hf | Ta |  W | Re |
            22.59E0_WP,22.56E0_WP,21.45E0_WP,19.3E0_WP,13.534E0_WP,             & ! Os | Ir | Pt | Au | Hg |
            11.85E0_WP,11.34E0_WP,9.78E0_WP,9.196E0_WP,0.0E0_WP,                & ! Tl | Pb | Bi | Po | At |
            9.73E0_WP,0.0E0_WP,5.E0_WP,10.07E0_WP,11.724E0_WP,                  & ! Rn | Fr | Ra | Ac | Th |
            15.37E0_WP,19.05E0_WP,20.45E0_WP,19.816E0_WP,13.67E0_WP,            & ! Pa |  U | Np | Pu | Am |
            13.51E0_WP,14.78E0_WP,15.1E0_WP,0.0E0_WP,0.0E0_WP,                  & ! Cm | Bk | Cf | Es | Fm |
            0.0E0_WP,0.0E0_WP,0.0E0_WP,0.0E0_WP,0.0E0_WP                        & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
!  4) Atomic radius --> Angstroem
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  ATOM_RD = (/                   & !
           0.000289E0_WP,0.53E0_WP,0.31E0_WP,1.67E0_WP,1.12E0_WP,0.87E0_WP,     & ! ES |  H | He | Li | Be |  B |
                    0.67E0_WP,0.56E0_WP,0.48E0_WP,0.42E0_WP,0.38E0_WP,          & !  C |  N |  O |  F | Ne |
                    1.90E0_WP,1.45E0_WP,1.18E0_WP,1.11E0_WP,0.98E0_WP,          & ! Na | Mg | Al | Si |  P |
                    0.88E0_WP,0.79E0_WP,0.71E0_WP,2.43E0_WP,1.94E0_WP,          & !  S | Cl | Ar |  K | Ca |
                    1.84E0_WP,1.76E0_WP,1.71E0_WP,1.66E0_WP,1.61E0_WP,          & ! Sc | Ti |  V | Cr | Mn |
                    1.56E0_WP,1.52E0_WP,1.49E0_WP,1.45E0_WP,1.42E0_WP,          & ! Fe | Co | Ni | Cu | Zn |
                    1.36E0_WP,1.25E0_WP,1.14E0_WP,1.03E0_WP,0.94E0_WP,          & ! Ga | Ge | As | Se | Br |
                    0.88E0_WP,2.65E0_WP,2.19E0_WP,2.12E0_WP,2.06E0_WP,          & ! Kr | Rb | Sr |  Y | Zr |
                    1.98E0_WP,1.90E0_WP,1.83E0_WP,1.78E0_WP,1.73E0_WP,          & ! Nb | Mo | Tc | Ru | Rh |
                    1.69E0_WP,1.65E0_WP,1.61E0_WP,1.56E0_WP,1.45E0_WP,          & ! Pd | Ag | Cd | In | Sn |
                    1.33E0_WP,1.23E0_WP,1.15E0_WP,1.08E0_WP,2.98E0_WP,          & ! Sb | Te |  I | Xe | Cs |
                    2.53E0_WP,1.95E0_WP,1.85E0_WP,2.47E0_WP,2.06E0_WP,          & ! Ba | La | Ce | Pr | Nd |
                    2.05E0_WP,2.38E0_WP,2.31E0_WP,2.33E0_WP,2.25E0_WP,          & ! Pm | Sm | Eu | Gd | Tb |
                    2.28E0_WP,1.75E0_WP,2.26E0_WP,2.22E0_WP,2.22E0_WP,          & ! Dy | Ho | Er | Tm | Yb |
                    2.17E0_WP,2.08E0_WP,2.00E0_WP,1.93E0_WP,1.88E0_WP,          & ! Lu | Hf | Ta |  W | Re |
                    1.85E0_WP,1.80E0_WP,1.77E0_WP,1.74E0_WP,1.71E0_WP,          & ! Os | Ir | Pt | Au | Hg |
                    1.56E0_WP,1.54E0_WP,1.43E0_WP,1.35E0_WP,1.38E0_WP,          & ! Tl | Pb | Bi | Po | At |
                    1.20E0_WP,3.48E0_WP,2.15E0_WP,1.95E0_WP,1.80E0_WP,          & ! Rn | Fr | Ra | Ac | Th |
                    1.80E0_WP,1.75E0_WP,1.75E0_WP,1.75E0_WP,1.75E0_WP,          & ! Pa |  U | Np | Pu | Am |
                    0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,          & ! Cm | Bk | Cf | Es | Fm |
                    0.00E0_WP,0.00E0_WP,0.00E0_WP,1.31E0_WP,1.26E0_WP           & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
!  5) Ionic radius --> Angstroem
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  IONI_RD = (/                   & !
                0.00084E0_WP,1.54E0_WP,0.00E0_WP,0.76E0_WP,0.34E0_WP,0.23E0_WP, & ! ES |  H | He | Li | Be |  B |
                    2.60E0_WP,1.46E0_WP,1.40E0_WP,1.33E0_WP,1.12E0_WP,          & !  C |  N |  O |  F | Ne |
                    1.02E0_WP,0.72E0_WP,0.54E0_WP,2.71E0_WP,2.12E0_WP,          & ! Na | Mg | Al | Si |  P |
                    1.84E0_WP,1.81E0_WP,1.54E0_WP,1.38E0_WP,1.00E0_WP,          & !  S | Cl | Ar |  K | Ca |
                    0.83E0_WP,0.80E0_WP,0.72E0_WP,0.84E0_WP,0.91E0_WP,          & ! Sc | Ti |  V | Cr | Mn |
                    0.82E0_WP,0.82E0_WP,0.78E0_WP,0.96E0_WP,0.83E0_WP,          & ! Fe | Co | Ni | Cu | Zn |
                    1.13E0_WP,0.90E0_WP,0.69E0_WP,0.69E0_WP,1.96E0_WP,          & ! Ga | Ge | As | Se | Br |
                    1.69E0_WP,1.52E0_WP,1.18E0_WP,1.06E0_WP,1.09E0_WP,          & ! Kr | Rb | Sr |  Y | Zr |
                    0.74E0_WP,0.92E0_WP,0.95E0_WP,0.77E0_WP,0.86E0_WP,          & ! Nb | Mo | Tc | Ru | Rh |
                    0.86E0_WP,1.13E0_WP,1.14E0_WP,1.32E0_WP,0.93E0_WP,          & ! Pd | Ag | Cd | In | Sn |
                    0.89E0_WP,2.11E0_WP,2.10E0_WP,1.90E0_WP,1.65E0_WP,          & ! Sb | Te |  I | Xe | Cs |
                    1.43E0_WP,1.22E0_WP,1.07E0_WP,1.06E0_WP,1.04E0_WP,          & ! Ba | La | Ce | Pr | Nd |
                    1.06E0_WP,1.11E0_WP,1.12E0_WP,0.97E0_WP,0.93E0_WP,          & ! Pm | Sm | Eu | Gd | Tb |
                    0.91E0_WP,0.89E0_WP,0.89E0_WP,0.87E0_WP,1.13E0_WP,          & ! Dy | Ho | Er | Tm | Yb |
                    0.85E0_WP,0.84E0_WP,0.72E0_WP,0.68E0_WP,0.72E0_WP,          & ! Lu | Hf | Ta |  W | Re |
                    0.89E0_WP,0.89E0_WP,0.85E0_WP,1.37E0_WP,1.27E0_WP,          & ! Os | Ir | Pt | Au | Hg |
                    1.49E0_WP,1.32E0_WP,0.96E0_WP,0.65E0_WP,2.27E0_WP,          & ! Tl | Pb | Bi | Po | At |
                    0.00E0_WP,1.80E0_WP,1.52E0_WP,1.18E0_WP,1.01E0_WP,          & ! Rn | Fr | Ra | Ac | Th |
                    1.13E0_WP,1.03E0_WP,1.10E0_WP,1.08E0_WP,1.07E0_WP,          & ! Pa |  U | Np | Pu | Am |
                    1.19E0_WP,1.18E0_WP,1.17E0_WP,0.00E0_WP,0.00E0_WP,          & ! Cm | Bk | Cf | Es | Fm |
                    0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP           & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
!  6) Covalent radius --> Angstroem
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  COVA_RD = (/                   & !
            0.00017E0_WP,0.31E0_WP,0.28E0_WP,1.28E0_WP,0.96E0_WP,0.83E0_WP,     & ! ES |  H | He | Li | Be |  B |
                    0.76E0_WP,0.71E0_WP,0.66E0_WP,0.57E0_WP,0.58E0_WP,          & !  C |  N |  O |  F | Ne |
                    1.66E0_WP,1.41E0_WP,1.21E0_WP,1.11E0_WP,1.07E0_WP,          & ! Na | Mg | Al | Si |  P |
                    1.05E0_WP,1.02E0_WP,1.06E0_WP,2.03E0_WP,1.76E0_WP,          & !  S | Cl | Ar |  K | Ca |
                    1.70E0_WP,1.60E0_WP,1.53E0_WP,1.39E0_WP,1.39E0_WP,          & ! Sc | Ti |  V | Cr | Mn |
                    1.32E0_WP,1.26E0_WP,1.24E0_WP,1.32E0_WP,1.22E0_WP,          & ! Fe | Co | Ni | Cu | Zn |
                    1.22E0_WP,1.20E0_WP,1.19E0_WP,1.20E0_WP,1.20E0_WP,          & ! Ga | Ge | As | Se | Br |
                    1.16E0_WP,2.10E0_WP,1.95E0_WP,1.90E0_WP,1.75E0_WP,          & ! Kr | Rb | Sr |  Y | Zr |
                    1.64E0_WP,1.54E0_WP,1.47E0_WP,1.46E0_WP,1.42E0_WP,          & ! Nb | Mo | Tc | Ru | Rh |
                    1.39E0_WP,1.45E0_WP,1.44E0_WP,1.42E0_WP,1.39E0_WP,          & ! Pd | Ag | Cd | In | Sn |
                    1.39E0_WP,1.38E0_WP,1.39E0_WP,1.40E0_WP,2.44E0_WP,          & ! Sb | Te |  I | Xe | Cs |
                    2.15E0_WP,2.07E0_WP,2.04E0_WP,2.03E0_WP,2.01E0_WP,          & ! Ba | La | Ce | Pr | Nd |
                    1.99E0_WP,1.98E0_WP,1.98E0_WP,1.96E0_WP,1.94E0_WP,          & ! Pm | Sm | Eu | Gd | Tb |
                    1.92E0_WP,1.92E0_WP,1.89E0_WP,1.90E0_WP,1.87E0_WP,          & ! Dy | Ho | Er | Tm | Yb |
                    1.87E0_WP,1.75E0_WP,1.70E0_WP,1.62E0_WP,1.51E0_WP,          & ! Lu | Hf | Ta |  W | Re |
                    1.44E0_WP,1.41E0_WP,1.36E0_WP,1.36E0_WP,1.32E0_WP,          & ! Os | Ir | Pt | Au | Hg |
                    1.45E0_WP,1.46E0_WP,1.48E0_WP,1.40E0_WP,1.50E0_WP,          & ! Tl | Pb | Bi | Po | At |
                    1.50E0_WP,2.60E0_WP,2.12E0_WP,2.15E0_WP,2.06E0_WP,          & ! Rn | Fr | Ra | Ac | Th |
                    2.00E0_WP,1.96E0_WP,1.90E0_WP,1.87E0_WP,1.80E0_WP,          & ! Pa |  U | Np | Pu | Am |
                    1.69E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,          & ! Cm | Bk | Cf | Es | Fm |
                    0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP           & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
!  7) Wigner-Seitz radius --> Angstroem
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  WISE_RD = (/                   & !
            0.00E0_WP,1.61E0_WP,0.00E0_WP,1.72E0_WP,0.99E0_WP,0.00E0_WP,        & ! ES |  H | He | Li | Be |  B |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & !  C |  N |  O |  F | Ne |
                2.08E0_WP,1.41E0_WP,1.10E0_WP,0.00E0_WP,0.00E0_WP,              & ! Na | Mg | Al | Si |  P |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,2.57E0_WP,1.73E0_WP,              & !  S | Cl | Ar |  K | Ca |
                1.25E0_WP,1.02E0_WP,0.94E0_WP,0.98E0_WP,1.13E0_WP,              & ! Sc | Ti |  V | Cr | Mn |
                1.12E0_WP,1.10E0_WP,0.95E0_WP,1.41E0_WP,1.22E0_WP,              & ! Fe | Co | Ni | Cu | Zn |
                1.16E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Ga | Ge | As | Se | Br |
                0.00E0_WP,2.75E0_WP,1.89E0_WP,1.38E0_WP,1.16E0_WP,              & ! Kr | Rb | Sr |  Y | Zr |
                1.63E0_WP,0.85E0_WP,0.95E0_WP,0.93E0_WP,1.03E0_WP,              & ! Nb | Mo | Tc | Ru | Rh |
                1.05E0_WP,1.60E0_WP,1.37E0_WP,1.27E0_WP,1.17E0_WP,              & ! Pd | Ag | Cd | In | Sn |
                1.34E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,2.98E0_WP,              & ! Sb | Te |  I | Xe | Cs |
                1.96E0_WP,1.41E0_WP,1.40E0_WP,1.40E0_WP,1.40E0_WP,              & ! Ba | La | Ce | Pr | Nd |
                0.00E0_WP,1.38E0_WP,1.79E0_WP,1.38E0_WP,1.36E0_WP,              & ! Pm | Sm | Eu | Gd | Tb |
                1.36E0_WP,1.35E0_WP,1.34E0_WP,1.34E0_WP,1.70E0_WP,              & ! Dy | Ho | Er | Tm | Yb |
                1.32E0_WP,1.10E0_WP,0.95E0_WP,0.86E0_WP,0.84E0_WP,              & ! Lu | Hf | Ta |  W | Re |
                0.83E0_WP,0.94E0_WP,1.06E0_WP,1.59E0_WP,1.40E0_WP,              & ! Os | Ir | Pt | Au | Hg |
                1.31E0_WP,1.22E0_WP,1.19E0_WP,0.00E0_WP,0.00E0_WP,              & ! Tl | Pb | Bi | Po | At |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,1.44E0_WP,1.25E0_WP,              & ! Rn | Fr | Ra | Ac | Th |
                1.14E0_WP,1.07E0_WP,1.95E0_WP,1.06E0_WP,0.00E0_WP,              & ! Pa |  U | Np | Pu | Am |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Cm | Bk | Cf | Es | Fm |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
!  8) Debye temperature at 0 K --> Kelvin
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  DEBY_TE = (/                   & !
             0.00E0_WP, 122.E0_WP,  26.E0_WP, 344.E0_WP,1481.E0_WP,1480.E0_WP,  & ! ES |  H | He | Li | Be |  B |
            2230.E0_WP,  68.E0_WP,  91.E0_WP, 0.00E0_WP,  75.E0_WP,             & !  C |  N |  O |  F | Ne |
             157.E0_WP, 403.E0_WP, 433.E0_WP, 645.E0_WP, 193.E0_WP,             & ! Na | Mg | Al | Si |  P |
             250.E0_WP, 115.E0_WP,  92.E0_WP,  91.E0_WP, 229.E0_WP,             & !  S | Cl | Ar |  K | Ca |
             346.E0_WP, 420.E0_WP, 399.E0_WP, 606.E0_WP, 409.E0_WP,             & ! Sc | Ti |  V | Cr | Mn |
             477.E0_WP, 460.E0_WP, 477.E0_WP, 347.E0_WP, 329.E0_WP,             & ! Fe | Co | Ni | Cu | Zn |
             325.E0_WP, 373.E0_WP, 282.E0_WP, 153.E0_WP, 0.00E0_WP,             & ! Ga | Ge | As | Se | Br |
              72.E0_WP,  57.E0_WP, 147.E0_WP, 248.E0_WP, 290.E0_WP,             & ! Kr | Rb | Sr |  Y | Zr |
             276.E0_WP, 423.E0_WP, 454.E0_WP, 555.E0_WP, 512.E0_WP,             & ! Nb | Mo | Tc | Ru | Rh |
             271.E0_WP, 227.E0_WP, 210.E0_WP, 112.E0_WP, 199.E0_WP,             & ! Pd | Ag | Cd | In | Sn |
             220.E0_WP, 152.E0_WP, 109.E0_WP,  64.E0_WP,  40.E0_WP,             & ! Sb | Te |  I | Xe | Cs |
             111.E0_WP, 150.E0_WP, 179.E0_WP, 152.E0_WP, 163.E0_WP,             & ! Ba | La | Ce | Pr | Nd |
             158.E0_WP, 169.E0_WP, 118.E0_WP, 182.E0_WP, 176.E0_WP,             & ! Pm | Sm | Eu | Gd | Tb |
             183.E0_WP, 190.E0_WP, 188.E0_WP, 200.E0_WP, 118.E0_WP,             & ! Dy | Ho | Er | Tm | Yb |
             183.E0_WP, 272.E0_WP, 246.E0_WP, 383.E0_WP, 416.E0_WP,             & ! Lu | Hf | Ta |  W | Re |
             467.E0_WP, 420.E0_WP, 237.E0_WP, 162.E0_WP,  72.E0_WP,             & ! Os | Ir | Pt | Au | Hg |
              79.E0_WP, 105.E0_WP, 120.E0_WP,  81.E0_WP, 0.00E0_WP,             & ! Tl | Pb | Bi | Po | At |
             0.00E0_WP,  39.E0_WP,  89.E0_WP, 124.E0_WP, 160.E0_WP,             & ! Rn | Fr | Ra | Ac | Th |
             185.E0_WP, 248.E0_WP, 259.E0_WP, 206.E0_WP, 121.E0_WP,             & ! Pa |  U | Np | Pu | Am |
             123.E0_WP, 0.00E0_WP, 0.00E0_WP, 0.00E0_WP, 0.00E0_WP,             & ! Cm | Bk | Cf | Es | Fm |
             0.00E0_WP, 0.00E0_WP, 0.00E0_WP, 0.00E0_WP, 0.00E0_WP              & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
!  9) Bulk modulus --> GPa
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  BULK_MD = (/                   & !
      0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,11.0E0_WP ,130.E0_WP ,320.E0_WP ,        & ! ES |  H | He | Li | Be |  B |
            33.0E0_WP ,0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,             & !  C |  N |  O |  F | Ne |
            6.30E0_WP ,45.0E0_WP ,76.0E0_WP ,100.E0_WP ,11.0E0_WP ,             & ! Na | Mg | Al | Si |  P |
            7.70E0_WP ,1.10E0_WP ,0.00E0_WP ,3.10E0_WP ,17.0E0_WP ,             & !  S | Cl | Ar |  K | Ca |
            57.0E0_WP ,110.E0_WP ,160.E0_WP ,160.E0_WP ,120.E0_WP ,             & ! Sc | Ti |  V | Cr | Mn |
            170.E0_WP ,180.E0_WP ,180.E0_WP ,140.E0_WP ,70.0E0_WP ,             & ! Fe | Co | Ni | Cu | Zn |
            0.00E0_WP ,0.00E0_WP ,22.0E0_WP ,8.30E0_WP ,1.90E0_WP ,             & ! Ga | Ge | As | Se | Br |
            0.00E0_WP ,2.50E0_WP ,0.00E0_WP ,41.0E0_WP ,0.00E0_WP ,             & ! Kr | Rb | Sr |  Y | Zr |
            170.E0_WP ,230.E0_WP ,0.00E0_WP ,220.E0_WP ,380.E0_WP ,             & ! Nb | Mo | Tc | Ru | Rh |
            180.E0_WP ,100.E0_WP ,42.0E0_WP ,0.00E0_WP ,58.0E0_WP ,             & ! Pd | Ag | Cd | In | Sn |
            42.0E0_WP ,64.0E0_WP ,7.70E0_WP ,0.00E0_WP ,1.60E0_WP ,             & ! Sb | Te |  I | Xe | Cs |
            9.40E0_WP ,28.0E0_WP ,22.0E0_WP ,29.0E0_WP ,32.0E0_WP ,             & ! Ba | La | Ce | Pr | Nd |
            33.0E0_WP ,38.0E0_WP ,8.30E0_WP ,38.0E0_WP ,38.7E0_WP ,             & ! Pm | Sm | Eu | Gd | Tb |
            41.0E0_WP ,40.0E0_WP ,44.0E0_WP ,45.0E0_WP ,31.0E0_WP ,             & ! Dy | Ho | Er | Tm | Yb |
            48.0E0_WP ,110.E0_WP ,200.E0_WP ,310.E0_WP ,370.E0_WP ,             & ! Lu | Hf | Ta |  W | Re |
            0.00E0_WP ,320.E0_WP ,230.E0_WP ,220.E0_WP ,25.0E0_WP ,             & ! Os | Ir | Pt | Au | Hg |
            43.0E0_WP ,46.0E0_WP ,31.0E0_WP ,0.00E0_WP ,0.00E0_WP ,             & ! Tl | Pb | Bi | Po | At |
            0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,90.0E0_WP ,             & ! Rn | Fr | Ra | Ac | Th |
            0.00E0_WP ,100.E0_WP ,0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,             & ! Pa |  U | Np | Pu | Am |
            0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,             & ! Cm | Bk | Cf | Es | Fm |
            0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,0.00E0_WP ,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 10) Shear modulus --> GPa
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  SHEA_MD = (/                   & !
            0.0E0_WP,0.00E0_WP,0.00E0_WP,4.20E0_WP,132.E0_WP,0.00E0_WP,         & ! ES |  H | He | Li | Be |  B |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & !  C |  N |  O |  F | Ne |
                3.30E0_WP,17.0E0_WP,26.0E0_WP,0.00E0_WP,0.00E0_WP,              & ! Na | Mg | Al | Si |  P |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,1.30E0_WP,7.40E0_WP,              & !  S | Cl | Ar |  K | Ca |
                29.0E0_WP,44.0E0_WP,47.0E0_WP,115.E0_WP,0.00E0_WP,              & ! Sc | Ti |  V | Cr | Mn |
                82.0E0_WP,76.0E0_WP,76.0E0_WP,48.0E0_WP,43.0E0_WP,              & ! Fe | Co | Ni | Cu | Zn |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,3.70E0_WP,0.00E0_WP,              & ! Ga | Ge | As | Se | Br |
                0.00E0_WP,0.00E0_WP,6.10E0_WP,26.0E0_WP,33.0E0_WP,              & ! Kr | Rb | Sr |  Y | Zr |
                38.0E0_WP,20.0E0_WP,0.00E0_WP,173.E0_WP,155.E0_WP,              & ! Nb | Mo | Tc | Ru | Rh |
                44.0E0_WP,30.0E0_WP,19.0E0_WP,0.00E0_WP,18.0E0_WP,              & ! Pd | Ag | Cd | In | Sn |
                20.0E0_WP,16.0E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Sb | Te |  I | Xe | Cs |
                4.90E0_WP,14.0E0_WP,14.0E0_WP,15.0E0_WP,16.0E0_WP,              & ! Ba | La | Ce | Pr | Nd |
                18.0E0_WP,20.0E0_WP,7.90E0_WP,22.0E0_WP,22.0E0_WP,              & ! Pm | Sm | Eu | Gd | Tb |
                25.0E0_WP,26.0E0_WP,28.0E0_WP,31.0E0_WP,10.0E0_WP,              & ! Dy | Ho | Er | Tm | Yb |
                27.0E0_WP,30.0E0_WP,67.0E0_WP,161.E0_WP,178.E0_WP,              & ! Lu | Hf | Ta |  W | Re |
                222.E0_WP,210.E0_WP,61.0E0_WP,27.0E0_WP,0.00E0_WP,              & ! Os | Ir | Pt | Au | Hg |
                2.80E0_WP,5.60E0_WP,12.0E0_WP,0.00E0_WP,0.00E0_WP,              & ! Tl | Pb | Bi | Po | At |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,31.0E0_WP,              & ! Rn | Fr | Ra | Ac | Th |
                0.00E0_WP,111.E0_WP,0.00E0_WP,43.0E0_WP,0.00E0_WP,              & ! Pa |  U | Np | Pu | Am |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Cm | Bk | Cf | Es | Fm |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 11) Young modulus --> GPa
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  YOUN_MD = (/                   & !
           0.00E0_WP,0.00E0_WP,0.00E0_WP,4.90E0_WP,287.E0_WP,0.00E0_WP,         & ! ES |  H | He | Li | Be |  B |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & !  C |  N |  O |  F | Ne |
                10.0E0_WP,45.0E0_WP,70.0E0_WP,47.0E0_WP,0.00E0_WP,              & ! Na | Mg | Al | Si |  P |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,20.0E0_WP,              & !  S | Cl | Ar |  K | Ca |
                74.0E0_WP,116.E0_WP,128.E0_WP,279.E0_WP,198.E0_WP,              & ! Sc | Ti |  V | Cr | Mn |
                211.E0_WP,209.E0_WP,200.E0_WP,130.E0_WP,108.E0_WP,              & ! Fe | Co | Ni | Cu | Zn |
                0.00E0_WP,0.00E0_WP,8.00E0_WP,10.0E0_WP,0.00E0_WP,              & ! Ga | Ge | As | Se | Br |
                0.00E0_WP,2.40E0_WP,0.00E0_WP,64.0E0_WP,67.0E0_WP,              & ! Kr | Rb | Sr |  Y | Zr |
                105.E0_WP,329.E0_WP,0.00E0_WP,447.E0_WP,275.E0_WP,              & ! Nb | Mo | Tc | Ru | Rh |
                121.E0_WP,85.0E0_WP,50.0E0_WP,11.0E0_WP,50.0E0_WP,              & ! Pd | Ag | Cd | In | Sn |
                55.0E0_WP,43.0E0_WP,0.00E0_WP,0.00E0_WP,1.70E0_WP,              & ! Sb | Te |  I | Xe | Cs |
                13.0E0_WP,37.0E0_WP,34.0E0_WP,37.0E0_WP,41.0E0_WP,              & ! Ba | La | Ce | Pr | Nd |
                46.0E0_WP,50.0E0_WP,18.0E0_WP,55.0E0_WP,56.0E0_WP,              & ! Pm | Sm | Eu | Gd | Tb |
                61.0E0_WP,64.0E0_WP,70.0E0_WP,74.0E0_WP,24.0E0_WP,              & ! Dy | Ho | Er | Tm | Yb |
                67.0E0_WP,78.0E0_WP,186.E0_WP,411.E0_WP,463.E0_WP,              & ! Lu | Hf | Ta |  W | Re |
                0.00E0_WP,528.E0_WP,168.E0_WP,78.0E0_WP,0.00E0_WP,              & ! Os | Ir | Pt | Au | Hg |
                8.00E0_WP,16.0E0_WP,32.0E0_WP,0.00E0_WP,0.00E0_WP,              & ! Tl | Pb | Bi | Po | At |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,79.0E0_WP,              & ! Rn | Fr | Ra | Ac | Th |
                0.00E0_WP,208.E0_WP,0.00E0_WP,96.0E0_WP,0.00E0_WP,              & ! Pa |  U | Np | Pu | Am |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Cm | Bk | Cf | Es | Fm |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 12) Poisson ratio --> dimensionless
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  POIS_RT = (/                   & !
           0.00E0_WP,0.00E0_WP,0.00E0_WP,0.36E0_WP,0.03E0_WP,0.00E0_WP,         & ! ES |  H | He | Li | Be |  B |
                0.20E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & !  C |  N |  O |  F | Ne |
                0.00E0_WP,0.29E0_WP,0.35E0_WP,0.27E0_WP,0.00E0_WP,              & ! Na | Mg | Al | Si |  P |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & !  S | Cl | Ar |  K | Ca |
                0.28E0_WP,0.32E0_WP,0.37E0_WP,0.21E0_WP,0.31E0_WP,              & ! Sc | Ti |  V | Cr | Mn |
                0.29E0_WP,0.31E0_WP,0.31E0_WP,0.34E0_WP,0.25E0_WP,              & ! Fe | Co | Ni | Cu | Zn |
                0.47E0_WP,0.28E0_WP,0.00E0_WP,0.33E0_WP,0.00E0_WP,              & ! Ga | Ge | As | Se | Br |
                0.00E0_WP,0.00E0_WP,0.28E0_WP,0.24E0_WP,0.34E0_WP,              & ! Kr | Rb | Sr |  Y | Zr |
                0.40E0_WP,0.31E0_WP,0.00E0_WP,0.30E0_WP,0.26E0_WP,              & ! Nb | Mo | Tc | Ru | Rh |
                0.34E0_WP,0.37E0_WP,0.30E0_WP,0.26E0_WP,0.36E0_WP,              & ! Pd | Ag | Cd | In | Sn |
                0.00E0_WP,0.33E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Sb | Te |  I | Xe | Cs |
                0.00E0_WP,0.28E0_WP,0.24E0_WP,0.28E0_WP,0.28E0_WP,              & ! Ba | La | Ce | Pr | Nd |
                0.28E0_WP,0.27E0_WP,0.15E0_WP,0.00E0_WP,0.26E0_WP,              & ! Pm | Sm | Eu | Gd | Tb |
                0.25E0_WP,0.23E0_WP,0.24E0_WP,0.21E0_WP,0.21E0_WP,              & ! Dy | Ho | Er | Tm | Yb |
                0.26E0_WP,0.37E0_WP,0.34E0_WP,0.28E0_WP,0.30E0_WP,              & ! Lu | Hf | Ta |  W | Re |
                0.25E0_WP,0.26E0_WP,0.38E0_WP,0.44E0_WP,0.00E0_WP,              & ! Os | Ir | Pt | Au | Hg |
                0.45E0_WP,0.44E0_WP,0.33E0_WP,0.00E0_WP,0.00E0_WP,              & ! Tl | Pb | Bi | Po | At |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.27E0_WP,              & ! Rn | Fr | Ra | Ac | Th |
                0.00E0_WP,0.23E0_WP,0.00E0_WP,0.21E0_WP,0.00E0_WP,              & ! Pa |  U | Np | Pu | Am |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Cm | Bk | Cf | Es | Fm |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 13) Work function --> eV
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  WORK_FC = (/                   & !
           0.00E0_WP,0.00E0_WP,0.00E0_WP,2.93E0_WP,4.98E0_WP,4.45E0_WP,         & ! ES |  H | He | Li | Be |  B |
                5.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & !  C |  N |  O |  F | Ne |
                2.36E0_WP,3.66E0_WP,4.17E0_WP,4.79E0_WP,0.00E0_WP,              & ! Na | Mg | Al | Si |  P |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,2.29E0_WP,2.87E0_WP,              & !  S | Cl | Ar |  K | Ca |
                3.50E0_WP,4.33E0_WP,4.20E0_WP,4.50E0_WP,4.10E0_WP,              & ! Sc | Ti |  V | Cr | Mn |
                4.74E0_WP,5.00E0_WP,5.20E0_WP,4.76E0_WP,4.25E0_WP,              & ! Fe | Co | Ni | Cu | Zn |
                4.32E0_WP,5.00E0_WP,3.75E0_WP,5.90E0_WP,0.00E0_WP,              & ! Ga | Ge | As | Se | Br |
                0.00E0_WP,2.26E0_WP,2.59E0_WP,3.10E0_WP,4.05E0_WP,              & ! Kr | Rb | Sr |  Y | Zr |
                4.33E0_WP,4.57E0_WP,0.00E0_WP,4.71E0_WP,4.98E0_WP,              & ! Nb | Mo | Tc | Ru | Rh |
                5.41E0_WP,4.63E0_WP,4.08E0_WP,4.09E0_WP,4.42E0_WP,              & ! Pd | Ag | Cd | In | Sn |
                4.63E0_WP,4.95E0_WP,0.00E0_WP,0.00E0_WP,1.95E0_WP,              & ! Sb | Te |  I | Xe | Cs |
                2.52E0_WP,3.50E0_WP,2.90E0_WP,0.00E0_WP,3.20E0_WP,              & ! Ba | La | Ce | Pr | Nd |
                0.00E0_WP,2.70E0_WP,2.50E0_WP,2.90E0_WP,3.00E0_WP,              & ! Pm | Sm | Eu | Gd | Tb |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Dy | Ho | Er | Tm | Yb |
                3.30E0_WP,3.90E0_WP,4.30E0_WP,4.61E0_WP,4.72E0_WP,              & ! Lu | Hf | Ta |  W | Re |
                5.93E0_WP,5.46E0_WP,5.55E0_WP,5.38E0_WP,4.47E0_WP,              & ! Os | Ir | Pt | Au | Hg |
                3.84E0_WP,4.25E0_WP,4.34E0_WP,5.00E0_WP,0.00E0_WP,              & ! Tl | Pb | Bi | Po | At |
                0.00E0_WP,2.10E0_WP,2.80E0_WP,3.20E0_WP,3.40E0_WP,              & ! Rn | Fr | Ra | Ac | Th |
                3.70E0_WP,3.73E0_WP,3.90E0_WP,3.60E0_WP,3.70E0_WP,              & ! Pa |  U | Np | Pu | Am |
                3.90E0_WP,3.80E0_WP,4.00E0_WP,3.30E0_WP,0.00E0_WP,              & ! Cm | Bk | Cf | Es | Fm |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 14) Valence --> dimensionless
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  VALE_IN = (/                   & !
           0.00E0_WP,1.00E0_WP,0.00E0_WP,1.00E0_WP,2.00E0_WP,3.00E0_WP,         & ! ES |  H | He | Li | Be |  B |
                4.00E0_WP,5.00E0_WP,2.00E0_WP,1.00E0_WP,0.00E0_WP,              & !  C |  N |  O |  F | Ne |
                1.00E0_WP,2.00E0_WP,3.00E0_WP,4.00E0_WP,5.00E0_WP,              & ! Na | Mg | Al | Si |  P |
                6.00E0_WP,7.00E0_WP,2.00E0_WP,1.00E0_WP,2.00E0_WP,              & !  S | Cl | Ar |  K | Ca |
                3.00E0_WP,4.00E0_WP,5.00E0_WP,6.00E0_WP,7.00E0_WP,              & ! Sc | Ti |  V | Cr | Mn |
                6.00E0_WP,5.00E0_WP,4.00E0_WP,4.00E0_WP,2.00E0_WP,              & ! Fe | Co | Ni | Cu | Zn |
                3.00E0_WP,4.00E0_WP,5.00E0_WP,6.00E0_WP,7.00E0_WP,              & ! Ga | Ge | As | Se | Br |
                2.00E0_WP,1.00E0_WP,2.00E0_WP,3.00E0_WP,4.00E0_WP,              & ! Kr | Rb | Sr |  Y | Zr |
                5.00E0_WP,6.00E0_WP,7.00E0_WP,8.00E0_WP,6.00E0_WP,              & ! Nb | Mo | Tc | Ru | Rh |
                4.00E0_WP,4.00E0_WP,2.00E0_WP,3.00E0_WP,4.00E0_WP,              & ! Pd | Ag | Cd | In | Sn |
                5.00E0_WP,6.00E0_WP,7.00E0_WP,8.00E0_WP,1.00E0_WP,              & ! Sb | Te |  I | Xe | Cs |
                2.00E0_WP,3.00E0_WP,4.00E0_WP,4.00E0_WP,3.00E0_WP,              & ! Ba | La | Ce | Pr | Nd |
                3.00E0_WP,3.00E0_WP,3.00E0_WP,3.00E0_WP,4.00E0_WP,              & ! Pm | Sm | Eu | Gd | Tb |
                3.00E0_WP,3.00E0_WP,3.00E0_WP,3.00E0_WP,3.00E0_WP,              & ! Dy | Ho | Er | Tm | Yb |
                3.00E0_WP,4.00E0_WP,5.00E0_WP,6.00E0_WP,7.00E0_WP,              & ! Lu | Hf | Ta |  W | Re |
                8.00E0_WP,8.00E0_WP,6.00E0_WP,5.00E0_WP,4.00E0_WP,              & ! Os | Ir | Pt | Au | Hg |
                3.00E0_WP,4.00E0_WP,5.00E0_WP,6.00E0_WP,7.00E0_WP,              & ! Tl | Pb | Bi | Po | At |
                8.00E0_WP,1.00E0_WP,2.00E0_WP,3.00E0_WP,4.00E0_WP,              & ! Rn | Fr | Ra | Ac | Th |
                5.00E0_WP,6.00E0_WP,7.00E0_WP,8.00E0_WP,6.00E0_WP,              & ! Pa |  U | Np | Pu | Am |
                4.00E0_WP,4.00E0_WP,4.00E0_WP,3.00E0_WP,3.00E0_WP,              & ! Cm | Bk | Cf | Es | Fm |
                3.00E0_WP,3.00E0_WP,3.00E0_WP,0.00E0_WP,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 15) Noninteger bonding valence --> dimensionless
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  VALE_NI = (/                   & !
           0.00E0_WP,0.00E0_WP,0.00E0_WP,1.09E0_WP,1.99E0_WP,0.00E0_WP,         & ! ES |  H | He | Li | Be |  B |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & !  C |  N |  O |  F | Ne |
                1.11E0_WP,2.08E0_WP,2.76E0_WP,0.00E0_WP,0.00E0_WP,              & ! Na | Mg | Al | Si |  P |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,1.21E0_WP,2.22E0_WP,              & !  S | Cl | Ar |  K | Ca |
                2.85E0_WP,3.20E0_WP,3.45E0_WP,3.53E0_WP,3.41E0_WP,              & ! Sc | Ti |  V | Cr | Mn |
                3.32E0_WP,3.09E0_WP,2.83E0_WP,2.57E0_WP,2.40E0_WP,              & ! Fe | Co | Ni | Cu | Zn |
                2.43E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Ga | Ge | As | Se | Br |
                0.00E0_WP,1.24E0_WP,2.32E0_WP,3.21E0_WP,3.15E0_WP,              & ! Kr | Rb | Sr |  Y | Zr |
                4.14E0_WP,4.42E0_WP,4.24E0_WP,4.05E0_WP,3.67E0_WP,              & ! Nb | Mo | Tc | Ru | Rh |
                3.15E0_WP,2.70E0_WP,2.48E0_WP,2.51E0_WP,0.00E0_WP,              & ! Pd | Ag | Cd | In | Sn |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,1.28E0_WP,              & ! Sb | Te |  I | Xe | Cs |
                2.58E0_WP,3.50E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Ba | La | Ce | Pr | Nd |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Pm | Sm | Eu | Gd | Tb |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Dy | Ho | Er | Tm | Yb |
                0.00E0_WP,3.97E0_WP,4.51E0_WP,4.29E0_WP,4.79E0_WP,              & ! Lu | Hf | Ta |  W | Re |
                4.72E0_WP,4.36E0_WP,3.90E0_WP,3.26E0_WP,2.52E0_WP,              & ! Os | Ir | Pt | Au | Hg |
                2.38E0_WP,3.50E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Tl | Pb | Bi | Po | At |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Rn | Fr | Ra | Ac | Th |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Pa |  U | Np | Pu | Am |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Cm | Bk | Cf | Es | Fm |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 16) -log_10 (electric resisitivity) --> m Ohm
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  ELEC_RE = (/                   & !
           0.00E0_WP,0.00E0_WP,0.00E0_WP,7.02E0_WP,7.39E0_WP,-4.00E0_WP,        & ! ES |  H | He | Li | Be |  B |
                5.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & !  C |  N |  O |  F | Ne |
                7.33E0_WP,7.36E0_WP,7.58E0_WP,3.00E0_WP,7.00E0_WP,              & ! Na | Mg | Al | Si |  P |
              -15.00E0_WP,-2.00E0_WP,0.00E0_WP,7.15E0_WP,7.47E0_WP,             & !  S | Cl | Ar |  K | Ca |
                6.26E0_WP,6.39E0_WP,6.69E0_WP,6.89E0_WP,5.79E0_WP,              & ! Sc | Ti |  V | Cr | Mn |
                7.01E0_WP,7.22E0_WP,7.15E0_WP,7.77E0_WP,7.23E0_WP,              & ! Fe | Co | Ni | Cu | Zn |
                6.85E0_WP,3.30E0_WP,6.52E0_WP,0.00E0_WP,-10.00E0_WP,            & ! Ga | Ge | As | Se | Br |
                0.00E0_WP,6.92E0_WP,6.89E0_WP,6.24E0_WP,6.38E0_WP,              & ! Kr | Rb | Sr |  Y | Zr |
                6.82E0_WP,7.30E0_WP,6.69E0_WP,7.15E0_WP,7.37E0_WP,              & ! Nb | Mo | Tc | Ru | Rh |
                7.00E0_WP,7.79E0_WP,7.15E0_WP,7.09E0_WP,6.96E0_WP,              & ! Pd | Ag | Cd | In | Sn |
                6.39E0_WP,4.00E0_WP,-7.00E0_WP,0.00E0_WP,6.69E0_WP,             & ! Sb | Te |  I | Xe | Cs |
                6.46E0_WP,6.21E0_WP,6.12E0_WP,6.15E0_WP,6.19E0_WP,              & ! Ba | La | Ce | Pr | Nd |
                6.12E0_WP,6.03E0_WP,6.04E0_WP,5.89E0_WP,5.92E0_WP,              & ! Pm | Sm | Eu | Gd | Tb |
                6.04E0_WP,6.03E0_WP,6.06E0_WP,6.15E0_WP,6.55E0_WP,              & ! Dy | Ho | Er | Tm | Yb |
                6.24E0_WP,6.52E0_WP,6.88E0_WP,7.30E0_WP,6.74E0_WP,              & ! Lu | Hf | Ta |  W | Re |
                7.09E0_WP,7.33E0_WP,6.96E0_WP,7.66E0_WP,6.02E0_WP,              & ! Os | Ir | Pt | Au | Hg |
                6.82E0_WP,6.68E0_WP,5.89E0_WP,6.37E0_WP,0.00E0_WP,              & ! Tl | Pb | Bi | Po | At |
                0.00E0_WP,0.00E0_WP,6.00E0_WP,0.00E0_WP,6.77E0_WP,              & ! Rn | Fr | Ra | Ac | Th |
                6.74E0_WP,6.55E0_WP,5.92E0_WP,5.82E0_WP,0.00E0_WP,              & ! Pa |  U | Np | Pu | Am |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,              & ! Cm | Bk | Cf | Es | Fm |
                0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP,0.00E0_WP               & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 17) lattice parameter a --> Angstroem
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  ALAT_PA = (/                   & !
  00.0000E0_WP, 4.7000E0_WP, 4.2420E0_WP, 3.5100E0_WP, 2.2858E0_WP, 5.0600E0_WP,& ! ES |  H | He | Li | Be |  B |
            2.4640E0_WP, 3.8610E0_WP, 5.4030E0_WP, 5.5000E0_WP, 4.4290E0_WP,    & !  C |  N |  O |  F | Ne |
            4.2906E0_WP, 3.2094E0_WP, 4.0495E0_WP, 5.4309E0_WP,11.4500E0_WP,    & ! Na | Mg | Al | Si |  P |
           10.4370E0_WP, 6.2235E0_WP, 5.2560E0_WP, 5.3280E0_WP, 5.5884E0_WP,    & !  S | Cl | Ar |  K | Ca |
            3.3090E0_WP, 2.9508E0_WP, 3.0300E0_WP, 2.9100E0_WP, 8.9125E0_WP,    & ! Sc | Ti |  V | Cr | Mn |
            2.8665E0_WP, 2.5071E0_WP, 3.5240E0_WP, 3.6149E0_WP, 2.6649E0_WP,    & ! Fe | Co | Ni | Cu | Zn |
            4.5197E0_WP, 5.6575E0_WP, 3.7598E0_WP, 9.0540E0_WP, 6.7265E0_WP,    & ! Ga | Ge | As | Se | Br |
            5.7060E0_WP, 5.5850E0_WP, 6.0849E0_WP, 3.6474E0_WP, 3.2320E0_WP,    & ! Kr | Rb | Sr |  Y | Zr |
            3.3004E0_WP, 3.1470E0_WP, 2.7350E0_WP, 2.7059E0_WP, 3.8034E0_WP,    & ! Nb | Mo | Tc | Ru | Rh |
            3.8907E0_WP, 4.0853E0_WP, 2.9794E0_WP, 3.2523E0_WP, 5.8318E0_WP,    & ! Pd | Ag | Cd | In | Sn |
            4.3070E0_WP, 4.4572E0_WP, 7.1802E0_WP, 6.2023E0_WP, 6.1410E0_WP,    & ! Sb | Te |  I | Xe | Cs |
            5.0280E0_WP, 3.7720E0_WP, 6.1410E0_WP, 3.6725E0_WP, 3.6580E0_WP,    & ! Ba | La | Ce | Pr | Nd |
           00.0000E0_WP, 3.6210E0_WP, 4.5810E0_WP, 3.6360E0_WP, 3.6010E0_WP,    & ! Pm | Sm | Eu | Gd | Tb |
            3.5930E0_WP, 3.5773E0_WP, 3.5588E0_WP, 3.5375E0_WP, 5.4847E0_WP,    & ! Dy | Ho | Er | Tm | Yb |
            3.5031E0_WP, 3.1964E0_WP, 3.3013E0_WP, 3.1652E0_WP, 2.7610E0_WP,    & ! Lu | Hf | Ta |  W | Re |
            2.7344E0_WP, 3.8390E0_WP, 3.9242E0_WP, 4.0782E0_WP, 3.0050E0_WP,    & ! Os | Ir | Pt | Au | Hg |
            3.4566E0_WP, 4.9508E0_WP, 6.6740E0_WP, 3.3590E0_WP,00.0000E0_WP,    & ! Tl | Pb | Bi | Po | At |
           00.0000E0_WP,00.0000E0_WP, 5.1480E0_WP, 5.6700E0_WP, 5.0842E0_WP,    & ! Rn | Fr | Ra | Ac | Th |
            3.9250E0_WP, 2.8537E0_WP, 6.6630E0_WP, 6.1830E0_WP, 3.4681E0_WP,    & ! Pa |  U | Np | Pu | Am |
            3.4960E0_WP, 3.4160E0_WP, 3.3800E0_WP,00.0000E0_WP,00.0000E0_WP,    & ! Cm | Bk | Cf | Es | Fm |
           00.0000E0_WP,00.0000E0_WP,00.0000E0_WP,00.0000E0_WP,00.0000E0_WP     & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 18) lattice parameter B --> Angstroem
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  BLAT_PA = (/                   & !
  00.0000E0_WP, 4.7000E0_WP, 4.2420E0_WP, 3.5100E0_WP, 2.2858E0_WP, 5.0600E0_WP,& ! ES |  H | He | Li | Be |  B |
            2.4640E0_WP, 3.8610E0_WP, 3.4290E0_WP, 3.2800E0_WP, 4.4290E0_WP,    & !  C |  N |  O |  F | Ne |
            4.2906E0_WP, 3.2094E0_WP, 4.0495E0_WP, 5.4309E0_WP, 5.5030E0_WP,    & ! Na | Mg | Al | Si |  P |
           12.8450E0_WP, 4.4561E0_WP, 5.2560E0_WP, 5.3280E0_WP, 5.5884E0_WP,    & !  S | Cl | Ar |  K | Ca |
            3.3090E0_WP, 2.9508E0_WP, 3.0300E0_WP, 2.9100E0_WP, 8.9125E0_WP,    & ! Sc | Ti |  V | Cr | Mn |
            2.8665E0_WP, 2.5071E0_WP, 3.5240E0_WP, 3.6149E0_WP, 2.6649E0_WP,    & ! Fe | Co | Ni | Cu | Zn |
            7.6633E0_WP, 5.6575E0_WP, 3.7598E0_WP, 9.0830E0_WP, 4.6451E0_WP,    & ! Ga | Ge | As | Se | Br |
            5.7060E0_WP, 5.5850E0_WP, 6.0849E0_WP, 3.6474E0_WP, 3.2320E0_WP,    & ! Kr | Rb | Sr |  Y | Zr |
            3.3004E0_WP, 3.1470E0_WP, 2.7350E0_WP, 2.7059E0_WP, 3.8034E0_WP,    & ! Nb | Mo | Tc | Ru | Rh |
            3.8907E0_WP, 4.0853E0_WP, 2.9794E0_WP, 3.2523E0_WP, 5.8318E0_WP,    & ! Pd | Ag | Cd | In | Sn |
            4.3070E0_WP, 4.4572E0_WP, 4.7102E0_WP, 6.2023E0_WP, 6.1410E0_WP,    & ! Sb | Te |  I | Xe | Cs |
            5.0280E0_WP, 3.7720E0_WP, 6.1410E0_WP, 3.6725E0_WP, 3.6580E0_WP,    & ! Ba | La | Ce | Pr | Nd |
           00.0000E0_WP, 3.6210E0_WP, 4.5810E0_WP, 3.6360E0_WP, 3.6010E0_WP,    & ! Pm | Sm | Eu | Gd | Tb |
            3.5930E0_WP, 3.5773E0_WP, 3.5588E0_WP, 3.5375E0_WP, 5.4847E0_WP,    & ! Dy | Ho | Er | Tm | Yb |
            3.5031E0_WP, 3.1964E0_WP, 3.3013E0_WP, 3.1652E0_WP, 2.7610E0_WP,    & ! Lu | Hf | Ta |  W | Re |
            2.7344E0_WP, 3.8390E0_WP, 3.9242E0_WP, 4.0782E0_WP, 3.0050E0_WP,    & ! Os | Ir | Pt | Au | Hg |
            3.4566E0_WP, 4.9508E0_WP, 6.1170E0_WP, 3.3590E0_WP,00.0000E0_WP,    & ! Tl | Pb | Bi | Po | At |
           00.0000E0_WP,00.0000E0_WP, 5.1480E0_WP, 5.6700E0_WP, 5.0842E0_WP,    & ! Rn | Fr | Ra | Ac | Th |
            3.9250E0_WP, 5.8695E0_WP, 4.7230E0_WP, 4.8220E0_WP, 3.4681E0_WP,    & ! Pa |  U | Np | Pu | Am |
            3.4960E0_WP, 3.4160E0_WP, 3.3800E0_WP,00.0000E0_WP,00.0000E0_WP,    & ! Cm | Bk | Cf | Es | Fm |
           00.0000E0_WP,00.0000E0_WP,00.0000E0_WP,00.0000E0_WP,00.0000E0_WP     & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 19) lattice parameter C --> Angstroem
!
      REAL (WP), DIMENSION(0:105), PARAMETER ::  CLAT_PA = (/                   & !
  00.0000E0_WP, 3.4000E0_WP, 4.2420E0_WP, 3.5100E0_WP, 3.5843E0_WP, 5.0600E0_WP,& ! ES |  H | He | Li | Be |  B |
            6.7110E0_WP, 6.2650E0_WP, 5.0860E0_WP, 7.2800E0_WP, 4.4290E0_WP,    & !  C |  N |  O |  F | Ne |
            4.2906E0_WP, 5.2108E0_WP, 4.0495E0_WP, 5.4309E0_WP,11.2610E0_WP,    & ! Na | Mg | Al | Si |  P |
           24.3690E0_WP, 8.1785E0_WP, 5.2560E0_WP, 5.3280E0_WP, 5.5884E0_WP,    & !  S | Cl | Ar |  K | Ca |
            5.2733E0_WP, 4.6855E0_WP, 3.0300E0_WP, 2.9100E0_WP, 8.9125E0_WP,    & ! Sc | Ti |  V | Cr | Mn |
            2.8665E0_WP, 4.0695E0_WP, 3.5240E0_WP, 3.6149E0_WP, 4.9468E0_WP,    & ! Fe | Co | Ni | Cu | Zn |
            4.5260E0_WP, 5.6575E0_WP,10.5475E0_WP,11.6010E0_WP, 8.7023E0_WP,    & ! Ga | Ge | As | Se | Br |
            5.7060E0_WP, 5.5850E0_WP, 6.0849E0_WP, 5.7306E0_WP, 5.1470E0_WP,    & ! Kr | Rb | Sr |  Y | Zr |
            3.3004E0_WP, 3.1470E0_WP, 4.3880E0_WP, 4.2815E0_WP, 3.8034E0_WP,    & ! Nb | Mo | Tc | Ru | Rh |
            3.8907E0_WP, 4.0853E0_WP, 5.6186E0_WP, 4.9461E0_WP, 3.1819E0_WP,    & ! Pd | Ag | Cd | In | Sn |
           11.2730E0_WP, 5.9290E0_WP, 9.8103E0_WP, 6.2023E0_WP, 6.1410E0_WP,    & ! Sb | Te |  I | Xe | Cs |
            5.0280E0_WP,12.1440E0_WP, 6.1410E0_WP,11.8354E0_WP,11.7990E0_WP,    & ! Ba | La | Ce | Pr | Nd |
           00.0000E0_WP, 2.6250E0_WP, 4.5810E0_WP, 5.7826E0_WP, 5.6936E0_WP,    & ! Pm | Sm | Eu | Gd | Tb |
            5.6537E0_WP, 5.6158E0_WP, 5.5874E0_WP, 5.5546E0_WP, 5.4847E0_WP,    & ! Dy | Ho | Er | Tm | Yb |
            5.5509E0_WP, 5.0511E0_WP, 3.3013E0_WP, 3.1652E0_WP, 4.4560E0_WP,    & ! Lu | Hf | Ta |  W | Re |
            4.3173E0_WP, 3.8390E0_WP, 3.9242E0_WP, 4.0782E0_WP, 3.0050E0_WP,    & ! Os | Ir | Pt | Au | Hg |
            5.5248E0_WP, 4.9508E0_WP, 3.3040E0_WP, 3.3590E0_WP,00.0000E0_WP,    & ! Tl | Pb | Bi | Po | At |
           00.0000E0_WP,00.0000E0_WP, 5.1480E0_WP, 5.6700E0_WP, 5.0842E0_WP,    & ! Rn | Fr | Ra | Ac | Th |
            3.2380E0_WP, 4.9548E0_WP, 4.8870E0_WP,10.9630E0_WP,11.2410E0_WP,    & ! Pa |  U | Np | Pu | Am |
           11.3310E0_WP,11.0690E0_WP,11.0250E0_WP,00.0000E0_WP,00.0000E0_WP,    & ! Cm | Bk | Cf | Es | Fm |
           00.0000E0_WP,00.0000E0_WP,00.0000E0_WP,00.0000E0_WP,00.0000E0_WP     & ! Md | No | Lr | Rf | Db |
                                                           /)                     !
!
! 20) Crystal structure
!
!                           CUB : simple cubic               BCC : body-centered cubic     
!                           FCC : face-centered cubic        HEX : simple hexagonal
!                           TRG : simple trigonal            BCM : based-centered monoclinic
!                           TEP : tetrahedral packing        TRC : simple triclinic
!                           FCO : face-centered orthorombic  BOR : base orthorombic
!                           MON : simple monoclinic          CTE : centered tetragonal
!                           ORT : simple orthorombic
!
!
      CHARACTER (LEN = 3), DIMENSION(0:105), PARAMETER  ::  CRYS_ST = (/        & !
                                    '   ','HEX','FCC','BCC','HEX','TRG',        & ! ES |  H | He | Li | Be |  B |
                                          'HEX','HEX','BCM','BCM','FCC',        & !  C |  N |  O |  F | Ne |
                                          'BCC','HEX','FCC','TEP','TRC',        & ! Na | Mg | Al | Si |  P |
                                          'FCO','BOR','FCC','BCC','FCC',        & !  S | Cl | Ar |  K | Ca |
                                          'HEX','HEX','BCC','BCC','BCC',        & ! Sc | Ti |  V | Cr | Mn |
                                          'BCC','HEX','FCC','FCC','HEX',        & ! Fe | Co | Ni | Cu | Zn |
                                          'BOR','FCC','TRG','MON','BOR',        & ! Ga | Ge | As | Se | Br |
                                          'FCC','BCC','FCC','HEX','HEX',        & ! Kr | Rb | Sr |  Y | Zr |
                                          'BCC','BCC','HEX','HEX','FCC',        & ! Nb | Mo | Tc | Ru | Rh |
                                          'FCC','FCC','HEX','CTE','CTE',        & ! Pd | Ag | Cd | In | Sn |
                                          'TRG','TRG','BOR','FCC','BCC',        & ! Sb | Te |  I | Xe | Cs |
                                          'BCC','HEX','HEX','HEX','HEX',        & ! Ba | La | Ce | Pr | Nd |
                                          '   ','TRG','BCC','HEX','HEX',        & ! Pm | Sm | Eu | Gd | Tb |
                                          'HEX','HEX','HEX','HEX','FCC',        & ! Dy | Ho | Er | Tm | Yb |
                                          'HEX','HEX','BCC','BCC','HEX',        & ! Lu | Hf | Ta |  W | Re |
                                          'HEX','FCC','FCC','FCC','TRG',        & ! Os | Ir | Pt | Au | Hg |
                                          'HEX','FCC','BCM','CUB','   ',        & ! Tl | Pb | Bi | Po | At |
                                          '   ','   ','BCC','FCC','FCC',        & ! Rn | Fr | Ra | Ac | Th |
                                          'CTE','BOR','ORT','MON','HEX',        & ! Pa |  U | Np | Pu | Am |
                                          'HEX','HEX','HEX','   ','   ',        & ! Cm | Bk | Cf | Es | Fm |
                                          '   ','   ','   ','   ','   '         & ! Md | No | Lr | Rf | Db |
                                                                      /)          !
!
! 21) Magnetic type
!
!
!                           DIA : diamagnetic               PAR : paramagnetic
!                           FER : ferromagnetic             AFM : antiferromagnetic
!
!
      CHARACTER (LEN = 3), DIMENSION(0:105), PARAMETER  ::  MAGN_TY = (/        & !
                                          '   ','DIA','DIA','PAR','DIA','DIA',  & ! ES |  H | He | Li | Be |  B |
                                           'DIA','DIA','PAR','   ','DIA',       & !  C |  N |  O |  F | Ne |
                                           'PAR','PAR','PAR','DIA','DIA',       & ! Na | Mg | Al | Si |  P |
                                           'DIA','DIA','DIA','PAR','PAR',       & !  S | Cl | Ar |  K | Ca |
                                           'PAR','PAR','PAR','AFM','PAR',       & ! Sc | Ti |  V | Cr | Mn |
                                           'FER','FER','FER','DIA','DIA',       & ! Fe | Co | Ni | Cu | Zn |
                                           'DIA','DIA','DIA','DIA','DIA',       & ! Ga | Ge | As | Se | Br |
                                           'DIA','PAR','PAR','PAR','PAR',       & ! Kr | Rb | Sr |  Y | Zr |
                                           'PAR','PAR','PAR','PAR','PAR',       & ! Nb | Mo | Tc | Ru | Rh |
                                           'PAR','DIA','DIA','DIA','DIA',       & ! Pd | Ag | Cd | In | Sn |
                                           'DIA','DIA','DIA','DIA','PAR',       & ! Sb | Te |  I | Xe | Cs |
                                           'PAR','PAR','PAR','PAR','PAR',       & ! Ba | La | Ce | Pr | Nd |
                                           '   ','PAR','PAR','FER','PAR',       & ! Pm | Sm | Eu | Gd | Tb |
                                           'PAR','PAR','PAR','PAR','PAR',       & ! Dy | Ho | Er | Tm | Yb |
                                           'PAR','PAR','PAR','PAR','PAR',       & ! Lu | Hf | Ta |  W | Re |
                                           'PAR','PAR','PAR','DIA','DIA',       & ! Os | Ir | Pt | Au | Hg |
                                           'DIA','DIA','DIA','   ','   ',       & ! Tl | Pb | Bi | Po | At |
                                           '   ','   ','   ','   ','PAR',       & ! Rn | Fr | Ra | Ac | Th |
                                           'PAR','PAR','   ','PAR','PAR',       & ! Pa |  U | Np | Pu | Am |
                                           '   ','   ','   ','   ','   ',       & ! Cm | Bk | Cf | Es | Fm |
                                           '   ','   ','   ','   ','   '        & ! Md | No | Lr | Rf | Db |
                                                                      /)          !
!
END MODULE ATOMIC_PROPERTIES
