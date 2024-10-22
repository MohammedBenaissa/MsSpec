!
!=======================================================================
!
!
!
!           ************************************************************
!           * ******************************************************** *
!           * *                                                      * *
!           * *         MULTIPLE-SCATTERING SPIN-INDEPENDENT         * *
!           * *                   SPECTROSCOPY CODE                  * *
!           * *                                                      * *
!           * ******************************************************** *
!           ************************************************************
!
!
!
!
!  Written by D. Sebilleau, Groupe Theorie,
!                           Departement Materiaux-Nanosciences,
!                           Institut de Physique de Rennes,
!                           UMR CNRS-Universite 6251,
!                           Universite de Rennes-1,
!                           35042 Rennes-Cedex,
!                           France
!
!  Contributors : M. Gavaza, H. Zhao, K. Hatada, J. Xu
!
!-----------------------------------------------------------------------
!
!     As a general rule in this code, although there might be a few
!     exceptions (...), a variable whose name starts with a 'I' is a
!     switch, with a 'J' is a loop index and with a 'N' is a number.
!
!-----------------------------------------------------------------------
!
!     The main subroutines are :
!
!                * PED_MS    : computes the photoelectron diffraction
!                              cross-section
!
!                * LED_MS    : computes the low-energy electron
!                              diffraction cross-section
!
!                * XAS_MS    : computes the X-ray absorption (EXAFS/XANES)
!                              cross-section
!
!                * AED_MS    : computes the Auger electron diffraction
!                              cross-section
!
!                * RES_MS    : computes the resonant elastic scattering
!                              cross-section
!
!                * ELS_MS    : computes the electron energy loss
!                              cross-section
!
!                * E2E_MS    : computes the (e-,2e-) coincidence
!                              cross-section
!
!                * E3E_MS    : computes the (e-,3e-) coincidence
!                              cross-section
!
!                * PLS_MS    : computes the photoemission energy loss
!                              spectroscopy cross-section (photoelectron
!                              diffraction from a plasmon peak)
!
!                * EIG_MS    : computes the spectral radius of
!                              the multiple scattering matrix (eigenvalues)
!
!                * CALC_TAU  : computes the matrix elements of the
!                              scattering path operator
!
!                * BL_AMP    : computes the multiple scattering amplitude
!                              on a given atom
!
!                * FINDPATHS : generates all the multiple scattering
!                              paths the electron can follow
!
!                * PATHOP    : calculates the contribution of a given
!                              path to the scattering path operator
!
!                * MATDIF    : computes the Rehr-Albers scattering
!                              matrices
!
!
!-----------------------------------------------------------------------
!
!     Always remember, when changing the input data file, to keep the
!     format. The rule here is that the last digit of any integer or
!     character data must correspond to the tab (+) while for real data,
!     the tab precedes the point.
!
!     Do not forget, before submitting a calculation, to check the
!     consistency of the input data with the corresponding maximal
!     values in the include file.
!
!     Warning: this version is only consistent with phagen_scf.f, everything
!              related to the formerly used mufpot.f has been suppressed.
!
!-----------------------------------------------------------------------
!
!     Please report any bug or problem at :
!
!                      didier.sebilleau@univ-rennes1.fr
!
!
!
!                                            Last modified : 27 May 2021
!
!=======================================================================
!
      PROGRAM SPEC
!
!  This routine reads the various input files and calls the subroutine
!   computing the cross-section. This is the driving routine of MsSpec.
!
      USE ACCURACY_REAL
!
      USE CLUSTER
      USE CLUS_ELEC
      USE CLUSTER_LIMITS
      USE EXP_TYPE
      USE INDAT
      USE INFILES
      USE INUNITS
      USE OUTUNITS
      USE TESTS
!
      USE CHECK
      USE CROSS_SECTIONS
      USE DATA_TREATMENT
      USE INITIALIZE_CALC
      USE INPUT_DATA
      USE READ_EXT_FILES
      USE TREATMENT_PED
!
      IMPLICIT NONE
!
      INTEGER               ::  JF,J_FILE
      INTEGER               ::  ICOM,N_INFILES
      INTEGER               ::  IPHA,ITRTL,NP
!
      REAL (WP)             ::  ZEM
!
!  Reading information from the launching script procspec
!
!      READ(*,10) N_INFILES                                          !
!      READ(*,10) ICOM                                               !
!
!  List of the input data files
!
!      DO JF = 1, N_INFILES                                          !
!        READ(*,20) INDATA(JF)                                       !
!      END DO                                                        !
      ICOM      = 5                                                 !
      N_INFILES = 1                                                 !
!
      INDATA(1) = 'data/spec14.dat'                                 !
!
!..........  Loop on the data files  ..........
!
      DO J_FILE = 1, N_INFILES                                      !
!
        OPEN(UNIT=ICOM, FILE=INDATA(J_FILE), STATUS='OLD')          !
!
!  Reading the input data file INDATA(J_FILE)
!
        CALL READ_DATA(ICOM,N_INFILES,J_FILE,ITRTL)                 !
!
!..........  Reading of the TL and radial matrix elements files  ..........
!..........  (photon-electron, Coulomb and no excitation case)  ..........
!
!
!  1) Reading the T-matrix elements for the various electron beams:
!     incoming, excited (i.e. not detected), outgoing 1,
!     outgoing 2 and outgoing 3
!
        IF(INC == 1) THEN                                           !
          CALL READ_TL(1,IUI4,INFILE4,IUO1,NAT,A)                   !
        END IF                                                      !
        IF(EXC == 1) THEN                                           !
          CALL READ_TL(2,IUI7,INFILE7,IUO1,NAT,A)                   !
        END IF                                                      !
        IF(OUT1 == 1) THEN                                          !
          CALL READ_TL(3,IUI9,INFILE9,IUO1,NAT,A)                   !
        END IF                                                      !
        IF(OUT2 == 1) THEN                                          !
          CALL READ_TL(4,IUI12,INFILE12,IUO1,NAT,A)                 !
        END IF                                                      !
        IF(OUT3 == 1) THEN                                          !
          CALL READ_TL(5,IUI15,INFILE15,IUO1,NAT,A)                 !
        END IF                                                      !
!
!  2) Reading the regular (RHOR), irregular (RHOI), Auger (RHOA),
!     Coulomb delocalized-localized (RHOL) and Coulomb delocalized-delocalized
!     (RHOD) radial integrals for the various electron beams
!
!     Auger corresponds to a localized-localized Coulomb integral
!
        IF((INC == 1) .AND. (EXCITATION(1) /= 'NOINTER')) THEN      !
          CALL READ_RADIAL(1,IUI5,INFILE5,IUO1)                     !
        END IF                                                      !
        IF((EXC == 1) .AND. (EXCITATION(2) /= 'NOINTER')) THEN      !
          CALL READ_RADIAL(2,IUI8,INFILE8,IUO1)                     !
        END IF                                                      !
        IF((OUT1 == 1) .AND. (EXCITATION(3) /= 'NOINTER')) THEN     !
          CALL READ_RADIAL(3,IUI10,INFILE10,IUO1)                   !
        END IF                                                      !
        IF((OUT2 == 1) .AND. (EXCITATION(4) /= 'NOINTER')) THEN     !
          CALL READ_RADIAL(4,IUI13,INFILE13,IUO1)                   !
        END IF                                                      !
        IF((OUT3 == 1) .AND. (EXCITATION(5) /= 'NOINTER')) THEN     !
          CALL READ_RADIAL(5,IUI16,INFILE16,IUO1)                   !
        END IF                                                      !
!
!  Reading the core state wave function
!
        IF(SPECTRO == 'PLS') THEN                                   !
          CALL READ_CORE_WF(IUI3,IUO1,INFILE3)                      !
        END IF                                                      !
!
!  Reading the external cluster
!
        IF(CLU == 1) THEN                                           !
          CALL READ_CLUSTER(IUO1,IUI1,INFILE1,IPRINT,IPHA,ZEM)      !
        END IF                                                      !
!
!  Checking surface-like atoms for mean square displacements
!                    calculations
!
        CALL CHECK_VIB(NAT)                                         !
!
!  Set up of the variables used for an internal calculation
!    of the mean free path and/or of mean square displacements
!
        CALL SET_INT_CALC(NAT)                                      !
!
!  Set up of the mean square displacements
!
        CALL INIT_MSD                                               !
!
!  Checking for overlayers of empty spheres
!
        CALL CHECK_EMPTY_SPHERES(NPLAN)                             !
!
!  Caculation of the cross-section
!
        CALL CROSS_SECTION(N_INFILES,J_FILE,IPHA,ZEM,NP)            !
!
!..........  End of the loop on the data files  ..........
!
      END DO                                                        !
!
!  Treatment of the data file to prepare for plotting
!
      CALL TREAT_DATA(N_INFILES,NP)                                 !
!
!  Closing the output files
!
      IF(ISOM == 0) THEN                                            !
        CLOSE(IUO2)                                                 !
      END IF                                                        !
      IF((ISOM == 0) .AND. (N_INFILES /= 1)) THEN                   !
        CLOSE(IUO1)                                                 !
      END IF                                                        !
!
!  Read formats:
!
  10  FORMAT(I2)
  20  FORMAT(A24)
!
      END PROGRAM SPEC
