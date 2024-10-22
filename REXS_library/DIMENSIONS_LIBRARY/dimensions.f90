!
!=======================================================================
!
MODULE DIMENSION_CODE
!
!  This module contains the dimensioning of the epsilon.f90 code
!
!
      IMPLICIT NONE
!
      INTEGER               ::  NATP_M                              ! max. number of
!                                                                   ! prototypical
      PARAMETER (NATP_M = 5)                                        ! atoms
!                                                                   !
      INTEGER               ::  NATCLU_M                            ! max. number of
!                                                                   ! atoms
      PARAMETER (NATCLU_M = 14)                                     ! in the cluster
!                                                                   !
      INTEGER               ::  NAT_EQ_M                            ! max. number of
!                                                                   ! equivalent
      PARAMETER (NAT_EQ_M = 8)                                      ! atoms in one class
!                                                                   !
      INTEGER               ::  N_CL_L_M                            ! max. number of
!                                                                   ! local
      PARAMETER (N_CL_L_M = 1)                                      ! atomic classes
!                                                                   !
      INTEGER               ::  NE_M                                ! max. number of
!                                                                   ! energy
      PARAMETER (NE_M = 1)                                          ! points
!                                                                   !
      INTEGER               ::  NL_M                                ! max. number of
!                                                                   ! basis
      PARAMETER (NL_M = 12)                                         ! functions used
!                                                                   !
      INTEGER               ::  LI_M                                ! max. l value
!                                                                   ! of initial
      PARAMETER (LI_M = 3)                                          ! core state
!                                                                   !
      INTEGER               ::  NEMET_M                             ! max. number
!                                                                   ! of non equivalent
      PARAMETER (NEMET_M = 1)                                       ! absorbers
!                                                                   !
      INTEGER               ::  NO_ST_M                             ! max. value
!                                                                   ! of Rehr-Albers
      PARAMETER (NO_ST_M = 2)                                       ! cut-off
!                                                                   !
      INTEGER               ::  N_SCAT_M                            ! max. value
!                                                                   ! for the
      PARAMETER (N_SCAT_M = 10)                                     ! scattering order
!                                                                   !
      INTEGER               ::  NSO_M                               ! max. number
!                                                                   ! for spin-orbit
      PARAMETER (NSO_M = 2)                                         ! components
!                                                                   !
      INTEGER               ::  NTEMP_M                             ! max. number
!                                                                   ! of temperature
      PARAMETER (NTEMP_M = 1)                                       ! points
!                                                                   !
      INTEGER               ::  NODES_EX_M                          ! max. number
!                                                                   ! of computer
      PARAMETER (NODES_EX_M = 1)                                    ! nodes
!                                                                   !
      INTEGER               ::  NSPIN_M                             ! max. number
!                                                                   ! of spin
      PARAMETER (NSPIN_M = 1)                                       ! orientations
!                                                                   !
      INTEGER               ::  NATBL_M                             ! max. number of AL
!                                                                   ! for which MS amplitudes BL
      PARAMETER (NATBL_M = 1)                                       ! are calculated
!                                                                   !
      INTEGER               ::  N_MESH_M                            ! max. number
!                                                                   ! of mesh points
      PARAMETER (N_MESH_M = 1500)                                   ! for radial grid
!                                                                   !
      INTEGER               ::  N_TH_M                              ! max. number
!                                                                   ! of theta
      PARAMETER (N_TH_M = 2000)                                     ! values
!                                                                   !
      INTEGER               ::  N_PH_M                              ! max. number
!                                                                   ! of phi
      PARAMETER (N_PH_M = 2000)                                     ! values
!                                                                   !
      INTEGER               ::  NDIM_M                              ! max. number
!                                                                   ! of lines for the TREAT_XXX
      PARAMETER (NDIM_M = 100000)                                   ! subroutine to read
!                                                                   !
      INTEGER               ::  N_TILT_M                            ! max. number
!                                                                   ! of
      PARAMETER (N_TILT_M = 11)                                     ! tilt angles
!                                                                   !
      INTEGER               ::  N_ORD_M                             ! max. number
!                                                                   ! of iterations
      PARAMETER (N_ORD_M = 200)                                     ! for the power method
!                                                                   !
      INTEGER               ::  L_MAX           !
      INTEGER               ::  NLP_M,NLA_M,N_MU_M,N_NU_M           !
      INTEGER               ::  NATM,LINMAX,LINMAXA                 !
      INTEGER               ::  LINIMAX,LINFMAX                     !
      INTEGER               ::  NLAMBDA_M,NSPIN2_M                  !
      INTEGER               ::  NT_M,NCG_M                          !
!
      PARAMETER (L_MAX = NL_M-1)                                    !
      PARAMETER (NLP_M = NL_M, NLA_M = NL_M)                        !
      PARAMETER (N_MU_M = NO_ST_M, N_NU_M = NO_ST_M/2)              !
      PARAMETER (NATM = NATP_M+3)                                   !
      PARAMETER (LINMAX = NLP_M*NLP_M, LINMAXA = NLA_M*NLA_M)       !
      PARAMETER (LINIMAX = (LI_M+1) * (LI_M+1))                     !
      PARAMETER (LINFMAX = (LI_M+4) * (LI_M+4))                     !
      PARAMETER (NLAMBDA_M = (NO_ST_M+2) * (NO_ST_M+1)/2)           !
      PARAMETER (NSPIN2_M = 3*NSPIN_M-2)                            !
      PARAMETER (NT_M = (NL_M-1) * (1 + (NSPIN_M-1) * NL_M))        !
      PARAMETER (NCG_M = 4*LI_M +2)                                 !
!                                                                   !
      INTEGER               ::  N_BESS                              ! max. value
!                                                                   ! for arrays
      PARAMETER (N_BESS = 100*NL_M)                                 ! in Bessel routine
!                                                                   !
      INTEGER               ::  NPATH_M                             ! max. number
!                                                                   ! of paths
      PARAMETER (NPATH_M = 500)                                     ! with contribution printed
!                                                                   !
      INTEGER               ::  N_GAUNT                             ! max. number n
!                                                                   ! for n!
      PARAMETER (N_GAUNT = 5*NL_M)                                  ! in nj symbols
!                                                                   !
      INTEGER               ::  NGR_M                               ! max. truncation
!                                                                   ! number
      PARAMETER (NGR_M = 10)                                        ! for correlation expansion
!                                                                   ! number
      INTEGER               ::  NLTWO,NLMM                          !
!                                                                   !
      PARAMETER (NLTWO = 2*NL_M, NLMM = LINMAX*NGR_M)               !
!
!
!          NATP_M    : MAXIMAL NUMBER OF ATOMS IN THE UNIT CELL
!                       OF THE SUBSTRATE (OR MAXIMAL NUMBER OF
!                       PROTOTYPICAL ATOMS)
!          NATCLU_M : MAXIMAL NUMBER OF ATOMS IN THE CLUSTER
!          NAT_EQ_M : MAXIMAL NUMBER OF EQUIVALENT ATOMS IN A  CLASS
!          N_CL_L_M : MAXIMAL NUMBER OF LOCAL ATOMIC CLASSES
!          NE_M     : MAXIMAL NUMBER OF ENERGY POINTS
!          NL_M     : MAXIMAL NUMBER OF BASIS FUNCIONS USED
!                       IN THE EXPANSIONS FOR  ANY ELECTRON
!          NLP_M    : MAXIMAL NUMBER OF BASIS FUNCIONS USED
!                       IN THE EXPANSIONS FOR THE (PHOTO)ELECTRON
!          NLA_M    : MAXIMAL NUMBER OF BASIS FUNCIONS USED
!                       IN THE EXPANSIONS FOR THE AUGER ELECTRON
!          LI_M     : MAXIMAL L VALUE OF THE INITIAL CORE STATE(S)
!          NEMET_M  : MAXIMAL NUMBER OF NON EQUIVALENT ABSORBERS
!          NO_ST_M  : MAXIMAL VALUE OF THE REHR-ALBERS CUT-OFF
!                       "NO" FOR THE STORAGE OF THE SCATTERING
!                       MATRICES
!          NDIF_M   : MAXIMAL VALUE FOR THE SCATTERING ORDER
!          NSO_M    : MAXIMAL NUMBER OF SPIN-ORBIT COMPONENTS
!          NTEMP_M  : MAXIMAL NUMBER OF TEMPERATURE POINTS
!          NODES_EX_M : (MAXIMAL NUMBER OF NODES - 1) AVAILABLE
!                          FOR PARALLEL CALCULATIONS
!          N_MESH_M : MAXIMAL NUMBER OF MESH POINTS FOR THE RADIAL GRID
!          N_MU_M   : MAXIMAL NUMBER OF MU INDICES IN THE R-A EXPANSIONS
!          N_NU_M   : MAXIMAL NUMBER OF NU INDICES IN THE R-A EXPANSIONS
!          NATM     : MAXIMAL NUMBER OF ATOMS IN THE UNIT CELL
!                       (SUBSTRATE+ADSORBATES)
!          LINMAX   : MAXIMAL VALUE OF THE LINEAR INDEX (L,M)
!          LINMAXA  : MAXIMAL VALUE OF THE LINEAR INDEX (AUGER CASE)
!          LINIMAX  : MAXIMAL VALUE OF THE LINEAR INDEX (LI,MI)
!          LINFMAX  : MAXIMAL VALUE OF THE LINEAR INDEX (LF,MF)
!          NLAMBDA_M: MAXIMAL NUMBER OF VALUES OF LAMBDA STORED
!          NSPIN_M  : MAXIMAL NUMBER OF SPIN ORIENTATIONS
!                         SPIN-INDEPENDENT CASE    : NSPIN_M = 1
!                         SPIN-DEPENDENT CASE      : NSPIN_M = 2
!          NATBL_M  : MAXIMAL NUMBER OF AL FOR WHICH MS AMPLITUDES
!                       BL ARE CALCULATED
!          NSPIN2_M : NUMBER OF SPIN DEPENDENT T MATRICES :
!                         SPIN-INDEPENDENT CASE    : NSPIN2_M = 1
!                         SPIN-DEPENDENT CASE      : NSPIN2_M = 4
!          NT_M     : MAXIMAL SIZE IN (L,M) OF THE T MATRICES :
!                         SPIN-INDEPENDENT CASE : NT_M = NL_M-1
!                         SPIN-DEPENDENT CASE   : NT_M = (NL_M-1)*(NL_M+1)
!          NTH_M    : MAXIMAL NUMBER OF THETA VALUES
!          NPH_M    : MAXIMAL NUMBER OF PHI VALUES
!          NDIM_M   : MAXIMAL NUMBER OF LINES FOR THE TREAT_XXX
!                         SUBROUTINE TO READ
!          N_TILT_M : MAXIMAL NUMBER OF TILT ANGLE VALUES (GENERATED
!                         BY THE ext_dir.f CODE)
!          NCG_M    : MAXIMAL VALUE FOR THE STORAGE INDEX OF
!                         CLEBSCH-GORDON COEFFICIENT (AUGER ELECTRON ONLY)
!          N_BESS   : MAXIMAL VALUE FOR THE ARRAYS IN THE BESSEL ROUTINE
!          NPATH_M  : MAXIMAL NUMBER OF PATHS WITH CONTRIBUTION PRINTED
!          N_GAUNT  : MAXIMAL NUMBER N FOR N! IN NJ SYMBOLS
!          NGR_M    : MAXIMAL TRUNCATION ORDER FOR CORRELATION EXPANSION
!          N_ORD_M  : MAXIMAL NUMBER OF ITERATIONS FOR THE POWER METHOD
!                         USED TO APPROXIMATE THE SPECTRAL RADIUS
!
!
!
!   **********************************************************************
!   **********                                                  **********
!   **********  WARNING : ALWAYS CHECK THE CONSISTENCY BETWEEN  **********
!   **********            THE PARAMETERS IN THE INPUT FILE AND  **********
!   **********               IN THE PHASE SHIFT AND CLUSTER     **********
!   **********                  FILES AND THE DIMENSIONS        **********
!   **********                        DEFINED ABOVE             **********
!   **********                                                  **********
!   **********                NL_M      >=  LMAX+1              **********
!   **********                NE_M      >=  NE                  **********
!   **********                NO_ST_M   >=  NO                  **********
!   **********                NATP_M    >=  NAT                 **********
!   **********                NEMET_M   >=  NEMET               **********
!   **********                L_INIT_M  >=  LI                  **********
!   **********                NATCLU_M  >=  Number of atoms     **********
!   **********                NDIF_M    >=  NDIF                **********
!   **********                NPATH_M   >=  NPATH               **********
!   **********                                                  **********
!   **********************************************************************
!
!
END MODULE DIMENSION_CODE

