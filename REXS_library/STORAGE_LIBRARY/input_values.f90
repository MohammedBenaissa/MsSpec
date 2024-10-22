!
!=======================================================================
!
MODULE ABSORBER
!
!  This module contains the parameters of the absorbing atom
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  NTYPEM,NEMET
!
      REAL (WP)             ::  EMIT(3)
!
END MODULE ABSORBER
!
!=======================================================================
!
MODULE ALGORITHMS
!
!  This module contains the parameters for the MS algorithms
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  2)  ::  ALGO_IN,ALGO_EX
      CHARACTER (LEN =  2)  ::  ALGO_O1,ALGO_O2,ALGO_O3
!
END MODULE ALGORITHMS
!
!=======================================================================
!
MODULE ANA_DIR
!
!  This module contains the direction of the analyzer
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)             ::  DIRANA(3,49),ANADIR(3,49)
      REAL (WP)             ::  RTH_IN(49),RPH_IN(49)
!
END MODULE ANA_DIR
!
!=======================================================================
!
MODULE APPROX_CS
!
!  This module contains variables controlling the calculation
!    of the cross-section
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  4)  ::  CS_TYPE
!
      INTEGER               ::  ISELRULE
!
END MODULE APPROX_CS
!
!=======================================================================
!
MODULE ATOMIC_INDEX
!
!  This module contains the atomic index
!
      IMPLICIT NONE
!
      INTEGER               ::  I_AT
!
END MODULE ATOMIC_INDEX
!
!=======================================================================
!
MODULE ATOMIC_MASS
!
!  This module contains the atomic mass as a function
!            of the atom number
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATM
!
      IMPLICIT NONE
!
      REAL (WP)             ::  XMT(NATM)
!
END MODULE ATOMIC_MASS
!
!=======================================================================
!
MODULE ATOMS
!
!  This module contains the content of the cluster
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATCLU_M
!
      IMPLICIT NONE
!
!
      CHARACTER (LEN =  2)  ::  CHEM(NATCLU_M)
!
      INTEGER               ::  I_GR,I_INV,NZAT(NATCLU_M)
!
      REAL (WP)             ::  VALZ(NATCLU_M),VALZ_MAX
!
END MODULE ATOMS
!
!=======================================================================
!
MODULE ATOMS_TIP
!
!  This module contains the content of the tip cluster
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATCLU_M
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  2)  ::  CHEM_TI(NATCLU_M)
!
      INTEGER               ::  I_GR_TI,I_INV_TI,NZ_TIAT(NATCLU_M)
!
      REAL (WP)             ::  VALZ_TI(NATCLU_M),VALZ_MAX_TI
!
END MODULE ATOMS_TIP
!
!=======================================================================
!
MODULE AUGER_PARAM
!
!  This module contains the information on the Auger process
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  6)  ::  AUGER
!
      INTEGER               ::  LE_MIN,LE_MAX,L_BOUNDS(0:20,2)
!
END MODULE AUGER_PARAM
!
!=======================================================================
!
MODULE AVERAGING
!
!  This module contains the averaging parameters
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : N_TH_M,N_PH_M
!
      IMPLICIT NONE
!
      INTEGER               ::  I_SET
!
      REAL (WP)             ::  TH_0(N_TH_M),PH_0(N_PH_M)
!
END MODULE AVERAGING
!
!=======================================================================
!
MODULE BEAMS_ALGO
!
!  This module stores the algorithm
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  2)  ::  ALGO(5)
!
END MODULE BEAMS_ALGO
!
!=======================================================================
!
MODULE CALC_TYPE
!
!  This module contains the switches for the type of calculation
!
      IMPLICIT NONE
!
      INTEGER               ::  IPHI,IE,ITHETA
      INTEGER               ::  IMOD,IPOL
      INTEGER               ::  I_EXT,I_TEST
!
END MODULE CALC_TYPE
!
!=======================================================================
!
MODULE CLUS_ELEC
!
!  This module contains the switches concerning the cluster and
!    on the electrons
!
      IMPLICIT NONE
!
      INTEGER               ::  NEL,NCL,NTAB_EL(5)
      INTEGER               ::  INC,EXC,OUT1,OUT2,OUT3
      INTEGER               ::  TIP,CLU,EIG,PLA
!
END MODULE CLUS_ELEC
!
!=======================================================================
!
MODULE CLUSTER
!
!  This module contains parameters of the cluster
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  UNIT
!
      INTEGER               ::  NAT
!
      REAL (WP)             ::  A,VINT
!
END MODULE CLUSTER
!
!=======================================================================
!
MODULE CLUSTER_COORD
!
!  This module contains parameters of the cluster coordinates
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATM,NATP_M,NAT_EQ_M,NATCLU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  NATCLU,N_PROT,NATYP(NATM),NCHTYP(NATP_M)
      INTEGER               ::  NCORR(NAT_EQ_M,NATP_M),INEW_AT(NATCLU_M)
!
      REAL (WP)             ::  SYM_AT(3,NATCLU_M)
!
END MODULE CLUSTER_COORD
!
!=======================================================================
!
MODULE CLUSTER_LIMITS
!
!  This module contains parameters of the cluster limits
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATCLU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  NPLAN
!
      REAL (WP)             ::  X_MAX(NATCLU_M),X_MIN(NATCLU_M)
      REAL (WP)             ::  Y_MAX(NATCLU_M),Y_MIN(NATCLU_M)
      REAL (WP)             ::  VAL(NATCLU_M)
!
END MODULE CLUSTER_LIMITS
!
!=======================================================================
!
MODULE CLUSTER_TIP
!
!  This module contains parameters of the tip
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  UNIT_TI
!
      INTEGER               ::  NAT_TI
!
      REAL (WP)             ::  A_TI,VINT_TI
!
END MODULE CLUSTER_TIP
!
!=======================================================================
!
MODULE CONV_ACC
!
!  This module contains the parameters for the
!     convergence acceleration
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_XN,I_VA,I_GN,I_WN,LEVIN
!
      COMPLEX (WP)          ::  ALPHA,BETA
!
END MODULE CONV_ACC
!
!=======================================================================
!
MODULE CONV_TYPE
!
!  This module contains the other parameters for the
!     convergence acceleration
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  4)  ::  METHOD
!
      INTEGER               ::  I_PWM,I_ACC,N_ONE,N_MAX
      INTEGER               ::  N_ITER,N_TABLE
!
      REAL (WP)             ::  SHIFT,ACC,EXPO
!
END MODULE CONV_TYPE
!
!=======================================================================
!
MODULE DAMPING
!
!  This module contains the parameters for the damping
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               :: IATT,IATT_TI
!
END MODULE DAMPING
!
!=======================================================================
!
MODULE DEB_WAL_CLU
!
!  This module contains the parameters for the calculation
!    of the cluster Debye-Waller factors
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATM
!
      IMPLICIT NONE
!
      INTEGER               ::  IMSD,IDWSPH
!
      REAL (WP)             ::  TD,T,RSJ,UJ2(NATM)
!
END MODULE DEB_WAL_CLU
!
!=======================================================================
!
MODULE DEB_WAL_TIP
!
!  This module contains the parameters for the calculation
!    of the tip Debye-Waller factors
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATM
!
      IMPLICIT NONE
!
      INTEGER               ::  IMSD_TI,IDWSPH_TI
!
      REAL (WP)             ::  TD_TI,T_TI,RSJ_TI,UJ2_TI(NATM)
!
END MODULE DEB_WAL_TIP
!
!=======================================================================
!
MODULE DICHROISM
!
!  This module contains the parameters for dichroic calculations
!
      IMPLICIT NONE
!
      INTEGER               ::  I_DICHR
!
END MODULE DICHROISM
!
!=======================================================================
!
MODULE EIGEN
!
!  This module contains the parameters for the calculation
!    of eigenvalues
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  NE_EIG,I_VIB,I_MFP,I_DAMP
!
      REAL (WP)             ::  EIG_INI,EIG_FIN
!
END MODULE EIGEN
!
!=======================================================================
!
MODULE ELEC_TYPE
!
!  This module contains the parameters for electron type
!
      IMPLICIT NONE
!
      INTEGER               ::  NEL,NTAB_EL(5)
!
END MODULE ELEC_TYPE
!
!=======================================================================
!
MODULE EMPTY_SPHERES
!
!  This module contains the parameters for the calculation
!    of empty spheres
!
      USE DIMENSION_CODE,       ONLY : NATCLU_M
!
      IMPLICIT NONE
!
      INTEGER               ::  N_ESL,I_ES(NATCLU_M)
!
END MODULE EMPTY_SPHERES
!
!=======================================================================
!
MODULE EXP_TYPE
!
!  This module contains the switches for the type of experiment
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  SPECTRO,STEREO,MODE
      CHARACTER (LEN =  7)  ::  EXCITATION(5)
!
      INTEGER               ::  IMOD
!
END MODULE EXP_TYPE
!
!=======================================================================
!
MODULE HEADER
!
!  This module contains variables controlling the headers
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  6)  ::  AUGER1
      CHARACTER (LEN =  1)  ::  NLI,EDGE
!
      INTEGER               ::  NI,NEDGE
!
END MODULE HEADER
!
!=======================================================================
!
MODULE INDAT
!
!  This module contains the names of the input data files
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 24)  ::  INDATA(100)
!
END MODULE INDAT
!
!=======================================================================
!
MODULE INFILES
!
!  This module contains the names of the input data files
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 24)  ::  INFILE1,INFILE2,INFILE3,INFILE4
      CHARACTER (LEN = 24)  ::  INFILE5,INFILE6,INFILE7,INFILE8
      CHARACTER (LEN = 24)  ::  INFILE9,INFILE10,INFILE11,INFILE12
      CHARACTER (LEN = 24)  ::  INFILE13,INFILE14,INFILE15,INFILE16
      CHARACTER (LEN = 24)  ::  INFILE17,INFILE18,INFILE19,INFILE20
!
END MODULE INFILES
!
!=======================================================================
!
MODULE INUNITS
!
!  This module contains the units of the input data files
!
      IMPLICIT NONE
!
      INTEGER               ::  IUI1,IUI2,IUI3,IUI4
      INTEGER               ::  IUI5,IUI6,IUI7,IUI8
      INTEGER               ::  IUI9,IUI10,IUI11,IUI12
      INTEGER               ::  IUI13,IUI14,IUI15,IUI16
      INTEGER               ::  IUI17,IUI18,IUI19,IUI20
!
END MODULE INUNITS
!
!=======================================================================
!
MODULE INIT_A
!
!  This module contains parameters controlling the Auger line
!
!
      INTEGER               ::  LI_C,LI_I,LI_A
!
END MODULE INIT_A
!
!=======================================================================
!
MODULE INIT_J
!
!  This module contains the values to initialize the calculation
!    of the angular momenta j
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  S_O
!
      INTEGER               ::  JF1,JF2,I_SO
!
END MODULE INIT_J
!
!=======================================================================
!
MODULE INIT_L
!
!  This module contains the values to initialize the calculation
!    of the angular momenta l
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  5)  ::  SELRULE
!
      INTEGER               ::  LI,NNL,LF1,LF2,ISTEP_LF,I_0
!
END MODULE INIT_L
!
!=======================================================================
!
MODULE INIT_L_I
!
!  This module contains the values to initialize the calculation
!    of the angular momenta l (incoming)
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  5)  ::  SELRULE_I
!
      INTEGER               ::  LI_I,NNL1,LF1_I,LF2_I
      INTEGER               ::  ISTEP_LF_I,I_0_I
!
END MODULE INIT_L_I
!
!=======================================================================
!
MODULE INIT_L_O
!
!  This module contains the values to initialize the calculation
!    of the angular momenta l (outgoing)
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  5)  ::  SELRULE_O
!
      INTEGER               ::  LI_O,NNL2,LF1_O,LF2_O
      INTEGER               ::  ISTEP_LF_O,I_0_O
!
END MODULE INIT_L_O
!
!=======================================================================
!
MODULE INIT_M
!
!  This module contains the values to initialize the calculation
!    of the angular momenta m
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  MULTIPLET
!
      INTEGER               ::  I_SHELL,I_MULT,L_MUL,J_MUL,S_MUL
!
END MODULE INIT_M
!
!=======================================================================
!
MODULE M_RULE
!
!  This module contains the information for the m selection rules
!
      IMPLICIT NONE
!
      INTEGER               ::  M_R
!
END MODULE M_RULE
!
!=======================================================================
!
MODULE NON_BULK
!
!  This module contains the information of the non bulk planes
!
      IMPLICIT NONE
!
      INTEGER               ::  NONVOL(100)
!
END MODULE NON_BULK
!
!=======================================================================
!
MODULE ONE_CHANNEL
!
!  This module contains the type of matrix elements (one-channel case)
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  2)  ::  CHANNEL,TYPCHAN
!
END MODULE ONE_CHANNEL
!
!=======================================================================
!
MODULE ONE_CHANNEL_I
!
!  This module contains the type of matrix elements (incoming)
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  2)  ::  CHANNEL_I,TYPCHAN_I
!
END MODULE ONE_CHANNEL_I
!
!=======================================================================
!
MODULE ONE_CHANNEL_O
!
!  This module contains the type of matrix elements (outgoing)
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  2)  ::  CHANNEL_O,TYPCHAN_O
!
END MODULE ONE_CHANNEL_O
!
!=======================================================================
!
MODULE OPTICAL_ELTS
!
!  This module contains the information for the selection rules
!
      IMPLICIT NONE
!
      INTEGER               ::  I_C1,I_E1,I_E2
      INTEGER               ::  I_E3,I_M1,I_M2
!
END MODULE OPTICAL_ELTS
!
!=======================================================================
!
MODULE OPTICAL_ELTS_I
!
!  This module contains the information for the selection rules (incoming)
!
      IMPLICIT NONE
!
      INTEGER               ::  I_C1_I,I_E1_I,I_E2_I
      INTEGER               ::  I_E3_I,I_M1_I,I_M2_I
!
END MODULE OPTICAL_ELTS_I
!
!=======================================================================
!
MODULE OPTICAL_ELTS_O
!
!  This module contains the information for the selection rules (outgoing)
!
      IMPLICIT NONE
!
      INTEGER               ::  I_C1_O,I_E1_O,I_E2_O
      INTEGER               ::  I_E3_O,I_M1_O,I_M2_O
!
END MODULE OPTICAL_ELTS_O
!
!=======================================================================
!
MODULE OUTFILES
!
!  This module contains the names of the output data files
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 24)  ::  OUTFILE1,OUTFILE2,OUTFILE3,OUTFILE4
!
END MODULE OUTFILES
!
!=======================================================================
!
MODULE OUTUNITS
!
!  This module contains the units of the output data files
!
      IMPLICIT NONE
!
      INTEGER               ::  IUO1,IUO2,IUO3,IUO4
      INTEGER               ::  IUSCR,IUSCR2
!
END MODULE OUTUNITS
!
!=======================================================================
!
MODULE PLAS_CAL
!
!  This module contains parameters controlling the plasmon calculation
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  I_ANL
!
      REAL (WP)             ::  V_FL,WF_C,DISP
!
END MODULE PLAS_CAL
!
!=======================================================================
!
MODULE PLAS_EXP
!
!  This module contains parameters controlling the plasmon experiment
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  N_PLA,I_BLK
!
      REAL (WP)             ::  E_PLA(4)
!
END MODULE PLAS_EXP
!
!=======================================================================
!
MODULE REXS_EXP
!
!  This module contains parameters controlling the REXS calculation
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  5)  ::  SELRULE_IN,SELRULE_O1
!
END MODULE REXS_EXP
!
!=======================================================================
!
MODULE SELECT_RULE
!
!  This module contains the selection rule for incoming/outgoing beams
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  5)  ::  SELRULE_I,SELRULE_O
!
END MODULE SELECT_RULE
!
!=======================================================================
!
MODULE SPECTRUM
!
!  This module contains the eigenvalue spectra
!
      USE DIMENSION_CODE,       ONLY : NE_M
!
      IMPLICIT NONE
!
      INTEGER               ::  NE
      INTEGER               ::  I_SPECTRUM(NE_M)
!
END MODULE SPECTRUM
!
!=======================================================================
!
MODULE SPIN_PARAM
!
!  This module contains the spin parameters
!
      IMPLICIT NONE
!
      INTEGER               ::  I_SPIN,NSPIN,NSPIN2,ISFLIP
      INTEGER               ::  IR_DIA,NSTEP
!
END MODULE SPIN_PARAM
!
!=======================================================================
!
MODULE TEMP_LOOP
!
!  This module contains temperature loop variables
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      INTEGER               ::  ITEMP,NTEMP
!
      REAL (WP)             ::  TEMP0,TEMP1
!
END MODULE TEMP_LOOP
!
!=======================================================================
!
MODULE TEMP01
!
!  This module contains temporary variables
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  1)  ::  EDGE_C,EDGE_I,EDGE_A,MULT
!
      INTEGER               ::  IM1,IM2,NEDGE_C,NEDGE_I,NEDGE_A
!
END MODULE TEMP01
!
!=======================================================================
!
MODULE TESTS
!
!  This module contains parameters for testing purposes
!
      IMPLICIT NONE
!
      INTEGER               ::  IPRINT,ISORT1,NPATHP,ISOM
!
END MODULE TESTS
!
!=======================================================================
!
MODULE TESTS_TI
!
!  This module contains parameters for testing purposes
!
      IMPLICIT NONE
!
      INTEGER               ::  IPRINT_TI
      INTEGER               ::  NPATHP_TI
!
END MODULE TESTS_TI
!
!=======================================================================
!
MODULE TEXT_CALC
!
!  This module contains text for the description of calculations
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 48)  ::  TEXTC(10)
!
END MODULE TEXT_CALC
!
!=======================================================================
!
MODULE TMP_01
!
!  This module contains temporary storage
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      REAL (WP)             ::  THFWD_IN(NATP_M),THFWD_EX(NATP_M)
      REAL (WP)             ::  THFWD_O1(NATP_M),THFWD_O2(NATP_M)
      REAL (WP)             ::  THFWD_O3(NATP_M)
!
END MODULE TMP_01
!
!=======================================================================
!
MODULE TMP_02
!
!  This module contains temporary storage
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      REAL (WP)             ::  THBWD_IN(NATP_M),THBWD_EX(NATP_M)
      REAL (WP)             ::  THBWD_O1(NATP_M),THBWD_O2(NATP_M)
      REAL (WP)             ::  THBWD_O3(NATP_M)
!
END MODULE TMP_02
!
!=======================================================================
!
MODULE TMP_03
!
!  This module contains text for the description of calculations
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 3)  ::  UNLENGTH_EX,UNLENGTH_IN
      CHARACTER (LEN = 3)  ::  UNLENGTH_O1,UNLENGTH_O2
      CHARACTER (LEN = 3)  ::  UNLENGTH_O3
!
END MODULE TMP_03
!
!=======================================================================
!
MODULE TMP_04
!
!  This module contains parameters for testing purposes
!
      IMPLICIT NONE
!
      INTEGER               ::  ITRTL_IN,ITRTL_EX
      INTEGER               ::  ITRTL_O1,ITRTL_O2
      INTEGER               ::  ITRTL_O3
!
END MODULE TMP_04
!
!=======================================================================
!
MODULE TMP_05
!
!  This module contains text for the description of calculations
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)             ::  ALPHAR,BETAR,RACC
!
END MODULE TMP_05
!
!=======================================================================
!
MODULE XM_RHO
!
!  This module contains the atomic mass and the density for all elements
!
      USE ACCURACY_REAL
!
      IMPLICIT NONE
!
      REAL (WP)             ::  XM_AT(0:105),RHO_AT(0:105)
!
END MODULE XM_RHO
!
!=======================================================================
!
MODULE VIBR_TYPE
!
!  This module contains parameters for the calculation of vibrations
!
      USE DIMENSION_CODE,       ONLY : NATP_M
!
      IMPLICIT NONE
!
      INTEGER               ::  I_FREE(NATP_M)
!
END MODULE VIBR_TYPE










