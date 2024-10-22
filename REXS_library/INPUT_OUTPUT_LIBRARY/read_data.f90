!
!=======================================================================
!
MODULE INPUT_DATA
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE READ_DATA(ICOM,N_INFILES,JFILE,ITRTL)
!
!  This subroutine reads the input data file of the MsSpec code,
!    and either stores them into modules or pass them as arguments
!
!
!  Input parameters:
!
!       * ICOM     : Fortran unit number for the input data file
!       * N_INFILES: number of input data files to be read
!       * JFILE    : current input data file number
!
!
!  Output parameters:
!
!       * ITRTL    : switch for the number of non zero digits above which
!                      T-matrix elements will be considered to be null. by
!                      default, these elements are read with nine decimal digits
!                          ---> 1    : only tl values >= +/- 0.1 will be kept
!                          ---> 2    : only tl values >= +/- 0.01 will be kept
!                          ...
!
!
!   Author :  D. Sébilleau
!
!                                            Last modified : 28 May 2021
!
!  Modules storing the data
!
!
      USE REAL_NUMBERS
!
      USE CLUSTER
      USE CLUSTER_TIP
      USE CLUS_ELEC
      USE DEB_WAL_CLU
      USE DEB_WAL_TIP
      USE DICHROISM
      USE EIGEN
      USE EXP_TYPE
      USE INDAT
      USE INIT_J
      USE INIT_M
      USE NON_BULK
      USE OUTFILES
      USE OUTUNITS
      USE SPIN_PARAM
      USE TESTS
!
      USE AMPLI_IN
      USE AMPLI_EX
      USE AMPLI_O1
      USE AMPLI_O2
      USE AMPLI_O3
!
      USE APPROXIMATIONS_IN
      USE APPROXIMATIONS_EX
      USE APPROXIMATIONS_O1
      USE APPROXIMATIONS_O2
      USE APPROXIMATIONS_O3
!
      USE ENERGY_IN
      USE ENERGY_EX
      USE ENERGY_O1
      USE ENERGY_O2
      USE ENERGY_O3
!
      USE PHI_IN
      USE PHI_O1
      USE PHI_O2
      USE PHI_O3
!
      USE THETA_IN
      USE THETA_O1
      USE THETA_O2
      USE THETA_O3
!
      USE POLAR_IN
      USE TEST_IN
!
      USE INITIALIZE_CALC
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 40)  ::  TITLE
      CHARACTER (LEN = 24)  ::  TEXT
      CHARACTER (LEN =  7)  ::  TESLEC,DUMMY,INTERACT_I,INTERACT_O
      CHARACTER (LEN =  7)  ::  INTERACT_P,INTERACT_E
      CHARACTER (LEN =  5)  ::  SELRULE
      CHARACTER (LEN =  4)  ::  STRING,SPECNAME(15)
      CHARACTER (LEN =  4)  ::  PARTICLE,T_CAL
      CHARACTER (LEN =  2)  ::  BEAM
!
      INTEGER, INTENT(IN)   ::  ICOM,N_INFILES,JFILE
!
      INTEGER, INTENT(OUT)  ::  ITRTL
!
      INTEGER               ::  I_SPEC,I_AMP
      INTEGER               ::  IRET,IBAS
      INTEGER               ::  NE,NPHI,NTHETA
      INTEGER               ::  I_REL
      INTEGER               ::  JLINE,NLINE
      INTEGER               ::  I_DWSPH,I_SPHER
      INTEGER               ::  N_CALL,NPOINT
!
      INTEGER               ::  MOD
!
      REAL (WP)             ::  ABS,MAX
!
      DATA SPECNAME    / 'PED ','LEED','EXAF','AED ','APEC','(E,2', &
                         '(E,3','EELS','REXS','STM ','PEEL','BEEM', &
                         'EIGE','    ','    '/
!
      NLINE = 2000                                                  !
!
      IRET    = 0                                                   !
      I_SO    = 0                                                   !
      STEREO  = ' NO'                                               !
      IBAS    = 0                                                   !
!
      IDWSPH     = 0                                                !
      IDWSPH_TI  = 0                                                !
!
      I_BASIS_IN = 0                                                !
      I_BASIS_EX = 0                                                !
      I_BASIS_O1 = 0                                                !
      I_BASIS_O2 = 0                                                !
      I_BASIS_O3 = 0
!
!..........   Reading of the input data in unit ICOM   ..........
!
      READ(ICOM,1) DUMMY                                            !
      READ(ICOM,2) TITLE                                            !
      READ(ICOM,1) DUMMY                                            !
!
      DO JLINE = 1, NLINE                                           !
        READ(ICOM,49) STRING                                        !
        IF(STRING == 'TYPE') THEN                                   !
          BACKSPACE ICOM                                            !
          GO TO 600                                                 !
        END IF                                                      !
      END DO                                                        !
!
!================================================================
!                 Type of calculation
!================================================================
!
      READ(ICOM,1) DUMMY                                            !
 600  READ(ICOM,3) TEXT                                            !
      READ(ICOM,1) DUMMY                                            !
!
      READ(ICOM,13) SPECTRO,I_REL,I_SPIN,I_DICHR                    !
!
      IF( (I_DICHR == 2) .AND. (I_SPIN == 0)) THEN                  !
        PRINT 514                                                   !
        STOP                                                        !
      END IF                                                        !
!
!  Setting the spectroscopy index I_SPEC
!
      IF(SPECTRO == 'PED') THEN                                     !
        I_SPEC = 1                                                  !
      ELSE IF(SPECTRO == 'LED') THEN                                !
        I_SPEC = 2                                                  !
      ELSE IF(SPECTRO == 'XAS') THEN                                !
        I_SPEC = 3                                                  !
      ELSE IF(SPECTRO == 'AED') THEN                                !
        I_SPEC = 4                                                  !
      ELSE IF(SPECTRO == 'APC') THEN                                !
        I_SPEC = 5                                                  !
      ELSE IF(SPECTRO == 'E2E') THEN                                !
        I_SPEC = 6                                                  !
      ELSE IF(SPECTRO == 'E3E') THEN                                !
        I_SPEC = 7                                                  !
      ELSE IF(SPECTRO == 'ELS') THEN                                !
        I_SPEC = 8                                                  !
      ELSE IF(SPECTRO == 'RES') THEN                                !
        I_SPEC = 9                                                  !
      ELSE IF(SPECTRO == 'STM') THEN                                !
        I_SPEC = 10                                                 !
      ELSE IF(SPECTRO == 'PLS') THEN                                !
        I_SPEC = 11                                                 !
      ELSE IF(SPECTRO == 'BEM') THEN                                !
        I_SPEC = 12                                                 !
      ELSE IF(SPECTRO == 'EIG') THEN                                !
        I_SPEC = 13                                                 !
      ELSE                                                          !
        I_SPEC = 14                                                 !
      END IF                                                        !
!
!  Checking for the experimental parameters to read
!  The name of the spectroscopy is used as an anchor point
!
      DO JLINE = 1, NLINE                                           !
        READ(ICOM,49) STRING                                        !
        IF(STRING == SPECNAME(I_SPEC)) THEN                         !
          GO TO 602                                                 !
         END IF                                                     !
      END DO                                                        !
!
  602 BACKSPACE ICOM                                                !
      BACKSPACE ICOM                                                !
!
!================================================================
!           Experimental parameters of the spectroscopy
!================================================================
!
      CALL READ_SPECTRO_INFO(ICOM,IRET,*605,*888)                   !
!
!================================================================
!           Experimental parameters of the beams
!================================================================
!                                                                   !
!
!  Setting up the beam information for storage in modules
!
      CALL BEAM_INI(SPECTRO,MODE)                                   !
!
      IF(SPECTRO == 'PED') THEN                                     ! one photon in, one electron out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'LED') THEN                                ! one electron in, one electron out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'XAS') THEN                                ! one photon in, one photon out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'AED') THEN                                ! nothing in, one electron out
!                                                                   !
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'APC') THEN                                ! one photon in, two electrons out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O2'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'E2E') THEN                                ! one electron in, two electrons out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O2'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'E3E') THEN                                ! one electron in, three electrons out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O2'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O3'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'ELS') THEN                                ! one electron in, two electrons out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'RES') THEN                                ! one photon in, one photon out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!                                                                   !
      ELSE IF(SPECTRO == 'PLS') THEN                                ! one photon in, one electron out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL READ_BEAM_INFO(ICOM,BEAM,PARTICLE)                     !
!
      END IF                                                        !
!
!  Setting various spectroscopy-dependent parameters
!
      CALL INIT_SPECTRO(SPECTRO,SELRULE,STEREO,I_SPEC,           &  !
                        EXCITATION,INTERACT_I,INTERACT_O)           !
!
!================================================================
!           Calculation parameters
!================================================================
!
      IF(CLU == 1) THEN                                             !
        T_CAL = 'CLUS'                                              !
        CALL READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*605)    !
      END IF                                                        !
!
      IF(EIG == 1) THEN                                             !
        T_CAL = 'EIGE'                                              !
        CALL READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*605)    !
      END IF                                                        !
!
      IF(INC == 1) THEN                                             !
        T_CAL = 'INCO'                                              !
        CALL READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*605)    !
      END IF                                                        !
!
      IF(EXC == 1) THEN                                             !
        T_CAL = 'EXCI'                                              !
        CALL READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*605)    !
      END IF                                                        !
!
      IF(OUT1 == 1) THEN                                            !
        T_CAL = 'OUTG'                                              !
        CALL READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*605)    !
      END IF                                                        !
!
      IF(OUT2 == 1) THEN                                            !
        T_CAL = '2nd '                                              !
        CALL READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*605)    !
      END IF                                                        !
!
      IF(OUT3 == 1) THEN                                            !
        T_CAL = '3rd '                                              !
        CALL READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*605)    !
      END IF                                                        !
!
      IF(TIP == 1) THEN                                             !
        T_CAL = 'TIP '                                              !
        CALL READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*605)    !
      END IF                                                        !
!
!================================================================
!                 Reading the input/output files
!================================================================
!
      CALL READ_IN_OUT_FILES(ICOM,IRET,*605)                        !
!
      GO TO 601                                                     !
!
!   Missing input data file(s)
!
 605  REWIND ICOM                                                   !
      DO JLINE = 1, NLINE                                           !
        READ(ICOM,5) TESLEC                                         !
        IF(TESLEC == 'CONTROL') THEN                                !
          BACKSPACE ICOM                                            !
          READ(ICOM,34) OUTFILE1,IUO1                               !
          GO TO 601                                                 !
        END IF                                                      !
      END DO                                                        !
!
      IF(NAT > NATP_M) IRET = 2                                     !
      IF(TIP == 1) THEN                                             !
        IF(NAT_TI > NATP_M) IRET = 2                                !
      END IF                                                        !
!
      NE     = MAX(NE_IN,NE_EX,NE_O1,NE_O2,NE_O3)                   !
      NPHI   = MAX(NPHI_IN,NPHI_O1,NPHI_O2,NPHI_O3)                 !
      NTHETA = MAX(NTHETA_IN,NTHETA_O1,NTHETA_O2,NTHETA_O3)         !
!
      IF(NE > NE_M) IRET = 2                                        !
!
!    Opening of the control file for printing
!
 601  IF((JFILE == 1) .OR. (ISOM == 0)) THEN                        !
        OPEN(UNIT=IUO1, FILE=OUTFILE1, STATUS='UNKNOWN')            !
      END IF                                                        !
      IF((N_INFILES > 1) .AND. (ISOM /= 0)) THEN                    !
         WRITE(IUO1,105) INDATA(JFILE)                              !
      END IF                                                        !
!
!    Error messages
!
 888  CALL WRITE_READ_ERRORS(IRET,IUO1)                             !
!
!    Initialization of variables for test calculations
!
      CALL INIT_TEST_CALC(SPECTRO,EXCITATION,INTERACT_P,INTERACT_E) !
!
      IF((N_INFILES == 1).OR.(IBAS == 1)) ISOM = 0                  !
!
      IF((IPOL_IN == 0) .AND. (I_DICHR > 0)) THEN                   !
        PRINT 513                                                   !
        STOP                                                        !
      END IF                                                        !
!
!  Set up of the switch controlling external
!    reading of the detector directions and
!    averaging over them for an undetected electron
!
      CALL EXTERNAL_ANGLES(SPECTRO,IRET)                            !
!
!
!..........   Writing of the input data in unit IUO1   ..........
!
!
      WRITE(IUO1,100)                                               !
      WRITE(IUO1,101)                                               !
      WRITE(IUO1,101)                                               !
      WRITE(IUO1,102) TITLE                                         !
      WRITE(IUO1,101)                                               !
      WRITE(IUO1,101)                                               !
      WRITE(IUO1,203)                                               !
!
      IF(I_TEST_IN == 2) THEN                                       !
        IF(ABS(IPOL_IN) == 1) THEN                                  !
          WRITE(IUO1,525)                                           !
        ELSE IF(ABS(IPOL_IN) == 2) THEN                             !
          WRITE(IUO1,526)                                           !
        END IF                                                      !
      END IF                                                        !
!
      WRITE(IUO1,103) TEXT                                          !
      WRITE(IUO1,104) SPECTRO,I_REL,I_SPIN,I_DICHR                  !
!
!================================================================
!              Writing the experimental parameters
!================================================================
!
      CALL WRITE_SPECTRO_INFO(IUO1)                                 !
!
!================================================================
!              Writing the beam(s) parameters
!================================================================
!
      IF(SPECTRO == 'PED') THEN                                     ! one photon in, one electron out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'LED') THEN                                ! one electron in, one electron out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'XAS') THEN                                ! one photon in, one photon out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'AED') THEN                                ! nothing in, one electron out
!                                                                   !
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'APC') THEN                                ! one photon in, two electrons out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O2'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'E2E') THEN                                ! one electron in, two electrons out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O2'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'E3E') THEN                                ! one electron in, three electrons out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O2'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O3'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'ELS') THEN                                ! one electron in, two electrons out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'IN'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'RES') THEN                                ! one photon in, one photon out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!                                                                   !
      ELSE IF(SPECTRO == 'PLS') THEN                                ! one photon in, one electron out
!                                                                   !
        BEAM     = 'IN'                                             !
        PARTICLE = 'PHOT'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
        BEAM     = 'O1'                                             !
        PARTICLE = 'ELEC'                                           !
        CALL WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)                    !
!
      END IF                                                        !
!
!
!================================================================
!              Writing the calculation parameters
!================================================================
!
      CALL WRITE_CALC_PARAM(IUO1,JFILE,SPECTRO)                     !
!
      I_AMP =  MAX(I_AMP_IN,I_AMP_EX,I_AMP_O1,I_AMP_O2,I_AMP_O3)    !
!
      IF(I_AMP == 1) WRITE(6,534)                                   !
!
!  Setting up of various parameters
!
      IF(SPECTRO == 'EIG') THEN                                     !
!
!  Switch for including vibrational damping into the MS matrix
!
!           I_VIB = 0 : no vibrations included
!           I_VIB = 1 : vibrations included
!
!         and mean free path-like damping
!
!           I_MFP = 0 : no Im(k) damping included
!           I_MFP = 1 : Im(k) damping included
!
        I_VIB = MOD(I_DAMP,2)                                       !
        IF(I_VIB == 1) THEN                                         !
          IDWSPH = 1                                                !
        ELSE                                                        !
          IDWSPH = 0                                                !
        END IF                                                      !
        IF(I_DAMP <= 1) THEN                                        !
          I_MFP = 0                                                 !
        ELSE                                                        !
          I_MFP = 1                                                 !
        END IF                                                      !
      END IF                                                        !
!
      NPOINT = NPHI * NE * NTHETA                                   !
      ISORT1 = 0                                                    !
      IF(NPOINT > 250) THEN                                         !
        ISORT1 = 1                                                  !
        WRITE(IUO1,510)                                             !
      END IF                                                        !
!
      I_SPHER = MAX(I_BASIS_IN,I_BASIS_EX,I_BASIS_O1,I_BASIS_O2)    !
      I_DWSPH = MAX(IDWSPH,IDWSPH_TI)                               !
!
!  Computes the the logarithm of the Gamma function
!
      CALL CALC_LOG_GAMMA(I_SPHER)                                  !
!
!  Setting up of the initial and final values for the angular momenta
!
!
      IF(SPECTRO /= 'RES') THEN                                     !
!
!  Non REXS case
!
        N_CALL = 1                                                  !
        CALL CHOOSE_COUMAT(SPECTRO,N_CALL)                          !
!
      ELSE                                                          !
!
!  REXS incoming beam
!
        N_CALL = 1                                                  !
        CALL CHOOSE_COUMAT(SPECTRO,N_CALL)                          !
!
!  REXS outgoing beam
!
        N_CALL = 2                                                  !
        CALL CHOOSE_COUMAT(SPECTRO,N_CALL)                          !
!
      END IF                                                        !
!
!  Initialization of the spin-orbit parameters
!
      CALL INIT_SPIN_ORBIT(I_SPIN,IRET,*888)                        !
!
!  Initialization of multiplet variables (Clebsch-Gordan ...)
!
      IF(I_MULT == 1) THEN                                          !
        CALL INIT_MULTIPLET(SPECTRO)                                !
      END IF                                                        !
!
!..........  Setting up various parameters  ..........
!..........        for consistency          ..........
!
      CALL INIT_PARAM_CONS(SPECTRO,ITRTL)                           !
!
!..........  Check of the dimensioning in the Gaussian case  ..........
!
      CALL STOP_EXT(SPECTRO)                                        !
!
!....................   Read FORMATs   ....................
!
!
   1  FORMAT(A7)
   2  FORMAT(21X,A40)
   3  FORMAT(27X,A24)
   5  FORMAT(49X,A7)
  13  FORMAT(7X,A3,9X,I1,9X,I1,9X,I1)
  34  FORMAT(9X,A24,5X,I2)
  49  FORMAT(27X,A4)
!
!....................   Write FORMATs   ....................
!
!
 100  FORMAT(//////////,'******************************',               &
             '****************************************************')
 101  FORMAT('*********************',40X,'*********************')
 102  FORMAT('*********************',A40,'*********************')
 103  FORMAT(21X,A24,/)
 104  FORMAT(10X,A3,9X,I1,9X,I1,9X,I1,9X,'SPECTRO,I_REL,I_SPIN,',       &
             'I_DICHR')
 105  FORMAT(///,'ooooooooooooooooooooooooooooooooooooooooo',           &
             'ooooooooooooooooooooooooooooooooooooooooo',/,             &
             'oooooooooooooooo',50X,'oooooooooooooooo',/,               &
             'oooooooooooooooo   INPUT DATA FILE   :  ',A24,            &
             '  oooooooooooooooo',/,'oooooooooooooooo',50X,             &
             'oooooooooooooooo',/,'oooooooooooooooooooooooooooo',       &
             'ooooooooooooooooooooooooooooooooooooooooooooooooo',       &
             'ooooo',///)
!
 203  FORMAT('**************************************************',      &
             '********************************',//////////)
 510  FORMAT(///,4X,' <<<<<<<<<<  AS THE CALCULATION HAS MORE THAN ',   &
             '250 POINTS, SOME OUTPUTS HAVE BEEN SUPRESSED  ',          &
             '>>>>>>>>>>',                                              &
             ///)
 513  FORMAT(///,15X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE IPOL = 0 AND ',   &
             'I_DICHR > 0  >>>>>>>>>>')
 514  FORMAT(///,15X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE I_DICHR = 2 ',    &
             'AND I_SPIN = 0  >>>>>>>>>>')
 525  FORMAT(///,14X,'ATOMIC CALCULATION : Z AXIS ALONG POLARIZATION ', &
                'DIRECTION',/,'  ',/,'   ',/,'  ')
 526  FORMAT(///,18X,'ATOMIC CALCULATION : Z AXIS ALONG LIGHT ',        &
                'DIRECTION',/,'  ',/,'   ',/,'  ')
 534  FORMAT(//,20X,'THIS CALCULATION OUTPUTS ALSO THE AMPLITUDES')
!
      END SUBROUTINE READ_DATA
!
!=======================================================================
!
      SUBROUTINE READ_SPECTRO_INFO(ICOM,IRET,*,*)
!
!  This subroutine reads the information about the selected spectroscopy
!    and stores it in COMMON blocks for further use
!
!
!  Input variables :
!
!
!                       ICOM      :  Fortran index of the input data file
!
!
!  Output variables :
!
!
!                       IRET      :  error return code
!                       *         :  label return for error in input data file
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified :  5 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO
!
      USE APPROX_CS
      USE EXP_TYPE
      USE HEADER
      USE INIT_A
      USE INIT_J
      USE INIT_L
      USE INIT_M
      USE PLAS_EXP
      USE REXS_EXP
      USE SPIN_PARAM
      USE TEMP01
!
      USE INITIALIZE_CALC
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  7)  ::  DUMMY
      CHARACTER (LEN =  6)  ::  AUGER
!
      INTEGER, INTENT(IN)   ::  ICOM
      INTEGER, INTENT(OUT)  ::  IRET
!
      INTEGER               ::  J,J_PLA
!
!  Skipping the three title lines
!
      READ(ICOM,1) DUMMY                                            !
      READ(ICOM,1) DUMMY                                            !
      READ(ICOM,1) DUMMY                                            !
!
!  Photoelectron diffraction case:
!
      IF(SPECTRO == 'PED') THEN                                     !
!
        READ(ICOM,10) NI,NLI,S_O,SELRULE,I_SO                       !
        READ(ICOM,11) IMOD                                          !
!
!........Getting core-state related parameters (LI, ...)
!
        CALL INIT_CORE_STATE(NLI,I_SO,S_O,LI,IRET,*900)             !
!
      ELSE IF(SPECTRO == 'LED') THEN                                !
!
        READ(ICOM,11) IMOD                                          !
        SELRULE = '00000'                                           !
!
      ELSE IF(SPECTRO == 'XAS') THEN                                !
!
        READ(ICOM,12) EDGE,NEDGE,SELRULE                            !
        READ(ICOM,11) IMOD                                          !
!
        LI = NEDGE / 2                                              !
        IF(NEDGE > 1) I_SO = 1                                      !
!
      ELSE IF(SPECTRO == 'AED') THEN                                !
!
        READ(ICOM,13) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A  !
        READ(ICOM,14) I_MULT,IM1,MULT,IM2                           !
        READ(ICOM,11) IMOD                                          !
!
        LI_C = NEDGE_C / 2                                          !
        LI_I = NEDGE_I / 2                                          !
        LI_A = NEDGE_A / 2                                          !
!
        IF((EDGE_I == EDGE_A) .AND. (LI_I == LI_A)) THEN            !
          I_SHELL = 1                                               !
        ELSE                                                        !
          I_SHELL = 0                                               !
        END IF                                                      !
!
        IF(EDGE_C == 'K') THEN                                      !
          AUGER = ' '//EDGE_C//EDGE_I//CHAR(48+NEDGE_I)//    &      !
                  EDGE_A//CHAR(48+NEDGE_A)                          !
        ELSE                                                        !
          AUGER = EDGE_C//CHAR(48+NEDGE_C)//EDGE_I//         &      !
                  CHAR(48+NEDGE_I)//EDGE_A//CHAR(48+NEDGE_A)        !
        END IF                                                      !
        AUGER1 = AUGER                                              !
!
      ELSE IF(SPECTRO == 'APC') THEN                                !
!
        READ(ICOM,10) NI,NLI,S_O,SELRULE,I_SO                       !
        READ(ICOM,13) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A  !
        READ(ICOM,14) I_MULT,IM1,MULT,IM2                           !
        READ(ICOM,11) IMOD                                          !
!
      ELSE IF(SPECTRO == 'E2E') THEN                                !
!
        READ(ICOM,18) NI,NLI,S_O,ISELRULE,I_SO                      !
        READ(ICOM,15) MODE,CS_TYPE                                  !
        READ(ICOM,11) IMOD                                          !
!
!........Getting core-state related parameters (LI, ...)
!
        CALL INIT_CORE_STATE(NLI,I_SO,S_O,LI,IRET,*900)             !
!
      ELSE IF(SPECTRO == 'E3E') THEN                                !
!
        READ(ICOM,18) NI,NLI,S_O,ISELRULE,I_SO                      !
        READ(ICOM,15) MODE,CS_TYPE                                  !
        READ(ICOM,11) IMOD                                          !
!
!........Getting core-state related parameters (LI, ...)
!
        CALL INIT_CORE_STATE(NLI,I_SO,S_O,LI,IRET,*900)             !
!
      ELSE IF(SPECTRO == 'ELS') THEN                                !
!
        READ(ICOM,18) NI,NLI,S_O,ISELRULE,I_SO                      !
        READ(ICOM,15) MODE,CS_TYPE                                  !
        READ(ICOM,11) IMOD                                          !
!
!........Getting core-state related parameters (LI, ...)
!
        CALL INIT_CORE_STATE(NLI,I_SO,S_O,LI,IRET,*900)             !
!
      ELSE IF(SPECTRO == 'RES') THEN                                !
!
        READ(ICOM,12) EDGE,NEDGE,SELRULE_IN,SELRULE_O1              !
        READ(ICOM,11) IMOD                                          !
!
      ELSE IF(SPECTRO == 'PLS') THEN                                !
!
        READ(ICOM,10) NI,NLI,S_O,SELRULE,I_SO                       !
        READ(ICOM,16) N_PLA,I_BLK                                   !
        READ(ICOM,17) (E_PLA(J_PLA), J_PLA=1,4)                     !
        READ(ICOM,11) IMOD                                          !
!
      ELSE IF(SPECTRO == 'STM') THEN                                !
!
        CONTINUE                                                    !
!
      ELSE IF(SPECTRO == 'BEM') THEN                                !
!
        CONTINUE                                                    !
!
      END IF                                                        !
!
!  Initialization of the spin-orbit parameters
!
      CALL INIT_SPIN_ORBIT(I_SPIN,IRET,*901)                        !
!
!  Initialization of multiplet variables (Clebsch-Gordan ...)
!
      IF(I_MULT == 1) THEN                                          !
        CALL INIT_MULTIPLET(SPECTRO)                                !
      END IF                                                        !
!
      GO TO 800
!
 800  RETURN                                                        !
 900  RETURN 1                                                      !
 901  RETURN 2                                                      !
!
!  Formats:
!
   1  FORMAT(A7)
!
  10  FORMAT(8X,I1,A1,8X,A3,4X,A5,8X,I2)
  11  FORMAT(9X,I1)
  12  FORMAT(8X,A1,I1,4X,A5)
  13  FORMAT(8X,A1,I1,8X,A1,I1,8X,A1,I1)
  14  FORMAT(9X,I1,8X,I1,A1,I1)
  15  FORMAT(7X,A3,6X,A4)
  16  FORMAT(9X,I1,9X,I1)
  17  FORMAT(8X,F5.2,5X,F5.2,5X,F5.2,5X,F5.2)
  18  FORMAT(8X,I1,A1,8X,A3,8X,I2,8X,I2)
!
      END SUBROUTINE READ_SPECTRO_INFO
!
!=======================================================================
!
      SUBROUTINE READ_BEAM_INFO(ICOM,BEAM,PARTICLE)
!
!  This subroutine reads the information about any beam (electron/photon)
!    and stores it in COMMON blocks for further use
!
!
!  Input variables :
!
!
!                       ICOM      :  Fortran index of the input data file
!                       BEAM      :  type of beam used
!                                            --> 'IN'   : incoming
!                                            --> 'EX'   : excited
!                                            --> 'O1'   : 1st outgoing
!                                            --> 'O2'   : 2nd outgoing
!                                            --> 'O3'   : 3rd outgoing
!                       PARTICLE  :  type of particles in the beam
!                                            --> 'ELEC' : electrons
!                                            --> 'PHOT' : photons
!
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 27 May 2021
!
!
      USE REAL_NUMBERS,       ONLY : ZERO
!
      USE AMPLI_IN
      USE AVER_IN
      USE ENERGY_IN
      USE INTEGR_IN
      USE PHI_IN
      USE POLAR_IN
      USE THETA_IN
!
      USE AMPLI_EX
      USE ENERGY_EX
!
      USE AMPLI_O1
      USE AVER_O1
      USE ENERGY_O1
      USE INTEGR_O1
      USE PHI_O1
      USE POLAR_O1
      USE THETA_O1
!
      USE AMPLI_O2
      USE AVER_O2
      USE ENERGY_O2
      USE INTEGR_O2
      USE PHI_O2
      USE POLAR_O2
      USE THETA_O2
!
      USE AMPLI_O3
      USE AVER_O3
      USE ENERGY_O3
      USE INTEGR_O3
      USE PHI_O3
      USE POLAR_O3
      USE THETA_O3
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  4)  ::  PARTICLE
      CHARACTER (LEN =  2)  ::  BEAM
!
      CHARACTER (LEN =  4)  ::  STRING1,STRING2
      CHARACTER (LEN =  4)  ::  B_READ(5),DUMMY
!
      INTEGER, INTENT(IN)   ::  ICOM
!
      INTEGER               ::  JLINE,NLINE
      INTEGER               ::  I_BEAM
      INTEGER               ::  IPHI,ITHETA,IE
      INTEGER               ::  IPOL,I_INT,I_AMP
      INTEGER               ::  I_AVER,ICHKDIR
      INTEGER               ::  NPHI,NTHETA,NE
!
      REAL (WP)             ::  PHI0,PHI1
      REAL (WP)             ::  THETA0,THETA1
      REAL (WP)             ::  E0,E1,ACCEPT
!
      DATA   B_READ      /'INCO','EXCI','OUTG','2nd ','3rd '/
!
      NLINE = 100                                                  !
!
!  Initialization of I_BEAM
!
      IF(BEAM == 'IN') THEN                                         !
        I_BEAM = 1                                                  !
      ELSE IF(BEAM == 'EX') THEN                                    !
        I_BEAM = 2                                                  !
      ELSE IF(BEAM == 'O1') THEN                                    !
        I_BEAM = 3                                                  !
      ELSE IF(BEAM == 'O2') THEN                                    !
        I_BEAM = 4                                                  !
      ELSE IF(BEAM == 'O3') THEN                                    !
        I_BEAM = 5                                                  !
      END IF                                                        !
!
!  Going to the selected beam
!
      DO JLINE = 1, NLINE                                           !
        READ(ICOM,10) STRING1,STRING2                               !
        IF(STRING1 == B_READ(I_BEAM)) THEN                          !
          IF(STRING2 == PARTICLE) THEN                              !
            READ(ICOM,10) DUMMY                                     !
            GO TO 100                                               !
          END IF                                                    !
        END IF                                                      !
      END DO                                                        !
!
!  Electron case:
!
 100  IF(PARTICLE == 'ELEC') THEN                                   !
!
        IF(I_BEAM /= 2) READ(ICOM,11) IPHI,NPHI,PHI0,PHI1           !
        IF(I_BEAM /= 2) READ(ICOM,11) ITHETA,NTHETA,THETA0,THETA1   !
        READ(ICOM,11) IE,NE,E0,E1                                   !
        IF(I_BEAM /= 2) READ(ICOM,11) I_INT,I_AMP                   !
        IF(I_BEAM /= 2) READ(ICOM,12) I_AVER,ACCEPT,ICHKDIR         !
        IF(E0 == ZERO) E0 = 0.0001E0_WP                             !
        GO TO 200                                                   !
!
!  Photon case:
!
      ELSE IF(PARTICLE == 'PHOT') THEN                              !
!
        READ(ICOM,11) IPHI,NPHI,PHI0,PHI1                           !
        READ(ICOM,11) ITHETA,NTHETA,THETA0,THETA1                   !
        READ(ICOM,11) IE,NE,E0,E1                                   !
        READ(ICOM,11) IPOL                                          !
        IF(E0 == ZERO) E0 = 0.0001E0_WP                             !
        GO TO 300                                                   !
!
      END IF                                                        !
!
!  Electron case
!
 200  IF(BEAM == 'IN') THEN                                         !
!
        IPHI_IN    = IPHI                                           !
        ITHETA_IN  = ITHETA                                         !
        IE_IN      = IE                                             !
        I_INT_IN   = I_INT                                          !
        I_AMP_IN   = I_AMP                                          !
        I_AVER_IN  = I_AVER                                         !
        ICHKDIR_IN = ICHKDIR                                        !
!
        NPHI_IN    = NPHI                                           !
        NTHETA_IN  = NTHETA                                         !
        NE_IN      = NE                                             !
!
        PHI0_IN    = PHI0                                           !
        PHI1_IN    = PHI1                                           !
        THETA0_IN  = THETA0                                         !
        THETA1_IN  = THETA1                                         !
        E0_IN      = E0                                             !
        E1_IN      = E1                                             !
        ACCEPT_IN  = ACCEPT                                         !
!
      ELSE IF(BEAM == 'EX') THEN                                    !
!
        I_AMP_EX   = I_AMP                                          !
        IE_EX      = IE                                             !
!
        NE_EX      = NE                                             !
!
        E0_EX      = E0                                             !
        E1_EX      = E1                                             !
!
      ELSE IF(BEAM == 'O1') THEN                                    !
!
        IPHI_O1    = IPHI                                           !
        ITHETA_O1  = ITHETA                                         !
        IE_O1      = IE                                             !
        I_INT_O1   = I_INT                                          !
        I_AMP_O1   = I_AMP                                          !
        I_AVER_O1  = I_AVER                                         !
        ICHKDIR_O1 = ICHKDIR                                        !
!
        NPHI_O1    = NPHI                                           !
        NTHETA_O1  = NTHETA                                         !
        NE_O1      = NE                                             !
!
        PHI0_O1    = PHI0                                           !
        PHI1_O1    = PHI1                                           !
        THETA0_O1  = THETA0                                         !
        THETA1_O1  = THETA1                                         !
        E0_O1      = E0                                             !
        E1_O1      = E1                                             !
        ACCEPT_O1  = ACCEPT                                         !
!
      ELSE IF(BEAM == 'O2') THEN                                    !
!
        IPHI_O2    = IPHI                                           !
        ITHETA_O2  = ITHETA                                         !
        IE_O2      = IE                                             !
        I_INT_O2   = I_INT                                          !
        I_AMP_O2   = I_AMP                                          !
        I_AVER_O2  = I_AVER                                         !
        ICHKDIR_O2 = ICHKDIR                                        !
!
        NPHI_O2    = NPHI                                           !
        NTHETA_O2  = NTHETA                                         !
        NE_O2      = NE                                             !
!
        PHI0_O2    = PHI0                                           !
        PHI1_O2    = PHI1                                           !
        THETA0_O2  = THETA0                                         !
        THETA1_O2  = THETA1                                         !
        E0_O2      = E0                                             !
        E1_O2      = E1                                             !
        ACCEPT_O2  = ACCEPT                                         !
!
      ELSE IF(BEAM == 'O3') THEN                                    !
!
        IPHI_O3    = IPHI                                           !
        ITHETA_O3  = ITHETA                                         !
        IE_O3      = IE                                             !
        I_INT_O3   = I_INT                                          !
        I_AMP_O3   = I_AMP                                          !
        I_AVER_O3  = I_AVER                                         !
        ICHKDIR_O3 = ICHKDIR                                        !
!
        NPHI_O3    = NPHI                                           !
        NTHETA_O3  = NTHETA                                         !
        NE_O3      = NE                                             !
!
        PHI0_O3    = PHI0                                           !
        PHI1_O3    = PHI1                                           !
        THETA0_O3  = THETA0                                         !
        THETA1_O3  = THETA1                                         !
        E0_O3      = E0                                             !
        E1_O3      = E1                                             !
        ACCEPT_O3  = ACCEPT                                         !
!
      END IF                                                        !
      GO TO 400                                                     !
!
!  Photon case
!
 300  IF(BEAM == 'IN') THEN                                         !
!
        IPHI_IN    = IPHI                                           !
        ITHETA_IN  = ITHETA                                         !
        IE_IN      = IE                                             !
        IPOL_IN    = IPOL                                           !
!
        NPHI_IN    = NPHI                                           !
        NTHETA_IN  = NTHETA                                         !
        NE_IN      = NE                                             !
!
        PHI0_IN    = PHI0                                           !
        PHI1_IN    = PHI1                                           !
        THETA0_IN  = THETA0                                         !
        THETA1_IN  = THETA1                                         !
        E0_IN      = E0                                             !
        E1_IN      = E1                                             !
!
      ELSE IF(BEAM == 'O1') THEN                                    !
!
        IPHI_O1    = IPHI                                           !
        ITHETA_O1  = ITHETA                                         !
        IE_O1      = IE                                             !
        IPOL_O1    = IPOL                                           !
!
        NPHI_O1    = NPHI                                           !
        NTHETA_O1  = NTHETA                                         !
        NE_O1      = NE                                             !
!
        PHI0_O1    = PHI0                                           !
        PHI1_O1    = PHI1                                           !
        THETA0_O1  = THETA0                                         !
        THETA1_O1  = THETA1                                         !
        E0_O1      = E0                                             !
        E1_O1      = E1                                             !
!
      END IF                                                        !
!
 400  REWIND ICOM                                                   !
!
!  Formats
!
  10  FORMAT(48X,A4,9X,A4)
  11  FORMAT(8X,I2,6X,I4,6X,F7.2,3X,F7.2)
  12  FORMAT(9X,I1,8X,F5.2,6X,I1)
!
      END SUBROUTINE READ_BEAM_INFO
!
!=======================================================================
!
      SUBROUTINE READ_CALC_PARAM(ICOM,JFILE,T_CAL,SPECTRO,IRET,*)
!
!  This subroutine reads the calculation parameters and stores them
!                       for further use
!
!
!  Input variables :
!
!
!                       ICOM      :  Fortran index of the input data file
!                       JFILE     :  index of the input data file
!                       T_CAL     :  type of calculation block
!                                            --> 'CLUS'   : cluster
!                                            --> 'EIGE'   : eigenvalue
!                                            --> 'INCO'   : incoming electron
!                                            --> 'EXCI'   : excited electron
!                                            --> 'OUTG'   : 1st outgoing electron
!                                            --> '2nd '   : 2nd outgoing electron
!                                            --> '3rd '   : 3rd outgoing electron
!                                            --> 'TIP '   : tip
!                       SPECTRO   :  spectroscopy condidered
!                                            --> 'PED'   : photoelectron diffraction
!                                            --> 'LED'   : low-energy electron diffraction
!                                            --> 'XAS'   : X-ray absorption spectroscopy
!                                            --> 'AED'   : Auger electron diffraction
!                                            --> 'APC'   : Auger photoelectron coincidence spectroscopy
!                                            --> 'E2E'   : (e,2e) spectroscopy
!                                            --> 'E3E'   : (e,3e) spectroscopy
!                                            --> 'ELS'   : electron energy loss spectroscopy
!                                            --> 'RES'   : resonant elastic X-ray scattering
!                                            --> 'STM'   : scanning tunneling microscopy
!                                            --> 'PLS'   : photoelectron energy loss spectroscopy
!                                            --> 'BEM'   : ballistic energy electron microscopy
!                                            --> 'EIG'   : eigenvalue calculation
!
!
!  Output variables :
!
!
!                       IRET      :  error code on exit
!                       *         :  return with error code IRET
!
!
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 14 Jun 2021
!
!
      USE REAL_NUMBERS,    ONLY : FOUR
      USE CONSTANTS_P1,    ONLY : BOHR
!
      USE ALGORITHMS
!
      USE APPROXIMATIONS_IN
      USE APPROXIMATIONS_EX
      USE APPROXIMATIONS_O1
      USE APPROXIMATIONS_O2
      USE APPROXIMATIONS_O3
!
      USE F_TH_IN
      USE F_TH_EX
      USE F_TH_O1
      USE F_TH_O2
      USE F_TH_O3
!
      USE LINLBD_IN
      USE LINLBD_EX
      USE LINLBD_O1
      USE LINLBD_O2
      USE LINLBD_O3
!
      USE MFP_IN
      USE MFP_EX
      USE MFP_O1
      USE MFP_O2
      USE MFP_O3
!
      USE RA_IN
      USE RA_EX
      USE RA_O1
      USE RA_O2
      USE RA_O3
!
      USE RENORM_IN
      USE RENORM_EX
      USE RENORM_O1
      USE RENORM_O2
      USE RENORM_O3
!
      USE SPIN_IN
      USE SPIN_EX
      USE SPIN_O1
      USE SPIN_O2
      USE SPIN_O3
!
      USE TEST_IN
      USE TEST_EX
      USE TEST_O1
      USE TEST_O2
      USE TEST_O3
!
      USE ATOMS
      USE ATOMS_TIP
      USE CLUSTER
      USE CLUSTER_TIP
      USE CONV_ACC
      USE CONV_TYPE
      USE DAMPING
      USE DEB_WAL_CLU
      USE DEB_WAL_TIP
      USE EIGEN
      USE INIT_L
      USE NON_BULK
      USE PLAS_CAL
      USE SPECTRUM
      USE TEMP_LOOP
      USE TESTS
      USE TESTS_TI
      USE TEXT_CALC
      USE TMP_01
      USE TMP_02
      USE TMP_03
      USE TMP_04
      USE TMP_05
!
      USE INITIALIZE_CALC
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 48)  ::  TITLE
      CHARACTER (LEN =  7)  ::  TESLEC
      CHARACTER (LEN =  4)  ::  T_CAL
      CHARACTER (LEN =  4)  ::  DUMMY,STRING1,STRING2
      CHARACTER (LEN =  3)  ::  SPECTRO
!
      INTEGER, INTENT(IN)   ::  ICOM,JFILE
      INTEGER, INTENT(OUT)  ::  IRET
!
      INTEGER               ::  JLINE,NLINE
      INTEGER               ::  I,J,K
      INTEGER               ::  NLEC,NDEB,NFIN
      INTEGER               ::  N_LAST,N_LINE_E,NAT_OLD
      INTEGER               ::  JAT
!
      INTEGER               ::  INT
!
      REAL (WP)             ::  SMALL
!
      REAL (WP)             ::  MIN,MAX,FLOAT,ABS
!
      COMPLEX (WP)          ::  CMPLX
!
      DATA SMALL   / 0.0001E0_WP /
!
      NLINE = 2000                                                  !
!                                                          !
!  Going to the selected calculation block                 ! STRING1 must be equal to 'CALC'
!                                                          ! and STRING2 to T_CAL
      DO JLINE = 1, NLINE                                           !
        READ(ICOM,5) STRING1,STRING2                                !
        IF(STRING1 == 'CALC') THEN                                  !
          IF(STRING2 == T_CAL) THEN                                 !
            BACKSPACE ICOM                                          !
            READ(ICOM,7) TITLE                                      !
            READ(ICOM,5) DUMMY                                      !
            GO TO 100                                               !
          END IF                                                    !
        END IF                                                      !
      END DO                                                        !
!
!================================================================
!              Reading the calculation parameters
!                           (cluster)
!================================================================
!
  100 IF(T_CAL == 'CLUS') THEN                                      !
!
        READ(ICOM,10) NAT,A,UNIT,I_GR                               !
!
        IF(UNIT == 'ATU') THEN                                      !
          A = A * BOHR                                              !
        END IF                                                      !
!
        READ(ICOM,11) ISOM,NONVOL(JFILE),NPATHP,VINT                !
        READ(ICOM,12) IDWSPH,IATT,IPRINT                            !
!
        READ(ICOM,13) IMSD,TD,T,RSJ                                 !
!
        NLEC = INT((NAT - SMALL) / 4) + 1                           !
!
        READ(ICOM,14) ITEMP,NTEMP,TEMP0,TEMP1                       !
!
        DO I = 1, NLEC                                              !
          NDEB = 4 * (I - 1) + 1                                    !
          NFIN = MIN(4 * I,NAT)                                     !
          READ(ICOM,15) (UJ2(J),J = NDEB,NFIN)                      !
        END DO                                                      !
!
        TEXTC(1) = TITLE                                            !
!
!================================================================
!              Reading the calculation parameters
!                      (eigenvalues)
!================================================================
!
      ELSE IF(T_CAL == 'EIGE') THEN                                  !
!
        READ(ICOM,20) NE_EIG,EIG_INI,EIG_FIN,I_DAMP                 !
!
        NE       = NE_EIG                                           !
        N_LINE_E = INT((FLOAT(NE_EIG) - SMALL) / FOUR) + 1          !
        N_LAST   = 4 - (4 * N_LINE_E - NE_EIG)                      !
!
        IF(N_LINE_E > 1) THEN                                       !
          DO JLINE = 1, N_LINE_E - 1                                !
             J = (JLINE -1 ) * 4                                    !
             READ(ICOM,21) I_SPECTRUM(J+1),I_SPECTRUM(J+2),       & !
                           I_SPECTRUM(J+3),I_SPECTRUM(J+4)          !
          END DO                                                    !
        END IF                                                      !
!
        J = 4 * (N_LINE_E - 1)                                      !
!
        READ(ICOM,21) (I_SPECTRUM(J+K), K=1,N_LAST)                 !
!
!........Reading convergence acceleration parameters
!
        READ(ICOM,22) I_PWM,METHOD,RACC,EXPO                        !
        READ(ICOM,23) N_MAX,N_ITER,N_TABLE,SHIFT                    !
        READ(ICOM,24) I_XN,I_VA,I_GN,I_WN                           !
        READ(ICOM,25) LEVIN,ALPHAR,BETAR                            !
!
        ACC = RACC                                                  !
        IF(ABS(I_PWM) <= 2) THEN                                    !
          I_ACC  = 0                                                !
          N_ITER = N_MAX                                            !
        ELSE IF(I_PWM == 3) THEN                                    !
          I_ACC  = 1                                                !
          N_ITER = N_MAX                                            !
        ELSE IF(I_PWM == -3) THEN                                   !
          I_ACC  = - 1                                              !
          N_ITER = N_MAX                                            !
        ELSE IF(I_PWM == 4) THEN                                    !
          I_ACC  = 2                                                !
        ELSE IF(I_PWM == -4) THEN                                   !
          I_ACC  = - 2                                              !
        END IF                                                      !
        IF(N_MAX < N_ITER) N_ITER = N_MAX                           !
!
        ALPHA = CMPLX(ALPHAR)                                       !
        BETA  = CMPLX(BETAR)                                        !
!
        TEXTC(2) = TITLE                                            !
!
!================================================================
!              Reading the calculation parameters
!                     (incoming electron)
!================================================================
!
      ELSE IF(T_CAL == 'INCO') THEN                                 !
!
        IF(SPECTRO == 'STM') THEN                                   !
          NAT_OLD = NAT                                             !
          NAT     = NAT_TI                                          !
        END IF                                                      !
!
        READ(ICOM,30) NO_IN,N_SCAT_IN,I_BASIS_IN,ALGO_IN            !
        READ(ICOM,31) I_REN_IN,N_REN_IN,REN_R_IN,REN_I_IN           !
!
        IF(I_BASIS_IN == 0) THEN                                    !
          IDWSPH = 0                                                !
          NO_IN  = 0                                                !
        END IF                                                      !
        IF(NO_IN < 0) NO_IN = 8                                     !
        NUMAX_IN(1) = NO_IN / 2                                     !
!
        READ(ICOM,32) ISFLIP_IN,IR_DIA_IN,ITRTL_IN,I_TEST_IN        !
!
        READ(ICOM,33) IFWD_IN,NTHOUT_IN,I_NO_IN,I_RA_IN             !
!
        IF(NTHOUT_IN == N_SCAT_IN - 1) IFWD_IN = 0                  !
!
        IF(I_RA_IN == 1) NO_IN = 0                                  !
        DO JAT = 1, NAT
          READ(ICOM,34) N_RA_IN(JAT),THFWD_IN(JAT),               & !
                        IBWD_IN(JAT),THBWD_IN(JAT)                  !
          IF(I_RA_IN == 0) THEN                                     !
            N_RA_IN(JAT)  = NO_IN                                   !
            NUMAX_IN(JAT) = NO_IN / 2                               !
          ELSE IF(I_RA_IN == 1) THEN                                !
            NUMAX_IN(JAT) = N_RA_IN(JAT) / 2                        !
            NO_IN         = MAX(N_RA_IN(JAT),NO_IN)                 !
          END IF                                                    !
        END DO                                                      !
        IF(NO_IN > NO_ST_M) IRET = 14                               !
!
        READ(ICOM,6) TESLEC                                         !
        IF(TESLEC == 'IPW,NCU') THEN                                !
          BACKSPACE ICOM                                            !
        ELSE                                                        !
          IRET = 8                                                  !
          RETURN 1                                                  !
        END IF                                                      !
!
        READ(ICOM,35) IPW_IN,NCUT_IN,PCTINT_IN,IPP_IN               !
        READ(ICOM,36) ILENGTH_IN,RLENGTH_IN,UNLENGTH_IN             !
        READ(ICOM,37) IMFP_IN,XMFP0_IN                              !
        READ(ICOM,38) IFTHET_IN,NFTHET_IN,R0_IN,R1_IN               !
!
!   Storing atomic data for mean square displacements and
!          mean free path inner calculation
!
        IF((IMSD == 1) .OR. (IMFP_IN == 1)) THEN                    !
          CALL ATDATA                                               !
        END IF                                                      !
        IF(SPECTRO == 'STM') THEN                                   !
          NAT = NAT_OLD                                             !
        END IF                                                      !
!
        TEXTC(3) = TITLE                                            !
!
!================================================================
!              Reading the calculation parameters
!                     (excited electron)
!================================================================
!
      ELSE IF(T_CAL == 'EXCI') THEN                                 !
!
        READ(ICOM,30) NO_EX,N_SCAT_EX,I_BASIS_EX,ALGO_IN            !
        READ(ICOM,31) I_REN_EX,N_REN_EX,REN_R_EX,REN_I_EX           !
!
        IF(I_BASIS_EX == 0) THEN                                    !
          IDWSPH = 0                                                !
          NO_EX  = 0                                                !
        END IF                                                      !
        IF(NO_EX < 0) NO_EX = 8                                     !
        NUMAX_EX(1) = NO_EX / 2                                     !
!
        READ(ICOM,32) ISFLIP_EX,IR_DIA_EX,ITRTL_EX,I_TEST_EX        !
!
        READ(ICOM,33) IFWD_EX,NTHOUT_EX,I_NO_EX,I_RA_EX             !
!
        IF(NTHOUT_EX == N_SCAT_EX-1) IFWD_EX = 0                    !
!
        IF(I_RA_EX == 1) NO_EX = 0                                  !
        DO JAT = 1, NAT                                             !
          READ(ICOM,34) N_RA_EX(JAT),THFWD_EX(JAT),               & !
                        IBWD_EX(JAT),THBWD_EX(JAT)                  !
          IF(I_RA_EX == 0) THEN                                     !
            N_RA_EX(JAT)  = NO_EX                                   !
            NUMAX_EX(JAT) = NO_EX / 2                               !
          ELSE IF(I_RA_EX == 1) THEN                                !
            NUMAX_EX(JAT) = N_RA_EX(JAT) / 2                        !
            NO_EX         = MAX(N_RA_EX(JAT),NO_EX)                 !
          END IF                                                    !
        END DO                                                      !
        IF(NO_EX > NO_ST_M) IRET = 14                               !
!
        READ(ICOM,6) TESLEC                                         !
        IF(TESLEC == 'IPW,NCU') THEN                                !
          BACKSPACE ICOM                                            !
        ELSE                                                        !
          IRET = 8                                                  !
          RETURN 1                                                  !
        END IF                                                      !
!
        READ(ICOM,35) IPW_EX,NCUT_EX,PCTINT_EX,IPP_EX               !
        READ(ICOM,36) ILENGTH_EX,RLENGTH_EX,UNLENGTH_EX             !
        READ(ICOM,37) IMFP_EX,XMFP0_EX                              !
        READ(ICOM,38) IFTHET_EX,NFTHET_EX,R0_EX,R1_EX               !
!
!   Storing atomic data for mean square displacements and
!          mean free path inner calculation
!
        IF((IMSD == 1) .OR. (IMFP_EX == 1)) THEN                    !
          CALL ATDATA                                               !
        END IF                                                      !
!
        TEXTC(4) = TITLE                                            !
!
!================================================================
!              Reading the calculation parameters
!                     (1st outgoing electron)
!================================================================
!
      ELSE IF(T_CAL == 'OUTG') THEN                                 !
!
        READ(ICOM,30) NO_O1,N_SCAT_O1,I_BASIS_O1,ALGO_O1            !
        READ(ICOM,31) I_REN_O1,N_REN_O1,REN_R_O1,REN_I_O1           !
!
        IF(I_BASIS_O1 == 0) THEN                                    !
          IDWSPH = 0                                                !
          NO_O1  = 0                                                !
        END IF                                                      !
        IF(NO_O1 < 0) NO_O1 = 8                                    !
        NUMAX_O1(1) = NO_O1 / 2                                     !
!
        READ(ICOM,32) ISFLIP_O1,IR_DIA_O1,ITRTL_O1,I_TEST_O1        !
!
        IF(SELRULE == '00000') THEN                                 !
          I_TEST_O1 = 1                                             ! setting the excitation
        END IF                                                      ! matrix elements to 1
!
        READ(ICOM,33) IFWD_O1,NTHOUT_O1,I_NO_O1,I_RA_O1             !
!
        IF(NTHOUT_O1 == N_SCAT_O1 - 1) IFWD_O1 = 0                  !
!
        IF(I_RA_O1 == 1) NO_O1 = 0                                  !
        DO JAT = 1, NAT                                             !
          READ(ICOM,34) N_RA_O1(JAT),THFWD_O1(JAT),               & !
                        IBWD_O1(JAT),THBWD_O1(JAT)                  !
          IF(I_RA_O1 == 0) THEN                                     !
            N_RA_O1(JAT)  = NO_O1                                   !
            NUMAX_O1(JAT) = NO_O1 / 2                               !
          ELSE IF(I_RA_O1 == 1) THEN                                !
            NUMAX_O1(JAT) = N_RA_O1(JAT) / 2                        !
            NO_O1         = MAX(N_RA_O1(JAT),NO_O1)                 !
          END IF                                                    !
        END DO                                                      !
        IF(NO_O1 > NO_ST_M) IRET = 14                               !
!
        READ(ICOM,6) TESLEC                                         !
        IF(TESLEC == 'IPW,NCU') THEN                                !
          BACKSPACE ICOM                                            !
        ELSE                                                        !
          IRET = 8                                                  !
          RETURN 1                                                  !
        END IF                                                      !
!
        READ(ICOM,35) IPW_O1,NCUT_O1,PCTINT_O1,IPP_O1               !
        READ(ICOM,36) ILENGTH_O1,RLENGTH_O1,UNLENGTH_O1             !
        READ(ICOM,37) IMFP_O1,XMFP0_O1                              !
        READ(ICOM,38) IFTHET_O1,NFTHET_O1,R0_O1,R1_O1               !
!
!   Storing atomic data for mean square displacements and
!          mean free path inner calculation
!
        IF((IMSD == 1) .OR. (IMFP_O1 == 1)) THEN                    !
          CALL ATDATA                                               !
        END IF                                                      !
!
        TEXTC(5) = TITLE                                            !
!
!================================================================
!              Reading the calculation parameters
!                     (2nd outgoing electron)
!================================================================
!
      ELSE IF(T_CAL == '2nd ') THEN                                 !
!
        READ(ICOM,30) NO_O2,N_SCAT_O2,I_BASIS_O2,ALGO_O2            !
        READ(ICOM,31) I_REN_O2,N_REN_O2,REN_R_O2,REN_I_O2           !
!
        IF(I_BASIS_O2 == 0) THEN                                    !
          IDWSPH = 0                                                !
          NO_O2  = 0                                                !
        END IF                                                      !
        IF(NO_O2 < 0) NO_O2 = 8                                     !
        NUMAX_O2(1) = NO_O2 / 2                                     !
!
        READ(ICOM,32) ISFLIP_O2,IR_DIA_O2,ITRTL_O2,I_TEST_O2        !
!
        READ(ICOM,33) IFWD_O2,NTHOUT_O2,I_NO_O2,I_RA_O2             !
!
        IF(NTHOUT_O2 == N_SCAT_O2 - 1) IFWD_O2 = 0                  !
!
        IF(I_RA_O2 == 1) NO_O2 = 0                                  !
        DO JAT = 1, NAT                                             !
          READ(ICOM,34) N_RA_O2(JAT),THFWD_O2(JAT),               & !
                        IBWD_O2(JAT),THBWD_O2(JAT)                  !
          IF(I_RA_O2 == 0) THEN                                     !
            N_RA_O2(JAT)  = NO_O2                                   !
            NUMAX_O2(JAT) = NO_O2 / 2                               !
          ELSE IF(I_RA_O2 == 1) THEN                                 !
            NUMAX_O2(JAT) = N_RA_O2(JAT) / 2                        !
            NO_O2         = MAX(N_RA_O2(JAT),NO_O2)                 !
          END IF                                                    !
        END DO                                                      !
        IF(NO_O2 > NO_ST_M) IRET = 14                               !
!
        READ(ICOM,6) TESLEC                                         !
        IF(TESLEC == 'IPW,NCU') THEN                                !
          BACKSPACE ICOM                                            !
        ELSE                                                        !
          IRET = 8                                                  !
          RETURN 1                                                  !
        END IF                                                      !
!
        READ(ICOM,35) IPW_O2,NCUT_O2,PCTINT_O2,IPP_O2               !
        READ(ICOM,36) ILENGTH_O2,RLENGTH_O2,UNLENGTH_O2             !
        READ(ICOM,37) IMFP_O2,XMFP0_O2                              !
        READ(ICOM,38) IFTHET_O2,NFTHET_O2,R0_O2,R1_O2               !
!
!   Storing atomic data for mean square displacements and
!          mean free path inner calculation
!
        IF((IMSD == 1) .OR. (IMFP_O2 == 1)) THEN                    !
          CALL ATDATA                                               !
        END IF                                                      !
!
        TEXTC(6) = TITLE                                            !
!
!================================================================
!              Reading the calculation parameters
!                     (3rd outgoing electron)
!================================================================
!
      ELSE IF(T_CAL == '3rd ') THEN                                 !
!
        READ(ICOM,30) NO_O3,N_SCAT_O3,I_BASIS_O3,ALGO_O3            !
        READ(ICOM,31) I_REN_O3,N_REN_O3,REN_R_O3,REN_I_O3           !
!
        IF(I_BASIS_O3 == 0) THEN                                    !
          IDWSPH = 0                                                !
          NO_O3  = 0                                                !
        END IF                                                      !
        IF(NO_O3 < 0) NO_O3 = 8                                     !
        NUMAX_O3(1) = NO_O3 / 2                                     !
!
        READ(ICOM,32) ISFLIP_O3,IR_DIA_O3,ITRTL_O3,I_TEST_O3        !
!
        READ(ICOM,33) IFWD_O3,NTHOUT_O3,I_NO_O3,I_RA_O3             !
!
        IF(NTHOUT_O3 == N_SCAT_O3 - 1) IFWD_O3 = 0                  !
!
        IF(I_RA_O3 == 1) NO_O3 = 0                                  !
        DO JAT = 1, NAT                                             !
          READ(ICOM,34) N_RA_O3(JAT),THFWD_O3(JAT),               & !
                        IBWD_O3(JAT),THBWD_O3(JAT)                  !
          IF(I_RA_O3 == 0) THEN                                     !
            N_RA_O3(JAT)  = NO_O3                                   !
            NUMAX_O3(JAT) = NO_O3 / 2                               !
          ELSE IF(I_RA_O3 == 1) THEN                                !
            NUMAX_O3(JAT) = N_RA_O3(JAT) / 2                        !
            NO_O3         = MAX(N_RA_O3(JAT),NO_O3)                 !
          END IF                                                    !
        END DO                                                      !
        IF(NO_O3 > NO_ST_M) IRET = 14                               !
!
        READ(ICOM,6) TESLEC                                         !
        IF(TESLEC == 'IPW,NCU') THEN                                !
          BACKSPACE ICOM                                            !
        ELSE                                                        !
          IRET = 8                                                  !
          RETURN 1                                                  !
        END IF                                                      !
!
        READ(ICOM,35) IPW_O3,NCUT_O3,PCTINT_O3,IPP_O3               !
        READ(ICOM,36) ILENGTH_O3,RLENGTH_O3,UNLENGTH_O3             !
        READ(ICOM,37) IMFP_O3,XMFP0_O3                              !
        READ(ICOM,38) IFTHET_O3,NFTHET_O3,R0_O3,R1_O3               !
!
!   Storing atomic data for mean square displacements and
!          mean free path inner calculation
!
        IF((IMSD == 1) .OR. (IMFP_O3 == 1)) THEN                    !
          CALL ATDATA                                               !
        END IF                                                      !
!
        TEXTC(7) = TITLE                                            !
!
!================================================================
!              Reading the calculation parameters
!                     for the plasmons
!================================================================
!
      ELSE IF(T_CAL == 'PLAS') THEN                                 !
!
        READ(ICOM,40) I_ANL,V_FL,WF_C,DISP                          !
!
        TEXTC(8) = TITLE                                            !
!
!================================================================
!              Reading the calculation parameters
!                     (tip related)
!================================================================
!
      ELSE IF(T_CAL == 'TIP ') THEN                                 !
!
        READ(ICOM,50) NAT_TI,A_TI,UNIT_TI,I_GR_TI                   !
!
        IF(UNIT_TI == 'ATU') THEN                                   !
          A_TI = A_TI * BOHR                                        !
        END IF                                                      !
!
        READ(ICOM,51) NPATHP_TI,VINT_TI                             !
        READ(ICOM,52) IDWSPH_TI,IATT_TI,IPRINT_TI                   !
!
        READ(ICOM,53) IMSD_TI,TD_TI,T_TI,RSJ_TI                     !
!
        NLEC = INT((NAT - SMALL) / 4) + 1                           !
!
        DO I = 1, NLEC                                              !
          NDEB = 4 * (I - 1) + 1                                    !
          NFIN = MIN(4 * I,NAT_TI)                                  !
          READ(ICOM,54) (UJ2_TI(J),J=NDEB,NFIN)                     !
        END DO                                                      !
!
        TEXTC(9) = TITLE                                            !
!
      END IF                                                        !
!
!  Rewinding the input data file for further use
!
      REWIND ICOM                                                   !
!
!  Storing IATT into the beam common for consistency
!              with previous version
!
      IATT_IN = IATT                                                !
      IATT_EX = IATT                                                !
      IATT_O1 = IATT                                                !
      IATT_O2 = IATT                                                !
      IATT_O3 = IATT                                                !
!
!  Read formats
!
   5  FORMAT(23X,A4,21X,A4)                                !
   6  FORMAT(49X,A7)                                       ! strings formats
   7  FORMAT(23X,A48)                                      !
!
  10  FORMAT(6X,I4,8X,F6.4,3X,A3,9X,I1)                    !
  11  FORMAT(9X,I1,8X,I2,6X,I4,8X,F6.2)                    !
  12  FORMAT(9X,I1,9X,I1,9X,I1)                            ! cluster formats
  13  FORMAT(9X,I1,6X,F8.3,2X,F8.3,5X,F4.2)                !
  14  FORMAT(8X,I2,6X,I4,6X,F7.2,3X,F7.2)                  !
  15  FORMAT(8X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5)              !
!
  20  FORMAT(7X,I3,6X,F7.2,3X,F7.2,6X,I1)                  !
  21  FORMAT(8X,I2,8X,I2,8X,I2,8X,I2)                      !
  22  FORMAT(8X,I2,6X,A4,9X,F7.5,2X,F6.3)                  ! eigenvalues formats
  23  FORMAT(5X,I5,6X,I4,6X,I4,8X,F6.3)                    !
  24  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)                      !
  25  FORMAT(8X,I2,6X,F7.2,3X,F7.2)                        !
!
  30  FORMAT(8X,I2,8X,I2,9X,I1,8X,A2)                      !
  31  FORMAT(9X,I1,9X,I1,6X,F8.3,2X,F8.3)                  !
  32  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)                      !
  33  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)                      !
  34  FORMAT(9X,I1,6X,F6.2,7X,I1,7X,F6.2)                  ! beams formats
  35  FORMAT(9X,I1,9X,I1,7X,F8.4,4X,I1)                    !
  36  FORMAT(9X,I1,7X,F6.2,4X,A3)                          !
  37  FORMAT(9X,I1,7X,F6.2)                                !
  38  FORMAT(9X,I1,6X,I4,8X,F6.3,5X,F6.3)                  !
!
  40  FORMAT(9X,I1,8X,A2,7X,A3,7X,A3)                      ! plasmons format
!
  50  FORMAT(6X,I4,8X,F6.4,3X,A3,9X,I1)                    !
  51  FORMAT(6X,I4,8X,F6.2)                                !
  52  FORMAT(9X,I1,9X,I1,9X,I1)                            ! tip formats
  53  FORMAT(9X,I1,6X,F7.2,3X,F7.2,6X,F4.2)                !
  54  FORMAT(8X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5)              !
!
      END SUBROUTINE READ_CALC_PARAM
!
!=======================================================================
!
      SUBROUTINE READ_IN_OUT_FILES(ICOM,IRET,*)
!
!  This subroutine reads the information about the selected spectroscopy
!    and stores it in COMMON blocks for further use
!
!
!  Input variables :
!
!
!                       ICOM      :  Fortran index of the input data file
!
!
!
!
!  Output variables :
!
!
!                       IRET      :  error code on exit
!                       *         :  return with error code IRET
!
!
!
!   Author :  D. Sébilleau
!
!
!                                           Last modified :  6 May 2021
!
      USE CLUS_ELEC
      USE INFILES
      USE INUNITS
      USE OUTFILES
      USE OUTUNITS
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  7)  ::  DUMMY
      CHARACTER (LEN =  4)  ::  STRING
!
      INTEGER, INTENT(IN)   ::  ICOM
      INTEGER, INTENT(OUT)  ::  IRET
!
      INTEGER               ::  JLINE,NLINE
!
      REAL (WP)             ::  MAX
!
      NLINE = 2000                                                  !
!
!================================================================
!                 Reading the cluster files
!================================================================
!
      DO JLINE = 1, NLINE                                           !
        READ(ICOM,5) STRING                                         !
        IF(STRING == 'SAMP') THEN                                   !
          READ(ICOM,6) DUMMY                                        !
          READ(ICOM,6) DUMMY                                        !
          READ(ICOM,6) DUMMY                                        !
          READ(ICOM,10) INFILE1,IUI1                                ! cluster file
          READ(ICOM,10) INFILE2,IUI2                                ! uj2 file
          READ(ICOM,10) INFILE3,IUI3                                ! core wave function file
          READ(ICOM,10) INFILE4,IUI4                                ! screened Coulomb file
          GO TO 613                                                 !
        END IF                                                      !
      END DO                                                        !
      IRET = 9                                                      !
      RETURN 1                                                      !
!
!================================================================
!             Reading the incoming electron files
!================================================================
!
 613  IF(INC == 1) THEN                                             !
        DO JLINE = 1, NLINE                                         !
          READ(ICOM,5) STRING                                       !
          IF(STRING == 'INCO') THEN                                 !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,10) INFILE5,IUI5                              ! T-matrix file
            READ(ICOM,10) INFILE6,IUI6                              ! K-directions file
            GO TO 614                                               !
          END IF                                                    !
        END DO                                                      !
        IRET = 9                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
!================================================================
!             Reading the excited electron files
!================================================================
!
 614  IF(EXC == 1) THEN                                             !
        DO JLINE = 1, NLINE                                         !
          READ(ICOM,5) STRING                                       !
          IF(STRING == 'EXCI') THEN                                 !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,10) INFILE7,IUI7                              ! T-matrix file
            READ(ICOM,10) INFILE8,IUI8                              ! radial integrals file
            GO TO 615                                               !
          END IF                                                    !
        END DO                                                      !
        IRET = 9                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
!================================================================
!             Reading the outgoing electron files
!================================================================
!
 615  IF(OUT1 == 1) THEN                                            !
        DO JLINE =1 , NLINE                                         !
          READ(ICOM,5) STRING                                       !
          IF(STRING == 'OUTG') THEN                                 !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,10) INFILE9,IUI9                              ! T-matrix file
            READ(ICOM,10) INFILE10,IUI10                            ! radial integrals file
            READ(ICOM,10) INFILE11,IUI11                            ! K-directions file
            GO TO 616                                               !
          END IF                                                    !
        END DO                                                      !
        IRET = 9                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
!================================================================
!             Reading the 2nd outgoing electron files
!================================================================
!
 616  IF(OUT2 == 1) THEN                                            !
        DO JLINE = 1, NLINE                                         !
          READ(ICOM,5) STRING                                       !
          IF(STRING == '2nd ') THEN                                 !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,10) INFILE12,IUI12                            ! T-matrix file
            READ(ICOM,10) INFILE13,IUI13                            ! radial integrals file
            READ(ICOM,10) INFILE14,IUI14                            ! K-directions file
            GO TO 617                                               !
          END IF                                                    !
        END DO                                                      !
        IRET = 9                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
!================================================================
!             Reading the 3rd outgoing electron files
!================================================================
!
 617  IF(OUT3 == 1) THEN                                            !
        DO JLINE = 1, NLINE                                         !
          READ(ICOM,5) STRING                                       !
          IF(STRING == '3rd ') THEN                                 !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,10) INFILE15,IUI15                            ! T-matrix file
            READ(ICOM,10) INFILE16,IUI16                            ! radial integrals file
            READ(ICOM,10) INFILE17,IUI17                            ! K-directions file
            GO TO 618                                               !
          END IF                                                    !
        END DO                                                      !
        IRET = 9                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
!================================================================
!                 Reading the tip file
!================================================================
!
 618  IF(TIP == 1) THEN                                             !
        DO JLINE = 1, NLINE                                         !
          READ(ICOM,5) STRING                                       !
          IF(STRING == 'TIP ') THEN                                 !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,6) DUMMY                                      !
            READ(ICOM,10) INFILE18,IUI18                            ! tip cluster file
            READ(ICOM,10) INFILE19,IUI19                            ! tip uj2 file
            GO TO 619                                               !
          END IF                                                    !
        END DO                                                      !
        IRET = 9                                                    !
        RETURN 1                                                    !
      END IF                                                        !
!
!================================================================
!                 Reading the output files
!================================================================
!
 619  DO JLINE = 1, NLINE                                           !
        READ(ICOM,5) STRING                                         !
        IF(STRING == 'PUT ') THEN                                   !
          READ(ICOM,6) DUMMY                                        !
          READ(ICOM,6) DUMMY                                        !
          READ(ICOM,6) DUMMY                                        !
          READ(ICOM,10) OUTFILE1,IUO1                               ! control file
          READ(ICOM,10) OUTFILE2,IUO2                               ! result file
          READ(ICOM,10) OUTFILE3,IUO3                               ! scattering amplitude file
          READ(ICOM,10) OUTFILE4,IUO4                               ! augmented clister file
          GO TO 620                                                 !
        END IF                                                      !
      END DO                                                        !
      IRET = 9                                                      !
      RETURN 1                                                      !
!
!    Setting up the unit number of the scratch files
!
 620  IUSCR = MAX(ICOM,IUI1,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,   & !
                  IUI9,IUI10,IUI11,IUI12,IUI13,IUI14,IUI15,       & !
                  IUI16,IUI17,IUI18,IUI19,IUO1,IUO2,IUO3,IUO4) + 1  !
      IUSCR2 = IUSCR + 1                                            !
!
!  Formats:
!
   5  FORMAT(37X,A4)
   6  FORMAT(A7)
!
  10  FORMAT(9X,A24,5X,I2)
!
      END SUBROUTINE READ_IN_OUT_FILES
!
!=======================================================================
!
      SUBROUTINE WRITE_READ_ERRORS(IRET,IUO1)
!
!  This subroutine prints errors encountered during the read of the input
!    data file and stops the code
!
!  Input variables :
!
!                       IRET       :  error code
!                       IUO1       :  Fortran unit number of the check file
!                       LE_MAX     :
!
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified :  6 May 2021
!
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  IRET,IUO1
!
!  Checking the error code
!
      IF(IRET == 0)  GO TO 888                                      !
      IF(IRET == 1)  GO TO 1                                        !
      IF(IRET == 2)  GO TO 2                                        !
      IF(IRET == 3)  GO TO 3                                        !
      IF(IRET == 4)  GO TO 4                                        !
      IF(IRET == 5)  GO TO 5                                        !
      IF(IRET == 6)  GO TO 6                                        !
      IF(IRET == 7)  GO TO 7                                        !
      IF(IRET == 8)  GO TO 8                                        !
      IF(IRET == 9)  GO TO 9                                        !
      IF(IRET == 10) GO TO 10                                       !
      IF(IRET == 11) GO TO 11                                       !
      IF(IRET == 12) GO TO 12                                       !
      IF(IRET == 13) GO TO 13                                       !
      IF(IRET == 14) GO TO 14                                       !
!
!  Stops
!
   1  WRITE(IUO1,101)                                               !
      STOP                                                          !
   2  WRITE(IUO1,102)                                               !
      STOP                                                          !
   3  WRITE(IUO1,103)                                               !
      STOP                                                          !
   4  WRITE(IUO1,104)                                               !
      STOP                                                          !
   5  WRITE(IUO1,105)                                               !
      STOP                                                          !
   6  WRITE(IUO1,106)                                               !
      STOP                                                          !
   7  WRITE(IUO1,107)                                               !
      STOP                                                          !
   8  WRITE(IUO1,108)                                               !
      STOP                                                          !
   9  WRITE(IUO1,109)                                               !
      STOP                                                          !
  10  WRITE(IUO1,110)                                               !
      STOP                                                          !
  11  WRITE(IUO1,111)                                               !
      STOP                                                          !
  12  WRITE(IUO1,112)                                               !
      STOP                                                          !
  13  WRITE(IUO1,113)                                               !
      STOP                                                          !
  14  WRITE(IUO1,114)                                               !
      STOP                                                          !
!
!  Formats
!
 101  FORMAT(///,'<<<<<<<<<<  (NAT,NE) > (NATP_M,NE_M)',                &
             ' - CHECK THE DIMENSIONING >>>>>>>>>>')
 102  FORMAT(///,22X,' <<<<<<<<<<  THIS STRUCTURE DOES NOT EXIST ',     &
      '  >>>>>>>>>>')
 103  FORMAT(///,4X,' <<<<<<<<<<  ONLY ONE OF THE VALUES IPHI,ITHETA ', &
      'AND IE CAN BE EQUAL TO 1  >>>>>>>>>>')
 104  FORMAT(///,8X,' <<<<<<<<<<  CHANGE THE DIMENSIONING OF PCREL ',   &
      'IN MAIN AND READ_DATA  >>>>>>>>>>')
 105  FORMAT(///,8X,'<<<<<<<<<<  WRONG NAME FOR THE INITIAL STATE',     &
            '  >>>>>>>>>>')
 106  FORMAT(///,'<<<<<<<<<<  LI IS LARGER THAN LI_M - ',               &
             'CHECK THE DIMENSIONING  >>>>>>>>>>')
 107  FORMAT(///,'<<<<<<<<<<  SPIN-ORBIT COMPONENT NOT CONSISTENT WITH',&
             ' THE VALUE OF LI  >>>>>>>>>>')
 108  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES THFWD DOES NOT ',     &
             'CORRESPOND TO NAT  >>>>>>>>>>')
 109  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES UJ2 DOES NOT ',       &
             'CORRESPOND TO NAT  >>>>>>>>>>')
 110  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES ATBAS DOES NOT ',     &
             'CORRESPOND TO NAT  >>>>>>>>>>')
 111  FORMAT(///,'<<<<<<<<<<  LI OR IMOD NOT CONSISTENT BETWEEN ',      &
             'PED AND AED FOR COINCIDENCE CALCULATION  >>>>>>>>>>')
 112  FORMAT(///,'<<<<<<<<<<  THE EXTERNAL DIRECTIONS FILE IS ',        &
             'NOT CONSISTENT WITH THE INPUT DATA FILE  >>>>>>>>>>')
 113  FORMAT(///,'<<<<<<<<<<      EXCURSIONS OF ANGLES SHOULD ',        &
             ' BE IDENTICAL      >>>>>>>>>>',/,'<<<<<<<<<<     ',       &
             'FOR BOTH ELECTRONS IN CLUSTER ROTATION MODE',             &
              '     >>>>>>>>>>',//)
 114  FORMAT(///,'<<<<<<<<<<  NO_ST_M IS TOO SMALL IN THE .inc FILE ',  &
             '>>>>>>>>>>',//)
!
 888  RETURN                                                        !
!
      END SUBROUTINE WRITE_READ_ERRORS
!
!=======================================================================
!
      SUBROUTINE WRITE_CALC_PARAM(IUO1,JFILE,SPECTRO)
!
!  This subroutine writes the calculation parameters
!
!
!  Input variables :
!
!
!                       IUO1      :  Fortran index of the control file
!                       JFILE     :  index of the input data file
!                       SPECTRO   :  spectroscopy condidered
!                                            --> 'PED'   : photoelectron diffraction
!                                            --> 'LED'   : low-energy electron diffraction
!                                            --> 'XAS'   : X-ray absorption spectroscopy
!                                            --> 'AED'   : Auger electron diffraction
!                                            --> 'APC'   : Auger photoelectron coincidence spectroscopy
!                                            --> 'E2E'   : (e,2e) spectroscopy
!                                            --> 'E3E'   : (e,3e) spectroscopy
!                                            --> 'ELS'   : electron energy loss spectroscopy
!                                            --> 'RES'   : resonant elastic X-ray scattering
!                                            --> 'STM'   : scanning tunneling microscopy
!                                            --> 'PLS'   : photoelectron energy loss spectroscopy
!                                            --> 'BEM'   : ballistic energy electron microscopy
!                                            --> 'EIG'   : eigenvalue calculation
!
!
!
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 14  Jun 2021
!
      USE REAL_NUMBERS,    ONLY : FOUR
      USE CONSTANTS_P1,    ONLY : BOHR
!
      USE ALGORITHMS
!
      USE AMPLI_IN
      USE AMPLI_EX
      USE AMPLI_O1
      USE AMPLI_O2
      USE AMPLI_O3
!
      USE APPROXIMATIONS_IN
      USE APPROXIMATIONS_EX
      USE APPROXIMATIONS_O1
      USE APPROXIMATIONS_O2
      USE APPROXIMATIONS_O3
!
      USE LINLBD_IN
      USE LINLBD_EX
      USE LINLBD_O1
      USE LINLBD_O2
      USE LINLBD_O3
!
      USE MFP_IN
      USE MFP_EX
      USE MFP_O1
      USE MFP_O2
      USE MFP_O3
!
      USE RA_IN
      USE RA_EX
      USE RA_O1
      USE RA_O2
      USE RA_O3
!
      USE RENORM_IN
      USE RENORM_EX
      USE RENORM_O1
      USE RENORM_O2
      USE RENORM_O3
!
      USE SPIN_IN
      USE SPIN_EX
      USE SPIN_O1
      USE SPIN_O2
      USE SPIN_O3
!
      USE TEST_IN
      USE TEST_EX
      USE TEST_O1
      USE TEST_O2
      USE TEST_O3
!
      USE ATOMS
      USE ATOMS_TIP
      USE CLUSTER
      USE CLUSTER_TIP
      USE CLUS_ELEC
      USE CONV_ACC
      USE CONV_TYPE
      USE DAMPING
      USE DEB_WAL_CLU
      USE DEB_WAL_TIP
      USE EIGEN
      USE INIT_L
      USE NON_BULK
      USE PLAS_CAL
      USE SPECTRUM
      USE TEMP_LOOP
      USE TESTS
      USE TESTS_TI
      USE TEXT_CALC
      USE TMP_01
      USE TMP_02
      USE TMP_03
      USE TMP_04
      USE TMP_05
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 3)   ::  SPECTRO
!
      INTEGER, INTENT(IN)   ::  IUO1,JFILE
!
      INTEGER               ::  I,J,K
      INTEGER               ::  NLEC,NDEB,NFIN
      INTEGER               ::  NUJ,NAT_OLD,N_LAST
      INTEGER               ::  JAT,JLINE
      INTEGER               ::  N_LINE_E,I_AMP
!
      INTEGER               ::  INT
!
      REAL (WP)             ::  SMALL,PIS180
      REAL (WP)             ::  A_OLD
!
      REAL (WP)             ::  FLOAT
!
      DATA SMALL      / 0.0001E0_WP /
      DATA PIS180     / 0.01745329251994329576923690768488612713E0_WP /
!
      I_AMP =  MAX(I_AMP_IN,I_AMP_EX,I_AMP_O1,I_AMP_O2,I_AMP_O3)    !
!
!================================================================
!              Writing the calculation parameters
!                           (cluster)
!================================================================
!
      IF(CLU == 1) THEN                                             !
!
        WRITE(IUO1,5)  TEXTC(1)                                     !
        WRITE(IUO1,10) NAT,A,UNIT,I_GR                              !
        WRITE(IUO1,11) ISOM,NONVOL(JFILE),NPATHP,VINT               !
        WRITE(IUO1,12) IDWSPH,IATT,IPRINT                           !
!
        WRITE(IUO1,13) IMSD,TD,T,RSJ                                !
        WRITE(IUO1,14) ITEMP,NTEMP,TEMP0,TEMP1                      !
!
        NLEC = INT((NAT - SMALL) / 4) + 1                           !
!
        DO I  = 1, NLEC                                             !
          NDEB = 4 * (I - 1) + 1                                    !
          NFIN = 4 * I                                              !
          IF(I == NLEC) NFIN = NAT                                  !
          NUJ = NFIN - NDEB + 1                                     !
          IF(NUJ == 1) WRITE(IUO1,15) (UJ2(J),J=NDEB,NFIN)          !
          IF(NUJ == 2) WRITE(IUO1,16) (UJ2(J),J=NDEB,NFIN)          !
          IF(NUJ == 3) WRITE(IUO1,17) (UJ2(J),J=NDEB,NFIN)          !
          IF(NUJ == 4) WRITE(IUO1,18) (UJ2(J),J=NDEB,NFIN)          !
        END DO                                                      !
!
        DO J = 1, NATM                                              !
          UJ2(J) = UJ2(J) / (A * A)                                 !
        END DO                                                      !
!
        IF(I_AMP == 1) WRITE(6,19)                                  !
!
      END IF                                                        !
!
!================================================================
!              Writing the calculation parameters
!                      (eigenvalues)
!================================================================
!
      IF(EIG == 1) THEN                                             !
!
        WRITE(IUO1,20)                                              !
        WRITE(IUO1,5) TEXTC(2)                                      !
        WRITE(IUO1,21) NE_EIG,EIG_INI,EIG_FIN,I_DAMP                !
!
        N_LINE_E = INT((FLOAT(NE_EIG) - SMALL) / FOUR) + 1          !
        DO JLINE = 1, N_LINE_E - 1                                  !
          J = (JLINE -1 ) * 4                                       !
          WRITE(IUO1,22) I_SPECTRUM(J+1),I_SPECTRUM(J+2),         & !
                         I_SPECTRUM(J+3),I_SPECTRUM(J+4)            !
        END DO                                                      !
        N_LAST = 4 - (4 * N_LINE_E - NE_EIG)                        !
        J      = 4 * (N_LINE_E - 1)                                 !
        WRITE(IUO1,22) (I_SPECTRUM(J+K),K=1,N_LAST)                 !
!
        WRITE(IUO1,23) I_PWM,METHOD,RACC,EXPO                       !
        WRITE(IUO1,24) N_MAX,N_ITER,N_TABLE,SHIFT                   !
        WRITE(IUO1,25) I_XN,I_VA,I_GN,I_WN                          !
        WRITE(IUO1,26) LEVIN,ALPHAR,BETAR                           !
!
      END IF                                                        !
!
!================================================================
!              Writing the calculation parameters
!                     (incoming electron)
!================================================================
!
      IF(INC == 1) THEN                                             !
!
        WRITE(IUO1,5) TEXTC(3)                                      !
!
        WRITE(IUO1,30) NO_IN,N_SCAT_IN,I_BASIS_IN,ALGO_IN           !
        WRITE(IUO1,31) I_REN_IN,N_REN_IN,REN_R_IN,REN_I_IN          !
!
        IF(SPECTRO == 'STM') THEN                                   !
          NAT_OLD = NAT                                             !
          A_OLD   = A                                               !
          NAT     = NAT_TI                                          !
          A       = A_TI                                            !
        END IF                                                      !
!
        WRITE(IUO1, 32) ISFLIP_IN,IR_DIA_IN,ITRTL_IN,I_TEST_IN      !
!
        IF(ISFLIP_IN == 0) THEN                                     !
          NSTEP_IN = 3                                              !
        ELSE                                                        !
          NSTEP_IN = 1                                              !
        END IF                                                      !
!
        WRITE(IUO1,33) IFWD_IN,NTHOUT_IN,I_NO_IN,I_RA_IN            !
        DO JAT = 1, NAT                                             !
          WRITE(IUO1,34) N_RA_IN(JAT),THFWD_IN(JAT),IBWD_IN(JAT), & !
                         THBWD_IN(JAT)                              !
          RTHFWD_IN(JAT) = THFWD_IN(JAT) * PIS180                   !
          RTHBWD_IN(JAT) = THBWD_IN(JAT) * PIS180                   !
        END DO                                                      !
        WRITE(IUO1,35) IPW_IN,NCUT_IN,PCTINT_IN,IPP_IN              !
        WRITE(IUO1,36) ILENGTH_IN,RLENGTH_IN,UNLENGTH_IN            !
        WRITE(IUO1,37) IMFP_IN,XMFP0_IN                             !
!
        IF(UNLENGTH_IN == 'ATU') RLENGTH_IN = RLENGTH_IN * BOHR / A !
        IF(UNLENGTH_IN == 'ANG') RLENGTH_IN = RLENGTH_IN / A        !
!
        IF(SPECTRO == 'STM') THEN                                   !
          NAT = NAT_OLD                                             !
          A   = A_OLD                                               !
        END IF                                                      !
!
      END IF                                                        !
!
!================================================================
!              Writing the calculation parameters
!                     (excited electron)
!================================================================
!
      IF(EXC == 1) THEN                                             !
!
        WRITE(IUO1,5) TEXTC(4)                                      !
        IF(SPECTRO /= 'EIG') THEN                                   !
!
          WRITE(IUO1,30) NO_EX,N_SCAT_EX,I_BASIS_EX,ALGO_EX         !
          WRITE(IUO1,31) I_REN_EX,N_REN_EX,REN_R_EX,REN_I_EX        !
!
          IF(SPECTRO == 'XAS') N_SCAT_EX = N_SCAT_EX + 1            !
!
          WRITE(IUO1,32) ISFLIP_EX,IR_DIA_EX,ITRTL_EX,I_TEST_EX     !
!
          IF(ISFLIP_EX == 0) THEN                                   !
            NSTEP_EX = 3                                            !
          ELSE                                                      !
            NSTEP_EX = 1                                            !
          END IF                                                    !
!
          WRITE(IUO1,33) IFWD_EX,NTHOUT_EX,I_NO_EX,I_RA_EX          !
          DO JAT = 1, NAT                                           !
            WRITE(IUO1,34) N_RA_EX(JAT),THFWD_EX(JAT),IBWD_EX(JAT),&!
                           THBWD_EX(JAT)                            !
            RTHFWD_EX(JAT) = THFWD_EX(JAT) * PIS180                 !
            RTHBWD_EX(JAT) = THBWD_EX(JAT) * PIS180                 !
          END DO                                                    !
          WRITE(IUO1,35) IPW_EX,NCUT_EX,PCTINT_EX,IPP_EX            !
          WRITE(IUO1,36) ILENGTH_EX,RLENGTH_EX,UNLENGTH_EX          !
        END IF                                                      !
        WRITE(IUO1,37) IMFP_EX,XMFP0_EX                             !
!
        IF(SPECTRO /= 'EIG') THEN                                   !
          IF(UNLENGTH_EX == 'ATU') RLENGTH_EX = RLENGTH_EX * BOHR / A
          IF(UNLENGTH_EX == 'ANG') RLENGTH_EX = RLENGTH_EX / A      !
        END IF                                                      !
!
      END IF                                                        !
!
!================================================================
!              Writing the calculation parameters
!                     (1st outgoing electron)
!================================================================
!
      IF(OUT1 == 1) THEN                                            !
!
        WRITE(IUO1,5) TEXTC(5)                                      !
!
        WRITE(IUO1,30) NO_O1,N_SCAT_O1,I_BASIS_O1,ALGO_O1           !
        WRITE(IUO1,31) I_REN_O1,N_REN_O1,REN_R_O1,REN_I_O1          !
!
        WRITE(IUO1,32) ISFLIP_O1,IR_DIA_O1,ITRTL_O1,I_TEST_O1       !
!
        IF(ISFLIP_O1 == 0) THEN                                     !
          NSTEP_O1 = 3                                              !
        ELSE                                                        !
          NSTEP_O1 = 1                                              !
        END IF                                                      !
!
        WRITE(IUO1,33) IFWD_O1,NTHOUT_O1,I_NO_O1,I_RA_O1            !
        DO JAT = 1, NAT                                             !
          WRITE(IUO1,34) N_RA_O1(JAT),THFWD_O1(JAT),IBWD_O1(JAT), & !
                         THBWD_O1(JAT)                              !
          RTHFWD_O1(JAT) = THFWD_O1(JAT) * PIS180                   !
          RTHBWD_O1(JAT) = THBWD_O1(JAT) * PIS180                   !
        END DO                                                      !
        WRITE(IUO1,35) IPW_O1,NCUT_O1,PCTINT_O1,IPP_O1              !
        WRITE(IUO1,36) ILENGTH_O1,RLENGTH_O1,UNLENGTH_O1            !
        WRITE(IUO1,37) IMFP_O1,XMFP0_O1                             !
!
        IF(UNLENGTH_O1 == 'ATU') RLENGTH_O1 = RLENGTH_O1 * BOHR / A !
        IF(UNLENGTH_O1 == 'ANG') RLENGTH_O1 = RLENGTH_O1 / A        !
!
      END IF                                                        !
!
!================================================================
!              Writing the calculation parameters
!                     (2nd outgoing electron)
!================================================================
!
      IF(OUT2 == 1) THEN                                            !
!
        WRITE(IUO1,5) TEXTC(6)                                      !
!
        WRITE(IUO1,30) NO_O2,N_SCAT_O2,I_BASIS_O2,ALGO_O2           !
        WRITE(IUO1,31) I_REN_O2,N_REN_O2,REN_R_O2,REN_I_O2          !
!
        WRITE(IUO1,32) ISFLIP_O2,IR_DIA_O2,ITRTL_O2,I_TEST_O2       !
!
        IF(ISFLIP_O2 == 0) THEN                                     !
          NSTEP_O2 = 3                                              !
        ELSE                                                        !
          NSTEP_O2 = 1                                              !
        END IF                                                      !
!
        WRITE(IUO1,33) IFWD_O2,NTHOUT_O2,I_NO_O2,I_RA_O2            !
        DO JAT = 1, NAT                                             !
          WRITE(IUO1,34) N_RA_O2(JAT),THFWD_O2(JAT),IBWD_O2(JAT), & !
                         THBWD_O2(JAT)                              !
          RTHFWD_O2(JAT) = THFWD_O2(JAT) * PIS180                   !
          RTHBWD_O2(JAT) = THBWD_O2(JAT) * PIS180                   !
        END DO                                                      !
        WRITE(IUO1,35) IPW_O2,NCUT_O2,PCTINT_O2,IPP_O2              !
        WRITE(IUO1,36) ILENGTH_O2,RLENGTH_O2,UNLENGTH_O2            !
        WRITE(IUO1,37) IMFP_O2,XMFP0_O2                             !
!
        IF(UNLENGTH_O2 == 'ATU') RLENGTH_O2 = RLENGTH_O2 *BOHR / A  !
        IF(UNLENGTH_O2 == 'ANG') RLENGTH_O2 = RLENGTH_O2 / A        !
!
      END IF                                                        !
!
!================================================================
!              Writing the calculation parameters
!                     (3rd outgoing electron)
!================================================================
!
      IF(OUT3 == 1) THEN                                            !
!
        WRITE(IUO1,5) TEXTC(7)                                      !
!
        WRITE(IUO1,30) NO_O3,N_SCAT_O3,I_BASIS_O3,ALGO_O3           !
        WRITE(IUO1,31) I_REN_O3,N_REN_O3,REN_R_O3,REN_I_O3          !
!
        WRITE(IUO1,32) ISFLIP_O3,IR_DIA_O3,ITRTL_O3,I_TEST_O3       !
!
        IF(ISFLIP_O3 == 0) THEN                                     !
          NSTEP_O3 = 3                                              !
        ELSE                                                        !
          NSTEP_O3 = 1                                              !
        END IF                                                      !
!
        WRITE(IUO1,33) IFWD_O3,NTHOUT_O3,I_NO_O3,I_RA_O3            !
        DO JAT = 1, NAT                                             !
          WRITE(IUO1,34) N_RA_O3(JAT),THFWD_O3(JAT),IBWD_O3(JAT), & !
                         THBWD_O3(JAT)                              !
          RTHFWD_O3(JAT) = THFWD_O3(JAT) * PIS180                   !
          RTHBWD_O3(JAT) = THBWD_O3(JAT) * PIS180                   !
        END DO                                                      !
        WRITE(IUO1,35) IPW_O3,NCUT_O3,PCTINT_O3,IPP_O3              !
        WRITE(IUO1,36) ILENGTH_O3,RLENGTH_O3,UNLENGTH_O3            !
        WRITE(IUO1,37) IMFP_O3,XMFP0_O3                             !
!
        IF(UNLENGTH_O3 == 'ATU') RLENGTH_O3 = RLENGTH_O3 * BOHR / A !
        IF(UNLENGTH_O3 == 'ANG') RLENGTH_O3 = RLENGTH_O3 / A        !
!
      END IF                                                        !
!
!================================================================
!              Writing the calculation parameters
!                     for the plasmons
!================================================================
!
      IF(PLA == 1) THEN                                             !
!
        WRITE(IUO1,5) TEXTC(8)                                      !
!
        WRITE(IUO1,40) I_ANL,V_FL,WF_C,DISP                         !
!
      END IF                                                        !
!
!================================================================
!              Writing the calculation parameters
!                     (tip related)
!================================================================
!
      IF(TIP == 1) THEN                                             !
!
        WRITE(IUO1,5) TEXTC(9)                                      !
!
        WRITE(IUO1,50) NAT_TI,A_TI,UNIT_TI,I_GR_TI                  !
        WRITE(IUO1,51) NPATHP_TI,VINT_TI                            !
        WRITE(IUO1,52) IDWSPH_TI,IATT_TI,IPRINT_TI                  !
        WRITE(IUO1,53) IMSD_TI,TD_TI,T_TI,RSJ_TI                    !
        DO I = 1, NLEC                                              !
          NDEB = 4 * (I - 1) + 1                                    !
          NFIN = 4 * I                                              !
          IF(I == NLEC) NFIN = NAT                                  !
          NUJ = NFIN - NDEB + 1                                     !
          IF(NUJ == 1) WRITE(IUO1,54) (UJ2_TI(J),J=NDEB,NFIN)       !
          IF(NUJ == 2) WRITE(IUO1,55) (UJ2_TI(J),J=NDEB,NFIN)       !
          IF(NUJ == 3) WRITE(IUO1,56) (UJ2_TI(J),J=NDEB,NFIN)       !
          IF(NUJ == 4) WRITE(IUO1,57) (UJ2_TI(J),J=NDEB,NFIN)       !
        END DO                                                      !
!
        DO J = 1, NATM                                              !
          UJ2_TI(J) = UJ2_TI(J) / (A_TI * A_TI)                     !
        END DO                                                      !
!
      END IF                                                        !
!
!  Write formats:
!
   5  FORMAT(///,21X,A48,/)                                            ! string format
!
  10  FORMAT(9X,I4,8X,F6.4,3X,A3,9X,I1,9X,'NAT,A,UNIT,I_GR')           !
  11  FORMAT(12X,I1,8X,I2,6X,I4,7X,F6.2,6X,'ISOM,NONVOL,NPATH,VINT')   !
  12  FORMAT(12X,I1,9X,I1,9X,I1,19X,'IDWSPH,IATT,IPRINT')              !
  13  FORMAT(12X,I1,6X,F8.3,2X,F8.3,5X,F4.2,6X,'IMSD,TD,T,RSJ')        !
  14  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'ITEMP,NTEMP,TEMP1,TEMP2')! cluster formats
  15  FORMAT(11X,F8.5,33X,'UJ2(NAT)  : SUBSTRATE')                     !
  16  FORMAT(11X,F8.5,2X,F8.5,23X,'UJ2(NAT)  : SUBSTRATE')             !
  17  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,13X,'UJ2(NAT)  : SUBSTRATE')     !
  18  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5,3X,'UJ2(NAT)  : ',     & !
             'SUBSTRATE')                                              !
  19  FORMAT(//,20X,'THIS CALCULATION OUTPUTS ALSO THE AMPLITUDES')    !
!
  20  FORMAT(///,22X,'TYPE OF CALCULATION : EIGENVALUE ANALYSIS')      !
  21  FORMAT(10X,I3,6X,F7.2,3X,F7.2,6X,I1,9X,'NE_EIG,EIG_INI,',      & !
             'EIG_FIN,I_DAMP')                                         !
  22  FORMAT(11X,I2,8X,I2,8X,I2,8X,I2,9X,'I_SPECTRUM(NE)')             ! eigenvalues formats
  23  FORMAT(11X,I2,6X,A4,9X,F7.5,2X,F6.3,5X,'I_PWM,METHOD,ACC,EXPO')  !
  24  FORMAT(8X,I5,6X,I4,6X,I4,8X,F6.3,5X,'N_MAX,N_ITER,N_TABLE,SHIFT')!
  25  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'I_XN,I_VA,I_GN,I_WN')        !
  26  FORMAT(11X,I2,6X,F7.2,3X,F7.2,16X,'L,ALPHA,BETA')                !
!
  30  FORMAT(11X,I2,8X,I2,9X,I1,8X,A2,9X,'NO,N_SCAT,I_BASIS,ALGO')     !
  31  FORMAT(12X,I1,9X,I1,6X,F8.3,2X,F8.3,5X,'I_REN,N_REN,REN_R,REN_I')!
  32  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'ISFLIP,IR_DIA,ITRTL,I_TEST') !
  33  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'IFWD,NTHOUT,I_NO,I_RA')      ! electron beams formats
  34  FORMAT(12X,I1,7X,F6.2,6X,I1,7X,F6.2,6X,'N_RA(NAT),THFWD(NAT)', & !
             ',IBWD(NAT),THBWD(NAT)')                                  !
  35  FORMAT(12X,I1,9X,I1,7X,F8.4,4X,I1,9X,'IPW,NCUT,PCTINT,IPP')      !
  36  FORMAT(12X,I1,7X,F6.2,4X,A3,19X,'ILENGTH,RLENGTH,UNLENGTH')      !
  37  FORMAT(12X,I1,7X,F6.2,26X,'IMFP,XMFP0')
!
  40  FORMAT(12X,I1,8X,A2,7X,A3,7X,A3,9X,'I_ANL,V_FL,WF_C,DISP')       ! plasmons format
!
  50  FORMAT(9X,I4,8X,F6.4,3X,A3,9X,I1,9X,'NAT,A,UNIT,I_GR')           !
  51  FORMAT(9X,I4,7X,F6.2,26X,' NPATH,VINT')                          !
  52  FORMAT(12X,I1,9X,I1,9X,I1,19X,'IDWSPH,IATT,IPRINT')              !
  53  FORMAT(12X,I1,6X,F7.2,3X,F7.2,6X,F4.2,6X,'IMSD,TD,T,RSJ')        ! tip formats
  54  FORMAT(11X,F8.5,33X,'UJ2(NAT)  : TIP')                           !
  55  FORMAT(11X,F8.5,2X,F8.5,23X,'UJ2(NAT)  : TIP')                   !
  56  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,13X,'UJ2(NAT)  : TIP')           !
  57  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5,3X,'UJ2(NAT)  : TIP')    !
!
      END SUBROUTINE WRITE_CALC_PARAM
!
!=======================================================================
!
      SUBROUTINE WRITE_SPECTRO_INFO(IUO1)
!
!  This subroutine reads the information about the selected spectroscopy
!    and stores it in COMMON blocks for further use
!
!
!  Input variables :
!
!
!                       IUO1      :  Fortran index of the control file
!
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 27 May 2021
!
!
      USE AVER_IN
      USE ENERGY_IN
      USE PHI_IN
      USE POLAR_IN
      USE THETA_IN
!
      USE ENERGY_EX
!
      USE AVER_O1
      USE ENERGY_O1
      USE PHI_O1
      USE POLAR_O1
      USE THETA_O1
!
      USE AVER_O2
      USE ENERGY_O2
      USE PHI_O2
      USE POLAR_O2
      USE THETA_O2
!
      USE AVER_O3
      USE ENERGY_O3
      USE PHI_O3
      USE POLAR_O3
      USE THETA_O3
!
      USE APPROX_CS
      USE EXP_TYPE
      USE HEADER
      USE INIT_J
      USE INIT_L
      USE INIT_M
      USE PLAS_EXP
      USE REXS_EXP
      USE TEMP01
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  IUO1
!
      INTEGER               ::  J_PLA
!
      REAL (WP)             ::  ABS
!
!  Photoelectron diffraction case:
!
      IF(SPECTRO == 'PED') THEN                                     !
!
        IF(IPHI_O1 == 1) THEN                                       !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,50)                                          !
          ELSE                                                      !
            WRITE(IUO1,51)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(IE_O1 == 1) WRITE(IUO1,52)                               !
        IF(ITHETA_O1 == 1) WRITE(IUO1,53)                           !
!
        WRITE(IUO1,10) NI,NLI,S_O,SELRULE,I_SO                      !
        WRITE(IUO1,11) IMOD                                         !
!
        IF(ABS(THETA0_IN) > 90.0E0_WP) THEN                         !
          WRITE(IUO1,200)                                           !
          STOP                                                      !
        END IF                                                      !
!
        IF(ITHETA_O1 == 1) THEN                                     !
          IF( (THETA0_O1 < -90.0E0_WP) .OR.                       & !
              (THETA1_O1 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,201)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O1 direction to account for detector entrance slit
!
        IF(I_AVER_O1 > 3)  I_AVER_O1 = 3                            !
        IF(I_AVER_O1 < 0)  I_AVER_O1 = 0                            !
        IF(I_AVER_O1 == 0) NDIR_O1   = 1                            !
        IF(I_AVER_O1 == 1) NDIR_O1   = 5                            !
        IF(I_AVER_O1 == 2) NDIR_O1   = 13                           !
        IF(I_AVER_O1 == 3) NDIR_O1   = 49                           !
!
      ELSE IF(SPECTRO == 'LED') THEN                                !
!
        IF(IPHI_IN == 1) THEN                                       !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,60)                                          !
          ELSE                                                      !
            WRITE(IUO1,61)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(IE_IN == 1) WRITE(IUO1,62)                               !
        IF(ITHETA_IN == 1) WRITE(IUO1,63)                           !
!
        WRITE(IUO1,11) IMOD                                         !
!
        IF(ITHETA_IN == 1) THEN                                     !
          IF( (THETA0_IN < -90.0E0_WP) .OR.                       & !
              (THETA1_IN > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,200)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
        IF(ITHETA_O1 == 1) THEN                                     !
          IF( (THETA0_O1 < -90.0E0_WP) .OR.                       & !
              (THETA1_O1 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,201)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O1 direction to account for detector entrance slit
!
        IF(I_AVER_O1 > 3)  I_AVER_O1 = 3                            !
        IF(I_AVER_O1 < 0)  I_AVER_O1 = 0                            !
        IF(I_AVER_O1 == 0) NDIR_O1   = 1                            !
        IF(I_AVER_O1 == 1) NDIR_O1   = 5                            !
        IF(I_AVER_O1 == 2) NDIR_O1   = 13                           !
        IF(I_AVER_O1 == 3) NDIR_O1   = 49                           !
!
      ELSE IF(SPECTRO == 'XAS') THEN                                !
!
        WRITE(IUO1,70)                                              !
        WRITE(IUO1,12) EDGE,NEDGE,SELRULE                           !
        WRITE(IUO1,11) IMOD                                         !
!
      ELSE IF(SPECTRO == 'AED') THEN                                !
!
        IF(IPHI_O1 == 1) THEN                                       !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,80)                                          !
          ELSE                                                      !
            WRITE(IUO1,81)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(ITHETA_O1 == 1) WRITE(IUO1,82)                           !
!
        WRITE(IUO1,13) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A !
        WRITE(IUO1,14) I_MULT,IM1,MULT,IM2                          !
        WRITE(IUO1,11) IMOD                                         !
!
!  Averaging over O1 direction to account for detector entrance slit
!
        IF(I_AVER_O1 > 3)  I_AVER_O1 = 3                            !
        IF(I_AVER_O1 < 0)  I_AVER_O1 = 0                            !
        IF(I_AVER_O1 == 0) NDIR_O1   = 1                            !
        IF(I_AVER_O1 == 1) NDIR_O1   = 5                            !
        IF(I_AVER_O1 == 2) NDIR_O1   = 13                           !
        IF(I_AVER_O1 == 3) NDIR_O1   = 49                           !
!
      ELSE IF(SPECTRO == 'APC') THEN                                !
!
        IF(IPHI_O1 == 1) WRITE(IUO1,90)                             !
        IF(ITHETA_O1 == 1) WRITE(IUO1,91)                           !
!
        WRITE(IUO1,92)                                              !
        WRITE(IUO1,10) NI,NLI,S_O,SELRULE,I_SO                      !
        WRITE(IUO1,93)                                              !
        WRITE(IUO1,13) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A !
        WRITE(IUO1,14) I_MULT,IM1,MULT,IM2                          !
        WRITE(IUO1,11) IMOD                                         !
        WRITE(IUO1,94)                                              !
!
        IF(ABS(THETA0_IN) > 90.0E0_WP) THEN                         !
          WRITE(IUO1,200)                                           !
          STOP                                                      !
        END IF                                                      !
!
        IF(ITHETA_O1 == 1) THEN                                     !
          IF( (THETA0_O1 < -90.0E0_WP) .OR.                       & !
              (THETA1_O1 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,201)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O1 direction to account for detector entrance slit
!
        IF(I_AVER_O1 > 3)  I_AVER_O1 = 3                            !
        IF(I_AVER_O1 < 0)  I_AVER_O1 = 0                            !
        IF(I_AVER_O1 == 0) NDIR_O1   = 1                            !
        IF(I_AVER_O1 == 1) NDIR_O1   = 5                            !
        IF(I_AVER_O1 == 2) NDIR_O1   = 13                           !
        IF(I_AVER_O1 == 3) NDIR_O1   = 49                           !
!
        IF(ITHETA_O2 == 1) THEN                                     !
          IF( (THETA0_O2 < -90.0E0_WP) .OR.                       & !
              (THETA1_O2 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,202)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O2 direction to account for detector entrance slit
!
        IF(I_AVER_O2 > 3)  I_AVER_O2 = 3                            !
        IF(I_AVER_O2 < 0)  I_AVER_O2 = 0                            !
        IF(I_AVER_O2 == 0) NDIR_O2   = 1                            !
        IF(I_AVER_O2 == 1) NDIR_O2   = 5                            !
        IF(I_AVER_O2 == 2) NDIR_O2   = 13                           !
        IF(I_AVER_O2 == 3) NDIR_O2   = 49                           !
!
        WRITE(IUO1,94)                                              !
!
      ELSE IF(SPECTRO == 'E2E') THEN                                !
!
        IF(IPHI_O1 == 1) WRITE(IUO1,100)                            !
        IF(ITHETA_O1 == 1) WRITE(IUO1,101)                          !
        IF(IE_O1 == 1) WRITE(IUO1,102)                              !
!
        WRITE(IUO1,18) NI,NLI,S_O,ISELRULE,I_SO                     !
        WRITE(IUO1,15) MODE,CS_TYPE                                 !
        WRITE(IUO1,11) IMOD                                         !
!
        IF(ABS(THETA0_IN) > 90.0E0_WP) THEN                         !
          WRITE(IUO1,200)                                           !
          STOP                                                      !
        END IF                                                      !
!
        IF(ITHETA_O1 == 1) THEN                                     !
          IF( (THETA0_O1 < -90.0E0_WP) .OR.                       & !
              (THETA1_O1 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,201)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O1 direction to account for detector entrance slit
!
        IF(I_AVER_O1 > 3)  I_AVER_O1 = 3                            !
        IF(I_AVER_O1 < 0)  I_AVER_O1 = 0                            !
        IF(I_AVER_O1 == 0) NDIR_O1   = 1                            !
        IF(I_AVER_O1 == 1) NDIR_O1   = 5                            !
        IF(I_AVER_O1 == 2) NDIR_O1   = 13                           !
        IF(I_AVER_O1 == 3) NDIR_O1   = 49                           !
!
        IF(ITHETA_O2 == 1) THEN                                     !
          IF( (THETA0_O2 < -90.0E0_WP) .OR.                       & !
              (THETA1_O2 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,202)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O2 direction to account for detector entrance slit
!
        IF(I_AVER_O2 > 3)  I_AVER_O2 = 3                            !
        IF(I_AVER_O2 < 0)  I_AVER_O2 = 0                            !
        IF(I_AVER_O2 == 0) NDIR_O2   = 1                            !
        IF(I_AVER_O2 == 1) NDIR_O2   = 5                            !
        IF(I_AVER_O2 == 2) NDIR_O2   = 13                           !
        IF(I_AVER_O2 == 3) NDIR_O2   = 49                           !
!
      ELSE IF(SPECTRO == 'E3E') THEN                                !
!
        IF(IPHI_O1 == 1) WRITE(IUO1,110)                            !
        IF(ITHETA_O1 == 1) WRITE(IUO1,111)                          !
        IF(IE_O1 == 1) WRITE(IUO1,112)                              !
!
        WRITE(IUO1,18) NI,NLI,S_O,ISELRULE,I_SO                     !
        WRITE(IUO1,15) MODE,CS_TYPE                                 !
        WRITE(IUO1,11) IMOD                                         !
!
        IF(ABS(THETA0_IN) > 90.0E0_WP) THEN                         !
          WRITE(IUO1,200)                                           !
          STOP                                                      !
        END IF                                                      !
!
        IF(ITHETA_O1 == 1) THEN                                     !
          IF( (THETA0_O1 < -90.0E0_WP) .OR.                       & !
              (THETA1_O1 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,201)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O1 direction to account for detector entrance slit
!
        IF(I_AVER_O1 > 3)  I_AVER_O1 = 3                            !
        IF(I_AVER_O1 < 0)  I_AVER_O1 = 0                            !
        IF(I_AVER_O1 == 0) NDIR_O1   = 1                            !
        IF(I_AVER_O1 == 1) NDIR_O1   = 5                            !
        IF(I_AVER_O1 == 2) NDIR_O1   = 13                           !
        IF(I_AVER_O1 == 3) NDIR_O1   = 49                           !
!
        IF(ITHETA_O2 == 1) THEN                                     !
          IF( (THETA0_O2 < -90.0E0_WP) .OR.                       & !
              (THETA1_O2 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,202)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O2 direction to account for detector entrance slit
!
        IF(I_AVER_O2 > 3)  I_AVER_O2 = 3                            !
        IF(I_AVER_O2 < 0)  I_AVER_O2 = 0                            !
        IF(I_AVER_O2 == 0) NDIR_O2   = 1                            !
        IF(I_AVER_O2 == 1) NDIR_O2   = 5                            !
        IF(I_AVER_O2 == 2) NDIR_O2   = 13                           !
        IF(I_AVER_O2 == 3) NDIR_O2   = 49                           !
!
        IF(ITHETA_O3 == 1) THEN                                     !
          IF( (THETA0_O3 < -90.0E0_WP) .OR.                       & !
              (THETA1_O3 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,203)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O3 direction to account for detector entrance slit
!
        IF(I_AVER_O3 > 3)  I_AVER_O3 = 3                            !
        IF(I_AVER_O3 < 0)  I_AVER_O3 = 0                            !
        IF(I_AVER_O3 == 0) NDIR_O3   = 1                            !
        IF(I_AVER_O3 == 1) NDIR_O3   = 5                            !
        IF(I_AVER_O3 == 2) NDIR_O3   = 13                           !
        IF(I_AVER_O3 == 3) NDIR_O3   = 49                           !
!
      ELSE IF(SPECTRO == 'ELS') THEN                                !
!
        IF(IPHI_O1 == 1) WRITE(IUO1,120)                            !
        IF(ITHETA_O1 == 1) WRITE(IUO1,121)                          !
        IF(IE_O1 == 1) WRITE(IUO1,122)                              !
!
        WRITE(IUO1,18) NI,NLI,S_O,ISELRULE,I_SO                     !
        WRITE(IUO1,15) MODE,CS_TYPE                                 !
        WRITE(IUO1,11) IMOD                                         !
!
        IF(ABS(THETA0_IN) > 90.0E0_WP) THEN                         !
          WRITE(IUO1,200)                                           !
          STOP                                                      !
        END IF                                                      !
!
        IF(ITHETA_O1 == 1) THEN                                     !
          IF( (THETA0_O1 < -90.0E0_WP) .OR.                       & !
              (THETA1_O1 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,201)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O1 direction to account for detector entrance slit
!
        IF(I_AVER_O1 > 3)  I_AVER_O1 = 3                            !
        IF(I_AVER_O1 < 0)  I_AVER_O1 = 0                            !
        IF(I_AVER_O1 == 0) NDIR_O1   = 1                            !
        IF(I_AVER_O1 == 1) NDIR_O1   = 5                            !
        IF(I_AVER_O1 == 2) NDIR_O1   = 13                           !
        IF(I_AVER_O1 == 3) NDIR_O1   = 49                           !
!
      ELSE IF(SPECTRO == 'RES') THEN                                !
!
        IF(IPHI_O1 == 1) WRITE(IUO1,130)                            !
        IF(ITHETA_O1 == 1) WRITE(IUO1,131)                          !
        IF(IE_O1 == 1) WRITE(IUO1,132)                              !
!
        WRITE(IUO1,12) EDGE,NEDGE,SELRULE_IN,SELRULE_O1             !
        WRITE(IUO1,11) IMOD                                         !
!
      ELSE IF(SPECTRO == 'PLS') THEN                                !
!
        IF(IPHI_O1 == 1) THEN                                       !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,140)                                         !
          ELSE                                                      !
            WRITE(IUO1,141)                                         !
          END IF                                                    !
        END IF                                                      !
        IF(IE_O1 == 1) WRITE(IUO1,142)                              !
        IF(ITHETA_O1 == 1) WRITE(IUO1,143)                          !
!
        WRITE(IUO1,10) NI,NLI,S_O,SELRULE,I_SO                      !
        WRITE(IUO1,16) N_PLA,I_BLK                                  !
        WRITE(IUO1,17) (E_PLA(J_PLA), J_PLA=1,4)                    !
        WRITE(IUO1,11) IMOD                                         !
!
        IF(ABS(THETA0_IN) > 90.0E0_WP) THEN                         !
          WRITE(IUO1,200)                                           !
          STOP                                                      !
        END IF                                                      !
!
        IF(ITHETA_O1 == 1) THEN                                     !
          IF( (THETA0_O1 < -90.0E0_WP) .OR.                       & !
              (THETA1_O1 > 90.0E0_WP) ) THEN                        !
            WRITE(IUO1,201)                                         !
            STOP                                                    !
          END IF                                                    !
        END IF                                                      !
!
!  Averaging over O1 direction to account for detector entrance slit
!
        IF(I_AVER_O1 > 3)  I_AVER_O1 = 3                            !
        IF(I_AVER_O1 < 0)  I_AVER_O1 = 0                            !
        IF(I_AVER_O1 == 0) NDIR_O1   = 1                            !
        IF(I_AVER_O1 == 1) NDIR_O1   = 5                            !
        IF(I_AVER_O1 == 2) NDIR_O1   = 13                           !
        IF(I_AVER_O1 == 3) NDIR_O1   = 49                           !
!
      ELSE IF(SPECTRO == 'STM') THEN                                !
!
        WRITE(IUO1,150)                                             !
!
      ELSE IF(SPECTRO == 'BEM') THEN                                !
!
        WRITE(IUO1,160)                                             !
!
      ELSE IF(SPECTRO == 'EIG') THEN                                !
!
        WRITE(IUO1,170)                                             !
!
      END IF                                                        !
!
!  Formats:
!
  10  FORMAT(11X,I1,A1,8X,A3,4X,A5,8X,I2,9X,'LI,S-O,SELRULE,I_SO')
  11  FORMAT(11X,I2,39X,'IMOD')
  12  FORMAT(11X,A1,I1,5X,A5,29X,'EDGE,SELRULE')
  13  FORMAT(11X,A1,I1,8X,A1,I1,8X,A1,I1,19X,'EDGE_C,EDGE_I,',        &
             'EDGE_A')
  14  FORMAT(12X,I1,8X,I1,A1,I1,28X,'I_MULT,MULT')
  15  FORMAT(10X,A3,6X,A4,29X,'MODE,CS_TYPE')
  16  FORMAT(12X,I1,9X,I1,29X,'N_PLA,I_BLK')
  17  FORMAT(11X,F5.2,5X,F5.2,5X,F5.2,5X,F5.2,6X,'E_PLA(J_PLA)')
  18  FORMAT(11X,I1,A1,8X,A3,8X,I2,8X,I2,9X,'LI,S-O,ISELRULE,I_SO')
!
  50  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL PHOTOELECTRON', &
             ' DIFFRACTION',/)
  51  FORMAT(///,21X,'TYPE OF CALCULATION : FULL HEMISPHERE',         &
             ' PHOTOELECTRON DIFFRACTION',/)
  52  FORMAT(///,21X,'TYPE OF CALCULATION : FINE STRUCTURE ',         &
             'OSCILLATIONS',/)
  53  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR PHOTOELECTRON',     &
             ' DIFFRACTION',/)
!
  60  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL LEED',          &
             ' VARIATIONS',/)
  61  FORMAT(///,21X,'TYPE OF CALCULATION : FULL HEMISPHERE',         &
             ' LEED',/)
  62  FORMAT(///,21X,'TYPE OF CALCULATION : LEED ENERGY ',            &
             'VARIATIONS',/)
  63  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR LEED',              &
             ' VARIATIONS',/)
!
  70  FORMAT(///,21X,'TYPE OF CALCULATION : XAS',/)
!
  80  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL AUGER ELECTRON',&
             ' DIFFRACTION',/)
  81  FORMAT(///,21X,'TYPE OF CALCULATION : FULL HEMISPHERE',         &
             ' AUGER ELECTRON DIFFRACTION',/)
  82  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR AUGER ELECTRON',    &
             ' DIFFRACTION',/)
!
  90  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL AUGER ',        &
             'PHOTOELECTRON COINCIDENCE SPECTROSCOPY,/')
  91  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR AUGER ',            &
             'PHOTOELECTRON COINCIDENCE SPECTROSCOPY',/)
  92  FORMAT(///,9X,'------------------------  FIRST ELECTRON : ',    &
             '------------------------')
  93  FORMAT(///,9X,'------------------------ SECOND ELECTRON : ',    &
             '------------------------')
  94  FORMAT(///,9X,'----------------------------------------------', &
             '----------------------')
!
 100  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL(e,2e) IMPACT ', &
             'SPECTROSCOPY',/)
 101  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR (e,2e) IMPACT ',    &
             'SPECTROSCOPY',/)
 102  FORMAT(///,21X,'TYPE OF CALCULATION : ENERGY (e,2e) IMPACT ',   &
             'SPECTROSCOPY',/)
!
 110  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL(e,3e) IMPACT ', &
             'SPECTROSCOPY',/)
 111  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR (e,3e) IMPACT ',    &
             'SPECTROSCOPY',/)
 112  FORMAT(///,21X,'TYPE OF CALCULATION : ENERGY (e,3e) IMPACT ',   &
             'SPECTROSCOPY',/)
!
 120  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL ELECTRON ',     &
             'ENERGY LOSS SPECTROSCOPY',/)
 121  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR ELECTRON ',         &
             'ENERGY LOSS SPECTROSCOPY',/)
 122  FORMAT(///,21X,'TYPE OF CALCULATION : ENERGY ELECTRON ',        &
             'ENERGY LOSS SPECTROSCOPY',/)
!
 130  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL RESONANT ',     &
             'ELASTIC X-RAY SCATTERING',/)
 131  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR RESONANT ',         &
             'ELASTIC X-RAY SCATTERING',/)
 132  FORMAT(///,21X,'TYPE OF CALCULATION : ENERGY RESONANT ',        &
             'ELASTIC X-RAY SCATTERING',/)
!
 140  FORMAT(///,21X,'TYPE OF CALCULATION : AZIMUTHAL PHOTOEMISSION', &
             ' LOSS SPECTROSCOPY',/)
 141  FORMAT(///,21X,'TYPE OF CALCULATION : FULL HEMISPHERE',         &
             ' PHOTOEMISSION ENERGY LOSS SPECTROSCOPY',/)
 142  FORMAT(///,21X,'TYPE OF CALCULATION : ENERGY PHOTOEMISSION ',   &
             'LOSS SPECTROSCOPY',/)
 143  FORMAT(///,21X,'TYPE OF CALCULATION : POLAR PHOTOEMISSION',     &
             ' LOSS SPECTROSCOPY',/)
!
 150  FORMAT(///,21X,'TYPE OF CALCULATION : SCANNING TUNNELING ',     &
                     'MICROSCOPY,/')
!
 160  FORMAT(///,21X,'TYPE OF CALCULATION : BALLISTIC ENERGY ',       &
                     'ELECTRON MICROSCOPY',/)
!
 170  FORMAT(///,21X,'TYPE OF CALCULATION : EIGENVALUE',              &
             ' ANALYSIS',/)
!
 200  FORMAT(///,2X,' <<<<<<<<<<  THE IN THETA VARIATION EXCEEDS THE ',&
             'PHYSICAL LIMITS (-90,+90)  >>>>>>>>>>',///)
 201  FORMAT(///,2X,' <<<<<<<<<<  THE O1 THETA VARIATION EXCEEDS THE ',&
             'PHYSICAL LIMITS (-90,+90)  >>>>>>>>>>',///)
 202  FORMAT(///,2X,' <<<<<<<<<<  THE O2 THETA VARIATION EXCEEDS THE ',&
             'PHYSICAL LIMITS (-90,+90)  >>>>>>>>>>',///)
 203  FORMAT(///,2X,' <<<<<<<<<<  THE O3 THETA VARIATION EXCEEDS THE ',&
             'PHYSICAL LIMITS (-90,+90)  >>>>>>>>>>',///)
!
      END SUBROUTINE WRITE_SPECTRO_INFO
!
!=======================================================================
!
      SUBROUTINE WRITE_BEAM_INFO(IUO1,BEAM,PARTICLE)
!
!  This subroutine WRITES the information about any beam (electron/photon)
!
!
!  Input variables :
!
!
!                       IUO1      :  Fortran index of the control file
!                       BEAM      :  type of beam used
!                                            --> 'IN'   : incoming
!                                            --> 'EX'   : excited
!                                            --> 'O1'   : 1st outgoing
!                                            --> 'O2'   : 2nd outgoing
!                                            --> 'O3'   : 3rd outgoing
!                       PARTICLE  :  type of particles in the beam
!                                            --> 'ELEC' : electrons
!                                            --> 'PHOT' : photons
!
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 27 May 2021
!
!
      USE AMPLI_IN
      USE AVER_IN
      USE ENERGY_IN
      USE INTEGR_IN
      USE PHI_IN
      USE POLAR_IN
      USE THETA_IN
!
      USE AMPLI_EX
      USE ENERGY_EX
!
      USE AMPLI_O1
      USE AVER_O1
      USE ENERGY_O1
      USE INTEGR_O1
      USE PHI_O1
      USE POLAR_O1
      USE THETA_O1
!
      USE AMPLI_O2
      USE AVER_O2
      USE ENERGY_O2
      USE INTEGR_O2
      USE PHI_O2
      USE POLAR_O2
      USE THETA_O2
!
      USE AMPLI_O3
      USE AVER_O3
      USE ENERGY_O3
      USE INTEGR_O3
      USE PHI_O3
      USE POLAR_O3
      USE THETA_O3
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 38)  ::  TITLE1
      CHARACTER (LEN =  8)  ::  TITLE2
      CHARACTER (LEN =  4)  ::  PARTICLE
      CHARACTER (LEN =  2)  ::  BEAM
!
      INTEGER, INTENT(IN)   ::  IUO1
!
      INTEGER               ::  I_BEAM
      INTEGER               ::  IPHI,ITHETA,IE
      INTEGER               ::  IPOL,I_INT,I_AMP
      INTEGER               ::  I_AVER,ICHKDIR
      INTEGER               ::  NPHI,NTHETA,NE
!
      REAL (WP)             ::  PHI0,PHI1
      REAL (WP)             ::  THETA0,THETA1
      REAL (WP)             ::  E0,E1,ACCEPT
!
!  Initialization of I_BEAM, TITLE1 and TITLE2
!
      IF(BEAM == 'IN') THEN                                         !
        I_BEAM = 1                                                  !
        TITLE1 = 'EXPERIMENTAL PARAMETERS: INCOMING ... '           !
      ELSE IF(BEAM == 'EX') THEN                                    !
        I_BEAM = 2                                                  !
        TITLE1 = 'EXPERIMENTAL PARAMETERS: EXCITED .... '           !
      ELSE IF(BEAM == 'O1') THEN                                    !
        I_BEAM = 3                                                  !
        TITLE1 = 'EXPERIMENTAL PARAMETERS: OUTGOING ... '           !
      ELSE IF(BEAM == 'O2') THEN                                    !
        I_BEAM = 4                                                  !
        TITLE1 = 'EXPERIMENTAL PARAMETERS: 2nd OUTGOING '           !
      ELSE IF(BEAM == 'O3') THEN                                    !
        I_BEAM = 5                                                  !
        TITLE1 = 'EXPERIMENTAL PARAMETERS: 3rd OUTGOING '           !
      END IF                                                        !
!
      IF(PARTICLE == 'ELEC') THEN                                   !
        TITLE2 = 'ELECTRON'                                         !
      ELSE IF(PARTICLE == 'PHOT') THEN                              !
        TITLE2 = 'PHOTON  '                                         !
      END IF                                                        !
!
!  Selection of beam
!
      IF(BEAM == 'IN') THEN                                         !
!
        IPHI    = IPHI_IN                                           !
        ITHETA  = ITHETA_IN                                         !
        IE      = IE_IN                                             !
        IPOL    = IPOL_IN                                           !
        I_INT   = I_INT_IN                                          !
        I_AMP   = I_AMP_IN                                          !
        I_AVER  = I_AVER_IN                                         !
        ICHKDIR = ICHKDIR_IN                                        !
!
        NPHI    = NPHI_IN                                           !
        NTHETA  = NTHETA_IN                                         !
        NE      = NE_IN                                             !
!
        PHI0    = PHI0_IN                                           !
        PHI1    = PHI1_IN                                           !
        THETA0  = THETA0_IN                                         !
        THETA1  = THETA1_IN                                         !
        E0      = E0_IN                                             !
        E1      = E1_IN                                             !
        ACCEPT  = ACCEPT_IN                                         !
!
      ELSE IF(BEAM == 'EX') THEN                                    !
!
        IE      = IE_EX                                             !
!
        NE      = NE_EX                                             !
!
        E0      = E0_EX                                             !
        E1      = E1_EX                                             !
!
      ELSE IF(BEAM == 'O1') THEN                                    !
!
        IPHI    = IPHI_O1                                           !
        ITHETA  = ITHETA_O1                                         !
        IE      = IE_O1                                             !
        IPOL    = IPOL_O1                                           !
        I_INT   = I_INT_O1                                          !
        I_AMP   = I_AMP_O1                                          !
        I_AVER  = I_AVER_O1                                         !
        ICHKDIR = ICHKDIR_O1                                        !
!
        NPHI    = NPHI_O1                                           !
        NTHETA  = NTHETA_O1                                         !
        NE      = NE_O1                                             !
!
        PHI0    = PHI0_O1                                           !
        PHI1    = PHI1_O1                                           !
        THETA0  = THETA0_O1                                         !
        THETA1  = THETA1_O1                                         !
        E0      = E0_O1                                             !
        E1      = E1_O1                                             !
        ACCEPT  = ACCEPT_O1                                         !
!
      ELSE IF(BEAM == 'O2') THEN                                    !
!
        IPHI    = IPHI_O2                                           !
        ITHETA  = ITHETA_O2                                         !
        IE      = IE_O2                                             !
        IPOL    = IPOL_O2                                           !
        I_INT   = I_INT_O2                                          !
        I_AMP   = I_AMP_O2                                          !
        I_AVER  = I_AVER_O2                                         !
        ICHKDIR = ICHKDIR_O2                                        !
!
        NPHI    = NPHI_O2                                           !
        NTHETA  = NTHETA_O2                                         !
        NE      = NE_O2                                             !
!
        PHI0    = PHI0_O2                                           !
        PHI1    = PHI1_O2                                           !
        THETA0  = THETA0_O2                                         !
        THETA1  = THETA1_O2                                         !
        E0      = E0_O2                                             !
        E1      = E1_O2                                             !
        ACCEPT  = ACCEPT_O2                                         !
!
      ELSE IF(BEAM == 'O3') THEN                                    !
!
        IPHI    = IPHI_O3                                           !
        ITHETA  = ITHETA_O3                                         !
        IE      = IE_O3                                             !
        IPOL    = IPOL_O3                                           !
        I_INT   = I_INT_O3                                          !
        I_AMP   = I_AMP_O3                                          !
        I_AVER  = I_AVER_O3                                         !
        ICHKDIR = ICHKDIR_O3                                        !
!
        NPHI    = NPHI_O3                                           !
        NTHETA  = NTHETA_O3                                         !
        NE      = NE_O3                                             !
!
        PHI0    = PHI0_O3                                           !
        PHI1    = PHI1_O3                                           !
        THETA0  = THETA0_O3                                         !
        THETA1  = THETA1_O3                                         !
        E0      = E0_O3                                             !
        E1      = E1_O3                                             !
        ACCEPT  = ACCEPT_O3                                         !
!
      END IF                                                        !
!
      WRITE(IUO1,5) TITLE1,TITLE2                                   !
!
!  Electron case:
!
      IF(PARTICLE == 'ELEC') THEN                                   !
!
        IF(I_BEAM /= 2) WRITE(IUO1,10) IPHI,NPHI,PHI0,PHI1          !
        IF(I_BEAM /= 2) WRITE(IUO1,11) ITHETA,NTHETA,THETA0,THETA1  !
        WRITE(IUO1,12) IE,NE,E0,E1                                  !
        IF(I_BEAM /= 2) WRITE(IUO1,13) I_INT,I_AMP                  !
        IF(I_BEAM /= 2) WRITE(IUO1,14) I_AVER,ACCEPT,ICHKDIR        !
!
!  Photon case:
!
      ELSE IF(PARTICLE == 'PHOT') THEN                              !
!
        WRITE(IUO1,20) IPHI,NPHI,PHI0,PHI1                          !
        WRITE(IUO1,21) ITHETA,NTHETA,THETA0,THETA1                  !
        WRITE(IUO1,22) IE,NE,E0,E1                                  !
        WRITE(IUO1,23) IPOL                                         !
!
      END IF                                                        !
!
!  Formats:
!
   5  FORMAT(///,21X,A38,A8,/)
!
  10  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IPHI,NPHI,PHI0,PHI1')
  11  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'ITHETA,NTHETA,THETA0,',  &
             'THETA1')
  12  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IE,NE,E0,E1')
  13  FORMAT(11X,I2,8X,I2,29X,'I_INT,I_AMP')
  14  FORMAT(12X,I1,8X,F5.2,6X,I1,19X,'I_AVER,ACCEPT,ICHKDIR')
!
  20  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IPHI,NPHI,PHI0,PHI1')
  21  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'ITHETA,NTHETA,THETA0,',  &
             'THETA1')
  22  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IE,NE,E0,E1')
  23  FORMAT(11X,I2,39X,'IPOL')
!
      END SUBROUTINE WRITE_BEAM_INFO
!
!=======================================================================
!
      SUBROUTINE WRITE_INFO_START(IUO1,J_EL,IPOTC,IPRINT,INV)
!
!  This routine writes the information on the spectroscopy selected into the check file
!     at the beginning of the calculation
!
!  Input variables :
!
!                       IUO1      :  Fortran unit for the check file *.lis
!                       J_EL      :  electron number
!                       IPOTC     :  switch indicating whether the potential is complex or not
!                       IPRINT    :  switch to control the printing
!                       INV       :  switch indicating whether matrix inversion is used or not
!                                    ---> 0 : series expansion
!                                    ---> 1 : matrix inversion
!                                    ---> 2 : correlation expansion
!                                    ---> 3 : renormalized expansion
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 27 May 2021
!
!
      USE ATOMIC_INDEX
      USE AUGER_MULTIPLET
      USE AUGER_PARAM
      USE CURRENT_BEAM
      USE CURRENT_CALC
      USE CURRENT_INIT_VAL
      USE CURRENT_SPIN
      USE DICHROISM
      USE EXP_TYPE
      USE INIT_M
      USE OPTICAL_ELTS
!
      USE ENERGY_IN
      USE PHI_IN
      USE THETA_IN
!
      USE ENERGY_O1
      USE PHI_O1
      USE THETA_O1
!
      USE ENERGY_O2
      USE PHI_O2
      USE THETA_O2
!
      USE ENERGY_O3
      USE PHI_O3
      USE THETA_O3
!
      IMPLICIT NONE
!
!
      INTEGER, INTENT(IN)        ::  IUO1,J_EL,IPOTC,IPRINT,INV
!
!----------  Photoelectron diffraction case (PED)  ----------
!
      IF((SPECTRO == 'PED') .OR. (SPECTRO == 'APC')) THEN           !
        IF(SPECTRO == 'PED') THEN                                   !
          IF(IPHI == 1) THEN                                        !
            IF(STEREO == ' NO') THEN                                !
              WRITE(IUO1,10)                                        !
            ELSE                                                    !
              WRITE(IUO1,11)                                        !
            END IF                                                  !
          END IF                                                    !
          IF(ITHETA == 1) WRITE(IUO1,12)                            !
          IF(IE == 1) WRITE(IUO1,13)                                !
          IF(I_TEST == 1) WRITE(IUO1,200)                           !
        END IF                                                      !
        IF((SPECTRO == 'APC') .AND. (J_EL == 1)) THEN               !
          WRITE(IUO1,401)                                           !
          WRITE(IUO1,201)                                           !
        END IF                                                      !
        IF(J_EL == 2) GO TO 888                                     !
        IF(IPRINT > 0) THEN                                         !
          WRITE(IUO1,202)                                           !
          WRITE(IUO1,203)                                           !
          IF(ISPIN == 0) THEN                                       !
            WRITE(IUO1,204)                                         !
          ELSE                                                      !
            WRITE(IUO1,206)                                         !
          END IF                                                    !
          WRITE(IUO1,203)                                           !
          IF(IPOTC == 0) THEN                                       !
            WRITE(IUO1,207)                                         !
          ELSE                                                      !
            WRITE(IUO1,208)                                         !
          END IF                                                    !
          WRITE(IUO1,203)                                           !
          IF(I_C1 == 0) THEN                                        !
            WRITE(IUO1,209)                                         !
            WRITE(IUO1,203)                                         !
            IF(IPOL == 0) THEN                                      !
              WRITE(IUO1,210)                                       !
            ELSE IF(ABS(IPOL) == 1) THEN                            !
              WRITE(IUO1,211)                                       !
            ELSE IF(IPOL == 2) THEN                                 !
              WRITE(IUO1,212)                                       !
            END IF                                                  !
            WRITE(IUO1,203)                                         !
            IF(I_DICHR > 0) THEN                                    !
              WRITE(IUO1,213)                                       !
            END IF                                                  !
            WRITE(IUO1,203)                                         !
            WRITE(IUO1,202)                                         !
            WRITE(IUO1,302)                                         !
            WRITE(IUO1,303) THETA0_IN,PHI0_IN                       !
            IF((SPECTRO == 'PED') .AND. (IMOD == 1)) THEN           !
              WRITE(IUO1,304)                                       !
            END IF                                                  !
          ELSE                                                      !
            WRITE(IUO1,214)                                         !
            WRITE(IUO1,203)                                         !
            WRITE(IUO1,215)                                         !
            WRITE(IUO1,203)                                         !
            WRITE(IUO1,202)                                         !
            WRITE(IUO1,302)                                         !
            WRITE(IUO1,303) THETA0_IN,PHI0_IN                       !
            IF((SPECTRO == 'PED') .AND. (IMOD == 1)) THEN           !
              WRITE(IUO1,304)                                       !
            END IF                                                  !
          END IF                                                    !
!
          IF(I_AT == 0) THEN                                        !
            IF(INV == 0) THEN                                       !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,219) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,220) N_SCAT                              !
              END IF                                                !
            ELSE IF(INV == 1) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,221)                                     !
              ELSE                                                  !
                WRITE(IUO1,222)                                     !
              END IF                                                !
            ELSE IF(INV == 2) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,223) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,224) N_SCAT                              !
              END IF                                                !
            ELSE IF(INV == 3) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,225) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,226) N_SCAT                              !
              END IF                                                !
            END IF                                                  !
          ELSE                                                      !
            IF(I_BASIS == 0) THEN                                   !
              WRITE(IUO1,217)                                       !
            ELSE                                                    !
              WRITE(IUO1,218)                                       !
            END IF                                                  !
          END IF                                                    !
!
       END IF                                                       !
 888   CONTINUE                                                     !
      END IF                                                        !
!
!----------  LEED case (LED)  ----------
!
      IF(SPECTRO == 'LED') THEN                                     !
        IF(IPHI == 1) THEN                                          !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,20)                                          !
          ELSE                                                      !
            WRITE(IUO1,21)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(ITHETA == 1) WRITE(IUO1,22)                              !
        IF(IE == 1) WRITE(IUO1,23)                                  !
        IF(IPRINT > 0) THEN                                         !
          WRITE(IUO1,202)                                           !
          WRITE(IUO1,203)                                           !
          IF(ISPIN == 0) THEN                                       !
            WRITE(IUO1,204)                                         !
          ELSE                                                      !
            WRITE(IUO1,206)                                         !
          END IF                                                    !
          WRITE(IUO1,203)                                           !
          IF(IPOTC == 0) THEN                                       !
            WRITE(IUO1,207)                                         !
          ELSE                                                      !
            WRITE(IUO1,208)                                         !
          END IF                                                    !
          WRITE(IUO1,203)                                           !
          WRITE(IUO1,202)                                           !
          WRITE(IUO1,300)                                           !
          WRITE(IUO1,301) THETA0_IN,PHI0_IN                         !
          IF((SPECTRO == 'LED') .AND. (IMOD == 1)) THEN             !
            WRITE(IUO1,304)                                         !
          END IF                                                    !
!
          IF(I_AT == 0) THEN                                        !
            IF(INV == 0) THEN                                       !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,219) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,220) N_SCAT                              !
              END IF                                                !
            ELSE IF(INV == 1) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,221)                                     !
              ELSE                                                  !
                WRITE(IUO1,222)                                     !
              END IF                                                !
            ELSE IF(INV == 2) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,223) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,224) N_SCAT                              !
              END IF                                                !
            ELSE IF(INV == 3) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,225) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,226) N_SCAT                              !
              END IF                                                !
            END IF                                                  !
          ELSE                                                      !
            IF(I_BASIS == 0) THEN                                   !
              WRITE(IUO1,217)                                       !
            ELSE                                                    !
              WRITE(IUO1,218)                                       !
            END IF                                                  !
          END IF                                                    !
!
        END IF                                                      !
      END IF                                                        !
!
!----------  Auger diffraction case (AED)  ----------
!
      IF((SPECTRO == 'AED') .OR. (SPECTRO == 'APC')) THEN           !
        IF(SPECTRO == 'AED') THEN                                   !
          IF(IPHI == 1) THEN                                        !
            IF(STEREO == ' NO') THEN                                !
              WRITE(IUO1,30)                                        !
            ELSE                                                    !
              WRITE(IUO1,31)                                        !
            END IF                                                  !
          END IF                                                    !
          IF(ITHETA == 1) WRITE(IUO1,32)                            !
          IF(I_TEST == 1) WRITE(IUO1,200)                           !
        END IF                                                      !
        IF((SPECTRO == 'APC') .AND. (J_EL == 2)) THEN               !
          WRITE(IUO1,402)                                           !
          WRITE(IUO1,201)                                           !
        END IF                                                      !
        IF((SPECTRO == 'AED') .OR. (J_EL == 2)) THEN                !
          IF(IPRINT > 0) THEN                                       !
            WRITE(IUO1,202)                                         !
            WRITE(IUO1,203)                                         !
            IF(ISPIN == 0) THEN                                     !
              WRITE(IUO1,204)                                       !
            ELSE                                                    !
              WRITE(IUO1,206)                                       !
            END IF                                                  !
            WRITE(IUO1,203)                                         !
            IF(IPOTC == 0) THEN                                     !
              WRITE(IUO1,207)                                       !
            ELSE                                                    !
              WRITE(IUO1,208)                                       !
            END IF                                                  !
            WRITE(IUO1,203)                                         !
            WRITE(IUO1,202)                                         !
            WRITE(IUO1,33) AUGER                                    !
            CALL AUGER_MULT                                         !
            IF(I_MULT == 0) THEN                                    !
              WRITE(IUO1,34)                                        !
            ELSE                                                    !
              WRITE(IUO1,35) MULTIPLET                              !
            END IF                                                  !
!
          IF(I_AT == 0) THEN                                        !
            IF(INV == 0) THEN                                       !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,219) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,220) N_SCAT                              !
              END IF                                                !
            ELSE IF(INV == 1) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,221)                                     !
              ELSE                                                  !
                WRITE(IUO1,222)                                     !
              END IF                                                !
            ELSE IF(INV == 2) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,223) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,224) N_SCAT                              !
              END IF                                                !
            ELSE IF(INV == 3) THEN                                  !
              IF(I_BASIS == 0) THEN                                 !
                WRITE(IUO1,225) N_SCAT                              !
              ELSE                                                  !
                WRITE(IUO1,226) N_SCAT                              !
              END IF                                                !
            END IF                                                  !
          ELSE                                                      !
            IF(I_BASIS == 0) THEN                                   !
              WRITE(IUO1,217)                                       !
            ELSE                                                    !
              WRITE(IUO1,218)                                       !
            END IF                                                  !
          END IF                                                    !
!
          END IF                                                    !
        END IF                                                      !
      END IF                                                        !
!
!----------  APECS case
!
      IF((SPECTRO == 'APC') .AND. (J_EL == 1)) THEN                 !
        IF(IPHI == 1) THEN                                          !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,40)                                          !
          ELSE                                                      !
            WRITE(IUO1,41)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(ITHETA == 1) WRITE(IUO1,42)                              !
        IF(I_TEST == 1) WRITE(IUO1,200)                             !
      END IF                                                        !
!
!----------  X-ray absorption case (XAS)  ----------
!
      IF(SPECTRO == 'XAS') THEN                                     !
        WRITE(IUO1,50)                                              !
      END IF                                                        !
!
!  Formats
!
!  PED
!
  10  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE AZIMUTHAL PHOTOELECTRON DIFFRACTION CALCULATION #####',  &
      '#####',/////)
  11  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE FULL ANGLE PHOTOELECTRON DIFFRACTION CALCULATION ',      &
      '##########',/////)
  12  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE POLAR PHOTOELECTRON DIFFRACTION CALCULATION #####',      &
      '#####',/////)
  13  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE FINE STRUCTURE OSCILLATIONS CALCULATION #####',          &
      '#####',/////)
!
!  LEED
!
  20  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE AZIMUTHAL LEED CALCULATION #####'                        &
      '#####',/////)
  21  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE FULL ANGLE LEED CALCULATION ',                           &
      '##########',/////)
  22  FORMAT(/////,6X,'########## BEGINNING ',                         &
      'OF THE POLAR LEED CALCULATION #####',                           &
      '#####',/////)
  23  FORMAT(/////,5X,'########## BEGINNING ',                         &
      'OF THE ENERGY LEED CALCULATION #####',                          &
      '#####',/////)
!
!  AED
!
  30  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE AZIMUTHAL AUGER DIFFRACTION CALCULATION #####',          &
      '#####',/////)
  31  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE FULL ANGLE AUGER DIFFRACTION CALCULATION ',              &
      '##########',/////)
  32  FORMAT(/////,6X,'########## BEGINNING ',                         &
      'OF THE POLAR AUGER DIFFRACTION CALCULATION #####',              &
      '#####',/////)
  33  FORMAT(////,31X,'AUGER LINE :',A6,//)
  34  FORMAT(///,20X,'CALCULATION MADE FOR THE FULL AUGER LINE',       &
             ' ',/,' ',/,' ')
  35  FORMAT(///,20X,'CALCULATION MADE FOR THE ',A3,' MULTIPLET ',     &
            'LINE',' ',/,' ',/,' ')
!
!  APECS
!
  40  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE AZIMUTHAL APECS DIFFRACTION CALCULATION #####',          &
      '#####',/////)
  41  FORMAT(/////,'########## BEGINNING ',                            &
      'OF THE FULL ANGLE APECS DIFFRACTION CALCULATION ',              &
      '##########',/////)
  42  FORMAT(/////,6X,'########## BEGINNING ',                         &
      'OF THE POLAR APECS DIFFRACTION CALCULATION #####',              &
      '#####',/////)
!
!  XAS
!
  50  FORMAT(/////,'##########              BEGINNING ',               &
      'OF THE EXAFS CALCULATION              ##########',/////)
!
!  General information
!
 200  FORMAT('  ----->  TEST CALCULATION : NO EXCITATION ',            &
                 'MATRIX ELEMENTS TAKEN INTO ACCOUNT  <-----',///)
 201  FORMAT('            ',/)
 202  FORMAT(24X,'+++++++++++++++++++++++++++++++++++++')
 203  FORMAT(24X,'+',35X,'+')
 204  FORMAT(24X,'+              STANDARD             +')
 205  FORMAT(24X,'+            RELATIVISTIC           +')
 206  FORMAT(24X,'+           SPIN-POLARIZED          +')
 207  FORMAT(24X,'+    REAL POTENTIAL CALCULATION     +')
 208  FORMAT(24X,'+   COMPLEX POTENTIAL CALCULATION   +')
 209  FORMAT(24X,'+                WITH               +')
 210  FORMAT(24X,'+        NON POLARIZED LIGHT        +')
 211  FORMAT(24X,'+     LINEARLY POLARIZED LIGHT      +')
 212  FORMAT(24X,'+    CIRCULARLY POLARIZED LIGHT     +')
 213  FORMAT(24X,'+         IN DICHROIC MODE          +')
 214  FORMAT(24X,'+                FOR                +')
 215  FORMAT(24X,'+          A SINGLE CHANNEL         +')
 217  FORMAT(//,26X,'(PLANE WAVES - ATOMIC CASE)')
 218  FORMAT(//,24X,'(SPHERICAL WAVES - ATOMIC CASE)')
 219  FORMAT(///,19X,'(PLANE WAVES MULTIPLE SCATTERING - ORDER ',      &
             ')')
 220  FORMAT(///,17X,'(SPHERICAL WAVES MULTIPLE SCATTERING - ',        &
                     'ORDER ',I1,')')
 221  FORMAT(///,17X,'(PLANE WAVES MULTIPLE SCATTERING - MATRIX ',     &
             'INVERSION)')
 222  FORMAT(///,15X,'(SPHERICAL WAVES MULTIPLE SCATTERING - MATRIX ', &
             'INVERSION)')
 223  FORMAT(///,19X,'(PLANE WAVES CORRELATION SCATTERING - ',         &
                     'ORDER ',I1, ')')
 224  FORMAT(///,17X,'(SPHERICAL WAVES CORRELATION SCATTERING - ',     &
                     'ORDER ',I1,')')
 225  FORMAT(///,19X,'(PLANE WAVES RENORMALIZED SCATTERING - ',        &
                     'ORDER ',I1,')')
 226  FORMAT(///,17X,'(SPHERICAL WAVES RENORMALIZED SCATTERING - ',    &
                     'ORDER ', I1,')')
!
!  Position of beams: incoming/outgoing(s), light/electron(s)
!
 300  FORMAT(////,31X,'POSITION OF THE INITIAL BEAM :',/)
 301  FORMAT(14X,'TH_BEAM = ',F6.2,' DEGREES',5X,'PHI_BEAM = ',        &
             F6.2,' DEGREES')
 302  FORMAT(////,31X,'POSITION OF THE LIGHT :',/)
 303  FORMAT(14X,'TH_LIGHT = ',F6.2,' DEGREES',5X,'PHI_LIGHT = ',      &
             F6.2,' DEGREES')
 304  FORMAT(14X,' (WHEN THE DETECTOR IS ALONG ',                      &
             'THE NORMAL TO THE SURFACE)')
!
 400  FORMAT(///,9X,'------------------------   INCOMING ELECTRON : ', &
             '------------------------')
 401  FORMAT(///,9X,'------------------------  FIRST OUT ELECTRON : ', &
             '------------------------')
 402  FORMAT(///,9X,'------------------------ SECOND OUT ELECTRON : ', &
             '------------------------')
 403  FORMAT(///,9X,'------------------------  THIRD OUT ELECTRON : ', &
             '------------------------')
!
      END SUBROUTINE WRITE_INFO_START
!
!=======================================================================
!
      SUBROUTINE WRITE_INFO_END(IUO1)
!
!  This routine writes the information on the spectroscopy selected into the check file
!     at the end of the calculation
!
!  Input variables :
!
!                       IUO1      :  Fortran unit for the check file *.lis
!
!
!
!   Author :  D. Sébilleau
!
!
!                                         Last modified : 27 May 2021
!
      USE CURRENT_CALC
      USE EXP_TYPE
!
      IMPLICIT  NONE
!
      INTEGER, INTENT(IN)        ::  IUO1
!
      IF(SPECTRO == 'PED') THEN                                     !
!
        IF(IPHI == 1) THEN                                          !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,10)                                          !
          ELSE                                                      !
            WRITE(IUO1,11)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(ITHETA == 1) WRITE(IUO1,12)                              !
        IF(IE == 1)     WRITE(IUO1,13)                              !
!
      ELSE IF(SPECTRO == 'LED') THEN                                !
!
        IF(IPHI == 1) THEN                                          !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,20)                                          !
          ELSE                                                      !
            WRITE(IUO1,21)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(ITHETA == 1) WRITE(IUO1,22)                              !
        IF(IE == 1)     WRITE(IUO1,23)                              !
!
      ELSE IF(SPECTRO == 'XAS') THEN                                !
!
        WRITE(IUO1,51)                                              !
!
      ELSE IF(SPECTRO == 'AED') THEN                                !
!
        IF(IPHI == 1) THEN                                          !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,30)                                          !
          ELSE                                                      !
            WRITE(IUO1,31)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(ITHETA == 1) WRITE(IUO1,32)                              !
!
      ELSE IF(SPECTRO == 'APC') THEN                                !
!
        IF(IPHI == 1) THEN                                          !
          IF(STEREO == ' NO') THEN                                  !
            WRITE(IUO1,40)                                          !
          ELSE                                                      !
            WRITE(IUO1,41)                                          !
          END IF                                                    !
        END IF                                                      !
        IF(ITHETA == 1) WRITE(IUO1,42)                              !
!
      END IF                                                        !
!
!  Formats!
!
!
!  PED
!
  10  FORMAT(/////,'########## END OF THE ',                      &
      'AZIMUTHAL PHOTOELECTRON DIFFRACTION CALCULATION #####',    &
      '#####')
  11  FORMAT(/////,'########## END OF THE ',                      &
      'FULL ANGLE PHOTOELECTRON DIFFRACTION CALCULATION #####',   &
      '#####')
  12  FORMAT(/////,'########## END OF THE ',                      &
      'POLAR PHOTOELECTRON DIFFRACTION CALCULATION ##########')
  13  FORMAT(/////,'########## END OF THE ',                      &
      'FINE STRUCTURE OSCILLATIONS CALCULATION #####',            &
      '#####')
!
!  LEED
!
  20  FORMAT(/////,'########## END ',                             &
      'OF THE AZIMUTHAL LEED CALCULATION #####',                  &
      '#####',/////)
  21  FORMAT(/////,'########## END OF THE ',                      &
      'FULL ANGLE LEED CALCULATION #####',                        &
      '#####')
  22  FORMAT(/////,6X,'########## END ',                          &
      'OF THE POLAR LEED CALCULATION #####',                      &
      '#####',/////)
  23  FORMAT(/////,5X,'########## END ',                          &
      'OF THE ENERGY LEED CALCULATION #####',                     &
      '#####',/////)
!
!  AED
!
  30  FORMAT(/////,'########## END ',                             &
      'OF THE AZIMUTHAL AUGER DIFFRACTION CALCULATION #####',     &
      '#####',/////)
  31  FORMAT(/////,'########## END ',                             &
      'OF THE FULL ANGLE AUGER DIFFRACTION CALCULATION #####',    &
      '#####',/////)
  32  FORMAT(/////,6X,'########## END ',                          &
      'OF THE POLAR AUGER DIFFRACTION CALCULATION #####',         &
      '#####',/////)
!
!  APECS
!
  40  FORMAT(/////,'########## END ',                             &
      'OF THE AZIMUTHAL APECS DIFFRACTION CALCULATION #####',     &
      '#####',/////)
  41  FORMAT(/////,'########## END ',                             &
      'OF THE FULL ANGLE APECS DIFFRACTION CALCULATION #####',    &
      '#####',/////)
  42  FORMAT(/////,6X,'########## END ',                          &
      'OF THE POLAR APECS DIFFRACTION CALCULATION #####',         &
      '#####',/////)
!
!  XAS
!
  51  FORMAT(/////,'##########                END OF THE ',       &
      'EXAFS CALCULATION                ##########')
!
      END SUBROUTINE WRITE_INFO_END
!
!=======================================================================
!
      SUBROUTINE HEADERS(IUO2)
!
!   This subroutine writes headers containing the main parameters
!       of the calculation in the result file. The number of
!       lines written depends of the spectroscopy
!
!
!   Author: Didier Sébilleau
!
!                                        Last modified : 26 May 2021
!
!
      USE ALGORITHMS
!
      USE ENERGY_IN
      USE PHI_IN
      USE THETA_IN
!
      USE ENERGY_O1
      USE PHI_O1
      USE THETA_O1
!
      USE ENERGY_O2
      USE PHI_O2
      USE THETA_O2
!
      USE ENERGY_O3
      USE PHI_O3
      USE THETA_O3
!
      USE AUGER_PARAM
      USE CLUSTER
      USE CURRENT_BEAM
      USE CURRENT_FINA_VAL
      USE CURRENT_INIT_VAL
      USE EXP_TYPE
      USE HEADER
      USE INDAT
      USE INFILES
      USE INIT_J
      USE INIT_L
      USE INIT_M
      USE XAS
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 25)  ::  DA_FILE,TL_FILE
      CHARACTER (LEN = 25)  ::  RD_FILE,CL_FILE
!
      INTEGER, INTENT(IN)   ::  IUO2
!
      INTEGER, PARAMETER    ::  MAX_LENGTH = 50
!
      INTEGER               ::  J_CHAR
      INTEGER               ::  N_CHAR1,N_CHAR2,N_CHAR3,N_CHAR4
!
      WRITE(IUO2,1)                                                 !
      WRITE(IUO2,2)                                                 !
!
!   Input files section:
!
!      Finding the filename
!
      IF(SPECTRO == 'PED') THEN                                     !
        DA_FILE = INDATA(1)                                         !
        CL_FILE = INFILE1                                           !
        TL_FILE = INFILE9                                           !
        RD_FILE = INFILE10                                          !
      ELSE IF(SPECTRO == 'XAS') THEN                                !
        CL_FILE = INDATA(1)                                         !
        CL_FILE = INFILE1                                           !
        TL_FILE = INFILE9                                           !
        RD_FILE = INFILE10                                          !
      ELSE IF(SPECTRO == 'LED') THEN                                !
        CL_FILE = INDATA(1)                                         !
        CL_FILE = INFILE1                                           !
        TL_FILE = INFILE9                                           !
        RD_FILE = INFILE10                                          !
      ELSE IF(SPECTRO == 'AED') THEN                                !
        CL_FILE = INDATA(1)                                         !
        CL_FILE = INFILE1                                           !
        TL_FILE = INFILE9                                           !
        RD_FILE = INFILE10                                          !
      ELSE IF(SPECTRO == 'APC') THEN                                !
        CL_FILE = INDATA(1)                                         !
        CL_FILE = INFILE1                                           !
        TL_FILE = INFILE9                                           !
        RD_FILE = INFILE10                                          !
      ELSE IF(SPECTRO == 'RES') THEN                                !
        CL_FILE = INDATA(1)                                         !
        CL_FILE = INFILE1                                           !
        TL_FILE = INFILE9                                           !
        RD_FILE = INFILE10                                          !
      ELSE IF(SPECTRO == 'ELS') THEN                                !
        CL_FILE = INDATA(1)                                         !
        CL_FILE = INFILE1                                           !
        TL_FILE = INFILE9                                           !
        RD_FILE = INFILE10                                          !
      ELSE IF(SPECTRO == 'EIG') THEN                                !
        CL_FILE = INDATA(1)                                         !
        CL_FILE = INFILE1                                           !
        TL_FILE = INFILE9                                           !
        RD_FILE = INFILE10                                          !
      END IF                                                        !
!
!      Checking the size of filenames
!
      N_CHAR1 = 0                                                   !
      DO J_CHAR = 1 , MAX_LENGTH                                    !
        IF(DA_FILE(J_CHAR:J_CHAR) == ' ') GO TO 500                 !
        N_CHAR1 = N_CHAR1 + 1                                       !
      END DO                                                        !
 500  CONTINUE                                                      !
!
      N_CHAR2 = 0                                                   !
      DO J_CHAR = 1, MAX_LENGTH                                     !
        IF(CL_FILE(J_CHAR:J_CHAR) == ' ') GO TO 501                 !
        N_CHAR2 = N_CHAR2 + 1                                       !
      END DO                                                        !
 501  CONTINUE                                                      !
!
      N_CHAR3 = 0                                                   !
      DO J_CHAR = 1 , MAX_LENGTH                                    !
        IF(RD_FILE(J_CHAR:J_CHAR) == ' ') GO TO 502                 !
        N_CHAR3 = N_CHAR3 + 1                                       !
      END DO                                                        !
 502  CONTINUE                                                      !
!
      N_CHAR4 = 0                                                   !
      DO J_CHAR = 1 , MAX_LENGTH                                    !
        IF(TL_FILE(J_CHAR:J_CHAR) == ' ') GO TO 503                 !
        N_CHAR4 = N_CHAR4 + 1                                       !
      END DO                                                        !
 503  CONTINUE                                                      !
!
      WRITE(IUO2,3) DA_FILE(6:N_CHAR1)                              !
      WRITE(IUO2,4) TL_FILE(4:N_CHAR4)                              !
      IF(SPECTRO /= 'LED') THEN                                     !
        WRITE(IUO2,5) RD_FILE(5:N_CHAR3)                            !
      END IF                                                        !
      WRITE(IUO2,6) CL_FILE(6:N_CHAR2)                              !
      WRITE(IUO2,2)                                                 !
!
!   Type of calculation
!
      WRITE(IUO2,2)                                                 !
!
      IF(SPECTRO == 'PED') THEN                                     !
        WRITE(IUO2,11) SPECTRO,ALGO_O1                              !
        IF(ALGO_O1 == 'SE') THEN                                    !
          WRITE(IUO2,12) NO,N_SCAT,IFWD,IPW,ILENGTH                 !
        ELSE IF(ALGO_O1 == 'CE') THEN                               !
          WRITE(IUO2,13) N_SCAT                                     !
        END IF                                                      !
        WRITE(IUO2,14) VINT                                         !
      ELSE IF(SPECTRO == 'XAS') THEN                                !
        WRITE(IUO2,11) SPECTRO,ALGO_IN                              !
        IF(ALGO_IN == 'SE') THEN                                    !
          WRITE(IUO2,12) NO,N_SCAT,IFWD,IPW,ILENGTH                 !
        ELSE IF(ALGO_IN == 'CE') THEN                               !
          WRITE(IUO2,13) N_SCAT                                     !
        END IF                                                      !
        WRITE(IUO2,14) VINT                                         !
      ELSE IF(SPECTRO == 'LED') THEN                                !
        WRITE(IUO2,11) SPECTRO,ALGO_O1                              !
        IF(ALGO_O1 == 'SE') THEN                                    !
          WRITE(IUO2,12) NO,N_SCAT,IFWD,IPW,ILENGTH                 !
        ELSE IF(ALGO_O1 == 'CE') THEN                               !
          WRITE(IUO2,13) N_SCAT                                     !
        END IF                                                      !
        WRITE(IUO2,14) VINT                                         !
      ELSE IF(SPECTRO == 'AED') THEN                                !
        WRITE(IUO2,11) SPECTRO,ALGO_O1                              !
        IF(ALGO_O1 == 'SE') THEN                                    !
          WRITE(IUO2,12) NO,N_SCAT,IFWD,IPW,ILENGTH                 !
        ELSE IF(ALGO_O1 == 'CE') THEN                               !
          WRITE(IUO2,13) N_SCAT                                     !
        END IF                                                      !
        WRITE(IUO2,14) VINT                                         !
      ELSE IF(SPECTRO == 'APC') THEN                                !
        WRITE(IUO2,15) SPECTRO,ALGO_O1,ALGO_O2                      !
        WRITE(IUO2,14) VINT                                         !
      ELSE IF(SPECTRO == 'EIG') THEN                                !
        WRITE(IUO2,11) SPECTRO,ALGO_O1                              !
      ELSE IF(SPECTRO == 'RES') THEN                                !
        CONTINUE                                                    !
      ELSE IF(SPECTRO == 'ELS') THEN                                !
        CONTINUE                                                    !
      END IF                                                        !
!
      WRITE(IUO2,2)                                                 !
!
!   Initial state parameters
!
      IF(SPECTRO == 'PED') THEN                                     !
        WRITE(IUO2,21) NI,NLI,S_O,SELRULE                           !
      ELSE IF(SPECTRO == 'XAS') THEN                                !
        WRITE(IUO2,22) EDGE,NEDGE,SELRULE                           !
      ELSE IF(SPECTRO == 'LED') THEN                                !
        CONTINUE                                                    !
      ELSE IF(SPECTRO == 'AED') THEN                                !
        WRITE(IUO2,24) AUGER,MULTIPLET                              !
      ELSE IF(SPECTRO == 'APC') THEN                                !
        WRITE(IUO2,21) NI,NLI,S_O,SELRULE                           !
        WRITE(IUO2,24) AUGER,MULTIPLET                              !
      ELSE IF(SPECTRO == 'RES') THEN                                !
        CONTINUE                                                    !
      ELSE IF(SPECTRO == 'ELS') THEN                                !
        CONTINUE                                                    !
      END IF                                                        !
!
      WRITE(IUO2,2)                                                 !
!
!   Angular and energy parameters
!
      IF(SPECTRO == 'PED') THEN                                     !
        WRITE(IUO2,35)                                              !
        WRITE(IUO2,34) THETA0_IN,PHI0_IN,E0_IN                      !
        WRITE(IUO2,2)                                               !
        WRITE(IUO2,36) 1                                            !
        WRITE(IUO2,31) THETA0_O1,THETA1_O1                          !
        WRITE(IUO2,32) PHI0_O1,PHI1_O1                              !
        WRITE(IUO2,33) E0_O1,E1_O1                                  !
      ELSE IF(SPECTRO == 'XAS') THEN                                !
        WRITE(IUO2,35)                                              !
        WRITE(IUO2,33) EK_INI,EK_FIN                                !
        WRITE(IUO2,34) THETA0_IN,PHI0_IN,E0_IN                      !
      ELSE IF(SPECTRO == 'LED') THEN                                !
        WRITE(IUO2,35)                                              !
        WRITE(IUO2,31) THETA0_IN,PHI0_IN                            !
        WRITE(IUO2,2)                                               !
        WRITE(IUO2,36) 1                                            !
        WRITE(IUO2,31) THETA0_O1,THETA1_O1                          !
        WRITE(IUO2,32) PHI0_O1,PHI1_O1                              !
        WRITE(IUO2,2)                                               !
        WRITE(IUO2,33) E0,E1                                        !
      ELSE IF(SPECTRO == 'AED') THEN                                !
        WRITE(IUO2,36) 1                                            !
        WRITE(IUO2,31) THETA0_O1,THETA1_O1                          !
        WRITE(IUO2,32) PHI0_O1,PHI1_O1                              !
      ELSE IF(SPECTRO == 'APC') THEN                                !
        WRITE(IUO2,35)                                              !
        WRITE(IUO2,34) THLUM,PHILUM,ELUM                            !
        WRITE(IUO2,2)                                               !
        WRITE(IUO2,37)                                              !
        WRITE(IUO2,31) THETA0_O1,THETA1_O1                          !
        WRITE(IUO2,32) PHI0_O1,PHI1_O1                              !
        WRITE(IUO2,33) E0_O1,E1_O1                                  !
        WRITE(IUO2,2)                                               !
        WRITE(IUO2,38)                                              !
        WRITE(IUO2,31) THETA0_O2,THETA1_O2                          !
        WRITE(IUO2,32) PHI0_O2,PHI1_O2                              !
      ELSE IF(SPECTRO == 'EIG') THEN                                !
        WRITE(IUO2,33) EK_INI,EK_FIN                                !
      ELSE IF(SPECTRO == 'RES') THEN                                !
        CONTINUE                                                    !
      ELSE IF(SPECTRO == 'ELS') THEN                                !
        CONTINUE                                                    !
      END IF                                                        !
!
!   End of headers
!
      WRITE(IUO2,2)                                                 !
      WRITE(IUO2,1)                                                 !
      WRITE(IUO2,39)                                                !
!
!   Formats
!
   1  FORMAT('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', &
             '!!!!!!!!!!!!!!!!')
   2  FORMAT('!',69X,'!')
   3  FORMAT('!',10X,'data file        :  ',A19,20X,'!')
   4  FORMAT('!',10X,'t-matrix file    :    ',A17,20X,'!')
   5  FORMAT('!',10X,'rad integral file: ',A20,20X,'!')
   6  FORMAT('!',10X,'cluster file     :  ',A19,20X,'!')
!
  11  FORMAT('!',10X,'spectroscopy     :      ',A3,8X,'algorithm : ',   &
             A2,10X,'!')
  12  FORMAT('!',15X,'NO = ',I1,'  N_SCAT = ',I2,'  IFWD = ',I1,        &
            '  IPW = ',I1,'  ILENGTH = ',I1,5X,'!')
  13  FORMAT('!',15X,'N_SCAT = ',I2,45X,'!')
  14  FORMAT('!',10X,'inner potential  :    ',F6.2,' eV',28X,'!')
  15  FORMAT('!',10X,'spectroscopy: ',A3,10X,'algorithm (photo): ',A2,  &
              11X,'!',/,'!',37X,'algorithm (auger): ',A2, 11X,'!')
!
  21  FORMAT('!',10X,'initial state    : ',I1,A1,1X,A3,                 &
             ' selection rules:',' SELRULE = ',A5,' !')
  22  FORMAT('!',10X,'initial state    : ',A1,I1,2X,' selection rules:',&
             ' SELRULE = ',A5,4X,'!')
  24  FORMAT('!',10X,'initial state    : ',A6,2X,' multiplet: ',A3,     &
             17X,'!')
!
  31  FORMAT('!',10X,'THETA_INI: ',F8.2,6X,'THETA_FIN: ',F8.2,15X,'!')
  32  FORMAT('!',10X,'PHI_INI  : ',F8.2,6X,'PHI_FIN  : ',F8.2,15X,'!')
  33  FORMAT('!',10X,'E_INI    : ',F8.2,' eV',3X,'E_FIN    : ',F8.2,    &
             ' eV',12X,'!')
  34  FORMAT('!',10X,'THETA_LUM: ',F8.2,2X,'PHI_LUM: ',F8.2,2X,         &
             'E_LUM: ',F8.2,' eV !')
  35  FORMAT('!',10X,'incoming beam    : ',40X,'!')
  36  FORMAT('!',10X,'outgoing beam ',I1,'  : ',40X,'!')
  37  FORMAT('!',10X,'photoelectron beam:',40X,'!')
  38  FORMAT('!',10X,'auger beam        :',40X,'!')
  39  FORMAT(71X)
!
      END SUBROUTINE HEADERS
!
!=======================================================================
!
      SUBROUTINE STOP_EXT(SPECTRO)
!
!  This routine stops the code when the dimension N_TILT_M in the
!      spec.inc file is insufficient for the number of values to
!      Gaussian average over (as generated by the ext_dir.f code)
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 27 May 2021
!
!
      USE INFILES
      USE INUNITS
      USE OUTUNITS
!
      USE EXTERNAL_O1
      USE EXTERNAL_O2
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  SPECTRO
!
      INTEGER               ::  NSET_O1,NSET_O2
!
      NSET_O1 = 1                                                   !
      NSET_O2 = 1                                                   !
!
      IF((SPECTRO == 'PED').OR.(SPECTRO == 'AED')) THEN             !
        IF(I_EXT_O1 == -1) THEN                                     !
          OPEN(UNIT=IUI11, FILE = INFILE11, STATUS='OLD')           !
          READ(IUI11,15) NSET_O1                                    !
          CLOSE(IUI11)                                              !
        END IF                                                      !
        IF(I_EXT_O2 == -1) THEN                                     !
          OPEN(UNIT=IUI11, FILE = INFILE11, STATUS='OLD')           !
          READ(IUI11,15) NSET_O2                                    !
          CLOSE(IUI11)                                              !
        END IF                                                      !
      END IF                                                        !
      IF(SPECTRO == 'APC') THEN                                     !
        IF(I_EXT_O1 == -1) THEN                                     !
          OPEN(UNIT=IUI11, FILE = INFILE11, STATUS='OLD')           !
          READ(IUI11,15) NSET_O1                                    !
          CLOSE(IUI11)                                              !
        END IF                                                      !
        IF(I_EXT_O2 == -1) THEN                                     !
          OPEN(UNIT=IUI9, FILE = INFILE9, STATUS='OLD')             !
          READ(IUI9,15) NSET_O2                                     !
          CLOSE(IUI9)                                               !
        END IF                                                      !
      END IF                                                        !
!
      IF(MAX(NSET_O1,NSET_O2) > N_TILT_M) THEN                      !
        WRITE(IUO1,10) MAX(NSET_O2,NSET_O1)                         !
        STOP                                                        !
      END IF                                                        !
!
!  Formats:
!
  10  FORMAT(///,16X,'<<<<<<<<<<  N_TILT_M SHOULD BE AT LEAST ', &
             I3,'  >>>>>>>>>>')
  15  FORMAT(8X,I3)
!
      END SUBROUTINE STOP_EXT
!
END MODULE INPUT_DATA
