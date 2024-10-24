!
!=======================================================================
!
MODULE CALC_PED_MS
!
!  This module perform the calculation of the MS cross-section of PED
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE PED_MS(NPLAN,ZEM,IPHA,N_INFILES,J_FILE,NP)
!
!   This subroutine computes the PED formula in the spin-independent case
!        from a non spin-orbit resolved initial core state LI.
!
!   Alternatively, it can compute the PED amplitude for the APECS process.
!
!   The calculation is performed using a series expansion for the
!           expression of the scattering path operator
!
!
!   Authors: Didier Sébilleau and Xu Junqing
!
!
!                                        Last modified : 17 Jun 2021
!
!
      USE REAL_NUMBERS,        ONLY : ZERO,ONE,TWO,EIGHT,HALF
      USE COMPLEX_NUMBERS,     ONLY : ZEROC,IC
      USE PI_ETC,              ONLY : PI
      USE CONSTANTS_P1,        ONLY : BOHR
      USE CONSTANTS_P2,        ONLY : ALPHA
      USE ENE_CHANGE,          ONLY : RYD
!
      USE ABSORBER
      USE ANA_DIR
      USE AVERAGING
      USE BEAMS_ALGO
      USE CLUSTER
      USE CLUSTER_LIMITS,      ONLY : VAL
      USE CURRENT_AVER
      USE CURRENT_CALC
      USE CLUSTER_COORD
      USE CURRENT_COEF_RENORM, ONLY : C_REN
      USE CURRENT_BEAM,        ONLY : N_SCAT
      USE CURRENT_EXT_DIR
      USE CURRENT_FINA_VAL
      USE CURRENT_FIXSCAN
      USE CURRENT_INIT_VAL,    ONLY : PHI0,THETA0,E0
      USE CURRENT_RENORM,      ONLY : I_REN
      USE CURRENT_SPIN
      USE CURRENT_T_MATRIX
      USE DEB_WAL_CLU
      USE DICHROISM
      USE EXP_TYPE
      USE EXTREMES
      USE INFILES
      USE INUNITS
      USE INIT_L
      USE INIT_J
      USE M_RULE
      USE OUTUNITS
      USE TESTS
!
      USE ENERGY_IN          !  \
      USE PHI_IN             !   > incoming photon beam
      USE THETA_IN           !  /
!
      USE ALGORITHMS
      USE ARC_SIN
      USE BEAM_AMPLITUDE
      USE CALC_TAU0J_CE
      USE CALC_TAU0J_MI
      USE CALC_TAU0J_RE
      USE CALC_TAU0J_SE
      USE CALC_TAU0J_SM

      USE DIRECTION_ANA
      USE ELECTRON_CHOICE
      USE ELECTRON_PHOTON
      USE HEADER
      USE INPUT_DATA
      USE INITIALIZE_CALC
      USE MEAN_FREE_PATH
      USE RENORMALIZATION
      USE TREATMENT_PED
      USE VECTOR
      USE VIBRATIONS
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)   ::  TAU_TYPE
      CHARACTER (LEN =  7)   ::  STAT
      CHARACTER (LEN = 13)   ::  OUTDATA1
      CHARACTER (LEN = 24)   ::  OUTFILE
!
      INTEGER, INTENT(IN)    ::  NPLAN,IPHA,N_INFILES,J_FILE
      INTEGER, INTENT(OUT)   ::  NP
!
      INTEGER                ::  JEMET,JNEM,JEL
      INTEGER                ::  I_CALL,IS
      INTEGER                ::  JPLAN,JE,JAT,LAT
      INTEGER                ::  JATLEM,J_FIXED,J_SCAN
      INTEGER                ::  J_SET,J,K,JEPS,JLINE
      INTEGER                ::  ISAUT,I_DUM1
      INTEGER                ::  N_DUM1,N_DUM2,N_DUM3
      INTEGER                ::  JREF,JSET,JDIR,JPES,JJ
      INTEGER                ::  MI,LF,MF,ILF,INDF
      INTEGER                ::  IOUT,IOUT2,I_MI
      INTEGER                ::  N
!
      INTEGER                ::  INT
!
      REAL (WP), INTENT(IN)  ::  ZEM
!
      REAL (WP)              ::  LUM(3),POL(3),PV(3)
      REAL (WP)              ::  DIRLUM(3),DIRPOL(3,2)
      REAL (WP)              ::  E_PH(NE_M)
      REAL (WP)              ::  PIS180,SMALL
      REAL (WP)              ::  ZSURF,ZSURFE,Z_PLANE
      REAL (WP)              ::  ECIN,CFM,Q,XMFP1,GAMMA
      REAL (WP)              ::  XK2UJ2,FIX_STEP,XINCRF,XINCRS
      REAL (WP)              ::  DPHI,RPHI,DTHETA,RTHETA
      REAL (WP)              ::  RTH,RPH,THD,PED
      REAL (WP)              ::  SSETDIR_1,SSETDIF_1
      REAL (WP)              ::  SSETDIR_2,SSETDIF_2
      REAL (WP)              ::  SSET2DIR_1,SSET2DIF_1
      REAL (WP)              ::  SSET2DIR_2,SSET2DIF_2
      REAL (WP)              ::  DTHETAP,DPHIP
      REAL (WP)              ::  SEPSDIR_1,SEPSDIF_1
      REAL (WP)              ::  SEPSDIR_2,SEPSDIF_2
      REAL (WP)              ::  SRDIF_1,SRDIR_1,SRDIF_2,SRDIR_2
      REAL (WP)              ::  SMIDIR_1,SMIDIF_1,SMIDIR_2,SMIDIF_2
      REAL (WP)              ::  ANALUM,DZZEM,ELUM,THLUM,PHILUM,THETALUM
      REAL (WP)              ::  VKR,W,XMFP
!
      REAL (WP)              ::  FLOAT,SIN,COS,ACOS
!
      COMPLEX (WP), ALLOCATABLE  ::  TAU(:,:,:)
!
      COMPLEX (WP)           ::  COEF
      COMPLEX (WP)           ::  TLT(0:NT_M,4,NATM,NE_M)
      COMPLEX (WP)           ::  AL(LINMAX),BL(LINMAX)
      COMPLEX (WP)           ::  SLF_1,SLF_2,SJDIR_1,SJDIR_2
      COMPLEX (WP)           ::  MLFLI_E1(2,-LI_M:LI_M,-LI_M-1:LI_M+1,0:LI_M+1,3)
      COMPLEX (WP)           ::  MLFLI_E2(2,-LI_M:LI_M,-LI_M-2:LI_M+2,0:LI_M+2,3,3)
      COMPLEX (WP)           ::  MLFLI_E3(2,-LI_M:LI_M,-LI_M-3:LI_M+3,0:LI_M+3,3,3)
      COMPLEX (WP)           ::  MLFLI_M1(2,-LI_M:LI_M,-LI_M-1:LI_M+1,0:LI_M,3)
      COMPLEX (WP)           ::  MLFLI_M2(2,-LI_M:LI_M,-LI_M-2:LI_M+2,0:LI_M+1,3,3)
      COMPLEX (WP)           ::  MLIL0(2,-LI_M:LI_M,LINFMAX)
!
!
      DATA PIS180     / 0.01745329251994329576923690768488612713E0_WP /
      DATA SMALL      / 0.0001E0_WP /
!
      TAU_TYPE = 'OUT'                                              !
!
      JEMET = NTYPEM                                                !
      JNEM  = NTYPEM                                                !
!
      NSET     = 1                                                  !
      JEL      = 3                                                  !
      OUTDATA1 = 'CROSS-SECTION'                                    !
!
!  Computation of the renormalization coeffcients
!
      IF(ALGO(3) == 'RE') THEN                                      !
        CALL COEF_RENORM(N_SCAT)                                    !
      END IF                                                        !
!
!  Parameters of the light
!
      ELUM   = E0_IN                                                !
      THLUM  = THETA0_IN                                            !
      PHILUM = PHI0_IN                                              !
!
!  Writing the beginning of the amplitude file (whenever necessary)
!
      I_CALL = 1                                                    !
      CALL WRITE_PED_AMPLITUDE(J_FILE,1,1,1,1,1,1,1,1,ZEROC,ZEROC,& !
                               ZEROC,ZEROC,I_CALL,IOUT,IOUT2,I_MI,& !
                               STAT,OUTFILE)                        !
!
!  Writing spectroscopy-related information in the check file unit IUO1
!
      CALL WRITE_INFO_START(IUO1,1,IPOTC,IPRINT,0)                  !
!
!  Position of the light when the analyzer is along the z axis :
!                  (X_LUM_Z,Y_LUM_Z,Z_LUM_Z)
!
      IS = 1                                                        !
      CALL INIT_INC_BEAM(IS,THLUM,PHILUM,THETA0,PHI0,IMOD,I_EXT,  & !
                         IPH_1,RTH,RPH,DIRLUM,LUM)                  !
!
      ZSURF = VAL(1)                                                !
!
!  Opening the cross-section result file
!
      IF((ISOM == 0) .OR. (J_FILE == 1)) THEN                       !
        OPEN(UNIT=IOUT, FILE=OUTFILE, STATUS=STAT)                  !
      END IF                                                        !
!
!  Writing the headers in the output file
!
      CALL HEADERS(IOUT)                                            !
!
      IF((ISOM == 0) .OR. ((ISOM > 0) .AND. (J_FILE == 1))) THEN    !
        WRITE(IOUT,12) SPECTRO,OUTDATA1                             !
        WRITE(IOUT,9) ISPIN,I_DICHR,I_SO,ISFLIP,ICHKDIR,IPHI,     & !
                      ITHETA,IE,IPH_1,I_EXT                         !
      END IF                                                        !
!
!  Writing the second part of the beginning of the amplitude file (whenever necessary)
!
      I_CALL= 2                                                     !
      CALL WRITE_PED_AMPLITUDE(J_FILE,1,1,1,1,1,1,1,1,ZEROC,ZEROC,& !
                               ZEROC,ZEROC,I_CALL,IOUT,IOUT2,I_MI,& !
                               STAT,OUTFILE)                        !
!
      IF(ISOM == 0) THEN                                            !
        WRITE(IOUT,79) NPLAN,NEMET,NTHETA,NPHI,NE                   !
      ELSE IF((ISOM /= 0) .AND. (J_FILE == 1)) THEN                 !
        WRITE(IOUT,11) NTHETA,NPHI,NE                               !
      END IF                                                        !
!
!   Loop over the planes
!
      DO JPLAN = 1, NPLAN                                           !
        Z_PLANE = VAL(JPLAN)                                        !
!
        IF((IPHA == 1) .OR. (IPHA == 2)) THEN                       !
          DZZEM = ABS(Z_PLANE - ZEM)                                !
          IF(DZZEM < SMALL) GO TO 10                                ! <-- absorber on plane JPLAN
          GO TO 1                                                   ! <-- absorber not on plane JPLAN
        END IF                                                      !
!
  10   CONTINUE                                                     !
!
       NTYPEM = 1                                                   !
       IF((ISORT1 == 0) .AND. (IPRINT > 0)) THEN                    !
         IF(I_TEST /=  2) THEN                                      !
          WRITE(IUO1,52) JPLAN,EMIT(1),EMIT(1),EMIT(1),NTYPEM       !
         END IF                                                     !
       END IF                                                       !
       IF(ISOM == 1) N = JPLAN                                      !
       ZSURFE = VAL(1) - EMIT(3)                                    !
!
!   Loop over the energies
!
        DO JE = 1, NE                                               !
         FMIN(0) = ONE                                              !
         FMAX(0) = ONE                                              !
         IF(NE > 1) THEN                                            !
           ECIN     = E0   + FLOAT(JE-1) * (E1 - E0) / FLOAT(NE-1)  !
           E_PH(JE) = ELUM + FLOAT(JE-1) * (E1 - E0) / FLOAT(NE-1)  !
         ELSE IF(NE == 1) THEN                                      !
           ECIN     = E0                                            !
           E_PH(JE) = ELUM                                          !
         END IF                                                     !
!
!   Coefficient outside the outer sums
!
         IF(I_TEST /=  1) THEN                                      !
           CFM = EIGHT * PI * E_PH(JE) * ALPHA                      !
         ELSE                                                       !
           CFM = ONE                                                !
         END IF                                                     !
!
!   Photon wave number in atomic units
!
         Q = E_PH(JE) * ALPHA * HALF / RYD                          !
!
!   Inelastic mean free path
!
         CALL MFP(ECIN,XMFP,*6)                                     !
         XMFP1 = XMFP / A                                           !
         IF(IPRINT > 0) WRITE(IUO1,56) A,XMFP1                      !
!
         GAMMA = ONE / (TWO * XMFP1)                                !
         IF(IPOTC == 0) THEN                                        !
           VK(JE) = VK(JE) + IC * GAMMA                             !
         END IF                                                     !
         IF(I_TEST /=  1) THEN                                      !
           VKR = REAL(VK(JE),KIND=WP)                               !
         ELSE                                                       !
           VKR = ONE                                                !
         END IF                                                     !
         IF(I_MI == 1) THEN                                         !
           WRITE(IOUT2,21) ECIN,VKR*CFM                             !
         END IF                                                     !
!
!  Averaged T-matrix elements taking into account vibrational effects
!
         IF(IDWSPH == 1) THEN                                       !
           IF(IMSD >= 1) WRITE(IUO1,22)                             !
           DO JAT = 1, N_PROT                                       !
             IF(IMSD == 0) THEN                                     !
               XK2UJ2 = VK2(JE) * UJ2(JAT)                          !
             ELSE                                                   !
               XK2UJ2 = VK2(JE) * UJ_SQ(JAT)                        !
               WRITE(IUO1,23) JAT,UJ_SQ(JAT)*A*A                    !
             END IF                                                 !
             CALL AV_T_MATRIX(JAT,JE,XK2UJ2,TLT)                    !
             DO LAT = 0, LMAX(JAT,JE)                               !
               TL(LAT,1,JAT,JE ) = TLT(LAT,1,JAT,JE)                !
             END DO                                                 !
           END DO                                                   !
         END IF                                                     !
!
         IF(ABS(I_EXT) >= 1) THEN                                   !
           OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')            !
           READ(IUI11,13) IDIR,NSET,N_DUM1                          !
           READ(IUI11,14) I_DUM1,N_DUM2,N_DUM3                      !
         END IF                                                     !
!
!   Storage of the coupling matrix elements MLFLI along the basis
!                       directions X,Y and Z
!
         CALL CALC_MLFLI(JE,E_PH(JE),JEL,Q,MLFLI_E1,MLFLI_E2,     & !
                         MLFLI_E3,MLFLI_M1,MLFLI_M2)                !
!
!   Calculation of the scattering path operator TAU
!
         IF(ALLOCATED(TAU)) THEN                                    !
           DEALLOCATE(TAU)                                          !
         END IF                                                     !
         ALLOCATE(TAU(LINMAX,LINFMAX,NATCLU_M))                     !
!
         JATLEM = 1                                                 !
!
         IF(ALGO(3) == 'CE') THEN                                   !
           CALL CALC_TAU_CE(TAU_TYPE,JATLEM,JE,TAU)                 !
         ELSE IF(ALGO(3) == 'MI') THEN                              !
           CALL CALC_TAU_MI(TAU_TYPE,JATLEM,JE,TAU)                 !
         ELSE IF(ALGO(3) == 'RE') THEN                              !
           CALL CALC_TAU_RE(TAU_TYPE,JATLEM,JE,TAU)                 !
         ELSE IF(ALGO(3) == 'SE') THEN                              !
           CALL CALC_TAU_SE(TAU_TYPE,JATLEM,JE,TAU)                 !
         ELSE IF(ALGO(3) == 'SM') THEN                              !
           CALL CALC_TAU_SM(TAU_TYPE,JATLEM,JE,TAU)                 !
         END IF                                                     !
!
!   Calculation of the photoelectron diffraction formula
!
!
!    Loop over the 'fixed' angle (i.e. the angle with the lower number of values)
!
         DO J_FIXED = 1, N_FIXED                                    !
          IF(N_FIXED > 1) THEN                                      !
            IF(I_EXT == 0) THEN                                     !
              FIX_STEP = (FIX1 - FIX0) / FLOAT(N_FIXED - 1)         !
              XINCRF   = FLOAT(J_FIXED - 1) * FIX_STEP              !
            ELSE                                                    !
              XINCRF   = ZERO                                       !
            END IF                                                  !
          ELSE IF(N_FIXED == 1) THEN                                !
            XINCRF     = ZERO                                       !
          END IF                                                    !
!
          IF(ABS(I_EXT) >= 1) THEN                                  !
            READ(IUI11,86) JSET,JLINE,THD,PED                       !
            IF(I_EXT == -1) BACKSPACE IUI11                         !
            THETA0 = THD                                            !
            PHI0   = PED                                            !
          END IF                                                    !
          IF(IPH_1 == 1) THEN                                       !
             IF(I_EXT == 0) THEN                                    !
               DPHI = PHI0 + XINCRF                                 !
             ELSE                                                   !
               DPHI = PED                                           !
             END IF                                                 !
             RPHI = DPHI * PIS180                                   !
             IF(IPRINT > 0) WRITE(IUO1,66) DPHI                     !
          ELSE                                                      !
             ISAUT = 0                                              !
             IF(I_EXT == 0) THEN                                    !
               DTHETA = THETA0 + XINCRF                             !
             ELSE                                                   !
               DTHETA = THD                                         !
             END IF                                                 !
             RTHETA = DTHETA * PIS180                               !
             IF(ABS(DTHETA) > 90.0E0_WP) ISAUT = ISAUT + 1          !
             IF(I_EXT >= 1)              ISAUT = 0                  !
             IF(I_TEST == 2)             ISAUT = 0                  !
             IF(ISAUT > 0) GO TO 8                                  !
             IF(IPRINT > 0) WRITE(IUO1,65) DTHETA                   !
             IF((IPRINT >  0) .AND. (I_TEST /=  2)) WRITE(IUO1,59)  !
             IF((IPRINT == 1) .AND. (I_TEST /=  2)) WRITE(IUO1,60)  !
!
!  THETA-dependent number of PHI points for stereographic
!     representation (to obtain a uniform sampling density).
!     (Courtesy of J. Osterwalder - University of Zurich)
!
             IF(STEREO == 'YES') THEN                               !
               N_SCAN = INT( (SCAN1 - SCAN0) *                    & !
                              SIN(RTHETA) / FIX_STEP + SMALL ) + 1  !
             END IF                                                 !
!
          END IF                                                    !
!
          IF((N_FIXED > 1) .AND. (IMOD == 1)) THEN                  !
!
!  When there are several sets of scans (N_FIXED > 1),
!   the initial position LUM of the light is recalculated
!      for each initial position (RTH,RPH) of the analyzer
!
            IS = 2                                                  !
            CALL INIT_INC_BEAM(IS,DTHETA,DPHI,THETA0,PHI0,IMOD,   & !
                               I_EXT,IPH_1,RTH,RPH,DIRLUM,LUM)      !
!
          END IF                                                    !
!
!    Loop over the scanned angle (i.e. the angle with the higher number of values)
!
          DO J_SCAN = 1, N_SCAN                                     !
           IF(N_SCAN > 1) THEN                                      !
             XINCRS = FLOAT(J_SCAN - 1) * (SCAN1 - SCAN0) /       & !
                      FLOAT(N_SCAN - 1)                             !
           ELSE IF(N_SCAN == 1) THEN                                !
             XINCRS = ZERO                                          !
           END IF                                                   !
           IF(I_EXT == -1) THEN                                     !
             READ(IUI11,86) JSET,JLINE,THD,PED                      !
             BACKSPACE IUI11                                        !
           END IF                                                   !
           IF(IPH_1 == 1) THEN                                      !
             ISAUT = 0                                              !
             IF(I_EXT == 0) THEN                                    !
               DTHETA = THETA0 + XINCRS                             !
             ELSE                                                   !
               DTHETA = THD                                         !
             END IF
             RTHETA = DTHETA * PIS180                               !
             IF(ABS(DTHETA) > 90.0E0_WP) ISAUT = ISAUT + 1          !
             IF(I_EXT >= 1)              ISAUT = 0                  !
             IF(I_TEST == 2)             ISAUT = 0                  !
             IF(ISAUT > 0) GO TO 8                                  !
             IF(IPRINT > 0) WRITE(IUO1,65) DTHETA                   !
             IF((IPRINT >  0) .AND. (I_TEST /=  2)) WRITE(IUO1,59)  !
             IF((IPRINT == 1) .AND. (I_TEST /=  2)) WRITE(IUO1,60)  !
           ELSE                                                     !
             IF(I_EXT == 0) THEN                                    !
               DPHI = PHI0 + XINCRS                                 !
             ELSE                                                   !
               DPHI = PED                                           !
             END IF                                                 !
             RPHI = DPHI * PIS180                                   !
             IF(IPRINT > 0) WRITE(IUO1,66) DPHI                     !
           END IF                                                   !
!
!  Loop over the sets of directions to average over (for Gaussian average)
!
!
           SSETDIR_1 = ZERO                                         !
           SSETDIF_1 = ZERO                                         !
           SSETDIR_2 = ZERO                                         !
           SSETDIF_2 = ZERO                                         !
!
           SSET2DIR_1 = ZERO                                        !
           SSET2DIF_1 = ZERO                                        !
           SSET2DIR_2 = ZERO                                        !
           SSET2DIF_2 = ZERO                                        !
!
           IF(I_EXT == -1) THEN                                     !
             JREF = INT(NSET) / 2 + 1                               !
           ELSE                                                     !
             JREF = 1                                               !
           END IF                                                   !
!
           DO J_SET = 1, NSET                                       !
           IF(I_EXT == -1) THEN                                     !
             READ(IUI11,86) JSET,JLINE,THD,PED,W                    !
             DTHETA = THD                                           !
             DPHI   = PED                                           !
             RTHETA = DTHETA * PIS180                               !
             RPHI   = DPHI   * PIS180                               !
!
!  Here, there are several sets of scans (NSET > 1), so
!   the initial position LUM of the light must be
!   recalculated for each initial position of the analyzer
!
             IS = 3                                                 !
             CALL INIT_INC_BEAM(IS,DTHETA,DPHI,TH_0(JSET),        & !
                                PH_0(JSET),IMOD,I_EXT,IPH_1,      & !
                                RTH,RPH,DIRLUM,LUM)                 !
           ELSE                                                     !
             W = ONE                                                !
           END IF                                                   !
!
           IF(I_EXT == -1) PRINT 89                                 !
!
           CALL DIR_ANA(VINT,ECIN,RTHETA,RPHI)                      !
!
           IF(J_SET == JREF) THEN                                   !
             DTHETAP = DTHETA                                       !
             DPHIP   = DPHI                                         !
           END IF                                                   !
!
           IF(I_EXT == -1) THEN                                     !
             WRITE(IUO1,88) DTHETA,DPHI                             !
           END IF                                                   !
!
!     ..........          Case IMOD=1 only          ..........
!
!  Calculation of the position of the light when the analyzer is at
!   (THETA,PHI). DIRLUM is the direction of the light and its initial
!   value (at (THETA0,PHI0)) is LUM. AXE is the direction of the theta
!   rotation axis and EPS is defined so that (AXE,DIRLUM,EPS) is a
!   direct orthonormal basis. The transform of a vector R by a rotation
!   of OMEGA about AXE is then given by
!
!     R' = R COS(OMEGA) + (AXE.R)(1-COS(OMEGA)) AXE + (AXE^R) SIN(OMEGA)
!
!   Here, DIRANA is the internal direction of the analyzer and ANADIR
!                      its external position
!
!   Note that when the initial position of the analyzer is (RTH,RPH)
!    which coincides with (RTH0,RPH0) only for the first fixed angle
!
           IF(IMOD == 1) THEN                                       !
             CALL ROTATE_INC_BEAM(RTHETA,RPHI,RTH,RPH,ITHETA,IPHI,& !
                                  J_SCAN,LUM,DIRLUM)                !
           END IF                                                   !
!
           IF(DIRLUM(3) >   ONE) DIRLUM(3) =   ONE                  !
           IF(DIRLUM(3) < - ONE) DIRLUM(3) = - ONE                  !
           THETALUM = ACOS(DIRLUM(3))                               !
           IF(I_TEST == 2) THETALUM = - THETALUM                    !
           COEF = DIRLUM(1) + IC * DIRLUM(2)                        !
!
           CALL ARCSIN(COEF,DIRLUM(3),PHILUM)                       !
           ANALUM = ANADIR(1,1) * DIRLUM(1) +                     & !
                    ANADIR(2,1) * DIRLUM(2) +                     & !
                    ANADIR(3,1) * DIRLUM(3)                         !
!
           SEPSDIR_1 = ZERO                                         !
           SEPSDIF_1 = ZERO                                         !
           SEPSDIR_2 = ZERO                                         !
           SEPSDIF_2 = ZERO                                         !
!
!   Loop over the directions of polarization
!
           DO JEPS = 1, NEPS                                        !
            IF((JEPS == 1) .AND. (IPOL >= 0)) THEN                  !
              DIRPOL(1,JEPS) =   COS(THETALUM) * COS(PHILUM)        !
              DIRPOL(2,JEPS) =   COS(THETALUM) * SIN(PHILUM)        !
              DIRPOL(3,JEPS) = - SIN(THETALUM)                      !
            ELSE                                                    !
              DIRPOL(1,JEPS) = - SIN(PHILUM)                        !
              DIRPOL(2,JEPS) =   COS(PHILUM)                        !
              DIRPOL(3,JEPS) =   ZERO                               !
            END IF                                                  !
            IF(ABS(IPOL) == 1) THEN                                 !
              IF(IPRINT > 0) THEN                                   !
                WRITE(IUO1,61) (DIRANA(J,1),J=1,3),               & !
                               (DIRLUM(K),  K=1,3),               & !
                               (DIRPOL(K,1),K=1,3),               & !
                                ANALUM                              !
              END IF                                                !
            ELSE
              IF((JEPS == 1) .AND. (IPRINT > 0)) THEN               !
                WRITE(IUO1,63) (DIRANA(J,1),J=1,3),               & !
                               (DIRLUM(K),  K=1,3),ANALUM           !
              END IF                                                !
            END IF                                                  !
            IF((JEPS == 1) .AND. (I_EXT == -1)) PRINT 89            !
!
!  Vector product LIGHT x POLARIZATION for magnetic dipole
!
            DO JJ = 1, 3                                            !
              POL(JJ) = DIRPOL(JJ,JEPS)                             !
            END DO                                                  !
            CALL PRVECT(DIRLUM,POL,PV,ONE)                          !
!
!   Calculation of the coupling matrix MLIL0 (dipole, quadrupole and octupole terms)
!

            CALL CALC_MLIL0(JEPS,I_TEST,DIRPOL,DIRLUM,MLFLI_E1,   & !
                            MLFLI_E2,MLFLI_E3,MLFLI_M1,MLFLI_M2,  & !
                            MLIL0)                                  !
!
            SRDIF_1 = ZERO                                          !
            SRDIR_1 = ZERO                                          !
            SRDIF_2 = ZERO                                          !
            SRDIR_2 = ZERO                                          !

!
!   Loop over the different directions of the analyzer contained in a cone
!                         (for angular averaging)
!
            DO JDIR = 1, NDIR                                       !
!
!   Computation of the multiple scattering amplitudes BL
!
             CALL BL_AMP(TAU_TYPE,ALGO_O1,JE,ZSURF,NATYP,NTYPEM,  & !
                         NTYPEM,JATLEM,JDIR,ZSURFE,GAMMA,TAU,AL,BL) !
!
             SMIDIR_1 = ZERO                                        !
             SMIDIF_1 = ZERO                                        !
             SMIDIR_2 = ZERO                                        !
             SMIDIF_2 = ZERO                                        !
!
!   Loop over the equiprobable azimuthal quantum numbers MI corresponding
!                        to the initial core state LI
!
             DO MI = - LI, LI                                       !
              SJDIR_1 = ZEROC                                       !
              SJDIR_2 = ZEROC                                       !
!
!   Calculation of the atomic emission (used a a reference for the output)
!
              DO LF = LF1, LF2, ISTEP_LF                            !
                ILF = LF * LF + LF + 1                              !
                DO MF = - LF, LF                                    !
                  INDF = ILF + MF                                   !
                  IF((MF < MI-M_R) .OR. (MF > MI+M_R)) GO TO 444    !
                  IF( (SELRULE == '00000') .AND.                  & !
                      (MF /=  MI) ) GO TO 444                       !
!
                  IF(I_REN >= 1) THEN                               !
                    AL(INDF) = C_REN(0) * AL(INDF)                  !
                  END IF                                            !
!
                  SJDIR_1 = SJDIR_1 + AL(INDF) * MLIL0(1,MI,INDF)   !
                  IF(I_DICHR >= 1) THEN                             !
                    SJDIR_2 = SJDIR_2 + AL(INDF) * MLIL0(2,MI,INDF) !
                  END IF                                            !
!
 444              CONTINUE                                          !
                END DO                                              !
              END DO                                                !
!
!   Contribution of the atoms through the BL amplitudes
!
              SLF_1 = ZEROC                                         !
              SLF_2 = ZEROC                                         !
!
              IF(I_TEST == 2) GO TO 111                             !
!
              DO LF = LF1, LF2, ISTEP_LF                            !
               ILF = LF * LF + LF + 1                               !
               DO MF = -LF, LF                                      !
                INDF = ILF + MF                                     !
                IF((MF < MI-M_R) .OR. (MF > MI+M_R)) GO TO 555      !
                IF( (SELRULE == '00000') .AND.                    & !
                    (MF /=  MI)) GO TO 555                          !
!
                SLF_1 = SLF_1 + BL(INDF) * MLIL0(1,MI,INDF)         !
                IF(I_DICHR >= 1) THEN                               !
                  SLF_2 = SLF_2 + BL(INDF) * MLIL0(2,MI,INDF)       !
                END IF                                              !
!
 555            CONTINUE                                            !
               END DO                                               !
              END DO                                                !
!
!   Writing the amplitudes (whenever necessary)
!
 111          I_CALL = 2
              CALL WRITE_PED_AMPLITUDE(J_FILE,JPLAN,JEMET,JE,     & !
                                       J_FIXED,J_SCAN,JEPS,JDIR,  & !
                                       MI,SJDIR_1,SLF_1,SJDIR_2,  & !
                                       SLF_2,I_CALL,IOUT,IOUT2,   & !
                                       I_MI,STAT,OUTFILE)           !
!
!   Computing the square modulus
!
              IF(SPECTRO == 'PED') THEN
                SMIDIF_1 = SMIDIF_1 + ABS(SLF_1)   * ABS(SLF_1)     !
                SMIDIR_1 = SMIDIR_1 + ABS(SJDIR_1) * ABS(SJDIR_1)   !
                IF(I_DICHR >= 1) THEN                               !
                  SMIDIF_2 = SMIDIF_2 + ABS(SLF_2)   * ABS(SLF_2)   !
                  SMIDIR_2 = SMIDIR_2 + ABS(SJDIR_2) * ABS(SJDIR_2) !
                END IF                                              !
              END IF                                                !
!
!   End of the loop over MI
!
             END DO                                                 !
!
             IF(SPECTRO == 'APC') GO TO 220                         !

             SRDIR_1 = SRDIR_1 + SMIDIR_1                           !
             SRDIF_1 = SRDIF_1 + SMIDIF_1                           !
             IF(I_DICHR >= 1) THEN                                  !
               SRDIR_2 = SRDIR_2 + SMIDIR_2                         !
               SRDIF_2 = SRDIF_2 + SMIDIF_2                         !
             END IF                                                 !
!
 220         CONTINUE                                               !
!
!   End of the loop on the directions of the analyzer
!                         +
!   Multiplication by (a_0)^2 to go back to SI units
!
            END DO                                                  !
!
            IF(SPECTRO == 'APC') GO TO 221                          !
!
!            SEPSDIF_1 = SEPSDIF_1 +                               & !
!                        SRDIF_1 * VKR * CFM * BOHR * BOHR / NDIR    !
!            SEPSDIR_1 = SEPSDIR_1 +                               & !
!                        SRDIR_1 * VKR * CFM * BOHR * BOHR / NDIR    !
            SEPSDIF_1 = SEPSDIF_1 + SRDIF_1 * VKR * CFM / NDIR      !
            SEPSDIR_1 = SEPSDIR_1 + SRDIR_1 * VKR * CFM / NDIR      !
            IF(I_DICHR >= 1) THEN                                   !
              SEPSDIF_2 = SEPSDIF_2 +                             & !
                          SRDIF_2 * VKR * CFM * BOHR * BOHR / NDIR  !
              SEPSDIR_2 = SEPSDIR_2 +                             & !
                          SRDIR_2 * VKR * CFM * BOHR * BOHR / NDIR  !
            END IF                                                  !
!
 221        CONTINUE                                                !
!
!   End of the loop on the polarization
!
           END DO                                                   !
!
           SSETDIR_1 = SSETDIR_1 + SEPSDIR_1 * W                    !
           SSETDIF_1 = SSETDIF_1 + SEPSDIF_1 * W                    !
           IF(ICHKDIR == 2) THEN                                    !
             IF(JSET == JREF) THEN                                  !
               SSET2DIR_1 = SEPSDIR_1                               !
               SSET2DIF_1 = SEPSDIF_1                               !
             END IF                                                 !
           END IF                                                   !
           IF(I_DICHR >= 1) THEN                                    !
             SSETDIR_2 = SSETDIR_2 + SEPSDIR_2 * W                  !
             SSETDIF_2 = SSETDIF_2 + SEPSDIF_2 * W                  !
             IF(ICHKDIR == 2) THEN                                  !
               IF(JSET == JREF) THEN                                !
                 SSET2DIR_2 = SEPSDIR_2                             !
                 SSET2DIF_2 = SEPSDIF_2                             !
               END IF                                               !
             END IF                                                 !
           END IF                                                   !
!
!   End of the loop on the set averaging
!
           END DO                                                   !
!
!   Writing the cross-section in the result file
!
           IF(SPECTRO == 'APC') GO TO 222                           !
!
           IF(I_DICHR == 0) THEN                                    !
             IF(ISOM == 2) THEN                                     !
              WRITE(IOUT,67) JPLAN,J_FILE,DTHETAP,DPHIP,ECIN,     & !
                             SSETDIR_1,SSETDIF_1                    !
              IF(ICHKDIR == 2) THEN                                 !
               WRITE(IOUT,67) JPLAN,J_FILE,DTHETAP,DPHIP,ECIN,    & !
                              SSET2DIR_1,SSET2DIF_1                 !
              END IF                                                !
             ELSE                                                   !
              WRITE(IOUT,67) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,      & !
                             SSETDIR_1,SSETDIF_1                    !
              IF(ICHKDIR == 2) THEN                                 !
               WRITE(IOUT,67) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,     & !
                              SSET2DIR_1,SSET2DIF_1                 !
              END IF                                                !
             END IF                                                 !
           ELSE                                                     !
             IF(ISOM == 2) THEN                                     !
              WRITE(IOUT,72) JPLAN,J_FILE,DTHETAP,DPHIP,ECIN,     & !
                             SSETDIR_1,SSETDIF_1,                 & !
                             SSETDIR_2,SSETDIF_2                    !
              IF(ICHKDIR == 2) THEN                                 !
               WRITE(IOUT,72) JPLAN,J_FILE,DTHETAP,DPHIP,ECIN,    & !
                              SSET2DIR_1,SSET2DIF_1,              & !
                              SSET2DIR_2,SSET2DIF_2                 !
              END IF                                                !
             ELSE                                                   !
              WRITE(IOUT,72) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,      & !
                             SSETDIR_1,SSETDIF_1,                 & !
                             SSETDIR_2,SSETDIF_2                    !
              IF(ICHKDIR == 2) THEN                                 !
               WRITE(IOUT,72) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,     & !
                              SSET2DIR_1,SSET2DIF_1,              & !
                              SSET2DIR_2,SSET2DIF_2                 !
              END IF                                                !
             END IF                                                 !
           END IF                                                   !
!
 222       CONTINUE                                                 !
!
!   End of the loop on the scanned angle
!
          END DO                                                    !
!
   8      CONTINUE                                                  !
!
!   End of the loop on the fixed angle
!
         END DO                                                     !
!
         DEALLOCATE(TAU)                                            !
!
!     End of the loop on the energy
!
         CLOSE(IUI11)                                               !

        END DO                                                      !
   1    CONTINUE                                                    !
!
!   End of the loop on the planes
!
      END DO                                                        !
!
      IF(ABS(I_EXT) >= 1) CLOSE(IUI11)                              !
      IF((ISOM == 0) .OR. (J_FILE == N_INFILES)) WRITE(IOUT,*)      !
      IF(SPECTRO == 'APC') CLOSE(IOUT)                              !
      IF(SPECTRO == 'APC') GO TO 7                                  !
      IF((NPLAN > 0) .AND. (ISOM == 0)) THEN                        !
        NP = 0                                                      !
        CALL TREAT_PED(ISOM,N_INFILES,J_FILE,NP)                    !
      END IF                                                        !
      IF(I_EXT == 2) THEN                                           !
         CALL WEIGHT_SUM(I_EXT,0,1)                                 !
      END IF                                                        !
!
!  Writing spectroscopy-related information in the check file unit IUO1
!
      CALL WRITE_INFO_END(IUO1)                                     !
!
      GO TO 7                                                       !
!
   6  WRITE(IUO1,55)                                                !
!
!  Formats:
!
   9  FORMAT(9(2X,I1),2X,I2)
  11  FORMAT(I4,2X,I4,2X,I4)
  12  FORMAT(2X,A3,11X,A13)
  13  FORMAT(6X,I1,1X,I3,2X,I4)
  14  FORMAT(6X,I1,1X,I3,3X,I3)
  21  FORMAT(10X,E12.6,3X,E12.6)
  22  FORMAT(16X,'INTERNAL CALCULATION OF MEAN SQUARE DISPLACEMENTS',/, &
             25X,' BY DEBYE UNCORRELATED MODEL:',/)
  23  FORMAT(21X,'ATOM TYPE ',I5,'  MSD = ',F8.6,' ANG**2')
  52  FORMAT(/////,2X,'******* PLANE NUMBER ',I3,' POSITION OF ',       &
      'THE ABSORBER : (',F6.3,',',F6.3,',',F6.3,') *******',/,2X,       &
      '******* ',19X,'THIS ABSORBER IS OF TYPE ',I2,20X,' *******')
  55  FORMAT(///,12X,' <<<<<<<<<<  THIS VALUE OF IMFP IS NOT',          &
      'AVAILABLE  >>>>>>>>>>')
  56  FORMAT(4X,'LATTICE PARAMETER A = ',F6.3,' ANGSTROEMS',4X,         &
      'MEAN FREE PATH = ',F6.3,' * A',//)
  59  FORMAT(//,15X,'THE SCATTERING DIRECTION IS GIVEN INSIDE ',        &
      'THE CRYSTAL')
  60  FORMAT(7X,'THE POSITIONS OF THE ATOMS ARE GIVEN WITH RESPECT ',   &
      'TO THE ABSORBER')
  61  FORMAT(///,4X,'..........  DIRECTION OF THE DETECTOR      : (',   &
            F6.3,',',F6.3,',',F6.3,                                     &
             ')  ..........',/,16X,'DIRECTION OF THE LIGHT    ',        &
             '     : (',F6.3,',',F6.3,',',F6.3,                         &
            ')',/,16X,'DIRECTION OF THE POLARIZATION  : (',             &
            F6.3,',',F6.3,',',F6.3,')',/,16X,'ANALYZER.LIGHT       ',   &
             '          :        ',F7.4)
  63  FORMAT(///,4X,'..........  DIRECTION OF THE DETECTOR      : (',   &
             F6.3,',',F6.3,',',F6.3,                                    &
             ')  ..........',/,16X,'DIRECTION OF THE LIGHT    ',        &
             '     : (',F6.3,',',F6.3,',',F6.3,')',/,16X,               &
             'ANALYZER.LIGHT               :        ',F7.4)
  65  FORMAT(////,3X,'++++++++++++++++++',9X,                           &
      'THETA = ',F6.2,' DEGREES',9X,'++++++++',                         &
      '++++++++++',///)
  66  FORMAT(////,3X,'++++++++++++++++++',9X,                           &
      'PHI = ',F6.2,' DEGREES',9X,'++++++++++',                         &
      '++++++++++',///)
  67  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,              &
             2X,E12.6)
  72  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,              &
             2X,E12.6,2X,E12.6,2X,E12.6)
  79  FORMAT(2X,I3,2X,I2,2X,I4,2X,I4,2X,I4)
  86  FORMAT(2X,I3,1X,I4,5X,F8.3,3X,F8.3,3X,E12.6)
  88  FORMAT(/,19X,'TILTED THETA =',F6.2,5X,'TILTED PHI =',             &
               F6.2)
  89  FORMAT(/,4X,'..........................................',         &
               '.....................................')
!
   7  RETURN                                                        !
!
      END SUBROUTINE PED_MS
!
!=======================================================================
!
      SUBROUTINE WRITE_PED_AMPLITUDE(J_FILE,JPLAN,JEMET,JE,J_FIXED,    &
                                     J_SCAN,JEPS,JDIR,MI,SJDIR_1,SLF_1,&
                                     SJDIR_2,SLF_2,I_CALL,IOUT,IOUT2,  &
                                     I_MI,STAT,OUTFILE)
!
!  This routine writes the amplitude of the signal(s) detected in a file,
!    either for APECS, or for further use in PED (orientated orbitals' case)
!
!
!  Input variables :
!
!                       J_FILE    :  input data file number
!                       JPLAN     :  plane number
!                       JEMET     :  absorber number (always 1 when using phagen_scf.f)
!                       JE        :  energy point
!                       J_FIXED   :  number of the "fixed" angle point
!                       J_SCAN    :  number of the "scanned" angle point
!                       JEPS      :  number of the polarization direction
!                       JDIR      :  number of the analyzer direction (for averaging)
!                       MI        :  azimuthal quantum number of the excited core state
!                       SJDIR_1   :  atomic amplitude
!                       SLF_1     :  multiple scattering amplitude
!                       SJDIR_2   :  second atomic amplitude (dichroic case)
!                       SLF_2     :  second multiple scattering amplitude (dichroic case)
!                       I_CALL    :  number of the call:
!                                       ---> 1 : writes the beginning of the file
!                                            2 : writes the second part of the beginning
!                                            3 : writes the amplitudes
!
!
!  Output variables :
!
!                       IOUT      :  cross-section file number
!                       IOUT2     :  amplitude file number
!                       I_MI      :  switch for amplitudes
!                       STAT      :  status of the output files
!                       OUTFILE   :  name of the cross-section file
!
!
!   Author: Didier Sébilleau
!
!                                        Last modified : 28 May 2021
!
!
      USE CURRENT_AVER
      USE CURRENT_CALC
      USE CURRENT_INIT_VAL
      USE CURRENT_FINA_VAL
      USE CURRENT_FIXSCAN
      USE CURRENT_SPIN
      USE DICHROISM
      USE EXP_TYPE
      USE INIT_J
      USE NON_BULK
      USE OUTFILES
      USE OUTUNITS
      USE TESTS
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 24)     ::  OUTFILE
      CHARACTER (LEN = 24)     ::  AMPFILE
      CHARACTER (LEN = 13)     ::  OUTDATA2
      CHARACTER (LEN =  7)     ::  STAT
!
      INTEGER, INTENT(IN)      ::  J_FILE,JPLAN,JEMET,JE
      INTEGER, INTENT(IN)      ::  J_FIXED,J_SCAN,JEPS,JDIR
      INTEGER, INTENT(IN )     ::  MI,I_CALL
      INTEGER, INTENT(OUT)     ::  IOUT,IOUT2,I_MI
!
      INTEGER                  ::  N_DOT,J_CHAR
!
      COMPLEX (WP), INTENT(IN) ::  SLF_1,SLF_2,SJDIR_1,SJDIR_2
!
!  First call to the subroutine
!
      IF(I_CALL == 1) THEN                                          !
!
!  Choosing the output filename
!
        IF(I_AMP == 1) THEN                                         !
          I_MI     = 1                                              !
          OUTDATA2 = 'MS AMPLITUDES'                                !
        ELSE                                                        !
          I_MI     = 0                                              !
        END IF                                                      !
!
        IF((SPECTRO == 'PED') .OR. (SPECTRO == 'AED')) THEN         !
          IOUT    = IUO2                                            !
          OUTFILE = OUTFILE2                                        !
          STAT    = 'UNKNOWN'                                       !
          IF(I_MI == 1) THEN                                        !
            IOUT2 = IUSCR2 + 1                                      !
            N_DOT = 1                                               !
            DO J_CHAR = 1, 24                                       !
              IF(OUTFILE(J_CHAR:J_CHAR) == '.') GO TO 888           !
              N_DOT = N_DOT + 1                                     !
            END DO                                                  !
  888       CONTINUE                                                !
            AMPFILE = OUTFILE(1:N_DOT)//'amp'                       !
            OPEN(UNIT=IOUT2, FILE=AMPFILE, STATUS=STAT)             !
          END IF                                                    !
        ELSE IF(SPECTRO == 'APC') THEN                              !
          IOUT    = IUSCR2 + 1                                      !
          OUTFILE = 'res/phot.amp'                                  !
          STAT    = 'UNKNOWN'                                       !
        END IF                                                      !
!
!  Second call to the subroutine
!
      ELSE IF(I_CALL == 2) THEN                                     !
!
        IF((ISOM == 0) .OR. ((ISOM > 0) .AND. (J_FILE == 1))) THEN  !
          IF(I_MI == 1) THEN                                        !
            WRITE(IOUT2,12) SPECTRO,OUTDATA2                        !
            WRITE(IOUT2,12) STEREO                                  !
            WRITE(IOUT2,19) ISPIN,I_DICHR,I_SO,ISFLIP,ICHKDIR,IPHI,&!
                            ITHETA,IE,IPH_1,I_EXT                   !
            WRITE(IOUT2,20) PHI0,THETA0,PHI1,THETA1,NONVOL(1)       !
          END IF                                                    !
        END IF                                                      !
!
!  Third call to the subroutine
!
      ELSE IF(I_CALL == 3) THEN                                     !
!
!  Writing the amplitudes
!
        IF(SPECTRO == 'APC') THEN                                   !
          WRITE(IOUT,87) J_FILE,JPLAN,JEMET,JE,J_FIXED,J_SCAN,    & !
                         JEPS,JDIR,MI,SJDIR_1,SLF_1                 !
          IF(I_DICHR >= 1) THEN                                     !
            WRITE(IOUT,87) J_FILE,JPLAN,JEMET,JE,J_FIXED,J_SCAN,  & !
                           JEPS,JDIR,MI,SJDIR_2,SLF_2               !
          END IF                                                    !
        ELSE                                                        !
          IF(I_MI == 1) THEN                                        !
            WRITE(IOUT2,87) J_FILE,JPLAN,JEMET,JE,J_FIXED,        & !
                            J_SCAN,JEPS,JDIR,MI,SJDIR_1,          & !
                            SLF_1                                   !
            IF(I_DICHR >= 1) THEN                                   !
              WRITE(IOUT2,87) J_FILE,JPLAN,JEMET,JE,J_FIXED,      & !
                              J_SCAN,JEPS,JDIR,MI,SJDIR_2,        & !
                              SLF_2                                 !
            END IF                                                  !
          END IF                                                    !
        END IF                                                      !
!
      END IF                                                        !
!
!  Formats:
!
  12  FORMAT(2X,A3,11X,A13)
  19  FORMAT(2(2X,I1),1X,I2,6(2X,I1),2X,I2)
  20  FORMAT(2(5X,F6.2,2X,F6.2),2X,I1)
  87  FORMAT(2X,I2,2X,I3,2X,I2,2X,I3,2X,I3,2X,I3,2X,I1,2X,I2,2X,I2, &
             2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
!
      END SUBROUTINE WRITE_PED_AMPLITUDE
!
END MODULE CALC_PED_MS
