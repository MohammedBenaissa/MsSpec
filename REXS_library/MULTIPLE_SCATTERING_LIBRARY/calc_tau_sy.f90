!
!======================================================================
!
MODULE CALC_TAU0J_SY
!
!  This module provides the subroutine to compute
!    Tau^{0 j} or Tau^{j 0} in the SYMMETRIZED Rehr-Albers
!    series expansion case
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE CALC_TAU_SY(TAU_TYPE,IATL,JE,TAU)
!
!  This subroutine computes the scattering path operator Tau^{0 0},
!    Tau^{0 j} or Tau^{j 0}
!
!  Input variables:
!
!     TAU_TYPE : scattering path operator type:
!                   'ABS' --> absorption case       Tau^{0 0}
!                   'INC' --> incoming beam case    Tau^{j 0}
!                   'OUT' --> outgoing beam case    Tau^{0 J}
!     IATL     : index of the atom for which BL is computed
!     JE       : energy point for which BL is computed
!
!  Output variable:
!
!     TAU     : scattering path operator
!                   Note: Tau^{j 0}, Tau^{0 j} and Tau^{0 0} are stored
!                         in the same way to avoid having 3 TAU arrays
!
!  This is the SYMMETRIZED series expansion version
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 22 Jun 2021
!
!
      USE REAL_NUMBERS,              ONLY : ZERO,ONE,TENTH,LARGE
      USE COMPLEX_NUMBERS,           ONLY : ZEROC,ONEC
!
      USE CALC_TYPE
      USE CLUSTER_COORD
      USE CLUSTER_LIMITS,            ONLY : VAL
      USE CURRENT_BEAM
      USE CURRENT_T_MATRIX,          ONLY : LMAX
      USE EXTREMES
      USE INIT_L
      USE OUTUNITS
      USE PATH_INFO
      USE PRINT_PATH_INFO
      USE TESTPA
      USE TESTPB
      USE TESTS
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)      ::  TAU_TYPE
!
      INTEGER, INTENT(INOUT)    ::  IATL
      INTEGER, INTENT(IN)       ::  JE
!
      INTEGER                   ::  NTYPEM,JNEM
      INTEGER                   ::  JPOS(0:N_SCAT_M,3),JPA(N_SCAT_M)
      INTEGER                   ::  NSCAT,N_SCA,I_CP,I_CP_PW,ND
      INTEGER                   ::  JATLEM,JAEM,JD,JP
      INTEGER                   ::  ITYP,INUM,NBTYPI
      INTEGER                   ::  JPT,JOP,JLINE,KD
      INTEGER                   ::  JOPA,NOPA,NPATOT
      INTEGER                   ::  DELTA1,DELTA2
      INTEGER                   ::  LMJMAX
!
      INTEGER                   ::  INT
!
      REAL (WP)                 ::  NPATH1(0:N_SCAT_M)
      REAL (WP)                 ::  XR(N_SCAT_M),YR(N_SCAT_M)
      REAL (WP)                 ::  ZR(N_SCAT_M),R(N_SCAT_M)
      REAL (WP)                 ::  ZSURF,THMI,PHMI,DIJ
      REAL (WP)                 ::  XPATOT,XINT0,DIST0
      REAL (WP)                 ::  FREF,XMAX
!
      REAL (WP)                 ::  FLOAT
!
      COMPLEX (WP), INTENT(OUT) ::  TAU(LINMAX,LINFMAX,NATCLU_M)
!
      COMPLEX (WP)              ::  PW(0:N_SCAT_M)
      COMPLEX (WP)              ::  RHOMI
      COMPLEX (WP)              ::  ZM,ZMF,ZLJ,ZMJ
      COMPLEX (WP)              ::  ZM1,ZM2,Z_M1,Z_M2
      COMPLEX (WP)              ::  Z_M_P,Z_MF1,Z_MF2
!
      NTYPEM = 1                                                    ! number of absorbers
      JNEM   = 1                                                    ! type of current absorber
!
      NSCAT = NATCLU - 1                                            !
!
      ZSURF = VAL(1)                                                !
!
      FMIN(0) = ONE                                                 !
      FMAX(0) = ONE                                                 !
!
!   Beam-dependent quantities
!
      IF(TAU_TYPE == 'ABS') THEN                                    !
        I_CP    = 1                                                 !
        I_CP_PW = 1                                                 !
      ELSE IF(TAU_TYPE == 'INC') THEN                               !
        I_CP    = 0                                                 !
        I_CP_PW = 1                                                 !
      ELSE IF(TAU_TYPE == 'OUT') THEN                               !
        I_CP    = 0                                                 !
        I_CP_PW = 0                                                 !
      END IF                                                        !
!
!   Initialization of TAU
!
      CALL INIT_TAU_M(JE,N_PROT,NATYP,LMAX,TAU,LMJMAX)              !
!
!   Initializations for the calculation of TAU
!
      IF(I_TEST == 2) GO TO 10                                      !
      IF(IPW == 0) THEN                                             !
!
        PW(0)  = ONEC                                               !
        PW(1)  = ONEC                                               !
        ND     = 0                                                  !
        TH01   = ZERO                                               !
        PHI01  = ZERO                                               !
        RHO01  = ZEROC                                              !
        THMI   = ZERO                                               !
        PHMI   = ZERO                                               !
        RHOMI  = ZEROC                                              !
        JATLEM = JNEM                                               !
        IF(NTYPEM > 1) THEN                                         !
          DO JAEM = NTYPEM-1, 1, -1                                 !
            JATLEM = JATLEM + NATYP(JAEM)                           !
          END DO                                                    !
        END IF                                                      !
        DO JD = 1, N_SCAT                                           !
          NPATH2(JD) = ZERO                                         !
          NPATH(JD)  = ZERO                                         !
          IT(JD)     = 0                                            !
          IN(JD)     = 0                                            !
          FMIN(JD)   = LARGE                                        !
          FMAX(JD)   = ZERO                                         !
        END DO                                                      !
        NTHOF = 0                                                   !
!
! Calculation of the maximal intensity for the paths of order NCUT
!   (plane waves). This will be taken as a reference for the IPW filter.

      ELSE IF(IPW == 1) THEN                                        !
        CALL PW_FILTER_INI(JE,N_PROT,NTYPEM,JNEM,NATYP,LMAX,      & !
                           I_CP_PW,VAL,FREF)                        !
      END IF                                                        !
!
!  Initialization of the paths
!
      IF(NPATHP > 0) THEN                                           !
        DO JP = 1, NPATHP-1                                         !
          FMN(JP)  = ZERO                                           !
          PATH(JP) = ZERO                                           !
          JON(JP)  = ZERO                                           !
        END DO                                                      !
        FMN(NPATHP)  = - ONE                                        !
        PATH(NPATHP) = ZERO                                         !
        JON(NPATHP)  = ZERO                                         !
      END IF                                                        !
      IREF = 0                                                      !
      IJ   = 1                                                      !
      IF(IPRINT == 3) THEN                                          !
        OPEN(UNIT=IUSCR, STATUS='SCRATCH')                          !
      END IF                                                        !
!
!  Calling the path generation
!
      IF(TAU_TYPE /= 'INC') THEN                                    !
        DIJ = ZERO                                                  !
        CALL FINDPATHS(ND,NTYPEM,JATLEM,I_CP,R,XR,YR,ZR,RHOMI,    & !
                       THMI,PHMI,ZSURF,JPOS,PW,JE,FREF,DIJ,TAU)     !
      ELSE                                                          !
        IREF = 0                                                    !
        IJ   = 1                                                    !
        DIJ  = ZERO                                                 !
        IATL = 0                                                    !
        DO ITYP = 1, N_PROT                                         !
          NBTYPI = NATYP(ITYP)                                      !
          DO INUM = 1, NBTYPI                                       !
            IATL      = IATL + 1                                    !  -->  IATL is redefined here ! CHECK
            JPOS(0,1) = ITYP                                        !
            JPOS(0,2) = INUM                                        !
            JPOS(0,3) = IATL                                        !
            LF1       = 0                                           !
            LF2       = LMAX(ITYP,JE)                               !
            ISTEP_LF  = 1                                           !
            CALL FINDPATHS(ND,ITYP,IATL,I_CP,R,XR,YR,ZR,RHOMI,    & !
                           THMI,PHMI,ZSURF,JPOS,PW,JE,FREF,DIJ,TAU) !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
!   Printing of the paths
!
      IF(NPATHP == 0) GO TO 10                                      !
      IF(NSCAT > 1) THEN                                            !
        XPATOT = ( FLOAT(NSCAT)**FLOAT(N_SCAT+1) - ONE ) /        & !
                 FLOAT(NSCAT-1)                                     !
      ELSE                                                          !
        XPATOT = FLOAT(N_SCAT + 1)                                  !
      END IF                                                        !
      IF(XPATOT < 2.14748E+09_WP) THEN                              !
        NPATOT = INT(XPATOT)                                        !
        IF(NPATOT < NPATHP) NPATHP = NPATOT - 1                     !
      END IF                                   !                    !
      WRITE(IUO1,200) NPATHP                                        !
      WRITE(IUO1,201)                                               !
      DO JPT = 1, NPATHP                                            !
        IF(PATH(NPATHP) > 2.14E+09_WP) THEN                         !
          WRITE(IUO1,202) JPT,JON(JPT),PATH(JPT),FMN(JPT),        & !
                          DMN(JPT),JNEM,(JPON(JPT,KD),KD=1,JON(JPT))!
        ELSE                                                        !
          WRITE(IUO1,203) JPT,JON(JPT),INT(PATH(JPT)),FMN(JPT),   & !
                          DMN(JPT),JNEM,(JPON(JPT,KD),KD=1,JON(JPT))!
        END IF                                                      !
      END DO                                                        !
      IF(IPRINT == 3) THEN                                          !
        IF(XPATOT > 2.14748E+09_WP) GO TO 20                        !
        WRITE(IUO1,204)                                             !
        WRITE(IUO1,205)                                             !
        NPATOT = INT(XPATOT)                                        !
        DO JOP = 0, N_SCAT                                          !
          IF(JOP == 0) THEN                                         !
            XINT0 = FMAX(0)                                         !
            DIST0 = ZERO                                            !
            WRITE(IUO1,206) JOP,JOP+1,XINT0,DIST0,JNEM              !
            GO TO 25                                                !
          END IF                                                    !
          WRITE(IUO1,207)                                           !
          DO JLINE = 1, NPATOT-1                                    !
            READ(IUSCR,100,ERR=25,END=25) JOPA,NOPA,XMAX,DIST0,   & !
                                          (JPA(KD),KD=1,JOPA)       !
            IF(JOPA == JOP) THEN                                    !
              IF(NOPA > 2.14E+09_WP) THEN                           !
               WRITE(IUO1,208) JOPA,NOPA,XMAX,DIST0,JNEM,         & !
                               (JPA(KD),KD=1,JOPA)                  !
              ELSE                                                  !
               WRITE(IUO1,206) JOPA,INT(NOPA),XMAX,DIST0,JNEM,    & !
                               (JPA(KD),KD=1,JOPA)                  !
              END IF                                                !
            END IF                                                  !
          END DO                                                    !
          IF(JOP == N_SCAT) WRITE(IUO1,209)                         !
  25      REWIND IUSCR                                              !
        END DO                                                      !
        GO TO 30
!
  20    WRITE(IUO1,210)                                             !
        CLOSE(IUSCR,STATUS='DELETE')                                !
!
  30  END IF                                                        !
!
      DO JD = 0, N_SCAT                                             !
        NPATH1(JD) = FLOAT(NSCAT)**FLOAT(JD)                        !
        IF(NPATH1(JD) > 2.14E+09_WP) THEN                           !
          IF(FMIN(JD) == 0.1E+21_WP) FMIN(JD) = ZERO                !
          WRITE(IUO1,211) JD,NPATH1(JD),NPATH2(JD),FMIN(JD),      & !
                          NPMI(JD),FMAX(JD),NPMA(JD)                !
          IF( (IPW == 1) .AND. (JD > NCUT) )  THEN                  !
              WRITE(IUO1,213) FREF * PCTINT                         !
          END IF                                                    !
        ELSE                                                        !
          IF(FMIN(JD) == 0.1E+21_WP) FMIN(JD) = ZERO                !
          WRITE(IUO1,212) JD,INT(NPATH1(JD) + TENTH),             & !
                             INT(NPATH2(JD) + TENTH),FMIN(JD),    & !
                             INT(NPMI(JD)   + TENTH),FMAX(JD),    & !
                             INT(NPMA(JD)   + TENTH)                !
          IF( (IPW == 1) .AND. (JD > NCUT) ) THEN                   !
             WRITE(IUO1,213) FREF * PCTINT                          !
          END IF                                                    !
        END IF                                                      !
      END DO                                                        !
 10   CONTINUE                                                      !
!
!   Calculation of the missing values of TAU(INDJ,INDF,JTYP) when
!     there is a relation between negative and positive matrix elements
!
         DO JTYP = 1, N_PROT                                        !
           LMJ  = LMAX(JTYP,JE)                                     !
           ILM  = I_LM(JTYP)                                        !
           I_LJ = ISTEP_L(JTYP)                                     !
           I_MJ = ISTEP_M(JTYP)                                     !
           IRL  = I_REL_MP(JTYP)                                    !
           IF((ILM == 1) .AND. (I_LJ == 2)) ILM = 0                 !
           ILMJ = MOD(LMJ,2)                                        !
!
           IF(IRL == 1) THEN                                        !
             RL = Z_L_P(JTYP)                                       !
             ZM = Z_M_P(JTYP)                                       !
           ELSE                                                     !
             GO TO 40                                               !
           END IF                                                   !
!
           DO LF = LF1, LF2, ISTEP_LF                               !
             ILF = LF * LF + LF + 1                                 !
!
             IF(I_LJ == 1) THEN                                     !
               LJ0 = 0                                              !
               LJ1 = LMJ                                            !
             ELSE                                                   !
               LJ0 = MOD(LF,2)                                      !
               LJ1 = LMJ - MOD(LF + ILMJ,2)                         !
             END IF                                                 !
!
             ZLF = RL**LF                                           !
!
             DO MF = -LF, LF                                        !
               INDF  = ILF + MF                                     !
               INDF_ = ILF - MF                                     !
!
               ZMF = ZLF * (ZM**MF)                                 !
!
               DO LJ = LJ0, LJ1, I_LJ                               !
                 ILJ = LJ * LJ + LJ + 1                             !
!
                 ZLJ = ZMF * (RL**LJ)                               !
!
                 IF(I_MJ == 1) THEN                                 !
                   MJ0 = - LJ                                       !
                   MJ1 =   LJ                                       !
                 ELSE                                               !
                   N2 = 0                                           !
                   IF(ILM == 1) THEN                                !
                     LJMLF = MOD(ABS(LJ - LF),2)                    !
                     IF(LJMLF == 1) THEN                            !
                       N2 = I_MJ / 2                                !
                     END IF                                         !
                   END IF                                           !
                   DELTA1 = MOD(ABS(LJ + MF + N2),I_MJ)             !
                   DELTA2 = MOD(ABS(LJ - MF + N2),I_MJ)             !
                   MJ0    = - LJ + DELTA1                           !
                   MJ1    =   LJ - DELTA2                           !
                 END IF                                             !
                 MJ1 = - MOD(ABS(MJ0),I_MJ)                         !
!
                 DO MJ = MJ0, MJ1, I_MJ                             !
                   INDJ  = ILJ + MJ                                 !
                   INDJ_ = ILJ - MJ                                 !
!
                   ZMJ = ZLJ / (ZM**MJ)                             !
!
                   IF(MJ /= 0) THEN                                 !
                     TAU(INDJ_,INDF_,JTYP) = TAU(INDJ,INDF,JTYP) *& !
                                             ZMJ                    !
                   END IF                                           !
                 END DO                                             !
               END DO                                               !
             END DO                                                 !
           END DO                                                   !
  40       CONTINUE                                                 !
         END DO                                                     !
!
!
      END SUBROUTINE CALC_TAU_SY
!
!======================================================================
!
      SUBROUTINE INIT_TAU_M(JE,N_PROT,NATYP,LMAX,TAU,LMJMAX)
!
!  This subroutine initializes to zero the TAU array for outgoing beam
!
!  Input variables:
!
!     JE       : energy point for which TAU is computed
!     N_PROT   : number of prototypical atoms
!     LMAX     : cut-off for the L-expansions
!     NATYP    : prototypical atoms indices array
!
!  Output variable:
!
!     TAU      : scattering path operator array
!     LMJMAX   : largest value of l_max
!
!
!
!  This is the SYMMETRIZED version
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 22 Jun 2021
!
!
      USE COMPLEX_NUMBERS,           ONLY : ZEROC
!
      USE INIT_L
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  JE,N_PROT
      INTEGER, INTENT(IN)       ::  LMAX(NATM,NE_M),NATYP(NATM)
      INTEGER, INTENT(OUT)      ::  LMJMAX
!
      INTEGER                   ::  JATL,JTYP,JNUM
      INTEGER                   ::  NBTYP,LMJ
      INTEGER                   ::  LJ,ILJ,MJ,INDJ
      INTEGER                   ::  LF,ILF,MF,INDF
!
      COMPLEX (WP), INTENT(OUT) ::  TAU(LINMAX,LINFMAX,NATP_M)
!
      LMJMAX = 0                                                    !
      JATL   = 0                                                    !
!
      DO JTYP = 1, N_PROT                                           !
!
        NBTYP  = NATYP(JTYP)                                        !
        LMJ    = LMAX(JTYP,JE)                                      !
        LMJMAX = MAX(LMJMAX,LMJ)                                    !
!
        DO LJ = 0, LMJ                                            !
!
          ILJ = LJ * LJ + LJ + 1                                  !
!
          DO MJ = -LJ, LJ                                         !
!
            INDJ = ILJ + MJ                                       !
!
            DO LF = LF1, LF2, ISTEP_LF                            !
!
              ILF = LF * LF + LF + 1                              !
!
              DO MF = -LF, LF                                     !
!
                INDF = ILF + MF                                   !
                TAU(INDJ,INDF,JTYP) = ZEROC                       !
!
              END DO                                              !
!
            END DO                                                !
!
          END DO                                                  !
!
        END DO                                                      !
!
      END DO                                                        !
!
      END SUBROUTINE INIT_TAU_M
!
!======================================================================
!
      SUBROUTINE PW_FILTER_INI(JE,N_PROT,NTYPEM,JNEM,NATYP,LMAX,  & !
                               I_CP,VAL,FREF)                       !
!
!  This subroutine computes the reference value for the IPW filter
!    (maximal intensity for the paths of order NCUT in plane waves)
!
!  Input variables:
!
!     JE       : energy point for which BL is computed
!     N_PROT   : number of prototypical atoms
!     NTYPEM   : index of the type of absorber
!     JNEM     : index of the absorber
!     NATYP    : prototypical atoms indices array
!     LMAX     : cut-off for the L-expansions
!     I_CP     : switch for close paths (=0) or open paths (=1)
!     VAL      : plane indices array
!
!  Output variable:
!
!     FREF     : reference value for the IPW filter
!
!  Note:
!
!     The paths are considered invariant by reversal of the direction.
!       Therefore, the FREF for Tau^{j 0} is taken as the FREF
!       for Tau^{0 j} (I_CP = 1)
!
!
!  This is the SYMMETRIZED version
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 22 Jun 2021
!
!
      USE REAL_NUMBERS,              ONLY : ZERO,ONE,LARGE
      USE COMPLEX_NUMBERS,           ONLY : ZEROC,ONEC
!
      USE CALC_TYPE
      USE CURRENT_BEAM
      USE EXTREMES
      USE PATH_INFO
      USE TESTPA
      USE TESTPB
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  JE,N_PROT,NTYPEM,JNEM
      INTEGER, INTENT(IN)       ::  LMAX(NATM,NE_M),NATYP(NATM)
      INTEGER, INTENT(INOUT)    ::  I_CP
!
      INTEGER                   ::  JPOS(0:N_SCAT_M,3)
      INTEGER                   ::  ND,JATLEM,JAEM,JD
      INTEGER                   ::  N_SCATOLD,NOOLD,I_BASISOLD
      INTEGER                   ::  LMJMAX
!
      REAL (WP), INTENT(IN)     ::  VAL(NATCLU_M)
      REAL (WP), INTENT(OUT)    ::  FREF
!
      REAL (WP)                 ::  XR(N_SCAT_M),YR(N_SCAT_M)
      REAL (WP)                 ::  ZR(N_SCAT_M),R(N_SCAT_M)
      REAL (WP)                 ::  ZSURF,THMI,PHMI
      REAL (WP)                 ::  DIJ
!
      COMPLEX (WP)              ::  TAU(LINMAX,LINFMAX,NATCLU_M)
      COMPLEX (WP)              ::  PW(0:N_SCAT_M)
      COMPLEX (WP)              ::  RHOMI
!
      ZSURF = VAL(1)                                                !
!
      PW(0)  = ONEC                                                 !
      PW(1)  = ONEC                                                 !
      ND     = 0                                                    !
      TH01   = ONE                                                  !
      PHI01  = ONE                                                  !
      RHO01  = ZEROC                                                !
      THMI   = ONE                                                  !
      PHMI   = ONE                                                  !
      RHOMI  = ZEROC                                                !
      JATLEM = JNEM                                                 !
!
      IF(NTYPEM > 1) THEN                                           !
        DO JAEM = NTYPEM-1, 1, -1                                   !
          JATLEM = JATLEM + NATYP(JAEM)                             !
        END DO                                                      !
      END IF                                                        !
!
      DO JD = 1, N_SCAT                                             !
        NPATH2(JD) = ZERO                                           !
        NPATH(JD)  = ZERO                                           !
        IT(JD)     = 0                                              !
        IN(JD)     = 0                                              !
        FMIN(JD)   = LARGE                                          !
        FMAX(JD)   = ZERO                                           !
      END DO                                                        !
!
      NTHOF = 0                                                     !
!
      N_SCATOLD    = N_SCAT                                         !
      NOOLD      = NO                                               !
      I_BASISOLD = I_BASIS                                          !
      N_SCAT       = NCUT                                           !
      NO         = 0                                                !
      I_BASIS    = 0                                                !
      IREF       = 1                                                !
      IPW        = 0                                                !
      IJ         = 0                                                !
      DIJ        = ZERO                                             !
      FREF       = ZERO                                             !
!
      CALL FINDPATHS(ND,NTYPEM,JATLEM,I_CP,R,XR,YR,ZR,RHOMI,      & !
                     THMI,PHMI,ZSURF,JPOS,PW,JE,FREF,DIJ,TAU)       !
!
      N_SCAT    = N_SCATOLD                                             !
      NO      = NOOLD                                               !
      I_BASIS = I_BASISOLD                                          !
      PW(0)   = ONEC                                                !
      PW(1)   = ONEC                                                !
      IPW     = 1                                                   !
      ND      = 0                                                   !
      TH01    = ZERO                                                !
      PHI01   = ZERO                                                !
      RHO01   = ZEROC                                               !
      THMI    = ZERO                                                !
      PHMI    = ZERO                                                !
      RHOMI   = ZEROC                                               !
      JATLEM  = JNEM                                                !
!
      IF(NTYPEM > 1) THEN                                           !
        DO JAEM = NTYPEM-1, 1, -1                                   !
          JATLEM = JATLEM + NATYP(JAEM)                             !
        END DO                                                      !
      END IF                                                        !
!
      DO JD = 1, N_SCAT                                             !
        NPATH2(JD) = ZERO                                           !
        NPATH(JD)  = ZERO                                           !
        IT(JD)     = 0                                              !
        IN(JD)     = 0                                              !
        FMIN(JD)   = LARGE                                          !
        FMAX(JD)   = ZERO                                           !
      END DO                                                        !
!
      NTHOF = 0                                                     !
!
!   New initialization of TAU after the PW calculation
!
      CALL INIT_TAU_M(JE,N_PROT,NATYP,LMAX,TAU,LMJMAX)              !
!
      END SUBROUTINE PW_FILTER_INI
!
!=======================================================================
!
      RECURSIVE SUBROUTINE FINDPATHS(ND,ITYP,IATL,I_CP,R,XR,YR,ZR,  &
                                     RHOMI,THMI,PHIMI,ZSURF,JPOS,   &
                                     PW,JE,FREF,DIJ,TAU)
!
!  This routine generates all the paths and filters them according to the
!      criteria given in the input data file (IFSPH,IFWD,IPW,ILENGTH).
!      It corresponds to the spin-independent case from a non spin-orbit
!                    resolved initial core state LI
!
!  Starting from two input atoms M and I in the path, it will
!      generate atoms J and K
!
!        M        J -- K
!          \    /
!           \  /
!             I
!       \_______/ \_____/
!           |        |
!           v        |
!         input      |->   constructed here
!
!  Input variables:
!
!     ND       : position of previous atom in the path
!     ITYP     : type of the atom from which the subroutine is called
!     IATL     : index of the atom from which the subroutine is called
!     I_CP     : closed path index
!                      I_CP = 0 : all open paths generated
!                      I_CP = 1 : only closed paths generated
!     RHOMI    : \
!     THMI     :  > k * r, theta and phi for incoming direction M-I
!     PHIMI    : /
!     ZSURF    : z-position of absorbing atom
!     JPOS     : position of atom ND in the path
!     PW       :
!     JE       : current energy point
!     FREF     : reference value for path selection
!     DIJ      : initial path length
!
!  Output variable:
!
!     ND       : position of atom K in the path
!     DIJ      : current path length
!     TAU      : scattering path operator of the paths generated
!
!
!      This is the SYMMETRIZED version.
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 22 Jun 2021
!
!
      USE REAL_NUMBERS,              ONLY : ONE
      USE COMPLEX_NUMBERS,           ONLY : IC
!
      USE CLUSTER_COORD
      USE CURRENT_BEAM
      USE CURRENT_T_MATRIX
      USE DEB_WAL_CLU
      USE INIT_L
      USE PATH_INFO
      USE ROTATION
      USE TEMP_STORAGE_1
      USE TESTPA
      USE TESTPB
!
      USE ROT_CUB
      USE SYM_OP
      USE TAU_PROT
      USE TAU_SYM_OP
!
      USE EULER_ANGLES
      USE VIBRATIONS
      USE WIGNER_ROTATIONS
!
      INTEGER, INTENT(INOUT)     ::  ND
      INTEGER, INTENT(INOUT)     ::  JPOS(0:N_SCAT_M,3)
      INTEGER, INTENT(IN)        ::  ITYP,IATL,I_CP
      INTEGER, INTENT(IN)        ::  JE
!
      INTEGER                    ::  N_TYP,I_ABS,IEULER
      INTEGER                    ::  JTYP,NBTYPJ,JNUM,JATL
      INTEGER                    ::  KTYP,NBTYPK,KNUM,KATL
!
      REAL (WP), INTENT(INOUT)   ::  XR(N_SCAT_M),YR(N_SCAT_M)
      REAL (WP), INTENT(INOUT)   ::  ZR(N_SCAT_M),R(N_SCAT_M)
      REAL (WP), INTENT(IN)      ::  THMI,PHIMI,ZSURF
      REAL (WP), INTENT(INOUT)   ::  FREF
      REAL (WP), INTENT(INOUT)   ::  DIJ
!
      REAL (WP)                  ::  COSFWDI,COSBWDI
      REAL (WP)                  ::  COSFWDJ,COSBWDJ
      REAL (WP)                  ::  COSTHMIJ,COSTHIJK
      REAL (WP)                  ::  THIJ,PHIIJ,THJK,PHIJK
      REAL (WP)                  ::  AMIJ,BMIJ,CMIJ
      REAL (WP)                  ::  AIJK,BIJK,CIJK
!
      REAL (WP)                  ::  COS,SQRT
!
      COMPLEX (WP), INTENT(IN)   ::  RHOMI
      COMPLEX (WP), INTENT(INOUT)::  PW(0:N_SCAT_M)
      COMPLEX (WP), INTENT(OUT)  ::  TAU(LINMAX,LINFMAX,NATP_M)
!
      COMPLEX (WP)               ::  RHOIJ,RHOJK
!
      IEULER = 1                                                    !
!
      IF(IFWD == 1)       COSFWDI = COS(RTHFWD(ITYP))               !
      IF(IBWD(ITYP) == 1) COSBWDI = COS(RTHBWD(ITYP))               !
!
!    I_CP = 0 : all open paths generated
!    I_CP = 1 : only closed paths generated
!
      IF((I_CP == 0) .OR. ( ND /= N_SCAT-1)) THEN                   !
        N_TYP = N_PROT                                              !
      ELSE                                                          !
        N_TYP = 1                                                   !
      END IF                                                        !
!
      DO JTYP = 1, N_TYP                                            !
       IF(IFWD == 1)       COSFWDJ = COS(RTHFWD(JTYP))              !
       IF(IBWD(JTYP) == 1) COSBWDJ = COS(RTHBWD(JTYP))              !
       ND = ND + 1                                                  !
!
!     I_ABS = 0 : the atom before the scatterer is not the absorber
!     I_ABS = 1 : the atom before the scatterer is the absorber
!     I_ABS = 2 : the atom after the scatterer is the absorber (XAS only)
!
       IF(ND == 1) THEN                                             !
         I_ABS = 1                                                  !
       ELSE                                                         !
         I_ABS = 0                                                  !
       END IF                                                       !
!
       IF((I_CP == 0) .OR. (ND /= N_SCAT)) THEN                     !
         NBTYPJ = NATYP(JTYP)                                       !
       ELSE                                                         !
         NBTYPJ = 1                                                 !
       END IF                                                       !
!
       DO JNUM = 1, NBTYPJ                                          !
          JATL = NCORR(JNUM,JTYP)                                   !
!
          IF(JATL == IATL) GO TO 12                                 !
          IF(INEW_AT(JATL) == 1) GO TO 12                           !
!
          XR(ND) = SYM_AT(1,JATL) - SYM_AT(1,IATL)                  !
          YR(ND) = SYM_AT(2,JATL) - SYM_AT(2,IATL)                  !
          ZR(ND) = SYM_AT(3,JATL) - SYM_AT(3,IATL)                  !
          R(ND)  = SQRT( XR(ND) * XR(ND) +                        & !
                         YR(ND) * YR(ND) +                        & !
                         ZR(ND) * ZR(ND) )                          !
          DIJ    = DIJ + R(ND)                                      !
!
          IF((ILE == 1) .AND. (DIJ > RLE)) IT(ND-1) = 1             !
          IF((IT(ND-1) == 1) .AND. (ND > 1)) GO TO 42               !
!
          JPOS(ND,1) = JTYP                                         !
          JPOS(ND,2) = JNUM                                         !
          JPOS(ND,3) = JATL                                         !
          NPATH(ND) = NPATH(ND) + ONE                               !
!
          IF(ND > 1) THEN                                           !
            COSTHMIJ=( XR(ND) * XR(ND-1) +                        & !
                       YR(ND) * YR(ND-1) +                        & !
                       ZR(ND) * ZR(ND-1) ) /                      & !
                       (R(ND) * R(ND-1))                            !
!
            IF(IFWD == 1) THEN                                      !
!
!  Forward scattering filter switched on
!
              CALL FWD_FILTER(ND,ITYP,COSFWDI,COSBWDI,COSTHMIJ)     !
!
            END IF                                                  !
!
          END IF                                                    !
!
          IF((IT(ND-1) == 1) .AND. (ND > 1)) GO TO 42               !
!
          RHOIJ = VK(JE) * R(ND)                                    !
          CALL CALC_THETA_PHI(XR(ND),YR(ND),ZR(ND),THIJ,PHIIJ)      !
!
          IF((ND > 1) .AND. ((ND-1) < N_SCAT)) THEN                 !
            IF(IDWSPH == 1) THEN                                    !
              DW(ND-1) = ONE                                        !
            ELSE                                                    !
              DW(ND-1) = EXP( - VK2(JE) * UJ2(ITYP) *              & !
                               (ONE - COSTHMIJ) )                   !
            END IF                                                  !
          END IF                                                    !
          IF(ND == 1) THEN                                          !
            RHO01 = RHOIJ                                           !
            TH01  = THIJ                                            !
            PHI01 = PHIIJ                                           !
            CALL DJMN2(TH01,RLM01,LF2,2)                            !
            GO TO 30                                                !
          END IF                                                    !
!
          IF(IPW == 1) THEN                                         !
!
!  Plane wave filter switched on
!
            CALL PW_FILTER(JE,JPOS,JTYP,ITYP,ND,THIJ,COSTHMIJ,    & !
                           FREF,R,PW)                               !
!
          END IF                                                    !
!
          IF((IT(ND-1) == 1) .OR. (IT(ND) == 1)) GO TO 42           !
!
          CALL EULER(THIJ,PHIIJ,THMI,PHIMI,AMIJ,BMIJ,CMIJ,IEULER)   !
          CALL MATDIF(NO,ND-1,LF2,ITYP,JTYP,JE,I_ABS,I_BASIS,     & !
                      AMIJ,BMIJ,CMIJ,RHOMI,RHOIJ)                   !
  30      CEX(ND)   = EXP(IC * RHOIJ) / R(ND)                       !
          CEXDW(ND) = CEX(ND) * DW(ND-1)                            !
!
          IF((IJ == 1) .OR. (ND == NCUT)) THEN                      !
            IF((I_CP == 0) .OR. (JATL == 1) .AND. (JNUM == 1))) THEN!
              CALL PATHOP(JPOS,ND,JE,I_CP,RHO01,PHI01,RHOIJ,      & !
                          THIJ,PHIIJ,FREF,IJ,DIJ,TAU)               !
              NPATH2(ND) = NPATH2(ND) + ONE                         !
            END IF                                                  !
          END IF                                                    !
          IF(ND == N_SCAT) GO TO 42                                 !
          I_ABS = 0                                                 !
!
          IF((I_CP == 0) .OR. (ND /= N_SCAT-1)) THEN                !
            N_TYP = N_PROT                                          !
          ELSE                                                      !
            N_TYP = 1                                               !
          END IF                                                    !
!
          DO KTYP =1, N_TYP                                         !
            ND = ND + 1                                             !
            IF(ND > N_SCAT) GO TO 20                                !
!
            IF((I_CP == 0) .OR. (ND /= N_SCAT)) THEN                !
              NBTYPK = NATYP(KTYP)                                  !
            ELSE                                                    !
              NBTYPK = 1                                            !
            END IF                                                  !
!
            DO KNUM = 1, NBTYPK                                     !
              KATL = NCORR(KNUM,KTYP)                               !
!
               IF(INEW_AT(KATL) == 1) GO TO 22                      !
               IF(KATL == JATL) GO TO 22                            !
!
              JPOS(ND,1) = KTYP                                     !
              JPOS(ND,2) = KNUM                                     !
              JPOS(ND,3) = KATL                                     !
              XR(ND)     = SYM_AT(1,KATL) - SYM_AT(1,JATL)          !
              YR(ND)     = SYM_AT(2,KATL) - SYM_AT(2,JATL)          !
              ZR(ND)     = SYM_AT(3,KATL) - SYM_AT(3,JATL)          !
              R(ND)      = SQRT( XR(ND) * XR(ND) +                & !
                                 YR(ND) * YR(ND) +                & !
                                 ZR(ND) * ZR(ND) )                  !
              DIJ        = DIJ + R(ND)                              !
!
              IF((ILE == 1) .AND. (DIJ > RLE)) IT(ND-1) = 1         !
              IF(IT(ND-1) == 1) GO TO 32                            !
!
              RHOJK     = R(ND) * VK(JE)                            !
              NPATH(ND) = NPATH(ND) + ONE                           !
              COSTHIJK  = ( XR(ND) * XR(ND-1) +                   & !
                            YR(ND) * YR(ND-1) +                   & !
                            ZR(ND) * ZR(ND-1) ) /                 & !
                          (R(ND)*R(ND-1))                           !
!
              IF(IFWD == 1) THEN                                    !
!
!  Forward scattering filter switched on
!
                CALL FWD_FILTER(ND,JTYP,COSFWDJ,COSBWDJ,COSTHIJK)   !
!
              END IF                                                !
!
              IF(IT(ND-1) == 1) GO TO 32                            !
!
              CALL CALC_THETA_PHI(XR(ND),YR(ND),ZR(ND),THJK,PHIJK)  !
!
              IF(ND-1 < N_SCAT) THEN                                !
                IF(IDWSPH == 1) THEN                                !
                  DW(ND-1) = ONE                                    !
                ELSE                                                !
                  DW(ND-1) = EXP( - VK2(JE) * UJ2(JTYP) *         & !
                                  (ONE - COSTHIJK) )                !
                END IF                                              !
              END IF                                                !
              IF(IPW == 1) THEN                                     !
!
!  Plane wave filter switched on
!
                CALL PW_FILTER(JE,JPOS,KTYP,JTYP,ND,THJK,         & !
                               COSTHIJK,FREF,R,PW)                  !
!
              END IF                                                !
!
              IF((IT(ND-1) == 1) .OR. (IT(ND) == 1)) GO TO 32       !
!
              CALL EULER(THJK,PHIJK,THIJ,PHIIJ,                   & !
                         AIJK,BIJK,CIJK,IEULER)                     !
              IF((I_CP == 1) .AND. (ND == N_SCAT)) I_ABS=2          !
              CALL MATDIF(NO,ND-1,LF2,JTYP,KTYP,JE,I_ABS,I_BASIS, & !
                          AIJK,BIJK,CIJK,RHOIJ,RHOJK)               !
              CEX(ND)   = EXP(IC * RHOJK) / R(ND)                   !
              CEXDW(ND) = CEX(ND) * DW(ND-1)                        !
              IF((IJ == 1).OR.(ND == NCUT)) THEN                    !
                IF((I_CP == 0) .OR. (KATL == 1) .AND.            & !
                                    (KNUM == 1)) THEN               !
                  CALL PATHOP(JPOS,ND,JE,I_CP,RHO01,PHI01,        & !
                              RHOJK,THJK,PHIJK,FREF,IJ,DIJ,TAU)     !
                  NPATH2(ND) = NPATH2(ND) + ONE                     !
                END IF                                              !
              END IF                                                !
              IF(ND == N_SCAT) GO TO 32                             !
              CALL FINDPATHS(ND,KTYP,KATL,I_CP,R,XR,YR,ZR,RHOJK,  & !
                             THJK,PHIJK,ZSURF,JPOS,PW,JE,FREF,    & !
                             DIJ,TAU)                               !
  32          DIJ = DIJ - R(ND)                                     !
!
  22          IF(IN(ND-1) == 1) NTHOF = NTHOF - 1                   !
!
              IT(ND-1) = 0                                          !
              IN(ND-1) = 0                                          !
            END DO                                                  !
  20        CONTINUE                                                !
            ND = ND - 1                                             !
          END DO                                                    !
  42      DIJ = DIJ - R(ND)                                         !
!
  12      IF(ND > 1) THEN                                           !
            IF(IN(ND-1) == 1) NTHOF = NTHOF - 1                     !
            IT(ND-1) = 0                                            !
            IN(ND-1) = 0                                            !
          END IF                                                    !
!
       END DO                                                       !
       ND = ND - 1                                                  !
      END DO                                                        !
!
      END SUBROUTINE FINDPATHS
!
!=======================================================================
!
      SUBROUTINE FWD_FILTER(ND,ITYP,COSFWDI,COSBWDI,COSTHMIJ)
!
!  This subroutine computes the forward scattering filter. It compares the
!    scattering angle with the forward/backward maximum scattering angle
!    as defined in the input data file.
!
!
!  Input variables :
!
!                       ND        :  index of the atom in the path
!                       ITYP      :  index of the atom I just before the atom J in the path
!                       COSFWDI   :  cosine of the maximum forward scattering angle
!                                    authorised on atom I
!                       COSBWDI   :  cosine of the maximum backward scattering angle
!                                    authorised on atom I
!                       COSTHMIJ  :  cosine of the scattering angle between directions
!                                    (KI) and (IJ) where K is the atom before I along the path
!
!
!  Output variable stored in module TESTPA :
!
!                       IT(ND-1)  :  switch array for selecting the path of scattering order ND-1
!                                     for rejection or not
!                                     ---> 0 : path computed in spherical waves
!                                     ---> 1 : path rejected for computation
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 14 May 2021
!
!
      USE REAL_NUMBERS,          ONLY : ZERO
!
      USE CURRENT_BEAM
      USE PATH_INFO
      USE TESTPA
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)        ::  ND,ITYP

      REAL (WP), INTENT(IN)      ::  COSFWDI,COSBWDI,COSTHMIJ
!
      REAL (WP)                  ::  SSMALL
!
      DATA SSMALL     / 0.0001E0_WP /

       IF(IBWD(ITYP) == 0) THEN                                     !
         IF(COSTHMIJ < COSFWDI) THEN                                !
           NTHOF    = NTHOF + 1                                     !
           IN(ND-1) = 1                                             !
           IF(NTHOF > NTHOUT) THEN                                  !
             IT(ND-1) = 1                                           !
           END IF                                                   !
         END IF                                                     !
       ELSE IF(IBWD(ITYP) == 1) THEN                                !
         IF( (COSTHMIJ > COSBWDI) .AND.                           & !
                     (COSTHMIJ < - SSMALL) ) THEN                   !
           NTHOF    = NTHOF + 1                                     !
           IN(ND-1) = 1                                             !
           IF(NTHOF > NTHOUT) THEN                                  !
             IT(ND-1) = 1                                           !
           END IF                                                   !
         END IF                                                     !
         IF((COSTHMIJ < COSFWDI) .AND. (COSTHMIJ >= ZERO)) THEN     !
           NTHOF    = NTHOF + 1                                     !
           IN(ND-1) = 1                                             !
           IF(NTHOF > NTHOUT) THEN                                  !
             IT(ND-1) = 1                                           !
           END IF                                                   !
         END IF                                                     !
      END IF                                                        !
!
      END SUBROUTINE FWD_FILTER
!
!=======================================================================
!
      SUBROUTINE CALC_THETA_PHI(X,Y,Z,THETA,PHI)
!
!   This subroutine computes the (theta,phi) angles from the
!     direction (X,Y,Z) of a vector
!
!  Input variables :
!
!                       X         :  x
!                       Y         :  y
!                       Z         :  z
!
!
!
!  Output variable
!
!                       THETA     :  theta
!                       PHI       :  phi
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified :  2 Jun 2021
!
!
      USE REAL_NUMBERS,          ONLY : ONE
!
      USE CURRENT_BEAM
      USE VECTOR
!
      IMPLICIT NONE
!
      REAL (WP), INTENT(IN)      ::  X,Y,Z
      REAL (WP), INTENT(OUT)     ::  THETA,PHI
!
      REAL (WP)                  ::  RR,TH,PH
!
!
      CALL XYZ_TO_R_THETA_PHI(X,Y,Z,RR,TH,PH)                       !
!
      THETA = TH                                                    !
      PHI   = PH                                                    !
!
      END SUBROUTINE CALC_THETA_PHI
!
!=======================================================================
!
      SUBROUTINE PW_FILTER(JE,JPOS,JTYP,ITYP,ND,THIJ,COSTHMIJ,     &
                           FREF,R,PW)
!
!  This subroutine computes the plane wave filter. It compares the amplitude
!    of a multiple scattering path computed in plane waves to a reference value FREF
!    which can be modulated by PCTINT for further flexibility
!
!
!  Input variables :
!
!                       JE        :  energy point
!                       JPOS      :  array containing the list of atoms along a path
!                       JTYP      :  index of the atom J just generated
!                       ITYP      :  index of the atom I just before the atom J in the path
!                       ND        :  index of the atom in the path
!                       THIJ      :  polar angle of the direction (I,J)
!                       COSTHMIJ  :  cosine of the scattering angle between directions
!                                    (MI) and (IJ) where M is the atom before I along the path
!                       FREF      :  reference intensity computed from subroutine CALL PW_FILTER_INI
!                       R         :  distance between atoms I and J
!
!
!  Output variables :
!
!                       PW        :  plane wave amplitude along the path
!
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 19 May 2021
!
!
      USE REAL_NUMBERS,           ONLY : ZERO
!
      USE CURRENT_BEAM
      USE CURRENT_T_MATRIX
      USE INIT_L
      USE TESTPA
      USE TESTPB
      USE TEMP_STORAGE_1
!
      USE SPHERICAL_HARMONICS,    ONLY : HARSPH2
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)        ::  JE,JTYP,ITYP,ND
      INTEGER, INTENT(IN)        ::  JPOS(0:N_SCAT_M,3)
!
      INTEGER                    ::  LMJ,LJ,LF
!
      REAL (WP), INTENT(IN)      ::  THIJ,COSTHMIJ,FREF
      REAL (WP), INTENT(IN)      ::  R(N_SCAT_M)
!
      REAL (WP)                  ::  XCOMP,PI4
      REAL (WP)                  ::  XMAXT
!
      REAL (WP)                  ::  MAX,ABS
!
      COMPLEX (WP), INTENT(OUT)  ::  PW(0:N_SCAT_M)
!
      COMPLEX (WP)               ::  PW1,PWI,FTHETA
      COMPLEX (WP)               ::  CTL,CTL2
      COMPLEX (WP)               ::  YLM1(0:NL_M,-NL_M:NL_M)
      COMPLEX (WP)               ::  YLM2(0:NL_M,-NL_M:NL_M)
!
      DATA XCOMP  / 1.0E-10_WP /
      DATA PI4    / 12.5663706143591729538505735331180E0_WP /          ! 4 pi
!
      CALL PW_SF(COSTHMIJ,JPOS(ND-1,1),JE,FTHETA)                   !
!
      PWI    = FTHETA * DW(ND-1) / R(ND)                            !
      PW(ND) = PW(ND-1) * PWI                                       !
      CTL2   = PI4 * PW(ND) * CEX(1) / VK(JE)                       !
      LMJ    = LMAX(ITYP,JE)                                        !
      IF(ND > NCUT) THEN                                            !
        IT(ND) = 1                                                  !
      ELSE                                                          !
        IT(ND) = 0                                                  !
      END IF                                                        !
!
      CALL HARSPH2(NL_M,TH01,YLM1,LF2)                              !
      CALL HARSPH2(NL_M,THIJ,YLM2,LMJ)                              !
!
      XMAXT = ZERO                                                  !
      DO LJ = 0, LMJ                                                !
        CTL = CTL2 * TL(LJ,1,JTYP,JE) * YLM2(LJ,0)                  !
        DO LF = LF1 , LF2, ISTEP_LF                                 !
          PW1   = CTL * YLM1(LF,0) * TL(LF,1,1,JE)                  !
          XMAXT = MAX(XMAXT,ABS(PW1))                               !
        END DO                                                      !
      END DO                                                        !
      IF((PCTINT * FREF - XMAXT < - XCOMP) .AND. (ND > NCUT)) THEN  !
        IT(ND) = 0                                                  !
      END IF                                                        !
!
      END SUBROUTINE PW_FILTER
!
!=======================================================================
!
      SUBROUTINE PW_SF(COSTH,JAT,JE,FTHETA)
!
!  This routine computes the plane wave scattering factor
!
!  Input variables :
!
!                       COSTH     :  cosine of the scattering angle
!                       JAT       :  index of the scattering atom
!                       JE        :  current energy point
!
!  Output variables :
!
!                       FTHETA    :  scattering amplitude
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 14 May 2021
!
!
      USE COMPLEX_NUMBERS,         ONLY : ZEROC
!
      USE CURRENT_T_MATRIX
!
      USE LEGENDRE_FUNCTIONS,      ONLY : POLLEG
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)        ::  JAT,JE
!
      INTEGER                    ::  NL,L
!
      REAL (WP), INTENT(IN)      ::  COSTH
!
      REAL (WP)                  ::  PL(0:150)
!
      REAL (WP)                  ::  FLOAT
!
      COMPLEX (WP), INTENT(OUT)  ::  FTHETA
!
      FTHETA = ZEROC                                                !
      NL     = LMAX(JAT,JE) + 1                                     !
!
!  Computing the Legendre polynomials
!
      CALL POLLEG(NL,COSTH,PL)                                      !
!
      DO L = 0, NL-1                                                !
!
        FTHETA = FTHETA + FLOAT(L + L + 1) * TL(L,1,JAT,JE) * PL(L) !
!
      END DO                                                        !
!
      FTHETA=FTHETA / VK(JE)                                        !
!
      END SUBROUTINE PW_SF
!
!=======================================================================
!
      SUBROUTINE MATDIF(NO,ND,LF,JTYP,KTYP,JE,I_ABS,I_BASIS,  &
                            A21,B21,C21,RHO1,RHO2)
!
!  This routine calculates the Rehr-Albers scattering matrix
!     F_{LAMBDA1,LAMBDA2}. The result is stored in the module
!     SCAT_MAT as F21(NSPIN2_M,NLAMBDA_M,NLAMBDA_M,N_SCAT_M).
!
!  The incoming value of rho = k * r is RHO1
!  The outgoing value of rho = k * r is RHO2
!  The scattering angle is characterized by the Euler angles (A21,B21,C21)
!
!  Input variables :
!
!                       NO        :  Rehr-Albers in dex of current beam
!                       ND        :  index of the scattering atom
!                       LF        :  final state anglar momentum l
!                       JTYP      :  type of atom J
!                       KTYP      :  type of atom K
!                       JE        :  current energy index
!                       I_ABS     :  absorber switch
!                                       I_ABS = 0  --> the atom before the scatterer is not the absorber
!                                       I_ABS = 1  --> the atom before the scatterer is the absorber
!                                       I_ABS = 2  --> the atom after the scatterer  is the absorber (XAS only)
!                       I_BASIS   :  approximation index
!                                       I_BASIS = 1  --> Hankel polynomials set to 1.0
!                                       I_BASIS > 1  --> standard calculation
!                       A21       :  \
!                       B21       :   > (alpha,beta,gamma) Euler scattering angles
!                       C21       :  /
!                       RHO1      :  k * r of incoming interatomic vector
!                       RHO2      :  k * r of outgoing interatomic vector
!
!  Output variables :
!
!                       F21       :  stored in module SCAT_MAT
!
!
!
!      This is the SYMMETRIZED version.
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 22 Jun 2021
!
!
!
      USE REAL_NUMBERS,              ONLY : ONE
      USE COMPLEX_NUMBERS
      USE PI_ETC,                    ONLY : PI
!
      USE CURRENT_LBDM_STORE
      USE CURRENT_LINLBD
      USE CURRENT_REHR_ALBERS
      USE CURRENT_T_MATRIX
      USE SCAT_MAT
      USE STORE_COEF,                ONLY : EXPF
!
      USE HANKEL_POLYNOMIALS
      USE WIGNER_ROTATIONS
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)        ::  NO,ND,LF,JTYP,KTYP,JE
      INTEGER, INTENT(IN)        ::  I_ABS,I_BASIS
!
      INTEGER                    ::  IB,LMJ,NN2
      INTEGER                    ::  NO1,NO2,NUMAX1,NUMAX2
      INTEGER                    ::  LAMBDA1,LAMBDA2
      INTEGER                    ::  LAMBDA1_1,LAMBDA1_2
      INTEGER                    ::  LAMBDA2_1,LAMBDA2_2
      INTEGER                    ::  LMIN,L
      INTEGER                    ::  NEQUAL
      INTEGER                    ::  NU1,NU2,MUMAX1,MUMAX2
      INTEGER                    ::  MU1,MU2,MU,MUP
!
      REAL (WP), INTENT(IN)      ::  A21,B21,C21
!
      REAL (WP)                  ::  SMALL
      REAL (WP)                  ::  RLM(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
      REAL (WP)                  ::  SIG1,SIG2,C1,CST1,CST2
!
      REAL (WP)                  ::  ABS,MAX,MIN
!
      COMPLEX (WP), INTENT(IN )  ::  RHO1,RHO2
!
      COMPLEX (WP)               ::  ONEOVK
      COMPLEX (WP)               ::  HLM1(0:NO_ST_M,0:NL_M-1)
      COMPLEX (WP)               ::  HLM2(0:NO_ST_M,0:NL_M-1)
      COMPLEX (WP)               ::  SL,SL_1_1,SL_1_2,SL_2_1,SL_2_2
      COMPLEX (WP)               ::  EXP1,EXP2,PROD1,PROD2
!
      DATA SMALL / 0.0001E0_WP /
!
      ONEOVK = ONE / VK(JE)                                         !
      IB     = 0                                                    !
      LMJ    = LMAX(JTYP,JE)                                        !
!
      IF(ABS(ABS(B21)-PI) < SMALL) IB = -1                          !
      IF(ABS(B21) < SMALL)         IB =  1                          !
      IF(NO == 8) THEN                                              !
        NN2 = LMAX(JTYP,JE) + 1                                     !
      ELSE                                                          !
        NN2 = NO                                                    !
      END IF                                                        !
!
!  NO is atom-dependent and is decreased with the rank of the scatterer
!    in the path when I_NO > 0. Here LAMBDA1 depends on the scatterer JTYP
!    while LAMBDA2 depends on the next atom (KTYP) in the path
!
      IF(I_NO == 0) THEN                                            !
        NO1 = N_RA(JTYP)                                            !
        NO2 = N_RA(KTYP)                                            !
      ELSE                                                          !
        NO1 = MAX(N_RA(JTYP) - (ND - 1) / I_NO,0)                   !
        NO2 = MAX(N_RA(KTYP) - ND       / I_NO,0)                   !
      END IF                                                        !
      IF(I_ABS == 0) THEN                                           !
        NUMAX1 = NO1 / 2                                            !
        NUMAX2 = NO2 / 2                                            !
      ELSE IF(I_ABS == 1) THEN                                      !
        NUMAX1 = MIN(LF,NO1 / 2)                                    !
        NUMAX2 = NO2 / 2                                            !
      ELSE IF(I_ABS == 2) THEN                                      !
        NUMAX1 = NO1 / 2                                            !
        NUMAX2 = MIN(LF,NO2 / 2)                                    !
      END IF                                                        !
      LBDM(1,ND) = (NO1 +1 ) * (NO1 + 2) / 2                        !
      LBDM(2,ND) = (NO2 +1 ) * (NO2 + 2) / 2                        !
!
      EXP2 = - EXP(- IC * A21)                                      !
      EXP1 =   EXP(- IC * C21)                                      !
!
      DO LAMBDA1 = 1, LBDMAX                                        !
        DO LAMBDA2 = 1, LBDMAX                                      !
          F21(1,LAMBDA2,LAMBDA1,ND) = ZEROC                         !
        END DO                                                      !
      END DO                                                        !
!
      IF(ABS(RHO1 - RHO2) > SMALL) THEN                             !
         CALL POLHAN(ISPHER,NUMAX1,LMJ,RHO1,HLM1)                   !
         CALL POLHAN(ISPHER,NN2,LMJ,RHO2,HLM2)                      !
         NEQUAL = 0                                                 !
      ELSE
         CALL POLHAN(ISPHER,NN2,LMJ,RHO1,HLM1)                      !
         NEQUAL = 1                                                 !
      END IF                                                        !
!
!  Calculation of the scattering matrix when the scattering angle
!              is different from 0 and pi
!
      IF(IB == 0) THEN                                              !
        CALL DJMN(B21,RLM,LMJ)                                      !
        DO NU1 = 0, NUMAX1                                          !
          MUMAX1 = NO1 - 2 * NU1                                    !
          IF(I_ABS == 1) MUMAX1 = MIN(LF - NU1,MUMAX1)              !
          DO NU2 = 0, NUMAX2                                        !
            MUMAX2 = NO2 - 2 * NU2                                  !
            PROD1  = ONEC                                           !
            DO MU1 = 0, MUMAX1                                      !
              LAMBDA1_1 = LBD(MU1,NU1)                              !
              IF(MU1 > 0) THEN                                      !
                LAMBDA1_2 = LBD(-MU1,NU1)                           !
                PROD1     = PROD1 * EXP1                            !
              ELSE                                                  !
                LAMBDA1_2 = LAMBDA1_1                               !
              END IF                                                !
              PROD2 = ONEC                                          !
              DO MU2 = 0, MUMAX2                                    !
                LAMBDA2_1 = LBD(MU2,NU2)                            !
                IF(MU2 > 0) THEN                                    !
                  LAMBDA2_2 = LBD(-MU2,NU2)                         !
                  PROD2 = PROD2 * EXP2                              !
                ELSE                                                !
                  LAMBDA2_2 = LAMBDA2_1                             !
                END IF                                              !
                LMIN = MAX(MU1,NU1,MU2+NU2)                         !
                SL_1_1 = ZEROC                                      !
                SL_1_2 = ZEROC                                      !
                SL_2_1 = ZEROC                                      !
                SL_2_2 = ZEROC                                      !
                DO L = LMIN, LMJ
                  IF(NEQUAL == 1) THEN
                    HLM2(MU2+NU2,L) = HLM1(MU2+NU2,L)
                  END IF
                  C1 = EXPF(MU1,L) / EXPF(MU2,L)
                  IF(ISPEED == 1) THEN
                    SL = FLOAT(L + L + 1) * C1 * TL(L,1,JTYP,JE)  & !
                                          * HLM1(NU1,L) *         & !
                                            HLM2(MU2+NU2,L)         !
                  ELSE
                    SL = FLOAT(L + L + 1) * C1 * TLT(L,1,JTYP,JE) & !
                                          * HLM1(NU1,L) *         & !
                                            HLM2(MU2+NU2,L)         !
                  END IF                                            !
                  IF(MU1 == 0) THEN                                 !
                    SL_1_1 = SL_1_1 + SL * RLM(-MU2,-MU1,L)         !
                    IF(MU2 > 0) THEN                                !
                      SL_2_1 = SL_2_1 + SL * RLM(MU2,-MU1,L)        !
                    END IF                                          !
                  ELSE IF(MU2 == 0) THEN                            !
                    SL_1_1 = SL_1_1 + SL * RLM(-MU2,-MU1,L)         !
                    IF(MU1 > 0) THEN                                !
                      SL_1_2 = SL_1_2 + SL * RLM(-MU2,MU1,L)        !
                    END IF                                          !
                  ELSE                                              !
                    SL_1_1 = SL_1_1 + SL * RLM(-MU2,-MU1,L)         !
                    SL_1_2 = SL_1_2 + SL * RLM(-MU2,MU1,L)          !
                    SL_2_1 = SL_2_1 + SL * RLM(MU2,-MU1,L)          !
                    SL_2_2 = SL_2_2 + SL * RLM(MU2,MU1,L)           !
                  END IF                                            !
                END DO                                              !
                IF(MU1 == 0) THEN                                   !
                  F21(1,LAMBDA2_1,LAMBDA1_1,ND) = SL_1_1 * PROD1 *& !
                                                  PROD2 * ONEOVK    !
                  IF(MU2 > 0) THEN
                    F21(1,LAMBDA2_2,LAMBDA1_1,ND) = SL_2_1 *      & !
                                                    PROD1 *       & !
                                                    ONEOVK / PROD2  !
                  END IF                                            !
                ELSE IF(MU2 == 0) THEN                              !
                  F21(1,LAMBDA2_1,LAMBDA1_1,ND) = SL_1_1 * PROD1 *& !
                                                  PROD2 * ONEOVK    !
                  IF(MU1 > 0) THEN                                  !
                    F21(1,LAMBDA2_1,LAMBDA1_2,ND) = SL_1_2 *      & !
                                                    ONEOVK /      & !
                                                   (PROD1 * PROD2)  !
                  END IF                                            !
                ELSE                                                !
                  F21(1,LAMBDA2_1,LAMBDA1_1,ND) = SL_1_1 *        & !
                                                  PROD1 * PROD2 * & !
                                                  ONEOVK            !
                  F21(1,LAMBDA2_2,LAMBDA1_1,ND) = SL_2_1 *        & !
                                                  PROD1 *         & !
                                                  ONEOVK / PROD2    !
                  F21(1,LAMBDA2_1,LAMBDA1_2,ND) = SL_1_2 *        & !
                                                  ONEOVK *        & !
                                                  PROD2 / PROD1     !
                  F21(1,LAMBDA2_2,LAMBDA1_2,ND) = SL_2_2 *        & !
                                                  ONEOVK /        & !
                                                  (PROD1 * PROD2)   !
                END IF                                              !
              END DO                                                !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
!
!  Calculation of the scattering matrix when the scattering angle
!     is equal to 0 (forward scattering) or pi (backscattering)
!
      ELSE                                                          !
        DO NU1 = 0, NUMAX1                                          !
          DO NU2 = 0, NUMAX2                                        !
            MUMAX1 = MIN(NO1 - 2 * NU1,NO1 - 2 * NU2)               !
            IF(I_ABS == 1) MUMAX1 = MIN(LF - NU1,MUMAX1)            !
            DO MU = 0, MUMAX1                                       !
              MUP       = MU * IB                                   !
              LAMBDA1_1 = LBD(MUP,NU1)                              !
              LAMBDA1_2 = LBD(-MUP,NU1)                             !
              LAMBDA2_1 = LBD(MU,NU2)                               !
              LAMBDA2_2 = LBD(-MU,NU2)                              !
              IF(MOD(MU,2) == 0) THEN                               !
                CST1 =   ONE                                        !
              ELSE                                                  !
                CST1 = - ONE                                        !
              END IF                                                !
              LMIN = MAX(NU1,MU + NU2)                              !
              SL = ZEROC                                            !
              DO L = LMIN, LMJ                                      !
                IF(NEQUAL == 1) THEN                                !
                   HLM2(MU+NU2,L) = HLM1(MU+NU2,L)                  !
                END IF                                              !
                IF(IB == -1) THEN                                   !
                  IF(MOD(L,2) == 0) THEN                            !
                    CST2 =   ONE                                    !
                  ELSE                                              !
                    CST2 = - ONE                                    !
                  END IF                                            !
                ELSE                                                !
                  CST2 = ONE                                        !
                END IF                                              !
                IF(ISPEED == 1) THEN                                !
                  SL = SL + FLOAT(L + L + 1) * CST2 *             & !
                             TL(L,1,JTYP,JE) * HLM1(NU1,L) *      & !
                             HLM2(MU+NU2,L)                         !
                ELSE
                  SL = SL + FLOAT(L + L + 1) * CST2 *             & !
                            TLT(L,1,JTYP,JE) * HLM1(NU1,L) *      & !
                            HLM2(MU+NU2,L)                          !
                END IF                                              !
              END DO                                                !
              F21(1,LAMBDA2_1,LAMBDA1_1,ND) = SL * CST1 * ONEOVK    !
              IF(MU > 0) THEN                                       !
                F21(1,LAMBDA2_2,LAMBDA1_2,ND) = SL * CST1 * ONEOVK  !
              END IF                                                !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
      END SUBROUTINE MATDIF
!
END MODULE CALC_TAU0J_SY
