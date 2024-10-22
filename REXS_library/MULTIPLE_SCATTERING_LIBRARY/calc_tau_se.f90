!
!======================================================================
!
MODULE CALC_TAU0J_SE
!
!  This module provides the subroutine to compute
!    Tau^{0 j} or Tau^{j 0} in the Rehr-Albers series expansion case
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE CALC_TAU_SE(TAU_TYPE,IATL,JE,TAU)
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
!  This is the series expansion version
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified :  7 Jun 2021
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
      CALL INIT_TAU_M(JE,N_PROT,NATYP,LMAX,TAU)                     !
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
!  Read formats:
!
 100  FORMAT(9X,I2,2X,E12.6,7X,E12.6,1X,F6.3,1X,10(I3,2X))
!
!  Write formats:
!
 200  FORMAT(/////,18X,'THE ',I3,' MORE INTENSE PATHS BY DECREASING', &
            ' ORDER :',/,24X,'(THE LENGTH IS GIVEN IN UNITS ',        &
            'OF A)')
 201  FORMAT(//,1X,'RANK',1X,'ORDER',4X,'No PATH',3X,                 &
             'INTENSITY',3X,'LENGTH',4X,'ABS',3X,                     &
             'ORDER OF THE SCATTERERS',/)
 202  FORMAT(I3,4X,I2,1X,E12.6,3X,E12.6,2X,F6.3,4X,I2,4X,10(I3,1X))
 203  FORMAT(I3,4X,I2,1X,I10,3X,E12.6,2X,F6.3,4X,I2,4X,10(I3,1X))
 204  FORMAT(/////,25X,' PATHS USED IN THE CALCULATION : ',           &
             /,24X,'(THE LENGTH IS GIVEN IN UNITS OF A)')
 205  FORMAT(//,1X,'JDIF',4X,'No OF THE PATH',2X,                     &
             'INTENSITY',3X,'LENGTH',4X,'ABSORBER',2X,                &
             'ORDER OF THE SCATTERERS',/)
 206  FORMAT(2X,I2,2X,I10,7X,E12.6,2X,F6.3,7X,I2,7X,10(I3,2X))
 207  FORMAT('   ')
 208  FORMAT(2X,I2,2X,E12.6,7X,E12.6,2X,F6.3,7X,I2,7X,10(I3,2X))
 209  FORMAT(///)
 210  FORMAT(10X,'<===== NUMBER OF PATHS TOO LARGE FOR PRINTING ',    &
             '=====>')
 211  FORMAT(//,2X,'ORDER ',I2,'  TOTAL NUMBER OF PATHS     : ',F15.1,&
                   /,10X,'  EFFECTIVE NUMBER OF PATHS : ',F15.1,      &
                   /,10X,'  MINIMAL INTENSITY         : ',E12.6,      &
                   2X,'No OF THE PATH : ',F15.1,                      &
                   /,10X,'  MAXIMAL INTENSITY         : ',E12.6,      &
                   2X,'No OF THE PATH : ',F15.1)
 212  FORMAT(//,2X,'ORDER ',I2,'  TOTAL NUMBER OF PATHS     : ',I10,  &
                   /,10X,'  EFFECTIVE NUMBER OF PATHS : ',I10,        &
                   /,10X,'  MINIMAL INTENSITY         : ',E12.6,      &
                   2X,'No OF THE PATH : ',I10,                        &
                   /,10X,'  MAXIMAL INTENSITY         : ',E12.6,      &
                   2X,'No OF THE PATH : ',I10)
 213  FORMAT(10X,'  CUT-OFF INTENSITY       : ',E12.6)
!
      END SUBROUTINE CALC_TAU_SE
!
!======================================================================
!
      SUBROUTINE INIT_TAU_M(JE,N_PROT,NATYP,LMAX,TAU)
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
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 14 May 2021
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
!
      INTEGER                   ::  JATL,JTYP,JNUM
      INTEGER                   ::  NBTYP,LMJ
      INTEGER                   ::  LJ,ILJ,MJ,INDJ
      INTEGER                   ::  LF,ILF,MF,INDF
!
      COMPLEX (WP), INTENT(OUT) ::  TAU(LINMAX,LINFMAX,NATCLU_M)
!
      JATL = 0                                                      !
!
      DO JTYP = 1, N_PROT                                           !
!
        NBTYP = NATYP(JTYP)                                         !
        LMJ   = LMAX(JTYP,JE)                                       !
!
        DO JNUM = 1, NBTYP                                          !
!
          JATL = JATL + 1                                           !
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
                  TAU(INDJ,INDF,JATL) = ZEROC                       !
!
                END DO                                              !
!
              END DO                                                !
!
            END DO                                                  !
!
          END DO                                                    !
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
!   Author :  D. Sébilleau
!
!                                          Last modified : 19 May 2021
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
      CALL INIT_TAU_M(JE,N_PROT,NATYP,LMAX,TAU)                     !
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
!
!   Author :  D. Sébilleau
!
!                                          Last modified :  2 Jun 2021
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
      COMPLEX (WP), INTENT(OUT)  ::  TAU(LINMAX,LINFMAX,NATCLU_M)
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
            IF((I_CP == 0) .OR. (JATL == 1)) THEN                   !
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
              IF(KATL == JATL) GO TO 22                             !
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
                IF((I_CP == 0) .OR. (KATL == 1)) THEN               !
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
!   Author :  D. Sébilleau
!
!                                          Last modified : 19 May 2021
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
      COMPLEX (WP)               ::  SL,SL_2_1,SL_2_2
      COMPLEX (WP)               ::  EXP1,EXP2,PROD1,PROD2
!
      DATA SMALL / 0.0001E0_WP /
!
      ONEOVK = ONE / VK(JE)                                         !
      IB     = 0                                                    !
      LMJ    = LMAX(JTYP,JE)                                        !
      IF(ABS(ABS(B21) - PI) < SMALL) IB = - 1                       !
      IF(ABS(B21) < SMALL)           IB =   1                       !
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
        NO1 = MAX(N_RA(JTYP) - (ND-1)/ I_NO,0)                      !
        NO2 = MAX(N_RA(KTYP) - ND / I_NO,0)                         !
      END IF                                                        !
      IF(I_ABS == 0) THEN                                           !
        NUMAX1 = NO1 / 2                                            !
        NUMAX2 = NO2 / 2                                            !
      ELSE IF(I_ABS == 1) THEN                                      !
        NUMAX1 = MIN(LF,NO1/2)                                      !
        NUMAX2 = NO2 / 2                                            !
      ELSE IF(I_ABS == 2) THEN                                      !
        NUMAX1 = NO1 / 2                                            !
        NUMAX2 = MIN(LF,NO2/2)                                      !
      END IF                                                        !
      LBDM(1,ND) = (NO1+1) * (NO1+2) / 2                            !
      LBDM(2,ND) = (NO2+1) * (NO2+2) / 2                            !
!
      EXP2 = - EXP(- IC * A21)                                      !
      EXP1 =   EXP(- IC * C21)                                      !
!
      DO LAMBDA1 = 1 , LBDMAX                                       !
        DO LAMBDA2 = 1, LBDMAX                                      !
          F21(1,LAMBDA2,LAMBDA1,ND) = ZEROC                         !
        END DO                                                      !
      END DO                                                        !
!
      IF(ABS(RHO1-RHO2) > SMALL) THEN                               !
         CALL POLHAN(I_BASIS,NUMAX1,LMJ,RHO1,HLM1)                  !
         CALL POLHAN(I_BASIS,NN2,LMJ,RHO2,HLM2)                     !
         NEQUAL = 0                                                 !
      ELSE                                                          !
         CALL POLHAN(I_BASIS,NN2,LMJ,RHO1,HLM1)                     !
         NEQUAL = 1                                                 !
      END IF                                                        !
!
!  Calculation of the scattering matrix when the scattering angle
!              is different from 0 and pi
!
      IF(IB == 0) THEN                                              !
        CALL DJMN2(B21,RLM,LMJ,0)                                   !
        DO NU1 = 0, NUMAX1                                          !
          MUMAX1  = NO1 - 2 * NU1                                   !
          IF(I_ABS == 1) MUMAX1 = MIN(LF-NU1,MUMAX1)                !
          DO NU2 = 0, NUMAX2                                        !
            MUMAX2 = NO2 - 2 * NU2                                  !
!
!  Case MU1 = 0
!
            LAMBDA1 = LBD(0,NU1)                                    !
!
!    Case MU2 = 0
!
            LAMBDA2 = LBD(0,NU2)                                    !
            LMIN = MAX(NU1,NU2)                                     !
            SL = ZEROC                                              !
            DO L = LMIN, LMJ                                        !
              IF(NEQUAL == 1) THEN                                  !
                HLM2(NU2,L) = HLM1(NU2,L)                           !
              END IF                                                !
              SL = SL + FLOAT(L+L+1) * RLM(0,0,L) *               & !
                        TL(L,1,JTYP,JE)           *               & !
                        HLM1(NU1,L) * HLM2(NU2,L)                   !
            END DO                                                  !
            F21(1,LAMBDA2,LAMBDA1,ND) = SL * ONEOVK                 !
!
!    Case MU2 > 0
!
            PROD2 = ONEC                                            !
            SIG2  = ONE                                             !
            DO MU2 = 1, MUMAX2                                      !
              LAMBDA2_1 = LBD(MU2,NU2)                              !
              LAMBDA2_2 = LBD(-MU2,NU2)                             !
              PROD2     = PROD2 * EXP2                              !
              SIG2      = - SIG2                                    !
              LMIN      = MAX(NU1,MU2+NU2)                          !
              SL        = ZEROC                                     !
              DO L = LMIN, LMJ                                      !
                IF(NEQUAL == 1) THEN                                !
                  HLM2(MU2+NU2,L) = HLM1(MU2+NU2,L)                 !
                END IF                                              !
                C1 = EXPF(0,L) / EXPF(MU2,L)                        !
                SL = SL + FLOAT(L+L+1) * RLM(MU2,0,L) * C1 *      & !
                          TL(L,1,JTYP,JE)                  *      & !
                          HLM1(NU1,L) * HLM2(MU2+NU2,L)             !
              END DO                                                !
              F21(1,LAMBDA2_1,LAMBDA1,ND) = SL * PROD2 * ONEOVK * SIG2
              F21(1,LAMBDA2_2,LAMBDA1,ND) = SL * ONEOVK / PROD2     !
            END DO                                                  !
!
!  Case MU1 > 0
!
            PROD1 = ONEC                                            !
            SIG1  = ONE                                             !
            DO MU1 = 1, MUMAX1                                      !
              LAMBDA1_1 = LBD(MU1,NU1)                              !
              LAMBDA1_2 = LBD(-MU1,NU1)                             !
              PROD1     = PROD1 * EXP1                              !
              SIG1      = - SIG1                                    !
!
!    Case MU2 = 0
!
              LAMBDA2 = LBD(0,NU2)                                  !
              LMIN    = MAX(MU1,NU1,NU2)                            !
              SL      = ZEROC                                       !
              DO L = LMIN, LMJ                                      !
                IF(NEQUAL == 1) THEN                                !
                  HLM2(NU2,L) = HLM1(NU2,L)                         !
                END IF
                C1 = EXPF(MU1,L) / EXPF(0,L)                        !
                SL = SL + FLOAT(L + L + 1) * RLM(0,MU1,L) * C1 *  & !
                          TL(L,1,JTYP,JE)                      *  & !
                          HLM1(NU1,L) * HLM2(NU2,L)                 !
              END DO                                                !
              F21(1,LAMBDA2,LAMBDA1_1,ND) = SL * PROD1 * ONEOVK * SIG1
              F21(1,LAMBDA2,LAMBDA1_2,ND) = SL * ONEOVK / PROD1     !
!
!    Case MU2 > 0
!
              PROD2 = ONEC                                          !
              SIG2  = SIG1                                          !
              DO MU2 = 1, MUMAX2                                    !
                LAMBDA2_1 = LBD(MU2,NU2)                            !
                LAMBDA2_2 = LBD(-MU2,NU2)                           !
                PROD2     = PROD2 * EXP2                            !
                SIG2      = - SIG2                                  !
                LMIN      = MAX(MU1,NU1,MU2+NU2)                    !
                SL_2_1    = ZEROC                                   !
                SL_2_2    = ZEROC                                   !
                DO L = LMIN, LMJ                                    !
                  IF(NEQUAL == 1) THEN                              !
                    HLM2(MU2+NU2,L) = HLM1(MU2+NU2,L)               !
                  END IF                                            !
                  C1 = EXPF(MU1,L) / EXPF(MU2,L)                    !
                  SL = FLOAT(L + L + 1) * C1 *                    & !
                       TL(L,1,JTYP,JE)       *                    & !
                       HLM1(NU1,L) * HLM2(MU2+NU2,L)                !
                  SL_2_1 = SL_2_1 + SL * RLM(MU2,-MU1,L)            !
                  SL_2_2 = SL_2_2 + SL * RLM(MU2,MU1,L)             !
                END DO
                F21(1,LAMBDA2_1,LAMBDA1_1,ND) = SL_2_2 * PROD1 *  & !
                                                PROD2 * ONEOVK * SIG2
                F21(1,LAMBDA2_2,LAMBDA1_1,ND) = SL_2_1 * PROD1 *  & !
                                                ONEOVK / PROD2      !
                F21(1,LAMBDA2_1,LAMBDA1_2,ND) = SL_2_1 * ONEOVK * & !
                                                PROD2 * SIG2 /    & !
                                                PROD1               !
                F21(1,LAMBDA2_2,LAMBDA1_2,ND) = SL_2_2 * ONEOVK / & !
                                                (PROD1 * PROD2)     !

              END DO                                                !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
!
!  Calculation of the scattering matrix when the scattering angle
!     is equal to 0 (forward scattering) or pi (backscattering)
!
      ELSE IF(IB == 1) THEN                                         !
        DO NU1 = 0, NUMAX1                                          !
          DO NU2 = 0, NUMAX2                                        !
            MUMAX1 = MIN(NO1-2*NU1,NO1-2*NU2)                       !
            IF(I_ABS == 1) MUMAX1 = MIN(LF-NU1,MUMAX1)              !
!
!  Case MU = 0
!
            LAMBDA1 = LBD(0,NU1)                                    !
            LAMBDA2 = LBD(0,NU2)                                    !
            LMIN    = MAX(NU1,NU2)                                  !
            SL      = ZEROC                                         !
            DO L = LMIN, LMJ                                        !
              IF(NEQUAL == 1) THEN                                  !
                 HLM2(NU2,L) = HLM1(NU2,L)                          !
              END IF                                                !
              SL = SL + FLOAT(L + L + 1)    *                     & !
                        TL(L,1,JTYP,JE)     *                     & !
                        HLM1(NU1,L) * HLM2(NU2,L)                   !
            END DO
            F21(1,LAMBDA2,LAMBDA1,ND) = SL * ONEOVK                 !
!
!  Case MU > 0
!
            CST1 = ONE                                              !
            DO MU = 1, MUMAX1                                       !
              LAMBDA1 = LBD(MU,NU2)                                 !
              LAMBDA2 = LBD(-MU,NU2)                                !
              CST1    = - CST1                                      !
              LMIN    = MAX(NU1,MU+NU2)                             !
              SL      = ZEROC                                       !
              DO L = LMIN, LMJ                                      !
                IF(NEQUAL == 1) THEN                                !
                   HLM2(MU+NU2,L) = HLM1(MU+NU2,L)                  !
                END IF                                              !
                SL = SL + FLOAT(L + L + 1) * CST1 *               & !
                         TL(L,1,JTYP,JE)         *                & !
                         HLM1(NU1,L) * HLM2(MU+NU2,L)               !
              END DO                                                !
              F21(1,LAMBDA1,LAMBDA1,ND) = SL * ONEOVK               !
              F21(1,LAMBDA2,LAMBDA2,ND) = SL * ONEOVK               !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
      ELSE IF(IB == -1) THEN                                        !
        DO NU1 = 0, NUMAX1                                          !
          DO NU2 = 0, NUMAX2                                        !
            MUMAX1 = MIN(NO1-2*NU1,NO1-2*NU2)                       !
            IF(I_ABS == 1) MUMAX1 = MIN(LF-NU1,MUMAX1)              !
!
!  Case MU = 0
!
            LAMBDA1 = LBD(0,NU1)                                    !
            LAMBDA2 = LBD(0,NU2)                                    !
            LMIN    = MAX(NU1,NU2)                                  !
            SL      = ZEROC                                         !
            DO L = LMIN, LMJ                                        !
              IF(NEQUAL == 1) THEN                                  !
                 HLM2(NU2,L) = HLM1(NU2,L)                          !
              END IF                                                !
              IF(MOD(L,2) == 0) THEN                                !
                CST2 =   ONE                                        !
              ELSE                                                  !
                CST2 = - ONE                                        !
              END IF                                                !
              SL = SL + FLOAT(L + L + 1) * CST2 *                 & !
                        TL(L,1,JTYP,JE)         *                 & !
                        HLM1(NU1,L) * HLM2(NU2,L)                   !
            END DO                                                  !
            F21(1,LAMBDA2,LAMBDA1,ND) = SL * ONEOVK                 !
!
!  Case MU > 0
!
            CST1 = ONE                                              !
            DO MU = 1, MUMAX1                                       !
              MUP       = - MU                                      !
              LAMBDA1_1 = LBD(MUP,NU1)                              !
              LAMBDA1_2 = LBD(-MUP,NU1)                             !
              LAMBDA2_1 = LBD(MU,NU2)                               !
              LAMBDA2_2 = LBD(-MU,NU2)                              !
              CST1      = - CST1                                    !
              LMIN      = MAX(NU1,MU+NU2)                           !
              SL        = ZEROC                                     !
              DO L = LMIN, LMJ                                      !
                IF(NEQUAL == 1) THEN                                !
                   HLM2(MU+NU2,L) = HLM1(MU+NU2,L)                  !
                END IF                                              !
                IF(MOD(L,2) == 0) THEN                              !
                  CST2 =   CST1                                     !
                ELSE                                                !
                  CST2 = - CST1                                     !
                END IF                                              !
                SL = SL + FLOAT(L + L + 1) * CST2 *               & !
                          TL(L,1,JTYP,JE)         *               & !
                          HLM1(NU1,L) * HLM2(MU+NU2,L)              !
              END DO                                                !
              F21(1,LAMBDA2_1,LAMBDA1_1,ND) = SL * ONEOVK           !
              F21(1,LAMBDA2_2,LAMBDA1_2,ND) = SL * ONEOVK           !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
      END SUBROUTINE MATDIF
!
!=======================================================================
!
      SUBROUTINE PATHOP(JPOS,JORDP,JE,I_CP,RHO01,PHI01,RHOIJ,      &
                        THIJ,PHIIJ,FREF,IJ,D,TAU)
!
!  This subroutine calculates the contribution of a given path to
!    the scattering path operator TAU.
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified :  3 Jun 2021
!
!
!
      USE REAL_NUMBERS,              ONLY : ZERO,ONE
      USE COMPLEX_NUMBERS
      USE PI_ETC,                    ONLY : PI
!
      USE CURRENT_BEAM
      USE EXTREMES
      USE INIT_L
      USE CURRENT_LBDM_STORE
      USE CURRENT_LINLBD
      USE CURRENT_REHR_ALBERS
      USE CURRENT_T_MATRIX
      USE LOCATE_MOD
      USE OUTUNITS
      USE PATH_INFO
      USE PRINT_PATH_INFO
      USE ROTATION
      USE SCAT_MAT
      USE TEMP_STORAGE_1
      USE TESTS
!
      USE HANKEL_POLYNOMIALS
      USE SORT1
      USE WIGNER_ROTATIONS
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)        ::  JE,I_CP,JORDP,IJ
      INTEGER, INTENT(IN)        ::  JPOS(0:N_SCAT_M,3)
!
      INTEGER                    ::  AMU1,JTYP,ITYP,JATL
      INTEGER                    ::  LMJ,NN,NN2
      INTEGER                    ::  NO1,NOJ,NUMX,NUMAXJ
      INTEGER                    ::  JSC,LBD1M1,LBD1M2,JPAT
      INTEGER                    ::  LAMBDA1,LAMBDA2,LAMBDA3,LAMBDAJ
      INTEGER                    ::  LBD2M,LBD3M
      INTEGER                    ::  LF,ILF,LJ,ILJ,INDJ
      INTEGER                    ::  NU1MAX1,NU1MAX,MU1MAX
      INTEGER                    ::  MU1,NU1
      INTEGER                    ::  MUJ,NUJ,MUJMAX,NUJMAX
      INTEGER                    ::  MJ,INDJP
      INTEGER                    ::  MF,INDF,INDFP
      INTEGER                    ::  KF,KD
      INTEGER                    ::  JMX
!
      REAL (WP), INTENT(IN)      ::  PHI01,THIJ,PHIIJ
      REAL (WP), INTENT(IN)      ::  D
      REAL (WP), INTENT(INOUT)   ::  FREF
!
      REAL (WP)                  ::  XCOMP,XMAX,FM1,FM2
      REAL (WP)                  ::  XINT,XINT1,XINT2,XINT3,XINT4
      REAL (WP)                  ::  RLMIJ(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
!
      COMPLEX (WP), INTENT(IN)   ::  RHO01,RHOIJ
      COMPLEX (WP), INTENT(OUT)  ::  TAU(LINMAX,LINFMAX,NATCLU_M)
!
      COMPLEX (WP)               ::  H(NLAMBDA_M,NLAMBDA_M)
      COMPLEX (WP)               ::  G(NLAMBDA_M,NLAMBDA_M)
      COMPLEX (WP)               ::  HLM01(0:NO_ST_M,0:NL_M-1)
      COMPLEX (WP)               ::  HLMIJ(0:NO_ST_M,0:NL_M-1)
      COMPLEX (WP)               ::  CF,CJ,OVK
      COMPLEX (WP)               ::  RLMF_0,RLMF_1
      COMPLEX (WP)               ::  SUM_NUJ_0,SUM_MUJ_0,SUM_NU1_0
      COMPLEX (WP)               ::  SUM_NUJ_1,SUM_MUJ_1,SUM_NU1_1
      COMPLEX (WP)               ::  SUM_NU1_2,SUM_NU1_3
      COMPLEX (WP)               ::  TL_J,EXP_J,EXP_F,SUM_1
      COMPLEX (WP)               ::  COEF,F
!
      COMPLEX (WP)               ::  CONJG
!
      DATA XCOMP / 1.0E-10_WP /
!
      OVK = ONE / VK(JE)                                            !
!
      IF(NPATHP > 0) THEN                                           !
        FM1  = FMIN(JORDP)                                          !
        XMAX = ZERO                                                 !
      END IF                                                        !
!
      EXP_J = EXP(- IC * (PHIIJ - PI))                              !
      EXP_F = EXP(  IC * PHI01)                                     !
!
      JTYP  = JPOS(JORDP,1)                                         !
      ITYP  = JPOS(1,1)                                             !
      JATL  = JPOS(JORDP,3)                                         !
!
      IF(I_CP == 0) THEN                                            !
        LMJ = LMAX(JTYP,JE)                                         !
      ELSE                                                          !
        LMJ = LF2                                                   !
      END IF                                                        !
      IF(NO == 8) THEN                                              !
        NN2 = LMJ + 1                                               !
      ELSE                                                          !
        NN2 = NO                                                    !
      END IF                                                        !
      IF(NO > LF2) THEN                                             !
        NN = LF2                                                    !
      ELSE                                                          !
        NN = NO                                                     !
      END IF                                                        !
!
!  NO is atom-dependent and is decreased with the rank of the scatterer
!    in the path when I_NO > 0 (except for the first scatterer ITYP for
!    which there is no such decrease)
!
      NO1 = N_RA(ITYP)                                              !
!
      IF(I_NO == 0) THEN                                            !
        IF(IJ == 1) THEN                                            !
          NOJ = N_RA(JTYP)                                          !
        ELSE                                                        !
          NOJ = 0                                                   !
        END IF                                                      !
      ELSE                                                          !
        IF(IJ == 1) THEN                                            !
          NOJ = MAX(N_RA(JTYP) - (JORDP-1) / I_NO,0)                !
        ELSE                                                        !
          NOJ = 0                                                   !
        END IF                                                      !
      END IF                                                        !
!
      NUMX   = NO1 / 2                                              !
      NUMAXJ = NOJ / 2                                              !
!
!   Calculation of the attenuation coefficients along the path
!
      COEF = CEX(1) * OVK                                           !
      DO JSC = 2, JORDP                                             !
        COEF = COEF * CEXDW(JSC)                                    !
      END DO                                                        !

!   Call of the subroutines used for the R-A termination matrix
!      This termination matrix is now merged into PATHOP
!
      CALL DJMN2(-THIJ,RLMIJ,LMJ,1)                                 !
      CALL POLHAN(I_BASIS,NN,LF2,RHO01,HLM01)                       !
      CALL POLHAN(I_BASIS,NN2,LMJ,RHOIJ,HLMIJ)                      !
!
      LBD1M1 = LBDM(1,1)                                            !
      LBD1M2 = LBDM(2,1)                                            !
!
!   Calculation of the L-independent part of TAU, called H
!
      IF(JORDP >= 3) THEN                                           !
        DO JPAT = 2, JORDP-1                                        !
          LBD2M = LBDM(1,JPAT)                                      !
          LBD3M = LBDM(2,JPAT)                                      !
          DO LAMBDA1 = 1 ,LBD1M1                                    !
            DO LAMBDA3 = 1, LBD3M                                   !
              SUM_1 = ZEROC                                         !
              DO LAMBDA2 = 1, LBD2M                                 !
                IF(JPAT > 2) THEN                                   !
                  SUM_1 = SUM_1 + H(LAMBDA2,LAMBDA1) *            & !
                                  F21(1,LAMBDA3,LAMBDA2,JPAT)       !
                ELSE                                                !
                  SUM_1 = SUM_1 + F21(1,LAMBDA2,LAMBDA1,1) *      & !
                                  F21(1,LAMBDA3,LAMBDA2,2)          !
                END IF                                              !
              END DO                                                !
              G(LAMBDA3,LAMBDA1) = SUM_1                            !
            END DO                                                  !
          END DO                                                    !
          DO LAMBDA1 = 1, LBD1M1                                    !
            DO LAMBDA2 = 1, LBD3M                                   !
              H(LAMBDA2,LAMBDA1) = G(LAMBDA2,LAMBDA1)               !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
      ELSE IF(JORDP == 2) THEN                                      !
        DO LAMBDA1 = 1, LBD1M1                                      !
          DO LAMBDA2 = 1, LBD1M2                                    !
            H(LAMBDA2,LAMBDA1) = F21(1,LAMBDA2,LAMBDA1,1)           !
          END DO                                                    !
        END DO                                                      !
      ELSE IF(JORDP == 1) THEN                                      !
        DO LAMBDA1 = 1 ,LBD1M1                                      !
          DO LAMBDA2 = 1, LBD1M1                                    !
            H(LAMBDA2,LAMBDA1) = ONEC                               !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
!   Calculation of the path operator TAU
!
      DO LF = LF1, LF2, ISTEP_LF                                    !
        ILF = LF * LF + LF + 1                                      !
!
        NU1MAX1 = MIN(LF,NUMX)                                      !
!
!  Case MF = 0
!
        DO LJ = 0, LMJ                                              !
          ILJ = LJ * LJ + LJ + 1                                    !
          NUJMAX = MIN(LJ,NUMAXJ)                                   !
          IF(JORDP == 1) THEN                                       !
            NU1MAX = MIN(NU1MAX1,LJ)                                !
          ELSE                                                      !
            NU1MAX = NU1MAX1                                        !
          END IF                                                    !
!
          TL_J = COEF * TL(LF,1,1,JE) * TL(LJ,1,JTYP,JE)            !
!
!    Case MJ = 0
!
          SUM_NU1_0 = ZEROC                                         !
!
          DO NU1 = 0, NU1MAX                                        !
            IF(JORDP > 1) THEN                                      !
              MU1MAX = MIN(LF-NU1,NO1-NU1-NU1)                      !
            ELSE                                                    !
              MU1MAX = MIN(LF-NU1,NO1-NU1-NU1,LJ)                   !
            END IF                                                  !
!
            DO MU1 = - MU1MAX, MU1MAX                               !
              LAMBDA1 = LBD(MU1,NU1)                                !
              AMU1    = ABS(MU1)                                    !
!
              RLMF_0 = HLM01(AMU1+NU1,LF) * RLM01(MU1,0,LF)         !
!
              SUM_NUJ_0 = ZEROC                                     !
!
              IF(JORDP > 1) THEN                                    !
                DO NUJ = 0, NUJMAX                                  !
                  MUJMAX = MIN(LJ,NOJ-NUJ-NUJ)                      !

                  SUM_MUJ_0 = ZEROC                                 !

                  DO MUJ = - MUJMAX, MUJMAX                         !

                    LAMBDAJ = LBD(MUJ,NUJ)                          !

                    SUM_MUJ_0 = SUM_MUJ_0 + H(LAMBDAJ,LAMBDA1) *  & !
                                            RLMIJ(MUJ,0,LJ)         !

                  END DO                                            !

                  SUM_NUJ_0 = SUM_NUJ_0 + SUM_MUJ_0 * HLMIJ(NUJ,LJ) !
!
                END DO                                              !
              ELSE                                                  !
                SUM_NUJ_0 = HLMIJ(NU1,LJ) * RLMIJ(MU1,0,LJ)         !
              END IF                                                !
!
              SUM_NU1_0 = SUM_NU1_0 + RLMF_0 * SUM_NUJ_0            !
!
            END DO                                                  !
!
          END DO                                                    !
!
          TAU(ILJ,ILF,JATL) = TAU(ILJ,ILF,JATL) +                 & !
                              TL_J * SUM_NU1_0                      !
!
          IF(NPATHP == 0) GO TO 35                                  !
!
          FM2         = FMAX(JORDP)                                 !
          XINT        = ABS(TL_J * SUM_NU1_0)                       !
          XMAX        = MAX(XINT,XMAX)                              !
          FMAX(JORDP) = MAX(FM2,XINT)                               !
          IF((FMAX(JORDP)-FM2) > XCOMP) THEN                        !
            NPMA(JORDP) = NPATH(JORDP)                              !
          END IF                                                    !
          IF((IREF == 1) .AND. (JORDP == NCUT)) THEN                !
            FREF        = FMAX(JORDP)                               !
          END IF                                                    !
 35       CONTINUE                                                  !
!
!    Case MJ > 0
!
          CJ = ONEC                                                 !
          DO MJ = 1, LJ                                             !
            INDJ  = ILJ + MJ                                        !
            INDJP = ILJ - MJ                                        !
            CJ    = CJ * EXP_J                                      !
!
            SUM_NU1_0 = ZEROC                                       !
            SUM_NU1_1 = ZEROC                                       !
!
            DO NU1 = 0, NU1MAX                                      !
              IF(JORDP > 1) THEN                                    !
                MU1MAX = MIN(LF-NU1,NO1-NU1-NU1)                    !
              ELSE                                                  !
                MU1MAX = MIN(LF-NU1,NO1-NU1-NU1,LJ)                 !
              END IF                                                !
!
              DO MU1 = - MU1MAX, MU1MAX                             !
                LAMBDA1 = LBD(MU1,NU1)                              !
                AMU1    = ABS(MU1)                                  !
!
                RLMF_0 = HLM01(AMU1+NU1,LF) * RLM01(MU1,0,LF)       !
!
                SUM_NUJ_0 = ZEROC                                   !
                SUM_NUJ_1 = ZEROC                                   !
!
                IF(JORDP > 1) THEN                                  !
                  DO NUJ = 0, NUJMAX                                !
                    MUJMAX = MIN(LJ,NOJ-NUJ-NUJ)                    !
!
                    SUM_MUJ_0 = ZEROC                               !
                    SUM_MUJ_1 = ZEROC                               !
!
                    DO MUJ = - MUJMAX, MUJMAX                       !
!
                      LAMBDAJ = LBD(MUJ,NUJ)                        !
!
                      SUM_MUJ_1 = SUM_MUJ_1 + H(LAMBDAJ,LAMBDA1) *& !
                                              RLMIJ(MUJ,-MJ,LJ)     !
                      SUM_MUJ_0 = SUM_MUJ_0 + H(LAMBDAJ,LAMBDA1) *& !
                                              RLMIJ(MUJ,MJ,LJ)      !
!
                    END DO                                          !
!
                    SUM_NUJ_0 = SUM_NUJ_0 + SUM_MUJ_0 * HLMIJ(NUJ,LJ)
                    SUM_NUJ_1 = SUM_NUJ_1 + SUM_MUJ_1 * HLMIJ(NUJ,LJ)
!
                  END DO                                            !
                ELSE                                                !
                  SUM_NUJ_1 = HLMIJ(NU1,LJ) * RLMIJ(MU1,-MJ,LJ)     !
                  SUM_NUJ_0 = HLMIJ(NU1,LJ) * RLMIJ(MU1,MJ,LJ)      !
                END IF                                              !
!
                SUM_NU1_0 = SUM_NU1_0 + RLMF_0 * SUM_NUJ_0          !
                SUM_NU1_1 = SUM_NU1_1 + RLMF_0 * SUM_NUJ_1          !
!
              END DO                                                !
!
            END DO                                                  !
!
            TAU(INDJP,ILF,JATL) = TAU(INDJP,ILF,JATL) +           & !
                                  CONJG(CJ) * TL_J * SUM_NU1_1      !
            TAU(INDJ,ILF,JATL)  = TAU(INDJ,ILF,JATL) +            & !
                                        CJ  * TL_J * SUM_NU1_0      !
!
            IF(NPATHP == 0) GO TO 45                                !
!
            FM2         = FMAX(JORDP)                               !
            XINT1       = ABS(CJ * TL_J * SUM_NU1_0)                !
            XINT2       = ABS(CONJG(CJ) * TL_J * SUM_NU1_1)         !
            XMAX        = MAX(XINT1,XINT2,XMAX)                     !
            FMAX(JORDP) = MAX(FM2,XINT1,XINT2)                      !
            IF((FMAX(JORDP)-FM2) > XCOMP)  THEN                     !
                NPMA(JORDP) = NPATH(JORDP)                          !
            END IF                                                  !
            IF((IREF == 1) .AND. (JORDP == NCUT)) THEN              !
                FREF    = FMAX(JORDP)                               !
            END IF                                                  !
 45         CONTINUE                                                !
          END DO                                                    !
        END DO                                                      !
!
!  Case MF > 0
!
        CF = ONEC                                                   !
        DO MF = 1, LF                                               !
          INDF  = ILF + MF                                          !
          INDFP = ILF - MF                                          !
          CF    = CF * EXP_F                                        !
!
          DO LJ = 0, LMJ                                            !
            ILJ = LJ * LJ + LJ + 1                                  !
            NUJMAX = MIN(LJ,NUMAXJ)                                 !
            IF(JORDP == 1) THEN                                     !
              NU1MAX = MIN(NU1MAX1,LJ)                              !
            ELSE                                                    !
              NU1MAX = NU1MAX1                                      !
            END IF                                                  !
!
            TL_J = COEF * TL(LF,1,1,JE) * TL(LJ,1,JTYP,JE)          !
!
!    Case MJ = 0
!
            SUM_NU1_0 = ZEROC                                       !
            SUM_NU1_1 = ZEROC                                       !
!
            DO NU1 = 0, NU1MAX                                      !
              IF(JORDP > 1) THEN                                    !
                MU1MAX = MIN(LF-NU1,NO1-NU1-NU1)                    !
              ELSE                                                  !
                MU1MAX = MIN(LF-NU1,NO1-NU1-NU1,LJ)                 !
              END IF                                                !
!
              DO MU1 = - MU1MAX, MU1MAX                             !
                LAMBDA1 = LBD(MU1,NU1)                              !
                AMU1    = ABS(MU1)                                  !
!
                RLMF_1 = HLM01(AMU1+NU1,LF) * RLM01(MU1,-MF,LF)     !
                RLMF_0 = HLM01(AMU1+NU1,LF) * RLM01(MU1,MF,LF)      !
!
                SUM_NUJ_0 = ZEROC                                   !
!
                IF(JORDP > 1) THEN                                  !
                  DO NUJ =0 , NUJMAX                                !
                    MUJMAX = MIN(LJ,NOJ-NUJ-NUJ)                    !
!
                    SUM_MUJ_0 = ZEROC                               !
!
                    DO MUJ = - MUJMAX, MUJMAX                       !
!
                      LAMBDAJ = LBD(MUJ,NUJ)                        !
!
                      SUM_MUJ_0 = SUM_MUJ_0 + H(LAMBDAJ,LAMBDA1) *& !
                                              RLMIJ(MUJ,0,LJ)       !
!
                    END DO                                          !
!
                    SUM_NUJ_0 = SUM_NUJ_0 + SUM_MUJ_0 * HLMIJ(NUJ,LJ)
!
                  END DO                                            !
                ELSE
                  SUM_NUJ_0 = HLMIJ(NU1,LJ) * RLMIJ(MU1,0,LJ)       !
                END IF
!
                SUM_NU1_0 = SUM_NU1_0 + RLMF_0 * SUM_NUJ_0          !
                SUM_NU1_1 = SUM_NU1_1 + RLMF_1 * SUM_NUJ_0          !
!
              END DO                                                !
!
            END DO                                                  !
!
            TAU(ILJ,INDF,JATL)  = TAU(ILJ,INDF,JATL)  +           & !
                                  CF * TL_J * SUM_NU1_0             !
            TAU(ILJ,INDFP,JATL) = TAU(ILJ,INDFP,JATL) +           & !
                                  CONJG(CF) * TL_J * SUM_NU1_1      !
!
            IF(NPATHP == 0) GO TO 25                                !
!
            FM2         = FMAX(JORDP)                               !
            XINT1       = ABS(CF * TL_J * SUM_NU1_0)                !
            XINT2       = ABS(CONJG(CF) * TL_J * SUM_NU1_1)         !
            XMAX        = MAX(XINT1,XINT2,XMAX)                     !
            FMAX(JORDP) = MAX(FM2,XINT1,XINT2)                      !
            IF((FMAX(JORDP)-FM2) > XCOMP) THEN                      !
                NPMA(JORDP) = NPATH(JORDP)                          !
            END IF                                                  !
            IF((IREF == 1) .AND. (JORDP == NCUT)) THEN              !
                FREF        = FMAX(JORDP)                           !
            END IF                                                  !
 25         CONTINUE                                                !
!
!    Case MJ > 0
!
            CJ = ONEC                                               !
            DO MJ = 1, LJ                                           !
              INDJ  = ILJ + MJ                                      !
              INDJP = ILJ - MJ                                      !
              CJ    = CJ * EXP_J                                    !
!
              SUM_NU1_0 = ZEROC                                     !
              SUM_NU1_1 = ZEROC                                     !
              SUM_NU1_2 = ZEROC                                     !
              SUM_NU1_3 = ZEROC                                     !
!
              DO NU1 = 0, NU1MAX                                    !
                IF(JORDP > 1) THEN                                  !
                  MU1MAX = MIN(LF-NU1,NO1-NU1-NU1)                  !
                ELSE                                                !
                  MU1MAX = MIN(LF-NU1,NO1-NU1-NU1,LJ)               !
                END IF                                              !
!
                DO MU1 = - MU1MAX, MU1MAX                           !
                  LAMBDA1 = LBD(MU1,NU1)                            !
                  AMU1    = ABS(MU1)                                !
!
                  RLMF_1 = HLM01(AMU1+NU1,LF) * RLM01(MU1,-MF,LF)   !
                  RLMF_0 = HLM01(AMU1+NU1,LF) * RLM01(MU1,MF,LF)    !
!
                  SUM_NUJ_0 = ZEROC                                 !
                  SUM_NUJ_1 = ZEROC                                 !
!
                  IF(JORDP > 1) THEN                                !
                    DO NUJ = 0, NUJMAX                              !
                      MUJMAX = MIN(LJ,NOJ-NUJ-NUJ)                  !
!
                      SUM_MUJ_0 = ZEROC                             !
                      SUM_MUJ_1 = ZEROC                             !
!
                      DO MUJ = - MUJMAX, MUJMAX                     !
!
                        LAMBDAJ = LBD(MUJ,NUJ)                      !
!
                        SUM_MUJ_1 = SUM_MUJ_1 + H(LAMBDAJ,LAMBDA1) *&
                                                RLMIJ(MUJ,-MJ,LJ)   !
                        SUM_MUJ_0 = SUM_MUJ_0 + H(LAMBDAJ,LAMBDA1) *&
                                                RLMIJ(MUJ,MJ,LJ)    !
!
                      END DO                                        !
!
                      SUM_NUJ_0 = SUM_NUJ_0 + SUM_MUJ_0 *         & !
                                              HLMIJ(NUJ,LJ)         !
                      SUM_NUJ_1 = SUM_NUJ_1 + SUM_MUJ_1 *         & !
                                              HLMIJ(NUJ,LJ)         !
!
                    END DO                                          !
                  ELSE                                              !
                    SUM_NUJ_1 = HLMIJ(NU1,LJ) * RLMIJ(MU1,-MJ,LJ)   !
                    SUM_NUJ_0 = HLMIJ(NU1,LJ) * RLMIJ(MU1,MJ,LJ)    !
                  END IF                                            !
!
                  SUM_NU1_0 = SUM_NU1_0 + RLMF_0 * SUM_NUJ_0        !
                  SUM_NU1_1 = SUM_NU1_1 + RLMF_0 * SUM_NUJ_1        !
                  SUM_NU1_2 = SUM_NU1_2 + RLMF_1 * SUM_NUJ_0        !
                  SUM_NU1_3 = SUM_NU1_3 + RLMF_1 * SUM_NUJ_1        !
!
                END DO                                              !
!
              END DO                                                !
!
              TAU(INDJP,INDF,JATL)  = TAU(INDJP,INDF,JATL)  +     & !
                                      CF * CONJG(CJ) * TL_J *     & !
                                      SUM_NU1_1                     !
              TAU(INDJP,INDFP,JATL) = TAU(INDJP,INDFP,JATL) +     & !
                                      CONJG(CF * CJ) * TL_J *     & !
                                      SUM_NU1_3                     !
              TAU(INDJ,INDF,JATL)   = TAU(INDJ,INDF,JATL)   +     & !
                                      CF * CJ * TL_J * SUM_NU1_0    !
              TAU(INDJ,INDFP,JATL)  = TAU(INDJ,INDFP,JATL)  +     & !
                                      CONJG(CF) * CJ * TL_J *     & !
                                      SUM_NU1_2                     !
!
              IF(NPATHP == 0) GO TO 15                              !
!
              FM2         = FMAX(JORDP)                             !
              XINT1       = ABS(CF * CJ * TL_J * SUM_NU1_0)         !
              XINT2       = ABS(CF * CONJG(CJ) * TL_J * SUM_NU1_1)  !
              XINT3       = ABS(CONJG(CF) * CJ * TL_J * SUM_NU1_2)  !
              XINT4       = ABS(CONJG(CF * CJ) * TL_J * SUM_NU1_3)  !
              XMAX        = MAX(XINT1,XINT2,XINT3,XINT4,XMAX)       ! <----------
              FMAX(JORDP) = MAX(FM2,XINT1,XINT2,XINT3,XINT4)        !            \
              IF((FMAX(JORDP)-FM2) > XCOMP) THEN                    !             \
                   NPMA(JORDP) = NPATH(JORDP)                       !              \
              END IF                                                !               \
              IF((IREF == 1) .AND. (JORDP == NCUT)) THEN            !                \
                  FREF    = FMAX(JORDP)                             !                 \
              END IF                                                !                  \
 15           CONTINUE                                              !                   \
            END DO                                                  !                    |
          END DO                                                    !                    |
        END DO                                                      !                    |
      END DO                                                        !                    |
!                                                                                        |
      IF(NPATHP == 0) GO TO 16                                      !                    |
!                                                                                        |
      FMIN(JORDP) = MIN(FM1,XMAX)                                   !                    |
!                                                                                        |
      IF(XMAX > FMN(NPATHP)) THEN                                   !                    |
!                                                                   ! LOCATE finds that XMAX
        CALL LOCATE(FMN,NPATHP,XMAX,JMX)                            ! is located between
!                                                                   ! FMN(JMX)and FMN(JMX+1)
        DO KF = NPATHP, JMX+2, -1                                   !
          FMN(KF)  = FMN(KF-1)                                      ! XMAX enters the list at
          JON(KF)  = JON(KF-1)                                      ! position JMX+1 and all
          PATH(KF) = PATH(KF-1)                                     ! values after are shifted up by 1
          DMN(KF)  = DMN(KF-1)                                      !
          DO KD = 1, 10                                             ! --> the NPATHP^th value disappears
            JPON(KF,KD) = JPON(KF-1,KD)                             ! and the (NPATH_P-1)^th one becomes
          END DO                                                    ! the final value
        END DO                                                      !
        FMN(JMX+1)  = XMAX                                          ! intensity of the new path
        JON(JMX+1)  = JORDP                                         ! order of the new path
        PATH(JMX+1) = NPATH(JORDP)                                  ! absolute number of the new path
        DMN(JMX+1)  = D                                             ! length of the new path
        DO KD = 1, JORDP                                            !
          JPON(JMX+1,KD) = JPOS(KD,3)                               ! list of the atoms of the new path
        END DO                                                      !
      END IF                                                        !
!
      IF((FMIN(JORDP)-FM1) < -XCOMP) NPMI(JORDP) = NPATH(JORDP)     !
      IF((IPRINT == 3) .AND. (IJ == 1)) THEN                        !
        WRITE(IUSCR,1) JORDP,NPATH(JORDP),XMAX,D,                 & !
                       (JPOS(KD,3),KD=1,JORDP)                      !
      END IF                                                        !
!
  16  RETURN                                                        !
!
!  Format:
!
   1  FORMAT(9X,I2,2X,E12.6,7X,E12.6,1X,F6.3,1X,10(I3,2X))
!
      END SUBROUTINE PATHOP
!
END MODULE CALC_TAU0J_SE
