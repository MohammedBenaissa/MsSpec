C
C======================================================================
C
C           ************************************************************
C           * ******************************************************** *
C           * *                                                      * *
C           * *         MULTIPLE-SCATTERING SPIN-INDEPENDENT         * *
C           * *        RESONANT ELASTIC X-RAY SCATTERING CODE        * *
C           * *              BASED ON MATRIX INVERSION               * *
C           * *                                                      * *
C           * ******************************************************** *
C           ************************************************************
C
C
C
C
C  Written by D. Sebilleau, Equipe de Physique des Surfaces et Interfaces,
C                           Institut de Physique de Rennes,
C                           UMR CNRS-Universite 6251,
C                           Universite de Rennes-1,
C                           35042 Rennes-Cedex,
C                           France
C
C  Contributions : M. Gavaza, H.-F. Zhao, K. Hatada
C
C-----------------------------------------------------------------------
C
C     As a general rule in this code, although there might be a few
C     exceptions (...), a variable whose name starts with a 'I' is a
C     switch, with a 'J' is a loop index and with a 'N' is a number.
C
C     The main subroutines are :
C
C                * PHDDIF    : computes the photoelectron diffraction
C                              formula
C
C                * LEDDIF    : computes the low-energy electron
C                              diffraction formula
C
C                * XASDIF    : computes the EXAFS or XANES formula
C                              depending on the energy
C
C                * AEDDIF    : computes the Auger electron diffraction
C                              formula
C
C                * FINDPATHS : generates the multiple scattering
C                              paths the electron will follow
C
C                * PATHOP    : calculates the contribution of a given
C                              path to the scattering path operator
C
C                * MATDIF    : computes the Rehr-Albers scattering
C                              matrices
C
C     A subroutine called NAME_A is the Auger equivalent of subroutine
C     NAME. The essentail difference between NAME and NAME_A is that
C     they do not contain the same arrays.
C
C     Always remember, when changing the input data file, to keep the
C     format. The rule here is that the last digit of any integer or
C     character data must correspond to the tab (+) while for real data,
C     the tab precedes the point.
C
C     Do not forget, before submitting a calculation, to check the
C     consistency of the input data with the corresponding maximal
C     values in the include file.
C
C-----------------------------------------------------------------------
C
C     Please report any bug or problem to me at :
C
C                      didier.sebilleau@univ-rennes1.fr
C
C
C
C  Last modified : 20 Jan 2014
C
C=======================================================================
C
      PROGRAM MAIN
C
C  This routine reads the various input files and calls the subroutine
C               performing the requested calculation
C
      INCLUDE 'spec.inc'
C
      COMMON /ADSORB/ IADS,NATA,NADS1,NADS2,NADS3,ADS(3,900),NCOUCH
      COMMON /APPROX/ NDIF,NO,ISPHER,IFWD,NTHOUT,RTHFWD(NATP_M),
     1                IBWD(NATP_M),RTHBWD(NATP_M),IPW,NCUT,PCTINT,IPP,
     2                ISPEED,IATTS,ILENGTH,RLENGTH
      COMMON /ATOMS/ VALZ(NATCLU_M),VALZ_MAX,NZAT(NATCLU_M),
     1               I_GR,I_INV,CHEM(NATCLU_M)
      COMMON /AUGER/ NLIN_A(0:20),L_BOUNDS(0:20,2),AUGER
      COMMON /BASES/ ATBAS(3*NATP_M),VECBAS(9)
      COMMON /CLUSLIM/ X_MAX(NATCLU_M),X_MIN(NATCLU_M),
     1                 Y_MAX(NATCLU_M),Y_MIN(NATCLU_M),
     2                 VAL(NATCLU_M),NPLAN
      COMMON /COOR/ NATCLU,N_PROT,NATYP(NATM),NCHTYP(NATP_M),
     1              NCORR(NAT_EQ_M,NATP_M),INEW_AT(NATCLU_M),
     2              SYM_AT(3,NATCLU_M)
      COMMON /DEBWAL/ IDCM,IDWSPH,TD,QD,TEMP,RSJ,UJ2(NATM)
      COMMON /INDAT/ INDATA(100)
      COMMON /INIT_A/ LI_C,LI_I,LI_A
      COMMON /INIT_L/ LI,INITL,NNL,LF1,LF2,ISTEP_LF
      COMMON /INIT_J/ JF1,JF2,I_SO,S_O
      COMMON /INIT_M/ I_SHELL,I_MULT,L_MUL,J_MUL,S_MUL,MULTIPLET
      COMMON /INFILES/ INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6,
     1                 INFILE7,INFILE8,INFILE9
      COMMON /INUNITS/ IUI1,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,IUI9
      COMMON /LIMAMA/ NIV,COUPUR
      COMMON /LPMOY/ ILPM,NZA,XMTA,RHOTA,XLPM0
      COMMON /MASSAT/ XMT(NATM)
      COMMON /MILLER/ IH,IK,II,IL,IVG0,IVN(3)
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /PARCAL/ NPHI0,NE0,NTHETA0,NFTHET0,NEPS0
      COMMON /PARCAL_A/ NPHI_A,NE_A,NTHETA_A,NFTHET_A
      COMMON /RELADS/ NRELA,PCRELA(3)
      COMMON /RELAX/ IREL,NREL,PCREL(10),OMEGA1,OMEGA2
      COMMON /RESEAU/ NCRIST,NCENTR,IBAS,NAT,A,BSURA,CSURA,UNIT
      COMMON /RXSGEN/ IE_RX,IPHI_RX,ITHETA_RX,NE,NPHI,NTHETA,
     1                E_INI_RX,E_FIN_RX,PHI_FIN_RX,THETA_FIN_RX
      COMMON /RXSINI/ THLUM_I,PHILUM_I,IPOL_I,NEPS_I,INTERACT_I
      COMMON /RXSFIN/ THLUM_O,PHILUM_O,IPOL_O,NEPS_O,INTERACT_O
      COMMON /SPIN/ ISPIN,IDICHR,NSPIN,NSPIN2,ISFLIP,IR_DIA,NSTEP
      COMMON /TESTS/ ITEST,IPRINT,ISORT1,NPATHP,ISOM
      COMMON /TRANS/ DLT(NE_M,NATM,0:18,5),TL(0:NT_M,4,NATM,NE_M),
     1               VK(NE_M),VK2(NE_M),IPOTC,ITL,LMAX(NATM,NE_M)
      COMMON /TL_AED/ DLT_A(NE_M,NATM,0:18,2),TL_A(0:NT_M,4,NATM,NE_M),
     1                VK_A(NE_M),VK2_A(NE_M),IPOTC_A,ITL_A,
     2                LMAX_A(NATM,NE_M)
      COMMON /TYPCAL/ IPHI,IE,ITHETA,IFTHET,IMOD,IPOL,I_CP,
     1                I_EXT,I_TEST
      COMMON /TYPCAL_A/ IPHI_A,IE_A,ITHETA_A,IFTHET_A,IMOD_A,I_CP_A,
     1                  I_EXT_A,I_TEST_A
      COMMON /TYPEM/ NEMET,IESURF,IEMET(NEMET_M)
      COMMON /TYPEXP/ SPECTRO,INTERACT,STEREO
      COMMON /VALIN/ PHI0,E0,THETA0,THLUM,PHILUM,ELUM,VINT,NONVOL(100)
      COMMON /XMRHO/ XMAT(0:99),RHOAT(0:99)
C
      DIMENSION VEC(3,3),VB1(3),VB2(3),VB3(3),VBS(3)
      DIMENSION ROT(3,3),EMET(3)
      DIMENSION VAL2(NATCLU_M)
      DIMENSION IRE(NATCLU_M,2)
      DIMENSION REL(NATCLU_M),RHOT(NATM)
      DIMENSION ATOME(3,NATCLU_M),COORD(3,NATCLU_M)
      DIMENSION NTYP(NATCLU_M),NATYP_OLD(NATM)
      DIMENSION LMAX_TMP(NATM,NE_M),DIST12(NATCLU_M,NATCLU_M)
      DIMENSION IBWD_TMP(NATP_M),RTHFWD_TMP(NATP_M),RTHBWD_TMP(NATP_M)
      DIMENSION UJ2_TMP(NATM),RHOT_TMP(NATM),XMT_TMP(NATM)
C
      COMPLEX VK,TL,DLT,TLSTAR,RHOK(NE_M,NATM,0:18,5,NSPIN2_M)
      COMPLEX VK_A,TL_A,DLT_A,TLSTAR_A
      COMPLEX RHOK_A(0:NT_M,NATM,0:40,2,NSPIN2_M),RAD_D,RAD_E
C
      INTEGER S_MUL,INV(2)
C
      CHARACTER RIEN
      CHARACTER*1 B
      CHARACTER*2 CHEM,R
      CHARACTER*3 UNIT,S_O,SPECTRO,MULTIPLET,STEREO
      CHARACTER*6 AUGER
      CHARACTER*7 INTERACT,INTERACT_I,INTERACT_O
      CHARACTER*24 INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6
      CHARACTER*24 INFILE7,INFILE8,INFILE9
      CHARACTER*24 INDATA
      CHARACTER*30 TUNIT,DUMMY
C
      DATA PI,BOHR,SMALL/3.141593,0.529177,0.001/
      DATA INV /1,0/
C
      LE_MAX=0
C
c      READ(*,776) NFICHLEC
c      READ(*,776) ICOM
c      DO JF=1,NFICHLEC
c        READ(*,777) INDATA(JF)
c      ENDDO
      nfichlec=1
      icom=5
      indata(1)='data/spec5.dat'
C
C..........  Loop on the data files  ..........
C
      DO JFICH=1,NFICHLEC
        OPEN(UNIT=ICOM, FILE=INDATA(JFICH), STATUS='OLD')
        CALL READ_DATA(ICOM,NFICHLEC,JFICH,ITRTL,*2,*1,*55,*74,*99,*504,
     1               *520,*540,*550,*570,*580,*590,*630)
C
C..........  Atomic case index  ..........
C
        I_AT=0
        IF((SPECTRO.EQ.'PHD').AND.(I_TEST.EQ.2)) I_AT=1
        IF((SPECTRO.EQ.'LED').AND.(I_TEST.EQ.2)) I_AT=1
        IF((SPECTRO.EQ.'AED').AND.(I_TEST_A.EQ.2)) I_AT=1
        IF((SPECTRO.EQ.'XAS').AND.(I_TEST.EQ.2)) I_AT=1
        IF((SPECTRO.EQ.'RES').AND.(I_TEST.EQ.2)) I_AT=1
        IF(SPECTRO.EQ.'APC') THEN
          IF((I_TEST.EQ.2).AND.(I_TEST_A.EQ.2)) I_AT=1
        ENDIF
C
        IF(IBAS.EQ.1) THEN
          IF(ITEST.EQ.0) THEN
            NEQ=(2*NIV+1)**3
          ELSE
            NEQ=(2*NIV+3)**3
          ENDIF
          IF(NEQ*NATP_M.GT.NATCLU_M) GOTO 518
        ENDIF
C
        IF(SPECTRO.EQ.'APC') THEN
          N_EL=2
        ELSE
          N_EL=1
        ENDIF
        IF((INTERACT.EQ.'COULOMB').OR.(INTERACT.EQ.'PHOCOUL')) THEN
          IF(I_MULT.EQ.0) THEN
            LE_MIN=ABS(LI_C-ABS(LI_I-LI_A))
            LE_MAX=LI_C+LI_A+LI_I
          ELSE
            LE_MIN=ABS(LI_C-L_MUL)
            LE_MAX=LI_C+L_MUL
          ENDIF
        ENDIF
C
C..........  Test of the dimensions against the input values  ..........
C
        IF(NO.GT.NO_ST_M) GOTO 600
        IF(LE_MAX.GT.LI_M) GOTO 620
C
        OPEN(UNIT=IUI2, FILE=INFILE2, STATUS='OLD')
        OPEN(UNIT=IUI3, FILE=INFILE3, STATUS='OLD')
        IF(INTERACT.EQ.'PHOCOUL') THEN
          OPEN(UNIT=IUI7, FILE=INFILE7, STATUS='OLD')
          OPEN(UNIT=IUI8, FILE=INFILE8, STATUS='OLD')
        ENDIF
C
C..........  Reading of the TL and radial matrix elements files  ..........
C..........      (dipolar excitation or no excitation case)       ..........
C
        IF(INTERACT.NE.'COULOMB') THEN
          IF(SPECTRO.EQ.'APC') WRITE(IUO1,418)
          READ(IUI2,3) NAT1,NE1,ITL,IPOTC,LMAX_MODE
          IF(ISPIN.EQ.0) THEN
            IF(NAT1.EQ.1) THEN
              WRITE(IUO1,561)
            ELSE
              WRITE(IUO1,560) NAT1
            ENDIF
          ENDIF
          IF((ITL.EQ.1).AND.(ISPIN.EQ.1)) THEN
            READ(IUI2,530) E_MIN,E_MAX,DE
          ENDIF
          IF((ISPIN.EQ.0).AND.(ITL.EQ.0)) THEN
            NLG=INT(NAT1-0.0001)/4 +1
            DO NN=1,NLG
              NRL=4*NN
              JD=4*(NN-1)+1
              IF(NN.EQ.NLG) NRL=NAT1
              READ(IUI2,555) (LMAX(JAT,1),JAT=JD,NRL)
              WRITE(IUO1,556) (LMAX(JAT,1),JAT=JD,NRL)
            ENDDO
C
C   Temporary storage of LMAX. Waiting for a version of PHAGEN
C           with LMAX dependent on the energy
C
            DO JE=1,NE
              DO JAT=1,NAT1
                LMAX(JAT,JE)=LMAX(JAT,1)
              ENDDO
            ENDDO
C
            NL1=1
            DO JAT=1,NAT1
              NL1=MAX0(NL1,LMAX(JAT,1)+1)
            ENDDO
            IF(NL1.GT.NL_M) GOTO 184
          ENDIF
          IF(ITL.EQ.0) READ(IUI3,101) NATR,NER
          IF(ISPIN.EQ.1) THEN
            READ(IUI3,106) L_IN,NATR,NER
            IF(LI.NE.L_IN) GOTO 606
          ENDIF
          NAT2=NAT+NATA
C
          IF((NAT1.NE.NAT2).OR.(NE1.NE.NE)) GOTO 180
          IF((ITL.EQ.0).AND.((NATR.NE.NAT2).OR.(NER.NE.NE))) GOTO 182
C
C..........  DL generated by MUFPOT and RHOK given  ..........
C..........     by S. M. Goldberg, C. S. Fadley     ..........
C..........     and S. Kono, J. Electron Spectr.    ..........
C..........      Relat. Phenom. 21, 285 (1981)      ..........
C
C
C       for Quadrupole approximation (in REXS case), RHOK and DLT should
C       be expanded from dimension 2 to 5
C
          IF(ITL.EQ.0) THEN
            DO JAT=1,NAT2
              IF((INITL.NE.0).AND.(IFTHET.NE.1)) THEN
                READ(IUI3,102) RIEN
                READ(IUI3,102) RIEN
                READ(IUI3,102) RIEN
              ENDIF
              DO JE=1,NE
                IF((IFTHET.EQ.1).OR.(INITL.EQ.0)) GOTO 121
                READ(IUI3,103) ENERGIE
                READ(IUI3,102) RIEN
                READ(IUI3,102) RIEN
                READ(IUI3,102) RIEN
 121            CONTINUE
                DO L=0,LMAX(JAT,JE)
                  READ(IUI2,7) VK(JE),TL(L,1,JAT,JE)
                  TL(L,1,JAT,JE)=CSIN(TL(L,1,JAT,JE))*CEXP((0.,1.)*
     1                         TL(L,1,JAT,JE))
                ENDDO
                IF((IFTHET.EQ.1).OR.(INITL.EQ.0)) GOTO 5
                DO LL=1,18
                  READ(IUI3,104) RH1,RH2,DEF1,DEF2
                  RHOK(JE,JAT,LL,1,1)=CMPLX(RH1)
                  RHOK(JE,JAT,LL,2,1)=CMPLX(RH2)
                  DLT(JE,JAT,LL,1)=CMPLX(DEF1)
                  DLT(JE,JAT,LL,2)=CMPLX(DEF2)
                ENDDO
  5             CONTINUE
              ENDDO
            ENDDO
          ELSE
C
C..........  TL and RHOK calculated by PHAGEN  ..........
C
            DO JE=1,NE
              NLG=INT(NAT2-0.0001)/4 +1
              IF(NE.GT.1) WRITE(IUO1,563) JE
              DO NN=1,NLG
                NRL=4*NN
                JD=4*(NN-1)+1
                IF(NN.EQ.NLG) NRL=NAT2
                READ(IUI2,555) (LMAX(JAT,JE),JAT=JD,NRL)
                WRITE(IUO1,556) (LMAX(JAT,JE),JAT=JD,NRL)
              ENDDO
              NL1=1
              DO JAT=1,NAT2
                NL1=MAX0(NL1,LMAX(JAT,1)+1)
              ENDDO
              IF(NL1.GT.NL_M) GOTO 184
              DO JAT=1,NAT2
                READ(IUI2,*) DUMMY
                DO L=0,LMAX(JAT,JE)
                  IF(LMAX_MODE.EQ.0) THEN
                    READ(IUI2,9) VK(JE),TLSTAR
                  ELSE
                    READ(IUI2,9) VK(JE),TLSTAR
                  ENDIF
                  TL(L,1,JAT,JE)=CONJG(TLSTAR)
                  VK(JE)=CONJG(VK(JE))
                ENDDO
              ENDDO
C
              IF((IFTHET.EQ.1).OR.(INITL.EQ.0)) GOTO 333
              IF(JE.EQ.1) THEN
                DO JDUM=1,4
                  READ(IUI3,102) RIEN
                ENDDO
              ENDIF
              DO JEMET=1,NEMET
                JM=IEMET(JEMET)
                READ(IUI3,105) (RHOK(JE,JM,NNL,IDQ,1),IDQ=1,5)
              ENDDO
  333         VK(JE)=VK(JE)*A
              VK2(JE)=CABS(VK(JE)*VK(JE))
            ENDDO
          ENDIF
C
          CLOSE(IUI2)
          CLOSE(IUI3)
C
C.......... Suppression of possible zeros in the TL array  ..........
C..........  (in case of the use of matrix inversion and   ..........
C..........             for energy variations)             ..........
C
          IF((ISPIN.EQ.0).AND.(ITL.EQ.1).AND.(LMAX_MODE.NE.0)) THEN
            CALL SUP_ZEROS(TL,LMAX,NE,NAT2,IUO1,ITRTL)
          ENDIF

        ENDIF
C
C..........  Reading of the TL and radial matrix elements files  ..........
C..........             (Coulomb excitation case)                ..........
C
        IF((INTERACT.EQ.'COULOMB').OR.(INTERACT.EQ.'PHOCOUL')) THEN
          IERR=0
          IF(INTERACT.EQ.'COULOMB') THEN
            IRD1=IUI2
            IRD2=IUI3
          ELSEIF(INTERACT.EQ.'PHOCOUL') THEN
            IRD1=IUI7
            IRD2=IUI8
          ENDIF
          IF(SPECTRO.EQ.'APC') WRITE(IUO1,419)
          READ(IRD1,3) NAT1_A,NE1_A,ITL_A,IPOTC_A,LMAX_MODE_A
          IF(ISPIN.EQ.0) THEN
            IF(NAT1_A.EQ.1) THEN
              WRITE(IUO1,561)
            ELSE
              WRITE(IUO1,560) NAT1_A
            ENDIF
          ENDIF
          IF((ITL_A.EQ.1).AND.(ISPIN.EQ.1)) THEN
            READ(IRD1,530) E_MIN_A,E_MAX_A,DE_A
          ENDIF
          IF(ITL_A.EQ.1) THEN
            READ(IRD2,107) LI_C2,LI_I2,LI_A2
            READ(IRD2,117) LE_MIN1,N_CHANNEL
            LE_MAX1=LE_MIN1+N_CHANNEL-1
            IF(I_TEST_A.NE.1) THEN
              IF((LE_MIN.NE.LE_MIN1).OR.(LE_MAX.NE.LE_MAX1)) GOTO 610
            ELSE
              LI_C2=0
              LI_I2=1
              LI_A2=0
              LE_MIN1=1
              N_CHANNEL=1
            ENDIF
          ENDIF
          IF((ISPIN.EQ.0).AND.(ITL_A.EQ.0)) THEN
            NLG=INT(NAT1_A-0.0001)/4 +1
            DO NN=1,NLG
              NRL=4*NN
              JD=4*(NN-1)+1
              IF(NN.EQ.NLG) NRL=NAT1_A
              READ(IRD1,555) (LMAX_A(JAT,1),JAT=JD,NRL)
              WRITE(IUO1,556) (LMAX_A(JAT,1),JAT=JD,NRL)
            ENDDO
C
C   Temporary storage of LMAX_A. Waiting for a version of PHAGEN
C           with LMAX_A dependent on the energy
C
            DO JE=1,NE1_A
              DO JAT=1,NAT1_A
                LMAX_A(JAT,JE)=LMAX_A(JAT,1)
              ENDDO
            ENDDO
C
            NL1_A=1
            DO JAT=1,NAT1_A
              NL1_A=MAX0(NL1_A,LMAX_A(JAT,1)+1)
            ENDDO
            IF(NL1_A.GT.NL_M) GOTO 184
          ENDIF
          IF(ITL_A.EQ.0) READ(IRD2,101) NATR_A,NER_A
          IF(ISPIN.EQ.1) THEN
            READ(IRD2,106) L_IN_A,NATR_A,NER_A
            IF(LI_C.NE.L_IN_A) GOTO 606
          ENDIF
          NAT2_A=NAT+NATA
          NAT2=NAT2_A
          IF((NAT1_A.NE.NAT2_A).OR.(NE1_A.NE.NE_A)) GOTO 180
          IF((ITL_A.EQ.0).AND.((NATR_A.NE.NAT2_A).OR.(NER_A.NE.NE)))
     1    GOTO 182
C
C..........  DL generated by MUFPOT and RHOK given  ..........
C..........     by S. M. Goldberg, C. S. Fadley     ..........
C..........     and S. Kono, J. Electron Spectr.    ..........
C..........      Relat. Phenom. 21, 285 (1981)      ..........
C
          IF(ITL_A.EQ.0) THEN
            CONTINUE
          ELSE
C
C..........  TL_A and RHOK_A calculated by PHAGEN  ..........
C
            DO JE=1,NE_A
              NLG=INT(NAT2_A-0.0001)/4 +1
              IF(NE_A.GT.1) WRITE(IUO1,563) JE
              DO NN=1,NLG
                NRL=4*NN
                JD=4*(NN-1)+1
                IF(NN.EQ.NLG) NRL=NAT2_A
                READ(IRD1,555) (LMAX_A(JAT,JE),JAT=JD,NRL)
                WRITE(IUO1,556) (LMAX_A(JAT,JE),JAT=JD,NRL)
              ENDDO
              DO JAT=1,NAT2_A
                READ(IRD1,*) DUMMY
                DO L=0,LMAX_A(JAT,JE)
                  IF(LMAX_MODE_A.EQ.0) THEN
                    READ(IRD1,9) VK_A(JE),TLSTAR
                  ELSE
                    READ(IRD1,7) VK_A(JE),TLSTAR
                  ENDIF
                  TL_A(L,1,JAT,JE)=CONJG(TLSTAR)
                  VK_A(JE)=CONJG(VK_A(JE))
                ENDDO
              ENDDO
C
              IF(IFTHET_A.EQ.1) GOTO 331
              DO LE=LE_MIN,LE_MAX
                DO JEMET=1,NEMET
                  JM=IEMET(JEMET)
                  READ(IRD2,109) L_E,LB_MIN,LB_MAX
                  IF(I_TEST_A.EQ.1) THEN
                    L_E=1
                    LB_MIN=0
                    LB_MAX=1
                  ENDIF
                  IF(LE.NE.L_E) IERR=1
                  L_BOUNDS(L_E,1)=LB_MIN
                  L_BOUNDS(L_E,2)=LB_MAX
                  DO LB=LB_MIN,LB_MAX
                    READ(IRD2,108) L_A,RAD_D,RAD_E
                    RHOK_A(LE,JM,L_A,1,1)=RAD_D
                    RHOK_A(LE,JM,L_A,2,1)=RAD_E
                    IF(I_TEST_A.EQ.1) THEN
                      IF(LB.EQ.LB_MIN) THEN
                        RHOK_A(LE,JM,L_A,1,1)=(0.0,0.0)
                        RHOK_A(LE,JM,L_A,2,1)=(1.0,0.0)
                      ELSEIF(LB.EQ.LB_MAX) THEN
                        RHOK_A(LE,JM,L_A,1,1)=(1.0,0.0)
                        RHOK_A(LE,JM,L_A,2,1)=(0.0,0.0)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
  331         VK_A(JE)=VK_A(JE)*A
              VK2_A(JE)=CABS(VK_A(JE)*VK_A(JE))
            ENDDO
          ENDIF
C
          CLOSE(IRD1)
          CLOSE(IRD2)
C
C.......... Suppression of possible zeros in the TL array  ..........
C..........  (in case of the use of matrix inversion and   ..........
C..........             for energy variations)             ..........
C
          IF((ISPIN.EQ.0).AND.(ITL_A.EQ.1).AND.(LMAX_MODE_A.NE.0)) THEN
            CALL SUP_ZEROS(TL_A,LMAX_A,NE_A,NAT2_A,IUO1,ITRTL)
          ENDIF
          IF(SPECTRO.EQ.'APC') WRITE(IUO1,420)
C
        ENDIF
C
C..........  Check of the consistency of the two TL and radial ..........
C..........           matrix elements for APECS                ..........
C
        IF(SPECTRO.EQ.'APC') THEN
C
          I_TL_FILE=0
          I_RD_FILE=0
C
          IF(NAT1.NE.NAT1_A) I_TL_FILE=1
          IF(NE1.NE.NE1_A) I_TL_FILE=1
          IF(ITL.NE.ITL_A) I_TL_FILE=1
          IF(IPOTC.NE.IPOTC_A) I_TL_FILE=1
C
          IF(LI_C.NE.LI_C2) I_RD_FILE=1
          IF(LI_I.NE.LI_I2) I_RD_FILE=1
          IF(LI_A.NE.LI_A2) I_RD_FILE=1
C
          IF(I_TL_FILE.EQ.1) GOTO 608
          IF(I_RD_FILE.EQ.1) GOTO 610
          IF(IERR.EQ.1) GOTO 610
C
        ENDIF
C
C..........  Calculation of the scattering factor (only)  ..........
C
        IF((IFTHET.EQ.0).AND.(IFTHET_A.EQ.0)) GO TO 8
        IF(IFTHET.EQ.1) THEN
          CALL PLOTFD(A,LMAX,ITL,NL1,NAT2,NE)
        ELSEIF(IFTHET_A.EQ.1) THEN
c        CALL PLOTFD_A(A,LMAX_A,ITL_A,NL1_A,NAT2_A,NE_A)
        ENDIF
        WRITE(IUO1,57)
        STOP
C
   8    IF(IBAS.EQ.0) THEN
C
C...............  Reading of an external cluster  ...............
C
C
C  Cluster originating from CLUSTER_NEW.F : IPHA=0
C  Cluster originating from PHAGEN_NEW.F  : IPHA=1 (atomic units), IPHA=2 (angstroems)
C  Other cluster                          : the first line must be text; then
C                                            free format : Atomic number,X,Y,Z,number
C                                            of the corresponding prototypical atom ;
C                                            All atoms corresponding to the same
C                                            prototypical atom must follow each other.
C                                            Moreover, the blocks of equivalent atoms
C                                            must be ordered by increasing number of
C                                            prototypical atom.
C
          VALZ_MIN=1000.0
          VALZ_MAX=-1000.0
C
          OPEN(UNIT=IUI4, FILE=INFILE4, STATUS='OLD')
          READ(IUI4,778,ERR=892) IPHA
          GOTO 893
  892     IPHA=3
          IF(UNIT.EQ.'ANG') THEN
            CUNIT=1./A
            TUNIT='ANGSTROEMS'
          ELSEIF(UNIT.EQ.'LPU') THEN
            CUNIT=1.
            TUNIT='UNITS OF THE LATTICE PARAMETER'
          ELSEIF(UNIT.EQ.'ATU') THEN
            CUNIT=BOHR/A
            TUNIT='ATOMIC UNITS'
          ELSE
            GOTO 890
          ENDIF
  893     NATCLU=0
          DO JAT=1,NAT2
            NATYP(JAT)=0
          ENDDO
          IF(IPHA.EQ.0) THEN
            CUNIT=1.
            TUNIT='UNITS OF THE LATTICE PARAMETER'
          ELSEIF(IPHA.EQ.1) THEN
            CUNIT=BOHR/A
            TUNIT='ATOMIC UNITS'
            IEMET(1)=1
          ELSEIF(IPHA.EQ.2) THEN
            CUNIT=1./A
            TUNIT='ANGSTROEMS'
            IEMET(1)=1
          ENDIF
          IF(IPRINT.EQ.2) THEN
            IF(I_AT.NE.1) THEN
              WRITE(IUO1,558) IUI4,TUNIT
              IF(IPHA.EQ.3) WRITE(IUO1,549)
            ENDIF
          ENDIF
          JATM=0
          DO JLINE=1,10000
            IF(IPHA.EQ.0) THEN
              READ(IUI4,125,END=780) R,NN,X,Y,Z,JAT
            ELSEIF(IPHA.EQ.1) THEN
              READ(IUI4,779,END=780) R,NN,X,Y,Z,JAT
            ELSEIF(IPHA.EQ.2) THEN
              READ(IUI4,779,END=780) R,NN,X,Y,Z,JAT
            ELSEIF(IPHA.EQ.3) THEN
              READ(IUI4,*,END=780) NN,X,Y,Z,JAT
            ENDIF
            JATM=MAX0(JAT,JATM)
            NATCLU=NATCLU+1
            IF(IPHA.NE.3) THEN
              CHEM(JAT)=R
            ELSE
              CHEM(JAT)='XX'
            ENDIF
            NZAT(JAT)=NN
            NATYP(JAT)=NATYP(JAT)+1
            COORD(1,NATCLU)=X*CUNIT
            COORD(2,NATCLU)=Y*CUNIT
            COORD(3,NATCLU)=Z*CUNIT
            VALZ(NATCLU)=Z*CUNIT
            IF((IPRINT.GE.2).AND.(I_AT.EQ.0)) THEN
              WRITE(IUO1,557) NATCLU,COORD(1,NATCLU),COORD(2,NATCLU),
     1                      COORD(3,NATCLU),JAT,NATYP(JAT),CHEM(JAT)
            ENDIF
          ENDDO
 780      NBZ=NATCLU
          IF(JATM.NE.NAT) GOTO 514
          CLOSE(IUI4)
C
          IF(NATCLU.GT.NATCLU_M) GOTO 510
          DO JA1=1,NATCLU
            DO JA2=1,NATCLU
              DIST12(JA1,JA2)=SQRT((COORD(1,JA1)-COORD(1,JA2))**2
     1                          +(COORD(2,JA1)-COORD(2,JA2))**2
     2                          +(COORD(3,JA1)-COORD(3,JA2))**2)
              IF((JA2.GT.JA1).AND.(DIST12(JA1,JA2).LT.0.001)) GOTO 895
            ENDDO
          ENDDO
C
          D_UP=VALZ_MAX-VALZ(1)
          D_DO=VALZ(1)-VALZ_MIN
          IF((D_DO.LE.D_UP).AND.(I_GR.EQ.2)) THEN
            I_INV=1
          ELSE
            I_INV=0
          ENDIF
        ELSE
C
C...............  Construction of an internal cluster  ...............
C
          CALL BASE
          CALL ROTBAS(ROT)
          IF(IVG0.EQ.2) THEN
            NMAX=NIV+1
          ELSE
            NMAX=(2*NIV+1)**3
          ENDIF
          IF((IPRINT.EQ.2).AND.(IVG0.LE.1)) THEN
            WRITE(IUO1,37)
            WRITE(IUO1,38) NIV
            DO NUM=1,NMAX
              CALL NUMAT(NUM,NIV,IA,IB,IC)
              WRITE(IUO1,17) NUM,IA,IB,IC
            ENDDO
            WRITE(IUO1,39)
          ENDIF
          CALL AMAS(NIV,ATOME,COORD,VALZ,IESURF,COUPUR,ROT,
     1            IRE,NATYP,NBZ,NAT2,NCOUCH,NMAX)
          IF((IREL.GE.1).OR.(NRELA.GT.0)) THEN
            CALL RELA(NBZ,NPLAN,NAT2,VALZ,VAL2,VAL,COORD,NATYP,REL,
     1              NCOUCH)
            IF(IREL.EQ.1) THEN
              DO JP=1,NPLAN
                VAL(JP)=VAL2(JP)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
C
C  Storage of the extremal values of x and y for each plane. They define
C    the exterior of the cluster when a new cluster has to be build to
C    support a point-group
C
        IF(I_GR.GE.1) THEN
          IF((IREL.EQ.0).OR.(IBAS.EQ.0)) THEN
            CALL ORDRE(NBZ,VALZ,NPLAN,VAL)
            WRITE(IUO1,50) NPLAN
            DO K=1,NPLAN
              WRITE(IUO1,29) K,VAL(K)
              X_MAX(K)=0.
              X_MIN(K)=0.
              Y_MAX(K)=0.
              Y_MIN(K)=0.
            ENDDO
          ENDIF
          DO JAT=1,NATCLU
            X=COORD(1,JAT)
            Y=COORD(2,JAT)
            Z=COORD(3,JAT)
            DO JPLAN=1,NPLAN
              IF(ABS(Z-VAL(JPLAN)).LT.SMALL) THEN
                X_MAX(JPLAN)=MAX(X,X_MAX(JPLAN))
                X_MIN(JPLAN)=MIN(X,X_MIN(JPLAN))
                Y_MAX(JPLAN)=MAX(Y,Y_MAX(JPLAN))
                Y_MIN(JPLAN)=MIN(Y,Y_MIN(JPLAN))
              ENDIF
            ENDDO
          ENDDO
        ENDIF
C
C  Instead of the symmetrization of the cluster (this version only)
C
        N_PROT=NAT
        NAT_ST=0
        DO JTYP=1,JATM
          NB_AT=NATYP(JTYP)
          IF(NB_AT.GT.NAT_EQ_M) GOTO 614
          DO JA=1,NB_AT
            NAT_ST=NAT_ST+1
            NCORR(JA,JTYP)=NAT_ST
          ENDDO
        ENDDO
        DO JC=1,3
          DO JA=1,NATCLU
            SYM_AT(JC,JA)=COORD(JC,JA)
          ENDDO
        ENDDO
C
C  Checking surface-like atoms for mean square displacements
C                    calculations
C
        CALL CHECK_VIB(NAT2)
C
C..........  Set up of the variables used for an internal  ..........
C..........   calculation of the mean free path and/or of  ..........
C..........         the mean square displacements          ..........
C
        IF((IDCM.EQ.1).OR.(ILPM.EQ.1)) THEN
          DO JTYP=1,NAT2
            XMT(JTYP)=XMAT(NZAT(JTYP))
            RHOT(JTYP)=RHOAT(NZAT(JTYP))
          ENDDO
          XMTA=XMT(1)
          RHOTA=RHOT(1)
          NZA=NZAT(1)
        ENDIF
        IF(IDCM.GT.0) THEN
          CALL CHNOT(3,VECBAS,VEC)
          DO J=1,3
            VB1(J)=VEC(J,1)
            VB2(J)=VEC(J,2)
            VB3(J)=VEC(J,3)
          ENDDO
          CPR=1.
          CALL PRVECT(VB2,VB3,VBS,CPR)
          VM=PRSCAL(VB1,VBS)
          QD=(6.*PI*PI*NAT/VM)**(1./3.)
        ENDIF
C
C..........  Writing of the contents of the cluster,  ..........
C..........  of the position of the different planes  ..........
C..........   and of their respective absorbers in    ..........
C..........          the control file IUO1            ..........
C
        IF(I_AT.EQ.1) GOTO 153
        IF((IPRINT.EQ.2).AND.(IBAS.GT.0)) THEN
          WRITE(IUO1,40)
          NCA=0
          DO J=1,NAT
            DO I=1,NMAX
              NCA=NCA+1
              WRITE(IUO1,20) J,I
              WRITE(IUO1,21) (ATOME(L,NCA),L=1,3)
              K=IRE(NCA,1)
              IF(K.EQ.0) THEN
                WRITE(IUO1,22)
              ELSE
                WRITE(IUO1,23) (COORD(L,K),L=1,3),IRE(NCA,2)
              ENDIF
            ENDDO
          ENDDO
          WRITE(IUO1,41)
        ENDIF
        IF(IBAS.EQ.1) THEN
          WRITE(IUO1,24)
          NATCLU=0
          DO I=1,NAT
            NN=NATYP(I)
            NATCLU=NATCLU+NATYP(I)
            WRITE(IUO1,26) NN,I
          ENDDO
          IF(IADS.EQ.1) NATCLU=NATCLU+NADS1+NADS2+NADS3
          WRITE(IUO1,782) NATCLU
          IF(NATCLU.GT.NATCLU_M) GOTO 516
          IF(IPRINT.EQ.3) WRITE(IUO1,559)
          IF(IPRINT.EQ.3) THEN
            NBTA=0
            DO JT=1,NAT2
              NBJT=NATYP(JT)
              DO JN=1,NBJT
                NBTA=NBTA+1
                WRITE(IUO1,557) NBTA,COORD(1,NBTA),COORD(2,NBTA),
     1                        COORD(3,NBTA),JT,JN,CHEM(JT)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
 153    IF((ITEST.EQ.1).AND.(IBAS.GT.0)) THEN
          CALL TEST(NIV,ROT,NATYP,NBZ,NAT2,IESURF,COUPUR,*56)
        ENDIF
        IF((IREL.EQ.0).OR.(IBAS.EQ.0)) THEN
          CALL ORDRE(NBZ,VALZ,NPLAN,VAL)
          IF(I_AT.EQ.0) WRITE(IUO1,50) NPLAN
          DO K=1,NPLAN
            IF(I_AT.EQ.0) WRITE(IUO1,29) K,VAL(K)
          ENDDO
        ENDIF
C
        IF(I_AT.EQ.0) WRITE(IUO1,30)
        IF((IPRINT.GT.0).AND.(I_AT.EQ.0)) THEN
          WRITE(IUO1,31) (IEMET(J),J=1,NEMET)
        ENDIF
        ZEM=1.E+20
        DO L=1,NPLAN
          Z=VAL(L)
          DO JEMED=1,NEMET
            CALL EMETT(JEMED,IEMET,Z,COORD,NATYP,EMET,NTEM,JNEM,*93)
            IF(I_AT.EQ.0) WRITE(IUO1,34) L,NTEM,EMET(1),EMET(2),EMET(3)
            IF((IPHA.EQ.1).OR.(IPHA.EQ.2)) ZEM=EMET(3)
            GO TO 33
  93        IF(I_AT.EQ.0) WRITE(IUO1,94) L,NTEM
  33        CONTINUE
          ENDDO
        ENDDO
C
C..........  Loop on the electrons involved in the  ..........
C..........  spectroscopy :  N_EL = 1 for PHD, XAS  ..........
C..........  LEED AED or RES and N_EL = 2 for APC   ..........
C
        DO J_EL=1,N_EL
C
C..........  Writing the information on the spectroscopies  ..........
C..........             in the control file IUO1            ..........
C
          IF(SPECTRO.EQ.'XAS') GOTO 566
          IF((SPECTRO.EQ.'APC').AND.(J_EL.EQ.1)) THEN
            IF(IPHI.EQ.1) THEN
              IF(STEREO.EQ.' NO') THEN
                WRITE(IUO1,236)
              ELSE
                WRITE(IUO1,248)
              ENDIF
            ENDIF
            IF(ITHETA.EQ.1) WRITE(IUO1,245)
            IF(I_TEST.EQ.1) WRITE(IUO1,234)
          ENDIF
C
C----------  Photoelectron diffraction case (PHD)  ----------
C
          IF((SPECTRO.EQ.'PHD').OR.(SPECTRO.EQ.'APC')) THEN
            IF(SPECTRO.EQ.'PHD') THEN
              IF(IPHI.EQ.1) THEN
                IF(STEREO.EQ.' NO') THEN
                  WRITE(IUO1,35)
                ELSE
                  WRITE(IUO1,246)
                ENDIF
              ENDIF
              IF(ITHETA.EQ.1) WRITE(IUO1,44)
              IF(IE.EQ.1) WRITE(IUO1,58)
              IF(INITL.EQ.0) WRITE(IUO1,118)
              IF(I_TEST.EQ.1) WRITE(IUO1,234)
            ENDIF
            IF((SPECTRO.EQ.'APC').AND.(J_EL.EQ.1)) THEN
              WRITE(IUO1,418)
              WRITE(IUO1,18)
            ENDIF
            IF(J_EL.EQ.2) GOTO 222
            IF(IPRINT.GT.0) THEN
              WRITE(IUO1,92)
              WRITE(IUO1,91)
              IF(ISPIN.EQ.0) THEN
                WRITE(IUO1,335)
              ELSE
                WRITE(IUO1,336)
              ENDIF
              WRITE(IUO1,91)
              IF(IPOTC.EQ.0) THEN
                WRITE(IUO1,339)
              ELSE
                WRITE(IUO1,334)
              ENDIF
              WRITE(IUO1,91)
              IF(INITL.NE.0) THEN
                WRITE(IUO1,337)
                WRITE(IUO1,91)
                IF(IPOL.EQ.0) THEN
                  WRITE(IUO1,88)
                ELSEIF(ABS(IPOL).EQ.1) THEN
                  WRITE(IUO1,87)
                ELSEIF(IPOL.EQ.2) THEN
                  WRITE(IUO1,89)
                ENDIF
                WRITE(IUO1,91)
                IF(IDICHR.GT.0) THEN
                  WRITE(IUO1,338)
                ENDIF
                WRITE(IUO1,91)
                WRITE(IUO1,92)
                WRITE(IUO1,90)
                WRITE(IUO1,43) THLUM,PHILUM
                IF((SPECTRO.EQ.'PHD').AND.(IMOD.EQ.1)) THEN
                  WRITE(IUO1,45)
                ENDIF
              ENDIF
C
              IF(INITL.EQ.2) THEN
                WRITE(IUO1,79) LI,LI-1,LI+1
                IF(I_SO.EQ.1) THEN
                  WRITE(IUO1,80) S_O
                ENDIF
                DO JE=1,NE
                  DO JEM=1,NEMET
                    JTE=IEMET(JEM)
                    IF(ISPIN.EQ.0) THEN
                      WRITE(IUO1,111)
     1                       JTE,(RHOK(JE,JTE,NNL,IDQ,1),IDQ=1,2)
                      IF(ITL.EQ.0) THEN
                        WRITE(IUO1,444)
     2                       JTE,(DLT(JE,JTE,NNL,IDQ),IDQ=1,2)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ELSEIF(INITL.EQ.-1) THEN
                WRITE(IUO1,82) LI,LI-1
                IF(I_SO.EQ.1) THEN
                  WRITE(IUO1,80) S_O
                ENDIF
                DO JE=1,NE
                  DO JEM=1,NEMET
                    JTE=IEMET(JEM)
                    IF(ISPIN.EQ.0) THEN
                      WRITE(IUO1,113) JTE,RHOK(JE,JTE,NNL,1,1)
                      IF(ITL.EQ.0) THEN
                        WRITE(IUO1,445) JTE,DLT(JE,JTE,NNL,1)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ELSEIF(INITL.EQ.1) THEN
                WRITE(IUO1,82) LI,LI+1
                IF(I_SO.EQ.1) THEN
                  WRITE(IUO1,80) S_O
                ENDIF
                DO JE=1,NE
                  DO JEM=1,NEMET
                    JTE=IEMET(JEM)
                    IF(ISPIN.EQ.0) THEN
                      WRITE(IUO1,113) JTE,RHOK(JE,JTE,NNL,2,1)
                      IF(ITL.EQ.0) THEN
                        WRITE(IUO1,445) JTE,DLT(JE,JTE,NNL,2)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
C
              IF(I_AT.EQ.0) THEN
                IF(INV(J_EL).EQ.0) THEN
                  IF(NDIF.EQ.1) THEN
                    IF(ISPHER.EQ.1) THEN
                      WRITE(IUO1,83)
                    ELSEIF(ISPHER.EQ.0) THEN
                      WRITE(IUO1,84)
                    ENDIF
                  ELSE
                    IF(ISPHER.EQ.0) THEN
                      WRITE(IUO1,97) NDIF
                    ELSE
                      WRITE(IUO1,98) NDIF
                    ENDIF
                  ENDIF
                ELSE
                  IF(ISPHER.EQ.0) THEN
                    WRITE(IUO1,122)
                  ELSE
                    WRITE(IUO1,120)
                  ENDIF
                ENDIF
              ELSE
                IF(ISPHER.EQ.0) THEN
                  WRITE(IUO1,85)
                ELSE
                  WRITE(IUO1,86)
                ENDIF
              ENDIF
C
            ENDIF
 222        CONTINUE
          ENDIF
C
C----------  LEED case (LED)  ----------
C
          IF(SPECTRO.EQ.'LED') THEN
            IF(IPHI.EQ.1) THEN
              IF(STEREO.EQ.' NO') THEN
                WRITE(IUO1,252)
              ELSE
                WRITE(IUO1,258)
              ENDIF
            ENDIF
            IF(ITHETA.EQ.1) WRITE(IUO1,254)
            IF(IE.EQ.1) WRITE(IUO1,256)
            IF(IPRINT.GT.0) THEN
              WRITE(IUO1,92)
              WRITE(IUO1,91)
              IF(ISPIN.EQ.0) THEN
                WRITE(IUO1,335)
              ELSE
                WRITE(IUO1,336)
              ENDIF
              WRITE(IUO1,91)
              IF(IPOTC.EQ.0) THEN
                WRITE(IUO1,339)
              ELSE
                WRITE(IUO1,334)
              ENDIF
              WRITE(IUO1,91)
              WRITE(IUO1,92)
              WRITE(IUO1,260)
              WRITE(IUO1,261) THLUM,PHILUM
              IF((SPECTRO.EQ.'LED').AND.(IMOD.EQ.1)) THEN
                WRITE(IUO1,45)
              ENDIF
C
              IF(I_AT.EQ.0) THEN
                IF(INV(J_EL).EQ.0) THEN
                  IF(NDIF.EQ.1) THEN
                    IF(ISPHER.EQ.1) THEN
                      WRITE(IUO1,83)
                    ELSEIF(ISPHER.EQ.0) THEN
                      WRITE(IUO1,84)
                    ENDIF
                  ELSE
                    IF(ISPHER.EQ.0) THEN
                      WRITE(IUO1,97) NDIF
                    ELSE
                      WRITE(IUO1,98) NDIF
                    ENDIF
                  ENDIF
                ELSE
                  IF(ISPHER.EQ.0) THEN
                    WRITE(IUO1,122)
                  ELSE
                    WRITE(IUO1,120)
                  ENDIF
                ENDIF
              ELSE
                IF(ISPHER.EQ.0) THEN
                  WRITE(IUO1,85)
                ELSE
                  WRITE(IUO1,86)
                ENDIF
              ENDIF
C
            ENDIF
          ENDIF
C
C----------  Auger diffraction case (AED)  ----------
C
          IF((SPECTRO.EQ.'AED').OR.(SPECTRO.EQ.'APC')) THEN
            IF(SPECTRO.EQ.'AED') THEN
              IF(IPHI_A.EQ.1) THEN
                IF(STEREO.EQ.' NO') THEN
                  WRITE(IUO1,235)
                ELSE
                  WRITE(IUO1,247)
                ENDIF
              ENDIF
              IF(ITHETA_A.EQ.1) WRITE(IUO1,244)
              IF(I_TEST_A.EQ.1) WRITE(IUO1,234)
            ENDIF
            IF((SPECTRO.EQ.'APC').AND.(J_EL.EQ.2)) THEN
              WRITE(IUO1,419)
              WRITE(IUO1,18)
            ENDIF
            IF((SPECTRO.EQ.'AED').OR.(J_EL.EQ.2)) THEN
              IF(IPRINT.GT.0) THEN
                WRITE(IUO1,92)
                WRITE(IUO1,91)
                IF(ISPIN.EQ.0) THEN
                  WRITE(IUO1,335)
                ELSE
                  WRITE(IUO1,336)
                ENDIF
                WRITE(IUO1,91)
                IF(IPOTC_A.EQ.0) THEN
                  WRITE(IUO1,339)
                ELSE
                  WRITE(IUO1,334)
                ENDIF
                WRITE(IUO1,91)
                WRITE(IUO1,92)
                WRITE(IUO1,95) AUGER
                CALL AUGER_MULT
                IF(I_MULT.EQ.0) THEN
                  WRITE(IUO1,154)
                ELSE
                  WRITE(IUO1,155) MULTIPLET
                ENDIF
C
                DO JEM=1,NEMET
                  JTE=IEMET(JEM)
                  WRITE(IUO1,112) JTE
                  DO LE=LE_MIN,LE_MAX
                    WRITE(IUO1,119) LE
                    LA_MIN=L_BOUNDS(LE,1)
                    LA_MAX=L_BOUNDS(LE,2)
                    DO LA=LA_MIN,LA_MAX
                      IF(ISPIN.EQ.0) THEN
                        WRITE(IUO1,115) LA,RHOK_A(LE,JTE,LA,1,1),
     1                                   RHOK_A(LE,JTE,LA,2,1)
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
C
                IF(I_AT.EQ.0) THEN
                  IF(INV(J_EL).EQ.0) THEN
                    IF(NDIF.EQ.1) THEN
                      IF(ISPHER.EQ.1) THEN
                        WRITE(IUO1,83)
                      ELSEIF(ISPHER.EQ.0) THEN
                        WRITE(IUO1,84)
                      ENDIF
                    ELSE
                      IF(ISPHER.EQ.0) THEN
                        WRITE(IUO1,97) NDIF
                      ELSE
                        WRITE(IUO1,98) NDIF
                      ENDIF
                    ENDIF
                  ELSE
                    IF(ISPHER.EQ.0) THEN
                      WRITE(IUO1,122)
                    ELSE
                      WRITE(IUO1,120)
                    ENDIF
                  ENDIF
                ELSE
                  IF(ISPHER.EQ.0) THEN
                    WRITE(IUO1,85)
                  ELSE
                    WRITE(IUO1,86)
                  ENDIF
                ENDIF
C
              ENDIF
            ENDIF
          ENDIF
C
C
C----------  Resonant Elastic X-ray Scattering case (RES)  ----------
C
          IF((SPECTRO.EQ.'RES')) THEN
            IF(IPHI.EQ.1) THEN
              IF(STEREO.EQ.' NO') THEN
                WRITE(IUO1,32)
              ELSE
                WRITE(IUO1,243)
              ENDIF
            ENDIF
            IF(ITHETA.EQ.1) WRITE(IUO1,46)
            IF(IE.EQ.1) WRITE(IUO1,58)
            IF(INITL.EQ.0) WRITE(IUO1,118)
            IF(I_TEST.EQ.1) WRITE(IUO1,234)
C
            IF(J_EL.EQ.2) GOTO 223
            IF(IPRINT.GT.0) THEN
              WRITE(IUO1,92)
              WRITE(IUO1,91)
              IF(ISPIN.EQ.0) THEN
                WRITE(IUO1,335)
              ELSE
                WRITE(IUO1,336)
              ENDIF
              WRITE(IUO1,91)
              IF(IPOTC.EQ.0) THEN
                WRITE(IUO1,339)
              ELSE
                WRITE(IUO1,334)
              ENDIF
              WRITE(IUO1,91)
C
C
C    INITL  = SWITCH TO SELECT THE CHANNELS IN THE FINAL STATE: 
C             (=-1) LF = LI - 1
C             (=0)  LF = LI (FOR THE PREPONDERANT AUGER CHANNEL)
C             (=1)  LF = LI + 1
C             (=2)  BOTH (-1) AND (1) (FULL DIPOLE SELECTION RULES)
C             (=4)  DIPOLE + QUADRUPOLE CHANNELS
C
              IF(INITL.NE.0) THEN
                WRITE(IUO1,337)
                WRITE(IUO1,91)
                IF(IPOL.EQ.0) THEN
                  WRITE(IUO1,88)
                ELSEIF(ABS(IPOL).EQ.1) THEN
                  WRITE(IUO1,87)
                ELSEIF(IPOL.EQ.2) THEN
                  WRITE(IUO1,89)
                ENDIF
                WRITE(IUO1,91)
                IF(IDICHR.GT.0) THEN
                  WRITE(IUO1,338)
                ENDIF
                WRITE(IUO1,91)
                WRITE(IUO1,92)
                WRITE(IUO1,90)
                WRITE(IUO1,43) THLUM,PHILUM
                IF((IMOD.EQ.1)) THEN
                  WRITE(IUO1,45)
                ENDIF
              ENDIF
C
              IF(INITL.EQ.2) THEN
                WRITE(IUO1,79) LI,LI-1,LI+1
                IF(I_SO.EQ.1) THEN
                  WRITE(IUO1,80) S_O
                ENDIF
                DO JE=1,NE
                  DO JEM=1,NEMET
                    JTE=IEMET(JEM)
                    IF(ISPIN.EQ.0) THEN
                      WRITE(IUO1,111)
     1                       JTE,(RHOK(JE,JTE,NNL,IDQ,1),IDQ=1,2)
                      IF(ITL.EQ.0) THEN
                        WRITE(IUO1,444)
     1                       JTE,(DLT(JE,JTE,NNL,IDQ),IDQ=1,2)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ELSEIF(INITL.EQ.-1) THEN
                WRITE(IUO1,82) LI,LI-1
                IF(I_SO.EQ.1) THEN
                  WRITE(IUO1,80) S_O
                ENDIF
                DO JE=1,NE
                  DO JEM=1,NEMET
                    JTE=IEMET(JEM)
                    IF(ISPIN.EQ.0) THEN
                      WRITE(IUO1,113) JTE,RHOK(JE,JTE,NNL,1,1)
                      IF(ITL.EQ.0) THEN
                        WRITE(IUO1,445) JTE,DLT(JE,JTE,NNL,1)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ELSEIF(INITL.EQ.1) THEN
                WRITE(IUO1,82) LI,LI+1
                IF(I_SO.EQ.1) THEN
                  WRITE(IUO1,80) S_O
                ENDIF
                DO JE=1,NE
                  DO JEM=1,NEMET
                    JTE=IEMET(JEM)
                    IF(ISPIN.EQ.0) THEN
                      WRITE(IUO1,113) JTE,RHOK(JE,JTE,NNL,2,1)
                      IF(ITL.EQ.0) THEN
                        WRITE(IUO1,445) JTE,DLT(JE,JTE,NNL,2)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ELSEIF(INITL.EQ.4) THEN
                WRITE(IUO1,78) LI,LI-2,LI-1,LI,LI+1,LI+2
                IF(I_SO.EQ.1) THEN
                  WRITE(IUO1,80) S_O
                ENDIF
                DO JE=1,NE
                  DO JEM=1,NEMET
                    JTE=IEMET(JEM)
                    IF(ISPIN.EQ.0) THEN
                      WRITE(IUO1,111)
     1                      JTE,(RHOK(JE,JTE,NNL,IDQ,1),IDQ=1,5)
                      IF(ITL.EQ.0) THEN
                        WRITE(IUO1,444)
     1                      JTE,(DLT(JE,JTE,NNL,1),IDQ=1,5)
                      ENDIF
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF
C
              IF(I_AT.EQ.0) THEN
                IF(INV(J_EL).EQ.0) THEN
                  IF(NDIF.EQ.1) THEN
                    IF(ISPHER.EQ.1) THEN
                      WRITE(IUO1,83)
                    ELSEIF(ISPHER.EQ.0) THEN
                      WRITE(IUO1,84)
                    ENDIF
                  ELSE
                    IF(ISPHER.EQ.0) THEN
                      WRITE(IUO1,97) NDIF
                    ELSE
                      WRITE(IUO1,98) NDIF
                    ENDIF
                  ENDIF
                ELSE
                  IF(ISPHER.EQ.0) THEN
                    WRITE(IUO1,122)
                  ELSE
                    WRITE(IUO1,120)
                  ENDIF
                ENDIF
              ELSE
                IF(ISPHER.EQ.0) THEN
                  WRITE(IUO1,85)
                ELSE
                  WRITE(IUO1,86)
                ENDIF
              ENDIF
C
            ENDIF
 223        CONTINUE
          ENDIF
C
C..........  Check of the dimensioning of the treatment routine  ..........
C
          CALL STOP_TREAT(NFICHLEC,NPLAN,NEMET,NE,NTHETA,NTHETA_A,
     1                  NPHI,NPHI_A,ISOM,I_EXT,I_EXT_A,SPECTRO)
C
C..........     Call of the subroutine performing either     ..........
C..........  the PhD, LEED, AED, EXAFS or APECS calculation  ..........
C
 566      IF(ISPIN.EQ.0) THEN
            IF(SPECTRO.EQ.'PHD') THEN
C              CALL PHDDIF_MI(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOK,
C     1                     NATCLU,NFICHLEC,JFICH,NP)
            ELSEIF(SPECTRO.EQ.'LED') THEN
c            CALL LEDDIF_MI(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOK,
c     1                     NATCLU,NFICHLEC,JFICH,NP)
            ELSEIF(SPECTRO.EQ.'AED') THEN
c            CALL AEDDIF_MI(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,
c     1                     RHOK_A,NATCLU,NFICHLEC,JFICH,NP,LE_MIN,
c     2                     LE_MAX)
            ELSEIF(SPECTRO.EQ.'XAS') THEN
c            CALL XASDIF_MI(NPLAN,VAL,ZEM,IPHA,RHOK,NFICHLEC,JFICH,NP)
            ELSEIF(SPECTRO.EQ.'APC') THEN
c            IF(J_EL.EQ.1) THEN
c              CALL PHDDIF_MI(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOK,
c     1                    NATCLU,NFICHLEC,JFICH,NP)
c            ELSEIF(J_EL.EQ.2) THEN
c              CALL AEDDIF_MI(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOK_A,
c     1                    NATCLU,NFICHLEC,JFICH,NP,LE_MIN,LE_MAX)
c            ENDIF
            ELSEIF(SPECTRO.EQ.'RES') THEN
              CALL RESDIF_MI(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOK,
     1                     NATCLU,NFICHLEC,JFICH,NP)
            ENDIF
          ELSEIF(ISPIN.EQ.1) THEN
c          IF(SPECTRO.EQ.'PHD') THEN
c            CALL PHDDIF_SP(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOK,
c     1                     NATCLU,NFICHLEC,JFICH,NP)
c          ELSEIF(SPECTRO.EQ.'AED') THEN
c            CALL AEDDIF_SP
c          ELSEIF(SPECTRO.EQ.'XAS') THEN
c            CALL XASDIF_SP
c          ENDIF
            continue
          ENDIF
C
C..........        End of the MS calculation :        ..........
C..........  direct exit or treatment of the results  ..........
C
C
C..........  End of the loop on the electrons  ..........
C
        ENDDO
C
        IF(SPECTRO.EQ.'PHD') THEN
          IF(IPHI.EQ.1) THEN
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,52)
            ELSE
              WRITE(IUO1,249)
            ENDIF
          ENDIF
          IF(ITHETA.EQ.1) WRITE(IUO1,49)
          IF(IE.EQ.1) WRITE(IUO1,59)
        ELSEIF(SPECTRO.EQ.'LED') THEN
          IF(IPHI.EQ.1) THEN
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,253)
            ELSE
              WRITE(IUO1,259)
            ENDIF
          ENDIF
          IF(ITHETA.EQ.1) WRITE(IUO1,255)
          IF(IE.EQ.1) WRITE(IUO1,257)
        ELSEIF(SPECTRO.EQ.'XAS') THEN
          WRITE(IUO1,51)
        ELSEIF(SPECTRO.EQ.'AED') THEN
          IF(IPHI_A.EQ.1) THEN
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,237)
            ELSE
              WRITE(IUO1,250)
            ENDIF
          ENDIF
          IF(ITHETA_A.EQ.1) WRITE(IUO1,238)
        ELSEIF(SPECTRO.EQ.'APC') THEN
          IF(IPHI.EQ.1) THEN
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,239)
            ELSE
              WRITE(IUO1,251)
            ENDIF
          ENDIF
          IF(ITHETA.EQ.1) WRITE(IUO1,240)
        ELSEIF(SPECTRO.EQ.'RES') THEN
          IF(IPHI.EQ.1) THEN
            IF(STEREO.EQ.' NO') THEN
              WRITE(IUO1,53)
            ELSE
              WRITE(IUO1,262)
            ENDIF
          ENDIF
          IF(ITHETA.EQ.1) WRITE(IUO1,48)
          IF(IE.EQ.1) WRITE(IUO1,59)
        ENDIF
C
        CLOSE(ICOM)
        IF((NFICHLEC.GT.1).AND.(ISOM.NE.0)) THEN
          WRITE(IUO1,562)
        ENDIF
        IF(ISOM.EQ.0) CLOSE(IUO2)
        IF((ISOM.EQ.0).AND.(NFICHLEC.NE.1)) CLOSE(IUO1)
C
C..........  End of the loop on the data files  ..........
C
      ENDDO
C
      IF(ISOM.NE.0) THEN
        JFF=1
        IF(ISPIN.EQ.0) THEN
          IF(SPECTRO.NE.'XAS') THEN
            CALL TREAT_PHD(ISOM,NFICHLEC,JFF,NP)
          ELSE
c            CALL TREAT_XAS(ISOM,NFICHLEC,NP)
          ENDIF
        ELSEIF(ISPIN.EQ.1) THEN
c          IF((SPECTRO.EQ.'PHD').OR.(SPECTRO.EQ.'AED')) THEN
c            CALL TREAT_PHD_SP(ISOM,NFICHLEC,JFF,NP)
c          ELSEIF(SPECTRO.EQ.'XAS') THEN
c            CALL TREAT_XAS_SP(ISOM,NFICHLEC,NP)
c          ENDIF
          continue
        ENDIF
      ENDIF
C
      IF((ISOM.NE.0).OR.(NFICHLEC.EQ.1)) CLOSE(IUO1)
      IF(ISOM.NE.0) CLOSE(IUO2)
      STOP
C
   1  WRITE(IUO1,60)
      STOP
   2  WRITE(IUO1,61)
      STOP
  55  WRITE(IUO1,65)
      STOP
  56  WRITE(IUO1,64)
      STOP
  74  WRITE(IUO1,75)
      STOP
  99  WRITE(IUO1,100)
      STOP
 180  WRITE(IUO1,181)
      STOP
 182  WRITE(IUO1,183)
      STOP
 184  WRITE(IUO1,185)
      STOP
 504  WRITE(IUO1,505)
      STOP
 510  WRITE(IUO1,511) IUI4
      STOP
 514  WRITE(IUO1,515)
      STOP
 516  WRITE(IUO1,517)
      STOP
 518  WRITE(IUO1,519)
      WRITE(IUO1,889)
      STOP
 520  WRITE(IUO1,521)
      STOP
 540  WRITE(IUO1,541)
      STOP
 550  WRITE(IUO1,551)
      STOP
 570  WRITE(IUO1,571)
      STOP
 580  WRITE(IUO1,581)
      STOP
 590  WRITE(IUO1,591)
      STOP
 600  WRITE(IUO1,601)
      STOP
 602  WRITE(IUO1,603)
      STOP
 604  WRITE(IUO1,605)
      STOP
 606  WRITE(IUO1,607)
      STOP
 608  WRITE(IUO1,609)
      STOP
 610  WRITE(IUO1,611)
      STOP
 614  WRITE(IUO1,615) NB_AT
      STOP
 620  WRITE(IUO1,621) LE_MAX
      STOP
 630  WRITE(IUO1,631)
      STOP
 890  WRITE(IUO1,891)
      STOP
 895  WRITE(IUO1,896) JA1,JA2
C
   3  FORMAT(5(5X,I4))
   7  FORMAT(3X,F9.4,1X,F9.4,5X,F12.9,5X,F12.9)
   9  FORMAT(3X,F9.4,1X,F9.4,5X,E12.6,5X,E12.6)
  17  FORMAT(12X,'ATOM NUMBER ',I4,10X,'CORRESPONDING TRANSLATIONS ',
     1': (',I3,',',I3,',',I3,')')
  18  FORMAT('            ',/)
  20  FORMAT(/,7X,'ATOM OF TYPE ',I2,' AND OF NUMBER ',I5)
  21  FORMAT(17X,'COORDINATES IN THE TOTAL CLUSTER : (',F7.3,',',
     1       F7.3,',',F7.3,')')
  22  FORMAT(22X,'THIS ATOM HAS BEEN SUPRESSED IN THE REDUCED CLUSTER')
  23  FORMAT(17X,'COORDINATES IN THE REDUCED CLUSTER :(',F7.3,',',
     1       F7.3,',',F7.3,')',5X,'NEW NUMBER : ',I4)
  24  FORMAT(///,29X,'CONTENTS OF THE REDUCED CLUSTER :',/)
  26  FORMAT(28X,I4,' ATOMS OF TYPE ',I2)
  29  FORMAT(/,20X,'THE Z POSITION OF PLANE ',I3,' IS : ',F6.3)
  30  FORMAT(///,23X,'THE ABSORBING ATOMS ARE OF TYPE :',/)
  31  FORMAT(38X,10(I2,3X),//)
  32  FORMAT(/////,'########## BEGINNING ',
     1'OF THE AZIMUTHAL RESONANT ELESTIC X-RAY SCATTERING CALCULATION',
     2'#####',/////)
  34  FORMAT(//,2X,'PLANE No ',I3,3X,'THE ABSORBER OF TYPE ',
     1I2,' IS POSITIONED AT (',F7.3,',',F7.3,',',F7.3,')')
  35  FORMAT(/////,'########## BEGINNING ',
     1'OF THE AZIMUTHAL PHOTOELECTRON DIFFRACTION CALCULATION #####',
     2'#####',/////)
  36  FORMAT(/////,'##########              BEGINNING ',
     1'OF THE EXAFS CALCULATION              ##########',/////)
  37  FORMAT(/////,'++++++++++++++++++++',
     1'     NUMBERING OF THE ATOMS GENERATED     +++++++++++++++++++')
  38  FORMAT(///,30X,'TRANSLATION LEVEL : ',I2,///)
  39  FORMAT(///,'++++++++++++++++++++++++++++++++++++++++++++++++',
     1'++++++++++++++++++++++++++++++++',/////)
  40  FORMAT(/////,'======================',
     1'     CONTENTS OF THE REDUCED CLUSTER     ======================',
     2      ///)
  41  FORMAT(///,'====================================================',
     1'============================',/////)
  43  FORMAT(14X,'TH_LIGHT = ',F6.2,' DEGREES',5X,'PHI_LIGHT = ',
     1       F6.2,' DEGREES')
  44  FORMAT(/////,'########## BEGINNING ',
     1'OF THE POLAR PHOTOELECTRON DIFFRACTION CALCULATION #####',
     2'#####',/////)
  45  FORMAT(14X,' (WHEN THE DETECTOR IS ALONG ',
     1       'THE NORMAL TO THE SURFACE)')
  46  FORMAT(/////,'########## BEGINNING ',
     1'OF THE POLAR RESONANT ELESTIC X-RAY SCATTERING CALCULATION ',
     2'##########',/////)
  48  FORMAT(/////,'########## END OF THE ',
     1'POLAR RESONANT ELESTIC X-RAY SCATTERING CALCULATION ##########')
  49  FORMAT(/////,'########## END OF THE ',
     1'POLAR PHOTOELECTRON DIFFRACTION CALCULATION ##########')
  50  FORMAT(///,22X,'THE CLUSTER IS COMPOSED OF ',I2,' PLANES :')
  51  FORMAT(/////,'##########                END OF THE ',
     1'EXAFS CALCULATION                ##########')
  52  FORMAT(/////,'########## END OF THE ',
     1'AZIMUTHAL PHOTOELECTRON DIFFRACTION CALCULATION #####',
     2'#####')
  53  FORMAT(/////,'########## END OF THE ',
     1'AZIMUTHAL RESONANT ELESTIC X-RAY SCATTERING CALCULATION #####',
     2'#####')
  57  FORMAT(///,27X,'CALCULATION OF THE SCATTERING FACTOR DONE')
  58  FORMAT(/////,'########## BEGINNING ',
     1'OF THE FINE STRUCTURE OSCILLATIONS CALCULATION #####',
     2'#####',/////)
  59  FORMAT(/////,'########## END OF THE ',
     1'FINE STRUCTURE OSCILLATIONS CALCULATION #####',
     2'#####')
  60  FORMAT(///,'<<<<<<<<<<  (NAT,NE,NEMET) > (NATP_M,NE_M,',
     1       'NEMET_M) - CHECK THE DIMENSIONING >>>>>>>>>>')
  61  FORMAT(///,22X,' <<<<<<<<<<  THIS STRUCTURE DOES NOT EXIST ',
     1'  >>>>>>>>>>')
  64  FORMAT(///,4X,' <<<<<<<<<<  NIV IS TOO SMALL, THE REDUCED ',
     1'CLUSTER HAS NOT CONVERGED YET  >>>>>>>>>>')
  65  FORMAT(///,4X,' <<<<<<<<<<  ONLY ONE OF THE VALUES IPHI,ITHETA ',
     1'ET IE CAN BE EQUAL TO 1  >>>>>>>>>>')
  75  FORMAT(///,8X,' <<<<<<<<<<  CHANGE THE DIMENSIONING OF PCREL ',
     1'IN MAIN ET READ_DATA  >>>>>>>>>>')
  78  FORMAT(//,18X,'INITIAL STATE L = ',I1,5X,'FINAL STATES L = ',
     1I1,4(',',I1),/)
  79  FORMAT(//,18X,'INITIAL STATE L = ',I1,5X,'FINAL STATES L = ',
     1I1,',',I1,/)
  80  FORMAT(15X,'(SPIN-ORBIT COMPONENT OF THE INITIAL CORE STATE : ',
     1          A3,')',//)
  81  FORMAT(18X,'(BOTH SPIN-ORBIT COMPONENTS TAKEN INTO ACCOUNT)')
  82  FORMAT(//,21X,'INITIAL STATE L = ',I1,5X,'FINAL STATE L = ',I1)
  83  FORMAT(//,32X,'(SPHERICAL WAVES)')
  84  FORMAT(//,34X,'(PLANE WAVES)')
  85  FORMAT(//,26X,'(PLANE WAVES - ATOMIC CASE)')
  86  FORMAT(//,24X,'(SPHERICAL WAVES - ATOMIC CASE)')
  87  FORMAT(24X,'+     LINEARLY POLARIZED LIGHT      +')
  88  FORMAT(24X,'+        NON POLARIZED LIGHT        +')
  89  FORMAT(24X,'+    CIRCULARLY POLARIZED LIGHT     +')
  90  FORMAT(////,31X,'POSITION OF THE LIGHT :',/)
  91  FORMAT(24X,'+',35X,'+')
  92  FORMAT(24X,'+++++++++++++++++++++++++++++++++++++')
  94  FORMAT(//,2X,'PLANE No ',I3,3X,'NO ABSORBER OF TYPE ',I2,
     1' IS PRESENT IN THIS PLANE')
  95  FORMAT(////,31X,'AUGER LINE :',A6,//)
  97  FORMAT(///,19X,'(PLANE WAVES MULTIPLE SCATTERING - ORDER ',I1,
     1       ')')
  98  FORMAT(///,17X,'(SPHERICAL WAVES MULTIPLE SCATTERING - ORDER ',
     1       I1,')')
 100  FORMAT(///,8X,'<<<<<<<<<<  WRONG NAME FOR THE INITIAL STATE',
     1      '  >>>>>>>>>>')
 101  FORMAT(24X,I3,24X,I3)
 102  FORMAT(A1)
 103  FORMAT(31X,F7.2)
 104  FORMAT(29X,F8.5,4X,F8.5,7X,F8.5,4X,F8.5)
 105  FORMAT(8X,f8.5,1X,f8.5,4X,f8.5,1X,f8.5,4X,f8.5,1X,f8.5,
     1       4X,f8.5,1X,f8.5,4X,f8.5,1X,f8.5)
 106  FORMAT(12X,I3,12X,I3,12X,I3)
 107  FORMAT(5X,I2,5X,I2,5X,I2)
 108  FORMAT(19X,I2,8X,F8.5,1X,F8.5,4X,F8.5,1X,F8.5)
 109  FORMAT(5X,I2,12X,I2,11X,I2)
 110  FORMAT(16X,'RADIAL MATRIX ELEMENTS FOR THE ABSORBER OF TYPE ',I2,
     1       ' :',/,22X,'(THE SPIN DOUBLET IS GIVEN AS : OUT/IN)',//)
 111  FORMAT(6X,'RADIAL MATRIX ELEMENTS FOR THE ABSORBER OF TYPE ',
     1       I2,' : (',F8.5,',',F8.5,')',4(/,59X,'(',F8.5,',',F8.5,')'))
 112  FORMAT(6X,'RADIAL MATRIX ELEMENTS FOR THE ABSORBER OF TYPE ',
     1       I2,' : ',/,8X,'(LE : ALLOWED VALUES FOR ESCAPING AUGER',
     2      ' ELECTRON)',/,
     2       8X,'(L  : INTERNAL VALUE THAT WILL BE SUMMED ON)',//)
 113  FORMAT(6X,'RADIAL MATRIX ELEMENT FOR THE ABSORBER OF ',
     *       'TYPE  ',I2,' : (',F8.5,',',F8.5,')')
 114  FORMAT(/)
 115  FORMAT(15X,'L = ',I2,5X,'(',F8.5,',',F8.5,')',5X,
     1       '(',F8.5,',',F8.5,')')
 117  FORMAT(12X,I2,5X,I2)
 118  FORMAT(/,37X,'AUGER ELECTRON DIFFRACTION',/)
 119  FORMAT(10X,'LE = ',I2,11X,'DIRECT INTEGRAL',8X,
     1       'EXCHANGE INTEGRAL')
 120  FORMAT(///,15X,'(SPHERICAL WAVES MULTIPLE SCATTERING - MATRIX ',
     1       'INVERSION)')
 122  FORMAT(///,17X,'(PLANE WAVES MULTIPLE SCATTERING - MATRIX ',
     1       'INVERSION)')
 125  FORMAT(11X,A2,5X,I2,3F10.4,12X,I4)
 154  FORMAT(///,20X,'CALCULATION MADE FOR THE FULL AUGER LINE',
     1       ' ',/,' ',/,' ')
 155  FORMAT(///,20X,'CALCULATION MADE FOR THE ',A3,' MULTIPLET ',
     1      'LINE',' ',/,' ',/,' ')
 181  FORMAT(///,'<<<<<<<<<<  NAT OR NE DIFFERENT BETWEEN THE INPUT ',
     1          'AND PHASE SHIFTS FILES  >>>>>>>>>>')
 183  FORMAT(///,'<<<<<<<<<<  NAT OR NE DIFFERENT BETWEEN THE INPUT ',
     1          'AND RADIAL MATRIX ELEMENTS FILES  >>>>>>>>>>')
 185  FORMAT(///,'<<<<<<<<<<  LMAX > NL_M-1 IN THE PHASE SHIFTS ',
     1          'FILE  >>>>>>>>>>')
 234  FORMAT('  ----->  TEST CALCULATION : NO EXCITATION ',
     1           'MATRIX ELEMENTS TAKEN INTO ACCOUNT  <-----',///)
 235  FORMAT(/////,'########## BEGINNING ',
     1'OF THE AZIMUTHAL AUGER DIFFRACTION CALCULATION #####',
     2'#####',/////)
 236  FORMAT(/////,'########## BEGINNING ',
     1'OF THE AZIMUTHAL APECS DIFFRACTION CALCULATION #####',
     2'#####',/////)
 237  FORMAT(/////,'########## END ',
     1'OF THE AZIMUTHAL AUGER DIFFRACTION CALCULATION #####',
     2'#####',/////)
 238  FORMAT(/////,6X,'########## END ',
     1'OF THE POLAR AUGER DIFFRACTION CALCULATION #####',
     2'#####',/////)
 239  FORMAT(/////,'########## END ',
     1'OF THE AZIMUTHAL APECS DIFFRACTION CALCULATION #####',
     2'#####',/////)
 240  FORMAT(/////,6X,'########## END ',
     1'OF THE POLAR APECS DIFFRACTION CALCULATION #####',
     2'#####',/////)
 243  FORMAT(/////,'########## BEGINNING ',
     1'OF THE FULL ANGLE RESONANT ELESTIC X-RAY SCATTERING CALCULATION',
     2' ##########',/////)
 244  FORMAT(/////,6X,'########## BEGINNING ',
     1'OF THE POLAR AUGER DIFFRACTION CALCULATION #####',
     2'#####',/////)
 245  FORMAT(/////,6X,'########## BEGINNING ',
     1'OF THE POLAR APECS DIFFRACTION CALCULATION #####',
     2'#####',/////)
 246  FORMAT(/////,'########## BEGINNING ',
     1'OF THE FULL ANGLE PHOTOELECTRON DIFFRACTION CALCULATION ',
     2'##########',/////)
 247  FORMAT(/////,'########## BEGINNING ',
     1'OF THE FULL ANGLE AUGER DIFFRACTION CALCULATION ',
     2'##########',/////)
 248  FORMAT(/////,'########## BEGINNING ',
     1'OF THE FULL ANGLE APECS DIFFRACTION CALCULATION ',
     2'##########',/////)
 249  FORMAT(/////,'########## END OF THE ',
     1'FULL ANGLE PHOTOELECTRON DIFFRACTION CALCULATION #####',
     2'#####')
 250  FORMAT(/////,'########## END ',
     1'OF THE FULL ANGLE AUGER DIFFRACTION CALCULATION #####',
     2'#####',/////)
 251  FORMAT(/////,'########## END ',
     1'OF THE FULL ANGLE APECS DIFFRACTION CALCULATION #####',
     2'#####',/////)
 252  FORMAT(/////,'########## BEGINNING ',
     1'OF THE AZIMUTHAL LEED CALCULATION #####',
     2'#####',/////)
 253  FORMAT(/////,'########## END ',
     1'OF THE AZIMUTHAL LEED CALCULATION #####',
     2'#####',/////)
 254  FORMAT(/////,6X,'########## BEGINNING ',
     1'OF THE POLAR LEED CALCULATION #####',
     2'#####',/////)
 255  FORMAT(/////,6X,'########## END ',
     1'OF THE POLAR LEED CALCULATION #####',
     2'#####',/////)
 256  FORMAT(/////,5X,'########## BEGINNING ',
     1'OF THE ENERGY LEED CALCULATION #####',
     2'#####',/////)
 257  FORMAT(/////,5X,'########## END ',
     1'OF THE ENERGY LEED CALCULATION #####',
     2'#####',/////)
 258  FORMAT(/////,'########## BEGINNING ',
     1'OF THE FULL ANGLE LEED CALCULATION ',
     2'##########',/////)
 259  FORMAT(/////,'########## END OF THE ',
     1'FULL ANGLE LEED CALCULATION #####',
     2'#####')
 260  FORMAT(////,31X,'POSITION OF THE INITIAL BEAM :',/)
 261  FORMAT(14X,'TH_BEAM = ',F6.2,' DEGREES',5X,'PHI_BEAM = ',
     1       F6.2,' DEGREES')
 262  FORMAT(/////,'########## END OF THE ',
     1'FULL ANGLE RESONANT ELESTIC X-RAY SCATTERING CALCULATION #####',
     2'#####')
 334  FORMAT(24X,'+   COMPLEX POTENTIAL CALCULATION   +')
 335  FORMAT(24X,'+              STANDARD             +')
 336  FORMAT(24X,'+           SPIN-POLARIZED          +')
 337  FORMAT(24X,'+                WITH               +')
 338  FORMAT(24X,'+         IN DICHROIC MODE          +')
 339  FORMAT(24X,'+    REAL POTENTIAL CALCULATION     +')
 418  FORMAT(///,9X,'------------------------  FIRST ELECTRON : ',
     1       '------------------------')
 419  FORMAT(///,9X,'------------------------ SECOND ELECTRON : ',
     1       '------------------------')
 420  FORMAT(///,9X,'----------------------------------------------',
     1       '----------------------')
 444  FORMAT(12X,'PHASE SHIFTS FOR THE ABSORBER OF TYPE  ',I2,' : ',
     1       '(',F8.5,',',F8.5,')',4(/,56X,'(',F8.5,',',F8.5,')'))
 445  FORMAT(12X,'PHASE SHIFT FOR THE ABSORBER OF TYPE  ',I2,' : (',
     1       F8.5,',',F8.5,')')
 505  FORMAT(///,'<<<<<<<<<<  LI IS LARGER THAN LI_M - ',
     1       'CHECK THE DIMENSIONING  >>>>>>>>>>')
 511  FORMAT(///,'<<<<<<<<<<  NATCLU_M IN THE .inc FILE IS NOT ',
     1       'CONSISTENT WITH THE NUMBER OF ATOMS READ FROM UNIT ',I2,
     2       '  >>>>>>>>>>')
 515  FORMAT(///,'<<<<<<<<<<  INCOMPATIBILITY BETWEEN THE VALUES OF ',
     1       'NAT IN THE DATA AND CLUSTER FILES  >>>>>>>>>>')
 517  FORMAT(///,'<<<<<<<<<<  THERE ARE MISSING VALUES FOR THFWD AND ',
     1       'IBWD  >>>>>>>>>>')
 519  FORMAT(///,'<<<<<<<<<<  NATCLU_M IN THE .inc FILE IS NOT',
     1       ' CONSISTENT WITH THE NUMBER OF ATOMS GENERATED BY THE ',
     2       'CODE  >>>>>>>>>>')
 521  FORMAT(///,'<<<<<<<<<<  SPIN-ORBIT COMPONENT NOT CONSISTENT WITH',
     1       ' THE VALUE OF LI  >>>>>>>>>>')
 530  FORMAT(3X,F9.4,3X,F9.4,3X,F9.4)
 535  FORMAT(29X,F8.5,1X,F8.5)
 541  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES THFWD DOES NOT ',
     1       'CORRESPOND TO NAT  >>>>>>>>>>')
 543  FORMAT(5X,F12.9,5X,F12.9)
 549  FORMAT(//,14X,' No ',10X,'COORDINATES',9X,'TYPE',2X,
     2       'SNo',2X,'SYM',/)
 551  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES UJ2 DOES NOT ',
     1       'CORRESPOND TO NAT  >>>>>>>>>>')
 555  FORMAT(4(7X,I2))
 556  FORMAT(28X,4(I2,5X))
 557  FORMAT(13X,I4,3X,'(',F7.3,',',F7.3,',',F7.3,')',2X,I4,2X,I4,
     1        3X,A2)
 558  FORMAT(/////,18X,'CONTENTS OF THE CLUSTER READ FROM UNIT ',
     1       I2,' : ',/,20X,'READ IN ',A30,//,15X,'No',13X,'(X,Y,Z)',
     2       10X,'CLASS',1X,'ATOM',/)
 559  FORMAT(/////,25X,'CONTENTS OF THE CLUSTER GENERATED : ',//,
     1       14X,' No ',10X,'COORDINATES',9X,'TYPE',2X,'SNo',2X,'SYM',/)
 560  FORMAT(////,12X,'MAXIMAL VALUES OF L FOR THE ',I3,
     1       ' PROTOTYPICAL ATOMS : ',//)
 561  FORMAT(////,18X,'MAXIMAL VALUE OF L FOR THE ',
     1       'PROTOTYPICAL ATOM : ',//)
 562  FORMAT(///,'oooooooooooooooo',12X,'END OF THE INPUT DATA FILE',
     1       13X,'oooooooooooooooo',///)
 563  FORMAT(//,20X,'ENERGY POINT No ',I3,' :',/)
 571  FORMAT(///,'<<<<<<<<<<  THE NUMBER OF LINES ATBAS DOES NOT ',
     1       'CORRESPOND TO NAT  >>>>>>>>>>')
 581  FORMAT(///,'<<<<<<<<<<  LI OR IMOD NOT CONSISTENT BETWEEN ',
     1       'PHD AND AED FOR COINCIDENCE CALCULATION  >>>>>>>>>>')
 591  FORMAT(///,'<<<<<<<<<<  THE EXTERNAL DIRECTIONS FILE IS ',
     1       'NOT CONSISTENT WITH THE INPUT DATA FILE  >>>>>>>>>>')
 601  FORMAT(///,'<<<<<<<<<<  NO_ST_M IS TOO SMALL IN THE .inc FILE ',
     1       '>>>>>>>>>>',//)
 603  FORMAT(///,'<<<<<<<<<<  NSPIN_M OR NSPIN2_M IS TOO SMALL IN THE ',
     1       '.inc FILE  >>>>>>>>>>',//)
 605  FORMAT(///,'<<<<<<<<<<  NT_M IS TOO SMALL IN THE .inc FILE ',
     1       '>>>>>>>>>>',//)
 607  FORMAT(///,'<<<<<<<<<<  THE INITIAL STATE LI IN THE INPUT DATA  ',
     1       'FILE IS DIFFERENT FROM THAT IN THE RADIAL MATRIX ',
     2       'ELEMENTS FILE  >>>>>>>>>>',//)
 609  FORMAT(///,'<<<<<<<<<<  THE TWO TL FILE ARE NOT COMPATIBLE  ',
     1           '>>>>>>>>>>',//)
 611  FORMAT(///,3X,'<<<<<<<<<<  THE RADIAL FILE FOR THE AUGER ',
     1           'ELECTRON IS NOT COMPATIBLE   >>>>>>>>>>',/,
     2           3X,'<<<<<<<<<<  ',17X,'WITH THE INPUT DATA FILE  ',
     3           16X,'>>>>>>>>>>',//)
 613  FORMAT(///,'<<<<<<<<<<  NATP_M SHOULD BE AT LEAST ',I3,'  IN ',
     1       'THE DIMENSIONNING FILE  >>>>>>>>>>',//)
 615  FORMAT(///,'<<<<<<<<<<  NAT_EQ_M SHOULD BE AT LEAST ',I3,'  IN ',
     1       'THE DIMENSIONNING FILE  >>>>>>>>>>',//)
 621  FORMAT(///,'<<<<<<<<<<  LI_M SHOULD BE AT LEAST ',I3,'  IN ',
     1       'THE DIMENSIONNING FILE  >>>>>>>>>>',//)
 631  FORMAT(///,'<<<<<<<<<<      EXCURSIONS OF ANGLES SHOULD ',
     1       ' BE IDENTICAL      >>>>>>>>>>',/,'<<<<<<<<<<     ',
     2       'FOR BOTH ELECTRONS IN CLUSTER ROTATION MODE',
     3        '     >>>>>>>>>>',//)
 776  FORMAT(I2)
 777  FORMAT(A24)
 778  FORMAT(30X,I1)
 779  FORMAT(11X,A2,5X,I2,3F10.4,I5)
 782  FORMAT(/////,22X,'THE CLUSTER GENERATED CONSISTS OF : ',I4,
     1       ' ATOMS')
 889  FORMAT(/////,'<<<<<<<<<<  DECREASE NIV OR INCREASE',
     1       ' NATCLU_M  >>>>>>>>>>')
 891  FORMAT(/////,'<<<<<<<<<<  WRONG NAME FOR THE COORDINATES ''',
     1       'UNITS  >>>>>>>>>>')
 896  FORMAT(///,10X,'<<<<<<<<<<  ERROR IN THE COORDINATES OF THE',
     1       '  ATOMS  >>>>>>>>>>',/,10X,'<<<<<<<<<<     ATOMS ',I4,
     2       ' AND ',I4,' ARE IDENTICAL    >>>>>>>>>>')
C
      END
C
C=======================================================================
C
      SUBROUTINE AMAS(NIV,ATOME,COORD,VALZ,ISURF,COUPUR,ROT,IRE,
     1                NATYP,NBZ,NAT2,NCOUCH,NMAX)
C
C  This routine generates a cluster from the knowledge of its
C                          lattice vectors
C
      INCLUDE 'spec.inc'
C
      COMMON /ADSORB/ IADS,NATA,NADS1,NADS2,NADS3,ADS(3,900),NCOUCH1
      COMMON /BASES/ ATBAS(3*NATP_M),VECBAS(9)
      COMMON /MILLER/ IM1,IM2,IM3,IM4,IVG0,IVN(3)
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /RESEAU/ NCRIST,NCENTR,IBAS,NAT,A,BSURA,CSURA,UNIT
C
      DIMENSION VALZ(NATCLU_M)
      DIMENSION ROT(3,3),IRE(NATCLU_M,2),NATYP(NATM),ITA(NATCLU_M)
      DIMENSION ATOME(3,NATCLU_M),ATRSU(3,NATCLU_M),COORD(3,NATCLU_M)
      DIMENSION ROTINV(3,3),XINIT(3,1),XFIN(3,1)
C
      CHARACTER*3 UNIT
C
      NCOUCH=0
      WRITE(IUO1,10) ISURF
  10  FORMAT(//,18X,'ATOM (0,0,0) ON THE SURFACE PLANE IS OF TYPE ',I2)
      NBZ=0
      CALL INVMAT(ROT,ROTINV)
      IF(IVG0.EQ.0) THEN
        CALL CHBASE(NATP_M,ATBAS)
      ENDIF
      NB1=0
      NB2=0
      DO NTYP=1,NAT
        NBAT=0
        DO NUM=1,NMAX
          NB1=NB1+1
          IRE(NB1,1)=0
          IRE(NB1,2)=0
          IF(IVG0.LE.1) THEN
            CALL NUMAT(NUM,NIV,IA,IB,IC)
          ELSE
            BSURA=1.
            CSURA=1.
          ENDIF
          IF(IVG0.LE.1) THEN
            XA=FLOAT(IA)
            XB=FLOAT(IB)
            XC=FLOAT(IC)
          ELSEIF(IVG0.EQ.2) THEN
            XA=FLOAT(NUM-1)
            XB=FLOAT(NUM-1)
            XC=FLOAT(NUM-1)
          ENDIF
          IF(IVG0.EQ.1) THEN
            IF(IVN(1).EQ.0) THEN
              ITA(NUM)=IA
            ELSEIF(IVN(2).EQ.0) THEN
              ITA(NUM)=IB
            ELSEIF(IVN(3).EQ.0) THEN
              ITA(NUM)=IC
            ENDIF
            IF((ITA(NUM).EQ.ITA(NUM-1)).AND.(NUM.GT.1)) GOTO 30
          ENDIF
          DO J=1,3
            K=J+3*(NTYP-1)
            O=ATBAS(K)
            ATOME(J,NB1)=O+XA*VECBAS(J)+XB*VECBAS(J+3)+XC*VECBAS(J+6)
          ENDDO
          DO I=1,3
            M=I+3*(ISURF-1)
            XINIT(I,1)=ATOME(I,NB1)-ATBAS(M)
          ENDDO
          CALL MULMAT(ROTINV,3,3,XINIT,3,1,XFIN)
          DO I=1,3
            ATRSU(I,NB1)=XFIN(I,1)
          ENDDO
          CALL TEST1(COUPUR,NB1,NB2,ATRSU,COORD,VALZ,NBAT,IRE,NBZ)
  30      CONTINUE
        ENDDO
        NATYP(NTYP)=NBAT
      ENDDO
      IF(IADS.GE.1) THEN
        N0=NBZ
        DO JADS=1,NADS1
          NB1=NB1+1
          DO I=1,3
            COORD(I,NB1)=ADS(I,JADS)
          ENDDO
          N1=0
          DO N=1,N0
            D=ABS(COORD(3,NB1)-VALZ(N))
            IF(D.LT.0.0001) N1=N1+1
          ENDDO
          IF(N1.EQ.0) THEN
            N0=N0+1
            VALZ(N0)=COORD(3,NB1)
          ENDIF
        ENDDO
        NANEW1=NADS1+NADS2
        NATYP(NAT+1)=NADS1
        IF(NANEW1.EQ.NADS1) GOTO 99
        DO JADS=NADS1+1,NANEW1
          NB1=NB1+1
          DO I=1,3
            COORD(I,NB1)=ADS(I,JADS)
          ENDDO
          N1=0
          DO N=1,N0
            D=ABS(COORD(3,NB1)-VALZ(N))
            IF(D.LT.0.0001) N1=N1+1
          ENDDO
          IF(N1.EQ.0) THEN
            N0=N0+1
            VALZ(N0)=COORD(3,NB1)
          ENDIF
        ENDDO
        NATYP(NAT+2)=NADS2
        NANEW2=NANEW1+NADS3
        IF(NANEW2.EQ.NANEW1) GOTO 99
        DO JADS=NANEW1+1,NANEW2
          NB1=NB1+1
          DO I=1,3
            COORD(I,NB1)=ADS(I,JADS)
          ENDDO
          N1=0
          DO N=1,N0
            D=ABS(COORD(3,NB1)-VALZ(N))
            IF(D.LT.0.0001) N1=N1+1
          ENDDO
          IF(N1.EQ.0) THEN
            N0=N0+1
            VALZ(N0)=COORD(3,NB1)
          ENDIF
        ENDDO
        NATYP(NAT+3)=NADS3
  99    CONTINUE
        NCOUCH=N0-NBZ
        NBZ=N0
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE BASE
C
C  This routine generates the lattice basis vectors for a given Bravais
C           lattice NCRIST centered according to NCENTR
C
      INCLUDE 'spec.inc'
C
      CHARACTER*15  BRAV(8),CENT(7)
      CHARACTER*31 RESEAU
C
      COMMON /BASES/ ATBAS(3*NATP_M),VECBAS(9)
      COMMON /CRANGL/ ALPHAD,BETAD,GAMMAD
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /RESEAU/ NCRIST,NCENTR,IBAS,NAT,A,BSURA,CSURA,UNIT
      COMMON /VECSYS/ ASYS(3),BSYS(3),CSYS(3)

C
      DIMENSION CUB(9),MNC(9),TCN(9),TRG(9),HEX(9)
C
      CHARACTER*3 UNIT
C
      DATA CUB /1.,0.,0.,   0.,1.,0.,   0.,0.,1./
      DATA MNC /1.,0.,1.,   0.,1.,0.,   0.,0.,1./
      DATA TCN /1.,0.,1.,   1.,1.,1.,   0.,0.,1./
      DATA TRG /0.,1.,1.,   -0.866025,-0.5,1.,   0.866025,-0.5,1./
      DATA HEX /1.,0.,0.,   -0.5,0.866025,0.,  0.,0.,1./
      DATA PIS180 /0.017453/
      DATA BRAV /'        CUBIQUE','     TETRAGONAL',' ORTHORHOMBIQUE',
     1           '   MONOCLINIQUE','    TRICLINIQUE','       TRIGONAL',
     2           '      HEXAGONAL','        EXTERNE'/
      DATA CENT /' ','CENTRE',' FACES CENTREES','(RHOMBOEDRIQUE)',
     1     ' FACE A CENTREE',' FACE B CENTREE',' FACE C CENTREE'/
C
      ALPHAR=ALPHAD*PIS180
      BETAR=BETAD*PIS180
      GAMMAR=GAMMAD*PIS180
      NAT3=NAT*3
      GO TO (1,1,1,2,3,4,5,6) NCRIST
C
   1  DO I=1,9
        VECBAS(I)=CUB(I)
      ENDDO
      IF(NCRIST.NE.1) THEN
        VECBAS(9)=CSURA
        IF(NCRIST.EQ.3) THEN
          VECBAS(5)=BSURA
        ENDIF
      ENDIF
      GO TO 6
C
   2  DO I=1,9
        VECBAS(I)=MNC(I)
      ENDDO
      VECBAS(1)=SIN(BETAR)
      VECBAS(3)=COS(BETAR)
      VECBAS(5)=BSURA
      VECBAS(9)=CSURA
      GO TO 6
C
   3  DO I=1,9
        VECBAS(I)=TCN(I)
      ENDDO
      VECBAS(1)=SIN(BETAR)
      VECBAS(3)=COS(BETAR)
      A2Y=(COS(GAMMAR)-COS(ALPHAR)*COS(BETAR))/SIN(BETAR)
      VECBAS(4)=BSURA*A2Y
      VECBAS(5)=BSURA*SQRT(SIN(ALPHAR)*SIN(ALPHAR)-A2Y*A2Y)
      VECBAS(6)=BSURA*COS(ALPHAR)
      VECBAS(9)=CSURA
      GO TO 6
C
   4  IF(((NCENTR.EQ.4).AND.(CSURA.NE.1.)).OR.(NCENTR.EQ.1)) GO TO 5
      ETA=-2.*SIN(ALPHAR/2.)/SQRT(3.)
      DZETA=SQRT(1.-ETA*ETA)
      DO I=1,3
        J=I+2*(I-1)
        J1=J+1
        J2=J+2
        VECBAS(J)=TRG(J)*ETA
        VECBAS(J1)=TRG(J1)*ETA
        VECBAS(J2)=TRG(J2)*DZETA
      ENDDO
      GO TO 6
C
   5  DO I=1,9
        VECBAS(I)=HEX(I)
      ENDDO
      VECBAS(9)=CSURA
C
   6  DO I=1,3
        ASYS(I)=VECBAS(I)
        BSYS(I)=VECBAS(I+3)
        CSYS(I)=VECBAS(I+6)
      ENDDO
      DCA=ABS(CSURA-1.)
      IF((NCRIST.EQ.6).AND.(DCA.LT.0.0001)) GO TO 8
      IF(NCRIST.EQ.8) GO TO 8
      IF(NCENTR.GT.1) THEN
        CALL CENTRE(VECBAS)
        IF(NCENTR.EQ.4) THEN
          DO I=1,9
            VECBAS(I)=VECBAS(I)*SQRT((1.-CSURA*CSURA)*3.)
          ENDDO
          DO I=1,3
            ASYS(I)=VECBAS(I)
            BSYS(I)=VECBAS(I+3)
            CSYS(I)=VECBAS(I+6)
          ENDDO
        ENDIF
      ENDIF
C
   8  RESEAU=BRAV(NCRIST)//' '//CENT(NCENTR)
      WRITE(IUO1,80) RESEAU,NAT
      WRITE(IUO1,81) (VECBAS(I),I=1,9)
      WRITE(IUO1,82)
      WRITE(IUO1,83) (ATBAS(I),I=1,NAT3)
C
  80  FORMAT(////,10X,'RESEAU CRISTALLIN DE TYPE : ',A29,/,16X,
     *       'CONTENANT',I3,' ATOMES DANS LA MAILLE ELEMENTAIRE',//)
  81  FORMAT(28X,'VECTEURS GENERATEURS :',//,26X,'A1 = (',F6.3,',',
     *F6.3,',',F6.3,')',/,26X,'A2 = (',F6.3,',',F6.3,',',F6.3,')',/,
     *26X,'A3 = (',F6.3,',',F6.3,',',F6.3,')')
  82  FORMAT(/,21X,'POSITIONS DES ATOMES DANS LA MAILLE :',/)
  83  FORMAT(29X,'(',F6.3,',',F6.3,',',F6.3,')')
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CENTRE(VECBAS)
C
C  This routine modifies the Bravais lattice basis vectors according to
C                 the way the lattice is centered
C
      COMMON /RESEAU/ NCRIST,NCENTR,IBAS,NAT,A,BSURA,CSURA,UNIT
C
      DIMENSION VECBAS(9),V1(9)
C
      CHARACTER*3 UNIT
C
      DO I=1,9
        V1(I)=VECBAS(I)
      ENDDO
      N1=NCENTR-1
      GO TO (2,3,4,5,6,7) N1
C
   2  DO I=1,3
        VECBAS(I)=-0.5*V1(I)+0.5*V1(I+3)+0.5*V1(I+6)
        VECBAS(I+3)=0.5*V1(I)-0.5*V1(I+3)+0.5*V1(I+6)
        VECBAS(I+6)=0.5*V1(I)+0.5*V1(I+3)-0.5*V1(I+6)
      ENDDO
      GO TO 8
C
   3  DO I=1,3
        VECBAS(I)=0.5*(V1(I+3)+V1(I+6))
        VECBAS(I+3)=0.5*(V1(I)+V1(I+6))
        VECBAS(I+6)=0.5*(V1(I)+V1(I+3))
      ENDDO
      GO TO 8
C
 4    DO I=1,3
        VECBAS(I)=(2./3.)*V1(I)+(1./3.)*V1(I+3)+(1./3.)*V1(I+6)
        VECBAS(I+3)=(-1./3.)*V1(I)+(1./3.)*V1(I+3)+(1./3.)*V1(I+6)
        VECBAS(I+6)=(-1./3.)*V1(I)-(2./3.)*V1(I+3)+(1./3.)*V1(I+6)
      ENDDO
      DO I=1,3
        VECBAS(3*I)=VECBAS(3*I)*SQRT(3./(1.-CSURA*CSURA))
      ENDDO
      GO TO 8
C
   5  DO I=1,3
        VECBAS(I+6)=0.5*(V1(I+3)+V1(I+6))
      ENDDO
      GO TO 8
C
   6  DO I=1,3
        VECBAS(I+6)=0.5*(V1(I)+V1(I+6))
      ENDDO
      GO TO 8
C
   7  DO I=1,3
        VECBAS(I+3)=0.5*(V1(I)+V1(I+3))
      ENDDO
C
   8  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CHBASE(NAT,ATBAS)
C
      COMMON /VECSYS/ ASYS(3),BSYS(3),CSYS(3)
C
      DIMENSION ATBAS(3*NAT),BASVEC(3,3),BAS1(1,3),BAS2(1,3)
C
      DO J=1,3
        BASVEC(1,J)=ASYS(J)
        BASVEC(2,J)=BSYS(J)
        BASVEC(3,J)=CSYS(J)
      ENDDO
C
      DO JAT=1,NAT
        DO J=1,3
          K=J+3*(JAT-1)
          BAS1(1,J)=ATBAS(K)
        ENDDO
        CALL MULMAT(BAS1,1,3,BASVEC,3,3,BAS2)
        DO J=1,3
          K=J+3*(JAT-1)
          ATBAS(K)=BAS2(1,J)
        ENDDO
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE CHNOT(NVEC,VEC1,VEC2)
C
C  This routine linearizes the storage of a two index array
C
      DIMENSION VEC1(3*NVEC),VEC2(3,NVEC)
C
      DO J=1,NVEC
        DO I=1,3
          VEC2(I,J)=VEC1(I+3*(J-1))
        ENDDO
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE INVMAT(B,BINV)
C
      DIMENSION B(3,3),BINV(3,3)
C
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
      A1=B(1,1)*B(2,2)*B(3,3)
      A2=B(2,1)*B(3,2)*B(1,3)
      A3=B(3,1)*B(1,2)*B(2,3)
      A4=B(1,1)*B(3,2)*B(2,3)
      A5=B(2,1)*B(1,2)*B(3,3)
      A6=B(3,1)*B(2,2)*B(1,3)
      DET=A1+A2+A3-A4-A5-A6
C
      IF(ABS(DET).LT.0.0001) GO TO 10
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            L=(I-J)*(I-K)*(J-K)
            IF(L.NE.0) THEN
              XNUM1=B(J,J)*B(K,K)-B(J,K)*B(K,J)
              XNUM2=B(I,K)*B(K,J)-B(K,K)*B(I,J)
              BINV(I,I)=XNUM1/DET
              BINV(I,J)=XNUM2/DET
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      GO TO 50
C
  10  WRITE(IUO1,60)
C
  60  FORMAT(5X,'NON INVERTIBLE MATRIX')
C
  50  CONTINUE
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE MULMAT(A1,IL1,IC1,A2,IL2,IC2,A3)
C
C  This routine performs the matrix multiplication of A1(IL1,IC1) by
C            A2(IL2,IC2) with the result stored in A3(IL1,IC2)
C
      DIMENSION A1(IL1,IC1),A2(IL2,IC2),A3(IL1,IC2)
C
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
      IF(IC1.NE.IL2) THEN
        WRITE(IUO1,10)
      ELSE
        DO I=1,IL1
          DO J=1,IC2
            A3(I,J)=0.
            DO K=1,IC1
              A3(I,J)=A3(I,J)+A1(I,K)*A2(K,J)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
  10  FORMAT(5X,'THESE MATRICES CANNOT BE MULTIPLIED')
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE NUMAT(NUM,NIVA,IL,IM,IN)
C
      DIMENSION I(100)
C
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
      L=2*NIVA+1
      IF(L.GT.100) THEN
        WRITE(IUO1,5)
        STOP
      ENDIF
      L1=NIVA+1
C
      DO K=1,L
        IF(K.LE.L1) THEN
          I(K)=K-1
        ELSE
          I(K)=L1-K
        ENDIF
      ENDDO
C
      Q1=FLOAT(NUM)/FLOAT(L*L)
      JR1=NUM-L*L*INT(Q1+0.0001)
      JS1=INT(Q1+0.9999)
      Q2=FLOAT(JR1)/FLOAT(L)
      JS2=INT(Q2+0.9999)
      IF(JR1.EQ.0) JS2=L
      Q3=FLOAT(NUM)/FLOAT(L)
      JR3=INT(Q3+0.0001)
      JS3=NUM-L*JR3
      IF(JS3.EQ.0) JS3=L
      IL=I(JS1)
      IM=I(JS2)
      IN=I(JS3)
C
   5  FORMAT(///,'<<<<<<<<<<  INCREASE THE SIZE OF I IN',
     1       '  THE NUMAT SUBROUTINE  >>>>>>>>>>')
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE RELA(NINI,NFIN,NAT,VALINI,VALFIN,VALFIN2,
     1                 COORD,NTYP,REL,L)
C
      INCLUDE 'spec.inc'
C
      COMMON /ADSORB/ I1,NATA,N1,N2,N3,ADS(3,900),NCOUCH
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /RELADS/ NRELA,PCRELA(3)
      COMMON /RELAX/ IREL,NREL,PCREL(10),OMEGA1,OMEGA2
C
      DIMENSION VALINI(NATCLU_M),VALFIN(NATCLU_M),REL(NATCLU_M)
      DIMENSION NTYP(NATM),COORD(3,NATCLU_M),LSP(2),DZA(2),DZB(2)
      DIMENSION DYA(2),DYB(2),VALFIN2(NATCLU_M),KZ(1000)
C
      DATA SMALL /0.0001/
C
      IF((IREL.EQ.1).OR.((IREL.EQ.0).AND.(NRELA.GT.0))) THEN
C
        CALL ORDRE(NINI,VALINI,NFIN,VALFIN)
        WRITE(IUO1,70) NFIN
        DO JPLAN=1,NFIN
          IF(JPLAN.LE.NRELA) THEN
            X1=1.
            X2=0.
            PCADS=PCRELA(JPLAN)
          ELSEIF((JPLAN.GT.NRELA).AND.(JPLAN.LE.L)) THEN
            X1=0.
            X2=0.
          ELSE
            X1=0.
            X2=1.
            PCSUBS=PCREL(JPLAN-L)
          ENDIF
          REL(JPLAN)=0.
          IF(JPLAN.GT.NREL+L) GO TO 20
          IF(JPLAN.EQ.NFIN) GO TO 20
          DPLAN=VALFIN(JPLAN)-VALFIN(JPLAN+1)
          REL(JPLAN)=DPLAN*(X1*PCADS+X2*PCSUBS)/100.
  20      DREL=VALFIN(JPLAN)+REL(JPLAN)
          WRITE(IUO1,30) JPLAN,VALFIN(JPLAN),DREL
        ENDDO
C
        NBR=0
        DO JTYP=1,NAT
          NBAT=NTYP(JTYP)
          DO NUM=1,NBAT
            NBR=NBR+1
            DO JPLAN=1,NFIN
              DIF=ABS(COORD(3,NBR)-VALFIN(JPLAN))
              IF(DIF.LT.SMALL) THEN
                COORD(3,NBR)=COORD(3,NBR)+REL(JPLAN)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
        DO JPLAN=1,NFIN
          VALFIN(JPLAN)=VALFIN(JPLAN)+REL(JPLAN)
        ENDDO
C
      ELSEIF(IREL.GE.2) THEN
C
        IP=0
        LSP(2)=0
        OMEGA=OMEGA1
  97    XN1=1.
        XN2=0.
        IP=IP+1
        CALL ORDRE(NINI,VALINI,NFIN,VALFIN)
        ZP=VALFIN(IP)
        CALL RZB110(OMEGA,DY1,DY2,DZ1,DZ2)
        DZA(IP)=DZ1
        DZB(IP)=DZ2
        DYA(IP)=DY1
        DYB(IP)=DY2
        IF(ABS(OMEGA).LT.SMALL) THEN
          LSP(IP)=1
        ELSE
          LSP(IP)=2
        ENDIF
        IF(LSP(IP).EQ.1) GOTO 95
        NBR=0
C
        DO JTYP=1,NAT-NATA
          NBAT=NTYP(JTYP)
          XN1=XN1+1.-FLOAT(JTYP)
          XN2=XN2-1.+FLOAT(JTYP)
          DO JNUM=1,NBAT
            NBR=NBR+1
            ZAT=COORD(3,NBR)-ZP
            IF(ABS(ZAT).LT.SMALL) THEN
              YAT=COORD(2,NBR)
              COORD(2,NBR)=YAT-XN1*DYA(IP)-XN2*DYB(IP)
              COORD(3,NBR)=ZAT+ZP+XN1*DZA(IP)+XN2*DZB(IP)
            ENDIF
          ENDDO
        ENDDO
C
  95    OMEGA=OMEGA2
        IF((IREL.EQ.3).AND.(IP.EQ.1)) GOTO 97
        LS=0
        DO I=1,IP
          LS=LS+LSP(I)
        ENDDO
        NBZ1=NFIN+LS-IP
        DO K=1,IP
          IF(LSP(K).EQ.2) THEN
            IF((K.EQ.2).AND.(LS.EQ.3)) THEN
              KN=K-1
            ELSE
              KN=K
            ENDIF
            VALINI(NBZ1-KN+1)=VALFIN(L+K)+DZB(K)
            REL(NBZ1-KN+1)=DZB(K)
          ELSE
            VALINI(NBZ1-K+1)=VALFIN(L+K)
            REL(NBZ1-K+1)=0.
          ENDIF
        ENDDO
C
        IL=0
        IR=0
        DO J=1,NFIN
          IS=0
          IF(J.LE.NRELA) THEN
            X1=1.
            X2=0.
            X3=0.
            PCADS=PCRELA(J)
            IS=1
          ELSEIF((J.GT.NRELA).AND.(J.LE.L)) THEN
            X1=0.
            X2=0.
            X3=0.
          ELSEIF((J.GT.L).AND.(J.LE.(L+IP))) THEN
            IR=IR+1
            IF(LSP(IR).EQ.1) THEN
              IF((IR.EQ.1).AND.(LSP(2).EQ.2)) GOTO 31
              X1=0.
              X2=1.
              X3=0.
              LT=MAX0(LSP(1),LSP(2))-1
              PCSUBS=PCREL(J-L-LT)
              IL=1
              IS=1
  31          CONTINUE
            ELSE
              X1=0.
              X2=0.
              X3=1.
            ENDIF
          ELSEIF((J.GT.(L+IP)).AND.(J.LE.(L+IP+NREL))) THEN
            X1=0.
            X2=1.
            X3=0.
            LT=MAX0(LSP(1),LSP(2))+IP-1
            PCSUBS=PCREL(J-L-LT+IL+1)
            IS=1
          ELSE
            X1=0.
            X2=0.
            X3=0.
          ENDIF
          DPLAN=VALFIN(J)-VALFIN(J+1)
          REL(J)=X3*DZA(IR)+DPLAN*(X1*PCADS+X2*PCSUBS)/100.
          VALINI(J)=VALFIN(J)+REL(J)
          IF(IS.EQ.1) THEN
            NBR=0
            DO JTYP=1,NAT
              NBAT=NTYP(JTYP)
              DO NUM=1,NBAT
                NBR=NBR+1
                DIF=ABS(COORD(3,NBR)-VALFIN(J))
                IF(DIF.LT.SMALL) THEN
                  COORD(3,NBR)=VALINI(J)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
C
        CALL ORDRE(NBZ1,VALINI,NFIN,VALFIN2)
        WRITE(IUO1,65) NFIN
        KZ(1)=0
        KZ(2)=LSP(1)
        KZ(3)=MAX0(LSP(1),LSP(2))
        DO KK=4,NFIN
          KZ(KK)=LS
        ENDDO
        DO JPLAN=1,NFIN
          IF(JPLAN.LE.L) THEN
            WRITE(IUO1,55) JPLAN,VALFIN(JPLAN),VALFIN2(JPLAN)
            VALINI(JPLAN)=VALFIN(JPLAN)
          ELSEIF((JPLAN.GT.L).AND.(JPLAN.LE.(L+LS))) THEN
            K=KZ(JPLAN-L) - INT((JPLAN-L)/2)
            IPLAN=JPLAN-K
            WRITE(IUO1,55) JPLAN,VALFIN(IPLAN),VALFIN2(JPLAN)
            VALINI(JPLAN)=VALFIN(IPLAN)
          ELSEIF(JPLAN.GT.(L+LS)) THEN
            IPLAN=JPLAN-LS+IP
            WRITE(IUO1,55) JPLAN,VALFIN(IPLAN),VALFIN2(JPLAN)
            VALINI(JPLAN)=VALFIN(IPLAN)
          ENDIF
        ENDDO
      ENDIF
C
  30  FORMAT(/,26X,'THE Z POSITION OF PLANE ',I3,' IS : ',F6.3,
     *      ' BEFORE RELAXATION AND : ',F6.3,' AFTER')
  55  FORMAT(/,26X,'THE Z  POSITION OF PLANE ',I3,' IS : ',F6.3,
     *        ' BEFORE RELAXATION AND : ',F6.3,' AFTER')
  65  FORMAT(//,44X,'THE SUMMATION IS PERFORMED OVER ',I2,' PLANES : ')
  70  FORMAT(//,44X,'THE SUMMATION IS PERFORMED OVER ',I2,' PLANES : ')
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE ROTBAS(ROT)
C
C  This routine calculates the basis vectors related to a surface
C         characterized by its Miller indices (IH,IK,II,IL)
C
      COMMON /MILLER/ IH,IK,II,IL,IVG0,IVN(3)
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /RESEAU/ NCRIST,NCENTR,IBAS,NAT,A,BSURA,CSURA,UNIT
      COMMON /VECSYS/ A1(3),A2(3),A3(3)
C
      DIMENSION ROT(3,3),VECT(3,3),A1STAR(3),A2STAR(3),A3STAR(3),B1(3)
      DIMENSION VECT1(3),XNORM(3),CHBASE(3,3),VECT2(3,3)
C
      CHARACTER*3 UNIT
C
      DATA PI /3.141593/
C
      IF((NCRIST.EQ.8).AND.(IVG0.GE.1)) GOTO 7
      XH=FLOAT(IH)
      XK=FLOAT(IK)
      XI=FLOAT(II)
      XL=FLOAT(IL)
      XI1=-XH-XK
      II1=INT(XI1)
      IF((NCRIST.EQ.7).AND.(XI.NE.XI1)) WRITE(IUO1,5) IH,IK,II1,IL
   5  FORMAT(5X,'THE SURFACE INDICES ARE NOT CORRECT,',/,5X,
     1       'FOR THE REST OF THE CALCULATION, THEY ARE TAKEN AS ',
     2       '(',I2,1X,I2,1X,I2,1X,I2,')')
      CPR=1.
      CALL PRVECT(A2,A3,B1,CPR)
      OMEGA=PRSCAL(A1,B1)/(2.*PI)
      CALL PRVECT(A2,A3,A1STAR,OMEGA)
      CALL PRVECT(A3,A1,A2STAR,OMEGA)
      CALL PRVECT(A1,A2,A3STAR,OMEGA)
      DO 10 I=1,3
        VECT1(I)=XH*A1STAR(I)+XK*A2STAR(I)+XL*A3STAR(I)
  10  CONTINUE
      DO 15 I=1,3
        ROT(I,3)=VECT1(I)/SQRT(PRSCAL(VECT1,VECT1))
  15  CONTINUE
      DO 20 I=1,3
        CHBASE(I,1)=A1(I)
        CHBASE(I,2)=A2(I)
        CHBASE(I,3)=A3(I)
        DO 25 J=1,3
          VECT(I,J)=0.
  25    CONTINUE
  20  CONTINUE
      XHKL=XH*XK*XL
      XHK=XH*XK
      XHL=XH*XL
      XKL=XK*XL
      IF(XHKL.NE.0.) THEN
        VECT(1,1)=-1./XH
        VECT(2,1)=1./XK
        VECT(1,2)=-1./XH
        VECT(3,2)=1./XL
        VECT(2,3)=-1./XK
        VECT(3,3)=1./XL
      ELSEIF(XHK.NE.0.) THEN
        VECT(1,1)=-1./XH
        VECT(2,1)=1./XK
      ELSEIF(XHL.NE.0.) THEN
        VECT(1,2)=-1./XH
        VECT(3,2)=1./XL
      ELSEIF(XKL.NE.0.) THEN
        VECT(2,3)=-1./XK
        VECT(3,3)=1./XL
      ELSEIF(XH.NE.0.) THEN
        VECT(2,2)=1./XH
      ELSEIF(XK.NE.0.) THEN
        VECT(3,3)=1./XK
      ELSEIF(XL.NE.0.) THEN
        VECT(1,1)=1./XL
      ENDIF
      CALL MULMAT(CHBASE,3,3,VECT,3,3,VECT2)
      DO 35 I=1,3
        XNORM(I)=SQRT(VECT2(1,I)**2+VECT2(2,I)**2+VECT2(3,I)**2)
  35  CONTINUE
      XMIN=AMIN1(XNORM(1),XNORM(2),XNORM(3))
      XMAX=AMAX1(XNORM(1),XNORM(2),XNORM(3))
      DO 40 I=1,3
        IF(XHKL.NE.0.) THEN
          IF(ABS(XMIN-XNORM(I)).LT.0.0001) THEN
            DO 45 J=1,3
              ROT(J,1)=VECT2(J,I)/XNORM(I)
  45        CONTINUE
          ENDIF
        ELSE
          IF(ABS(XMAX-XNORM(I)).LT.0.0001) THEN
            DO 50 J=1,3
              ROT(J,1)=VECT2(J,I)/XNORM(I)
  50        CONTINUE
          ENDIF
        ENDIF
  40  CONTINUE
      ROT(1,2)=ROT(2,3)*ROT(3,1)-ROT(3,3)*ROT(2,1)
      ROT(2,2)=ROT(3,3)*ROT(1,1)-ROT(3,1)*ROT(1,3)
      ROT(3,2)=ROT(1,3)*ROT(2,1)-ROT(2,3)*ROT(1,1)
      IF(NCRIST.EQ.7) THEN
        WRITE(IUO1,85) IH,IK,II1,IL
      ELSE
        WRITE(IUO1,80) IH,IK,IL
      ENDIF
      WRITE(IUO1,65) ROT(1,1),ROT(2,1),ROT(3,1)
      WRITE(IUO1,70) ROT(1,2),ROT(2,2),ROT(3,2)
      WRITE(IUO1,75) ROT(1,3),ROT(2,3),ROT(3,3)
      GOTO 37
   7  DO 17 I=1,3
        DO 27 J=1,3
          ROT(I,J)=0.
          IF(I.EQ.J) ROT(I,J)=1.
  27    CONTINUE
  17  CONTINUE
      IF(IVG0.EQ.1) WRITE(IUO1,48)
      IF(IVG0.EQ.2) WRITE(IUO1,47)
  47  FORMAT(//,25X,'LINEAR CHAIN STUDY ')
  48  FORMAT(//,35X,'PLANE STUDY')
  65  FORMAT(26X,'ISURF = (',F6.3,',',F6.3,',',F6.3,')')
  70  FORMAT(26X,'JSURF = (',F6.3,',',F6.3,',',F6.3,')')
  75  FORMAT(26X,'KSURF = (',F6.3,',',F6.3,',',F6.3,')')
  80  FORMAT(//,18X,'BASIS VECTORS FOR THE SURFACE (',I2,1X,I2,1X,
     *I2,') :',/)
  85  FORMAT(//,18X,'BASIS VECTORS FOR THE SURFACE (',I2,1X,I2,1X,
     *I2,1X,I2,') :',/)
C
  37  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE RZB110(OMEGA,DY1,DY2,DZ1,DZ2)
C
      A1=COS(OMEGA)
      ALPHA=SIN(OMEGA)
      BETA=A1-3.
      GAMMA=SQRT(3.)*(5./3.-A1)
      DELTA=SQRT(SQRT(3.)*(1./3.+A1)/GAMMA)
      CSA=SQRT(3.)*(-BETA-ALPHA*DELTA)/6.
      SNA=SQRT(1.-CSA*CSA)
      CSB=-SQRT(3.)*BETA/3. -CSA
      SNB=-SQRT(3.)*ALPHA/3. +SNA
      DY1=(SQRT(3.)*CSB-1.)/4.
      DY2=(1.-SQRT(3.)*CSA)/4.
      DZ1=(SQRT(3.)*SNB-SQRT(2.))/4.
      DZ2=(SQRT(3.)*SNA-SQRT(2.))/4.
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE TEST1(COUPUR,NB1,NB2,ATOME,COORD,VAL,NBAT,IRE,NBZ)
C
      INCLUDE 'spec.inc'
C
      DIMENSION ATOME(3,NATCLU_M),COORD(3,NATCLU_M),VAL(NATCLU_M)
      DIMENSION IRE(NATCLU_M,2)
C
      DIST2=0.
      DO 10 I=1,3
        DIST2=DIST2+ATOME(I,NB1)*ATOME(I,NB1)
  10  CONTINUE
      DIST=SQRT(DIST2)
      V=0.0001
      IF((ATOME(3,NB1).LE.V).AND.(DIST.LE.COUPUR)) THEN
        NBAT=NBAT+1
        NB2=NB2+1
        IRE(NB1,1)=NB2
        IRE(NB1,2)=NBAT
        DO 20 I=1,3
          COORD(I,NB2)=ATOME(I,NB1)
  20    CONTINUE
        IF(NBZ.EQ.0) THEN
          NBZ=NBZ+1
          VAL(NBZ)=COORD(3,NB2)
        ELSE
          N1=0
          DO N=1,NBZ
            D=ABS(COORD(3,NB2)-VAL(N))
            IF(D.LT.0.0001) N1=N1+1
          ENDDO
          IF(N1.EQ.0) THEN
            NBZ=NBZ+1
            VAL(NBZ)=COORD(3,NB2)
          ENDIF
        ENDIF
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE TEST(NIV,ROT,NATYP,NBZ,NAT2,ISURF,COUP,*)
C
      INCLUDE 'spec.inc'
C
      DIMENSION ATOME1(3,NATCLU_M),COORD1(3,NATCLU_M)
      DIMENSION IRE1(NATCLU_M,2),NATYP(NATM)
      DIMENSION NATYP1(NATM),VALZ1(NATCLU_M),ROT(3,3)
C
      NMAX1=(2*NIV+3)**3
      NV1=NIV+1
      CALL AMAS(NV1,ATOME1,COORD1,VALZ1,ISURF,COUP,ROT,IRE1,
     1          NATYP1,NBZ,NAT2,NCOUCH,NMAX1)
      DO 10 I=1,NAT2
        IF(NATYP(I).NE.NATYP1(I)) RETURN 1
  10  CONTINUE
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE ARCSIN(U,CST,RANGLE)
C
C   For a given complex number U, this subroutine calculates its phase
C   Warning : it is implicitely assumed that U = sin(theta) exp(i*phi)
C   with theta > or = to 0 which is always the case when theta is obtained
C   from the coordinates of a given vector r by the ACOS intrinsic function.
C
C   When sin(theta) = 0, then phi = 0 if cos(theta) = 1 and pi if
C   cos(theta) = -1. Cos(theta) is the variable CST.
C
      COMPLEX U,CANGLE
C
      IF(CABS(U).LT.0.0001) THEN
        IF(CST.GT.0.) THEN
          RANGLE=0.
        ELSEIF(CST.LT.0.) THEN
          RANGLE=3.141593
        ENDIF
      ELSE
        CANGLE=(0.,-1.)*CLOG(U/CABS(U))
        RANGLE=REAL(CANGLE)
      ENDIF
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE ATDATA
C
C  This routine contains the atomic mass and the density of all the
C    elements,and the equivalence between their atomic number and
C    chemical symbol.
C
C  Value Z = 0 added for empty spheres. The values entered in this
C    case are arbitrary and set to the corresponding Z = 1 value
C    divided by 1836 (the ratio of the mass of the proton and electron).
C
C                                           Last modified : 25 Apr 2013
C
      REAL XMAT(0:99),RHOAT(0:99),XM_AT,RHO_AT
C
      COMMON /XMRHO/ XM_AT(0:99),RHO_AT(0:99)
C
      DATA XMAT/0.00055,1.00794,4.00260,6.941,9.01218,10.81,12.011,
     1          14.0067,15.9994,18.998403,20.179,22.98977,
     2          24.305,26.98154,28.0855,30.97376,32.06,35.453,
     3          39.948,39.0983,40.08,44.9559,47.88,50.9415,
     4          51.996,54.9380,55.847,58.9332,58.69,63.546,
     5          65.38,69.72,72.59,74.9216,78.96,79.904,83.80,
     6          85.4678,87.62,88.9059,91.22,92.9064,95.94,98.,
     7          101.07,102.9055,106.42,107.8682,112.41,114.82,
     8          118.69,121.75,127.60,126.9045,131.29,132.9054,
     9          137.33,138.9055,140.12,140.9077,144.24,145.,
     *          150.36,151.96,157.25,158.9254,162.50,164.9304,
     *          167.26,168.9342,173.04,174.967,178.49,180.9479,
     *          183.85,186.207,190.2,192.22,195.08,196.9665,
     *          200.59,204.383,207.2,208.9804,209.,210.,222.,
     *          223.,226.0254,227.0278,232.0381,231.0359,
     *          238.0289,237.0482,244.,243.,247.,247.,251.,252./
C
      DATA RHOAT/0.0007,0.0708,0.122,0.533,1.845,2.34,2.26,0.81,1.14,
     1           1.108,1.207,0.969,1.735,2.6941,2.32,1.82,2.07,
     2           1.56,1.40,0.860,1.55,2.980,4.53,6.10,7.18,7.43,
     3           7.860,8.9,8.876,8.94,7.112,5.877,5.307,5.72,
     4           4.78,3.11,2.6,1.529,2.54,4.456,6.494,8.55,10.20,
     5           11.48,12.39,12.39,12.00,10.48,8.63,7.30,7.30,
     6           6.679,6.23,4.92,3.52,1.870,3.5,6.127,6.637,
     7           6.761,6.994,7.20,7.51,5.228,7.8772,8.214,8.525,
     8           8.769,9.039,9.294,6.953,9.811,13.29,16.624,19.3,
     9           20.98,22.53,22.39,21.41,18.85,13.522,11.83,11.33,
     *           9.730,9.30,0.0,4.4,0.0,5.,10.05,11.70,15.34,
     *           18.92,20.21,19.80,13.64,13.49,14.,0.0,0.0/
C
      DO J=0,99
        XM_AT(J)=XMAT(J)
        RHO_AT(J)=RHOAT(J)
      ENDDO
C
      END
C
C=======================================================================
C
      SUBROUTINE AUGER_MULT
C
C  This subroutine computes all the possible multiplets that are
C     contained in a given Auger transition line. It assumes that
C     the atom has closed shells only.
C
C                                        Last modified :  9 March 2006
C
      COMMON /INIT_A/ LI,L2,L1
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
      CHARACTER*1 SC(0:1),LC(0:6),JC(0:7)
      CHARACTER*3 MULTIPLET(112)
C
      DATA SC /'1','3'/
      DATA LC /'S','P','D','F','G','H','I'/
      DATA JC /'0','1','2','3','4','5','6','7'/
C
      WRITE(IUO1,10)
      N_MULT=0
      DO NS=0,1
        DO L=ABS(L1-L2),L1+L2
          DO J=ABS(L-NS),L+NS
            N_MULT=N_MULT+1
            MULTIPLET(N_MULT)=SC(NS)//LC(L)//JC(J)
            WRITE(IUO1,20) MULTIPLET(N_MULT)
          ENDDO
        ENDDO
      ENDDO
C
  10  FORMAT(///,26X,'THE POSSIBLE MULTIPLETS ARE :',/,'  ')
  20  FORMAT(58X,A3)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE BESPHE(NL,IBES,X1,FL)
C
C  This routine computes the spherical Bessel functions for
C                       a real argument X1.
C
C       IBES=1 : Bessel function
C       IBES=2 : Neumann function
C       IBES=3 : Hankel function of the first kind
C       IBES=4 : Hankel function of the second kind
C       IBES=5 : Modified Bessel function
C       IBES=6 : Modified Neumann function
C       IBES=7 : Modified Hankel function
C
C                                         Last modified : 8 Nov 2006
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'spec.inc'
C
      COMPLEX*16 FL(0:2*NL_M),FLNN(0:N_BESS),GL(0:N_BESS),UN,I,C1,C2
      COMPLEX*16 ZERO,CNORM
C
      DOUBLE PRECISION SCALN(0:N_BESS)
C
      REAL X1
C
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
      ECH=37.D0
      COMP=1.D37
      COMM=1.D-37
      X=DBLE(X1)
      NX=INT(X1)
      NREC=5*MAX0(NL-1,NX)
      IF(NREC.GT.N_BESS) GOTO 16
      ITEST=0
      ZERO=(0.D0,0.D0)
      UN=(1.D0,0.D0)
      C1=UN
      I=(0.D0,1.D0)
      C2=I
      DEB=1.D0
      IF((IBES.EQ.3).OR.(IBES.EQ.4)) THEN
        IBES1=1
        IF(IBES.EQ.4) C2=-I
      ELSEIF(IBES.EQ.7) THEN
        IBES1=5
        C2=-UN
      ELSE
        IBES1=IBES
      ENDIF
C
C   Case where the argument is zero
C
      IF(DABS(X).LT.0.000001D0) THEN
        IF((IBES.EQ.1).OR.(IBES.EQ.5)) THEN
          FL(0)=UN
          DO 10 L=1,NL-1
            FL(L)=ZERO
  10      CONTINUE
          ITEST=1
        ELSE
          ITEST=-1
        ENDIF
      ENDIF
      IF(ITEST) 11,12,13
  11  WRITE(IUO1,14)
      STOP
  16  WRITE(IUO1,17) NREC
      STOP
  15  IBES1=IBES1+1
C
C   Initial values
C
  12  A=-1.D0
      B=1.D0
      IF(IBES1.EQ.1) THEN
        FL(0)=UN*DSIN(X)/X
        FLNN(NREC)=ZERO
        SCALN(NREC)=0.D0
        FLNN(NREC-1)=UN*DEB
        SCALN(NREC-1)=0.D0
      ELSEIF(IBES1.EQ.2) THEN
        GL(0)=-UN*DCOS(X)/X
        GL(1)=GL(0)/X -DSIN(X)/X
      ELSEIF(IBES1.EQ.5) THEN
        A=1.D0
        B=-1.D0
        FL(0)=UN*DSINH(X)/X
        FLNN(NREC)=ZERO
        SCALN(NREC)=0.D0
        FLNN(NREC-1)=UN*DEB
        SCALN(NREC-1)=0.D0
      ELSEIF(IBES1.EQ.6) THEN
        A=1.D0
        B=-1.D0
        GL(0)=UN*DCOSH(X)/X
        GL(1)=(DSINH(X)-GL(0))/X
      ENDIF
C
C   Downward reccurence for the spherical Bessel function
C
      IF((IBES1.EQ.1).OR.(IBES1.EQ.5)) THEN
        DO 30 L=NREC-1,1,-1
          ECHEL=0.D0
          SCALN(L-1)=SCALN(L)
          REN=DEXP(SCALN(L)-SCALN(L+1))
          FLNN(L-1)=A*(REN*FLNN(L+1)-B*DFLOAT(2*L+1)*FLNN(L)/X)
          IF(CDABS(FLNN(L-1)).GT.COMP) THEN
            ECHEL=-ECH
          ELSEIF(CDABS(FLNN(L-1)).LT.COMM) THEN
            ECHEL=ECH
          ENDIF
          IF(ECHEL.NE.0.D0 ) SCALN(L-1)=ECHEL+SCALN(L-1)
          FLNN(L-1)=FLNN(L-1)*DEXP(ECHEL)
  30    CONTINUE
        CNORM=FL(0)/FLNN(0)
        DO 40 L=1,NL-1
          FL(L)=CNORM*FLNN(L)*DEXP(SCALN(0)-SCALN(L))
  40    CONTINUE
      ELSE
C
C   Upward recurrence for the spherical Neumann function
C
        DO 20 L=1,NL-1
          IF(IBES.EQ.7) C1=(-UN)**(L+2)
          GL(L+1)=A*GL(L-1)+B*DFLOAT(2*L+1)*GL(L)/X
          IF(IBES1.NE.IBES) THEN
C
C   Calculation of the spherical Hankel function
C
            FL(L+1)=C1*(FL(L+1)+C2*GL(L+1))
          ELSE
            FL(L+1)=GL(L+1)
          ENDIF
  20    CONTINUE
        IF(IBES1.EQ.IBES) THEN
          FL(0)=GL(0)
          FL(1)=GL(1)
        ELSE
          FL(0)=C1*(FL(0)+C2*GL(0))
          FL(1)=C1*(FL(1)+C2*GL(1))
        ENDIF
        IBES1=IBES
      ENDIF
      IF(IBES.NE.IBES1) GOTO 15
C
  13  RETURN
C
  14  FORMAT(/////,3X,'<<<<<<<<<< THE ARGUMENT OF THE BESSEL ',
     1       'FUNCTIONS IS NUL >>>>>>>>>>')
  17  FORMAT(/////,3X,'<<<<<<<<<< THE DIMENSIONNING N_BESS ',
     1       'IS NOT CORRECT FOR SUBROUTINE BESPHE >>>>>>>>>>',//,
     2       15X,'<<<<<<<<<< IT SHOULD BE AT LEAST : ',I5,
     3       ' >>>>>>>>>>')
C
      END
C
C=======================================================================
C
      SUBROUTINE BESPHE2(NL,IBES,X,FL)
C
C  This routine computes the spherical Bessel functions for
C                       a real argument X1.
C
C       IBES=1 : Bessel function
C       IBES=2 : Neumann function
C       IBES=3 : Hankel function of the first kind
C       IBES=4 : Hankel function of the second kind
C       IBES=5 : Modified Bessel function
C       IBES=6 : Modified Neumann function
C       IBES=7 : Modified Hankel function
C
C                                         Last modified : 8 Nov 2006
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'spec.inc'
C
      COMPLEX*16 FL(0:2*NL_M),FLNN(0:N_BESS),GL(0:N_BESS),UN,I,C1,C2
      COMPLEX*16 ZERO,CNORM
C
      DOUBLE PRECISION SCALN(0:N_BESS)
C
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
      ECH=37.D0
      COMP=1.D37
      COMM=1.D-37
      NX=INT(X)
      NREC=5*MAX0(NL-1,NX)
      IF(NREC.GT.N_BESS) GOTO 16
      ITEST=0
      ZERO=(0.D0,0.D0)
      UN=(1.D0,0.D0)
      C1=UN
      I=(0.D0,1.D0)
      C2=I
      DEB=1.D0
      IF((IBES.EQ.3).OR.(IBES.EQ.4)) THEN
        IBES1=1
        IF(IBES.EQ.4) C2=-I
      ELSEIF(IBES.EQ.7) THEN
        IBES1=5
        C2=-UN
      ELSE
        IBES1=IBES
      ENDIF
C
C   Case where the argument is zero
C
      IF(DABS(X).LT.0.000001D0) THEN
        IF((IBES.EQ.1).OR.(IBES.EQ.5)) THEN
          FL(0)=UN
          DO 10 L=1,NL-1
            FL(L)=ZERO
  10      CONTINUE
          ITEST=1
        ELSE
          ITEST=-1
        ENDIF
      ENDIF
      IF(ITEST) 11,12,13
  11  WRITE(IUO1,14)
      STOP
  16  WRITE(IUO1,17) NREC
      STOP
  15  IBES1=IBES1+1
C
C   Initial values
C
  12  A=-1.D0
      B=1.D0
      IF(IBES1.EQ.1) THEN
        FL(0)=UN*DSIN(X)/X
        FLNN(NREC)=ZERO
        SCALN(NREC)=0.D0
        FLNN(NREC-1)=UN*DEB
        SCALN(NREC-1)=0.D0
      ELSEIF(IBES1.EQ.2) THEN
        GL(0)=-UN*DCOS(X)/X
        GL(1)=GL(0)/X -DSIN(X)/X
      ELSEIF(IBES1.EQ.5) THEN
        A=1.D0
        B=-1.D0
        FL(0)=UN*DSINH(X)/X
        FLNN(NREC)=ZERO
        SCALN(NREC)=0.D0
        FLNN(NREC-1)=UN*DEB
        SCALN(NREC-1)=0.D0
      ELSEIF(IBES1.EQ.6) THEN
        A=1.D0
        B=-1.D0
        GL(0)=UN*DCOSH(X)/X
        GL(1)=(DSINH(X)-GL(0))/X
      ENDIF
C
C   Downward reccurence for the spherical Bessel function
C
      IF((IBES1.EQ.1).OR.(IBES1.EQ.5)) THEN
        DO 30 L=NREC-1,1,-1
          ECHEL=0.D0
          SCALN(L-1)=SCALN(L)
          REN=DEXP(SCALN(L)-SCALN(L+1))
          FLNN(L-1)=A*(REN*FLNN(L+1)-B*DFLOAT(2*L+1)*FLNN(L)/X)
          IF(CDABS(FLNN(L-1)).GT.COMP) THEN
            ECHEL=-ECH
          ELSEIF(CDABS(FLNN(L-1)).LT.COMM) THEN
            ECHEL=ECH
          ENDIF
          IF(ECHEL.NE.0.D0 ) SCALN(L-1)=ECHEL+SCALN(L-1)
          FLNN(L-1)=FLNN(L-1)*DEXP(ECHEL)
  30    CONTINUE
        CNORM=FL(0)/FLNN(0)
        DO 40 L=1,NL-1
          FL(L)=CNORM*FLNN(L)*DEXP(SCALN(0)-SCALN(L))
  40    CONTINUE
      ELSE
C
C   Upward recurrence for the spherical Neumann function
C
        DO 20 L=1,NL-1
          IF(IBES.EQ.7) C1=(-UN)**(L+2)
          GL(L+1)=A*GL(L-1)+B*DFLOAT(2*L+1)*GL(L)/X
          IF(IBES1.NE.IBES) THEN
C
C   Calculation of the spherical Hankel function
C
            FL(L+1)=C1*(FL(L+1)+C2*GL(L+1))
          ELSE
            FL(L+1)=GL(L+1)
          ENDIF
  20    CONTINUE
        IF(IBES1.EQ.IBES) THEN
          FL(0)=GL(0)
          FL(1)=GL(1)
        ELSE
          FL(0)=C1*(FL(0)+C2*GL(0))
          FL(1)=C1*(FL(1)+C2*GL(1))
        ENDIF
        IBES1=IBES
      ENDIF
      IF(IBES.NE.IBES1) GOTO 15
C
  13  RETURN
C
  14  FORMAT(/////,3X,'<<<<<<<<<< THE ARGUMENT OF THE BESSEL ',
     1       'FUNCTIONS IS NUL >>>>>>>>>>')
  17  FORMAT(/////,3X,'<<<<<<<<<< THE DIMENSIONNING N_BESS ',
     1       'IS NOT CORRECT FOR SUBROUTINE BESPHE >>>>>>>>>>',//,
     2       15X,'<<<<<<<<<< IT SHOULD BE AT LEAST : ',I5,
     3       ' >>>>>>>>>>')
C
      END
C
C=======================================================================
C
      SUBROUTINE CHECK_VIB(NAT2)
C
C  This subroutines checks the geometrical environment of each atom
C     to identify those which can move "freely" in one direction, in
C     order to see whether the mean square displacement in this
C     direction is of bulk type or surface type
C
C  An atom is considered to move freely in one direction if no other
C     atom is present in the tetragonal cell of height ALENGTH * A
C     and base edge 2 * A, whose base is centered on the atom considered
C
C  Only prototypical atoms are considered as all equivalent atoms are
C     in the same geometrical environment
C
C  Surface-like atoms are then identified as having I_FREE = 1
C
C                                     Last modified : 24 Apr 2013
C
      INCLUDE 'spec.inc'
C
      COMMON /COOR/ NATCLU,N_PROT,NATYP(NATM),NCHTYP(NATP_M),
     1              NCORR(NAT_EQ_M,NATP_M),INEW_AT(NATCLU_M),
     2              COORD(3,NATCLU_M)
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /VIBRAT/ I_FREE(NATP_M)
C
      INTEGER NSUR(NATP_M)
C
      DATA SMALL /0.0001/
C
      ALENGTH=4.
C
C.................... Checking the z direction ....................
C
      WRITE(IUO1,11)
      N_SUR=0
C
C  Loop on the prototypical atoms
C
      DO JTYP=1,N_PROT
C
        I_FREE(JTYP)=0
        JAT0=NCORR(1,JTYP)
        XA=COORD(1,JAT0)
        YA=COORD(2,JAT0)
        ZA=COORD(3,JAT0)
C
C  Loop on the surrounding atoms
C
        I_ACC=0
C
        DO JAT=1,NAT2
C
          IF(JAT.EQ.JAT0) GOTO 10
C
          X=COORD(1,JAT)
          Y=COORD(2,JAT)
          Z=COORD(3,JAT)
C
C  Considering only atoms with Z > ZA
C
          IF(Z.LT.(ZA+SMALL)) GOTO 10
C
C  Lateral and vertical distances between the two atoms
C
          D_LAT=(X-XA)*(X-XA)+(Y-YA)*(Y-YA)
          D_VER=(Z-ZA)*(Z-ZA)
C
          IF(D_VER.LT.(ALENGTH+SMALL)) THEN
            IF(D_LAT.LT.(1.+SMALL)) THEN
              I_ACC=I_ACC+1
            ENDIF
          ENDIF
C
          IF(I_ACC.GE.1) GOTO 10
C
  10      CONTINUE
C
        ENDDO
C
        IF(I_ACC.EQ.0) THEN
          I_FREE(JTYP)=1
          N_SUR=N_SUR+1
          NSUR(N_SUR)=JTYP
        ENDIF
C
      ENDDO
C
      WRITE(IUO1,12) (NSUR(J),J=1,N_SUR)
C
  11  FORMAT(//,18X,'SURFACE-LIKE ATOMS FOR MSD CALCULATIONS: ',/)
  12  FORMAT(20X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE DJMN(RBETA,R,LMAX)
C
C  This routine calculates Wigner rotation matrices R^{L}_{M1 M2} up to
C             order LMAX, following Messiah's convention.
C                   They are stored as R(M2,M1,L).
C
C                                         Last modified : 20 Oct 2006
C
      INCLUDE 'spec.inc'
C
      COMMON /COEFRLM/ CF(0:2*NL_M-2,0:2*NL_M-2,0:2*NL_M-2)
      COMMON /EXPROT/ EXPR(0:2*NL_M-2,0:2*NL_M-2)
C
      INTEGER EPS0
C
      DIMENSION R(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
C
      DATA SMALL,SQR2 /0.001,1.4142136/
C
      C=COS(RBETA)*0.5
      S=SIN(RBETA)*0.5
      CC=C+C
      CMUL=-1.
      IF(ABS(S).LT.SMALL) THEN
        IF(C.GT.0.) EPS0=1
        IF(C.LT.0.) EPS0=-1
        DO L=0,LMAX
          DO M1=-L,L
            DO M2=-L,L
              IF(M1.NE.M2*EPS0) THEN
                R(M2,M1,L)=0.
              ELSE
                IF(EPS0.EQ.1) THEN
                  R(M2,M1,L)=1.
                ELSE
                  IF(MOD(L+M1,2).EQ.0) THEN
                    R(M2,M1,L)=1.
                  ELSE
                    R(M2,M1,L)=-1.
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        S1=S*SQR2
        C1=0.5+C
        R(0,0,0)=1.0
        R(-1,-1,1)=C1
        R(0,-1,1)=S1
        R(1,-1,1)=1.-C1
        R(-1,0,1)=-S1
        R(0,0,1)=CC
        R(1,0,1)=S1
        R(-1,1,1)=1.-C1
        R(0,1,1)=-S1
        R(1,1,1)=C1
C
        PRODL=-S
        COEF=-S/C1
        CL=-1.
        DO L=2,LMAX
          CL=-CL
          L1=L-1
          FLL1=CC*FLOAT(L+L1)
          FLL2=1./(FLOAT(L*L1)*CC)
          PRODL=-PRODL*S
C
C  Case M = 0
C
          R_1=EXPR(0,L)*PRODL
          R(-L,0,L)=R_1
C
          R(L,0,L)=R_1*CL
          R(0,-L,L)=R_1*CL
C
          R(0,L,L)=R_1
C
          CM2=CL
          DO M2=-L1,-1
            CM2=CM2*CMUL
            CF1=CF(L1,0,-M2)/FLL1
            CF2=FLL1/CF(L,0,-M2)
            IF(-M2.LT.L1) THEN
              R_A=CF2*(R(M2,0,L1)-R(M2,0,L-2)*CF1)
            ELSE
              R_A=CF2*R(M2,0,L1)
            ENDIF
C
            R(M2,0,L)=R_A
C
            R(-M2,0,L)=R_A*CM2
            R(0,M2,L)=R_A*CM2
C
            R(0,-M2,L)=R_A
C
          ENDDO
C
          R(0,0,L)=FLL1*R(0,0,L1)/CF(L,0,0)-
     1             R(0,0,L-2)*CF(L1,0,0)/CF(L,0,0)
C
C  Case M > 0
C
          PRODM=1.
          CM=CL
          FLLM=0.
          DO M=1,L1
            CM=-CM
            PRODM=PRODM*COEF
            FLLM=FLLM+FLL2
C
            R_1=EXPR(M,L)*PRODL*PRODM
            R_2=R_1/(PRODM*PRODM)
C
            R(-L,M,L)=R_1
            R(-L,-M,L)=R_2
C
            R(L,-M,L)=R_1*CM
            R(M,-L,L)=R_1*CM
            R(L,M,L)=R_2*CM
            R(-M,-L,L)=R_2*CM
C
            R(-M,L,L)=R_1
            R(M,L,L)=R_2
C
            CM2=CM
            DO M2=-L1,-M
              CM2=-CM2
              D0=FLOAT(M2)*FLLM
              CF1=CF(L1,M,-M2)/FLL1
              CF2=FLL1/CF(L,M,-M2)
              IF((M.LT.L1).AND.(-M2.LT.L1)) THEN
                R_A=CF2*((1.-D0)*R(M2,M,L1)-R(M2,M,L-2)*CF1)
                R_B=CF2*((1.+D0)*R(M2,-M,L1)-R(M2,-M,L-2)*CF1)
              ELSE
                R_A=CF2*(1.-D0)*R(M2,M,L1)
                R_B=CF2*(1.+D0)*R(M2,-M,L1)
              ENDIF
C
              R(M2,M,L)=R_A
              R(M2,-M,L)=R_B
C
              R(-M2,-M,L)=R_A*CM2
              R(M,M2,L)=R_A*CM2
              R(-M,M2,L)=R_B*CM2
              R(-M2,M,L)=R_B*CM2
C
              R(-M,-M2,L)=R_A
              R(M,-M2,L)=R_B
C
            ENDDO
          ENDDO
C
          PRODM=PRODM*COEF
          R_1=PRODL*PRODM
          R_2=PRODL/PRODM
          R(-L,L,L)=R_1
          R(L,-L,L)=R_1
          R(L,L,L)=R_2
          R(-L,-L,L)=R_2
C
        ENDDO
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE DJMN2(RBETA,R,LMAX,ISWITCH)
C
C  This routine calculates Wigner rotation matrices R^{L}_{M1 M2} up to
C             order LMAX, following Messiah's convention.
C  They are stored as R(M2,M1,L) and multiplied (ISWITCH=1) or divided
C             by EXPF.
C
C                                         Last modified : 20 Oct 2006
C
      INCLUDE 'spec.inc'
C
      COMMON /COEFRLM/ CF(0:2*NL_M-2,0:2*NL_M-2,0:2*NL_M-2)
      COMMON /EXPFAC/ EXPF(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /EXPROT/ EXPR(0:2*NL_M-2,0:2*NL_M-2)
C
      INTEGER EPS0
C
      DIMENSION R(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
C
      DATA SMALL,SQR2 /0.001,1.4142136/
C
      C=COS(RBETA)*0.5
      S=SIN(RBETA)*0.5
      CC=C+C
      CMUL=-1.
      IF(ABS(S).LT.SMALL) THEN
        IF(C.GT.0.) EPS0=1
        IF(C.LT.0.) EPS0=-1
        DO L=0,LMAX
          DO M1=-L,L
            DO M2=-L,L
              IF(M1.NE.M2*EPS0) THEN
                R(M2,M1,L)=0.
              ELSE
                IF(EPS0.EQ.1) THEN
                  R(M2,M1,L)=1.
                ELSE
                  IF(MOD(L+M1,2).EQ.0) THEN
                    R(M2,M1,L)=1.
                  ELSE
                    R(M2,M1,L)=-1.
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        S1=S*SQR2
        C1=0.5+C
        R(0,0,0)=1.0
        R(-1,-1,1)=C1
        R(0,-1,1)=S1
        R(1,-1,1)=1.-C1
        R(-1,0,1)=-S1
        R(0,0,1)=CC
        R(1,0,1)=S1
        R(-1,1,1)=1.-C1
        R(0,1,1)=-S1
        R(1,1,1)=C1
C
        PRODL=-S
        COEF=-S/C1
        CL=-1.
        DO L=2,LMAX
          CL=-CL
          L1=L-1
          FLL1=CC*FLOAT(L+L1)
          FLL2=1./(FLOAT(L*L1)*CC)
          PRODL=-PRODL*S
C
C  Case M = 0
C
          R_1=EXPR(0,L)*PRODL
          R(-L,0,L)=R_1
C
          R(L,0,L)=R_1*CL
          R(0,-L,L)=R_1*CL
C
          R(0,L,L)=R_1
C
          CM2=CL
          DO M2=-L1,-1
            CM2=CM2*CMUL
            CF1=CF(L1,0,-M2)/FLL1
            CF2=FLL1/CF(L,0,-M2)
            IF(-M2.LT.L1) THEN
              R_A=CF2*(R(M2,0,L1)-R(M2,0,L-2)*CF1)
            ELSE
              R_A=CF2*R(M2,0,L1)
            ENDIF
C
            R(M2,0,L)=R_A
C
            R(-M2,0,L)=R_A*CM2
            R(0,M2,L)=R_A*CM2
C
            R(0,-M2,L)=R_A
C
          ENDDO
C
          R(0,0,L)=FLL1*R(0,0,L1)/CF(L,0,0)-
     1             R(0,0,L-2)*CF(L1,0,0)/CF(L,0,0)
C
C  Case M > 0
C
          PRODM=1.
          CM=CL
          FLLM=0.
          DO M=1,L1
            CM=-CM
            PRODM=PRODM*COEF
            FLLM=FLLM+FLL2
C
            R_1=EXPR(M,L)*PRODL*PRODM
            R_2=R_1/(PRODM*PRODM)
C
            R(-L,M,L)=R_1
            R(-L,-M,L)=R_2
C
            R(L,-M,L)=R_1*CM
            R(M,-L,L)=R_1*CM
            R(L,M,L)=R_2*CM
            R(-M,-L,L)=R_2*CM
C
            R(-M,L,L)=R_1
            R(M,L,L)=R_2
C
            CM2=CM
            DO M2=-L1,-M
              CM2=-CM2
              D0=FLOAT(M2)*FLLM
              CF1=CF(L1,M,-M2)/FLL1
              CF2=FLL1/CF(L,M,-M2)
              IF((M.LT.L1).AND.(-M2.LT.L1)) THEN
                R_A=CF2*((1.-D0)*R(M2,M,L1)-R(M2,M,L-2)*CF1)
                R_B=CF2*((1.+D0)*R(M2,-M,L1)-R(M2,-M,L-2)*CF1)
              ELSE
                R_A=CF2*(1.-D0)*R(M2,M,L1)
                R_B=CF2*(1.+D0)*R(M2,-M,L1)
              ENDIF
C
              R(M2,M,L)=R_A
              R(M2,-M,L)=R_B
C
              R(-M2,-M,L)=R_A*CM2
              R(M,M2,L)=R_A*CM2
              R(-M,M2,L)=R_B*CM2
              R(-M2,M,L)=R_B*CM2
C
              R(-M,-M2,L)=R_A
              R(M,-M2,L)=R_B
C
            ENDDO
          ENDDO
C
          PRODM=PRODM*COEF
          R_1=PRODL*PRODM
          R_2=PRODL/PRODM
          R(-L,L,L)=R_1
          R(L,-L,L)=R_1
          R(L,L,L)=R_2
          R(-L,-L,L)=R_2
C
        ENDDO
      ENDIF
C
      IF(ISWITCH.EQ.1) THEN
        DO L=0,LMAX
          DO M1=-L,L
            DO M2=-L,L
              R(M2,M1,L)=SQRT(FLOAT(L+L+1))*R(M2,M1,L)*EXPF(ABS(M2),L)
            ENDDO
          ENDDO
        ENDDO
      ELSEIF(ISWITCH.EQ.2) THEN
        DO L=0,LMAX
          DO M1=-L,L
            DO M2=-L,L
              R(M2,M1,L)=SQRT(FLOAT(L+L+1))*R(M2,M1,L)/EXPF(ABS(M2),L)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE EMETT(JEM,IEMET,Z,COORD,NATYP,EMET,NM,
     1                 JNEM,*)
C
C  This routine looks for the position of an absorber of type IEMET(JEM)
C    situated in the plane at Z. The result is stored in EMET(3)
C
      INCLUDE 'spec.inc'
C
      DIMENSION IEMET(NEMET_M)
      DIMENSION EMET(3),DIST(NATCLU_M),COORD(3,NATCLU_M),NATYP(NATM)
C
      KEMET=0
      JNT=0
      IEM=IEMET(JEM)
      IF(IEM.GT.1) THEN
        DO JTP=1,IEM-1
          JNT=JNT+NATYP(JTP)
        ENDDO
      ENDIF
      NB=NATYP(IEM)
      XMIN=1000000.
C
      DO J=1,NB
        JN=J+JNT
        DELTAZ=ABS(COORD(3,JN)-Z)
        IF(DELTAZ.LT.0.0001) THEN
          XX=COORD(1,JN)
          XY=COORD(2,JN)
          XZ=COORD(3,JN)
          DIST(J)=SQRT(XX*XX+XY*XY+XZ*XZ)
          IF(DIST(J).LT.XMIN) THEN
            XMIN=DIST(J)
            NM=IEM
            JNEM=J
            DO I=1,3
              EMET(I)=COORD(I,JN)
            ENDDO
          ENDIF
          KEMET=KEMET+1
        ENDIF
      ENDDO
C
      IF(KEMET.EQ.0) THEN
        NM=IEM
        RETURN 1
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE EULER(RTHETA1,RPHI1,RTHETA2,RPHI2,RALPHA,RBETA,RGAMMA,
     1                 IROT)
C
C  This routine calculates the Euler angles RALPHA,RBETA,RGAMMA corresponding
C       to the rotation r1(RTHETA1,RPHI1) ----> r2(RTHETA2,RPHI2)
C
C       IROT=1 : r ---> z represented by (0,RTHETA,PI-RPHI)
C       IROT=0 : r ---> z represented by (0,-RTHETA,-RPHI)
C
C
      COMPLEX U1,U2
C
      DATA PI /3.141593/
C
      IF(IROT.EQ.1) THEN
        EPS=1
      ELSE
        EPS=-1
      ENDIF
      DPHI=RPHI2-RPHI1
      A1=SIN(RTHETA1)*COS(RTHETA2)
      A2=COS(RTHETA1)*SIN(RTHETA2)
      A3=COS(RTHETA1)*COS(RTHETA2)
      A4=SIN(RTHETA1)*SIN(RTHETA2)
      U1=A1-A2*COS(DPHI)-(0.,1.)*SIN(RTHETA2)*SIN(DPHI)
      U2=A1*COS(DPHI)-A2+(0.,1.)*SIN(RTHETA1)*SIN(DPHI)
      U3=A3+A4*COS(DPHI)
      IF(U3.GT.1.) U3=1.
      IF(U3.LT.-1.) U3=-1.
      RBETA=ACOS(U3)
      IF(ABS(SIN(RBETA)).GT.0.0001) THEN
        U1=EPS*U1/SIN(RBETA)
        U2=EPS*U2/SIN(RBETA)
        CALL ARCSIN(U1,U3,RALPHA)
        CALL ARCSIN(U2,U3,RGAMMA)
      ELSE
        RALPHA=0.
        IF(ABS(U3-1.0).LT.0.0001) THEN
          RGAMMA=0.
        ELSE
          RGAMMA=PI
        ENDIF
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE GAUNT(L2,M2,L1,M1,GNT)
C
C   This subroutine calculates the Gaunt coefficient G(L2,L3|L1)
C    using a downward recursion scheme due to Schulten and Gordon
C    for the Wigner's 3j symbols. The result is stored as GNT(L3),
C    making use of the selection rule M3 = M1 - M2.
C
C   Ref. :  K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
C
C                                          Last modified :  8 Dec 2008
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'spec.inc'
C
      REAL GNT(0:N_GAUNT)
C
      DOUBLE PRECISION F(0:N_GAUNT),G(0:N_GAUNT),A(0:N_GAUNT)
      DOUBLE PRECISION A1(0:N_GAUNT),B(0:N_GAUNT)
C
      COMMON /LOGAMAD/ GLD(0:N_GAUNT,2)
C
      DATA PI4/12.566370614359D0/
C
      L12=L1+L2
      K12=L1-L2
C
      DO J=1,N_GAUNT
        GNT(J)=0.
      ENDDO
C
      IF((ABS(M1).GT.L1).OR.(ABS(M2).GT.L2)) GOTO 10
C
      M3=M1-M2
      LM1=L1+M1
      LM2=L2+M2
      KM1=L1-M1
      KM2=L2-M2
C
      IF(MOD(M1,2).EQ.0) THEN
        COEF=DSQRT(DFLOAT((2*L1+1)*(2*L2+1))/PI4)
      ELSE
        COEF=-DSQRT(DFLOAT((2*L1+1)*(2*L2+1))/PI4)
      ENDIF
C
      F(L12+1)=0.D0
      G(L12+1)=0.D0
      A(L12+1)=0.D0
      A1(L12+1)=0.D0
      D1=GLD(2*L2+1,1)-GLD(2*L12+2,1)
      D2=GLD(2*L1+1,1)-GLD(LM2+1,1)
      D3=GLD(L12+M3+1,1)-GLD(KM2+1,1)
      D4=GLD(L12-M3+1,1)-GLD(LM1+1,1)-GLD(KM1+1,1)
C
      IF(MOD(KM1-KM2,2).EQ.0) THEN
        F(L12)=DSQRT(DEXP(D1+D2+D3+D4))
      ELSE
        F(L12)=-DSQRT(DEXP(D1+D2+D3+D4))
      ENDIF
C
      D5=0.5D0*(GLD(2*L1+1,1)+GLD(2*L2+1,1)-GLD(2*L12+2,1))
      D6=GLD(L12+1,1)-GLD(L1+1,1)-GLD(L2+1,1)
C
      IF(MOD(K12,2).EQ.0) THEN
        G(L12)=DEXP(D5+D6)
      ELSE
        G(L12)=-DEXP(D5+D6)
      ENDIF
C
      A(L12)=2.D0*DSQRT(DFLOAT(L1*L2*(1+2*L12)*(L12*L12-M3*M3)))
      B(L12)=-DFLOAT((2*L12+1)*((L2*L2-L1*L1-K12)*M3+L12*(L12+1)
     1        *(M2+M1)))
      A1(L12)=2.D0*DFLOAT(L12)*DSQRT(DFLOAT(L1*L2*(1+2*L12)))
C
      IF(ABS(M3).LE.L12) THEN
        GNT(L12)=SNGL(COEF*F(L12)*G(L12)*DSQRT(DFLOAT(2*L12+1)))
      ELSE
        GNT(L12)=0.
      ENDIF
C
      JMIN=MAX0(ABS(K12),ABS(M3))
C
      DO J=L12-1,JMIN,-1
        J1=J+1
        J2=J+2
        A(J)=DSQRT(DFLOAT((J*J-K12*K12))*DFLOAT((L12+1)*(L12+1)-J*J)
     1             *DFLOAT(J*J-M3*M3))
        B(J)=-DFLOAT((2*J+1)*
     1        (L2*(L2+1)*M3-L1*(L1+1)*M3+J*J1*(M2+M1)))
        A1(J)=DFLOAT(J)*DSQRT(DFLOAT((J*J-K12*K12)*
     1                              ((L12+1)*(L12+1)-J*J)))
        F(J)=-(DFLOAT(J1)*A(J2)*F(J2)+B(J1)*F(J1))/(DFLOAT(J2)*A(J1))
        G(J)=-(DFLOAT(J1)*A1(J2)*G(J2))/(DFLOAT(J2)*A1(J1))
        GND=COEF*F(J)*G(J)*DSQRT(DFLOAT(2*J+1))
C
        IF(ABS(M3).LE.J) THEN
          GNT(J)=SNGL(GND)
        ELSE
          GNT(J)=0.
        ENDIF
C
      ENDDO
C
  10  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE GAUNT2(L2,M2,L1,M1,GNT)
C
C   This subroutine calculates the Gaunt coefficient G(L2,L3|L1)
C    using a downward recursion scheme due to Schulten and Gordon
C    for the Wigner's 3j symbols. The result is stored as GNT(L3),
C    making use of the selection rule M3 = M1 - M2.
C
C   Ref. :  K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
C
C                                          Last modified :  8 Dec 2008
C   This is the double precision version
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'spec.inc'
C
      REAL*8 F(0:N_GAUNT),G(0:N_GAUNT),A(0:N_GAUNT),A1(0:N_GAUNT)
      REAL*8 B(0:N_GAUNT),GNT(0:N_GAUNT)
C
      COMMON /LOGAMAD/ GLD(0:N_GAUNT,2)
C
      DATA PI4/12.566370614359D0/
C
      L12=L1+L2
      K12=L1-L2
C
      DO J=1,N_GAUNT
        GNT(J)=0.D0
      ENDDO
C
      IF((ABS(M1).GT.L1).OR.(ABS(M2).GT.L2)) GOTO 10
C
      M3=M1-M2
      LM1=L1+M1
      LM2=L2+M2
      KM1=L1-M1
      KM2=L2-M2
C
      IF(MOD(M1,2).EQ.0) THEN
        COEF=DSQRT(DFLOAT((2*L1+1)*(2*L2+1))/PI4)
      ELSE
        COEF=-DSQRT(DFLOAT((2*L1+1)*(2*L2+1))/PI4)
      ENDIF
C
      F(L12+1)=0.D0
      G(L12+1)=0.D0
      A(L12+1)=0.D0
      A1(L12+1)=0.D0
      D1=GLD(2*L2+1,1)-GLD(2*L12+2,1)
      D2=GLD(2*L1+1,1)-GLD(LM2+1,1)
      D3=GLD(L12+M3+1,1)-GLD(KM2+1,1)
      D4=GLD(L12-M3+1,1)-GLD(LM1+1,1)-GLD(KM1+1,1)
C
      IF(MOD(KM1-KM2,2).EQ.0) THEN
        F(L12)=DSQRT(DEXP(D1+D2+D3+D4))
      ELSE
        F(L12)=-DSQRT(DEXP(D1+D2+D3+D4))
      ENDIF
C
      D5=0.5D0*(GLD(2*L1+1,1)+GLD(2*L2+1,1)-GLD(2*L12+2,1))
      D6=GLD(L12+1,1)-GLD(L1+1,1)-GLD(L2+1,1)
C
      IF(MOD(K12,2).EQ.0) THEN
        G(L12)=DEXP(D5+D6)
      ELSE
        G(L12)=-DEXP(D5+D6)
      ENDIF
C
      A(L12)=2.D0*DSQRT(DFLOAT(L1*L2*(1+2*L12)*(L12*L12-M3*M3)))
      B(L12)=-DFLOAT((2*L12+1)*((L2*L2-L1*L1-K12)*M3+L12*(L12+1)
     1        *(M2+M1)))
      A1(L12)=2.D0*DFLOAT(L12)*DSQRT(DFLOAT(L1*L2*(1+2*L12)))
C
      IF(ABS(M3).LE.L12) THEN
        GNT(L12)=COEF*F(L12)*G(L12)*DSQRT(DFLOAT(2*L12+1))
      ELSE
        GNT(L12)=0.D0
      ENDIF
C
      JMIN=MAX0(ABS(K12),ABS(M3))
C
      DO J=L12-1,JMIN,-1
        J1=J+1
        J2=J+2
        A(J)=DSQRT(DFLOAT((J*J-K12*K12))*DFLOAT((L12+1)*(L12+1)-J*J)
     1             *DFLOAT(J*J-M3*M3))
        B(J)=-DFLOAT((2*J+1)*
     1        (L2*(L2+1)*M3-L1*(L1+1)*M3+J*J1*(M2+M1)))
        A1(J)=DFLOAT(J)*DSQRT(DFLOAT((J*J-K12*K12)*
     1                              ((L12+1)*(L12+1)-J*J)))
        F(J)=-(DFLOAT(J1)*A(J2)*F(J2)+B(J1)*F(J1))/(DFLOAT(J2)*A(J1))
        G(J)=-(DFLOAT(J1)*A1(J2)*G(J2))/(DFLOAT(J2)*A1(J1))
        GND=COEF*F(J)*G(J)*DSQRT(DFLOAT(2*J+1))
C
        IF(ABS(M3).LE.J) THEN
          GNT(J)=GND
        ELSE
          GNT(J)=0.D0
        ENDIF
C
      ENDDO
C
  10  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE HARSPH(NL,THETA,PHI,YLM,NC)
C
C  This routine computes the complex spherical harmonics using Condon and
C                  Shortley phase convention.
C
      INCLUDE 'spec.inc'
C
      COMMON /EXPFAC2/ EXPF2(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /FACTSQ/ FSQ(0:2*NL_M-2)
C
      COMPLEX YLM(0:NL,-NL:NL),COEF,YMM,YMMP,C
C
      DATA SQ4PI_INV,SQR3_INV /0.282095,0.488602/
      DATA PI,SMALL /3.141593,0.0001/
C
      X=COS(THETA)
      IF(ABS(X).LT.SMALL) X=0.0
      IF(ABS(X+1.).LT.SMALL) X=-1.0
      IF(ABS(X-1.).LT.SMALL) X=1.0
C
      YLM(0,0)=CMPLX(SQ4PI_INV)
      YLM(1,0)=X*SQR3_INV
      DO L=2,NC
        Y=1./FLOAT(L)
        YLM(L,0)=X*SQRT(4.-Y*Y)*YLM(L-1,0) -
     1           (1.-Y)*SQRT(1.+2./(FLOAT(L)-1.5))*YLM(L-2,0)
      ENDDO
C
      C2=-1.
      IF((THETA.GE.0.).AND.(THETA.LE.PI)) THEN
        C=-0.5*SQRT(1.-X*X)*EXP((0.,1.)*PHI)
      ELSE
        C=0.5*SQRT(1.-X*X)*EXP((0.,1.)*PHI)
      ENDIF
C
      C1=1.
      COEF=(1.,0.)
      DO M=1,NC
        C1=C1*C2
        COEF=COEF*C
        YMM=SQ4PI_INV*COEF*FSQ(M)
        YLM(M,M)=YMM
        YLM(M,-M)=C1*CONJG(YMM)
        YMMP=X*SQRT(FLOAT(M+M+3))*YMM
        YLM(M+1,M)=YMMP
        YLM(M+1,-M)=C1*CONJG(YMMP)
        IF(M.LT.NC-1) THEN
          DO L=M+2,NC
            YLM(L,M)=(X*(L+L-1)*EXPF2(L-1,M)*YLM(L-1,M) -
     1                (L+M-1)*EXPF2(L-2,M)*YLM(L-2,M))/(EXPF2(L,M)*
     2                (L-M))
            YLM(L,-M)=C1*CONJG(YLM(L,M))
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE HARSPH2(NL,THETA,PHI,YLM,NC)
C
C  This routine computes the complex spherical harmonics using Condon and
C                  Shortley phase convention. This version for m=0 only
C
      INCLUDE 'spec.inc'
C
      COMMON /EXPFAC2/ EXPF2(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /FACTSQ/ FSQ(0:2*NL_M-2)
C
      COMPLEX YLM(0:NL,-NL:NL),COEF,YMM,YMMP,C
C
      DATA SQ4PI_INV,SQR3_INV /0.282095,0.488602/
      DATA PI,SMALL /3.141593,0.0001/
C
      X=COS(THETA)
      IF(ABS(X).LT.SMALL) X=0.0
      IF(ABS(X+1.).LT.SMALL) X=-1.0
      IF(ABS(X-1.).LT.SMALL) X=1.0
C
      YLM(0,0)=CMPLX(SQ4PI_INV)
      YLM(1,0)=X*SQR3_INV
      DO L=2,NC
        Y=1./FLOAT(L)
        YLM(L,0)=X*SQRT(4.-Y*Y)*YLM(L-1,0) -
     1           (1.-Y)*SQRT(1.+2./(FLOAT(L)-1.5))*YLM(L-2,0)
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE HARSPH3(NL,THETA,PHI,YLM2,NC)
C
C  This routine computes the complex spherical harmonics using Condon and
C                  Shortley phase convention.
C
      INCLUDE 'spec.inc'
C
      COMMON /EXPFAC2/ EXPF2(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /FACTSQ/ FSQ(0:2*NL_M-2)
C
      COMPLEX YLM(0:NL,-NL:NL),COEF,YMM,YMMP,C
      COMPLEX YLM2(LINMAX)
C
      DATA SQ4PI_INV,SQR3_INV /0.282095,0.488602/
      DATA PI,SMALL /3.141593,0.0001/
C
      X=COS(THETA)
      IF(ABS(X).LT.SMALL) X=0.0
      IF(ABS(X+1.).LT.SMALL) X=-1.0
      IF(ABS(X-1.).LT.SMALL) X=1.0
C
      YLM(0,0)=CMPLX(SQ4PI_INV)
      YLM(1,0)=X*SQR3_INV
      DO L=2,NC
        Y=1./FLOAT(L)
        YLM(L,0)=X*SQRT(4.-Y*Y)*YLM(L-1,0) -
     1           (1.-Y)*SQRT(1.+2./(FLOAT(L)-1.5))*YLM(L-2,0)
      ENDDO
C
      C2=-1.
      IF((THETA.GE.0.).AND.(THETA.LE.PI)) THEN
        C=-0.5*SQRT(1.-X*X)*EXP((0.,1.)*PHI)
      ELSE
        C=0.5*SQRT(1.-X*X)*EXP((0.,1.)*PHI)
      ENDIF
C
      C1=1.
      COEF=(1.,0.)
      DO M=1,NC
        C1=C1*C2
        COEF=COEF*C
        YMM=SQ4PI_INV*COEF*FSQ(M)
        YLM(M,M)=YMM
        YLM(M,-M)=C1*CONJG(YMM)
        YMMP=X*SQRT(FLOAT(M+M+3))*YMM
        YLM(M+1,M)=YMMP
        YLM(M+1,-M)=C1*CONJG(YMMP)
        IF(M.LT.NC-1) THEN
          DO L=M+2,NC
            YLM(L,M)=(X*(L+L-1)*EXPF2(L-1,M)*YLM(L-1,M) -
     1                (L+M-1)*EXPF2(L-2,M)*YLM(L-2,M))/(EXPF2(L,M)*
     2                (L-M))
            YLM(L,-M)=C1*CONJG(YLM(L,M))
          ENDDO
        ENDIF
      ENDDO
C
      DO L=0,NC
        IL=L*L+L+1
        DO M=-L,L
          IND=IL+M
          YLM2(IND)=YLM(L,M)
        ENDDO
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE HEADERS(IUO2)
C
C   This subroutine writes headers containing the main parameters
C       of the calculation in the result file. The number of
C       lines written depends of the spectroscopy
C
C                                          Last modified : 20 Jan 2014
C
      INCLUDE 'spec.inc'
C
      INTEGER S_MUL
C
      CHARACTER*24 INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6
      CHARACTER*24 INFILE7,INFILE8,INFILE9
      CHARACTER*7  INTERACT
      CHARACTER*7  INTERACT_I,INTERACT_O
      CHARACTER*6  AUGER
      CHARACTER*3  SPECTRO,STEREO,MULTIPLET,S_O
      CHARACTER*2  ALGO1,ALGO2,ALGO3,ALGO4
      CHARACTER*1  EDGE,NLI
C
      COMMON /ALGORITHM/ ALGO1,ALGO2,ALGO3,ALGO4
      COMMON /APPROX/ NDIF,NO,ISPHER,IFWD,NTHOUT,RTHFWD(NATP_M),
     1                IBWD(NATP_M),RTHBWD(NATP_M),IPW,NCUT,PCTINT,IPP,
     2                ISPEED,IATTS,ILENGTH,RLENGTH
      COMMON /EXAFS/ NE_X,EK_INI,EK_FIN,EPH_INI
      COMMON /HEADER/ NI,NLI,AUGER,EDGE,NEDGE
      COMMON /INFILES/ INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6,
     1                 INFILE7,INFILE8,INFILE9
      COMMON /INIT_A/ LI_C,LI_I,LI_A
      COMMON /INIT_J/ JF1,JF2,I_SO,S_O
      COMMON /INIT_L/ LI,INITL,NNL,LF1,LF2,ISTEP_LF
      COMMON /INIT_L_I/ LI_1,INITL_I,NNL1,LF1_I,LF2_I,ISTEP_LF_I
      COMMON /INIT_L_O/ LI_2,INITL_O,NNL2,LF1_O,LF2_O,ISTEP_LF_O
      COMMON /INIT_M/ I_SHELL,I_MULT,L_MUL,J_MUL,S_MUL,MULTIPLET
      COMMON /MOYEN/ IMOY,NDIR,ACCEPT,ICHKDIR
      COMMON /PARCAL/ NPHI,NE,NTHETA,NFTHET,NEPS
      COMMON /RXSGEN/ IE_RX,IPHI_RX,ITHETA_RX,NE_RX,NPHI_RX,NTHETA_RX,
     1                E_INI_RX,E_FIN_RX,PHI_FIN_RX,THETA_FIN_RX
      COMMON /RXSINI/ THLUM_I,PHILUM_I,IPOL_I,NEPS_I,INTERACT_I
      COMMON /RXSFIN/ THLUM_O,PHILUM_O,IPOL_O,NEPS_O,INTERACT_O
C
      COMMON /PARCAL_A/ NPHI_A,NE_A,NTHETA_A,NFTHET_A
      COMMON /TYPCAL/ IPHI,IE,ITHETA,IFTHET,IMOD,IPOL,I_CP,
     1                I_EXT,I_TEST
      COMMON /TYPCAL_A/ IPHI_A,IE_A,ITHETA_A,IFTHET_A,IMOD_A,I_CP_A,
     1                  I_EXT_A,I_TEST_A
      COMMON /TYPEXP/ SPECTRO,INTERACT,STEREO
      COMMON /VALIN/ PHI0,E0,THETA0,THLUM,PHILUM,ELUM,VINT,NONVOL(100)
      COMMON /VALIN_AV/ I_SET,TH_0(NTH_M),PH_0(NPH_M)
      COMMON /VALFIN/ PHI1,EFIN,THETA1
      COMMON /VALEX_A/ PHI0_A,THETA0_A,PHI1_A,THETA1_A
C
      WRITE(IUO2,1)
      WRITE(IUO2,2)
C
C   Input files section:
C
C      Checking the size of filenames
C
      N_CHAR1=0
      DO J_CHAR=1,24
        IF(INFILE1(J_CHAR:J_CHAR).EQ.' ') GOTO 500
        N_CHAR1=N_CHAR1+1
      ENDDO
 500  CONTINUE
C
      N_CHAR2=0
      DO J_CHAR=1,24
        IF(INFILE2(J_CHAR:J_CHAR).EQ.' ') GOTO 501
        N_CHAR2=N_CHAR2+1
      ENDDO
 501  CONTINUE
C
      N_CHAR3=0
      DO J_CHAR=1,24
        IF(INFILE3(J_CHAR:J_CHAR).EQ.' ') GOTO 502
        N_CHAR3=N_CHAR3+1
      ENDDO
 502  CONTINUE
C
      N_CHAR4=0
      DO J_CHAR=1,24
        IF(INFILE4(J_CHAR:J_CHAR).EQ.' ') GOTO 503
        N_CHAR4=N_CHAR4+1
      ENDDO
 503  CONTINUE
C
      WRITE(IUO2,3) INFILE1(6:N_CHAR1)
      WRITE(IUO2,4) INFILE2(4:N_CHAR2)
      IF(INTERACT.NE.'NOINTER') THEN
        WRITE(IUO2,5) INFILE3(5:N_CHAR3)
      ENDIF
      WRITE(IUO2,6) INFILE4(6:N_CHAR4)
      WRITE(IUO2,2)
C
C   Type of calculation
C
      WRITE(IUO2,2)
C
      IF(SPECTRO.EQ.'PHD') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'LED') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'AED') THEN
        WRITE(IUO2,11) SPECTRO,ALGO2
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'APC') THEN
        WRITE(IUO2,15) SPECTRO,ALGO1,ALGO2
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'EIG') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
      ELSEIF(SPECTRO.EQ.'RES') THEN
        WRITE(IUO2,11) SPECTRO,ALGO1
        IF(ALGO1.EQ.'SE') THEN
          WRITE(IUO2,12) NO,NDIF,IFWD,IPW,ILENGTH
        ELSEIF(ALGO1.EQ.'CE') THEN
          WRITE(IUO2,13) NDIF
        ENDIF
        WRITE(IUO2,14) VINT
      ELSEIF(SPECTRO.EQ.'ELS') THEN
        CONTINUE
      ENDIF
C
      WRITE(IUO2,2)
C
C   Initial state parameters
C
      IF(SPECTRO.EQ.'PHD') THEN
        WRITE(IUO2,21) NI,NLI,S_O,INITL
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        WRITE(IUO2,22) EDGE,NEDGE,INITL
      ELSEIF(SPECTRO.EQ.'LED') THEN
        CONTINUE
      ELSEIF(SPECTRO.EQ.'AED') THEN
        WRITE(IUO2,24) AUGER,MULTIPLET
      ELSEIF(SPECTRO.EQ.'APC') THEN
        WRITE(IUO2,21) NI,NLI,S_O,INITL
        WRITE(IUO2,24) AUGER,MULTIPLET
      ELSEIF(SPECTRO.EQ.'RES') THEN
        WRITE(IUO2,25) NI,NLI,S_O,INITL_I,INITL_O
      ELSEIF(SPECTRO.EQ.'ELS') THEN
        CONTINUE
      ENDIF
C
      WRITE(IUO2,2)
C
C   Angular and energy parameters
C
      IF(SPECTRO.EQ.'PHD') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,34) THLUM,PHILUM,ELUM
        WRITE(IUO2,2)
        WRITE(IUO2,36)
        WRITE(IUO2,31) THETA0,THETA1
        WRITE(IUO2,32) PHI0,PHI1
        WRITE(IUO2,33) E0,EFIN
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,33) EK_INI,EK_FIN
        WRITE(IUO2,34) THLUM,PHILUM,ELUM
      ELSEIF(SPECTRO.EQ.'LED') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,31) THLUM,PHILUM
        WRITE(IUO2,2)
        WRITE(IUO2,36)
        WRITE(IUO2,31) THETA0,THETA1
        WRITE(IUO2,32) PHI0,PHI1
        WRITE(IUO2,2)
        WRITE(IUO2,33) E0,EFIN
      ELSEIF(SPECTRO.EQ.'AED') THEN
        WRITE(IUO2,36)
        WRITE(IUO2,31) THETA0_A,THETA1_A
        WRITE(IUO2,32) PHI0_A,PHI1_A
      ELSEIF(SPECTRO.EQ.'APC') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,34) THLUM,PHILUM,ELUM
        WRITE(IUO2,2)
        WRITE(IUO2,37)
        WRITE(IUO2,31) THETA0,THETA1
        WRITE(IUO2,32) PHI0,PHI1
        WRITE(IUO2,33) E0,EFIN
        WRITE(IUO2,2)
        WRITE(IUO2,38)
        WRITE(IUO2,31) THETA0_A,THETA1_A
        WRITE(IUO2,32) PHI0_A,PHI1_A
      ELSEIF(SPECTRO.EQ.'EIG') THEN
        WRITE(IUO2,33) EK_INI,EK_FIN
      ELSEIF(SPECTRO.EQ.'RES') THEN
        WRITE(IUO2,35)
        WRITE(IUO2,34) THLUM_I,PHILUM_I,ELUM
        WRITE(IUO2,2)
        WRITE(IUO2,36)
        WRITE(IUO2,31) THETA0,THETA1
        WRITE(IUO2,32) PHI0,PHI1
        WRITE(IUO2,33) E0,EFIN
      ELSEIF(SPECTRO.EQ.'ELS') THEN
        CONTINUE
      ENDIF
C
C   End of headers
C
      WRITE(IUO2,2)
      WRITE(IUO2,1)
      WRITE(IUO2,39)
C
C   Formats
C
   1  FORMAT('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',
     1       '!!!!!!!!!!!!!!!!')
   2  FORMAT('!',69X,'!')
   3  FORMAT('!',10X,'data file        :  ',A19,20X,'!')
   4  FORMAT('!',10X,'t-matrix file    :    ',A17,20X,'!')
   5  FORMAT('!',10X,'rad integral file: ',A20,20X,'!')
   6  FORMAT('!',10X,'cluster file     :  ',A19,20X,'!')
C
  11  FORMAT('!',10X,'spectroscopy     :      ',A3,8X,'algorithm : ',
     1       A2,10X,'!')
  12  FORMAT('!',15X,'NO = ',I1,'  NDIF = ',I2,'  IFWD = ',I1,
     1      '  IPW = ',I1,'  ILENGTH = ',I1,5X,'!')
  13  FORMAT('!',15X,'NDIF = ',I2,45X,'!')
  14  FORMAT('!',10X,'inner potential  :    ',F6.2,' eV',28X,'!')
  15  FORMAT('!',10X,'spectroscopy: ',A3,10X,'algorithm (photo): ',A2,
     1        11X,'!',/,'!',37X,'algorithm (auger): ',A2, 11X,'!')
C
  21  FORMAT('!',10X,'initial state    : ',I1,A1,1X,A3,
     1       ' selection rules:',' INITL = ',I2,6X,'!')
  22  FORMAT('!',10X,'initial state    : ',A1,I1,2X,' selection rules:',
     1       ' INITL = ',I2,8X,'!')
  24  FORMAT('!',10X,'initial state    : ',A6,2X,' multiplet: ',A3,
     1       17X,'!')
  25  FORMAT('!',10X,'initial state    : ',I1,A1,1X,A3,
     1       ' selection rules:', ' INITL_I = ',I2,4X,'!',/,'!',51X,
     2       ' INITL_O = ',I2,4X,'!')
C
  31  FORMAT('!',10X,'THETA_INI: ',F8.2,6X,'THETA_FIN: ',F8.2,15X,'!')
  32  FORMAT('!',10X,'PHI_INI  : ',F8.2,6X,'PHI_FIN  : ',F8.2,15X,'!')
  33  FORMAT('!',10X,'E_INI    : ',F8.2,' eV',3X,'E_FIN    : ',F8.2,
     1       ' eV',12X,'!')
  34  FORMAT('!',10X,'THETA_LUM: ',F8.2,2X,'PHI_LUM: ',F8.2,2X,
     1       'E_LUM: ',F8.2,' eV !')
  35  FORMAT('!',10X,'incoming beam    : ',40X,'!')
  36  FORMAT('!',10X,'outgoing beam    : ',40X,'!')
  37  FORMAT('!',10X,'photoelectron beam:',40X,'!')
  38  FORMAT('!',10X,'auger beam        :',40X,'!')
  39  FORMAT(71X)
C
      RETURN
C
      END
C
C=======================================================================
C
      INTEGER FUNCTION IG(J)
C
C  This function is returns the value 1 if J is an integer
C   and 2 if it is a half-integer
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      REAL*8 J,JJ
C
      DATA SMALL /0.0001D0/
C
      JJ=ABS(J+J)
C
      LL=INT(JJ+SMALL)
C
      IF(MOD(LL,2).EQ.0) THEN
        IG=1
      ELSE
        IG=2
      ENDIF
C
      END

C
C=======================================================================
C
      SUBROUTINE LOCATE(XX,N,X,J)
C
C
C            This subroutine is taken from the book :
C          "Numerical Recipes : The Art of Scientific
C           Computing" par W.H. PRESS, B.P. FLANNERY,
C               S.A. TEUKOLSKY et W.T. VETTERLING
C               (Cambridge University Press 1992)
C
C    It performs a search in an ordered table using a bisection method.
C    Given a monotonic array XX(1:N) and a value X, it returns J such
C                  that X is between XX(J) and XX(J+1).
C
      INTEGER J,N
      INTEGER JL,JM,JU
C
      REAL X,XX(N)
C
      JL=0
      JU=N+1
  10  IF(JU-JL.GT.1)THEN
        JM=(JU+JL)/2
        IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
          JL=JM
        ELSE
          JU=JM
        ENDIF
        GOTO 10
      ENDIF
      J=JL
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE LPM(E,XLPM,*)
C
C  This routine generates the electron mean free path
C
C       ILPM=-1: XLPM is set to 1.E+30
C       ILPM=0 : XLPM is the value given in the input data file
C       ILPM=1 : XLPM computed from Tokutaka et al, Surf. Sci. 149,349 (1985)
C       ILPM=2 : XLPM computed from the Seah and Dench expression
C
C                                          Last modified : 15 Sep 2009
C
      COMMON /LPMOY/ ILPM,NZ,XMAT,RHO,XLPM0
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /TESTS/ ITEST,IPRINT,ISORT1,NPATHP,ISOM
      COMMON /VALIN/ PHI0,E0,THETA0,THLUM,PHILUM,ELUM,VINT,NONVOL(100)
C
      E=E+VINT
C
      IF(ILPM.EQ.-1) THEN
        XLPM=1.E+30
      ELSEIF(ILPM.EQ.0) THEN
        XLPM=XLPM0
      ELSEIF(ILPM.EQ.1) THEN
        Q=FLOAT(NZ)*RHO/XMAT
        CSTE1=ALOG(Q/4.50)/(ALOG(7.74/4.50))
        CSTE2=ALOG(Q/3.32)/(ALOG(7.74/3.32))
        CSTE3=ALOG(Q/3.32)/(ALOG(4.50/3.32))
        A1=0.7271+0.2595*ALOG(E)
        A2=-3.2563+0.9395*ALOG(E)
        A3=-2.5716+0.8226*ALOG(E)
        IF(E.GE.350.) GO TO 10
        XLN=CSTE1*(0.0107-0.0083*ALOG(E))+A1
        GO TO 20
  10    IF((NZ.GE.24).AND.(NZ.LE.74)) GO TO 30
        XLN=CSTE2*(1.6551-0.2890*ALOG(E))+A2
        GO TO 20
  30    IF(NZ.GE.42) GO TO 40
        XLN=CSTE3*(0.6847-0.1169*ALOG(E))+A2
        GO TO 20
  40    XLN=CSTE1*(0.9704-0.1721*ALOG(E))+A3
  20    XLPM=EXP(XLN)
      ELSEIF(ILPM.EQ.2) THEN
        XLPM=1430./(E**2)+0.54*SQRT(E)
      ELSE
        RETURN 1
      ENDIF
C
      E=E-VINT
      IF(IPRINT.GT.0) WRITE(IUO1,80) E,XLPM
C
  80  FORMAT(/////,2X,'=========  E = ',F7.2,' eV',5X,'MEAN',
     1' FREE PATH = ',F6.3,' ANGSTROEMS  ','=========')
C
      RETURN
C
      END


C
C=======================================================================
C
      SUBROUTINE N_J(J1,MJ1,J2,MJ2,MJ6,NJ,I_INT,N_IN)
C
C   This subroutine calculates Wigner's 3j and 6j coefficients
C    using a downward recursion scheme due to Schulten and Gordon.
C    The 3j are defined as (J1 J2 J) where in fact L1=MJ1, etc are
C                          (L1 L2 L)
C    azimuthal quantum numbers, and the 6j as {J1 J2 J} where now
C                                             {L1 L2 L}
C    J1, L1, etc are the same kind of orbital quantum numbers.
C    The result is stored as NJ(J).
C
C    The parameter N allows to choose between 3j and 6j calculation, and
C    Clebsch-Gordan. It can take the values :
C
C                    N = 2 ----> Clebsch-Gordan
C                    N = 3 ----> Wigner's 3j
C                    N = 6 ----> Wigner's 6j
C
C    The Clebsch-Gordan coefficients are related to Wigner's 3j through :
C
C       CG(J1,M1,J2,M2|J,MJ) = ( J1 J2  J  )*sqrt(2*J+1)*(-1)**(J1-J2+MJ)
C                              ( M1 M2 -MJ )
C    I_INT is a flag that returns 1 if the index J of the nj symbol
C      is integer and 0 if it is a half integer.
C
C   Note : For 3j, MJ6 is ignored while for 6j, we have :
C
C                J1=J1   MJ1=L1   J2=J2   MJ2=L2   MJ6=L
C
C   Ref. : K. Schulten and R. G. Gordon, J. Math. Phys. 16, 1961 (1975)
C
C                      Last modified :  8 Dec 2008 ----> D. Sebilleau
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'spec.inc'
C
      REAL*4 NJ(0:N_GAUNT)
C
      REAL*8 J1,J2,J,MJ1,MJ2,MJ,JP1,JP2
      REAL*8 F(0:N_GAUNT),A(0:N_GAUNT),B(0:N_GAUNT)
      REAL*8 JL12,JK12,MJ6,SIG,GLD
      REAL*8 JJ1,JJ2,JL1,JL2,JL3,JJ12,JJ_MIN,JJ_MAX
C
      COMMON /LOGAMAD/ GLD(0:N_GAUNT,2)
C
      DATA SMALL /0.0001D0/
C
      IS=0
C
      IF(N_IN.EQ.2) THEN
        N_OU=3
        I_CG=1
      ELSE
        N_OU=N_IN
        I_CG=0
      ENDIF
C
      IF(N_OU.EQ.3) THEN
C
C------------------------------  3j case  ---------------------------------
C
C
C  Test to check if J1 and J2 are integer or semi-integer
C
C     Integer      : IG=1
C     Half-integer : IG=2
C
C  Each angular momentum J is represented by the integer index L and
C    the corresponding MJ by M
C
        L1=INT(J1+SMALL)
        L2=INT(J2+SMALL)
        M1=INT(MJ1+SIGN(SMALL,MJ1))
        M2=INT(MJ2+SIGN(SMALL,MJ2))
        DIF1=J1-DFLOAT(L1)
        DIF2=J2-DFLOAT(L2)
C
C  IGx is a flag telling the code which case of Gamma function to use :
C
C     IGx = 1 : integer case
C     IGx = 2 : half-integer case
C
        IF(ABS(DIF1).LT.SMALL) THEN
          IG1=1
        ELSE
          IG1=2
        ENDIF
        IF(ABS(DIF2).LT.SMALL) THEN
          IG2=1
        ELSE
          IG2=2
        ENDIF
        IF(IG1.EQ.IG2) THEN
          IGG=1
          IF(IG1.EQ.2) IS=1
        ELSE
          IGG=2
        ENDIF
C
C  Here, we assume that (J1,J2) are both either integer or half-integer
C  If J is integer, the corresponding index is L = j (for loops or storage)
C    while if J is an half-integer, this index is L=  j - 1/2 = int(j)
C
C  Integer indices are used for loops and for storage while true values
C    are used for the initial values. When J1 and J2 are both half-integers,
C    the values of J are integer and L should be increased by 1
C
        JL12=J1+J2
        JK12=J1-J2
C
        L12=INT(JL12 + SIGN(SMALL,JL12))
        K12=INT(JK12 + SIGN(SMALL,JK12))
C
        LM1=INT(J1+MJ1 + SIGN(SMALL,J1+MJ1))
        LM2=INT(J2+MJ2 + SIGN(SMALL,J2+MJ2))
        KM1=INT(J1-MJ1 + SIGN(SMALL,J1-MJ1))
        KM2=INT(J2-MJ2 + SIGN(SMALL,J2-MJ2))
C
        MJ=-MJ1-MJ2
C
        M=INT(MJ+SIGN(SMALL,MJ))
        L12M=INT(JL12+MJ+SIGN(SMALL,JL12+MJ))
        K12M=INT(JL12-MJ+SIGN(SMALL,JL12-MJ))
        L1_2=INT(J1+J1+SIGN(SMALL,J1))
        L2_2=INT(J2+J2+SIGN(SMALL,J2))
        L12_2=INT(JL12+JL12+SIGN(SMALL,JL12))
C
        IF(IG(JL12).EQ.1) THEN
          I_INT=1
        ELSE
          I_INT=0
        ENDIF
C
C  Initialisation of the 3j symbol NJ(J) = (J1  J2   J)
C                                          (MJ1 MJ2 MJ)
C
        DO L=0,L12
          NJ(L)=0.
        ENDDO
C
        IF((ABS(MJ1).GT.J1).OR.(ABS(MJ2).GT.J2)) GOTO 10
C
C  Initial values (J1+J2+1) and (J1+J2) for J to be used in the downward
C     recursion scheme. This scheme writes as
C
C     J A(J+1) NJ(J+1) + B(J) NJ(J) + (J+1) A(J) NJ(J-1) = 0
C
        F(L12+1)=0.D0
        A(L12+1)=0.D0
        D1=GLD(L2_2+1,1)-GLD(L12_2+2,1)
        D2=GLD(L1_2+1,1)-GLD(LM2+1,1)
        D3=GLD(L12M+1,1)-GLD(KM2+1,1)
        D4=GLD(K12M+1,1)-GLD(LM1+1,1)-GLD(KM1+1,1)
C
        N12=INT(JK12-MJ + SIGN(SMALL,JK12-MJ))
C
        IF(I_CG.EQ.1) THEN
          IF(MOD(N12,2).EQ.0) THEN
            SIG=1.D0
          ELSE
            SIG=-1.D0
          ENDIF
        ENDIF
C
        IF(MOD(N12,2).EQ.0) THEN
          F(L12)=DSQRT(DEXP(D1+D2+D3+D4))
        ELSE
          F(L12)=-DSQRT(DEXP(D1+D2+D3+D4))
        ENDIF
C
        A(L12)=2.D0*DSQRT(J1*J2*(1.D0+2.D0*JL12)*(JL12*JL12-MJ*MJ))
        B(L12)=-(2.D0*JL12+1.D0)*((J1*J1-J2*J2+JK12)*MJ-JL12*
     1          (JL12+1.D0)*(MJ2-MJ1))
C
        IF(ABS(M).LE.L12) THEN
          IF(I_CG.EQ.0) THEN
            NJ(L12)=SNGL(F(L12))
          ELSE
            NJ(L12)=SNGL(F(L12)*SIG*DSQRT(JL12+JL12+1.D0))
          ENDIF
        ELSE
          NJ(L12)=0.
        ENDIF
C
        LMIN=MAX0(ABS(K12),ABS(M))
C
C  Downward recursion for NJ(J)
C
        DO L=L12-1,LMIN,-1
          LP1=L+1
          LP2=L+2
C
C  Value of the angular momentum J corresponding to the loop index L
C
          IF(IGG.EQ.1) THEN
            J=DFLOAT(L)
            JP1=DFLOAT(LP1)
            JP2=DFLOAT(LP2)
          ELSE
            J=DFLOAT(L) + 0.5D0
            JP1=DFLOAT(LP1) + 0.5D0
            JP2=DFLOAT(LP2) + 0.5D0
          ENDIF
C
          A(L)=DSQRT((J*J-JK12*JK12)*((JL12+1.D0)*(JL12+1.D0)-J*J)
     1               *(J*J-MJ*MJ))
          B(L)=-(2.D0*J+1.D0)*(J1*(J1+1.D0)*MJ-J2*(J2+1.D0)*MJ-
     1           J*JP1*(MJ2-MJ1))
          F(L)=-(JP1*A(LP2)*F(LP2)+B(LP1)*F(LP1))/(JP2*A(LP1))
C
          IF(ABS(MJ).LE.J) THEN
            IF(I_CG.EQ.0) THEN
              NJ(L)=SNGL(F(L))
            ELSE
              NJ(L)=SNGL(F(L)*SIG*DSQRT(J+J+1.D0))
            ENDIF
          ELSE
            NJ(L)=0.
          ENDIF
C
        ENDDO
C
  10    CONTINUE
C
      ELSEIF(N_OU.EQ.6) THEN
C
C------------------------------  6j case  ---------------------------------
C
C  Change of notation for greater readability ---> NJ(JJ)
C
C    True angular momentum value : begins with a J (JJn,JLn)
C    Corresponding integer storage and loop index : begins by L (LJn,LLn)
C
        JJ1=J1
        JJ2=J2
        JL1=MJ1
        JL2=MJ2
        JL3=MJ6
C
        LJ1=INT(JJ1+SIGN(SMALL,JJ1))
        LJ2=INT(JJ2+SIGN(SMALL,JJ2))
        LL1=INT(JL1+SIGN(SMALL,JL1))
        LL2=INT(JL2+SIGN(SMALL,JL2))
        LL3=INT(JL3+SIGN(SMALL,JL3))
C
        JJ12=JJ1-JJ2
        JL12=JL1-JL2
C
        LJ12=INT(JJ12+SIGN(SMALL,JJ12))
        LL12=INT(JL12+SIGN(SMALL,JL12))
C
        JJ_MIN=MAX(ABS(LJ12),ABS(LL12))
        JJ_MAX=MIN(JJ1+JJ2,JL1+JL2)
        LJJ_MIN=INT(JJ_MIN+SIGN(SMALL,JJ_MIN))
        LJJ_MAX=INT(JJ_MAX+SIGN(SMALL,JJ_MAX))
C
C  Initialisation of the 6j symbol NJ(J) = {J1  J2  J }
C                                          {L1  L2  L3}
C
        DO L=0,LJJ_MAX
          NJ(L)=0.
        ENDDO
C
C  Initial values (J1+J2+1) and (J1+J2) for J to be used in the downward
C     recursion scheme. This scheme writes as
C
C     J A(J+1) NJ(J+1) + B(J) NJ(J) + (J+1) A(J) NJ(J-1) = 0
C
C  There are two possible initial values as max(|J1-J2|,|L1-L2|) <= J <=
C    min(J1+J2,L1+L2) :
C
C    {J1  J2  L1+L2}  and  {J1  J2  J1+J2}  =  {L1 L2 J1+J2}
C    {L1  L2    L3 }       {L1  L2    L3 }     {J1 J2   L3 }
C
C    They can be calculated from equation (6.3.1) of Edmonds page 97
C
        F(LJJ_MAX+1)=0.D0
        A(LJJ_MAX+1)=0.D0
C
        IF(ABS(JJ_MAX-JL1-JL2).LT.SMALL) THEN
          F(LJJ_MAX)=SIXJ_IN(JJ1,JJ2,JL1,JL2,JL3)
        ELSE
          F(LJJ_MAX)=SIXJ_IN(JL1,JL2,JJ1,JJ2,JL3)
        ENDIF
        NJ(LJJ_MAX)=SNGL(F(LJJ_MAX))
C
        A(LJJ_MAX)=SQRT((JJ_MAX*JJ_MAX-(JJ1-JJ2)*(JJ1-JJ2))*
     1                  ((JJ1+JJ2+1.D0)*(JJ1+JJ2+1.D0)-JJ_MAX*JJ_MAX)*
     2                  (JJ_MAX*JJ_MAX-(JL1-JL2)*(JL1-JL2))*
     3                  ((JL1+JL2+1.D0)*(JL1+JL2+1.D0)-JJ_MAX*JJ_MAX))
        B(LJJ_MAX)=(JJ_MAX+JJ_MAX+1.D0)*(JJ_MAX*(JJ_MAX+1.D0)*
     1             (-JJ_MAX*(JJ_MAX+1.D0)+JJ1*(JJ1+1.D0)+
     2               JJ2*(JJ2+1.D0))+
     3             JL1*(JL1+1.D0)*(JJ_MAX*(JJ_MAX+1.D0)+
     4               JJ1*(JJ1+1.D0)-JJ2*(JJ2+1.D0))+
     5             JL2*(JL2+1.D0)*(JJ_MAX*(JJ_MAX+1.D0)-
     6               JJ1*(JJ1+1.D0)+JJ2*(JJ2+1.D0))-
     6             (JJ_MAX+JJ_MAX)*(JJ_MAX+1.D0)*JL3*(JL3+1.D0))
C
        IF(IG(JJ_MAX).EQ.1) THEN
          I_INT=1
        ELSE
          I_INT=0
        ENDIF
C
C  Downward recurrence relation
C
        DO L=LJJ_MAX-1,LJJ_MIN,-1
          LP1=L+1
          LP2=L+2
C
C  Value of the angular momentum J corresponding to the loop index L
C
          IF(IG(JJ_MAX).EQ.1) THEN
            J=DFLOAT(L)
            JP1=DFLOAT(LP1)
            JP2=DFLOAT(LP2)
          ELSE
            J=DFLOAT(L) + 0.5D0
            JP1=DFLOAT(LP1) + 0.5D0
            JP2=DFLOAT(LP2) + 0.5D0
          ENDIF
C
          A(L)=SQRT((J*J-(JJ1-JJ2)*(JJ1-JJ2))*
     1              ((JJ1+JJ2+1.D0)*(JJ1+JJ2+1.D0)-J*J)*
     2              (J*J-(JL1-JL2)*(JL1-JL2))*
     3              ((JL1+JL2+1.D0)*(JL1+JL2+1.D0)-J*J))
          B(L)=(J+J+1)*(J*JP1*(-J*JP1+JJ1*(JJ1+1.D0)+JJ2*(JJ2+1.D0))+
     1                  JL1*(JL1+1.D0)*(J*JP1+JJ1*(JJ1+1.D0)-
     2                    JJ2*(JJ2+1.D0))+
     3                  JL2*(JL2+1.D0)*(J*JP1-JJ1*(JJ1+1.D0)+
     4                    JJ2*(JJ2+1.D0))-
     5                  (J+J)*JP1*JL3*(JL3+1.D0))
C
          F(L)=-(JP1*A(LP2)*F(LP2)+B(LP1)*F(LP1))/(JP2*A(LP1))
          NJ(L)=SNGL(F(L))
C
        ENDDO
C
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE ORDRE2(NINI,VALINI,NFIN,VALFIN)
C
C  Given a set of **integer** numbers VALINI, this routine orders them
C       and suppresses the values appearing more than once. The remaining
C       values are stored in VALFIN.
C
C       VALINI(K+1).GT.VALINI(K) : decreasing order
C       VALINI(K+1).LT.VALINI(K) : increasing order
C
C
C
      INTEGER VALINI(NINI),VALFIN(NINI),R1
C
      LOGICAL BUBBLE
C
      DO J=1,NINI-1
        K=J
        BUBBLE=.TRUE.
150     IF(K.GE.1.AND.BUBBLE) THEN
          IF(VALINI(K+1).LT.VALINI(K)) THEN
            R1=VALINI(K)
            VALINI(K)=VALINI(K+1)
            VALINI(K+1)=R1
          ELSE
            BUBBLE=.FALSE.
          ENDIF
          K=K-1
          GOTO 150
        ENDIF
      ENDDO
C
      JFIN=1
      VALFIN(1)=VALINI(1)
      DO J=1,NINI-1
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)).GT.0) THEN
          JFIN=JFIN+1
          VALFIN(JFIN)=VALINI(J+1)
        ENDIF
      ENDDO
      NFIN=JFIN
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE ORDRE(NINI,VALINI,NFIN,VALFIN)
C
C  Given a set of **real** numbers VALINI, this routine orders them and
C       suppresses the values appearing more than once. The remaining
C       values are stored in VALFIN.
C
C       VALINI(K+1).GT.VALINI(K) : decreasing order
C       VALINI(K+1).LT.VALINI(K) : increasing order
C
C
      DIMENSION VALINI(NINI),VALFIN(NINI)
C
      LOGICAL BUBBLE
C
      DATA SMALL /0.0001/
C
      DO J=1,NINI-1
        K=J
        BUBBLE=.TRUE.
150     IF(K.GE.1.AND.BUBBLE) THEN
          IF(VALINI(K+1).GT.VALINI(K)) THEN
            R1=VALINI(K)
            VALINI(K)=VALINI(K+1)
            VALINI(K+1)=R1
          ELSE
            BUBBLE=.FALSE.
          END IF
          K=K-1
          GOTO 150
        ENDIF
      ENDDO
C
      JFIN=1
      VALFIN(1)=VALINI(1)
      DO J=1,NINI-1
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)).GT.SMALL) THEN
          JFIN=JFIN+1
          VALFIN(JFIN)=VALINI(J+1)
        ENDIF
      ENDDO
      NFIN=JFIN
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE PLM(X,PLMM,NC)
C
C  This routine computes the Legendre functions. It is a modified version
C       of that written by W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY and
C       W.T. VETTERLING in "Numerical Recipes : The Art of Scientific
C           Computing" (Cambridge University Press 1992).
C
      DIMENSION PLMM(0:100,0:100)
C
      PLMM(0,0)=1.
      PLMM(1,0)=X
      DO L=2,NC
        PLMM(L,0)=(X*(L+L-1)*PLMM(L-1,0)-(L-1)*PLMM(L-2,0))/L
      ENDDO
C
      DO M=1,NC
        PMM=1.
        FACT=1.
        SOMX2=SQRT(1.-X*X)
        FACT=1.
        DO I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.
        ENDDO
        PMMP1=X*FACT*PMM
        PLMM(M,M)=PMM
        PLMM(M+1,M)=PMMP1
        IF(M.LT.NC-1) THEN
          DO L=M+2,NC
            PLL=(X*(L+L-1)*PMMP1-(L+M-1)*PMM)/(L-M)
            PMM=PMMP1
            PMMP1=PLL
            PLMM(L,M)=PLL
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
C
      END
C
C=============================================================================
C
      SUBROUTINE POLHAN(ISPHER,NO,NC,RHO,HLM)
C
C  This routine calculates a function HLM(L,M), related to the the Hankel
C      polynomials and their derivatives with respect to z=1/ikr,
C      necessary for the Rehr-Albers expansion of the propagator.
C
      INCLUDE 'spec.inc'
C
      COMPLEX HLM(0:NO_ST_M,0:NL_M-1),RHO,Z,ONEC
C
      ONEC=(1.,0.)
C
      IF(ISPHER.GE.1) THEN
        Z=(0.,-1.)/RHO
C
C  Case M = 0
C
        HLM(0,0)=ONEC
        HLM(0,1)=ONEC-Z
        DO L=2,NC
          HLM(0,L)=HLM(0,L-2)-FLOAT(L+L-1)*Z*HLM(0,L-1)
        ENDDO
C
C  Case M > 0
C
        IF(NO.GE.1) THEN
          DO M=1,NO
            HLM(M,M)=-Z*HLM(M-1,M-1)*FLOAT(M+M-1)
            HLM(M,M+1)=HLM(M,M)*FLOAT(M+M+1)*(ONEC-Z*FLOAT(M+1))
            DO L=M+2,NC
              HLM(M,L)=HLM(M,L-2)-FLOAT(L+L-1)*Z*(HLM(M,L-1)+
     1                                            HLM(M-1,L-1))
            ENDDO
          ENDDO
        ENDIF
      ELSE
        DO M=0,NO
          DO L=M,NC
            HLM(M,L)=ONEC
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE POLLEG(NC,X,PL)
C
C  This routine computes the Legendre polynomials up to order NC-1
C
      DIMENSION PL(0:100)
C
      PL(0)=1.
      PL(1)=X
      DO 10 L=2,NC-1
        L1=L-1
        L2=L-2
        L3=2*L-1
        PL(L)=(X*FLOAT(L3)*PL(L1)-FLOAT(L1)*PL(L2))/FLOAT(L)
  10  CONTINUE
C
      RETURN
C
      END
C
C=======================================================================
C
      FUNCTION PRSCAL(A1,A2)
C
C  This function computes the dot product of the two vectors A1 and A2
C
      DIMENSION A1(3),A2(3)
C
      PRSCAL=A1(1)*A2(1)+A1(2)*A2(2)+A1(3)*A2(3)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE PRVECT(A1,A2,A3,C)
C
C  This function computes the vector product of the two vectors A1 and A2.
C            The result is A3; C is a scaling factor
C
      DIMENSION A1(3),A2(3),A3(3)
C
      A3(1)=(A1(2)*A2(3)-A1(3)*A2(2))/C
      A3(2)=(A1(3)*A2(1)-A1(1)*A2(3))/C
      A3(3)=(A1(1)*A2(2)-A1(2)*A2(1))/C
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE READ_DATA(ICOM,NFICHLEC,JFICH,ITRTL,*,*,*,*,*,*,*,*,
     1                     *,*,*,*,*)
C
C  This subroutine reads the input data from unit ICOM and writes
C     them in the control file IUO1. Then, it stores the data in 
C                     the various COMMON blocks
C
C                                          Last modified : 20 Jan 2014
C
      INCLUDE 'spec.inc'
C
      COMPLEX*16 ALPHA,BETA 
C
      REAL*8 J1,J2,MJ1,MJ2,MJ3,JJ,DEXPF2,DFSQ,DXDEN,DEXPF
      REAL*8 JJ_MIN,JJ_MAX,JJ12,JL12,SMALL,SQPI,ACC,EXPO,SHIFT
      REAL*8 GLD,FACT1L,FACT2L
C
      REAL THFWD_IN(NATP_M),THFWD_EX(NATP_M)
      REAL THFWD_O1(NATP_M),THFWD_O2(NATP_M)
      REAL THBWD_IN(NATP_M),THBWD_EX(NATP_M)
      REAL THBWD_O1(NATP_M),THBWD_O2(NATP_M)
      REAL GLG(0:N_GAUNT),NJ(0:N_GAUNT),PCRELA_TI(3)
      REAL ALPHAR,BETAR,RACC
C
      INTEGER S_MUL,ISPEC(2),EXC,OUT1,OUT2,TIP,NONVOL_TI(100)
C
      CHARACTER*24 INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6
      CHARACTER*24 INFILE7,INFILE8,INFILE9,INFILE10,INFILE11
      CHARACTER*24 INFILE12,INFILE13,INFILE14,INFILE15,INFILE16
      CHARACTER*24 INFILE17
      CHARACTER*24 OUTFILE1,OUTFILE2,OUTFILE3,OUTFILE4
      CHARACTER*24 INDATA
      CHARACTER*7 TESLEC,RIEN,INTERACT,INTERACT_I,INTERACT_O
      CHARACTER*6 AUGER,AUGER1
      CHARACTER*4 STRING,SPECNAME(10),METHOD
      CHARACTER*4 TEXTE1(10),TEXTE2(12),TEXTE3(12),TEXTE4(12),TEXTE5(12)
      CHARACTER*4 TEXTE6(12),TEXTE7(12),TEXTE7B(12),TEXTE7C(12)
      CHARACTER*4 TEXTE7D(12),TEXTE7E(12),TEXTE7F(12),TEXTE8(12)
      CHARACTER*4 TEXTE10(12)
      CHARACTER*3 CODRES(8),CODCTR(7),CRIST,CENTR,UNIT,UNIT_TI
      CHARACTER*3 UNLENGTH_IN,UNLENGTH_EX,UNLENGTH_O1,UNLENGTH_O2
      CHARACTER*3 S_O,SPECTRO,MULTIPLET,STEREO
      CHARACTER*2 CHEM
      CHARACTER*1 NLI,EDGE,EDGE_C,EDGE_I,EDGE_A,MULT
C
      COMMON /ADSORB/ IADS,NATA,NADS1,NADS2,NADS3,ADS(3,900),NCOUCH
      COMMON /ADSORB_TI/ IADS_TI,NATA_TI,NADS1_TI,NADS2_TI,NADS3_TI,
     1                ADS_TI(3,900),NCOUCH_TI
      COMMON /AMPLI/ I_AMP
      COMMON /APPROX_IN/ NDIF_IN,NO_IN,ISPHER_IN,IFWD_IN,NTHOUT_IN,
     1                RTHFWD_IN(NATP_M),IBWD_IN(NATP_M),
     2                RTHBWD_IN(NATP_M),IPW_IN,NCUT_IN,PCTINT_IN,IPP_IN,
     3                ISPEED_IN,IATTS_IN,ILENGTH_IN,RLENGTH_IN
      COMMON /APPROX_EX/ NDIF_EX,NO_EX,ISPHER_EX,IFWD_EX,NTHOUT_EX,
     1                RTHFWD_EX(NATP_M),IBWD_EX(NATP_M),
     2                RTHBWD_EX(NATP_M),IPW_EX,NCUT_EX,PCTINT_EX,IPP_EX,
     3                ISPEED_EX,IATTS_EX,ILENGTH_EX,RLENGTH_EX
      COMMON /APPROX_O1/ NDIF_O1,NO_O1,ISPHER_O1,IFWD_O1,NTHOUT_O1,
     1                RTHFWD_O1(NATP_M),IBWD_O1(NATP_M),
     2                RTHBWD_O1(NATP_M),IPW_O1,NCUT_O1,PCTINT_O1,IPP_O1,
     3                ISPEED_O1,IATTS_O1,ILENGTH_O1,RLENGTH_O1
      COMMON /APPROX_O2/ NDIF_O2,NO_O2,ISPHER_O2,IFWD_O2,NTHOUT_O2,
     1                RTHFWD_O2(NATP_M),IBWD_O2(NATP_M),
     2                RTHBWD_O2(NATP_M),IPW_O2,NCUT_O2,PCTINT_O2,IPP_O2,
     3                ISPEED_O2,IATTS_O2,ILENGTH_O2,RLENGTH_O2
      COMMON /ATOMS/ VALZ(NATCLU_M),VALZ_MAX,NZAT(NATCLU_M),
     1               I_GR,I_INV,CHEM(NATCLU_M)
      COMMON /ATOMS_TI/ VALZ_TI(NATCLU_M),VALZ_MAX_TI,NZAT_TI(NATCLU_M),
     1               I_GR_TI,I_INV_TI,CHEM_TI(NATCLU_M)
      COMMON /AUGER/ NLIN_A(0:20),L_BOUNDS(0:20,2),AUGER
      COMMON /BASES/ ATBAS(3*NATP_M),VECBAS(9)
      COMMON /COEFRLM/ CF(0:2*NL_M-2,0:2*NL_M-2,0:2*NL_M-2)
      COMMON /CONVACC/ ALPHA,BETA,I_XN,I_VA,I_GN,I_WN,LEVIN
      COMMON /CONVTYP/ SHIFT,ACC,EXPO,I_PWM,I_ACC,N_ONE,N_MAX,N_ITER,
     1                 N_TABLE,METHOD
      COMMON /C_G/ CG(-LI_M:LI_M+1,2,2)
      COMMON /C_G_A/ CGA(0:NCG_M,-NCG_M:NCG_M,0:NCG_M,-NCG_M:NCG_M,
     1                   0:2*NCG_M)
      COMMON /C_G_M/ CG_S(2,2,2)
      COMMON /CONVOL/ G_MAX,G_HOL,G_SLO,G_CEN
      COMMON /CLUELEC/ NEL,NCL,INC,EXC,OUT1,OUT2,TIP,CRYS
      COMMON /CRANGL/ ALPHAD,BETAD,GAMMAD
      COMMON /DEBWAL/ IDCM,IDWSPH,TD,QD,T,RSJ,UJ2(NATM)
      COMMON /DEBWAL_TI/ IDCM_TI,IDWSPH_TI,TD_TI,QD_TI,T_TI,
     1                   RSJ_TI,UJ2_TI(NATM)
      COMMON /DEXPFAC2/ DEXPF2(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /DFACTSQ/ DFSQ(0:2*NL_M-2)
      COMMON /EIGEN/ NE_EIG,EIG_INI,EIG_FIN,I_VIB,I_MFP
      COMMON /EXAFS/ NE_X,EK_INI,EK_FIN,EPH_INI
      COMMON /EXPFAC/ EXPF(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /EXPFAC2/ EXPF2(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /EXPROT/ EXPR(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /FACTSQ/ FSQ(0:2*NL_M-2)
      COMMON /FDIF/ R0,R1
      COMMON /FIXSCAN/ N_FIXED,N_SCAN,IPH_1,FIX0,FIX1,SCAN0,SCAN1
      COMMON /FIXSCAN_A/ N_FIXED_A,N_SCAN_A,IPH_1_A,FIX0_A,FIX1_A,
     1                   SCAN0_A,SCAN1_A
      COMMON /HEADER/ NI,NLI,AUGER1,EDGE,NEDGE
      COMMON /INDAT/ INDATA(100)
      COMMON /INFILES2/ INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6,
     1                  INFILE7,INFILE8,INFILE9,INFILE10,INFILE11,
     2                  INFILE12,INFILE13,INFILE14,INFILE15,INFILE16,
     2                  INFILE17
      COMMON /INUNITS2/ IUI1,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,IUI9,
     1                  IUI10,IUI11,IUI12,IUI13,IUI14,IUI15,IUI16,IUI17
      COMMON /INIT_A/ LI_C,LI_I,LI_A
      COMMON /INIT_J/ JF1,JF2,I_SO,S_O
      COMMON /INIT_L/ LI,INITL,NNL,LF1,LF2,ISTEP_LF
      COMMON /INIT_L_I/ LI_1,INITL_I,NNL1,LF1_I,LF2_I,ISTEP_LF_I
      COMMON /INIT_L_O/ LI_2,INITL_O,NNL2,LF1_O,LF2_O,ISTEP_LF_O
      COMMON /INIT_M/ I_SHELL,I_MULT,L_MUL,J_MUL,S_MUL,MULTIPLET
      COMMON /LIMAMA/ NIV,COUPUR
      COMMON /LINLBD_IN/ LBD_IN(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX_IN,
     1                   NUMAX_IN(NATP_M)
      COMMON /LINLBD_EX/ LBD_EX(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX_EX,
     1                   NUMAX_EX(NATP_M)
      COMMON /LINLBD_O1/ LBD_O1(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX_O1,
     1                   NUMAX_O1(NATP_M)
      COMMON /LINLBD_O2/ LBD_O2(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX_O2,
     2                   NUMAX_O2(NATP_M)
      COMMON /LOGAMAD/ GLD(0:N_GAUNT,2)
      COMMON /LOOP_TEMP/ ITEMP,NTEMP,TEMP0,TEMP1
      COMMON /LPMOY/ ILPM,NZA,XM,RH,XLPM0
      COMMON /MILLER/ IH,IK,II,IL,IVG0,IVN(3)
      COMMON /MOYEN/ IMOY,NDIR,ACCEPT,ICHKDIR
      COMMON /MOYEN_A/ IMOY_A,NDIR_A,ACCEPT_A,ICHKDIR_A
      COMMON /OUTFILES/ OUTFILE1,OUTFILE2,OUTFILE3,OUTFILE4
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /PARCAL/ NPHI,NE,NTHETA,NFTHET,NEPS
      COMMON /PARCAL_A/ NPHI_A,NE_A,NTHETA_A,NFTHET_A
      COMMON /RA_IN/ I_NO_IN,I_RA_IN,N_RA_IN(NATP_M)
      COMMON /RA_EX/ I_NO_EX,I_RA_EX,N_RA_EX(NATP_M)
      COMMON /RA_O1/ I_NO_O1,I_RA_O1,N_RA_O1(NATP_M)
      COMMON /RA_O2/ I_NO_O2,I_RA_O2,N_RA_O2(NATP_M)
      COMMON /RELADS/ NRELA,PCRELA(3)
      COMMON /RELAX/ IREL,NREL,PCREL(10),OMEGA1,OMEGA2
      COMMON /RESEAU/ NCRIST,NCENTR,IBAS,NAT,A,BSURA,CSURA,UNIT
      COMMON /RESEAU_TI/ NAT_TI,A_TI,UNIT_TI
      COMMON /RXSGEN/ IE_RX,IPHI_RX,ITHETA_RX,NE_RX,NPHI_RX,NTHETA_RX,
     1                E_INI_RX,E_FIN_RX,PHI_FIN_RX,THETA_FIN_RX
      COMMON /RXSINI/ THLUM_I,PHILUM_I,IPOL_I,NEPS_I,INTERACT_I
      COMMON /RXSFIN/ THLUM_O,PHILUM_O,IPOL_O,NEPS_O,INTERACT_O
      COMMON /SPECTRUM/ I_SPECTRUM(NE_M)
      COMMON /SPIN_IN/ ISPIN_IN,IDICHR_IN,NSPIN_IN,NSPIN2_IN,ISFLIP_IN,
     1                 IR_DIA_IN,NSTEP_IN
      COMMON /SPIN_EX/ ISPIN_EX,IDICHR_EX,NSPIN_EX,NSPIN2_EX,ISFLIP_EX,
     1                 IR_DIA_EX,NSTEP_EX
      COMMON /SPIN_O1/ ISPIN_O1,IDICHR_O1,NSPIN_O1,NSPIN2_O1,ISFLIP_O1,
     1                 IR_DIA_O1,NSTEP_O1
      COMMON /SPIN_O2/ ISPIN_O2,IDICHR_O2,NSPIN_O2,NSPIN2_O2,ISFLIP_O2,
     1                 IR_DIA_O2,NSTEP_O2
      COMMON /TESTS/ ITEST,IPRINT,ISORT1,NPATHP,ISOM
      COMMON /TESTS_TI/ ITEST_TI,IPRINT_TI,ISORT1_TI,NPATHP_TI,ISOM_TI
      COMMON /TL_TRUNC/ ITRTL_IN,ITRTL_EX,ITRTL_O1,ITRTL_O2
      COMMON /TYPCAL/ IPHI,IE,ITHETA,IFTHET,IMOD,IPOL,I_CP,
     1                I_EXT,I_TEST
      COMMON /TYPCAL_A/ IPHI_A,IE_A,ITHETA_A,IFTHET_A,IMOD_A,I_CP_A,
     1                  I_EXT_A,I_TEST_A
      COMMON /TYPEM/ NEMET,IESURF,IEMET(NEMET_M)
      COMMON /TYPEM_TI/ NEMET_TI,IESURF_TI,IEMET_TI(NEMET_M)
      COMMON /TYPEXP/ SPECTRO,INTERACT,STEREO
      COMMON /UJ2_T/  UJ2T(NATP_M,NTEMP_M)
      COMMON /VALIN/ PHI0,E0,THETA0,THLUM,PHILUM,ELUM,VINT,NONVOL(100)
      COMMON /VALIN_AV/ I_SET,TH_0(NTH_M),PH_0(NPH_M)
      COMMON /VALFIN/ PHI1,EFIN,THETA1
      COMMON /VALEX_A/ PHI0_A,THETA0_A,PHI1_A,THETA1_A
      COMMON /XMRHO/ XMAT(0:99),RHOAT(0:99)
C
      DATA CODRES /'CUB','TET','ORB','MNC','TCN','TRG','HEX','EXT'/
      DATA CODCTR /'P','I','F','R','A','B','C'/
      DATA PIS180,BOHR /0.017453,0.529177/  
      DATA SPECNAME /'PhD ','LEED','EXAF','AED ','    ','REXS','EELS',
     1               'STM ','    ','EIGE'/
      DATA SQPI,SMALL /1.772453850906D0,1.D-6/
C
      NLINE=1000
      I_EXT=0
      I_EXT_A=0
      IVG0=0
      IRET=0
      NCRIST=0
      NCENTR=0
      I_SO=0
      STEREO=' NO'
      ISPEC(2)=0
      IBAS=0
      IADS_TI=0
C
      DO J=1,10
        PCREL(J)=0.
      ENDDO
C
      IDWSPH=0
      IDWSPH_TI=0
      ISPHER_IN=0
      ISPHER_EX=0
      ISPHER_O1=0
      ISPHER_O2=0
C
C   Initialization of the switches giving the kind 
C     of electron involved (incoming, excited,
C     1st outgoing, 2nd outgoing) and the possible 
C     existence of a tip over the cluster
C
      INC=0
      EXC=0
      OUT1=0
      OUT2=0
      TIP=0
      CRYS=0
C
C
C..........   Reading of the input data in unit ICOM   ..........
C
C
      READ(ICOM,1) RIEN
      READ(ICOM,2) TEXTE1
      READ(ICOM,1) RIEN
C
C  Checking for the crystal structure if generated internally
C      (not used if potential generated by phagen)
C
      DO JLINE=1,NLINE
        READ(ICOM,49,END=619) STRING
        IF(STRING.EQ.'CRYS') THEN
          CRYS=1
          BACKSPACE ICOM
          BACKSPACE ICOM
          GOTO 612
        ENDIF
      ENDDO
  619 REWIND ICOM
      TEXTE2(1)='    '
      TEXTE2(2)='CRYS'
      TEXTE2(3)='TAL '
      TEXTE2(4)='STRU'
      TEXTE2(5)='CTUR'
      TEXTE2(6)='E (E'
      TEXTE2(7)='XTER'
      TEXTE2(8)='NAL)'
      TEXTE2(9)=' :  '
      TEXTE2(10)='    '
      TEXTE2(11)='    '
      TEXTE2(12)='    '
C
      DO JLINE=1,NLINE
        READ(ICOM,49) STRING
        IF(STRING.EQ.'TYPE') THEN
          BACKSPACE ICOM
          GOTO 600
        ENDIF
      ENDDO 
C
C================================================================
C                     Crystal Structure
C================================================================
C
  612 READ(ICOM,1) RIEN
      READ(ICOM,2) TEXTE2
      READ(ICOM,1) RIEN
C
      READ(ICOM,3) CRIST,CENTR,IBAS,NAT
      READ(ICOM,4) A,BSURA,CSURA,UNIT
C
      IF(IBAS.EQ.0) THEN
        DO JLINE=1,NLINE
          READ(ICOM,5) TESLEC
          IF(TESLEC.EQ.'SPECTRO') THEN
            BACKSPACE ICOM
            BACKSPACE ICOM
            BACKSPACE ICOM
            GOTO 600
          ENDIF
        ENDDO
      ENDIF
C
      READ(ICOM,6) ALPHAD,BETAD,GAMMAD
      READ(ICOM,7) IH,IK,II,IL
      READ(ICOM,8) NIV,COUPUR,ITEST,IESURF 
      IF(NAT.GT.1) THEN
        DO I=1,NAT
          J=3*(I-1)
          READ(ICOM,9) ATBAS(1+J),ATBAS(2+J),ATBAS(3+J),CHEM(I),
     1                 NZAT(I)
        ENDDO
      ELSE
        READ(ICOM,9) X1,Y1,Z1,CHEM(1),NZA
      ENDIF 
C
      READ(ICOM,5) TESLEC
      IF(TESLEC.EQ.'VECBAS ') THEN
        BACKSPACE ICOM
      ELSE
        IRET=10
        GOTO 605
      ENDIF
C
      DO I=1,8
        IF(CRIST.EQ.CODRES(I)) NCRIST=I
        IF(I.NE.8) THEN
          IF(CENTR.EQ.CODCTR(I)) NCENTR=I
        ENDIF
      ENDDO
      IF((NCRIST.EQ.0).OR.(NCENTR.EQ.0)) THEN
        IRET=1
        GOTO 605
      ENDIF
C
      IF(NCRIST.EQ.8) THEN
        DO I=1,3
          J=3*(I-1)
          IVN(I)=1
          READ(ICOM,9) VECBAS(1+J),VECBAS(2+J),VECBAS(3+J)
          IF(ABS(VECBAS(1+J)).LT.0.0001) THEN
            IF(ABS(VECBAS(2+J)).LT.0.0001) THEN
              IF(ABS(VECBAS(3+J)).LT.0.0001) THEN
                IVG0=IVG0+1
                IVN(I)=0
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSE
        READ(ICOM,9) X3,Y3,Z3
        READ(ICOM,9) X4,Y4,Z4
        READ(ICOM,9) X5,Y5,Z5
      ENDIF
      READ(ICOM,10) IREL,NREL,(PCREL(I),I=1,2)
      IF(IREL.EQ.1) THEN
        IF(NREL.GT.2) THEN
          NLIGNE=INT(FLOAT(NREL-2)/4.)+1
          DO J=1,NLIGNE
            READ(ICOM,11) (PCREL(I),I=1,4)
          ENDDO
        ENDIF
        IF(NREL.GT.10) THEN
          IRET=4
          GOTO 605
        ENDIF
      ELSEIF(IREL.EQ.0) THEN
        NREL=0
      ENDIF
      IF(NREL.EQ.0) THEN
        DO JREL=1,10
          PCREL(JREL)=0.
        ENDDO
      ENDIF
      READ(ICOM,12) OMEGAD1,OMEGAD2,IADS
C
C================================================================
C                 Type of calculation
C================================================================
C
      READ(ICOM,1) RIEN
 600  READ(ICOM,2) TEXTE3
      READ(ICOM,1) RIEN
C
      READ(ICOM,13) SPECTRO,ISPIN,IDICHR,I_AMP
C
      IF((IDICHR.EQ.2).AND.(ISPIN.EQ.0)) THEN
        PRINT 514
        STOP
      ENDIF
C
      IF(ISPIN.EQ.0) THEN
        NSPIN2=1
        NSPIN=1
      ELSEIF(ISPIN.EQ.1) THEN
        NSPIN2=4
        NSPIN=2
      ENDIF
C
C  Setting various spectroscopy-dependent parameters
C
C  'INTERACT' describes the interaction process. It takes the 
C             following values
C
C             = 'NOINTER'  :  no excitation of electron        (ex. LEED)
C             = 'PHOTONS'  :  excitation with photons          (ex. PED)
C             = 'COULOMB'  :  Coulomb excitation               (ex. AED)
C             = 'PHOCOUL'  :  photon + Coulomb - 2 electrons-  (ex. APECS)
C             = 'PHOPHOT'  :  photon for incoming and outgoing (ex. REXS)
C
      IF(SPECTRO.EQ.'PHD') THEN
        INTERACT='PHOTONS'
        ISPEC(1)=1
        NEL=1
        NCL=1
        OUT1=1
      ELSEIF(SPECTRO.EQ.'LED') THEN
        INTERACT='NOINTER'
        ISPEC(1)=2
        NEL=1
        NCL=1
        OUT1=1
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        INTERACT='PHOTONS'
        IF(IDICHR.GT.1) THEN
          PRINT 512
          STOP
        ENDIF
        ISPEC(1)=3
        NEL=1
        NCL=1
        EXC=1
      ELSEIF(SPECTRO.EQ.'AED') THEN
        INTERACT='COULOMB'
        ISPEC(1)=4
        NEL=1
        NCL=1
        OUT1=1
      ELSEIF(SPECTRO.EQ.'APC') THEN
        INTERACT='PHOCOUL'
        ISPEC(1)=1
        ISPEC(2)=4
        NEL=2
        NCL=1
        OUT1=1
        OUT2=1
      ELSEIF(SPECTRO.EQ.'RES') THEN
        INTERACT_I='PHOTONS'
        INTERACT_O='PHOTONS'
        INTERACT='PHOPHOT'
        ISPEC(1)=6
        NEL=1
        NCL=1
        EXC=1
      ELSEIF(SPECTRO.EQ.'ELS') THEN
        INTERACT='COULOMB'
        ISPEC(1)=7
        NEL=1
        NCL=1
        INC=1
        EXC=1
        OUT1=1
      ELSEIF(SPECTRO.EQ.'STM') THEN
        INTERACT='NOINTER'
        ISPEC(1)=8
        NEL=1
        NCL=2
        INC=1
        EXC=1
        OUT1=1
        TIP=1
      ELSEIF(SPECTRO.EQ.'EIG') THEN
        INTERACT='NOINTER'
        ISPEC(1)=10
        NEL=1
        NCL=1
        EXC=1
      ENDIF
C
C  Checking for the experimental parameters to read
C  The name of the spectroscopy is used as an anchor point
C
      DO JEL=1,NEL
        DO JLINE=1,100
          READ(ICOM,49) STRING
          IF(STRING.EQ.SPECNAME(ISPEC(JEL))) THEN
            GOTO 602
          ENDIF
        ENDDO
  602   BACKSPACE ICOM
        BACKSPACE ICOM
C
C================================================================
C                   Experimental parameters
C================================================================
C
C     Photoelectron experimental parameters 
C     (Photoelectron diffraction and APECS)
C
        IF(ISPEC(JEL).EQ.1) THEN
          READ(ICOM,1) RIEN
          READ(ICOM,2) TEXTE4
          READ(ICOM,1) RIEN
C
          READ(ICOM,20) NI,NLI,S_O,INITL,I_SO
C
          IF((NLI.EQ.'s').OR.(NLI.EQ.'S')) THEN
            LI=0
          ELSEIF((NLI.EQ.'p').OR.(NLI.EQ.'P')) THEN
            LI=1
          ELSEIF((NLI.EQ.'d').OR.(NLI.EQ.'D')) THEN
            LI=2
          ELSEIF((NLI.EQ.'f').OR.(NLI.EQ.'F')) THEN
            LI=3
          ELSEIF((NLI.EQ.'g').OR.(NLI.EQ.'G')) THEN
            LI=4
          ELSE
            IRET=5
            GOTO 605
          ENDIF
          IF(LI.GT.LI_M) THEN
            IRET=6
            GOTO 605
          ENDIF
          IF(I_SO.LE.0) THEN
            S_O='   '
          ELSEIF(I_SO.EQ.1) THEN
            IF(S_O.EQ.'1/2') THEN
              IF(LI.GT.1) IRET=7
            ELSEIF(S_O.EQ.'3/2') THEN
              IF((LI.LT.1).OR.(LI.GT.2)) IRET=7
            ELSEIF(S_O.EQ.'5/2') THEN
              IF((LI.LT.2).OR.(LI.GT.3)) IRET=7
            ELSEIF(S_O.EQ.'7/2') THEN
              IF((LI.LT.3).OR.(LI.GT.4)) IRET=7
            ELSEIF(S_O.EQ.'9/2') THEN
              IF(LI.NE.4) IRET=7
            ENDIF
          ELSEIF(I_SO.EQ.2) THEN
            S_O='   '
          ENDIF
C
          READ(ICOM,14) IPHI,NPHI,PHI0,PHI1
          READ(ICOM,14) ITHETA,NTHETA,THETA0,THETA1
          READ(ICOM,14) IE,NE,E0,E1
          READ(ICOM,14) IFTHET,NFTHET,R0,R1
          READ(ICOM,17) THLUM,PHILUM,ELUM,IPOL
          READ(ICOM,18) IMOD,IMOY,ACCEPT,ICHKDIR
C
C     LEED experimental parameters
C
        ELSEIF(ISPEC(JEL).EQ.2) THEN
          READ(ICOM,1) RIEN
          READ(ICOM,2) TEXTE4
          READ(ICOM,1) RIEN
C
          READ(ICOM,14) IPHI,NPHI,PHI0,PHI1
          READ(ICOM,14) ITHETA,NTHETA,THETA0,THETA1
          READ(ICOM,14) IE,NE,E0,E1
          READ(ICOM,14) IFTHET,NFTHET,R0,R1
          READ(ICOM,17) TH_INI,PHI_INI
          READ(ICOM,18) IMOD,IMOY,ACCEPT,ICHKDIR
C
          THLUM=TH_INI
          PHILUM=PHI_INI
          ELUM=0.
          IDICHR=0
          INITL=0
C
C
C     EXAFS experimental parameters
C
        ELSEIF(ISPEC(JEL).EQ.3) THEN
          READ(ICOM,1) RIEN
          READ(ICOM,2) TEXTE5
          READ(ICOM,1) RIEN
C
          READ(ICOM,39) EDGE,NEDGE,INITL
          READ(ICOM,43) THLUM,PHILUM,IPOL
          READ(ICOM,19) NE_X,EK_INI,EK_FIN,EPH_INI
C
          LI=NEDGE/2
          IF(NEDGE.GT.1) I_SO=1 
C
C     Auger electron experimental parameters
C
        ELSEIF(ISPEC(JEL).EQ.4) THEN
C
          READ(ICOM,1) RIEN
          READ(ICOM,2) TEXTE6
          READ(ICOM,1) RIEN
C
          READ(ICOM,40) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A
          READ(ICOM,42) I_MULT,IM1,MULT,IM2
          READ(ICOM,14) IPHI_A,NPHI_A,PHI0_A,PHI1_A
          READ(ICOM,14) ITHETA_A,NTHETA_A,THETA0_A,THETA1_A
          READ(ICOM,14) IE,NE,E0,E1
          READ(ICOM,14) IFTHET_A,NFTHET_A,R0_A,R1_A
          READ(ICOM,14) I_INT
          READ(ICOM,18) IMOD_A,IMOY_A,ACCEPT_A,ICHKDIR_A
C
          LI_C=NEDGE_C/2
          LI_I=NEDGE_I/2
          LI_A=NEDGE_A/2
C
          IF((EDGE_I.EQ.EDGE_A).AND.(LI_I.EQ.LI_A)) THEN
            I_SHELL=1
          ELSE
            I_SHELL=0
          ENDIF
C
          IE_A=0
          NE_A=1
          I_CP_A=0
C
          IF(EDGE_C.EQ.'K') THEN
            AUGER=' '//EDGE_C//EDGE_I//CHAR(48+NEDGE_I)//
     1            EDGE_A//CHAR(48+NEDGE_A)
          ELSE
            AUGER=EDGE_C//CHAR(48+NEDGE_C)//EDGE_I//
     1            CHAR(48+NEDGE_I)//EDGE_A//CHAR(48+NEDGE_A)
          ENDIF
          AUGER1=AUGER
C
C     Resonant elastic X-ray spectroscopy parameters
C
        ELSEIF(ISPEC(JEL).EQ.6) THEN
C
          READ(ICOM,1) RIEN
          READ(ICOM,2) TEXTE8
          READ(ICOM,1) RIEN
C
          READ(ICOM,44) EDGE,NEDGE,IE_RX,IPHI_RX,ITHETA_RX
          READ(ICOM,43) THLUM_I,PHILUM_I,IPOL_I,INITL_I
          READ(ICOM,43) THLUM_O,PHILUM_O,IPOL_O,INITL_O
          READ(ICOM,46) NE_RX,E_INI_RX,E_FIN_RX
          READ(ICOM,47) NPHI_RX,PHI_FIN_RX
          READ(ICOM,47) NTHETA_RX,THETA_FIN_RX
          READ(ICOM,48) G_MAX,G_HOL,G_SLO,G_CEN
C
C     Eigenvalue calculation parameters
C
        ELSEIF(ISPEC(JEL).EQ.10) THEN
C
          READ(ICOM,1) RIEN
          READ(ICOM,2) TEXTE10
          READ(ICOM,1) RIEN
C
          READ(ICOM,53) NE_EIG,EIG_INI,EIG_FIN,I_DAMP
C
          NE=NE_EIG
          N_LINE_E=INT((FLOAT(NE_EIG)-0.0001)/4.)+1
          N_LAST=4-(4*N_LINE_E-NE_EIG)
C
          IF(N_LINE_E.GT.1) THEN
            DO JLINE=1,N_LINE_E-1
               J=(JLINE-1)*4
               READ(ICOM,7) I_SPECTRUM(J+1),I_SPECTRUM(J+2),
     1                      I_SPECTRUM(J+3),I_SPECTRUM(J+4)
            ENDDO
          ENDIF
C
          J=4*(N_LINE_E-1)
C
          READ(ICOM,7) (I_SPECTRUM(J+K), K=1,N_LAST)
C
          READ(ICOM,55) I_PWM,METHOD,RACC,EXPO
          READ(ICOM,56) N_MAX,N_ITER,N_TABLE,SHIFT
          READ(ICOM,57) I_XN,I_VA,I_GN,I_WN
          READ(ICOM,58) LEVIN,ALPHAR,BETAR
C
          ACC=DBLE(RACC)
          IF(ABS(I_PWM).LE.2) THEN
            I_ACC=0
            N_ITER=N_MAX
          ELSEIF(I_PWM.EQ.3) THEN
            I_ACC=1
            N_ITER=N_MAX
          ELSEIF(I_PWM.EQ.-3) THEN
            I_ACC=-1
            N_ITER=N_MAX
          ELSEIF(I_PWM.EQ.4) THEN
            I_ACC=2
          ELSEIF(I_PWM.EQ.-4) THEN
            I_ACC=-2
          ENDIF
          IF(N_MAX.LT.N_ITER) N_ITER=N_MAX
C
          ALPHA=DCMPLX(ALPHAR)
          BETA=DCMPLX(BETAR)
C 
        ENDIF
C
        REWIND ICOM
C
C  End of the reading of the experimental parameters
C
      ENDDO
C
C  Non Auger electron case
C
      IF(ISPEC(1).LE.2) THEN
        IF(IPHI.EQ.-1) THEN
          IPHI=1
          I_EXT=0
          ICHKDIR=0
          STEREO='YES'
          IF(ABS(PHI1-PHI0).LT.0.0001) THEN
           PHI0=0.
            PHI1=360.
            NPHI=361
          ENDIF
          IF(ABS(THETA1-THETA0).LT.0.0001) THEN
            THETA0=0.
            THETA1=88.
            NTHETA=89
          ENDIF
        ELSEIF(IPHI.EQ.2) THEN
          IPHI=1
          I_EXT=1
        ELSEIF(IPHI.EQ.3) THEN
          IPHI=1
          I_EXT=-1
        ELSEIF(ITHETA.EQ.2) THEN
          ITHETA=1
          I_EXT=1
        ELSEIF(ITHETA.EQ.3) THEN
          ITHETA=1
          I_EXT=-1
        ELSEIF(IE.EQ.2) THEN
          IE=1
          I_EXT=1
        ELSEIF(IE.EQ.3) THEN
          IE=1
          I_EXT=-1
        ELSEIF(IE.EQ.4) THEN
          IF(SPECTRO.EQ.'PHD') THEN
            IE=1
            I_EXT=2
            IMOD=0
          ELSE
            IE=1
            I_EXT=1
          ENDIF
        ENDIF
      ENDIF
C
      ICALC=IPHI*IE+IPHI*ITHETA+IE*ITHETA
      IF((ICALC.NE.0).AND.(IFTHET.EQ.0)) IRET=3
C
C  When the direction of the analyzer might be experimentally 
C    inaccurate, the calculation will be done for nine
C    direction across the one given in the data file 
C    with an increment of one degree.
C  
      IF(ICHKDIR.EQ.1) THEN
        IF((ITHETA.EQ.1).AND.(IPHI.EQ.0)) THEN
          NPHI=9
          PHI0=PHI0-4.
          PHI1=PHI0+8.
        ELSEIF((IPHI.EQ.1).AND.(ITHETA.EQ.0)) THEN
          NTHETA=9
          THETA0=THETA0-4.
          THETA1=THETA0+8.
        ENDIF
      ENDIF 
C  
C  Initialization of the values for the scanned angle 
C                 and the "fixed" one
C  
      IF(IPHI.EQ.1) THEN
        N_FIXED=NTHETA
        N_SCAN=NPHI
        FIX0=THETA0
        FIX1=THETA1
        SCAN0=PHI0
        SCAN1=PHI1
        IPH_1=0
      ELSEIF(ITHETA.EQ.1) THEN
        N_FIXED=NPHI
        N_SCAN=NTHETA
        FIX0=PHI0
        FIX1=PHI1
        SCAN0=THETA0
        SCAN1=THETA1
        IPH_1=1
      ELSEIF(IE.EQ.1) THEN
        IF(NTHETA.GE.NPHI) THEN
          N_FIXED=NPHI
          N_SCAN=NTHETA
          FIX0=PHI0
          FIX1=PHI1
          SCAN0=THETA0
          SCAN1=THETA1
          IPH_1=1
        ELSE
          N_FIXED=NTHETA
          N_SCAN=NPHI
          FIX0=THETA0
          FIX1=THETA1
          SCAN0=PHI0
          SCAN1=PHI1
          IPH_1=0
        ENDIF
      ENDIF
C
C  Auger electron case
C
      IF((SPECTRO.EQ.'AED').OR.(SPECTRO.EQ.'APC')) THEN
        IF(IPHI_A.EQ.-1) THEN
          IPHI_A=1
          I_EXT_A=0
          ICHKDIR_A=0
          STEREO='YES'
          IF(ABS(PHI1_A-PHI0_A).LT.0.0001) THEN
            PHI0_A=0.
            PHI1_A=360.
            NPHI_A=361
          ENDIF
          IF(ABS(THETA1_A-THETA0_A).LT.0.0001) THEN
            THETA0_A=0.
            THETA1_A=88.
            NTHETA_A=89
          ENDIF
        ELSEIF(IPHI_A.EQ.2) THEN
          IPHI_A=1
          I_EXT_A=1
        ELSEIF(IPHI_A.EQ.3) THEN
          IPHI_A=1
          I_EXT_A=-1
        ELSEIF(ITHETA_A.EQ.2) THEN
          ITHETA_A=1
          I_EXT_A=1
        ELSEIF(ITHETA_A.EQ.3) THEN
          ITHETA_A=1
          I_EXT_A=-1
        ENDIF
C
C  Check for the consistency of the data for the two electrons in 
C     APECS, in particular when the sample is rotated (IMOD=1)
C
        IF(SPECTRO.EQ.'APC') THEN
          IF((LI_C.NE.LI).OR.(IMOD_A.NE.IMOD)) THEN
            IRET=11
            GOTO 605
          ENDIF
          DTH=THETA1-THETA0
          DTH_A=THETA1_A-THETA0_A
          DPH=PHI1-PHI0
          DPH_A=PHI1_A-PHI0_A
          IF((IMOD_A.EQ.1).AND.(IPHI_A.NE.IPHI)) IRET=13
          IF((IMOD_A.EQ.1).AND.(ITHETA_A.NE.ITHETA)) IRET=13
          IF((IMOD_A.EQ.1).AND.(NPHI_A.NE.NPHI)) IRET=13
          IF((IMOD_A.EQ.1).AND.(NTHETA_A.NE.NTHETA)) IRET=13
          IF(NPHI.GT.1) THEN
            IF((IMOD_A.EQ.1).AND.(DPH_A.NE.DPH)) IRET=13
          ENDIF
          IF(NTHETA.GT.1) THEN
            IF((IMOD_A.EQ.1).AND.(DTH_A.NE.DTH)) IRET=13
          ENDIF
        ENDIF
C
C  When the direction of the analyzer might be experimentally 
C    inaccurate, the calculation will be done for nine
C    direction across the one given in the data file 
C    with an increment of one degree.
C
        IF(ICHKDIR_A.EQ.1) THEN
          IF((ITHETA_A.EQ.1).AND.(IPHI_A.EQ.0)) THEN
            NPHI_A=9
            PHI0_A=PHI0_A-4.
            PHI1_A=PHI0_A+8.
          ELSEIF((IPHI_A.EQ.1).AND.(ITHETA_A.EQ.0)) THEN
            NTHETA_A=9
            THETA0_A=THETA0_A-4.
            THETA1_A=THETA0_A+8.
          ENDIF
        ENDIF 
C
C  Initialization of the values for the scanned angle 
C               and the "fixed" one
C
        IF(IPHI_A.EQ.1) THEN
          N_FIXED_A=NTHETA_A
          N_SCAN_A=NPHI_A
          FIX0_A=THETA0_A
          FIX1_A=THETA1_A
          SCAN0_A=PHI0_A
          SCAN1_A=PHI1_A
          IPH_1_A=0
        ELSEIF(ITHETA_A.EQ.1) THEN
          N_FIXED_A=NPHI_A
          N_SCAN_A=NTHETA_A
          FIX0_A=PHI0_A
          FIX1_A=PHI1_A
          SCAN0_A=THETA0_A
          SCAN1_A=THETA1_A
          IPH_1_A=1
        ENDIF
      ENDIF
C
      IF(SPECTRO.EQ.'XAS') THEN
        I_CP=1
        NE=NE_X
      ELSE
        I_CP=0
      ENDIF
C
C================================================================
C              Reading the calculation parameters
C         (cluster related or independent of energy)
C================================================================
C
      DO JLINE=1,NLINE
        READ(ICOM,50) STRING
        IF(STRING.EQ.'CLUS') THEN
          BACKSPACE ICOM
          BACKSPACE ICOM
          GOTO 603
        ENDIF
      ENDDO
C
  603 READ(ICOM,1) RIEN
      READ(ICOM,2) TEXTE7
      READ(ICOM,1) RIEN
C
      READ(ICOM,52) NAT,A,UNIT
C
      IF(UNIT.EQ.'ATU') THEN
        A=A*BOHR
      ENDIF
C
      READ(ICOM,23) NEMET
C
      BACKSPACE ICOM
      NLG=INT((NEMET-0.0001)/3) +1
      DO N=1,NLG
        NRL=3*N
        JD=3*(N-1)+1
        IF(N.EQ.NLG) NRL=NEMET
        READ(ICOM,24) NEMO,(IEMET(J), J=JD, NRL)
        IF(N.EQ.1) NEMET1=NEMO
      ENDDO
C
      READ(ICOM,25) ISOM,NONVOL(JFICH),NPATHP,VINT 
      READ(ICOM,30) IDWSPH,ISPEED,IATTS,IPRINT
C
      IF(IDWSPH.EQ.0) ISPEED=1
C
      READ(ICOM,31) IDCM,TD,T,RSJ
C
      NLEC=INT((NAT-0.0001)/4)+1
C
      READ(ICOM,14) ITEMP,NTEMP,TEMP0,TEMP1
C
      DO I=1,NLEC
        NDEB=4*(I-1) + 1
        NFIN=MIN0(4*I,NAT)
        READ(ICOM,33) (UJ2(J),J=NDEB,NFIN)
      ENDDO
C
C================================================================
C              Reading the calculation parameters
C                     (incoming electron)
C================================================================
C
      IF(INC.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,50) STRING
          IF(STRING.EQ.'INCO') THEN
            BACKSPACE ICOM
            BACKSPACE ICOM
            GOTO 607
          ENDIF
        ENDDO
  607   READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE7B
        READ(ICOM,1) RIEN
C
        IF(SPECTRO.EQ.'STM') THEN 
          NAT_OLD=NAT
          NAT=NAT_TI
        ENDIF
C
        READ(ICOM,21) NO_IN,NDIF_IN,ISPHER_IN,I_GR_IN
C
        IF(ISPHER_IN.EQ.0) THEN 
          IDWSPH=0
          NO_IN=0
        ENDIF
        IF(NO_IN.LT.0) NO_IN=8
        NUMAX_IN(1)=NO_IN/2
C
        READ(ICOM,22) ISFLIP_IN,IR_DIA_IN,ITRTL_IN,I_TEST_IN
C
C
        READ(ICOM,26) IFWD_IN,NTHOUT_IN,I_NO_IN,I_RA_IN
C
        IF(NTHOUT_IN.EQ.NDIF_IN-1) IFWD_IN=0
C
        IF(I_RA_IN.EQ.1) NO_IN=0
        DO JAT=1,NAT
          READ(ICOM,27) N_RA_IN(JAT),THFWD_IN(JAT),
     1                  IBWD_IN(JAT),THBWD_IN(JAT)
          IF(I_RA_IN.EQ.0) THEN
            N_RA_IN(JAT)=NO_IN
            NUMAX_IN(JAT)=NO_IN/2
          ELSEIF(I_RA_IN.EQ.1) THEN
            NUMAX_IN(JAT)=N_RA_IN(JAT)/2
            NO_IN=MAX(N_RA_IN(JAT),NO_IN)
          ENDIF
        ENDDO
C
        READ(ICOM,5) TESLEC
        IF(TESLEC.EQ.'IPW,NCU') THEN
          BACKSPACE ICOM
        ELSE
          IRET=8
          GOTO 605
        ENDIF
C
        READ(ICOM,28) IPW_IN,NCUT_IN,PCTINT_IN,IPP_IN
        READ(ICOM,29) ILENGTH_IN,RLENGTH_IN,UNLENGTH_IN
        READ(ICOM,32) ILPM_IN,XLPM0_IN
C
        IF((IDCM.EQ.1).OR.(ILPM_IN.EQ.1)) THEN
          CALL ATDATA
        ENDIF
        IF(SPECTRO.EQ.'STM') THEN 
          NAT=NAT_OLD
        ENDIF
      ENDIF
C
C================================================================
C              Reading the calculation parameters
C                     (excited electron)
C================================================================
C
      IF(EXC.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,50) STRING
          IF(STRING.EQ.'EXCI') THEN
            BACKSPACE ICOM
            BACKSPACE ICOM
            GOTO 608
          ENDIF
        ENDDO
  608   READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE7C
        READ(ICOM,1) RIEN
C
        READ(ICOM,21) NO_EX,NDIF_EX,ISPHER_EX,I_GR_EX
C
        IF(ISPHER_EX.EQ.0) THEN 
          IDWSPH=0
          NO_EX=0
        ENDIF
        IF(NO_EX.LT.0) NO_EX=8
        NUMAX_EX(1)=NO_EX/2
C
        READ(ICOM,22) ISFLIP_EX,IR_DIA_EX,ITRTL_EX,I_TEST_EX
C
C
        READ(ICOM,26) IFWD_EX,NTHOUT_EX,I_NO_EX,I_RA_EX
C
        IF(NTHOUT_EX.EQ.NDIF_EX-1) IFWD_EX=0
C
        IF(I_RA_EX.EQ.1) NO_EX=0
        DO JAT=1,NAT
          READ(ICOM,27) N_RA_EX(JAT),THFWD_EX(JAT),
     1                  IBWD_EX(JAT),THBWD_EX(JAT)
          IF(I_RA_EX.EQ.0) THEN
            N_RA_EX(JAT)=NO_EX
            NUMAX_EX(JAT)=NO_EX/2
          ELSEIF(I_RA_EX.EQ.1) THEN
            NUMAX_EX(JAT)=N_RA_EX(JAT)/2
            NO_EX=MAX(N_RA_EX(JAT),NO_EX)
          ENDIF
        ENDDO
C
        READ(ICOM,5) TESLEC
        IF(TESLEC.EQ.'IPW,NCU') THEN
          BACKSPACE ICOM
        ELSE
          IRET=8
          GOTO 605
        ENDIF
C
        READ(ICOM,28) IPW_EX,NCUT_EX,PCTINT_EX,IPP_EX
        READ(ICOM,29) ILENGTH_EX,RLENGTH_EX,UNLENGTH_EX
        READ(ICOM,32) ILPM_EX,XLPM0_EX
C
        IF((IDCM.EQ.1).OR.(ILPM_EX.EQ.1)) THEN
          CALL ATDATA
        ENDIF
      ENDIF
C
C================================================================
C              Reading the calculation parameters
C                     (1st outgoing electron)
C================================================================
C
      IF(OUT1.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,50) STRING
          IF(STRING.EQ.'OUTG') THEN
            BACKSPACE ICOM
            BACKSPACE ICOM
            GOTO 609
          ENDIF
        ENDDO
  609   READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE7D
        READ(ICOM,1) RIEN
C
        READ(ICOM,21) NO_O1,NDIF_O1,ISPHER_O1,I_GR_O1
C
        IF(ISPHER_O1.EQ.0) THEN 
          IDWSPH=0
          NO_O1=0
        ENDIF
        IF(NO_O1.LT.0) NO_O1=8
        NUMAX_O1(1)=NO_O1/2
C
        READ(ICOM,22) ISFLIP_O1,IR_DIA_O1,ITRTL_O1,I_TEST_O1
C
C
        READ(ICOM,26) IFWD_O1,NTHOUT_O1,I_NO_O1,I_RA_O1
C
        IF(NTHOUT_O1.EQ.NDIF_O1-1) IFWD_O1=0
C
        IF(I_RA_O1.EQ.1) NO_O1=0
        DO JAT=1,NAT
          READ(ICOM,27) N_RA_O1(JAT),THFWD_O1(JAT),
     1                  IBWD_O1(JAT),THBWD_O1(JAT)
          IF(I_RA_O1.EQ.0) THEN
            N_RA_O1(JAT)=NO_O1
            NUMAX_O1(JAT)=NO_O1/2
          ELSEIF(I_RA_O1.EQ.1) THEN
            NUMAX_O1(JAT)=N_RA_O1(JAT)/2
            NO_O1=MAX(N_RA_O1(JAT),NO_O1)
          ENDIF
        ENDDO
C
        READ(ICOM,5) TESLEC
        IF(TESLEC.EQ.'IPW,NCU') THEN
          BACKSPACE ICOM
        ELSE
          IRET=8
          GOTO 605
        ENDIF
C
        READ(ICOM,28) IPW_O1,NCUT_O1,PCTINT_O1,IPP_O1
        READ(ICOM,29) ILENGTH_O1,RLENGTH_O1,UNLENGTH_O1
        READ(ICOM,32) ILPM_O1,XLPM0_O1
C
        IF((IDCM.EQ.1).OR.(ILPM_O1.EQ.1)) THEN
          CALL ATDATA
        ENDIF
      ENDIF
C
C================================================================
C              Reading the calculation parameters
C                     (2nd outgoing electron)
C================================================================
C
      IF(OUT2.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,50) STRING
          IF(STRING.EQ.'2nd ') THEN
            BACKSPACE ICOM
            BACKSPACE ICOM
            GOTO 610
          ENDIF
        ENDDO
  610   READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE7E
        READ(ICOM,1) RIEN
C
        READ(ICOM,21) NO_O2,NDIF_O2,ISPHER_O2,I_GR_O2
C
        IF(ISPHER_O2.EQ.0) THEN 
          IDWSPH=0
          NO_O2=0
        ENDIF
        IF(NO_O2.LT.0) NO_O2=8
        NUMAX_O2(1)=NO_O2/2
C
        READ(ICOM,22) ISFLIP_O2,IR_DIA_O2,ITRTL_O2,I_TEST_O2
C
C
        READ(ICOM,26) IFWD_O2,NTHOUT_O2,I_NO_O2,I_RA_O2
C
        IF(NTHOUT_O2.EQ.NDIF_O2-1) IFWD_O2=0
C
        IF(I_RA_O2.EQ.1) NO_O2=0
        DO JAT=1,NAT
          READ(ICOM,27) N_RA_O2(JAT),THFWD_O2(JAT),
     1                  IBWD_O2(JAT),THBWD_O2(JAT)
          IF(I_RA_O2.EQ.0) THEN
            N_RA_O2(JAT)=NO_O2
            NUMAX_O2(JAT)=NO_O2/2
          ELSEIF(I_RA_O2.EQ.1) THEN
            NUMAX_O2(JAT)=N_RA_O2(JAT)/2
            NO_O2=MAX(N_RA_O2(JAT),NO_O2)
          ENDIF
        ENDDO
C
        READ(ICOM,5) TESLEC
        IF(TESLEC.EQ.'IPW,NCU') THEN
          BACKSPACE ICOM
        ELSE
          IRET=8
          GOTO 605
        ENDIF
C
        READ(ICOM,28) IPW_O2,NCUT_O2,PCTINT_O2,IPP_O2
        READ(ICOM,29) ILENGTH_O2,RLENGTH_O2,UNLENGTH_O2
        READ(ICOM,32) ILPM_O2,XLPM0_O2
C
        IF((IDCM.EQ.1).OR.(ILPM_O2.EQ.1)) THEN
          CALL ATDATA
        ENDIF
      ENDIF
C
C================================================================
C              Reading the calculation parameters
C                     (tip related)
C================================================================
C
      IF(TIP.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,50) STRING
          IF(STRING.EQ.'TIP ') THEN
            BACKSPACE ICOM
            BACKSPACE ICOM
            GOTO 611
          ENDIF
        ENDDO
C
  611   READ(ICOM,1) RIEN
        READ(ICOM,2) TEXTE7F
        READ(ICOM,1) RIEN
C
        READ(ICOM,52) NAT_TI,A_TI,UNIT_TI
C
        IF(UNIT_TI.EQ.'ATU') THEN
          A_TI=A_TI*BOHR
        ENDIF
C
        READ(ICOM,23) NEMET_TI
C
        BACKSPACE ICOM
        NLG=INT((NEMET_TI-0.0001)/3) +1
        DO N=1,NLG
          NRL=3*N
          JD=3*(N-1)+1
          IF(N.EQ.NLG) NRL=NEMET_TI
          READ(ICOM,24) NEMO_TI,(IEMET_TI(J), J=JD, NRL)
          IF(N.EQ.1) NEMET1_TI=NEMO_TI
        ENDDO
C
        READ(ICOM,25) ISOM_TI,NONVOL_TI(JFICH),NPATHP_TI,VINT_TI 
        READ(ICOM,30) IDWSPH_TI,ISPEED_TI,IATTS_TI,IPRINT_TI
C
        IF(IDWSPH_TI.EQ.0) ISPEED_TI=1
C
        READ(ICOM,31) IDCM_TI,TD_TI,T_TI,RSJ_TI
C
        NLEC=INT((NAT-0.0001)/4)+1
C
        DO I=1,NLEC
          NDEB=4*(I-1) + 1
          NFIN=MIN0(4*I,NAT_TI)
          READ(ICOM,33) (UJ2_TI(J),J=NDEB,NFIN)
        ENDDO
      ENDIF
C
C  Reading the names of the input/output files
C
      REWIND ICOM
C
C================================================================
C                 Reading the cluster files
C================================================================
C
      DO JLINE=1,NLINE
        READ(ICOM,51) STRING
        IF(STRING.EQ.'SAMP') THEN
          READ(ICOM,1) RIEN
          READ(ICOM,1) RIEN
          READ(ICOM,1) RIEN
          READ(ICOM,34) INFILE2,IUI2
          READ(ICOM,34) INFILE3,IUI3
          GOTO 613
        ENDIF
      ENDDO
      IRET=9
      GOTO 605
C
C================================================================
C             Reading the incoming electron files
C================================================================
C
 613  IF(INC.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,51) STRING
          IF(STRING.EQ.'INCO') THEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,34) INFILE4,IUI4
            READ(ICOM,34) INFILE5,IUI5
            READ(ICOM,34) INFILE6,IUI6
            GOTO 614
          ENDIF
        ENDDO
        IRET=9
        GOTO 605
      ENDIF
C
C================================================================
C             Reading the excited electron files
C================================================================
C
 614  IF(EXC.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,51) STRING
          IF(STRING.EQ.'EXCI') THEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,34) INFILE7,IUI7
            READ(ICOM,34) INFILE8,IUI8
            GOTO 615
          ENDIF
        ENDDO
        IRET=9
        GOTO 605
      ENDIF
C
C================================================================
C             Reading the outgoing electron files
C================================================================
C
 615  IF(OUT1.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,51) STRING
          IF(STRING.EQ.'OUTG') THEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,34) INFILE9,IUI9
            READ(ICOM,34) INFILE10,IUI10
            READ(ICOM,34) INFILE11,IUI11
            GOTO 616
          ENDIF
        ENDDO
        IRET=9
        GOTO 605
      ENDIF
C
C================================================================
C             Reading the 2nd outgoing electron files
C================================================================
C
 616  IF(OUT2.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,51) STRING
          IF(STRING.EQ.'2nd ') THEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,34) INFILE12,IUI12
            READ(ICOM,34) INFILE13,IUI13
            READ(ICOM,34) INFILE14,IUI14
            GOTO 617
          ENDIF
        ENDDO
        IRET=9
        GOTO 605
      ENDIF
C
C================================================================
C                 Reading the tip files
C================================================================
C
 617  IF(TIP.EQ.1) THEN
        DO JLINE=1,NLINE
          READ(ICOM,51) STRING
          IF(STRING.EQ.'TIP ') THEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,34) INFILE15,IUI15
            READ(ICOM,34) INFILE16,IUI16
            GOTO 618
          ENDIF
        ENDDO
        IRET=9
        GOTO 605
      ENDIF
C
C================================================================
C                 Reading the temperature dependent uj2 file
C================================================================
C
      IF((ITEMP.EQ.1).AND.(IDCM.EQ.0)) THEN
        DO JLINE=1,NLINE
          READ(ICOM,51) STRING
          IF(STRING.EQ.'TEM ') THEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,1) RIEN
            READ(ICOM,34) INFILE17,IUI17
            GOTO 618
          ENDIF
        ENDDO
        IRET=9
        GOTO 605
      ENDIF
C
C================================================================
C                 Reading the output files
C================================================================
C
C
 618  DO JLINE=1,NLINE
        READ(ICOM,51) STRING
        IF(STRING.EQ.'PUT ') THEN
          READ(ICOM,1) RIEN
          READ(ICOM,1) RIEN
          READ(ICOM,1) RIEN
          READ(ICOM,34) OUTFILE1,IUO1
          READ(ICOM,34) OUTFILE2,IUO2
          READ(ICOM,34) OUTFILE3,IUO3
          READ(ICOM,34) OUTFILE4,IUO4
          GOTO 620
        ENDIF
      ENDDO
      IRET=9
      GOTO 605
C
C    Setting up the unit number of the scratch files
C
 620  IUSCR=MAX0(ICOM,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,
     1           IUI9,IUI10,IUI11,IUI12,IUI13,IUI14,IUI15,
     2           IUI16,IUI17,IUO1,IUO2,IUO3,IUO4)+1
      IUSCR2=IUSCR+1
C
C================================================================
C            Reading the sample adsorbate file if any
C================================================================
C
C
      IF(IADS.GE.1) THEN
        OPEN(UNIT=IUI3, FILE=INFILE3, STATUS='OLD')
        READ(IUI3,1) RIEN
        READ(IUI3,12) NATA,NADS1,NADS2,NADS3
        IF(NATA.EQ.1) THEN
          NADS2=0
          NADS3=0
        ELSEIF(NATA.EQ.2) THEN
          NADS3=0
        ENDIF
        READ(IUI3,35) (NZAT(I),I=NAT+1,NAT+NATA)
        READ(IUI3,36) (CHEM(I),I=NAT+1,NAT+NATA)
        READ(IUI3,37) (UJ2(NAT+J),J=1,NATA)
        READ(IUI3,38) NRELA,(PCRELA(I),I=1,NRELA)
        IF(NRELA.EQ.0) THEN
          DO JRELA=1,3
            PCRELA(JRELA)=0.
          ENDDO
        ENDIF
        NADS=NADS1+NADS2+NADS3
        DO JADS=1,NADS
          READ(IUI3,9) (ADS(I,JADS),I=1,3)
        ENDDO
        CLOSE(IUI3)
      ELSE
        NATA=0
        NRELA=0
      ENDIF
C
C================================================================
C            Reading the tip adsorbate file if any
C================================================================
C
C
      IF(IADS.GE.1) THEN
        OPEN(UNIT=IUI16, FILE=INFILE16, STATUS='OLD')
        READ(IUI16,1) RIEN
        READ(IUI16,12) NATA_TI,NADS1_TI,NADS2_TI,NADS3_TI
        IF(NATA_TI.EQ.1) THEN
          NADS2_TI=0
          NADS3_TI=0
        ELSEIF(NATA_TI.EQ.2) THEN
          NADS3_TI=0
        ENDIF
        READ(IUI16,35) (NZAT_TI(I),I=NAT+1,NAT+NATA)
        READ(IUI16,36) (CHEM_TI(I),I=NAT+1,NAT+NATA)
        READ(IUI16,37) (UJ2_TI(NAT+J),J=1,NATA)
        READ(IUI16,38) NRELA_TI,(PCRELA_TI(I),I=1,NRELA)
        IF(NRELA_TI.EQ.0) THEN
          DO JRELA=1,3
            PCRELA_TI(JRELA)=0.
          ENDDO
        ENDIF
        NADS_TI=NADS1_TI+NADS2_TI+NADS3_TI
        DO JADS=1,NADS_TI
          READ(IUI16,9) (ADS_TI(I,JADS),I=1,3)
        ENDDO
        CLOSE(IUI16)
      ELSE
        NATA_TI=0
        NRELA_TI=0
      ENDIF
C
      GOTO 601
C
C   Missing input data file(s)
C
 605  REWIND ICOM
      DO JLINE=1,NLINE
        READ(ICOM,5) TESLEC 
        IF(TESLEC.EQ.'CONTROL') THEN
          BACKSPACE ICOM
          READ(ICOM,34) OUTFILE1,IUO1
          GOTO 601
        ENDIF
      ENDDO
C
C    Opening of the control file for printing
C
 601  IF((JFICH.EQ.1).OR.(ISOM.EQ.0)) THEN
        OPEN(UNIT=IUO1, FILE=OUTFILE1, STATUS='UNKNOWN')
      ENDIF
      IF((NFICHLEC.GT.1).AND.(ISOM.NE.0)) THEN
         WRITE(IUO1,105) INDATA(JFICH)
      ENDIF
C
C    Error messages
C
      IF(IRET.EQ.1) RETURN 1
      IF(IRET.EQ.3) RETURN 3
      IF(IRET.EQ.4) RETURN 4
      IF(IRET.EQ.5) RETURN 5
      IF(IRET.EQ.6) RETURN 6
      IF(IRET.EQ.7) RETURN 7
      IF(IRET.EQ.8) RETURN 8
      IF(IRET.EQ.9) RETURN 9
      IF(IRET.EQ.10) RETURN 10
      IF(IRET.EQ.11) RETURN 11
      IF(IRET.EQ.12) RETURN 12
      IF(IRET.EQ.13) RETURN 13
C
C    Initializations for particular cases
C
      IF((SPECTRO.EQ.'AED').OR.(SPECTRO.EQ.'APC')) THEN
        I_TEST=I_TEST_O1
        I_TEST_A=I_TEST
      ELSEIF(SPECTRO.EQ.'PHD') THEN
        I_TEST=I_TEST_O1
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        I_TEST=I_TEST_EX
      ELSEIF(SPECTRO.EQ.'LED') THEN
        I_TEST=I_TEST_O1
      ELSEIF(SPECTRO.EQ.'RES') THEN
        I_TEST=I_TEST_EX
      ENDIF
C
      IF(I_TEST.EQ.1) THEN
        IF(INTERACT.EQ.'DIPOLAR') THEN
          LI=0
          IF(SPECTRO.NE.'RES') THEN
            INITL=1
            IPOL=1
          ELSE
            INITL_I=1
            INITL_O=1
            IPOL_I=1
            IPOL_O=1
          ENDIF
        ELSEIF(INTERACT.EQ.'COULOMB') THEN
          LI_C=0
          LI_I=0
        ENDIF
      ENDIF
C
      IF(I_TEST.EQ.2) THEN
        IF(ABS(IPOL).EQ.1) THEN
          THLUM=-90.
          PHILUM=0.
        ELSEIF(ABS(IPOL).EQ.2) THEN
          THLUM=0.
          PHILUM=0.
        ENDIF
        IMOD=0
        VINT=0.
        VINT_TI=0.
        A=1.
      ENDIF
C
      IF((NFICHLEC.EQ.1).OR.(IBAS.EQ.1)) ISOM=0
C
      IF(SPECTRO.EQ.'RES') IPOL=IPOL_I*IPOL_O
      IF((IPOL.EQ.0).AND.(IDICHR.GT.0)) THEN
        PRINT 513
        STOP
      ENDIF
C
C  Set up of the switch controlling external 
C    reading of the detector directions and 
C    averaging over them for an undetected electron
C
      IF(SPECTRO.EQ.'APC') THEN
        IF((I_EXT.EQ.-1).OR.(I_EXT_A.EQ.-1)) THEN
          IF(I_EXT*I_EXT_A.EQ.0) THEN
            WRITE(IUO1,523)
            I_EXT=-1
            I_EXT_A=-1
            OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')
            OPEN(UNIT=IUI14, FILE=INFILE14, STATUS='OLD')
            READ(IUI11,713) IDIR,NSET
            READ(IUI14,713) IDIR_A,NSET_A
            IF(IDIR.EQ.2) THEN
              IF(NSET.NE.NSET_A) WRITE(IUO1,524) NSET,NSET_A
              STOP
            ENDIF
          ENDIF
        ENDIF
        IF(I_INT.EQ.1) THEN
          I_EXT=2
        ELSEIF(I_INT.EQ.2) THEN
          I_EXT_A=2
        ELSEIF(I_INT.EQ.3) THEN
          I_EXT=2
          I_EXT_A=2
        ENDIF
      ENDIF
C
      IF(I_EXT.EQ.-1) THEN
        OPEN(UNIT=IUI6, FILE=INFILE11, STATUS='OLD')
        READ(IUI11,701) IDIR,I_SET,N_POINTS
        READ(IUI11,702) I_PH,N_FIXED,N_SCAN
        DO JS=1,I_SET
          READ(IUI11,703) TH_0(JS),PH_0(JS)
        ENDDO
        CLOSE(IUI11)
        IF(IDIR.NE.2) IRET=12
        IF(I_PH.NE.IPH_1) IPH_1=I_PH
        IF((SPECTRO.EQ.'PHD').OR.(SPECTRO.EQ.'APC')) THEN
          IF(I_PH.EQ.0) THEN
            NTHETA=N_FIXED
            NPHI=N_SCAN
          ELSE
            NTHETA=N_SCAN
            NPHI=N_FIXED
          ENDIF
          ICHKDIR=2
        ENDIF
      ENDIF
      IF(I_EXT.GE.1) THEN
        OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')
        READ(IUI11,701) IDIR,I_SET,N_POINTS
        CLOSE(IUI11)
        IF((IDIR.NE.1).AND.(I_EXT.EQ.2)) IRET=12
        N_FIXED=N_POINTS
        N_SCAN=1
        NTHETA=N_POINTS
        NPHI=1
      ENDIF
      IF(I_EXT_A.GE.1) THEN
        IF(SPECTRO.EQ.'APC') THEN
          OPEN(UNIT=IUI14, FILE=INFILE14, STATUS='OLD')
          READ(IUI14,701) IDIR_A,I_SET_A,N_POINTS_A
          CLOSE(IUI14)
        ELSE
          OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')
          READ(IUI11,701) IDIR_A,I_SET_A,N_POINTS_A
          CLOSE(IUI11)
        ENDIF
        IF((IDIR_A.NE.1).AND.(I_EXT_A.EQ.2)) IRET=12
        N_FIXED_A=N_POINTS_A
        N_SCAN_A=1
        NTHETA_A=N_POINTS_A
        NPHI_A=1
      ENDIF
C
      IF(I_EXT_A.EQ.-1) THEN
        IF(SPECTRO.EQ.'APC') THEN
          OPEN(UNIT=IUI14, FILE=INFILE14, STATUS='OLD')
          READ(IUI14,701) IDIR_A,I_SET_A,N_POINTS_A
          READ(IUI14,702) I_PH_A,N_FIXED_A,N_SCAN_A
          CLOSE(IUI14)
        ELSE
          OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')
          READ(IUI11,701) IDIR_A,I_SET_A,N_POINTS_A
          READ(IUI11,702) I_PH_A,N_FIXED_A,N_SCAN_A
          CLOSE(IUI11)
        ENDIF
        IF(IDIR_A.NE.2) IRET=12
        IF(I_PH_A.EQ.0) THEN
          NTHETA_A=N_FIXED_A
          NPHI_A=N_SCAN_A
        ELSE
          NTHETA_A=N_SCAN_A
          NPHI_A=N_FIXED_A
        ENDIF
        ICHKDIR_A=2
      ENDIF
C
C
C..........   Writing of the input data in unit IUO1   ..........
C
C
      WRITE(IUO1,100)
      WRITE(IUO1,101)
      WRITE(IUO1,101)
      WRITE(IUO1,102) TEXTE1
      WRITE(IUO1,101)
      WRITE(IUO1,101)
      WRITE(IUO1,203)
C
      IF(I_TEST.NE.2) THEN
        WRITE(IUO1,201) TEXTE2
      ELSE
        IF(ABS(IPOL).EQ.1) THEN
          WRITE(IUO1,525)
        ELSEIF(ABS(IPOL).EQ.2) THEN
          WRITE(IUO1,526)
        ENDIF
      ENDIF
C
      IF(NAT.GT.NATP_M) RETURN 2 
      IF(NEMET.GT.NEMET_M) RETURN 2
      IF(TIP.EQ.1) THEN
        IF(NAT_TI.GT.NATP_M) RETURN 2 
        IF(NEMET_TI.GT.NEMET_M) RETURN 2
      ENDIF
      IF(NE.GT.NE_M) RETURN 2
C
      IF(I_TEST_IN.EQ.2) GOTO 606
      IF(I_TEST_EX.EQ.2) GOTO 606
      IF(I_TEST_O1.EQ.2) GOTO 606
      IF(I_TEST_O2.EQ.2) GOTO 606
C
C================================================================
C              Writing the crystal structure
C                (if generated internally)
C================================================================
C
      IF(IBAS.EQ.1) THEN
        WRITE(IUO1,103) CRIST,CENTR,IBAS,NAT
        IF(NCRIST.EQ.1) THEN
          BSURA=1.
          CSURA=1.
          WRITE(IUO1,304) A
        ELSEIF((NCRIST.EQ.2).OR.(NCRIST.EQ.7).OR.(NCRIST.EQ.6)) THEN
          BSURA=1.
          WRITE(IUO1,404) A,CSURA
          IF((NCRIST.EQ.6).AND.(CSURA.EQ.1.)) THEN
            WRITE(IUO1,206) ALPHAD
          ELSEIF(NCRIST.EQ.4) THEN
            WRITE(IUO1,306) BETAD
          ENDIF
        ELSEIF((NCRIST.EQ.3).OR.(NCRIST.EQ.5).OR.(NCRIST.EQ.8)) THEN
          WRITE(IUO1,104) A,BSURA,CSURA
          IF(NCRIST.NE.3) THEN
            WRITE(IUO1,106) ALPHAD,BETAD,GAMMAD
          ENDIF
        ENDIF
        IF(NCRIST.EQ.7) THEN
          WRITE(IUO1,107) IH,IK,II,IL
        ELSE
          WRITE(IUO1,207) IH,IK,IL
        ENDIF
        WRITE(IUO1,108) NIV,COUPUR,ITEST,IESURF
        IF(NAT.GT.1) THEN
          DO I=1,NAT
            J=3*(I-1)
            WRITE(IUO1,109) ATBAS(1+J),ATBAS(2+J),ATBAS(3+J),CHEM(I),
     1                      NZAT(I)
          ENDDO
        ENDIF
        IF(NCRIST.EQ.8) THEN
          DO I=1,3
            J=3*(I-1)
            WRITE(IUO1,209) VECBAS(1+J),VECBAS(2+J),VECBAS(3+J)
          ENDDO
        ENDIF
        IF(IREL.GE.1) THEN
          WRITE(IUO1,110) IREL,NREL,(PCREL(I),I=1,2)
          IF(NREL.GT.2) THEN
            NLIGNE=INT(FLOAT(NREL-2)/4.)+1
            DO J=1,NLIGNE
              WRITE(IUO1,210) (PCREL(I),I=1,4)
            ENDDO
          ENDIF
          IF(NREL.GT.10) RETURN 4
          WRITE(IUO1,112) OMEGAD1,OMEGAD2,IADS
        ENDIF
        IF((IREL.EQ.0).AND.(IADS.EQ.1)) WRITE(IUO1,212) IADS
        IF(IADS.GE.1) THEN
          WRITE(IUO1,501) 
          DO JADS=1,NADS
            IF(JADS.LE.NADS1) THEN 
              IF(JADS.EQ.1) WRITE(IUO1,303) NAT+1
              WRITE(IUO1,309) (ADS(I,JADS),I=1,3)
            ELSEIF((JADS.GT.NADS1).AND.(JADS.LE.(NADS1+NADS2))) THEN 
              IF(JADS.EQ.(NADS1+1)) WRITE(IUO1,303) NAT+2
              WRITE(IUO1,309) (ADS(I,JADS),I=1,3)
            ELSEIF(JADS.GT.(NADS1+NADS2)) THEN
              IF(JADS.EQ.(NADS2+1)) WRITE(IUO1,303) NAT+3
              WRITE(IUO1,309) (ADS(I,JADS),I=1,3)
            ENDIF
          ENDDO
        ENDIF
        IF((IREL.GT.0).OR.(NRELA.GT.0)) WRITE(IUO1,502)
        IF(NRELA.GT.0) THEN
          WRITE(IUO1,311) (PCRELA(I),I=1,NRELA)
        ENDIF
        IF(IBAS.GT.0) THEN
          OMEGA1=OMEGAD1*PIS180
          OMEGA2=OMEGAD2*PIS180
        ENDIF
C
        IF(IREL.GT.0) THEN
          WRITE(IUO1,211) (PCREL(I),I=1,NREL)
        ENDIF
      ENDIF
C
      IF(IBAS.EQ.0) THEN
        WRITE(IUO1,204) A,IBAS
      ENDIF
C
C
      IF(I_AMP.EQ.1) WRITE(6,534)
C
 606  IF(SPECTRO.EQ.'APC') WRITE(IUO1,517)
C
C================================================================
C              Writing the experimental parameters
C                         (PhD case)
C================================================================
C
      IF(SPECTRO.EQ.'PHD') THEN
C
        IF(IPHI.EQ.1) THEN
          IF(STEREO.EQ.' NO') THEN 
            WRITE(IUO1,503)
          ELSE
            WRITE(IUO1,527)
          ENDIF
        ENDIF
        IF(IE.EQ.1) WRITE(IUO1,504)
        IF(ITHETA.EQ.1) WRITE(IUO1,505)
        IF(IFTHET.EQ.1) WRITE(IUO1,506)
C
        WRITE(IUO1,201) TEXTE4
        WRITE(IUO1,113) ISPIN,IDICHR
        WRITE(IUO1,120) NI,NLI,S_O,INITL,I_SO
C
        WRITE(IUO1,114) IPHI,NPHI,PHI0,PHI1
        WRITE(IUO1,115) ITHETA,NTHETA,THETA0,THETA1
        WRITE(IUO1,116) IE,NE,E0,E1
        WRITE(IUO1,216) IFTHET,NFTHET,R0,R1
C
        IF((ITHETA.EQ.1).AND.(IFTHET.EQ.0)) THEN
          IF((THETA0.LT.-90.0).OR.(THETA1.GT.90.0)) THEN
            WRITE(IUO1,508)
            STOP
          ENDIF 
          IF(ABS(THLUM).GT.90.0) THEN
            WRITE(IUO1,509)
            STOP
          ENDIF
        ENDIF 
C 
        WRITE(IUO1,117) THLUM,PHILUM,ELUM,IPOL
        WRITE(IUO1,118) IMOD,IMOY,ACCEPT,ICHKDIR
C
        IF(IMOY.GT.3) IMOY=3
        IF(IMOY.LT.0) IMOY=0
        IF(IMOY.EQ.0) NDIR=1
        IF(IMOY.EQ.1) NDIR=5
        IF(IMOY.EQ.2) NDIR=13
        IF(IMOY.EQ.3) NDIR=49
        IF((LI.EQ.0).AND.(INITL.NE.0)) INITL=1
C
C================================================================
C              Writing the experimental parameters
C                         (LEED case)
C================================================================
C
      ELSEIF(SPECTRO.EQ.'LED') THEN
C
        IF(IPHI.EQ.1) THEN
          IF(STEREO.EQ.' NO') THEN 
            WRITE(IUO1,529)
          ELSE
            WRITE(IUO1,530)
          ENDIF
        ENDIF
        IF(IE.EQ.1) WRITE(IUO1,531)
        IF(ITHETA.EQ.1) WRITE(IUO1,532)
        IF(IFTHET.EQ.1) WRITE(IUO1,506)
C
        WRITE(IUO1,201) TEXTE4
        WRITE(IUO1,141) ISPIN
C
        WRITE(IUO1,114) IPHI,NPHI,PHI0,PHI1
        WRITE(IUO1,115) ITHETA,NTHETA,THETA0,THETA1
        WRITE(IUO1,116) IE,NE,E0,E1
        WRITE(IUO1,216) IFTHET,NFTHET,R0,R1
C
        IF((ITHETA.EQ.1).AND.(IFTHET.EQ.0)) THEN
          IF((THETA0.LT.-90.0).OR.(THETA1.GT.90.0)) THEN
            WRITE(IUO1,508)
            STOP
          ENDIF 
        ENDIF 
C 
        WRITE(IUO1,142) TH_INI,PHI_INI
        WRITE(IUO1,118) IMOD,IMOY,ACCEPT,ICHKDIR
C
        IF(IMOY.GT.3) IMOY=3
        IF(IMOY.LT.0) IMOY=0
        IF(IMOY.EQ.0) NDIR=1
        IF(IMOY.EQ.1) NDIR=5
        IF(IMOY.EQ.2) NDIR=13
        IF(IMOY.EQ.3) NDIR=49
C
C================================================================
C              Writing the experimental parameters
C                         (XAS case)
C================================================================
C
      ELSEIF(SPECTRO.EQ.'XAS') THEN
C
        WRITE(IUO1,507)
        WRITE(IUO1,201) TEXTE5
        WRITE(IUO1,113) ISPIN,IDICHR
        WRITE(IUO1,111) IPOL
        WRITE(IUO1,134) EDGE,NEDGE,INITL,THLUM,PHILUM
        WRITE(IUO1,119) NE_X,EK_INI,EK_FIN,EPH_INI
C
C================================================================
C              Writing the experimental parameters
C                         (AED case)
C================================================================
C
      ELSEIF(SPECTRO.EQ.'AED') THEN
C
        IF(IPHI_A.EQ.1) THEN
          IF(STEREO.EQ.' NO') THEN 
            WRITE(IUO1,515)
          ELSE
            WRITE(IUO1,528)
          ENDIF
        ENDIF
        IF(ITHETA_A.EQ.1) WRITE(IUO1,516)
        WRITE(IUO1,201) TEXTE6
        WRITE(IUO1,113) ISPIN,IDICHR
        WRITE(IUO1,135) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A
        WRITE(IUO1,140) I_MULT,IM1,MULT,IM2
        WRITE(IUO1,136) IPHI_A,NPHI_A,PHI0_A,PHI1_A
        WRITE(IUO1,137) ITHETA_A,NTHETA_A,THETA0_A,THETA1_A
        WRITE(IUO1,138) IFTHET_A,NFTHET_A,R0_A,R1_A
        WRITE(IUO1,139) I_INT
        WRITE(IUO1,118) IMOD_A,IMOY_A,ACCEPT_A,ICHKDIR_A
C
        IF(IMOY_A.GT.3) IMOY_A=3
        IF(IMOY_A.LT.0) IMOY_A=0
        IF(IMOY_A.EQ.0) NDIR_A=1
        IF(IMOY_A.EQ.1) NDIR_A=5
        IF(IMOY_A.EQ.2) NDIR_A=13
        IF(IMOY_A.EQ.3) NDIR_A=49
C
C================================================================
C              Writing the experimental parameters
C                         (APECS case)
C================================================================
C
      ELSEIF(SPECTRO.EQ.'APC') THEN
C
        WRITE(IUO1,518)
        IF(IPHI.EQ.1) WRITE(IUO1,503)
        IF(ITHETA.EQ.1) WRITE(IUO1,505)
        IF(IFTHET.EQ.1) WRITE(IUO1,506)
C
        WRITE(IUO1,201) TEXTE4
        WRITE(IUO1,113) ISPIN,IDICHR
        WRITE(IUO1,120) NI,NLI,S_O,INITL,I_SO
C
        WRITE(IUO1,114) IPHI,NPHI,PHI0,PHI1
        WRITE(IUO1,115) ITHETA,NTHETA,THETA0,THETA1
        WRITE(IUO1,116) IE,NE,E0,E1
        WRITE(IUO1,216) IFTHET,NFTHET,R0,R1
C
        IF((ITHETA.EQ.1).AND.(IFTHET.EQ.0)) THEN
          IF((THETA0.LT.-90.0).OR.(THETA1.GT.90.0)) THEN
            WRITE(IUO1,508)
            STOP
          ENDIF
          IF(ABS(THLUM).GT.90.0) THEN
            WRITE(IUO1,509)
            STOP
          ENDIF
        ENDIF 
C
        WRITE(IUO1,117) THLUM,PHILUM,ELUM,IPOL
        WRITE(IUO1,118) IMOD,IMOY,ACCEPT,ICHKDIR
C
        IF(IMOY.GT.3) IMOY=3
        IF(IMOY.LT.0) IMOY=0
        IF(IMOY.EQ.0) NDIR=1
        IF(IMOY.EQ.1) NDIR=5
        IF(IMOY.EQ.2) NDIR=13
        IF(IMOY.EQ.3) NDIR=49
        IF((LI.EQ.0).AND.(INITL.NE.0)) INITL=1
C
        WRITE(IUO1,519)
        IF(IPHI_A.EQ.1) WRITE(IUO1,515)
        IF(ITHETA_A.EQ.1) WRITE(IUO1,516)
        WRITE(IUO1,201) TEXTE6
        WRITE(IUO1,113) ISPIN,IDICHR
        WRITE(IUO1,135) EDGE_C,NEDGE_C,EDGE_I,NEDGE_I,EDGE_A,NEDGE_A
        WRITE(IUO1,140) I_MULT,IM1,MULT,IM2
        WRITE(IUO1,136) IPHI_A,NPHI_A,PHI0_A,PHI1_A
        WRITE(IUO1,137) ITHETA_A,NTHETA_A,THETA0_A,THETA1_A
        WRITE(IUO1,138) IFTHET_A,NFTHET_A,R0_A,R1_A
        WRITE(IUO1,139) I_INT
        WRITE(IUO1,118) IMOD_A,IMOY_A,ACCEPT_A,ICHKDIR_A
C
        IF(IMOY_A.GT.3) IMOY_A=3
        IF(IMOY_A.LT.0) IMOY_A=0
        IF(IMOY_A.EQ.0) NDIR_A=1
        IF(IMOY_A.EQ.1) NDIR_A=5
        IF(IMOY_A.EQ.2) NDIR_A=13
        IF(IMOY_A.EQ.3) NDIR_A=49
C
        WRITE(IUO1,520)
C
C================================================================
C              Writing the experimental parameters
C                         (REXS case)
C================================================================
C
      ELSEIF(SPECTRO.EQ.'RES') THEN
C
        WRITE(IUO1,533)
        WRITE(IUO1,201) TEXTE8
        WRITE(IUO1,144) EDGE,NEDGE,IE_RX,IPHI_RX
        WRITE(IUO1,143) THLUM_I,PHILUM_I,IPOL_I,INITL_I
        WRITE(IUO1,145) THLUM_O,PHILUM_O,IPOL_O,INITL_O
        WRITE(IUO1,146) NE_RX,E_INI_RX,E_FIN_RX
        WRITE(IUO1,147) NPHI_RX,PHI_FIN_RX
        WRITE(IUO1,159) NTHETA_RX,THETA_FIN_RX
        WRITE(IUO1,148) G_MAX,G_HOL,G_SLO,G_CEN
C
C================================================================
C              Writing the experimental parameters
C                         (EELS case)
C================================================================
C
      ELSEIF(SPECTRO.EQ.'ELS') THEN
C
        WRITE(IUO1,535)
C
C================================================================
C              Writing the calculation parameters
C                      (Eigenvalue case)
C================================================================
C
      ELSEIF(SPECTRO.EQ.'EIG') THEN
C
        WRITE(IUO1,536)
        WRITE(IUO1,201) TEXTE10
        WRITE(IUO1,153) NE_EIG,EIG_INI,EIG_FIN,I_DAMP
        DO JLINE=1,N_LINE_E-1
          J=(JLINE-1)*4
          WRITE(IUO1,154) I_SPECTRUM(J+1),I_SPECTRUM(J+2),
     1                    I_SPECTRUM(J+3),I_SPECTRUM(J+4)
        ENDDO
        J=4*(N_LINE_E-1)
        WRITE(IUO1,154) (I_SPECTRUM(J+K),K=1,N_LAST)
C
        WRITE(IUO1,155) I_PWM,METHOD,RACC,EXPO
        WRITE(IUO1,156) N_MAX,N_ITER,N_TABLE,SHIFT
        WRITE(IUO1,157) I_XN,I_VA,I_GN,I_WN
        WRITE(IUO1,158) LEVIN,ALPHAR,BETAR
C
      ENDIF
C
C================================================================
C              Writing the calculation parameters
C                     (sample cluster case)
C================================================================
C
      WRITE(IUO1,201) TEXTE7

      IF(SPECTRO.NE.'EIG') THEN
        DO N=1,NLG
          NRL=3*N
          JD=3*(N-1)+1
          IF(N.EQ.NLG) NRL=NEMET
          IF(N.EQ.1) NEMO=NEMET1
          IF(N.LT.NLG) THEN
            WRITE(IUO1,123) NEMO,(IEMET(J), J=JD, NRL)
          ELSE
            NTE=NEMET-JD+1
            IF(NTE.EQ.1) WRITE(IUO1,223) NEMO,(IEMET(J),J=JD,NEMET)
            IF(NTE.EQ.2) WRITE(IUO1,323) NEMO,(IEMET(J),J=JD,NEMET)
            IF(NTE.EQ.3) WRITE(IUO1,123) NEMO,(IEMET(J),J=JD,NEMET)
          ENDIF
        ENDDO
        WRITE(IUO1,124) ISOM,NONVOL(JFICH),NPATHP,VINT
        WRITE(IUO1,129) IDWSPH,ISPEED,IATTS,IPRINT
      ELSE
        WRITE(IUO1,149) VINT
      ENDIF
      WRITE(IUO1,130) IDCM,TD,T,RSJ
      WRITE(IUO1,150) ITEMP,NTEMP,TEMP0,TEMP1
      DO I=1,NLEC
        NDEB=4*(I-1) + 1
        NFIN=4*I
        IF(I.EQ.NLEC) NFIN=NAT
        NUJ=NFIN-NDEB+1
        IF(NUJ.EQ.1) WRITE(IUO1,132) (UJ2(J),J=NDEB,NFIN)
        IF(NUJ.EQ.2) WRITE(IUO1,232) (UJ2(J),J=NDEB,NFIN)
        IF(NUJ.EQ.3) WRITE(IUO1,332) (UJ2(J),J=NDEB,NFIN)
        IF(NUJ.EQ.4) WRITE(IUO1,432) (UJ2(J),J=NDEB,NFIN)
      ENDDO
      IF(IADS.EQ.1) THEN
        IF(NATA.EQ.1) WRITE(IUO1,133) (UJ2(J),J=NAT+1,NAT+NATA)
        IF(NATA.EQ.2) WRITE(IUO1,233) (UJ2(J),J=NAT+1,NAT+NATA)
        IF(NATA.EQ.3) WRITE(IUO1,333) (UJ2(J),J=NAT+1,NAT+NATA)
      ENDIF
C
      DO J=1,NATM
        UJ2(J)=UJ2(J)/(A*A)
      ENDDO
C
C================================================================
C              Writing the calculation parameters
C                     (incoming electron case)
C================================================================
C
      IF(INC.EQ.1) THEN
        WRITE(IUO1,201) TEXTE7B
C
        WRITE(IUO1,121) NO_IN,NDIF_IN,ISPHER_IN,I_GR_IN
C
        IF(SPECTRO.EQ.'STM') THEN
          NAT_OLD=NAT
          A_OLD=A
          NAT=NAT_TI
          A=A_TI
        ENDIF
C
        WRITE(IUO1,122) ISFLIP_IN,IR_DIA_IN,ITRTL_IN,I_TEST_IN
C
        IF(ISFLIP_IN.EQ.0) THEN
          NSTEP_IN=3
        ELSE
          NSTEP_IN=1
        ENDIF
C
        WRITE(IUO1,125) IFWD_IN,NTHOUT_IN,I_NO_IN,I_RA_IN
        DO JAT=1,NAT
          WRITE(IUO1,126) N_RA_IN(JAT),THFWD_IN(JAT),IBWD_IN(JAT),
     1                    THBWD_IN(JAT)
          RTHFWD_IN(JAT)=THFWD_IN(JAT)*PIS180
          RTHBWD_IN(JAT)=THBWD_IN(JAT)*PIS180
        ENDDO
        WRITE(IUO1,127) IPW_IN,NCUT_IN,PCTINT_IN,IPP_IN
        WRITE(IUO1,128) ILENGTH_IN,RLENGTH_IN,UNLENGTH_IN
        WRITE(IUO1,131) ILPM_IN,XLPM0_IN
C
        IF(UNLENGTH_IN.EQ.'ATU') RLENGTH_IN=RLENGTH_IN*BOHR/A
        IF(UNLENGTH_IN.EQ.'ANG') RLENGTH_IN=RLENGTH_IN/A
C
        IF(SPECTRO.EQ.'STM') THEN
          NAT=NAT_OLD
          A=A_OLD
        ENDIF
      ENDIF
C
C================================================================
C              Writing the calculation parameters
C                     (excited electron case)
C================================================================
C
      IF(EXC.EQ.1) THEN
        WRITE(IUO1,201) TEXTE7C
        IF(SPECTRO.NE.'EIG') THEN
C
          WRITE(IUO1,121) NO_EX,NDIF_EX,ISPHER_EX,I_GR_EX
C
          IF(SPECTRO.EQ.'XAS') NDIF_EX=NDIF_EX+1
C
          WRITE(IUO1,122) ISFLIP_EX,IR_DIA_EX,ITRTL_EX,I_TEST_EX
C
          IF(ISFLIP_EX.EQ.0) THEN
            NSTEP_EX=3
          ELSE
            NSTEP_EX=1
          ENDIF
C
          WRITE(IUO1,125) IFWD_EX,NTHOUT_EX,I_NO_EX,I_RA_EX
          DO JAT=1,NAT
            WRITE(IUO1,126) N_RA_EX(JAT),THFWD_EX(JAT),IBWD_EX(JAT),
     1                      THBWD_EX(JAT)
            RTHFWD_EX(JAT)=THFWD_EX(JAT)*PIS180
            RTHBWD_EX(JAT)=THBWD_EX(JAT)*PIS180
          ENDDO
          WRITE(IUO1,127) IPW_EX,NCUT_EX,PCTINT_EX,IPP_EX
          WRITE(IUO1,128) ILENGTH_EX,RLENGTH_EX,UNLENGTH_EX
        ENDIF
        WRITE(IUO1,131) ILPM_EX,XLPM0_EX
C
        IF(SPECTRO.NE.'EIG') THEN
          IF(UNLENGTH_EX.EQ.'ATU') RLENGTH_EX=RLENGTH_EX*BOHR/A
          IF(UNLENGTH_EX.EQ.'ANG') RLENGTH_EX=RLENGTH_EX/A
        ENDIF
      ENDIF
C
C================================================================
C              Writing the calculation parameters
C                     (outgoing electron case)
C================================================================
C
      IF(OUT1.EQ.1) THEN
        WRITE(IUO1,201) TEXTE7D
C
        WRITE(IUO1,121) NO_O1,NDIF_O1,ISPHER_O1,I_GR_O1
C
        WRITE(IUO1,122) ISFLIP_O1,IR_DIA_O1,ITRTL_O1,I_TEST_O1
C
        IF(ISFLIP_O1.EQ.0) THEN
          NSTEP_O1=3
        ELSE
          NSTEP_O1=1
        ENDIF
C
        WRITE(IUO1,125) IFWD_O1,NTHOUT_O1,I_NO_O1,I_RA_O1
        DO JAT=1,NAT
          WRITE(IUO1,126) N_RA_O1(JAT),THFWD_O1(JAT),IBWD_O1(JAT),
     1                    THBWD_O1(JAT)
          RTHFWD_O1(JAT)=THFWD_O1(JAT)*PIS180
          RTHBWD_O1(JAT)=THBWD_O1(JAT)*PIS180
        ENDDO
        WRITE(IUO1,127) IPW_O1,NCUT_O1,PCTINT_O1,IPP_O1
        WRITE(IUO1,128) ILENGTH_O1,RLENGTH_O1,UNLENGTH_O1
        WRITE(IUO1,131) ILPM_O1,XLPM0_O1
C
        IF(UNLENGTH_O1.EQ.'ATU') RLENGTH_O1=RLENGTH_O1*BOHR/A
        IF(UNLENGTH_O1.EQ.'ANG') RLENGTH_O1=RLENGTH_O1/A
      ENDIF
C
C================================================================
C              Writing the calculation parameters
C                   (2nd outgoing electron case)
C================================================================
C
      IF(OUT2.EQ.1) THEN
        WRITE(IUO1,201) TEXTE7E
C
        WRITE(IUO1,121) NO_O2,NDIF_O2,ISPHER_O2,I_GR_O2
C
        WRITE(IUO1,122) ISFLIP_O2,IR_DIA_O2,ITRTL_O2,I_TEST_O2
C
        IF(ISFLIP_O2.EQ.0) THEN
          NSTEP_O2=3
        ELSE
          NSTEP_O2=1
        ENDIF
C
        WRITE(IUO1,125) IFWD_O2,NTHOUT_O2,I_NO_O2,I_RA_O2
        DO JAT=1,NAT
          WRITE(IUO1,126) N_RA_O2(JAT),THFWD_O2(JAT),IBWD_O2(JAT),
     1                    THBWD_O2(JAT)
          RTHFWD_O2(JAT)=THFWD_O2(JAT)*PIS180
          RTHBWD_O2(JAT)=THBWD_O2(JAT)*PIS180
        ENDDO
        WRITE(IUO1,127) IPW_O2,NCUT_O2,PCTINT_O2,IPP_O2
        WRITE(IUO1,128) ILENGTH_O2,RLENGTH_O2,UNLENGTH_O2
        WRITE(IUO1,131) ILPM_O2,XLPM0_O2
C
        IF(UNLENGTH_O2.EQ.'ATU') RLENGTH_O2=RLENGTH_O2*BOHR/A
        IF(UNLENGTH_O2.EQ.'ANG') RLENGTH_O2=RLENGTH_O2/A
      ENDIF
C
C================================================================
C              Writing the calculation parameters
C                     (tip cluster case)
C================================================================
C
      IF(TIP.EQ.1) THEN
        WRITE(IUO1,201) TEXTE7F
C
        DO N=1,NLG
          NRL=3*N
          JD=3*(N-1)+1
          IF(N.EQ.NLG) NRL=NEMET_TI
          IF(N.EQ.1) NEMO_TI=NEMET1_TI
          IF(N.LT.NLG) THEN
            WRITE(IUO1,123) NEMO_TI,(IEMET_TI(J), J=JD, NRL)
          ELSE
            NTE=NEMET_TI-JD+1
            IF(NTE.EQ.1) WRITE(IUO1,223) NEMO_TI,
     1                   (IEMET_TI(J),J=JD,NEMET)
            IF(NTE.EQ.2) WRITE(IUO1,323) NEMO_TI,
     1                   (IEMET_TI(J),J=JD,NEMET)
            IF(NTE.EQ.3) WRITE(IUO1,123) NEMO_TI,
     1                   (IEMET_TI(J),J=JD,NEMET)
          ENDIF
        ENDDO
        WRITE(IUO1,124) ISOM_TI,NONVOL_TI(JFICH),NPATHP_TI,VINT_TI
        WRITE(IUO1,129) IDWSPH_TI,ISPEED_TI,IATTS_TI,IPRINT_TI
        WRITE(IUO1,130) IDCM_TI,TD_TI,T_TI,RSJ_TI
        DO I=1,NLEC
          NDEB=4*(I-1) + 1
          NFIN=4*I
          IF(I.EQ.NLEC) NFIN=NAT
          NUJ=NFIN-NDEB+1
          IF(NUJ.EQ.1) WRITE(IUO1,132) (UJ2_TI(J),J=NDEB,NFIN)
          IF(NUJ.EQ.2) WRITE(IUO1,232) (UJ2_TI(J),J=NDEB,NFIN)
          IF(NUJ.EQ.3) WRITE(IUO1,332) (UJ2_TI(J),J=NDEB,NFIN)
          IF(NUJ.EQ.4) WRITE(IUO1,432) (UJ2_TI(J),J=NDEB,NFIN)
        ENDDO
        IF(IADS.EQ.1) THEN
          IF(NATA_TI.EQ.1) WRITE(IUO1,133) 
     1                     (UJ2_TI(J),J=NAT_TI+1,NAT_TI+NATA_TI)
          IF(NATA_TI.EQ.2) WRITE(IUO1,233) 
     1                     (UJ2_TI(J),J=NAT_TI+1,NAT_TI+NATA_TI)
          IF(NATA_TI.EQ.3) WRITE(IUO1,333) 
     1                     (UJ2_TI(J),J=NAT_TI+1,NAT_TI+NATA_TI)
        ENDIF
C
        DO J=1,NATM
          UJ2_TI(J)=UJ2_TI(J)/(A_TI*A_TI)
        ENDDO
      ENDIF
C
C  Setting up of various parameters
C
      IF(SPECTRO.EQ.'EIG') THEN
C
C  Switch for including vibrational damping into the MS matrix
C
C           I_VIB = 0 : no vibrations included
C           I_VIB = 1 : vibrations included
C
C         and mean free path-like damping
C
C           I_MFP = 0 : no Im(k) damping included
C           I_MFP = 1 : Im(k) damping included
C
        I_VIB=MOD(I_DAMP,2)
        IF(I_VIB.EQ.1) THEN
          IDWSPH=1
        ELSE
          IDWSPH=0
        ENDIF
        IF(I_DAMP.LE.1) THEN
          I_MFP=0
        ELSE
          I_MFP=1
        ENDIF
      ENDIF
C
      I_SPHER=MAX(ISPHER_IN,ISPHER_EX,ISPHER_O1,ISPHER_O2)
      I_DWSPH=MAX(IDWSPH,IDWSPH_TI)
C
      QD=0.
      QD_TI=0.
      IF(E0.EQ.0.) E0=0.0001
      NPOINT=NPHI*NE*NTHETA
      ISORT1=0
      IF(NPOINT.GT.250) THEN
        ISORT1=1
        WRITE(IUO1,510)
      ENDIF
C
      IF(I_DWSPH.EQ.1) THEN
        NFAC=N_GAUNT
      ELSE
        NFAC=4*NL_M
      ENDIF
C
C  Storage of the logarithm of the Gamma function GLD(N+1,N_INT)
C  for integer (N_INT=1) and semi-integer (N_INT=2) values :
C
C    GLD(N+1,1)   =   Log(N!) for N integer
C    GLD(N+1/2,2) =   Log(N!) for N semi-integer
C
      IF((I_SPHER.GE.0).OR.(I_MULT.EQ.1)) THEN
        GLG(1)=0.0
        GLD(1,1)=0.D0
        GLD(1,2)=DLOG(SQPI/2.D0)
        DO I=2,NFAC
          J=I-1
          GLG(I)=GLG(J)+ALOG(FLOAT(J))
          GLD(I,1)=GLD(J,1)+DLOG(DFLOAT(J))
          GLD(I,2)=GLD(J,2)+DLOG(DFLOAT(J) +0.5D0)
        ENDDO
      ELSEIF((IFTHET.EQ.1).AND.(ITEST.EQ.1)) THEN
        GLG(1)=0.0
        DO I=2,NFAC
          J=I-1
          GLG(I)=GLG(J)+ALOG(FLOAT(J))
        ENDDO
      ENDIF 
      EXPF(0,0)=1.
      EXPR(0,0)=1.
      FACT1L=0.D0
      DO L=1,2*NL_M-2
        XDEN=1./SQRT(FLOAT(L+L+1))
        DXDEN=1.D0/DSQRT(DFLOAT(L+L+1))
        FACT1L=FACT1L+DLOG(DFLOAT(L))
        FACT2L=DLOG(DFLOAT(L+1))
        DO M1=0,L
          EXPF(M1,L)=EXP(0.5*(GLG(L+M1+1)-GLG(L-M1+1)))
          DEXPF=DEXP(0.5D0*(GLD(L+M1+1,1)-GLD(L-M1+1,1)))
          EXPR(M1,L)=EXP(0.5*(GLG(L+L+1)-GLG(L+M1+1)-GLG(L-M1+1)))
          EXPF2(L,M1)=EXPF(M1,L)*XDEN
          DEXPF2(L,M1)=DEXPF*DXDEN
          IF(M1.GT.0) THEN 
            FACT2L=FACT2L+DLOG(DFLOAT(1+L+M1))
          ENDIF
          IF(L.LT.NL_M) THEN
            DO M2=0,L
              CF(L,M1,M2)=SQRT(FLOAT((L*L-M1*M1)*(L*L-M2*M2)))/
     1                    FLOAT(L)
            ENDDO
          ENDIF
        ENDDO
        FSQ(L)=EXP(0.5*REAL(FACT2L-FACT1L))
        DFSQ(L)=DEXP(0.5D0*(FACT2L-FACT1L))
      ENDDO
C
C  Setting up of the initial and final values for the angular momenta
C            and the number of polarization directions
C
      IF(INITL.LT.-2) THEN                            
        INITL=1                                                         
        WRITE(IUO1,511)                                                        
      ENDIF 
      IF(INITL.GT.2) THEN
        IF((INITL.NE.20).OR.(INITL.NE.40)) THEN
          IF(INITL.NE.42) THEN
            INITL=1
            WRITE(IUO1,511)                                                        
          ENDIF
        ENDIF
      ENDIF
C     
      NEPS=2-ABS(IPOL)
      IF(SPECTRO.NE.'RES') THEN
        NEPS=2-ABS(IPOL)
      ELSE
        NEPS_I=2-ABS(IPOL_I)
        NEPS_O=2-ABS(IPOL_O)
      ENDIF
C
      IF(IDICHR.GE.1) THEN
        IF(SPECTRO.NE.'RES') THEN
          NEPS=1
        ELSE
          NEPS_I=1
          NEPS_O=1
        ENDIF
      ENDIF
C
      IF(SPECTRO.NE.'RES') THEN
C
C  Non REXS case
C
        IF((INITL.GE.-2).AND.(INITL.LE.2)) THEN      
          LF1=LI+INITL
          LF2=LF1
          ISTEP_LF=1
        ELSEIF(INITL.EQ.20) THEN
          LF1=LI-1
          LF2=LI+1
          ISTEP_LF=2
        ELSEIF(INITL.EQ.40) THEN
          LF1=LI-2
          LF2=LI+2
          ISTEP_LF=2
        ELSEIF(INITL.EQ.24) THEN
          LF1=LI-2
          LF2=LI+2
          ISTEP_LF=1
        ELSEIF(INITL.EQ.42) THEN
          LF1=LI-2
          LF2=LI+2
          ISTEP_LF=1
        ENDIF
      ELSE
C
C  REXS incoming beam
C
        IF((INITL_I.GE.-2).AND.(INITL_I.LE.2)) THEN      
          LF1_I=LI+INITL_I
          LF2_I=LF1_I
          ISTEP_LF_I=1
        ELSEIF(INITL_I.EQ.20) THEN
          LF1_I=LI-1
          LF2_I=LI+1
          ISTEP_LF_I=2
        ELSEIF(INITL_I.EQ.40) THEN
          LF1_I=LI-2
          LF2_I=LI+2
          ISTEP_LF_I=2
        ELSEIF(INITL.EQ.24) THEN
          LF1_I=LI-2
          LF2_I=LI+2
          ISTEP_LF_I=1
        ELSEIF(INITL_I.EQ.42) THEN
          LF1_I=LI-2
          LF2_I=LI+2
          ISTEP_LF_I=1
        ENDIF
C
C  REXS outgoing beam
C
        IF((INITL_O.GE.-2).AND.(INITL_O.LE.2)) THEN      
          LF1_O=LI+INITL_O
          LF2_O=LF1_O
          ISTEP_LF_O=1
        ELSEIF(INITL_O.EQ.20) THEN
          LF1_O=LI-1
          LF2_O=LI+1
          ISTEP_LF_O=2
        ELSEIF(INITL_O.EQ.40) THEN
          LF1_O=LI-2
          LF2_O=LI+2
          ISTEP_LF_O=2
        ELSEIF(INITL_O.EQ.24) THEN
          LF1_O=LI-2
          LF2_O=LI+2
          ISTEP_LF_O=1
        ELSEIF(INITL_O.EQ.42) THEN
          LF1_O=LI-2
          LF2_O=LI+2
          ISTEP_LF_O=1
        ENDIF
C
        LF1=MIN(LF1_I,LF1_O)
        LF2=MAX(LF1_I,LF1_O)
        ISTEP_LF=MIN(ISTEP_LF_I,ISTEP_LF_O)
C
      ENDIF
C
C  Initialization of the values of ji if spin-orbit is taken
C                  into account.
C
C  Here :   JI is the loop index going from JF1 to JF2 with :
C
C                   JI=1    : ji = li + 1/2
C                   JI=2    : ji = li - 1/2
C
      IF(I_SO.LE.0) THEN
        JF1=1
        JF2=2
      ELSEIF(I_SO.EQ.1) THEN
        IF(S_O.EQ.'1/2') THEN
          IF(LI.EQ.0) THEN
            JF1=1
            JF2=1
          ELSEIF(LI.EQ.1) THEN
            JF1=2
            JF2=2
          ENDIF
        ELSEIF(S_O.EQ.'3/2') THEN
          IF(LI.EQ.1) THEN
            JF1=1
            JF2=1
          ELSEIF(LI.EQ.2) THEN
            JF1=2
            JF2=2
          ENDIF
        ELSEIF(S_O.EQ.'5/2') THEN
          IF(LI.EQ.2) THEN
            JF1=1
            JF2=1
          ELSEIF(LI.EQ.3) THEN
            JF1=2
            JF2=2
          ENDIF
        ELSEIF(S_O.EQ.'7/2') THEN
          IF(LI.EQ.3) THEN
            JF1=1
            JF2=1
          ELSEIF(LI.EQ.4) THEN
            JF1=2
            JF2=2
          ENDIF
        ELSEIF(S_O.EQ.'9/2') THEN
          IF(LI.EQ.4) THEN
            JF1=1
            JF2=1
          ELSE
            RETURN 7
          ENDIF
        ELSE
          RETURN 7
        ENDIF
      ELSEIF(I_SO.EQ.2) THEN
        JF1=1
        JF2=2
      ELSE
        RETURN 7
      ENDIF
C
      IF(NI.LE.5) THEN
         NNL=NI*(NI-1)/2 +LI+1
      ELSEIF(NI.EQ.6) THEN
         NNL=NI*(NI-1)/2 +LI
      ELSEIF(NI.EQ.7) THEN
         NNL=NI*(NI-1)/2 +LI-3
      ENDIF
      NNL1=NNL
      NNL2=NNL
      LI_1=LI
      LI_2=LI
C
C  Storage of the Clebsch-Gordan coefficients for the spin-orbit 
C  dependent coupling matrix elements in the array CG(MJI,JI,JSPIN).
C
C  Here :           JI=1    : ji = li + 1/2
C                   JI=2    : ji = li - 1/2
C                   MJI     : mji + 1/2 
C                   JSPIN=1 : msi = +1/2
C                   JSPIN=2 : msi = -1/2
C
C              so that all indices remain integer
C
      IF((I_SO.GT.0).OR.(ISPIN.EQ.1)) THEN
        DO JS=1,2
          DO JI=1,2
            DO MJI=-LI,LI+1
              CG(MJI,JI,JS)=0.0
            ENDDO
          ENDDO
        ENDDO
        DO MJI=-LI,LI+1
          CG(MJI,1,1)=SQRT(FLOAT(LI+MJI)/FLOAT(LI+LI+1))
          CG(MJI,1,2)=SQRT(FLOAT(LI-MJI+1)/FLOAT(LI+LI+1))
          IF((MJI.GT.-LI).AND.(MJI.LT.LI+1)) THEN
            CG(MJI,2,1)=-SQRT(FLOAT(LI-MJI+1)/FLOAT(LI+LI+1))
            CG(MJI,2,2)=SQRT(FLOAT(LI+MJI)/FLOAT(LI+LI+1)) 
          ENDIF 
        ENDDO
      ENDIF
C
C
C  Storage of the Clebsch-Gordan coefficients for the Auger multiplet 
C                    dependent coupling matrix elements 
C                   in the array CGA(LJ1,MJ1,LJ2,MJ2,LJ).
C
C  Here :       LJ1  is an integer index related to J1 (LJ1=2*J1)
C               LMJ1  is an integer index related to MJ1 (LMJ1=2*MJ1)
C               LJ2  is an integer index related to J2 (LJ2=2*J2)
C               LMJ2  is an integer index related to MJ2 (LMJ2=2*MJ2)
C               LJ is an integer index related to J :
C                     J = FLOAT(LJ) for J integer
C                     J = FLOAT(LJ) + 0.5 for J half integer
C
C              so that all indices remain integer
C
      IF((SPECTRO.EQ.'AED').OR.(SPECTRO.EQ.'APC')) THEN
        IF(I_MULT.EQ.1) THEN
          N=3
          MJ3=0.D0
          LJ_MAX=2*(LI_I+LI_A+1)
          DO LJ1=0,LJ_MAX
            J1=DFLOAT(LJ1)/2.D0
            DO LMJ1=-LJ1,LJ1,2
              MJ1=DFLOAT(LMJ1)/2.D0
              DO LJ2=0,LJ_MAX
                J2=DFLOAT(LJ2)/2.D0
                DO LMJ2=-LJ2,LJ2,2
                  MJ2=DFLOAT(LMJ2)/2.D0
                  CALL N_J(J1,MJ1,J2,MJ2,MJ3,NJ,I_INT,N)
C
                  JJ12=J1-J2
                  JL12=MJ1-MJ2
C
                  LJ12=INT(JJ12+SIGN(SMALL,JJ12))
                  LL12=INT(JL12+SIGN(SMALL,JL12))
C
                  JJ_MIN=ABS(LJ12)
                  JJ_MAX=J1+J2
                  LJJ_MIN=INT(JJ_MIN+SIGN(SMALL,JJ_MIN))
                  LJJ_MAX=INT(JJ_MAX+SIGN(SMALL,JJ_MAX))
C
                  DO LJJ=LJJ_MIN,LJJ_MAX,1
                    IF(I_INT.EQ.1) THEN
                      JJ=DFLOAT(LJJ)
                    ELSE
                      JJ=DFLOAT(LJJ)+0.5D0
                    ENDIF 
                    L_EXP=INT(J1-J2+MJ1+MJ2)
                    IF(MOD(L_EXP,2).EQ.0) THEN  
                      CGA(LJ1,LMJ1,LJ2,LMJ2,LJJ)=NJ(LJJ)*
     1                                           SQRT(2.*REAL(JJ)+1.)  
                    ELSE
                      CGA(LJ1,LMJ1,LJ2,LMJ2,LJJ)=-NJ(LJJ)*
     1                                            SQRT(2.*REAL(JJ)+1.)   
                    ENDIF      
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C  Storage of another of the spin Clebsch-Gordan used 
C    when the Auger line is multiplet-resolved. It 
C    originates from the coupling of SA and SC, 
C    the spins of the Auger electron of the original 
C    core electron (which is supposed to be the same 
C    as that of the photoelectron).
C
C    CG_S(I,J,K) with : I = 1 ---> MSA = -1/2
C                       I = 2 ---> MSA =  1/2
C                       J = 1 ---> MSC = -1/2
C                       J = 2 ---> MSC =  1/2
C                       K = 1 ---> S   =  0
C                       K = 2 ---> S   =  1
C
C                       MS = MSA+MSC
C
      IF(I_MULT.EQ.1) THEN
        CG_S(1,1,1)=0.
        CG_S(1,1,2)=1.
        CG_S(1,2,1)=-0.707107
        CG_S(1,2,2)= 0.707107
        CG_S(2,1,1)= 0.707107
        CG_S(2,1,2)= 0.707107
        CG_S(2,2,1)= 0.
        CG_S(2,2,2)= 1.
      ENDIF
C
C  Initialization of the variables used when only one multiplet
C    is taken into account in the Auger peak
C
      IF(I_MULT.EQ.1) THEN
        MULTIPLET=CHAR(48+IM1)//MULT//CHAR(48+IM2)
        IF(MOD(IM1,2).EQ.0) THEN
          WRITE(IUO1,522) IM1
          STOP
        ENDIF
        S_MUL=(IM1-1)/2
        J_MUL=IM2
        IF(MULT.EQ.'S') THEN
          L_MUL=0
        ELSEIF(MULT.EQ.'P') THEN
          L_MUL=1
        ELSEIF(MULT.EQ.'D') THEN
          L_MUL=2
        ELSEIF(MULT.EQ.'F') THEN
          L_MUL=3
        ELSEIF(MULT.EQ.'G') THEN
          L_MUL=4
        ELSEIF(MULT.EQ.'H') THEN
          L_MUL=5
        ELSEIF(MULT.EQ.'I') THEN
          L_MUL=6
        ELSEIF(MULT.EQ.'K') THEN
          L_MUL=7
        ELSEIF(MULT.EQ.'L') THEN
          L_MUL=8
        ELSEIF(MULT.EQ.'M') THEN
          L_MUL=9
        ELSE
          WRITE(IUO1,521) MULTIPLET
          STOP
        ENDIF
      ENDIF
C
C..........  Setting up various parameters  ..........
C..........        for consistency          ..........
C
      IF(SPECTRO.EQ.'STM') THEN
        ISPEED_IN=ISPEED_TI
        IATTS_IN=IATTS_TI
        ISPEED_EX=ISPEED_TI
        IATTS_EX=IATTS_TI
        ISPEED_O1=ISPEED
        IATTS_O1=IATTS
        ISPEED_O2=ISPEED
        IATTS_O2=IATTS
      ELSE
        ISPEED_IN=ISPEED
        IATTS_IN=IATTS
        ISPEED_EX=ISPEED
        IATTS_EX=IATTS
        ISPEED_O1=ISPEED
        IATTS_O1=IATTS
        ISPEED_O2=ISPEED
        IATTS_O2=IATTS
      ENDIF
C
      ISPIN_IN=ISPIN
      IDICHR_IN=IDICHR
      NSPIN_IN=NSPIN
      NSPIN2_IN=NSPIN2
      ISPIN_EX=ISPIN
      IDICHR_EX=IDICHR
      NSPIN_EX=NSPIN
      NSPIN2_EX=NSPIN2
      ISPIN_O1=ISPIN
      IDICHR_O1=IDICHR
      NSPIN_O1=NSPIN
      NSPIN2_O1=NSPIN2
      ISPIN_O2=ISPIN
      IDICHR_O2=IDICHR
      NSPIN_O2=NSPIN
      NSPIN2_O2=NSPIN2
C
      IF(NEL.EQ.1) THEN
        IF(INC.EQ.1) THEN
          ITRTL=ITRTL_IN
        ELSEIF(EXC.EQ.1) THEN
          ITRTL=ITRTL_EX
        ELSEIF(OUT1.EQ.1) THEN
          ITRTL=ITRTL_O1
        ENDIF
      ELSEIF(NEL.EQ.2) THEN
        ITRTL=MAX(ITRTL_O1,ITRTL_O2)
      ENDIF
C
C..........  Check of the dimensioning in the Gaussian case  ..........
C
      CALL STOP_EXT(I_EXT,I_EXT_A,SPECTRO)
C
C..........  Redistributing the common blocks  ..........
C
      CALL CHG_COM
C
C....................   Read FORMAT   ....................
C
C
   1  FORMAT(A7)
   2  FORMAT(21X,12A4)
   3  FORMAT(7X,A3,9X,A1,9X,I1,6X,I4)
   4  FORMAT(8X,F6.3,4X,F6.3,4X,F6.3,3X,A3)
   5  FORMAT(49X,A7)
   6  FORMAT(7X,F6.2,4X,F6.2,4X,F6.2)
   7  FORMAT(8X,I2,8X,I2,8X,I2,8X,I2)
   8  FORMAT(8X,I2,8X,F6.3,3X,I3,9X,I1)
   9  FORMAT(8X,F9.6,1X,F9.6,1X,F9.6,2X,A2,2X,I2)
  10  FORMAT(9X,I1,8X,I2,7X,F5.1,5X,F5.1)
  11  FORMAT(7X,F5.1,3(5X,F5.1))
  12  FORMAT(7X,F6.2,4X,F6.2,6X,I1)
  13  FORMAT(7X,A3,9X,I1,9X,I1,9X,I1)
  14  FORMAT(8X,I2,6X,I4,6X,F7.2,3X,F7.2)
!  15  FORMAT(7X,I3,7X,I3,7X,I3,7X,I3)
!  16  FORMAT(6X,F7.2,3X,F7.2,3X,F7.2,5X,F6.3)
  17  FORMAT(6X,F7.2,3X,F7.2,3X,F7.2,5X,I2)
  18  FORMAT(9X,I1,9X,I1,8X,F5.2,6X,I1)
  19  FORMAT(7X,I3,6X,F7.2,3X,F7.2,3X,F7.2)
  20  FORMAT(8X,I1,A1,8X,A3,7X,I2,8X,I2)
  21  FORMAT(8X,I2,8X,I2,9X,I1,9X,I1)
  22  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)
  23  FORMAT(8X,I2)
  24  FORMAT(8X,I2,3(8X,I2))
  25  FORMAT(9X,I1,8X,I2,6X,I4,8X,F6.2)
  26  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)
  27  FORMAT(9X,I1,6X,F6.2,7X,I1,7X,F6.2)
  28  FORMAT(9X,I1,9X,I1,7X,F8.4,4X,I1)
  29  FORMAT(9X,I1,7X,F6.2,4X,A3)
  30  FORMAT(9X,I1,8X,I2,9X,I1,9X,I1)
  31  FORMAT(9X,I1,6X,F7.2,3X,F7.2,6X,F4.2)
  32  FORMAT(9X,I1,7X,F6.2)
  33  FORMAT(8X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5)
  34  FORMAT(9X,A24,5X,I2)
  35  FORMAT(18X,I2,8X,I2,8X,I2)
  36  FORMAT(18X,A2,8X,A2,8X,A2)
  37  FORMAT(18X,F8.5,2X,F8.5,2X,F8.5)
  38  FORMAT(9X,I1,7X,F5.1,5X,F5.1,5X,F5.1)
  39  FORMAT(8X,A1,I1,8X,I2)
  40  FORMAT(8X,A1,I1,8X,A1,I1,8X,A1,I1)
  42  FORMAT(9X,I1,8X,I1,A1,I1)
  43  FORMAT(6X,F7.2,3X,F7.2,5X,I2,8X,I2)
  44  FORMAT(8X,A1,I1,8X,I2,9X,I1,9X,I1,9X,I1)
  46  FORMAT(7X,I3,6X,F7.2,3X,F7.2)
  47  FORMAT(7X,I3,6X,F7.2)
  48  FORMAT(8X,F5.2,5X,F5.2,5X,F5.2,5X,F5.2)
  49  FORMAT(27X,A4)
  50  FORMAT(48X,A4)
  51  FORMAT(37X,A4)
  52  FORMAT(6X,I4,8X,F6.4,3X,A3)
  53  FORMAT(7X,I3,6X,F7.2,3X,F7.2,6X,I1)
  55  FORMAT(8X,I2,6X,A4,9X,F7.5,2X,F6.3)
  56  FORMAT(5X,I5,6X,I4,6X,I4,8X,F6.3)
  57  FORMAT(9X,I1,9X,I1,9X,I1,9X,I1)
  58  FORMAT(8X,I2,6X,F7.2,3X,F7.2)
C
C
C....................   Write FORMAT   ....................
C
C
 100  FORMAT(//////////,'******************************',
     1       '****************************************************')
 101  FORMAT('*********************',40X,'*********************')
 102  FORMAT('*********************',10A4,'*********************')
 103  FORMAT(10X,A3,9X,A1,9X,I1,6X,I4,9X,'CRIST,CENTR,IBAS,NAT')
 104  FORMAT(11X,F6.3,4X,F6.3,4X,F6.3,15X,'A,BSURA,CSURA')
 105  FORMAT(///,'ooooooooooooooooooooooooooooooooooooooooo',
     1       'ooooooooooooooooooooooooooooooooooooooooo',/,
     2       'oooooooooooooooo',50X,'oooooooooooooooo',/,
     3       'oooooooooooooooo   INPUT DATA FILE   :  ',A24,
     4       '  oooooooooooooooo',/,'oooooooooooooooo',50X,
     5       'oooooooooooooooo',/,'oooooooooooooooooooooooooooo',
     6       'ooooooooooooooooooooooooooooooooooooooooooooooooo',
     7       'ooooo',///)
 106  FORMAT(10X,F6.2,4X,F6.2,4X,F6.2,16X,'ALPHAD,BETAD,GAMMAD')
 107  FORMAT(11X,I2,8X,I2,8X,I2,8X,I2,9X,'H,K,I,L')
 108  FORMAT(12X,I1,8X,F6.3,3X,I3,9X,I1,9X,'NIV,COUPUR,ITEST,IESURF')
 109  FORMAT(11X,F9.6,1X,F9.6,1X,F9.6,2X,A2,2X,I2,4X,'ATBAS,CHEM(NAT)',
     1       ',NZAT(NAT)')
 110  FORMAT(12X,I1,8X,I2,7X,F5.1,5X,F5.1,7X,'IREL,NREL,PCREL(NREL)')
 111  FORMAT(12X,I1,39X,'IPOL')
 112  FORMAT(10X,F6.2,4X,F6.2,6X,I1,19X,'OMEGA1,OMEGA2,IADS')
 113  FORMAT(12X,I1,9X,I1,29X,'ISPIN,IDICHR')
 114  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IPHI,NPHI,PHI0,PHI1')
 115  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'ITHETA,NTHETA,THETA0,',
     1       'THETA1')
 116  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IE,NE,E0,E1')
 117  FORMAT(9X,F7.2,3X,F7.2,3X,F7.2,5X,I2,9X,'THLUM,PHILUM,ELUM,IPOL')
 118  FORMAT(12X,I1,9X,I1,8X,F5.2,6X,I1,9X,'IMOD,IMOY,ACCEPT,ICHKDIR')
 119  FORMAT(10X,I3,6X,F7.2,3X,F7.2,3X,F7.2,6X,'NE,EK_INI,',
     1       'EK_FIN,EPH_INI')
 120  FORMAT(11X,I1,A1,8X,A3,7X,I2,8X,I2,9X,'LI,S-O,INITL,I_SO')
 121  FORMAT(11X,I2,8X,I2,9X,I1,9X,I1,9X,'NO,NDIF,ISPHER,I_GR')
 122  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'ISFLIP,IR_DIA,ITRTL,I_TEST')
 123  FORMAT(11X,I2,3(8X,I2),9X,'NEMET,IEMET(NEMET)')
 124  FORMAT(12X,I1,8X,I2,6X,I4,7X,F6.2,6X,'ISOM,NONVOL,NPATH,VINT')
 125  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'IFWD,NTHOUT,I_NO,I_RA')
 126  FORMAT(12X,I1,7X,F6.2,6X,I1,7X,F6.2,6X,'N_RA(NAT),THFWD(NAT)',
     1       ',IBWD(NAT),THBWD(NAT)')
 127  FORMAT(12X,I1,9X,I1,7X,F8.4,4X,I1,9X,'IPW,NCUT,PCTINT,IPP')
 128  FORMAT(12X,I1,7X,F6.2,4X,A3,19X,'ILENGTH,RLENGTH,UNLENGTH')
 129  FORMAT(12X,I1,8X,I2,9X,I1,9X,I1,9X,'IDWSPH,ISPEED,IATT,IPRINT')
 130  FORMAT(12X,I1,6X,F7.2,3X,F7.2,6X,F4.2,6X,'IDCM,TD,T,RSJ')
 131  FORMAT(12X,I1,7X,F6.2,26X,'ILPM,XLPM0')
 132  FORMAT(11X,F8.5,33X,'UJ2(NAT)  : ',
     1       'SUBSTRATE')
 133  FORMAT(11X,F8.5,33X,'UJ2(NATA) : ',
     1       'ADSORBATES')
 134  FORMAT(11X,A1,I1,8X,I2,6X,F7.2,3X,F7.2,6X,'EDGE,INITL,THLUM,',
     1       'PHILUM')
 135  FORMAT(11X,A1,I1,8X,A1,I1,8X,A1,I1,19X,'EDGE_C,EDGE_I,',
     1       'EDGE_A')
 136  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IPHI_A,NPHI_A,PHI0_A,',
     1       'PHI1_A')
 137  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'ITHETA_A,NTHETA_A,',
     1       'THETA0_A,THETA1_A')
 138  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IFTHET_A,NFTHET_A,R0_A,',
     1       'R1_A')
 139  FORMAT(11X,I2,39X,'I_INT')
 140  FORMAT(12X,I1,8X,I1,A1,I1,28X,'I_MULT,MULT')
 141  FORMAT(12X,I1,39X,'ISPIN')
 142  FORMAT(9X,F7.2,3X,F7.2,26X,'TH_INI,PHI_INI')
 143  FORMAT(9X,F7.2,3X,F7.2,6X,I1,8X,I2,9X,'THLUM_I,PHILUM_I',
     1       ',IPOL_I,ISR_I')
 144  FORMAT(11X,A1,I1,9X,I1,9X,I1,9X,I1,9X,'EDGE,IE,IPHI,ITHETA')
 145  FORMAT(9X,F7.2,3X,F7.2,6X,I1,8X,I2,9X,'THLUM_O,PHILUM_O',
     1       ',IPOL_O,ISR_O')
 146  FORMAT(10X,I3,6X,F7.2,3X,F7.2,16X,'NE,E_INI,E_FIN')
 147  FORMAT(10X,I3,6X,F7.2,26X,'NPHI,PHI_FIN')
 148  FORMAT(11X,F5.2,5X,F5.2,5X,F5.2,5X,F5.2,6X,'G_MAX,G_HOL,G_SLO,',
     1       'G_CEN')
 149  FORMAT(10X,F6.2,36X,'VINT')
 150  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'ITEMP,NTEMP,TEMP1,TEMP2')
 153  FORMAT(10X,I3,6X,F7.2,3X,F7.2,6X,I1,9X,'NE_EIG,EIG_INI,',
     1       'EIG_FIN,I_DAMP')
 154  FORMAT(11X,I2,8X,I2,8X,I2,8X,I2,9X,'I_SPECTRUM(NE)') 
 155  FORMAT(11X,I2,6X,A4,9X,F7.5,2X,F6.3,5X,'I_PWM,METHOD,ACC,EXPO')
 156  FORMAT(8X,I5,6X,I4,6X,I4,8X,F6.3,5X,'N_MAX,N_ITER,N_TABLE,SHIFT')
 157  FORMAT(12X,I1,9X,I1,9X,I1,9X,I1,9X,'I_XN,I_VA,I_GN,I_WN')
 158  FORMAT(11X,I2,6X,F7.2,3X,F7.2,16X,'L,ALPHA,BETA')
 159  FORMAT(10X,I3,6X,F7.2,26X,'NTHETA,THETA_FIN')
C
 201  FORMAT(///,21X,12A4,////)
 203  FORMAT('**************************************************',
     1       '********************************',//////////)
 204  FORMAT(11X,F6.3,5X,I1,29X,'A,IBAS')
 206  FORMAT(10X,F6.2,36X,'ALPHAD')
 207  FORMAT(11X,I2,8X,I2,8X,I2,19X,'H,K,L')
 209  FORMAT(11X,F9.6,1X,F9.6,1X,F9.6,12X,'VECBAS')
 210  FORMAT(10X,F5.1,3(5X,F5.1),7X,'PCREL(NREL)')
 211  FORMAT(20X,'SUBSTRATE : ',10(F5.1,','))
 212  FORMAT(32X,I1,19X,'IADS')
 216  FORMAT(11X,I2,6X,I4,6X,F7.2,3X,F7.2,6X,'IFTHET,NFTHET,R0,R1')
 223  FORMAT(11X,I2,1(8X,I2),29X,'NEMET,IEMET(NEMET)')
 232  FORMAT(11X,F8.5,2X,F8.5,23X,'UJ2(NAT)  : ',
     1       'SUBSTRATE')
 233  FORMAT(11X,F8.5,2X,F8.5,23X,'UJ2(NATA) : ',
     1       'ADSORBATES')
C
 303  FORMAT(/,33X,'ATOMS OF TYPE ',I1,' :',/)
 304  FORMAT(11X,F6.3,35X,'A')
 306  FORMAT(10X,F6.2,36X,'BETAD')
 309  FORMAT(11X,F9.6,1X,F9.6,1X,F9.6,12X,'XADS,YADS,ZADS')
 311  FORMAT(20X,'ADSORBATE : ',3(F5.1,','))
 323  FORMAT(11X,I2,2(8X,I2),19X,'NEMET,IEMET(NEMET)')
 332  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,13X,'UJ2(NAT)  : ',
     1       'SUBSTRATE')
 333  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,13X,'UJ2(NATA) : ',
     1       'ADSORBATES')
C
 404  FORMAT(11X,F6.3,4X,F6.3,25X,'A,CSURA')
 432  FORMAT(11X,F8.5,2X,F8.5,2X,F8.5,2X,F8.5,3X,'UJ2(NAT)  : ',
     1       'SUBSTRATE')
C
 501  FORMAT(//,30X,'POSITION OF THE ADSORBATES :')
 502  FORMAT(///,25X,'VALUE OF THE RELAXATIONS :',/)
 503  FORMAT(///,14X,'TYPE OF CALCULATION : AZIMUTHAL PHOTOELECTRON',
     1       ' DIFFRACTION')
 504  FORMAT(///,18X,'TYPE OF CALCULATION : FINE STRUCTURE ',
     1       'OSCILLATIONS')
 505  FORMAT(///,16X,'TYPE OF CALCULATION : POLAR PHOTOELECTRON',
     1       ' DIFFRACTION')
 506  FORMAT(///,23X,'TYPE OF CALCULATION : SCATTERING FACTOR')
 507  FORMAT(///,28X,'TYPE OF CALCULATION : EXAFS')
 508  FORMAT(///,2X,' <<<<<<<<<<  THE THETA VARIATION EXCEEDS THE ',
     1       'PHYSICAL LIMITS (-90,+90)  >>>>>>>>>>',///)
 509  FORMAT(///,2X,' <<<<<<<<<<  THE THLUM VARIATION EXCEEDS THE ',
     1       'PHYSICAL LIMITS (-90,+90)  >>>>>>>>>>',///)
 510  FORMAT(///,4X,' <<<<<<<<<<  AS THE CALCULATION HAS MORE THAN ',
     1       '250 POINTS, SOME OUTPUTS HAVE BEEN SUPRESSED  >>>>>>>>>>',
     2       ///)
 511  FORMAT(///,4X,' <<<<<<<<<<  INCORRECT VALUE OF INITL, THE ',
     1       'CALCULATION IS PERFORMED WITH INITL = 1  >>>>>>>>>>')
 512  FORMAT(///,4X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE A SPIN RESOLVED ',
     1       'EXAFS EXPERIMENT : DECREASE IDICHR  >>>>>>>>>>')
 513  FORMAT(///,15X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE IPOL = 0 AND ',
     1       'IDICHR > 0  >>>>>>>>>>')
 514  FORMAT(///,15X,' <<<<<<<<<<  IMPOSSIBLE TO HAVE IDICHR = 2 AND ',
     1       'ISPIN = 0  >>>>>>>>>>')
 515  FORMAT(///,12X,'TYPE OF CALCULATION : AZIMUTHAL AUGER ELECTRON',
     1       ' DIFFRACTION')
 516  FORMAT(///,16X,'TYPE OF CALCULATION : POLAR AUGER ELECTRON',
     1       ' DIFFRACTION')
 517  FORMAT(///,10X,'TYPE OF CALCULATION : AUGER PHOTOELECTRON ',
     1       'COINCIDENCE SPECTROSCOPY')
 518  FORMAT(///,9X,'------------------------  FIRST ELECTRON : ',
     1       '------------------------')
 519  FORMAT(///,9X,'------------------------ SECOND ELECTRON : ',
     1       '------------------------')
 520  FORMAT(///,9X,'----------------------------------------------',
     1       '----------------------')
 521  FORMAT(///,4X,' <<<<<<<<<<  ',A3,' IS NOT IMPLEMENTED IN THIS ',
     1       'VERSION  >>>>>>>>>>')
 522  FORMAT(///,4X,' <<<<<<<<<<  WRONG NAME FOR THE MULTIPLET',
     1       '  >>>>>>>>>>',/,4X,' <<<<<<<<<<  ODD NUMBER ',
     2       'EXPECTED INSTEAD OF',I2,'  >>>>>>>>>>')
 523  FORMAT(///,4X,' <<<<<<<<<<  BOTH DETECTOR DIRECTIONS MUST BE ',
     1       'EITHER INTERNAL OR EXTERNAL  >>>>>>>>>>',/,8X,
     2              ' -----> PROCEEDING WITH EXTERNAL DIRECTIONS',/)
 524  FORMAT(///,4X,' <<<<<<<<<<  AVERAGING OVER ',I3,' DOMAINS ',
     1       'FOR PHOTOELECTRON   >>>>>>>>>>',/,4X,
     2       ' <<<<<<<<<<  AVERAGING OVER ',I3,' DOMAINS ',
     3       'FOR AUGER ELECTRON   >>>>>>>>>>',/,8X,
     4       ' -----> IMPOSSIBLE : CHECK INPUT FILES !')
 525  FORMAT(///,14X,'ATOMIC CALCULATION : Z AXIS ALONG POLARIZATION ',
     1          'DIRECTION',/,'  ',/,'   ',/,'  ')
 526  FORMAT(///,18X,'ATOMIC CALCULATION : Z AXIS ALONG LIGHT ',
     1          'DIRECTION',/,'  ',/,'   ',/,'  ')
 527  FORMAT(///,11X,'TYPE OF CALCULATION : FULL HEMISPHERE',
     1       ' PHOTOELECTRON DIFFRACTION')
 528  FORMAT(///,10X,'TYPE OF CALCULATION : FULL HEMISPHERE',
     1       ' AUGER ELECTRON DIFFRACTION')
 529  FORMAT(///,14X,'TYPE OF CALCULATION : AZIMUTHAL LEED',
     1       ' VARIATIONS')
 530  FORMAT(///,11X,'TYPE OF CALCULATION : FULL HEMISPHERE',
     1       ' LEED')
 531  FORMAT(///,18X,'TYPE OF CALCULATION : LEED ENERGY ',
     1       'VARIATIONS')
 532  FORMAT(///,16X,'TYPE OF CALCULATION : POLAR LEED',
     1       ' VARIATIONS')
 533  FORMAT(///,16X,'TYPE OF CALCULATION : RESONANT ELASTIC X-RAY',
     1       ' SCATTERING')
 534  FORMAT(//,20X,'THIS CALCULATION OUTPUTS ALSO THE AMPLITUDES')
 535  FORMAT(///,16X,'TYPE OF CALCULATION : ELECTRON ENERGY LOSS',
     1       ' SPECTROSCOPY')
 536  FORMAT(///,22X,'TYPE OF CALCULATION : EIGENVALUE',
     1       ' ANALYSIS')
C
 701  FORMAT(6X,I1,1X,I3,2X,I4)
 702  FORMAT(6X,I1,1X,I3,3X,I3)
 703  FORMAT(15X,F8.3,3X,F8.3)
 713  FORMAT(6X,I1,1X,I3)
C
      RETURN
C
      END 
C
C=======================================================================
C
      SUBROUTINE CHG_COM
C
C   This subroutine redistributes some common blocks depending 
C     on the number of electrons involved in the spectroscopy.
C     It allows in particular to simplify the structure when 
C     only one electron is involved in the process.
C
C                                         Last modified : 13 Oct 2010
C
      INCLUDE 'spec.inc'
C
      COMMON /APPROX/ NDIF,NO,ISPHER,IFWD,NTHOUT,RTHFWD(NATP_M),
     1                IBWD(NATP_M),RTHBWD(NATP_M),IPW,NCUT,PCTINT,IPP,
     2                ISPEED,IATTS,ILENGTH,RLENGTH
      COMMON /APPROX_IN/ NDIF_IN,NO_IN,ISPHER_IN,IFWD_IN,NTHOUT_IN,
     1                RTHFWD_IN(NATP_M),IBWD_IN(NATP_M),
     2                RTHBWD_IN(NATP_M),IPW_IN,NCUT_IN,PCTINT_IN,IPP_IN,
     3                ISPEED_IN,IATTS_IN,ILENGTH_IN,RLENGTH_IN
      COMMON /APPROX_EX/ NDIF_EX,NO_EX,ISPHER_EX,IFWD_EX,NTHOUT_EX,
     1                RTHFWD_EX(NATP_M),IBWD_EX(NATP_M),
     2                RTHBWD_EX(NATP_M),IPW_EX,NCUT_EX,PCTINT_EX,IPP_EX,
     3                ISPEED_EX,IATTS_EX,ILENGTH_EX,RLENGTH_EX
      COMMON /APPROX_O1/ NDIF_O1,NO_O1,ISPHER_O1,IFWD_O1,NTHOUT_O1,
     1                RTHFWD_O1(NATP_M),IBWD_O1(NATP_M),
     2                RTHBWD_O1(NATP_M),IPW_O1,NCUT_O1,PCTINT_O1,IPP_O1,
     3                ISPEED_O1,IATTS_O1,ILENGTH_O1,RLENGTH_O1
      COMMON /APPROX_O2/ NDIF_O2,NO_O2,ISPHER_O2,IFWD_O2,NTHOUT_O2,
     1                RTHFWD_O2(NATP_M),IBWD_O2(NATP_M),
     2                RTHBWD_O2(NATP_M),IPW_O2,NCUT_O2,PCTINT_O2,IPP_O2,
     3                ISPEED_O2,IATTS_O2,ILENGTH_O2,RLENGTH_O2
      COMMON /CLUELEC/ NEL,NCL,INC,EXC,OUT1,OUT2,TIP,CRYS
      COMMON /INFILES/ INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6,
     1                 INFILE7,INFILE8,INFILE9
      COMMON /INFILES2/ INFILE1B,INFILE2B,INFILE3B,INFILE4B,INFILE5B,
     1                  INFILE6B,INFILE7B,INFILE8B,INFILE9B,INFILE10B,
     2                  INFILE11B,INFILE12B,INFILE13B,INFILE14B,
     3                  INFILE15B,INFILE16B,INFILE17B
      COMMON /INUNITS/ IUI1,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,IUI9
      COMMON /INUNITS2/ IUI1B,IUI2B,IUI3B,IUI4B,IUI5B,IUI6B,IUI7B,IUI8B,
     1                  IUI9B,IUI10B,IUI11B,IUI12B,IUI13B,IUI14B,IUI15B,
     2                  IUI16B,IUI17B
      COMMON /LINLBD/ LBD(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX,NUMAX(NATP_M)
      COMMON /LINLBD_IN/ LBD_IN(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX_IN,
     1                   NUMAX_IN(NATP_M)
      COMMON /LINLBD_EX/ LBD_EX(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX_EX,
     1                   NUMAX_EX(NATP_M)
      COMMON /LINLBD_O1/ LBD_O1(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX_O1,
     1                   NUMAX_O1(NATP_M)
      COMMON /LINLBD_O2/ LBD_O2(-N_MU_M:N_MU_M,0:N_NU_M),LBDMAX_O2,
     2                   NUMAX_O2(NATP_M)
      COMMON /LPMOY/ ILPM,NZA,XM,RH,XLMP0
      COMMON /LPMOY_TMP/ NZA_TMP,XM0_TMP,RH0_TMP
      COMMON /LPMOY_IN/ ILPM_IN,XLPM0_IN
      COMMON /LPMOY_EX/ ILPM_EX,XLPM0_EX
      COMMON /LPMOY_O1/ ILPM_O1,XLPM0_O1
      COMMON /LPMOY_O2/ ILPM_O2,XLPM0_O2
      COMMON /RA/ I_NO,I_RA,N_RA(NATP_M)
      COMMON /RA_IN/ I_NO_IN,I_RA_IN,N_RA_IN(NATP_M)
      COMMON /RA_EX/ I_NO_EX,I_RA_EX,N_RA_EX(NATP_M)
      COMMON /RA_O1/ I_NO_O1,I_RA_O1,N_RA_O1(NATP_M)
      COMMON /RA_O2/ I_NO_O2,I_RA_O2,N_RA_O2(NATP_M)
      COMMON /RESEAU/ NCRIST,NCENTR,IBAS,NAT,A,BSURA,CSURA,UNIT
      COMMON /SPIN/ ISPIN,IDICHR,NSPIN,NSPIN2,ISFLIP,IR_DIA,NSTEP
      COMMON /SPIN_IN/ ISPIN_IN,IDICHR_IN,NSPIN_IN,NSPIN2_IN,ISFLIP_IN,
     1                 IR_DIA_IN,NSTEP_IN
      COMMON /SPIN_EX/ ISPIN_EX,IDICHR_EX,NSPIN_EX,NSPIN2_EX,ISFLIP_EX,
     1                 IR_DIA_EX,NSTEP_EX
      COMMON /SPIN_O1/ ISPIN_O1,IDICHR_O1,NSPIN_O1,NSPIN2_O1,ISFLIP_O1,
     1                 IR_DIA_O1,NSTEP_O1
      COMMON /SPIN_O2/ ISPIN_O2,IDICHR_O2,NSPIN_O2,NSPIN2_O2,ISFLIP_O2,
     1                 IR_DIA_O2,NSTEP_O2
      COMMON /TYPEXP/ SPECTRO,INTERACT,STEREO
C
      CHARACTER*24 INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6
      CHARACTER*24 INFILE7,INFILE8,INFILE9
      CHARACTER*24 INFILE1B,INFILE2B,INFILE3B,INFILE4B,INFILE5B
      CHARACTER*24 INFILE6B,INFILE7B,INFILE8B,INFILE9B,INFILE10B
      CHARACTER*24 INFILE11B,INFILE12B,INFILE13B,INFILE14B,INFILE15B
      CHARACTER*24 INFILE16B,INFILE17B
      CHARACTER*7 INTERACT
      CHARACTER*3 SPECTRO,STEREO,UNIT
C
      INTEGER EXC,OUT1,OUT2,TIP
C
C  One electron spectroscopies (PhD,LEED,XAS,AED,REXS,...)
C
      IF(NEL.EQ.1) THEN
C
        IF(INC.EQ.1) THEN
          NDIF=NDIF_IN
          NO=NO_IN
          ISPHER=ISPHER_IN
          IFWD=IFWD_IN
          NTHOUT=NTHOUT_IN
          IPW=IPW_IN
          NCUT=NCUT_IN
          PCTINT=PCTINT_IN
          IPP=IPP_IN
          ISPEED=ISPEED_IN
          IATTS=IATTS_IN
          ILENGTH=ILENGTH_IN
          RLENGTH=RLENGTH_IN
          I_NO=I_NO_IN
          I_RA=I_RA_IN
          ISPIN=ISPIN_IN
          IDICHR=IDICHR_IN
          NSPIN=NSPIN_IN
          NSPIN2=NSPIN2_IN
          ISFLIP=ISFLIP_IN
          IR_DIA=IR_DIA_IN
          NSTEP=NSTEP_IN
          ISPEED=ISPEED_IN
          IATTS=IATTS_IN
          ILPM=ILPM_IN
          XLMP0=XLPM0_IN
          DO JAT=1,NAT
             RTHFWD(JAT)=RTHFWD_IN(JAT)
             RTHBWD(JAT)=RTHBWD_IN(JAT)
             IBWD(JAT)=IBWD_IN(JAT)
             NUMAX(JAT)=NUMAX_IN(JAT)
             N_RA(JAT)=N_RA_IN(JAT)
          ENDDO
C
        ELSEIF(EXC.EQ.1) THEN
          NDIF=NDIF_EX
          NO=NO_EX
          ISPHER=ISPHER_EX
          IFWD=IFWD_EX
          NTHOUT=NTHOUT_EX
          IPW=IPW_EX
          NCUT=NCUT_EX
          PCTINT=PCTINT_EX
          IPP=IPP_EX
          ISPEED=ISPEED_EX
          IATTS=IATTS_EX
          ILENGTH=ILENGTH_EX
          RLENGTH=RLENGTH_EX
          I_NO=I_NO_EX
          I_RA=I_RA_EX
          ISPIN=ISPIN_EX
          IDICHR=IDICHR_EX
          NSPIN=NSPIN_EX
          NSPIN2=NSPIN2_EX
          ISFLIP=ISFLIP_EX
          IR_DIA=IR_DIA_EX
          NSTEP=NSTEP_EX
          ISPEED=ISPEED_EX
          IATTS=IATTS_EX
          ILPM=ILPM_EX
          XLMP0=XLPM0_EX
          DO JAT=1,NAT
             RTHFWD(JAT)=RTHFWD_EX(JAT)
             RTHBWD(JAT)=RTHBWD_EX(JAT)
             IBWD(JAT)=IBWD_EX(JAT)
             NUMAX(JAT)=NUMAX_EX(JAT)
             N_RA(JAT)=N_RA_EX(JAT)
          ENDDO
C
        ELSEIF(OUT1.EQ.1) THEN
          NDIF=NDIF_O1
          NO=NO_O1
          ISPHER=ISPHER_O1
          IFWD=IFWD_O1
          NTHOUT=NTHOUT_O1
          IPW=IPW_O1
          NCUT=NCUT_O1
          PCTINT=PCTINT_O1
          IPP=IPP_O1
          ISPEED=ISPEED_O1
          IATTS=IATTS_O1
          ILENGTH=ILENGTH_O1
          RLENGTH=RLENGTH_O1
          I_NO=I_NO_O1
          I_RA=I_RA_O1
          ISPIN=ISPIN_O1
          IDICHR=IDICHR_O1
          NSPIN=NSPIN_O1
          NSPIN2=NSPIN2_O1
          ISFLIP=ISFLIP_O1
          IR_DIA=IR_DIA_O1
          NSTEP=NSTEP_O1
          ISPEED=ISPEED_O1
          IATTS=IATTS_O1
          ILPM=ILPM_O1
          XLMP0=XLPM0_O1
          DO JAT=1,NAT
             RTHFWD(JAT)=RTHFWD_O1(JAT)
             RTHBWD(JAT)=RTHBWD_O1(JAT)
             IBWD(JAT)=IBWD_O1(JAT)
             NUMAX(JAT)=NUMAX_O1(JAT)
             N_RA(JAT)=N_RA_O1(JAT)
          ENDDO
        ENDIF
      ENDIF
C
C  Consistency of input files numbers with existing version
C
C    Cluster and adsorbate files
C
      
      IF(NEL.EQ.1) THEN
C
C  One electron spectroscopies (PhD,LEED,XAS,AED,REXS,...)
C
        INFILE4=INFILE2B
        INFILE5=INFILE3B
C
        IUI4=IUI2B
        IUI5=IUI3B
C
C    tl, rad and external directions files
C
        IF(INC.EQ.1) THEN
C
          INFILE2=INFILE4B
          INFILE3=INFILE5B
          INFILE6=INFILE6B
C
          IUI2=IUI4B
          IUI3=IUI5B
          IUI6=IUI6B
C
        ELSEIF(EXC.EQ.1) THEN
C
          INFILE2=INFILE7B
          INFILE3=INFILE8B
C
          IUI2=IUI7B
          IUI3=IUI8B
C
        ELSEIF(OUT1.EQ.1) THEN
C
          INFILE2=INFILE9B
          INFILE3=INFILE10B
          INFILE6=INFILE11B
C
          IUI2=IUI9B
          IUI3=IUI10B
          IUI6=IUI11B
C
        ENDIF
C
      ELSEIF(NEL.EQ.2) THEN
C
C  Case of APECS
C
        INFILE2=INFILE9B
        INFILE3=INFILE10B
        INFILE6=INFILE11B
C
        IUI2=IUI9B
        IUI3=IUI10B
        IUI6=IUI11B
C
        INFILE7=INFILE12B
        INFILE8=INFILE13B
        INFILE9=INFILE14B
C
        IUI7=IUI12B
        IUI8=IUI13B
        IUI9=IUI14B
C
      ENDIF
C
      RETURN
C
      END

C
C=======================================================================
C
      FUNCTION SIG2(RJ,JTYP)
C
C  This routine evaluates the mean square displacements.
C
      INCLUDE 'spec.inc'
C
      COMMON /DEBWAL/ IDCM,IDWSPH,TD,QD,T,RSJ,UJ2(NATM)
      COMMON /MASSAT/ XM(NATM)
      COMMON /RESEAU/ N1,N2,N3,N4,A0,R1,R2,UN
C
      REAL MJ
C
      CHARACTER*3 UN
C
      DATA COEF/145.52539/
      DATA RZ2,RZ4,RZ6/1.644934,1.082323,1.017343/
C
      A=TD/T
      BJ=QD*RJ
      U=BJ/A
      MJ=XM(JTYP)
      C=COEF/(2.*MJ*TD)
      COMP=RZ2-U*U*RZ4+U*U*U*U*RZ6
      X1=0.
      X2=0.
      X3=0.
      X4=0.
      DO 10 N=1,8
        Z=FLOAT(N)
        X1=X1+EXP(-Z*A)*((A/Z)+(1./(Z*Z)))
        X2=X2+1./(Z**8+U*U*(Z**6))
        X3=X3+EXP(-Z*A)*Z/(Z*Z+U*U)
        X4=X4+EXP(-Z*A)/(Z*Z+U*U)
  10  CONTINUE
      P1=1.+4.*(RZ2-X1)/(A*A)
      P2=-2.*(1.-COS(BJ))/(BJ*BJ)
      P3=-4.*(COMP-(U**6)*X2)/(A*A)
      P4=4.*SIN(BJ)*X3/(A*BJ)
      P5=4.*COS(BJ)*X4/(A*A)
      SIG2=C*(P1+P2+P3+P4+P5)/(A0*A0)
C
      RETURN
C
      END
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION SIXJ_IN(J1,J2,L1,L2,L3)
C
C  This function calculates the initial value {J1 J2 L1+L2}
C                                             {L1 L2   L3 }
C
C  A 6j symbol {J1 J2 J3} is non zero only if
C              {J4 J5 J6}
C
C   (J1,J2,J3),(J4,J5,J3),(J2,J4,J6) and (J1,J5,J6) satisfy the triangular inequality :
C
C       (a,b,c) non zero if |a-b| <= c <= (a+b) . This means also that (a+b) and c must
C       have the same nature (integer or half-integer).
C
C   (J1,J2,J3) and (J4,J5,J3) are taken care of by the bounds of J3, JJ_MIN and JJ_MAX,
C       as chosen in the N_J routine. Here we check the two last ones.
C
C                                        Last modified :  8 Dec 2008
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'spec.inc'
C
      REAL*8 J1,J2,L1,L2,L3
C
      COMMON /LOGAMAD/ GLD(0:N_GAUNT,2)
C
      DATA SMALL /0.0001/
C
      IZERO=0
C
C  Check for unphysical values of L3
C
      IF(DABS(J2-L1).GT.L3) IZERO=1
      IF(J2+L1.LT.L3) IZERO=1
      IF(IG(J2+L1).NE.IG(L3)) IZERO=1
      IF(DABS(J1-L2).GT.L3) IZERO=1
      IF(J1+L2.LT.L3) IZERO=1
      IF(IG(J1+L2).NE.IG(L3)) IZERO=1
C
      IF(IZERO.EQ.1) THEN
        SIXJ_IN=0.D0
      ELSE
C
C  Storage indices of the angular momenta.
C
        LJ1=INT(J1+SIGN(SMALL,J1))
        LJ2=INT(J2+SIGN(SMALL,J2))
        LL1=INT(L1+SIGN(SMALL,L1))
        LL2=INT(L2+SIGN(SMALL,L2))
        LL3=INT(L3+SIGN(SMALL,L3))
        LL1_2=INT(L1+L1+SIGN(SMALL,L1))
        LL2_2=INT(L2+L2+SIGN(SMALL,L2))
C
        MSIGN=INT(J1+J2+L1+L2+SIGN(SMALL,J1+J2+L1+L2))
        IF(MOD(MSIGN,2).EQ.0) THEN
          SIGNE=1.D0
        ELSE
          SIGNE=-1.D0
        ENDIF
C
        D1=GLD(LL1_2+1,1) + GLD(LL2_2+1,1) - GLD(LL1_2+LL2_2+2,1)
        D2=GLD(INT(J1+J2+L1+L2)+2,IG(J1+J2+L1+L2)) -
     1   GLD(INT(J1+J2-L1-L2)+1,IG(J1+J2-L1-L2))
        D3=GLD(INT(J1-J2+L1+L2)+1,IG(J1-J2+L1+L2)) -
     1     GLD(INT(J1+L2-L3)+1,IG(J1+L2-L3))
        D4=GLD(INT(J2-J1+L1+L2)+1,IG(J2-J1+L1+L2)) -
     1     GLD(INT(-J1+L2+L3)+1,IG(-J1+L2+L3))
        D5=GLD(INT(J1-L2+L3)+1,IG(J1-L2+L3)) -
     1     GLD(INT(J1+L2+L3)+2,IG(J1+L2+L3))
        D6=GLD(INT(J2+L3-L1)+1,IG(J2+L3-L1)) -
     1     GLD(INT(J2-L3+L1)+1,IG(J2-L3+L1))
        D7=GLD(INT(L1+L3-J2)+1,IG(L1+L3-J2)) +
     1     GLD(INT(L1+J2+L3)+2,IG(L1+J2+L3))
C
        SIXJ_IN=SIGNE*DSQRT(DEXP(D1+D2+D3+D4+D5+D6-D7))
C
      ENDIF
C
      END
C
C=======================================================================
C
      SUBROUTINE SPH_HAR(NL,X,CF,YLM,NC)
C
C  This routine computes the complex spherical harmonics using Condon and
C                  Shortley phase convention.
C
C  If the angular direction R=(THETAR,PHIR)  is given in cartesian
C      coordinates by (XR,YR,ZR), the arguments of the subroutine are :
C
C                    X  = ZR         = cos(THETAR)
C                    CF = XR + i YR  = sin(THETAR)*exp(i PHIR)
C
C          NL is the dimensioning of the YLM array and NC is
C                   the maximum l value to be computed.
C
      INCLUDE 'spec.inc'
C
      COMMON /EXPFAC2/ EXPF2(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /FACTSQ/ FSQ(0:2*NL_M-2)
C
      COMPLEX YLM(0:NL,-NL:NL),COEF,YMM,YMMP,C,CF
C
      DATA SQ4PI_INV,SQR3_INV /0.282095,0.488602/
C
C
      YLM(0,0)=CMPLX(SQ4PI_INV)
      YLM(1,0)=X*SQR3_INV
      DO L=2,NC
        Y=1./FLOAT(L)
        YLM(L,0)=X*SQRT(4.-Y*Y)*YLM(L-1,0) -
     1           (1.-Y)*SQRT(1.+2./(FLOAT(L)-1.5))*YLM(L-2,0)
      ENDDO
C
      C2=-1.
      C=-0.5*CF
C
      C1=1.
      COEF=(1.,0.)
      DO M=1,NC
        C1=C1*C2
        COEF=COEF*C
        YMM=SQ4PI_INV*COEF*FSQ(M)
        YLM(M,M)=YMM
        YLM(M,-M)=C1*CONJG(YMM)
        YMMP=X*SQRT(FLOAT(M+M+3))*YMM
        YLM(M+1,M)=YMMP
        YLM(M+1,-M)=C1*CONJG(YMMP)
        IF(M.LT.NC-1) THEN
          DO L=M+2,NC
            YLM(L,M)=(X*(L+L-1)*EXPF2(L-1,M)*YLM(L-1,M) -
     1                (L+M-1)*EXPF2(L-2,M)*YLM(L-2,M))/(EXPF2(L,M)*
     2                (L-M))
            YLM(L,-M)=C1*CONJG(YLM(L,M))
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE SPH_HAR2(NL,X,CF,YLM,NC)
C
C  This routine computes the complex spherical harmonics using Condon and
C                  Shortley phase convention.
C
C  If the angular direction R=(THETAR,PHIR)  is given in cartesian
C      coordinates by (XR,YR,ZR), the arguments of the subroutine are :
C
C                    X  = ZR         = cos(THETAR)
C                    CF = XR + i YR  = sin(THETAR)*exp(i PHIR)
C
C          NL is the dimensioning of the YLM array and NC is
C                   the maximum l value to be computed.
C
C  This is the double precision version of sph_har.f
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INCLUDE 'spec.inc'
C
      COMMON /DEXPFAC2/ DEXPF2(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /DFACTSQ/ DFSQ(0:2*NL_M-2)
C
      REAL*8 DEXPF2,DFSQ
C
      COMPLEX*16 YLM(0:NL,-NL:NL),COEF,YMM,YMMP,C,CF
C
      DATA SQ4PI_INV,SQR3_INV /0.282094791774D0,0.488602511903D0/
C
C
      YLM(0,0)=DCMPLX(SQ4PI_INV)
      YLM(1,0)=X*SQR3_INV
      DO L=2,NC
        Y=1.D0/DFLOAT(L)
        YLM(L,0)=X*DSQRT(4.D0-Y*Y)*YLM(L-1,0) -
     1           (1.D0-Y)*DSQRT(1.D0+2.D0/(DFLOAT(L)-1.5D0))*YLM(L-2,0)
      ENDDO
C
      C2=-1.D0
      C=-0.5D0*CF
C
      C1=1.D0
      COEF=(1.D0,0.D0)
      DO M=1,NC
        C1=C1*C2
        COEF=COEF*C
        YMM=SQ4PI_INV*COEF*DFSQ(M)
        YLM(M,M)=YMM
        YLM(M,-M)=C1*DCONJG(YMM)
        YMMP=X*DSQRT(DFLOAT(M+M+3))*YMM
        YLM(M+1,M)=YMMP
        YLM(M+1,-M)=C1*DCONJG(YMMP)
        IF(M.LT.NC-1) THEN
          DO L=M+2,NC
            YLM(L,M)=(X*DFLOAT(L+L-1)*DEXPF2(L-1,M)*YLM(L-1,M) -
     1                DFLOAT(L+M-1)*DEXPF2(L-2,M)*YLM(L-2,M))/
     2                (DEXPF2(L,M)*DFLOAT(L-M))
            YLM(L,-M)=C1*DCONJG(YLM(L,M))
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE STOP_EXT(I_EXT,I_EXT_A,SPECTRO)
C
C  This routine stops the code when the dimension N_TILT_M in the
C      spec.inc file is insufficient for the number of values to
C      Gaussian average over (as generated by the ext_dir.f code)
C
      INCLUDE 'spec.inc'
C
      COMMON /INFILES/ INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6,
     1                 INFILE7,INFILE8,INFILE9
      COMMON /INUNITS/ IUI1,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,IUI9
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
      CHARACTER*24 INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6
      CHARACTER*24 INFILE7,INFILE8,INFILE9
      CHARACTER*3 SPECTRO
C
      NSET=1
      NSET_A=1
C
      IF((SPECTRO.EQ.'PHD').OR.(SPECTRO.EQ.'AED')) THEN
        IF(I_EXT.EQ.-1) THEN
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          READ(IUI6,15) IDIR,NSET
          CLOSE(IUI6)
        ENDIF
        IF(I_EXT_A.EQ.-1) THEN
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          READ(IUI6,15) IDIR,NSET_A
          CLOSE(IUI6)
        ENDIF
      ENDIF
      IF(SPECTRO.EQ.'APC') THEN
        IF(I_EXT.EQ.-1) THEN
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          READ(IUI6,15) IDIR,NSET
          CLOSE(IUI6)
        ENDIF
        IF(I_EXT_A.EQ.-1) THEN
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
          READ(IUI9,15) IDIR,NSET_A
          CLOSE(IUI9)
        ENDIF
      ENDIF
C
      IF(MAX(NSET,NSET_A).GT.N_TILT_M) THEN
        WRITE(IUO1,10) MAX(NSET,NSET_A)
        STOP
      ENDIF
C
  10  FORMAT(///,16X,'<<<<<<<<<<  N_TILT_M SHOULD BE AT LEAST ',
     1       I3,'  >>>>>>>>>>')
  15  FORMAT(6X,I1,1X,I3)
C
      RETURN
C
      END

C
C=======================================================================
C
      SUBROUTINE STOP_TREAT(NFICHLEC,NPLAN,NEMET,NE,NTHETA,NTHETA_A,
     1                      NPHI,NPHI_A,ISOM,I_EXT,I_EXT_A,SPECTRO)
C
C   This subroutine stops the code before the long MS calculations
C     when the dimensioning NDIM_M of the treatment routines
C     (treat_aed,treat_apc,treat_phd,treat_xas) is insufficient.
C
C
C                                         Last modified : 06 Oct 2006
C
      INCLUDE 'spec.inc'
C
      CHARACTER*3 SPECTRO
C
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
      IF(ISOM.EQ.0) THEN
C
C   Photoelectron diffraction case
C
        IF(SPECTRO.EQ.'PHD') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            NDP=NEMET*NTHETA*NPHI*NE
          ELSEIF(I_EXT.EQ.-1) THEN
            NDP=NEMET*NTHETA*NPHI*NE*2
          ELSEIF(I_EXT.EQ.2) THEN
            NDP=NEMET*NTHETA*NE
          ENDIF
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 10
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
C
C   Auger electron diffraction case
C
        ELSEIF(SPECTRO.EQ.'AED') THEN
          IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
            NDP=NEMET*NTHETA_A*NPHI_A*NE
          ELSEIF(I_EXT_A.EQ.-1) THEN
            NDP=NEMET*NTHETA_A*NPHI_A*NE*2
          ELSEIF(I_EXT_A.EQ.2) THEN
            NDP=NEMET*NTHETA_A*NE
          ENDIF
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 20
          IF((NTHETA_A.GT.NTH_M).OR.(NPHI_A.GT.NPH_M)) GOTO 50
C
C   X-ray absorption case
C
        ELSEIF(SPECTRO.EQ.'XAS') THEN
          NDP=NEMET*NE
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 30
C
C   Auger Photoelectron coincidence spectroscopy case
C
        ELSEIF(SPECTRO.EQ.'APC') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NEMET*NTHETA*NPHI*NE*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NEMET*NTHETA*NPHI*NE*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NEMET*NTHETA*NPHI*NE*NTHETA_A
            ENDIF
          ELSEIF(I_EXT.EQ.-1) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NEMET*NTHETA*NPHI*NE*2*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NEMET*NTHETA*NPHI*NE*2*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NEMET*NTHETA*NPHI*NE*2*NTHETA_A
            ENDIF
          ELSEIF(I_EXT.EQ.2) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NEMET*NTHETA*NE*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NEMET*NTHETA*NE*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NEMET*NTHETA*NE*NTHETA_A
            ENDIF
          ENDIF
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 40
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
          IF((NTHETA_A.GT.NTH_M).OR.(NPHI_A.GT.NPH_M)) GOTO 50
C
C   Resonant elestic X-ray scattering case
C
        ELSEIF(SPECTRO.EQ.'RES') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            NDP=NEMET*NTHETA*NPHI*NE
          ELSEIF(I_EXT.EQ.-1) THEN
            NDP=NEMET*NTHETA*NPHI*NE*2
          ELSEIF(I_EXT.EQ.2) THEN
            NDP=NEMET*NTHETA*NE
          ENDIF
          NTT=NPLAN*NDP
          IF(NTT.GT.NDIM_M) GOTO 10
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
        ENDIF
C
      ELSE
C
C   Photoelectron diffraction case
C
        IF(SPECTRO.EQ.'PHD') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            NDP=NTHETA*NPHI*NE
          ELSEIF(I_EXT.EQ.-1) THEN
            NDP=NTHETA*NPHI*NE*2
          ELSEIF(I_EXT.EQ.2) THEN
            NDP=NTHETA*NE
          ENDIF
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 10
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
C
C   Auger electron diffraction case
C
        ELSEIF(SPECTRO.EQ.'AED') THEN
          IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
            NDP=NTHETA_A*NPHI_A*NE
          ELSEIF(I_EXT_A.EQ.-1) THEN
            NDP=NTHETA_A*NPHI_A*NE*2
          ELSEIF(I_EXT_A.EQ.2) THEN
            NDP=NTHETA_A*NE
          ENDIF
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 20
          IF((NTHETA_A.GT.NTH_M).OR.(NPHI_A.GT.NPH_M)) GOTO 50
C
C   X-ray absorption case
C
        ELSEIF(SPECTRO.EQ.'XAS') THEN
          NDP=NE
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 30
C
C   Auger Photoelectron coincidence spectroscopy case
C
        ELSEIF(SPECTRO.EQ.'APC') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NTHETA*NPHI*NE*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NTHETA*NPHI*NE*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NTHETA*NPHI*NE*NTHETA_A
            ENDIF
          ELSEIF(I_EXT.EQ.-1) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NTHETA*NPHI*NE*2*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NTHETA*NPHI*NE*2*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NTHETA*NPHI*NE*2*NTHETA_A
            ENDIF
          ELSEIF(I_EXT.EQ.2) THEN
            IF((I_EXT_A.EQ.0).OR.(I_EXT_A.EQ.1)) THEN
              NDP=NTHETA*NE*NTHETA_A*NPHI_A
            ELSEIF(I_EXT_A.EQ.-1) THEN
              NDP=NTHETA*NE*NTHETA_A*NPHI_A*2
            ELSEIF(I_EXT_A.EQ.2) THEN
              NDP=NTHETA*NE*NTHETA_A
            ENDIF
          ENDIF
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 40
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
          IF((NTHETA_A.GT.NTH_M).OR.(NPHI_A.GT.NPH_M)) GOTO 50
C
C   Resonant elestic X-ray scattering case
C
        ELSEIF(SPECTRO.EQ.'RES') THEN
          IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
            NDP=NTHETA*NPHI*NE
          ELSEIF(I_EXT.EQ.-1) THEN
            NDP=NTHETA*NPHI*NE*2
          ELSEIF(I_EXT.EQ.2) THEN
            NDP=NTHETA*NE
          ENDIF
          NTT=NFICHLEC*NDP
          IF(NTT.GT.NDIM_M) GOTO 10
          IF((NTHETA.GT.NTH_M).OR.(NPHI.GT.NPH_M)) GOTO 50
        ENDIF
      ENDIF
C
      GOTO 5
C
  10  WRITE(IUO1,11) NTT
      STOP
  20  WRITE(IUO1,21) NTT
      STOP
  30  WRITE(IUO1,31) NTT
      STOP
  40  WRITE(IUO1,41) NTT
      STOP
  50  WRITE(IUO1,51) MAX(NTHETA,NPHI,NTHETA_A,NPHI_A)
      STOP
C
  11  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',
     1       ' IN THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ',
     2       'FOR THE TREAT_PHD SUBROUTINE   >>>>>>>>>>',/,8X,
     3       '<<<<<<<<<<          NDIM_M BE AT LEAST ',I8,
     4       '         >>>>>>>>>>')
  21  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',
     1       ' IN THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ',
     2       'FOR THE TREAT_AED SUBROUTINE   >>>>>>>>>>',/,8X,
     3       '<<<<<<<<<<          NDIM_M BE AT LEAST ',I8,
     4       '         >>>>>>>>>>')
  31  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',
     1       ' IN THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ',
     2       'FOR THE TREAT_XAS SUBROUTINE   >>>>>>>>>>',/,8X,
     3       '<<<<<<<<<<          NDIM_M BE AT LEAST ',I8,
     4       '         >>>>>>>>>>')
  41  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',
     1       ' IN THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ',
     2       'FOR THE TREAT_APC SUBROUTINE   >>>>>>>>>>',/,8X,
     3       '<<<<<<<<<<          NDIM_M BE AT LEAST ',I8,
     4       '         >>>>>>>>>>')
  51  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF NTH_M OR NPH_M TOO SMALL',
     1       'IN THE INCLUDE FILE - SHOULD BE AT LEAST ',I6,
     2       '  >>>>>>>>>>')
C
   5  RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE SUP_ZEROS(TL,LMAX,NE,NAT,IUO1,ITRTL)
C
C   This routine suppresses possible zeros in the TL arrays so that
C     the code runs faster because of lower values of LMAX. Actually,
C     the TL array is not modified, it is just the LMAX array that is
C     altered. This is particularly useful for energy variations or
C     for matrix inversion
C
      INCLUDE 'spec.inc'
C
      COMPLEX TL_,TL(0:NT_M,4,NATM,NE_M)
C
      INTEGER LMAX(NATM,NE_M)
C
      IF(ITRTL.EQ.1) THEN
        SMALL=0.1
      ELSEIF(ITRTL.EQ.2) THEN
        SMALL=0.01
      ELSEIF(ITRTL.EQ.3) THEN
        SMALL=0.001
      ELSEIF(ITRTL.EQ.4) THEN
        SMALL=0.0001
      ELSEIF(ITRTL.EQ.5) THEN
        SMALL=0.00001
      ELSEIF(ITRTL.EQ.6) THEN
        SMALL=0.000001
      ELSEIF(ITRTL.EQ.7) THEN
        SMALL=0.0000001
      ELSEIF(ITRTL.EQ.8) THEN
        SMALL=0.00000001
      ELSE
        ITRTL=9
        SMALL=0.000000001
      ENDIF
C
      WRITE(IUO1,10)
      WRITE(IUO1,15) ITRTL
C
      DO JE=1,NE
        WRITE(IUO1,20) JE
        DO JAT=1,NAT
          NONZERO=0
          LM=LMAX(JAT,JE)
          DO L=0,LM
            TL_=TL(L,1,JAT,JE)
            IF((ABS(REAL(TL_)).GE.SMALL).OR.
     1         (ABS(AIMAG(TL_)).GE.SMALL)) THEN
              NONZERO=NONZERO+1
            ENDIF
          ENDDO
          LMAX(JAT,JE)=NONZERO-1
          WRITE(IUO1,30) JAT,LM,NONZERO-1
        ENDDO
      ENDDO
C
      WRITE(IUO1,40)
C
  10  FORMAT(//,'   ---> CHECK FOR ZEROS IN THE TL FILE TO REDUCE',
     1          ' THE AMOUNT OF COMPUTING :',/)
  15  FORMAT(/,'  (ONLY THE MATRIX ELEMENTS NON ZERO ',
     1         'TO THE FIRST ',I1,' DECIMAL DIGITS ARE KEPT)',/)
  20  FORMAT(/,15X,'ENERGY POINT No ',I3,/)
  30  FORMAT(8X,'PROTOTYPICAL ATOM No ',I5,'  INITIAL LMAX = ',I2,
     1            '   FINAL LMAX = ',I2)
  40  FORMAT(//)
C
      RETURN
C
      END
C
C=======================================================================
C
      FUNCTION UJ_SQ(JTYP)
C
C  This routine evaluates the mean square displacements UJ_SQ,
C    first along une direction (x, y or z): UJ2 within the Debye model,
C              using the Debye function formulation
C
C  X1 is the Debye function phi_1
C  UJ_SQ is given in unit of the square of the lattice parameter A0
C  Temperatures are expressed in Kelvin
C
C  The coefficient COEF equals:
C
C           3 hbar^{2} N_A 10^{3} / (4 k_B)
C
C    where N_A is the Avogadro number, k_B is Boltzmann's constant
C    and 10^3 arises from the fact that the atomic mass is
C    expressed in grams
C
C  Then UJ_SQ is obtained as UJ_SQ = (2 + RSJ) UJJ for surface atoms
C                            UJ_SQ = 3 UJJ for bulk atoms
C
C
C  For empty spheres, two possibilities are provided. By construction,
C    they are very light (their mass is taken as 1/1836 of the mass
C    of a H atom) and so they will vibrate a lot (IDCM = 1). When
C    setting IDCM = 2, their mean square displacement is set to a
C    tiny value so that they hardly vibrate (frozen empty spheres)
C
C                                        Last modified : 25 Apr 2013
C
      INCLUDE 'spec.inc'
C
      COMMON /DEBWAL/ IDCM,IDWSPH,TD,QD,T,RSJ,UJ2(NATM)
      COMMON /MASSAT/ XM(NATM)
      COMMON /RESEAU/ N1,N2,N3,N4,A0,R1,R2,UN
      COMMON /VIBRAT/ I_FREE(NATP_M)
C
      REAL MJ
C
      CHARACTER*3 UN
C
      DATA COEF /36.381551/
      DATA RZ2  /1.644934/
C
      N_MAX=20
C
C  Computation of the 1D mean square displacement UJ2
C
      A=TD/T
      MJ=XM(JTYP)
      C=COEF/(MJ*TD)
      X1=0.
C
      DO N=1,N_MAX
        Z=FLOAT(N)
        X1=X1+EXP(-Z*A)*(A+1./Z)/Z
      ENDDO
C
      P1=1.+4.*(RZ2-X1)/(A*A)
      UJJ=C*P1/(A0*A0)
C
C  3D mean square displacement UJ_SQ
C
      IF(IDCM.EQ.1) THEN
        UJ_SQ=(3.+FLOAT(I_FREE(JTYP))*(RSJ-1.))*UJJ
      ELSEIF(IDCM.EQ.2) THEN
        UJ_SQ=1.0E-20
      ENDIF
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE DIRAN(VINT,ECIN,J_EL)
C
C  This subroutine calculates the direction(s) of the analyzer with
C             or without an angular averaging.
C
C             DIRANA is the internal direction
C             ANADIR is the external direction
C
C             J_EL is the type of electron : 1 ---> photoelectron
C                                            2 ---> Auger electron
C
C                                            Last modified : 23/03/2006
C
      COMPLEX COEF,IC
C
      COMMON /DIRECT/ DIRANA(3,49),ANADIR(3,49),RTHEXT,
     1                RPHI,THETAR(49),PHIR(49)
      COMMON /DIRECT_A/ DIRANA_A(3,49),ANADIR_A(3,49),RTHEXT_A,
     1                  RPHI_A,THETAR_A(49),PHIR_A(49)
      COMMON /MOYEN/ IMOY,NDIR,ACCEPT,ICHKDIR
      COMMON /MOYEN_A/ IMOY_A,NDIR_A,ACCEPT_A,ICHKDIR_A
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /TESTS/ I1,IPRINT,I2,N1,I3
      COMMON /TYPCAL/ IPHI,IE,ITHETA,IFTHET,IMOD,IPOL,I_CP,I_EXT,I_TEST
      COMMON /TYPCAL_A/ IPHI_A,IE_A,ITHETA_A,IFTHET_A,IMOD_A,I_CP_A,
     1                  I_EXT_A,I_TEST_A
C
      DATA PI,PIS2,PIS180 /3.141593,1.570796,0.017453/
C
      IC=(0.,1.)
C
      IF(J_EL.EQ.1) THEN
        ANADIR(1,1)=SIN(RTHEXT)*COS(RPHI)
        ANADIR(2,1)=SIN(RTHEXT)*SIN(RPHI)
        ANADIR(3,1)=COS(RTHEXT)
        IF((ABS(I_EXT).LE.1).AND.(I_TEST.NE.2)) THEN
          CALL REFRAC(VINT,ECIN,RTHEXT,RTHINT)
        ELSE
          RTHINT=RTHEXT
        ENDIF
        IF((IPRINT.GT.0).AND.(I_EXT.NE.2)) THEN
          DTHEXT=RTHEXT/PIS180
          DTHINT=RTHINT/PIS180
          IF(I_TEST.NE.2) WRITE(IUO1,20) DTHEXT,DTHINT
        ENDIF
        DIRANA(1,1)=SIN(RTHINT)*COS(RPHI)
        DIRANA(2,1)=SIN(RTHINT)*SIN(RPHI)
        DIRANA(3,1)=COS(RTHINT)
        THETAR(1)=RTHINT
        PHIR(1)=RPHI
C
C The change in the definition below is necessary as RPHI is
C   used to define the rotation axis of the direction of the detector
C   when doing polar variations
C
        IF(ITHETA.EQ.1) THEN
          IF(RPHI.GT.PIS2) THEN
            RTHEXT=-RTHEXT
            RPHI=RPHI-PI
          ELSEIF(RPHI.LT.-PIS2) THEN
            RTHEXT=-RTHEXT
            RPHI=RPHI+PI
          ENDIF
        ENDIF
C
        IF(IMOY.GE.1) THEN
          N=2**(IMOY-1)
          S=SIN(ACCEPT*PI/180.)
          RN=FLOAT(N)
          J=1
          DO K1=-N,N
            RK1=FLOAT(K1)
            DO K2=-N,N
              RK2=FLOAT(K2)
              D=SQRT(RK1*RK1+RK2*RK2)
              IF((D-RN).GT.0.000001) GOTO 10
              IF((K1.EQ.0).AND.(K2.EQ.0)) GOTO 10
              C=SQRT(RN*RN-(RK1*RK1+RK2*RK2)*S*S)
              J=J+1
C
              ANADIR(1,J)=(RK1*S*COS(RTHEXT)*COS(RPHI)
     1                      -RK2*S*SIN(RPHI)+C*ANADIR(1,1))/RN
              ANADIR(2,J)=(RK1*S*COS(RTHEXT)*SIN(RPHI)
     1                      +RK2*S*COS(RPHI)+C*ANADIR(2,1))/RN
              ANADIR(3,J)=(-RK1*S*SIN(RTHEXT) +C*ANADIR(3,1))/RN
              THETA_R=ACOS(ANADIR(3,J))
              COEF=ANADIR(1,J)+IC*ANADIR(2,J)
              CALL ARCSIN(COEF,ANADIR(3,J),PHI_R)
              IF((ABS(I_EXT).LE.1).AND.(I_TEST.NE.2)) THEN
                CALL REFRAC(VINT,ECIN,THETA_R,THINT_R)
              ELSE
                THINT_R=THETA_R
              ENDIF
C
              DIRANA(1,J)=SIN(THINT_R)*COS(PHI_R)
              DIRANA(2,J)=SIN(THINT_R)*SIN(PHI_R)
              DIRANA(3,J)=COS(THINT_R)
C
              THETAR(J)=THINT_R
              PHIR(J)=PHI_R
  10          CONTINUE
            ENDDO
          ENDDO
        ENDIF
C
      ELSEIF(J_EL.EQ.2) THEN
        ANADIR_A(1,1)=SIN(RTHEXT_A)*COS(RPHI_A)
        ANADIR_A(2,1)=SIN(RTHEXT_A)*SIN(RPHI_A)
        ANADIR_A(3,1)=COS(RTHEXT_A)
        IF((ABS(I_EXT_A).LE.1).AND.(I_TEST_A.NE.2)) THEN
          CALL REFRAC(VINT,ECIN,RTHEXT_A,RTHINT_A)
        ELSE
          RTHINT_A=RTHEXT_A
        ENDIF
        IF((IPRINT.GT.0).AND.(I_EXT_A.NE.2)) THEN
          DTHEXT_A=RTHEXT_A/PIS180
          DTHINT_A=RTHINT_A/PIS180
          IF(I_TEST_A.NE.2) WRITE(IUO1,21) DTHEXT_A,DTHINT_A
        ENDIF
        DIRANA_A(1,1)=SIN(RTHINT_A)*COS(RPHI_A)
        DIRANA_A(2,1)=SIN(RTHINT_A)*SIN(RPHI_A)
        DIRANA_A(3,1)=COS(RTHINT_A)
        THETAR_A(1)=RTHINT_A
        PHIR_A(1)=RPHI_A
C
C The change in the definition below is necessary as RPHI is
C   used to define the rotation axis of the direction of the detector
C   when doing polar variations
C
        IF(ITHETA_A.EQ.1) THEN
          IF(RPHI_A.GT.PIS2) THEN
            RTHEXT_A=-RTHEXT_A
            RPHI_A=RPHI_A-PI
          ELSEIF(RPHI_A.LT.-PIS2) THEN
            RTHEXT_A=-RTHEXT_A
            RPHI_A=RPHI_A+PI
          ENDIF
        ENDIF
C
        IF(IMOY_A.GE.1) THEN
          N=2**(IMOY_A-1)
          S=SIN(ACCEPT_A*PI/180.)
          RN=FLOAT(N)
          J=1
          DO K1=-N,N
            RK1=FLOAT(K1)
            DO K2=-N,N
              RK2=FLOAT(K2)
              D=SQRT(RK1*RK1+RK2*RK2)
              IF((D-RN).GT.0.000001) GOTO 15
              IF((K1.EQ.0).AND.(K2.EQ.0)) GOTO 15
              C=SQRT(RN*RN-(RK1*RK1+RK2*RK2)*S*S)
              J=J+1
C
              ANADIR_A(1,J)=(RK1*S*COS(RTHEXT_A)*COS(RPHI_A)
     1                      -RK2*S*SIN(RPHI_A)+C*ANADIR_A(1,1))/RN
              ANADIR_A(2,J)=(RK1*S*COS(RTHEXT_A)*SIN(RPHI_A)
     1                      +RK2*S*COS(RPHI_A)+C*ANADIR_A(2,1))/RN
              ANADIR_A(3,J)=(-RK1*S*SIN(RTHEXT_A) +C*ANADIR_A(3,1))/RN
              THETA_R_A=ACOS(ANADIR_A(3,J))
              COEF=ANADIR_A(1,J)+IC*ANADIR_A(2,J)
              CALL ARCSIN(COEF,ANADIR_A(3,J),PHI_R_A)
              IF((ABS(I_EXT_A).LE.1).AND.(I_TEST_A.NE.2)) THEN
                CALL REFRAC(VINT,ECIN,THETA_R_A,THINT_R_A)
              ELSE
                THINT_R_A=THETA_R_A
              ENDIF
C
              DIRANA_A(1,J)=SIN(THINT_R_A)*COS(PHI_R_A)
              DIRANA_A(2,J)=SIN(THINT_R_A)*SIN(PHI_R_A)
              DIRANA_A(3,J)=COS(THINT_R_A)
C
              THETAR_A(J)=THINT_R_A
              PHIR_A(J)=PHI_R_A
  15          CONTINUE
            ENDDO
          ENDDO
        ENDIF
C
      ENDIF
C
  20  FORMAT(/,10X,'PHOTOELECTRON EXTERNAL THETA  =',F7.2,5X,
     1        'INTERNAL THETA =', F7.2)
  21  FORMAT(/,10X,'AUGER ELECTRON EXTERNAL THETA =',F7.2,5X,
     1        'INTERNAL THETA =', F7.2)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE REFRAC(VINT,EKIN,RTHETA,RTHINT)
C
C  This routine calculates the refraction of a plane wave beam induced
C     by the surface potential barrier VINT. EKIN is the kinetic energy
C     outside the crystal.
C
C                                          Last modified :  3 Dec 2008
C
      DATA PIS180,SMALL /0.017453,0.001/
C
      IF(VINT.LT.0.) VINT=ABS(VINT)
      IF(ABS(VINT).LT.SMALL) THEN
        RTHINT=RTHETA
      ELSE
        U=VINT/(EKIN+VINT)
        DTHETA=RTHETA/PIS180
        REFRA=SIN(RTHETA)*SIN(RTHETA)*(1.-U)
        RTHINT=ASIN(SQRT(REFRA))
        IF(DTHETA.LT.0.) THEN
          RTHINT=-RTHINT
        ENDIF
      ENDIF
C
      RETURN
C
      END

C
C=======================================================================
C
      SUBROUTINE WEIGHT_SUM(ISOM,I_EXT,I_EXT_A,JEL)
C
C   This subroutine performs a weighted sum of the results
C     corresponding to different directions of the detector.
C     The directions and weights are read from an external input file
C
C   JEL is the electron undetected (i.e. for which the outgoing
C     directions are integrated over the unit sphere). It is always
C     1 for one electron spectroscopies (PHD). For APECS, It can be
C     1 (photoelectron) or 2 (Auger electron) or even 0 (no electron
C     detected)
C
C                                         Last modified : 31 Jan 2007
C
      INCLUDE 'spec.inc'
C
      PARAMETER(N_MAX=5810,NPM=20)
C
      REAL*4 W(N_MAX),W_A(N_MAX),ECIN(NE_M)
      REAL*4 DTHETA(N_MAX),DPHI(N_MAX),DTHETAA(N_MAX),DPHIA(N_MAX)
      REAL*4 SR_1,SF_1,SR_2,SF_2
      REAL*4 SUMR_1(NPM,NE_M,N_MAX),SUMR_2(NPM,NE_M,N_MAX)
      REAL*4 SUMF_1(NPM,NE_M,N_MAX),SUMF_2(NPM,NE_M,N_MAX)
C
      CHARACTER*3 SPECTRO,SPECTRO2
      CHARACTER*5 LIKE
      CHARACTER*13 OUTDATA
      CHARACTER*24 INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6
      CHARACTER*24 INFILE7,INFILE8,INFILE9
C
      COMMON /INFILES/ INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6,
     1                 INFILE7,INFILE8,INFILE9
      COMMON /INUNITS/ IUI1,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,IUI9
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
C
C
      DATA JVOL,JTOT/0,-1/
      DATA LIKE /'-like'/
C
      REWIND IUO2
C
      READ(IUO2,15) SPECTRO,OUTDATA
      IF(SPECTRO.NE.'APC') THEN
        READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        READ(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM
        SPECTRO2='XAS'
      ELSE
        READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        READ(IUO2,9) ISPIN_A,IDICHR_A,I_SO_A,ISFLIP_A,ICHKDIR_A,IPHI_A,
     1               ITHETA_A,IE_A
        READ(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM
        READ(IUO2,8) NPHI_A,NTHETA_A
        IF(JEL.EQ.1) THEN
          SPECTRO2='AED'
        ELSEIF(JEL.EQ.2) THEN
          SPECTRO2='PHD'
        ELSEIF(JEL.EQ.0) THEN
          SPECTRO2='XAS'
        ENDIF
      ENDIF
C
      IF(NPLAN.GT.NPM) THEN
        WRITE(IUO1,4) NPLAN+2
        STOP
      ENDIF
C
C  Reading the number of angular points
C
      IF(SPECTRO.NE.'APC') THEN
        OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
        READ(IUI6,1) N_POINTS
        READ(IUI6,5) I_DIM,N_DUM1,N_DUM2
        N_POINTS_A=1
      ELSE
        IF(JEL.EQ.1) THEN
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          READ(IUI6,1) N_POINTS
          READ(IUI6,5) I_DIM,N_DUM1,N_DUM2
          IF(I_EXT_A.EQ.0) THEN
            N_POINTS_A=NTHETA_A*NPHI_A
          ELSE
            OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
            READ(IUI9,1) N_POINTS_A
            READ(IUI9,5) I_DIM,N_DUM1,N_DUM2
          ENDIF
          NTHETA0=NTHETA_A
          NPHI0=NPHI_A
        ELSEIF(JEL.EQ.2) THEN
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
          READ(IUI9,1) N_POINTS_A
          READ(IUI9,5) I_DIM,N_DUM1,N_DUM2
          IF(I_EXT.EQ.0) THEN
            N_POINTS=NTHETA*NPHI
          ELSE
            OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
            READ(IUI6,1) N_POINTS
            READ(IUI6,5) I_DIM,N_DUM1,N_DUM2
          ENDIF
          NTHETA0=NTHETA
          NPHI0=NPHI
        ELSEIF(JEL.EQ.0) THEN
          OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')
          READ(IUI6,1) N_POINTS
          READ(IUI9,1) N_POINTS_A
          READ(IUI6,5) I_DIM,N_DUM1,N_DUM2
          READ(IUI9,5) I_DIM,N_DUM1,N_DUM2
        ENDIF
      ENDIF
C
      IF(SPECTRO.NE.'APC') THEN
        NANGLE=1
      ELSE
        IF(JEL.EQ.1) THEN
          NANGLE=N_POINTS_A
        ELSEIF(JEL.EQ.2) THEN
          NANGLE=N_POINTS
        ELSEIF(JEL.EQ.0) THEN
          NANGLE=1
        ENDIF
      ENDIF
C
C  Initialization of the arrays
C
      DO JE=1,NE
        DO JANGLE=1,NANGLE
          DO JPLAN=1,NPLAN+2
            SUMR_1(JPLAN,JE,JANGLE)=0.
            SUMF_1(JPLAN,JE,JANGLE)=0.
            IF(IDICHR.GT.0) THEN
              SUMR_2(JPLAN,JE,JANGLE)=0.
              SUMF_2(JPLAN,JE,JANGLE)=0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
C  Reading of the data to be angle integrated
C
      DO JE=1,NE
C
        DO JANGLE=1,N_POINTS
          IF(I_EXT.NE.0) READ(IUI6,2) TH,PH,W(JANGLE)
          DO JANGLE_A=1,N_POINTS_A
            IF((I_EXT_A.NE.0).AND.(JANGLE.EQ.1)) THEN
              READ(IUI9,2) THA,PHA,W_A(JANGLE_A)
            ENDIF
C
            DO JPLAN=1,NPLAN+2
C
              IF(IDICHR.EQ.0) THEN
                IF(SPECTRO.NE.'APC') THEN
                  READ(IUO2,3) JDUM,DTHETA(JANGLE),DPHI(JANGLE),
     1                         ECIN(JE),SR_1,SF_1
                ELSE
                  READ(IUO2,13) JDUM,DTHETA(JANGLE),DPHI(JANGLE),
     1                          ECIN(JE),DTHETAA(JANGLE_A),
     2                          DPHIA(JANGLE_A),SR_1,SF_1
                ENDIF
              ELSE
                IF(SPECTRO.NE.'APC') THEN
                  READ(IUO2,23) JDUM,DTHETA(JANGLE),DPHI(JANGLE),
     1                          ECIN(JE),SR_1,SF_1,SR_2,SF_2
                ELSE
                  READ(IUO2,24) JDUM,DTHETA(JANGLE),DPHI(JANGLE),
     1                          ECIN(JE),DTHETAA(JANGLE_A),
     2                          DPHIA(JANGLE_A),SR_1,SF_1,SR_2,SF_2
                ENDIF
              ENDIF
C
              IF(JEL.EQ.1) THEN
                SUMR_1(JPLAN,JE,JANGLE_A)=SUMR_1(JPLAN,JE,JANGLE_A)+
     1                                    SR_1*W(JANGLE)
                SUMF_1(JPLAN,JE,JANGLE_A)=SUMF_1(JPLAN,JE,JANGLE_A)+
     1                                    SF_1*W(JANGLE)
              ELSEIF(JEL.EQ.2) THEN
                SUMR_1(JPLAN,JE,JANGLE)=SUMR_1(JPLAN,JE,JANGLE)+
     1                                  SR_1*W_A(JANGLE_A)
                SUMF_1(JPLAN,JE,JANGLE)=SUMF_1(JPLAN,JE,JANGLE)+
     1                                  SF_1*W_A(JANGLE_A)
              ELSEIF(JEL.EQ.0) THEN
                SUMR_1(JPLAN,JE,1)=SUMR_1(JPLAN,JE,1)+
     1                             SR_1*W(JANGLE)*W_A(JANGLE_A)
                SUMF_1(JPLAN,JE,1)=SUMF_1(JPLAN,JE,1)+
     1                             SF_1*W(JANGLE)*W_A(JANGLE_A)
              ENDIF
              IF(IDICHR.GT.0) THEN
                IF(JEL.EQ.1) THEN
                  SUMR_2(JPLAN,JE,JANGLE_A)=SUMR_2(JPLAN,JE,JANGLE_A)+
     1                                      SR_2*W(JANGLE)
                  SUMF_2(JPLAN,JE,JANGLE_A)=SUMF_2(JPLAN,JE,JANGLE_A)+
     1                                      SF_2*W(JANGLE)
                ELSEIF(JEL.EQ.2) THEN
                  SUMR_2(JPLAN,JE,JANGLE)=SUMR_2(JPLAN,JE,JANGLE)+
     1                                    SR_2*W_A(JANGLE_A)
                  SUMF_2(JPLAN,JE,JANGLE)=SUMF_2(JPLAN,JE,JANGLE)+
     1                                    SF_2*W_A(JANGLE_A)
                ELSEIF(JEL.EQ.0) THEN
                  SUMR_2(JPLAN,JE,1)=SUMR_2(JPLAN,JE,1)+
     1                               SR_2*W(JANGLE)*W_A(JANGLE_A)
                  SUMF_2(JPLAN,JE,1)=SUMF_2(JPLAN,JE,1)+
     1                               SF_2*W(JANGLE)*W_A(JANGLE_A)
                ENDIF
              ENDIF
C
            ENDDO
C
          ENDDO
          IF(I_EXT_A.NE.0) THEN
            REWIND IUI9
            READ(IUI9,1) NDUM
            READ(IUI9,1) NDUM
          ENDIF
        ENDDO
C
        IF(I_EXT.NE.0) THEN
          REWIND IUI6
          READ(IUI6,1) NDUM
          READ(IUI6,1) NDUM
        ENDIF
      ENDDO
C
      CLOSE(IUI6)
      CLOSE(IUI9)
      REWIND IUO2
C
      WRITE(IUO2,16) SPECTRO2,LIKE,SPECTRO,OUTDATA
      IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
        WRITE(IUO2,19) ISPIN,IDICHR,I_SO,ISFLIP
        WRITE(IUO2,18) NE,NPLAN,ISOM
      ELSEIF(JEL.EQ.1) THEN
        WRITE(IUO2,20) ISPIN_A,IDICHR_A,I_SO_A,ISFLIP_A,ICHKDIR_A,
     1                 IPHI_A,ITHETA_A,IE_A
        WRITE(IUO2,21) NPHI0,NTHETA0,NE,NPLAN,ISOM
      ELSEIF(JEL.EQ.2) THEN
        WRITE(IUO2,20) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        WRITE(IUO2,21) NPHI0,NTHETA0,NE,NPLAN,ISOM
      ENDIF
C
      DO JE=1,NE
        DO JANGLE=1,NANGLE
          IF(SPECTRO.EQ.'APC') THEN
            IF(JEL.EQ.1) THEN
              THETA=DTHETAA(JANGLE)
              PHI=DPHIA(JANGLE)
            ELSEIF(JEL.EQ.2) THEN
              THETA=DTHETA(JANGLE)
              PHI=DPHI(JANGLE)
            ENDIF
          ENDIF
C
          DO JPLAN=1,NPLAN
            IF(IDICHR.EQ.0) THEN
              IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
                WRITE(IUO2,33) JPLAN,ECIN(JE),SUMR_1(JPLAN,JE,JANGLE),
     1                         SUMF_1(JPLAN,JE,JANGLE)
              ELSE
                WRITE(IUO2,34) JPLAN,THETA,PHI,ECIN(JE),
     1                         SUMR_1(JPLAN,JE,JANGLE),
     2                         SUMF_1(JPLAN,JE,JANGLE)
              ENDIF
            ELSE
              IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
                WRITE(IUO2,43) JPLAN,ECIN(JE),SUMR_1(JPLAN,JE,JANGLE),
     1                         SUMF_1(JPLAN,JE,JANGLE),
     2                         SUMR_2(JPLAN,JE,JANGLE),
     3                         SUMF_2(JPLAN,JE,JANGLE)
              ELSE
                WRITE(IUO2,44) JPLAN,THETA,PHI,ECIN(JE),
     1                         SUMR_1(JPLAN,JE,JANGLE),
     2                         SUMF_1(JPLAN,JE,JANGLE),
     3                         SUMR_2(JPLAN,JE,JANGLE),
     4                         SUMF_2(JPLAN,JE,JANGLE)
              ENDIF
            ENDIF
          ENDDO
C
          IF(IDICHR.EQ.0) THEN
            IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
              WRITE(IUO2,33) JVOL,ECIN(JE),SUMR_1(NPLAN+1,JE,JANGLE),
     1                       SUMF_1(NPLAN+1,JE,JANGLE)
              WRITE(IUO2,33) JTOT,ECIN(JE),SUMR_1(NPLAN+2,JE,JANGLE),
     1                       SUMF_1(NPLAN+2,JE,JANGLE)
            ELSE
              WRITE(IUO2,34) JVOL,THETA,PHI,ECIN(JE),
     1                       SUMR_1(NPLAN+1,JE,JANGLE),
     2                       SUMF_1(NPLAN+1,JE,JANGLE)
              WRITE(IUO2,34) JTOT,THETA,PHI,ECIN(JE),
     1                       SUMR_1(NPLAN+2,JE,JANGLE),
     2                       SUMF_1(NPLAN+2,JE,JANGLE)
            ENDIF
          ELSE
            IF((SPECTRO.NE.'APC').OR.(JEL.EQ.0)) THEN
              WRITE(IUO2,43) JVOL,ECIN(JE),SUMR_1(NPLAN+1,JE,JANGLE),
     1                       SUMF_1(NPLAN+1,JE,JANGLE),
     2                       SUMR_2(NPLAN+1,JE,JANGLE),
     3                       SUMF_2(NPLAN+1,JE,JANGLE)
              WRITE(IUO2,43) JTOT,ECIN(JE),SUMR_1(NPLAN+2,JE,JANGLE),
     1                       SUMF_1(NPLAN+2,JE,JANGLE),
     2                       SUMR_2(NPLAN+2,JE,JANGLE),
     3                       SUMF_2(NPLAN+2,JE,JANGLE)
            ELSE
              WRITE(IUO2,44) JVOL,THETA,PHI,ECIN(JE),
     1                       SUMR_1(NPLAN+1,JE,JANGLE),
     2                       SUMF_1(NPLAN+1,JE,JANGLE),
     3                       SUMR_2(NPLAN+1,JE,JANGLE),
     4                       SUMF_2(NPLAN+1,JE,JANGLE)
              WRITE(IUO2,44) JTOT,THETA,PHI,ECIN(JE),
     1                       SUMR_1(NPLAN+2,JE,JANGLE),
     2                       SUMF_1(NPLAN+2,JE,JANGLE),
     3                       SUMR_2(NPLAN+2,JE,JANGLE),
     4                       SUMF_2(NPLAN+2,JE,JANGLE)
            ENDIF
          ENDIF
C
        ENDDO
      ENDDO
C
   1  FORMAT(13X,I4)
   2  FORMAT(15X,F8.3,3X,F8.3,3X,E12.6)
   3  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
   4  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL ',
     1       'IN THE WEIGHT_SUM SUBROUTINE - INCREASE NPM TO ',I3,
     2       '>>>>>>>>>>')
   5  FORMAT(6X,I1,1X,I3,3X,I3)
   8  FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
   9  FORMAT(9(2X,I1),2X,I2)
  13  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,
     1       F6.2,2X,E12.6,2X,E12.6)
  15  FORMAT(2X,A3,11X,A13)
  16  FORMAT(2X,A3,A5,1X,A3,2X,A13)
  18  FORMAT(I4,2X,I3,2X,I1)
  19  FORMAT(4(2X,I1))
  20  FORMAT(8(2X,I1))
  21  FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
  23  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,
     1       2X,E12.6,2X,E12.6)
  24  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,
     1       F6.2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
  33  FORMAT(2X,I3,2X,F8.2,2X,E12.6,2X,E12.6)
  34  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
  43  FORMAT(2X,I3,2X,F8.2,2X,E12.6,2X,E12.6,
     1       2X,E12.6,2X,E12.6)
  44  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,
     1       2X,E12.6,2X,E12.6)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE DWSPH(JTYP,JE,X,TLT,ISPEED)
C
C  This routine recomputes the T-matrix elements taking into account the
C                mean square displacements.
C
C  When the argument X is tiny, no vibrations are taken into account
C
C                                      Last modified : 25 Apr 2013
C
      INCLUDE 'spec.inc'
C
      COMMON /TRANS/ DLT(NE_M,NATM,0:18,5),TL(0:NT_M,4,NATM,NE_M),
     1               VK(NE_M),VK2(NE_M),IPOTC,ITL,LMAX(NATM,NE_M)
C
      DIMENSION GNT(0:N_GAUNT)
C
      COMPLEX TLT(0:NT_M,4,NATM,NE_M),TL,SL1,DLT,VK,ZEROC
C
      COMPLEX*16 FFL(0:2*NL_M)
C
      DATA PI4,EPS /12.566371,1.0E-10/
C
      ZEROC=(0.,0.)
C
      IF(X.GT.EPS) THEN
C
C  Standard case: vibrations
C
        IF(ISPEED.LT.0) THEN
          NSUM_LB=ABS(ISPEED)
        ENDIF
C
        COEF=PI4*EXP(-X)
        NL2=2*LMAX(JTYP,JE)+2
        IBESP=5
        MG1=0
        MG2=0
C
        CALL BESPHE(NL2,IBESP,X,FFL)
C
        DO L=0,LMAX(JTYP,JE)
          XL=FLOAT(L+L+1)
          SL1=ZEROC
C
          DO L1=0,LMAX(JTYP,JE)
            XL1=FLOAT(L1+L1+1)
            CALL GAUNT(L,MG1,L1,MG2,GNT)
            L2MIN=ABS(L1-L)
            IF(ISPEED.GE.0) THEN
              L2MAX=L1+L
            ELSEIF(ISPEED.LT.0) THEN
              L2MAX=L2MIN+2*(NSUM_LB-1)
            ENDIF
            SL2=0.
C
            DO L2=L2MIN,L2MAX,2
              XL2=FLOAT(L2+L2+1)
              C=SQRT(XL1*XL2/(PI4*XL))
              SL2=SL2+C*GNT(L2)*REAL(DREAL(FFL(L2)))
            ENDDO
C
            SL1=SL1+SL2*TL(L1,1,JTYP,JE)
          ENDDO
C
          TLT(L,1,JTYP,JE)=COEF*SL1
C
        ENDDO
C
      ELSE
C
C  Argument X tiny: no vibrations
C
        DO L=0,LMAX(JTYP,JE)
C
          TLT(L,1,JTYP,JE)=TL(L,1,JTYP,JE)
C
        ENDDO
C
      ENDIF
C
      RETURN
C
      END

C
C=======================================================================
C
      SUBROUTINE FACDIF1(VKE,RJ,RJK,THRJ,PHIRJ,BETA,GAMMA,L,M,
     1                   FSPH,JAT,JE,*)
C
C  This routine computes a spherical wave scattering factor
C
C                                        Last modified : 03/04/2006
C
      INCLUDE 'spec.inc'
C
      DIMENSION PLMM(0:100,0:100)
      DIMENSION D(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
C
      COMPLEX HLM(0:NO_ST_M,0:NL_M-1),HLN(0:NO_ST_M,0:NL_M-1),FSPH,RHOJ
      COMPLEX HLM1,HLM2,HLM3,HLM4,ALMU,BLMU,SLP,SNU,SMU,TL,VK,DLT,VKE
      COMPLEX RHOJK
C
      COMMON /APPROX/ NDIF,NO,ISPHER,IFWD,NTHOUT,RTHFWD(NATP_M),
     1                IBWD(NATP_M),RTHBWD(NATP_M),IPW,NCUT,PCTINT,IPP,
     2                ISPEED,IATTS,ILENGTH,RLENGTH
      COMMON /EXPFAC/ EXPF(0:2*NL_M-2,0:2*NL_M-2)
      COMMON /TRANS/ DLT(NE_M,NATM,0:18,5),TL(0:NT_M,4,NATM,NE_M),
     1               VK(NE_M),VK2(NE_M),IPOTC,ITL,LMAX(NATM,NE_M)
      COMMON /TYPCAL/ I2,I3,I4,IFTHET,I5,I6,I7,I8,I9
C
      DATA PI/3.141593/
C
      A=1.
      INTER=0
      IF(ITL.EQ.1) VKE=VK(JE)
      RHOJ=VKE*RJ
      RHOJK=VKE*RJK
      HLM1=(1.,0.)
      HLM2=(1.,0.)
      HLM3=(1.,0.)
      HLM4=(1.,0.)
      IEM=1
      CSTH=COS(BETA)
      IF((IFTHET.EQ.0).OR.(THRJ.LT.0.0001)) THEN
        INTER=1
        BLMU=SQRT(4.*PI/FLOAT(2*L+1))*CEXP((0.,-1.)*M*(PHIRJ-PI))
      ENDIF
      CALL PLM(CSTH,PLMM,LMAX(JAT,JE))
      IF(ISPHER.EQ.0) NO1=0
      IF(ISPHER.EQ.1) THEN
        IF(NO.EQ.8) THEN
          NO1=LMAX(JAT,JE)+1
        ELSE
          NO1=NO
        ENDIF
        CALL POLHAN(ISPHER,NO1,LMAX(JAT,JE),RHOJ,HLM)
        IF(IEM.EQ.0) THEN
          HLM4=HLM(0,L)
        ENDIF
        IF(RJK.GT.0.0001) THEN
          NDUM=0
          CALL POLHAN(ISPHER,NDUM,LMAX(JAT,JE),RHOJK,HLN)
        ENDIF
        CALL DJMN(THRJ,D,L)
        A1=ABS(D(0,M,L))
        IF(((A1.LT.0.0001).AND.(IFTHET.EQ.1)).AND.(INTER.EQ.0)) RETURN 1
      ENDIF
      MUMAX=MIN0(L,NO1)
      SMU=(0.,0.)
      DO 10 MU=0,MUMAX
        IF(MOD(MU,2).EQ.0) THEN
          B=1.
        ELSE
          B=-1.
          IF(SIN(BETA).LT.0.) THEN
            A=-1.
          ENDIF
        ENDIF
        IF(ISPHER.LE.1) THEN
          ALMU=(1.,0.)
          C=1.
        ENDIF
        IF(ISPHER.EQ.0) GOTO 40
        IF(INTER.EQ.0) BLMU=CMPLX(D(M,0,L))
        IF(MU.GT.0) THEN
          C=B*FLOAT(L+L+1)/EXPF(MU,L)
          ALMU=(D(M,MU,L)*CEXP((0.,-1.)*MU*GAMMA)+B*
     *          CEXP((0.,1.)*MU*GAMMA)*D(M,-MU,L))/BLMU
        ELSE
          C=1.
          ALMU=CMPLX(D(M,0,L))/BLMU
        ENDIF
  40    SNU=(0.,0.)
        NU1=INT(0.5*(NO1-MU)+0.0001)
        NUMAX=MIN0(NU1,L-MU)
        DO 20 NU=0,NUMAX
          SLP=(0.,0.)
          LPMIN=MAX0(MU,NU)
          DO 30 LP=LPMIN,LMAX(JAT,JE)
            IF(ISPHER.EQ.1) THEN
              HLM1=HLM(NU,LP)
              IF(RJK.GT.0.0001) HLM3=HLN(0,LP)
            ENDIF
            SLP=SLP+FLOAT(2*LP+1)*TL(LP,1,JAT,JE)*HLM1*PLMM(LP,MU)*HLM3
  30      CONTINUE
          IF(ISPHER.EQ.1) THEN
            HLM2=HLM(MU+NU,L)
          ENDIF
          SNU=SNU+SLP*HLM2
  20    CONTINUE
        SMU=SMU+SNU*C*ALMU*A*B
  10  CONTINUE
      FSPH=SMU/(VKE*HLM4)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE FACDIF(COSTH,JAT,JE,FTHETA)
C
C  This routine computes the plane wave scattering factor
C
      INCLUDE 'spec.inc'
C
      COMMON /TRANS/ DLT(NE_M,NATM,0:18,5),TL(0:NT_M,4,NATM,NE_M),
     1               VK(NE_M),VK2(NE_M),IPOTC,ITL,LMAX(NATM,NE_M)
C
      DIMENSION PL(0:100)
C
      COMPLEX FTHETA,DLT,TL,VK
C
      FTHETA=(0.,0.)
      NL=LMAX(JAT,JE)+1
      CALL POLLEG(NL,COSTH,PL)
      DO 20 L=0,NL-1
        FTHETA=FTHETA+(2*L+1)*TL(L,1,JAT,JE)*PL(L)
  20  CONTINUE
      FTHETA=FTHETA/VK(JE)
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE PLOTFD(A,LMX,ITL,NL,NAT,NE)
C
C  This routine prepares the output for a plot of the scattering factor
C
      INCLUDE 'spec.inc'
C
      CHARACTER*24 OUTFILE1,OUTFILE2,OUTFILE3,OUTFILE4
C
      DIMENSION LMX(NATM,NE_M)
C
      COMPLEX FSPH,VKE
C
      CHARACTER*3 S_O
C
      COMMON /APPROX/ NDIF,NO,ISPHER,IFWD,NTHOUT,RTHFWD(NATP_M),
     1                IBWD(NATP_M),RTHBWD(NATP_M),IPW,NCUT,PCTINT,IPP,
     2                ISPEED,IATTS,ILENGTH,RLENGTH
      COMMON /FDIF/ R1,R2
      COMMON /INIT_L/ L,I2,I3,I4,I5,I10
      COMMON /INIT_J/ JF1,JF2,I_SO,S_O
      COMMON /OUTFILES/ OUTFILE1,OUTFILE2,OUTFILE3,OUTFILE4
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /PARCAL/ N3,N4,N5,NFTHET,N6
      COMMON /TYPCAL/ IPHI,IE,ITHETA,I7,I8,I9,I12,I13,I14
      COMMON /VALIN/ PHI0,E0,THETA0,U1,U2,U3,VINT,N7(100)
      COMMON /VALFIN/ PHI1,EFIN,THETA1
C
      DATA PI,CONV/3.141593,0.512314/
C
      OPEN(UNIT=IUO3, FILE=OUTFILE3, STATUS='UNKNOWN')
      IF(ISPHER.EQ.0) THEN
        L=0
        LMAX=0
      ELSE
        LMAX=L
      ENDIF
      PHITOT=360.
      THTOT=360.*ITHETA*(1-IPHI)+180.*ITHETA*IPHI
      NPHI=(NFTHET+1)*IPHI+(1-IPHI)
      NTHT=(NFTHET+1)*ITHETA*(1-IPHI)+(NFTHET/2+1)*ITHETA*IPHI+
     *     (1-ITHETA)
      NE=NFTHET*IE + (1-IE)
      WRITE(IUO3,1) ISPHER,NL,NAT,L,NTHT,NPHI,NE,E0,EFIN
      DO 10 JT=1,NTHT
        DTHETA=THETA1+FLOAT(JT-1)*THTOT/FLOAT(MAX0(NTHT-1,1))
        RTHETA=DTHETA*PI/180.
        TEST=SIN(RTHETA)
        IF(TEST.GE.0.) THEN
          POZ=PI
          EPS=1.
        ELSE
          POZ=0.
          EPS=-1.
        ENDIF
        BETA=RTHETA*EPS
        IF(ABS(TEST).LT.0.0001) THEN
          NPHIM=1
        ELSE
          NPHIM=NPHI
        ENDIF
        DO 20 JP=1,NPHIM
          DPHI=PHI1+FLOAT(JP-1)*PHITOT/FLOAT(MAX0(NPHI-1,1))
          RPHI=DPHI*PI/180.
          GAMMA=POZ-RPHI
          DO 30 JE=1,NE
            IF(NE.EQ.1) THEN
              ECIN=E0
            ELSE
              ECIN=E0+FLOAT(JE-1)*(EFIN-E0)/FLOAT(NE-1)
            ENDIF
            IF(ITL.EQ.0) VKE=SQRT(ECIN-ABS(VINT))*CONV*A*(1.,0.)
            DO 40 JAT=1,NAT
              IF(L.GT.LMX(JAT,JE)) GOTO 90
              DO 50 M=-LMAX,LMAX
                CALL FACDIF1(VKE,R1,R2,THETA0,PHI0,BETA,GAMMA,L,M,
     1                       FSPH,JAT,JE,*60)
                GOTO 70
  60            WRITE(IUO1,80)
                STOP
  70            REFTH=REAL(FSPH)
                XIMFTH=AIMAG(FSPH)
                WRITE(IUO3,5) JE,JAT,L,M,REFTH,XIMFTH,DTHETA,DPHI,ECIN
  50          CONTINUE
              GOTO 40
  90          WRITE(IUO1,100) JAT
              STOP
  40        CONTINUE
  30      CONTINUE
  20    CONTINUE
  10  CONTINUE
      CLOSE(IUO3)
   1  FORMAT(5X,I1,2X,I2,2X,I4,2X,I2,2X,I3,2X,I3,2X,I3,2X,F8.2,2X,F8.2)
   5  FORMAT(1X,I3,1X,I4,1X,I2,1X,I3,1X,F6.3,1X,F6.3,1X,F6.2,
     1       1X,F6.2,1X,F8.2)
  80  FORMAT(15X,'<<<<<  WRONG VALUE OF THETA0 : THE DENOMINATOR ',
     1           'IS ZERO  >>>>>')
 100  FORMAT(15X,'<<<<<  THE VALUE OF L EST IS TOO LARGE FOR ATOM',
     1           ' : ',I2,'  >>>>>')
C
      RETURN
C
      END
C
C=======================================================================
C
      SUBROUTINE INV_MAT_MS(JE,TAU)
C
C   This subroutine stores the multiple scattering matrix and computes
C     the scattering path operator TAU^{j 0} exactly, without explicitely
C     using the inverse matrix.
C
C                           (Photoelectron case)
C
C                                         Last modified : 28 Mar 2007
C
      INCLUDE 'spec.inc'
C
      PARAMETER(NLTWO=2*NL_M)
C
      COMPLEX*16 HL1(0:NLTWO),SM(LINMAX*NATCLU_M,LINMAX*NATCLU_M)
      COMPLEX*16 IN(LINMAX*NATCLU_M,LINMAX)
      COMPLEX*16 SUM_L,ONEC,IC,ZEROC
      COMPLEX*16 YLM(0:NLTWO,-NLTWO:NLTWO),TLJ,TLK,EXPKJ
C
      COMPLEX TAU(LINMAX,LINFMAX,NATCLU_M)
      COMPLEX VK,DLT,TL
C
      REAL*8 PI,ATTKJ,GNT(0:N_GAUNT),XKJ,YKJ,ZKJ,RKJ,ZDKJ,KRKJ
C
      INTEGER IPIV(LINMAX*NATCLU_M)
C
      CHARACTER*1 CH
C
      COMMON /COOR/ NATCLU,N_PROT,NATYP(NATM),NCHTYP(NATP_M),
     1              NCORR(NAT_EQ_M,NATP_M),INEW_AT(NATCLU_M),
     2              SYM_AT(3,NATCLU_M)
      COMMON /INIT_L/ LI,INITL,NNL,LF1,LF2,ISTEP_LF
      COMMON /TRANS/ DLT(NE_M,NATM,0:18,5),TL(0:NT_M,4,NATM,NE_M),
     1               VK(NE_M),VK2(NE_M),IPOTC,ITL,LMAX(NATM,NE_M)
C
      DATA PI /3.1415926535898D0/
C
      ONEC=(1.D0,0.D0)
      IC=(0.D0,1.D0)
      ZEROC=(0.D0,0.D0)
      IBESS=3
      CH='N'
C
C  Construction of the multiple scattering matrix MS = (I-GoT).
C    Elements are stored using a linear index LINJ representing
C    (J,LJ)
C
      JLIN=0
      DO JTYP=1,N_PROT
        NBTYPJ=NATYP(JTYP)
        LMJ=LMAX(JTYP,JE)
        DO JNUM=1,NBTYPJ
          JATL=NCORR(JNUM,JTYP)
          XJ=SYM_AT(1,JATL)
          YJ=SYM_AT(2,JATL)
          ZJ=SYM_AT(3,JATL)
C
          DO LJ=0,LMJ
            ILJ=LJ*LJ+LJ+1
            TLJ=DCMPLX(TL(LJ,1,JTYP,JE))
            DO MJ=-LJ,LJ
              INDJ=ILJ+MJ
              JLIN=JLIN+1
C
              KLIN=0
              DO KTYP=1,N_PROT
                NBTYPK=NATYP(KTYP)
                LMK=LMAX(KTYP,JE)
                DO KNUM=1,NBTYPK
                  KATL=NCORR(KNUM,KTYP)
                  IF(KATL.NE.JATL) THEN
                    XKJ=DBLE(SYM_AT(1,KATL)-XJ)
                    YKJ=DBLE(SYM_AT(2,KATL)-YJ)
                    ZKJ=DBLE(SYM_AT(3,KATL)-ZJ)
                    RKJ=DSQRT(XKJ*XKJ+YKJ*YKJ+ZKJ*ZKJ)
                    KRKJ=DBLE(VK(JE))*RKJ
                    ATTKJ=DEXP(-DIMAG(DCMPLX(VK(JE)))*RKJ)
                    EXPKJ=(XKJ+IC*YKJ)/RKJ
                    ZDKJ=ZKJ/RKJ
                    CALL SPH_HAR2(2*NL_M,ZDKJ,EXPKJ,YLM,LMJ+LMK)
                    CALL BESPHE2(LMJ+LMK+1,IBESS,KRKJ,HL1)
                  ENDIF
C
                  DO LK=0,LMK
                    ILK=LK*LK+LK+1
                    L_MIN=ABS(LK-LJ)
                    L_MAX=LK+LJ
                    TLK=DCMPLX(TL(LK,1,KTYP,JE))
                    DO MK=-LK,LK
                      INDK=ILK+MK
                      KLIN=KLIN+1
                      SM(KLIN,JLIN)=ZEROC
                      SUM_L=ZEROC
                      IF(KATL.NE.JATL) THEN
                        CALL GAUNT2(LK,MK,LJ,MJ,GNT)
C
                        DO L=L_MIN,L_MAX,2
                          M=MJ-MK
                          IF(ABS(M).LE.L) THEN
                            SUM_L=SUM_L+(IC**L)*HL1(L)*
     1                                   YLM(L,M)*GNT(L)
                          ENDIF
                        ENDDO
                        SUM_L=SUM_L*ATTKJ*4.D0*PI*IC
                      ELSE
                        SUM_L=ZEROC
                      ENDIF
C
                      IF(KLIN.EQ.JLIN) THEN
                        SM(KLIN,JLIN)=ONEC-TLK*SUM_L
                        IF(JTYP.EQ.1) THEN
                          IN(KLIN,JLIN)=ONEC
                        ENDIF
                      ELSE
                        SM(KLIN,JLIN)=-TLK*SUM_L
                        IF(JTYP.EQ.1) THEN
                          IN(KLIN,JLIN)=ZEROC
                        ENDIF
                      ENDIF
C
                    ENDDO
                  ENDDO
C
                ENDDO
              ENDDO
C
            ENDDO
          ENDDO
C
        ENDDO
      ENDDO
C
      LW2=(LMAX(1,JE)+1)*(LMAX(1,JE)+1)
C
C   Partial inversion of the multiple scattering matrix MS and
C     multiplication by T : the LAPACK subroutine performing
C
C                          A * x = b
C
C     is used where b is the block column corresponding to
C     the absorber 0 in the identity matrix. x is then TAU^{j 0}.
C
      CALL ZGETRF(JLIN,JLIN,SM,LINMAX*NATCLU_M,IPIV,INFO1)
      IF(INFO1.NE.0) THEN
        WRITE(6,*) '     --->  INFO1 =',INFO1
      ELSE
        CALL ZGETRS(CH,JLIN,LW2,SM,LINMAX*NATCLU_M,IPIV,
     1              IN,LINMAX*NATCLU_M,INFO)
      ENDIF
C
C   Storage of the Tau matrix
C
      JLIN=0
      DO JTYP=1,N_PROT
        NBTYPJ=NATYP(JTYP)
        LMJ=LMAX(JTYP,JE)
        DO JNUM=1,NBTYPJ
          JATL=NCORR(JNUM,JTYP)
C
          DO LJ=0,LMJ
            ILJ=LJ*LJ+LJ+1
            TLJ=DCMPLX(TL(LJ,1,JTYP,JE))
            DO MJ=-LJ,LJ
              INDJ=ILJ+MJ
              JLIN=JLIN+1
C
              KLIN=0
              DO KTYP=1,N_PROT
                NBTYPK=NATYP(KTYP)
                LMK=LMAX(KTYP,JE)
                DO KNUM=1,NBTYPK
                  KATL=NCORR(KNUM,KTYP)
C
                  DO LK=0,LMK
                    ILK=LK*LK+LK+1
                    DO MK=-LK,LK
                      INDK=ILK+MK
                      KLIN=KLIN+1
                      IF((JATL.EQ.1).AND.(LJ.LE.LF2)) THEN
                        TAU(INDK,INDJ,KATL)=CMPLX(IN(KLIN,JLIN)*TLJ)
                      ENDIF
                    ENDDO
                  ENDDO
C
                ENDDO
              ENDDO
C
            ENDDO
          ENDDO
C
        ENDDO
      ENDDO
C
      RETURN
C
      END

C
C=======================================================================
C
C                         LAPACK Ax=b subroutines
C
C=======================================================================
C
C
C======================================================================
C
      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRS solves a system of linear equations
*     A * X = B,  A**T * X = B,  or  A**H * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by ZGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by ZGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASWP, ZTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
        INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
        INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
        INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZGETRS', -INFO )
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
        CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
        CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
        CALL ZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A**T * X = B  or A**H * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
        CALL ZTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
        CALL ZTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,
     $               LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
        CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of ZGETRS
*
      END
C
C
C======================================================================
C
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
        IEEECK = 0
        RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
        IEEECK = 0
        RETURN
      END IF
*
      RETURN
      END
C
C
C======================================================================
C
      INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     January 2007
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  ILAENV returns an INTEGER
*  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
*  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines (DEPRECATED)
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR method
*               for nonsymmetric eigenvalue problems (DEPRECATED)
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*          12 <= ISPEC <= 16:
*               xHSEQR or one of its subroutines,
*               see IPARMQ for detailed explanation
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ
      EXTERNAL           IEEECK, IPARMQ
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1: 1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
        IF( IC.GE.97 .AND. IC.LE.122 ) THEN
          SUBNAM( 1: 1 ) = CHAR( IC-32 )
          DO 20 I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   20     CONTINUE
        END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
        IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
          SUBNAM( 1: 1 ) = CHAR( IC+64 )
          DO 30 I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )SUBNAM( I:
     $             I ) = CHAR( IC+64 )
   30     CONTINUE
        END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
        IF( IC.GE.225 .AND. IC.LE.250 ) THEN
          SUBNAM( 1: 1 ) = CHAR( IC-32 )
          DO 40 I = 2, 6
            IC = ICHAR( SUBNAM( I: I ) )
            IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I: I ) = CHAR( IC-32 )
   40     CONTINUE
        END IF
      END IF
*
      C1 = SUBNAM( 1: 1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      C4 = C3( 2: 3 )
*
      GO TO ( 50, 60, 70 )ISPEC
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
          IF( SNAME ) THEN
            NB = 32
          ELSE
            NB = 32
          END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
          IF( SNAME ) THEN
            NB = 32
          ELSE
            NB = 32
          END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
          IF( SNAME ) THEN
            NB = 32
          ELSE
            NB = 32
          END IF
        ELSE IF( C3.EQ.'TRI' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        END IF
      ELSE IF( C2.EQ.'PO' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        END IF
      ELSE IF( C2.EQ.'SY' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
          NB = 32
        ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
          NB = 64
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          NB = 64
        ELSE IF( C3.EQ.'TRD' ) THEN
          NB = 32
        ELSE IF( C3.EQ.'GST' ) THEN
          NB = 64
        END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NB = 32
          END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NB = 32
          END IF
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NB = 32
          END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NB = 32
          END IF
        END IF
      ELSE IF( C2.EQ.'GB' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            IF( N4.LE.64 ) THEN
              NB = 1
            ELSE
              NB = 32
            END IF
          ELSE
            IF( N4.LE.64 ) THEN
              NB = 1
            ELSE
              NB = 32
            END IF
          END IF
        END IF
      ELSE IF( C2.EQ.'PB' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            IF( N2.LE.64 ) THEN
              NB = 1
            ELSE
              NB = 32
            END IF
          ELSE
            IF( N2.LE.64 ) THEN
              NB = 1
            ELSE
              NB = 32
            END IF
          END IF
        END IF
      ELSE IF( C2.EQ.'TR' ) THEN
        IF( C3.EQ.'TRI' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        END IF
      ELSE IF( C2.EQ.'LA' ) THEN
        IF( C3.EQ.'UUM' ) THEN
          IF( SNAME ) THEN
            NB = 64
          ELSE
            NB = 64
          END IF
        END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
        IF( C3.EQ.'EBZ' ) THEN
          NB = 1
        END IF
      END IF
      ILAENV = NB
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
          IF( SNAME ) THEN
            NBMIN = 2
          ELSE
            NBMIN = 2
          END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
          IF( SNAME ) THEN
            NBMIN = 2
          ELSE
            NBMIN = 2
          END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
          IF( SNAME ) THEN
            NBMIN = 2
          ELSE
            NBMIN = 2
          END IF
        ELSE IF( C3.EQ.'TRI' ) THEN
          IF( SNAME ) THEN
            NBMIN = 2
          ELSE
            NBMIN = 2
          END IF
        END IF
      ELSE IF( C2.EQ.'SY' ) THEN
        IF( C3.EQ.'TRF' ) THEN
          IF( SNAME ) THEN
            NBMIN = 8
          ELSE
            NBMIN = 8
          END IF
        ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
          NBMIN = 2
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRD' ) THEN
          NBMIN = 2
        END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NBMIN = 2
          END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NBMIN = 2
          END IF
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NBMIN = 2
          END IF
        ELSE IF( C3( 1: 1 ).EQ.'M' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NBMIN = 2
          END IF
        END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
        IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. C3.EQ.
     $       'QLF' ) THEN
          IF( SNAME ) THEN
            NX = 128
          ELSE
            NX = 128
          END IF
        ELSE IF( C3.EQ.'HRD' ) THEN
          IF( SNAME ) THEN
            NX = 128
          ELSE
            NX = 128
          END IF
        ELSE IF( C3.EQ.'BRD' ) THEN
          IF( SNAME ) THEN
            NX = 128
          ELSE
            NX = 128
          END IF
        END IF
      ELSE IF( C2.EQ.'SY' ) THEN
        IF( SNAME .AND. C3.EQ.'TRD' ) THEN
          NX = 32
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
        IF( C3.EQ.'TRD' ) THEN
          NX = 32
        END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NX = 128
          END IF
        END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
        IF( C3( 1: 1 ).EQ.'G' ) THEN
          IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. C4.EQ.
     $          'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. C4.EQ.'BR' )
     $           THEN
            NX = 128
          END IF
        END IF
      END IF
      ILAENV = NX
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
        ILAENV = IEEECK( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
        ILAENV = IEEECK( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines.
*
      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
      RETURN
*
*     End of ILAENV
*
      END
C
C
C======================================================================
C
      LOGICAL FUNCTION LSAME(CA,CB)
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER CA,CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
*     ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
*
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
        IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
        IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
*
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
        IF (INTA.GE.129 .AND. INTA.LE.137 .OR.
     +        INTA.GE.145 .AND. INTA.LE.153 .OR.
     +        INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
        IF (INTB.GE.129 .AND. INTB.LE.137 .OR.
     +        INTB.GE.145 .AND. INTB.LE.153 .OR.
     +        INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
*
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
        IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
        IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
      LSAME = INTA .EQ. INTB
*
*     RETURN
*
*     End of LSAME
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IZAMAX
      EXTERNAL           DLAMCH, IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZGETF2', -INFO )
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum
*
      SFMIN = DLAMCH('S')
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
        JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
        IPIV( J ) = JP
        IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
          IF( JP.NE.J )
     $         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
          IF( J.LT.M ) THEN
            IF( ABS(A( J, J )) .GE. SFMIN ) THEN
              CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
            ELSE
              DO 20 I = 1, M-J
                A( J+I, J ) = A( J+I, J ) / A( J, J )
   20         CONTINUE
            END IF
          END IF
*
        ELSE IF( INFO.EQ.0 ) THEN
*
          INFO = J
        END IF
*
        IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
          CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),
     $                  LDA, A( J+1, J+1 ), LDA )
        END IF
   10 CONTINUE
      RETURN
*
*     End of ZGETF2
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETF2, ZLASWP, ZTRSM
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'ZGETRF', -INFO )
        RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
        CALL ZGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
        DO 20 J = 1, MIN( M, N ), NB
          JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
          CALL ZGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
          IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
          DO 10 I = J, MIN( M, J+JB-1 )
            IPIV( I ) = J - 1 + IPIV( I )
   10     CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
          CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
          IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
            CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
            CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
            IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
              CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
            END IF
          END IF
   20   CONTINUE
      END IF
      RETURN
*
*     End of ZGETRF
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
        IX0 = K1
        I1 = K1
        I2 = K2
        INC = 1
      ELSE IF( INCX.LT.0 ) THEN
        IX0 = 1 + ( 1-K2 )*INCX
        I1 = K2
        I2 = K1
        INC = -1
      ELSE
        RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
        DO 30 J = 1, N32, 32
          IX = IX0
          DO 20 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
              DO 10 K = J, J + 31
                TEMP = A( I, K )
                A( I, K ) = A( IP, K )
                A( IP, K ) = TEMP
   10         CONTINUE
            END IF
            IX = IX + INCX
   20     CONTINUE
   30   CONTINUE
      END IF
      IF( N32.NE.N ) THEN
        N32 = N32 + 1
        IX = IX0
        DO 50 I = I1, I2, INC
          IP = IPIV( IX )
          IF( IP.NE.I ) THEN
            DO 40 K = N32, N
              TEMP = A( I, K )
              A( I, K ) = A( IP, K )
              A( IP, K ) = TEMP
   40       CONTINUE
          END IF
          IX = IX + INCX
   50   CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASWP
*
      END
C
C
C======================================================================
C
      SUBROUTINE XERBLA(SRNAME,INFO)
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER INFO
      CHARACTER*6 SRNAME
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE (*,FMT=9999) SRNAME,INFO
*
      STOP
*
 9999 FORMAT (' ** On entry to ',A6,' parameter number ',I2,' had ',
     +       'an illegal value')
*
*     End of XERBLA
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*
*  Purpose
*  =======
*
*  ZGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL CONJA,CONJB,NOTA,NOTB
*     ..
*     .. Parameters ..
      DOUBLE COMPLEX ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
*     and the number of rows of  B  respectively.
*
      NOTA = LSAME(TRANSA,'N')
      NOTB = LSAME(TRANSB,'N')
      CONJA = LSAME(TRANSA,'C')
      CONJB = LSAME(TRANSB,'C')
      IF (NOTA) THEN
        NROWA = M
        NCOLA = K
      ELSE
        NROWA = K
        NCOLA = M
      END IF
      IF (NOTB) THEN
        NROWB = K
      ELSE
        NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF ((.NOT.NOTA) .AND. (.NOT.CONJA) .AND.
     +    (.NOT.LSAME(TRANSA,'T'))) THEN
        INFO = 1
      ELSE IF ((.NOT.NOTB) .AND. (.NOT.CONJB) .AND.
     +         (.NOT.LSAME(TRANSB,'T'))) THEN
        INFO = 2
      ELSE IF (M.LT.0) THEN
        INFO = 3
      ELSE IF (N.LT.0) THEN
        INFO = 4
      ELSE IF (K.LT.0) THEN
        INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
        INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
        INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
        INFO = 13
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('ZGEMM ',INFO)
        RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
        IF (BETA.EQ.ZERO) THEN
          DO 20 J = 1,N
            DO 10 I = 1,M
              C(I,J) = ZERO
   10       CONTINUE
   20     CONTINUE
        ELSE
          DO 40 J = 1,N
            DO 30 I = 1,M
              C(I,J) = BETA*C(I,J)
   30       CONTINUE
   40     CONTINUE
        END IF
        RETURN
      END IF
*
*     Start the operations.
*
      IF (NOTB) THEN
        IF (NOTA) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
          DO 90 J = 1,N
            IF (BETA.EQ.ZERO) THEN
              DO 50 I = 1,M
                C(I,J) = ZERO
   50         CONTINUE
            ELSE IF (BETA.NE.ONE) THEN
              DO 60 I = 1,M
                C(I,J) = BETA*C(I,J)
   60         CONTINUE
            END IF
            DO 80 L = 1,K
              IF (B(L,J).NE.ZERO) THEN
                TEMP = ALPHA*B(L,J)
                DO 70 I = 1,M
                  C(I,J) = C(I,J) + TEMP*A(I,L)
   70           CONTINUE
              END IF
   80       CONTINUE
   90     CONTINUE
        ELSE IF (CONJA) THEN
*
*           Form  C := alpha*conjg( A' )*B + beta*C.
*
          DO 120 J = 1,N
            DO 110 I = 1,M
              TEMP = ZERO
              DO 100 L = 1,K
                TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
  100         CONTINUE
              IF (BETA.EQ.ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
  110       CONTINUE
  120     CONTINUE
        ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
          DO 150 J = 1,N
            DO 140 I = 1,M
              TEMP = ZERO
              DO 130 L = 1,K
                TEMP = TEMP + A(L,I)*B(L,J)
  130         CONTINUE
              IF (BETA.EQ.ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
  140       CONTINUE
  150     CONTINUE
        END IF
      ELSE IF (NOTA) THEN
        IF (CONJB) THEN
*
*           Form  C := alpha*A*conjg( B' ) + beta*C.
*
          DO 200 J = 1,N
            IF (BETA.EQ.ZERO) THEN
              DO 160 I = 1,M
                C(I,J) = ZERO
  160         CONTINUE
            ELSE IF (BETA.NE.ONE) THEN
              DO 170 I = 1,M
                C(I,J) = BETA*C(I,J)
  170         CONTINUE
            END IF
            DO 190 L = 1,K
              IF (B(J,L).NE.ZERO) THEN
                TEMP = ALPHA*DCONJG(B(J,L))
                DO 180 I = 1,M
                  C(I,J) = C(I,J) + TEMP*A(I,L)
  180           CONTINUE
              END IF
  190       CONTINUE
  200     CONTINUE
        ELSE
*
*           Form  C := alpha*A*B'          + beta*C
*
          DO 250 J = 1,N
            IF (BETA.EQ.ZERO) THEN
              DO 210 I = 1,M
                C(I,J) = ZERO
  210         CONTINUE
            ELSE IF (BETA.NE.ONE) THEN
              DO 220 I = 1,M
                C(I,J) = BETA*C(I,J)
  220         CONTINUE
            END IF
            DO 240 L = 1,K
              IF (B(J,L).NE.ZERO) THEN
                TEMP = ALPHA*B(J,L)
                DO 230 I = 1,M
                  C(I,J) = C(I,J) + TEMP*A(I,L)
  230           CONTINUE
              END IF
  240       CONTINUE
  250     CONTINUE
        END IF
      ELSE IF (CONJA) THEN
        IF (CONJB) THEN
*
*           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
*
          DO 280 J = 1,N
            DO 270 I = 1,M
              TEMP = ZERO
              DO 260 L = 1,K
                TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
  260         CONTINUE
              IF (BETA.EQ.ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
  270       CONTINUE
  280     CONTINUE
        ELSE
*
*           Form  C := alpha*conjg( A' )*B' + beta*C
*
          DO 310 J = 1,N
            DO 300 I = 1,M
              TEMP = ZERO
              DO 290 L = 1,K
                TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
  290         CONTINUE
              IF (BETA.EQ.ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
  300       CONTINUE
  310     CONTINUE
        END IF
      ELSE
        IF (CONJB) THEN
*
*           Form  C := alpha*A'*conjg( B' ) + beta*C
*
          DO 340 J = 1,N
            DO 330 I = 1,M
              TEMP = ZERO
              DO 320 L = 1,K
                TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
  320         CONTINUE
              IF (BETA.EQ.ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
  330       CONTINUE
  340     CONTINUE
        ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
          DO 370 J = 1,N
            DO 360 I = 1,M
              TEMP = ZERO
              DO 350 L = 1,K
                TEMP = TEMP + A(L,I)*B(J,L)
  350         CONTINUE
              IF (BETA.EQ.ZERO) THEN
                C(I,J) = ALPHA*TEMP
              ELSE
                C(I,J) = ALPHA*TEMP + BETA*C(I,J)
              END IF
  360       CONTINUE
  370     CONTINUE
        END IF
      END IF
*
      RETURN
*
*     End of ZGEMM .
*
      END
C
C
C======================================================================
C
      SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*  ZGERU  performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Arguments
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (M.LT.0) THEN
        INFO = 1
      ELSE IF (N.LT.0) THEN
        INFO = 2
      ELSE IF (INCX.EQ.0) THEN
        INFO = 5
      ELSE IF (INCY.EQ.0) THEN
        INFO = 7
      ELSE IF (LDA.LT.MAX(1,M)) THEN
        INFO = 9
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('ZGERU ',INFO)
        RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (INCY.GT.0) THEN
        JY = 1
      ELSE
        JY = 1 - (N-1)*INCY
      END IF
      IF (INCX.EQ.1) THEN
        DO 20 J = 1,N
          IF (Y(JY).NE.ZERO) THEN
            TEMP = ALPHA*Y(JY)
            DO 10 I = 1,M
              A(I,J) = A(I,J) + X(I)*TEMP
   10       CONTINUE
          END IF
          JY = JY + INCY
   20   CONTINUE
      ELSE
        IF (INCX.GT.0) THEN
          KX = 1
        ELSE
          KX = 1 - (M-1)*INCX
        END IF
        DO 40 J = 1,N
          IF (Y(JY).NE.ZERO) THEN
            TEMP = ALPHA*Y(JY)
            IX = KX
            DO 30 I = 1,M
              A(I,J) = A(I,J) + X(IX)*TEMP
              IX = IX + INCX
   30       CONTINUE
          END IF
          JY = JY + INCY
   40   CONTINUE
      END IF
*
      RETURN
*
*     End of ZGERU .
*
      END
C
C
      DOUBLE PRECISION FUNCTION DCABS1(Z)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
*     ..
*     ..
*  Purpose
*  =======
*
*  DCABS1 computes absolute value of a double complex number
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG
*
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
      END
C
C
C======================================================================
C
      SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ZA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*)
*     ..
*
*  Purpose
*  =======
*
*     scales a vector by a constant.
*     jack dongarra, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      DO 10 I = 1,N
        ZX(IX) = ZA*ZX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 DO 30 I = 1,N
        ZX(I) = ZA*ZX(I)
   30 CONTINUE
      RETURN
      END
C
C
C======================================================================
C
      SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*),ZY(*)
*     ..
*
*  Purpose
*  =======
*
*     interchanges two vectors.
*     jack dongarra, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE COMPLEX ZTEMP
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*       code for unequal increments or equal increments not equal
*         to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        ZTEMP = ZX(IX)
        ZX(IX) = ZY(IY)
        ZY(IY) = ZTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
*
*       code for both increments equal to 1
   20 DO 30 I = 1,N
        ZTEMP = ZX(I)
        ZX(I) = ZY(I)
        ZY(I) = ZTEMP
   30 CONTINUE
      RETURN
      END
C
C
C======================================================================
C
      SUBROUTINE ZTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*     .. Scalar Arguments ..
      DOUBLE COMPLEX ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,*),B(LDB,*)
*     ..
*
*  Purpose
*  =======
*
*  ZTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
*
*  The matrix X is overwritten on B.
*
*  Arguments
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DCONJG,MAX
*     ..
*     .. Local Scalars ..
      DOUBLE COMPLEX TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      DOUBLE COMPLEX ONE
      PARAMETER (ONE= (1.0D+0,0.0D+0))
      DOUBLE COMPLEX ZERO
      PARAMETER (ZERO= (0.0D+0,0.0D+0))
*     ..
*
*     Test the input parameters.
*
      LSIDE = LSAME(SIDE,'L')
      IF (LSIDE) THEN
        NROWA = M
      ELSE
        NROWA = N
      END IF
      NOCONJ = LSAME(TRANSA,'T')
      NOUNIT = LSAME(DIAG,'N')
      UPPER = LSAME(UPLO,'U')
*
      INFO = 0
      IF ((.NOT.LSIDE) .AND. (.NOT.LSAME(SIDE,'R'))) THEN
        INFO = 1
      ELSE IF ((.NOT.UPPER) .AND. (.NOT.LSAME(UPLO,'L'))) THEN
        INFO = 2
      ELSE IF ((.NOT.LSAME(TRANSA,'N')) .AND.
     +         (.NOT.LSAME(TRANSA,'T')) .AND.
     +         (.NOT.LSAME(TRANSA,'C'))) THEN
        INFO = 3
      ELSE IF ((.NOT.LSAME(DIAG,'U')) .AND. (.NOT.LSAME(DIAG,'N'))) THEN
        INFO = 4
      ELSE IF (M.LT.0) THEN
        INFO = 5
      ELSE IF (N.LT.0) THEN
        INFO = 6
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
        INFO = 9
      ELSE IF (LDB.LT.MAX(1,M)) THEN
        INFO = 11
      END IF
      IF (INFO.NE.0) THEN
        CALL XERBLA('ZTRSM ',INFO)
        RETURN
      END IF
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (ALPHA.EQ.ZERO) THEN
        DO 20 J = 1,N
          DO 10 I = 1,M
            B(I,J) = ZERO
   10     CONTINUE
   20   CONTINUE
        RETURN
      END IF
*
*     Start the operations.
*
      IF (LSIDE) THEN
        IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*inv( A )*B.
*
          IF (UPPER) THEN
            DO 60 J = 1,N
              IF (ALPHA.NE.ONE) THEN
                DO 30 I = 1,M
                  B(I,J) = ALPHA*B(I,J)
   30           CONTINUE
              END IF
              DO 50 K = M,1,-1
                IF (B(K,J).NE.ZERO) THEN
                  IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                  DO 40 I = 1,K - 1
                    B(I,J) = B(I,J) - B(K,J)*A(I,K)
   40             CONTINUE
                END IF
   50         CONTINUE
   60       CONTINUE
          ELSE
            DO 100 J = 1,N
              IF (ALPHA.NE.ONE) THEN
                DO 70 I = 1,M
                  B(I,J) = ALPHA*B(I,J)
   70           CONTINUE
              END IF
              DO 90 K = 1,M
                IF (B(K,J).NE.ZERO) THEN
                  IF (NOUNIT) B(K,J) = B(K,J)/A(K,K)
                  DO 80 I = K + 1,M
                    B(I,J) = B(I,J) - B(K,J)*A(I,K)
   80             CONTINUE
                END IF
   90         CONTINUE
  100       CONTINUE
          END IF
        ELSE
*
*           Form  B := alpha*inv( A' )*B
*           or    B := alpha*inv( conjg( A' ) )*B.
*
          IF (UPPER) THEN
            DO 140 J = 1,N
              DO 130 I = 1,M
                TEMP = ALPHA*B(I,J)
                IF (NOCONJ) THEN
                  DO 110 K = 1,I - 1
                    TEMP = TEMP - A(K,I)*B(K,J)
  110             CONTINUE
                  IF (NOUNIT) TEMP = TEMP/A(I,I)
                ELSE
                  DO 120 K = 1,I - 1
                    TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  120             CONTINUE
                  IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                END IF
                B(I,J) = TEMP
  130         CONTINUE
  140       CONTINUE
          ELSE
            DO 180 J = 1,N
              DO 170 I = M,1,-1
                TEMP = ALPHA*B(I,J)
                IF (NOCONJ) THEN
                  DO 150 K = I + 1,M
                    TEMP = TEMP - A(K,I)*B(K,J)
  150             CONTINUE
                  IF (NOUNIT) TEMP = TEMP/A(I,I)
                ELSE
                  DO 160 K = I + 1,M
                    TEMP = TEMP - DCONJG(A(K,I))*B(K,J)
  160             CONTINUE
                  IF (NOUNIT) TEMP = TEMP/DCONJG(A(I,I))
                END IF
                B(I,J) = TEMP
  170         CONTINUE
  180       CONTINUE
          END IF
        END IF
      ELSE
        IF (LSAME(TRANSA,'N')) THEN
*
*           Form  B := alpha*B*inv( A ).
*
          IF (UPPER) THEN
            DO 230 J = 1,N
              IF (ALPHA.NE.ONE) THEN
                DO 190 I = 1,M
                  B(I,J) = ALPHA*B(I,J)
  190           CONTINUE
              END IF
              DO 210 K = 1,J - 1
                IF (A(K,J).NE.ZERO) THEN
                  DO 200 I = 1,M
                    B(I,J) = B(I,J) - A(K,J)*B(I,K)
  200             CONTINUE
                END IF
  210         CONTINUE
              IF (NOUNIT) THEN
                TEMP = ONE/A(J,J)
                DO 220 I = 1,M
                  B(I,J) = TEMP*B(I,J)
  220           CONTINUE
              END IF
  230       CONTINUE
          ELSE
            DO 280 J = N,1,-1
              IF (ALPHA.NE.ONE) THEN
                DO 240 I = 1,M
                  B(I,J) = ALPHA*B(I,J)
  240           CONTINUE
              END IF
              DO 260 K = J + 1,N
                IF (A(K,J).NE.ZERO) THEN
                  DO 250 I = 1,M
                    B(I,J) = B(I,J) - A(K,J)*B(I,K)
  250             CONTINUE
                END IF
  260         CONTINUE
              IF (NOUNIT) THEN
                TEMP = ONE/A(J,J)
                DO 270 I = 1,M
                  B(I,J) = TEMP*B(I,J)
  270           CONTINUE
              END IF
  280       CONTINUE
          END IF
        ELSE
*
*           Form  B := alpha*B*inv( A' )
*           or    B := alpha*B*inv( conjg( A' ) ).
*
          IF (UPPER) THEN
            DO 330 K = N,1,-1
              IF (NOUNIT) THEN
                IF (NOCONJ) THEN
                  TEMP = ONE/A(K,K)
                ELSE
                  TEMP = ONE/DCONJG(A(K,K))
                END IF
                DO 290 I = 1,M
                  B(I,K) = TEMP*B(I,K)
  290           CONTINUE
              END IF
              DO 310 J = 1,K - 1
                IF (A(J,K).NE.ZERO) THEN
                  IF (NOCONJ) THEN
                    TEMP = A(J,K)
                  ELSE
                    TEMP = DCONJG(A(J,K))
                  END IF
                  DO 300 I = 1,M
                    B(I,J) = B(I,J) - TEMP*B(I,K)
  300             CONTINUE
                END IF
  310         CONTINUE
              IF (ALPHA.NE.ONE) THEN
                DO 320 I = 1,M
                  B(I,K) = ALPHA*B(I,K)
  320           CONTINUE
              END IF
  330       CONTINUE
          ELSE
            DO 380 K = 1,N
              IF (NOUNIT) THEN
                IF (NOCONJ) THEN
                  TEMP = ONE/A(K,K)
                ELSE
                  TEMP = ONE/DCONJG(A(K,K))
                END IF
                DO 340 I = 1,M
                  B(I,K) = TEMP*B(I,K)
  340           CONTINUE
              END IF
              DO 360 J = K + 1,N
                IF (A(J,K).NE.ZERO) THEN
                  IF (NOCONJ) THEN
                    TEMP = A(J,K)
                  ELSE
                    TEMP = DCONJG(A(J,K))
                  END IF
                  DO 350 I = 1,M
                    B(I,J) = B(I,J) - TEMP*B(I,K)
  350             CONTINUE
                END IF
  360         CONTINUE
              IF (ALPHA.NE.ONE) THEN
                DO 370 I = 1,M
                  B(I,K) = ALPHA*B(I,K)
  370           CONTINUE
              END IF
  380       CONTINUE
          END IF
        END IF
      END IF
*
      RETURN
*
*     End of ZTRSM .
*
      END
C
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
        CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
        BASE = BETA
        T = IT
        IF( LRND ) THEN
          RND = ONE
          EPS = ( BASE**( 1-IT ) ) / 2
        ELSE
          RND = ZERO
          EPS = BASE**( 1-IT )
        END IF
        PREC = EPS*BASE
        EMIN = IMIN
        EMAX = IMAX
        SFMIN = RMIN
        SMALL = ONE / RMAX
        IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
          SFMIN = SMALL*( ONE+EPS )
        END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
        RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
        RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
        RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
        RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
        RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
        RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
        RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
        RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
        RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
        RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      FIRST  = .FALSE.
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
        ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
        A = 1
        C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10   CONTINUE
        IF( C.EQ.ONE ) THEN
          A = 2*A
          C = DLAMC3( A, ONE )
          C = DLAMC3( C, -A )
          GO TO 10
        END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
        B = 1
        C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20   CONTINUE
        IF( C.EQ.A ) THEN
          B = 2*B
          C = DLAMC3( A, B )
          GO TO 20
        END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
        QTR = ONE / 4
        SAVEC = C
        C = DLAMC3( C, -A )
        LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
        B = LBETA
        F = DLAMC3( B / 2, -B / 100 )
        C = DLAMC3( F, A )
        IF( C.EQ.A ) THEN
          LRND = .TRUE.
        ELSE
          LRND = .FALSE.
        END IF
        F = DLAMC3( B / 2, B / 100 )
        C = DLAMC3( F, A )
        IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
        T1 = DLAMC3( B / 2, A )
        T2 = DLAMC3( B / 2, SAVEC )
        LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
        LT = 0
        A = 1
        C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30   CONTINUE
        IF( C.EQ.ONE ) THEN
          LT = LT + 1
          A = A*LBETA
          C = DLAMC3( A, ONE )
          C = DLAMC3( C, -A )
          GO TO 30
        END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      FIRST = .FALSE.
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
        ZERO = 0
        ONE = 1
        TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
        CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
        B = LBETA
        A = B**( -LT )
        LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
        B = TWO / 3
        HALF = ONE / 2
        SIXTH = DLAMC3( B, -HALF )
        THIRD = DLAMC3( SIXTH, SIXTH )
        B = DLAMC3( THIRD, -HALF )
        B = DLAMC3( B, SIXTH )
        B = ABS( B )
        IF( B.LT.LEPS )
     $      B = LEPS
*
        LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10   CONTINUE
        IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
          LEPS = B
          C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
          C = DLAMC3( HALF, -C )
          B = DLAMC3( HALF, C )
          C = DLAMC3( HALF, -B )
          B = DLAMC3( HALF, C )
          GO TO 10
        END IF
*+       END WHILE
*
        IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
        RBASE = ONE / LBETA
        SMALL = ONE
        DO 20 I = 1, 3
          SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20   CONTINUE
        A = DLAMC3( ONE, SMALL )
        CALL DLAMC4( NGPMIN, ONE, LBETA )
        CALL DLAMC4( NGNMIN, -ONE, LBETA )
        CALL DLAMC4( GPMIN, A, LBETA )
        CALL DLAMC4( GNMIN, -A, LBETA )
        IEEE = .FALSE.
*
        IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
          IF( NGPMIN.EQ.GPMIN ) THEN
            LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
          ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
            LEMIN = NGPMIN - 1 + LT
            IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
          ELSE
            LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
*
        ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
          IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
            LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
          ELSE
            LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
*
        ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
          IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
            LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
          ELSE
            LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
            IWARN = .TRUE.
          END IF
*
        ELSE
          LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
          IWARN = .TRUE.
        END IF
        FIRST = .FALSE.
***
* Comment out this if block if EMIN is ok
        IF( IWARN ) THEN
          FIRST = .TRUE.
          WRITE( 6, FMT = 9999 )LEMIN
        END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
        IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
        LRMIN = 1
        DO 30 I = 1, 1 - LEMIN
          LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30   CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
        CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
        EMIN = EMIN - 1
        A = B1
        B1 = DLAMC3( A / BASE, ZERO )
        C1 = DLAMC3( B1*BASE, ZERO )
        D1 = ZERO
        DO 20 I = 1, BASE
          D1 = D1 + B1
   20   CONTINUE
        B2 = DLAMC3( A*RBASE, ZERO )
        C2 = DLAMC3( B2 / RBASE, ZERO )
        D2 = ZERO
        DO 30 I = 1, BASE
          D2 = D2 + B2
   30   CONTINUE
        GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
        LEXP = TRY
        EXBITS = EXBITS + 1
        GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
        UEXP = LEXP
      ELSE
        UEXP = TRY
        EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
        EXPSUM = 2*LEXP
      ELSE
        EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
        EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
        EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
        Z = Z*RECBAS
        IF( Y.LT.ONE )
     $      OLDY = Y
        Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
        Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
C
C
C======================================================================
C
      INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
*
*  Purpose
*  =======
*
*       This program sets problem and machine dependent parameters
*       useful for xHSEQR and its subroutines. It is called whenever
*       ILAENV is called with 12 <= ISPEC <= 16
*
*  Arguments
*  =========
*
*       ISPEC  (input) integer scalar
*              ISPEC specifies which tunable parameter IPARMQ should
*              return.
*
*              ISPEC=12: (INMIN)  Matrices of order nmin or less
*                        are sent directly to xLAHQR, the implicit
*                        double shift QR algorithm.  NMIN must be
*                        at least 11.
*
*              ISPEC=13: (INWIN)  Size of the deflation window.
*                        This is best set greater than or equal to
*                        the number of simultaneous shifts NS.
*                        Larger matrices benefit from larger deflation
*                        windows.
*
*              ISPEC=14: (INIBL) Determines when to stop nibbling and
*                        invest in an (expensive) multi-shift QR sweep.
*                        If the aggressive early deflation subroutine
*                        finds LD converged eigenvalues from an order
*                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
*                        then the next QR sweep is skipped and early
*                        deflation is applied immediately to the
*                        remaining active diagonal block.  Setting
*                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
*                        multi-shift QR sweep whenever early deflation
*                        finds a converged eigenvalue.  Setting
*                        IPARMQ(ISPEC=14) greater than or equal to 100
*                        prevents TTQRE from skipping a multi-shift
*                        QR sweep.
*
*              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
*                        a multi-shift QR iteration.
*
*              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
*                        following meanings.
*                        0:  During the multi-shift QR sweep,
*                            xLAQR5 does not accumulate reflections and
*                            does not use matrix-matrix multiply to
*                            update the far-from-diagonal matrix
*                            entries.
*                        1:  During the multi-shift QR sweep,
*                            xLAQR5 and/or xLAQRaccumulates reflections and uses
*                            matrix-matrix multiply to update the
*                            far-from-diagonal matrix entries.
*                        2:  During the multi-shift QR sweep.
*                            xLAQR5 accumulates reflections and takes
*                            advantage of 2-by-2 block structure during
*                            matrix-matrix multiplies.
*                        (If xTRMM is slower than xGEMM, then
*                        IPARMQ(ISPEC=16)=1 may be more efficient than
*                        IPARMQ(ISPEC=16)=2 despite the greater level of
*                        arithmetic work implied by the latter choice.)
*
*       NAME    (input) character string
*               Name of the calling subroutine
*
*       OPTS    (input) character string
*               This is a concatenation of the string arguments to
*               TTQRE.
*
*       N       (input) integer scalar
*               N is the order of the Hessenberg matrix H.
*
*       ILO     (input) INTEGER
*       IHI     (input) INTEGER
*               It is assumed that H is already upper triangular
*               in rows and columns 1:ILO-1 and IHI+1:N.
*
*       LWORK   (input) integer scalar
*               The amount of workspace available.
*
*  Further Details
*  ===============
*
*       Little is known about how best to choose these parameters.
*       It is possible to use different values of the parameters
*       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
*
*       It is probably best to choose different parameters for
*       different matrices and different parameters at different
*       times during the iteration, but this has not been
*       implemented --- yet.
*
*
*       The best choices of most of the parameters depend
*       in an ill-understood way on the relative execution
*       rate of xLAQR3 and xLAQR5 and on the nature of each
*       particular eigenvalue problem.  Experiment may be the
*       only practical way to determine which choices are most
*       effective.
*
*       Following is a list of default values supplied by IPARMQ.
*       These defaults may be adjusted in order to attain better
*       performance in any particular computational environment.
*
*       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
*                        Default: 75. (Must be at least 11.)
*
*       IPARMQ(ISPEC=13) Recommended deflation window size.
*                        This depends on ILO, IHI and NS, the
*                        number of simultaneous shifts returned
*                        by IPARMQ(ISPEC=15).  The default for
*                        (IHI-ILO+1).LE.500 is NS.  The default
*                        for (IHI-ILO+1).GT.500 is 3*NS/2.
*
*       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
*
*       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
*                        a multi-shift QR iteration.
*
*                        If IHI-ILO+1 is ...
*
*                        greater than      ...but less    ... the
*                        or equal to ...      than        default is
*
*                                0               30       NS =   2+
*                               30               60       NS =   4+
*                               60              150       NS =  10
*                              150              590       NS =  **
*                              590             3000       NS =  64
*                             3000             6000       NS = 128
*                             6000             infinity   NS = 256
*
*                    (+)  By default matrices of this order are
*                         passed to the implicit double shift routine
*                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
*                         values of NS are used only in case of a rare
*                         xLAHQR failure.
*
*                    (**) The asterisks (**) indicate an ad-hoc
*                         function increasing from 10 to 64.
*
*       IPARMQ(ISPEC=16) Select structured matrix multiply.
*                        (See ISPEC=16 above for details.)
*                        Default: 3.
*
*     ================================================================
*     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 75, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            NH, NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
*     ..
*     .. Executable Statements ..
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
        NH = IHI - ILO + 1
        NS = 2
        IF( NH.GE.30 )
     $      NS = 4
        IF( NH.GE.60 )
     $      NS = 10
        IF( NH.GE.150 )
     $      NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )
        IF( NH.GE.590 )
     $      NS = 64
        IF( NH.GE.3000 )
     $      NS = 128
        IF( NH.GE.6000 )
     $      NS = 256
        NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
        IPARMQ = NMIN
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
        IPARMQ = NIBBLE
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
        IPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
        IF( NH.LE.KNWSWP ) THEN
          IPARMQ = NS
        ELSE
          IPARMQ = 3*NS / 2
        END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
        IPARMQ = 0
        IF( NS.GE.KACMIN )
     $      IPARMQ = 1
        IF( NS.GE.K22MIN )
     $      IPARMQ = 2
*
      ELSE
*        ===== invalid value of ispec =====
        IPARMQ = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
C
C
C======================================================================
C
      INTEGER FUNCTION IZAMAX(N,ZX,INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX ZX(*)
*     ..
*
*  Purpose
*  =======
*
*     finds the index of element having max. absolute value.
*     jack dongarra, 1/15/85.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION SMAX
      INTEGER I,IX
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL DCABS1
*     ..
      IZAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IZAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) GO TO 20
*
*        code for increment not equal to 1
*
      IX = 1
      SMAX = DCABS1(ZX(1))
      IX = IX + INCX
      DO 10 I = 2,N
        IF (DCABS1(ZX(IX)).LE.SMAX) GO TO 5
        IZAMAX = I
        SMAX = DCABS1(ZX(IX))
    5   IX = IX + INCX
   10 CONTINUE
      RETURN
*
*        code for increment equal to 1
*
   20 SMAX = DCABS1(ZX(1))
      DO 30 I = 2,N
        IF (DCABS1(ZX(I)).LE.SMAX) GO TO 30
        IZAMAX = I
        SMAX = DCABS1(ZX(I))
   30 CONTINUE
      RETURN
      END
C
C
C=======================================================================
C
      SUBROUTINE COUMAT(ITL,MI,LF,MF,
     1                  DIRIN,DIROU,QIN,QOU,
     2                  DELTA,RADIAL,MATRIX)
C
C       This routine calculates the spin-independent PhD optical matrix
C     elements for dipolar or quadrupole excitations. It is stored in
C     MATRIX(JDIR,JPOL,IO)
C
C  Here, the conventions are :
C
C             IPOL=1  : linearly polarized light
C             IPOL=2  : circularly polarized light
C
C             JPOL=1  : +/x polarization for circular/linear light
C             JPOL=2  : -/y polarization for circular/linear light
C
C  When IDICHR=0, JDIR = 1,2 and 3 correspond respectively to the x,y
C       and z directions for the linear polarization. But for IDICHR=1,
C       these basis directions are those of the position of the light.
C
C                                     Last modified : 20 Jan 2014
C
      INCLUDE 'spec.inc'
C
      COMMON /INIT_L/ LI,INITL,L2,L3,L4,L5
      COMMON /INIT_L_I/ LI_I,INITL_I,NNL1,LF1_I,LF2_I,ISTEP_LF_I
      COMMON /INIT_L_O/ LI_O,INITL_O,NNL2,LF1_O,LF2_O,ISTEP_LF_O
      COMMON /RXSINI/ THLUM_I,PHILUM_I,IPOL_I,NEPS_I,INTERACT_I
      COMMON /RXSFIN/ THLUM_O,PHILUM_O,IPOL_O,NEPS_O,INTERACT_O
      COMMON /SPIN/ I1,IDICHR,N1,N2,I2,I8,N3
      COMMON /TYPCAL/ I3,I4,I5,I6,I7,IPOL,I9,I10,I_TEST
C
      CHARACTER*7 INTERACT_I,INTERACT_O
C
      COMPLEX SUM_1,SUM_2,DELTA,YLM(3,-1:1),RADIAL
      COMPLEX MATRIX(3,2,2)
      COMPLEX ONEC,IC,IL,COEF,PROD
      COMPLEX PRODIN,PRODOU
C
      REAL RLM(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1),GNT(0:N_GAUNT)
      REAL THETA(3),PHI(3)
C
C
      REAL DIRIN(3),DIROU(3)
      REAL QIN,QOU
C
      INTEGER DQ
      REAL AIIN,AIOU,A
C
      COMPLEX AIFIN(3),AIFOU(3)
      COMPLEX SUM_1IN,SUM_1OU,SUM_2IN,SUM_2OU
C
C
      DATA PI4S3,C_LIN,SQR2 /4.188790,1.447202,1.414214/
      DATA PIS2 /1.570796/
      DATA SMALL /0.0001/
C
      ONEC=(1.,0.)
      IC=(0.,1.)
C
C
      DO JDIR=1,3
        DO I=1,2
          DO IO=1,2
            MATRIX(JDIR,I,IO)=(0.0,0.0)
          ENDDO
        ENDDO
      ENDDO
C
C
C
C          I do it here just for REXS test, it should be modified later
C                 by haifeng
C
      INITL=INITL_I
C
C
C
      IF(INITL.EQ.0) THEN
        DO JDIR=1,3
          DO J=1,2
            DO IO=1,2
              MATRIX(JDIR,J,IO)=ONEC
            ENDDO
          ENDDO
        ENDDO
      ELSEIF(ABS(INITL).EQ.1.OR.INITL.EQ.2) THEN
C
        M=MF-MI
C
        IF(MOD(LF,4).EQ.0) THEN
          IL=ONEC
        ELSEIF(MOD(LF,4).EQ.1) THEN
          IL=IC
        ELSEIF(MOD(LF,4).EQ.2) THEN
          IL=-ONEC
        ELSEIF(MOD(LF,4).EQ.3) THEN
          IL=-IC
        ENDIF
C
        CALL GAUNT(LI,MI,LF,MF,GNT)
C
        IF(ITL.EQ.0) THEN
c          COEF=CEXP(IC*DELTA)*CONJG(IL)
          COEF=CEXP(IC*DELTA)*IL
        ELSE
          IF(IDICHR.EQ.0) THEN
            COEF=PI4S3*IL
          ELSE
            COEF=C_LIN*IL
          ENDIF
        ENDIF
C
        PROD=COEF*RADIAL*GNT(1)
C
        IF(IDICHR.EQ.0) THEN
          YLM(1,-1)=(0.345494,0.)
          YLM(1,0)=(0.,0.)
          YLM(1,1)=(-0.345494,0.)
          YLM(2,-1)=(0.,-0.345494)
          YLM(2,0)=(0.,0.)
          YLM(2,1)=(0.,-0.345494)
          YLM(3,-1)=(0.,0.)
          YLM(3,0)=(0.488602,0.)
          YLM(3,1)=(0.,0.)
C
          DO JDIR=1,3
            DO IO=1,2
              MATRIX(JDIR,1,IO)=PROD*CONJG(YLM(JDIR,M))
            ENDDO
          ENDDO
C
        ELSEIF(IDICHR.GE.1) THEN
C
          THETA(1)=PIS2
          PHI(1)=0.
          THETA(2)=PIS2
          PHI(2)=PIS2
          THETA(3)=0.
          PHI(3)=0.
C
          DO JDIR=1,3
            CALL DJMN(THETA(JDIR),RLM,1)
            SUM_1=RLM(-1,M,1)*PROD*CEXP((0.,-1.)*M*PHI(JDIR))
            SUM_2=RLM(1,M,1)*PROD*CEXP((0.,-1.)*M*PHI(JDIR))
            IF(IPOL.EQ.2) THEN
              MATRIX(JDIR,1,1)=SQR2*SUM_1
              MATRIX(JDIR,2,1)=SQR2*SUM_2
              MATRIX(JDIR,1,2)=MATRIX(JDIR,1,1)
              MATRIX(JDIR,2,2)=MATRIX(JDIR,2,1)
            ELSEIF(ABS(IPOL).EQ.1) THEN
              MATRIX(JDIR,1,1)=(SUM_2-SUM_1)
              MATRIX(JDIR,2,1)=(SUM_2+SUM_1)*IC
              MATRIX(JDIR,1,2)=MATRIX(JDIR,1,1)
              MATRIX(JDIR,2,2)=MATRIX(JDIR,2,1)
            ENDIF
          ENDDO
        ENDIF
C
      ELSEIF(INITL.EQ.4) THEN
C
C
        M=MF-MI
C
        COEF=(1.0,0.0)
C
C
        DQ=ABS(LF-LI)
        CALL ANGINT_XYZ22(LI,MI,LF,MF,DIRIN,DQ,AIFIN)
        CALL ANGINT_XYZ22(LI,MI,LF,MF,DIROU,DQ,AIFOU)
        AIIN=0.0
        DO I=1,3
          A=ABS(AIFIN(I))
          AIIN=AIIN+A*A
        ENDDO
C
        AIOU=0.0
        DO I=1,3
          A=ABS(AIFOU(I))
          AIOU=AIOU+A*A
        ENDDO
C
        IF(AIIN.GT.SMALL.OR.AIOU.GT.SMALL) THEN

          PRODIN=COEF*RADIAL
          PRODOU=COEF*RADIAL
C
          IF (DQ.EQ.2) THEN
            PRODIN=PRODIN*IC*QIN
            PRODOU=PRODOU*IC*QOU
          ENDIF
C
          IF(IDICHR.EQ.0) THEN
            DO JDIR=1,3
              MATRIX(JDIR,1,1)=PRODIN*AIFIN(JDIR)
              MATRIX(JDIR,1,2)=PRODOU*AIFOU(JDIR)
              MATRIX(JDIR,1,2)=CONJG(MATRIX(JDIR,1,2))
            ENDDO
C
          ELSEIF(IDICHR.GE.1) THEN
C
            THETA(1)=PIS2
            PHI(1)=0.
            THETA(2)=PIS2
            PHI(2)=PIS2
            THETA(3)=0.
            PHI(3)=0.
C
            DO JDIR=1,3
              CALL DJMN(THETA(JDIR),RLM,1)
              SUM_1IN=RLM(-1,M,1)*PRODIN*CEXP((0.,-1.)*M*PHI(JDIR))
              SUM_1OU=RLM(-1,M,1)*PRODOU*CEXP((0.,-1.)*M*PHI(JDIR))
              SUM_2IN=RLM(1,M,1)*PRODIN*CEXP((0.,-1.)*M*PHI(JDIR))
              SUM_2OU=RLM(1,M,1)*PRODOU*CEXP((0.,-1.)*M*PHI(JDIR))
              IF(IPOL.EQ.2) THEN
                MATRIX(JDIR,1,1)=SQR2*SUM_1IN
                MATRIX(JDIR,1,2)=SQR2*SUM_1OU
                MATRIX(JDIR,2,1)=SQR2*SUM_2IN
                MATRIX(JDIR,2,2)=SQR2*SUM_2OU
              ELSEIF(ABS(IPOL).EQ.1) THEN
                MATRIX(JDIR,1,1)=(SUM_2IN-SUM_1IN)
                MATRIX(JDIR,1,2)=(SUM_2OU-SUM_1OU)
                MATRIX(JDIR,2,1)=(SUM_2IN+SUM_1IN)*IC
                MATRIX(JDIR,2,2)=(SUM_2OU+SUM_1OU)*IC
              ENDIF
              DO I=1,2
                MATRIX(JDIR,I,2)=CONJG(MATRIX(JDIR,I,2))
              ENDDO
            ENDDO
          ENDIF
C
        ENDIF !IF AI.GT.SMALL/
C
      ELSE
        DO JDIR=1,3
          DO J=1,2
            DO IO=1,2
              MATRIX(JDIR,J,IO)=ONEC
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
C
      END
C
C
C======================================================================
C
      SUBROUTINE ANGINT_XYZ22(L1,M1,L2,M2,DIR,DQ0,F)
C
C       This subroutine calculates the angular integral for radial
C     matrix with dipole approximation (DQ=1) or quadrupole
C     approximation (DQ=2).
C       note that Spherical harmonic with (L2,M2) has conjugation.
C
C         DIR(1)=SIN(THETA2)*SIN(PHI2)
C         DIR(2)=SIN(THETA2)*SIN(PHI2)
C         DIR(3)=COS(THETA2)
C
C       when dipole approximation is used, THETA2 and PHI2 are not used.
C
C                                     last modified 2014.12.20
C
C
      INCLUDE 'spec.inc'
C
      INTEGER L1,M1,L2,M2,DQ0
      REAL DIR(3)
      COMPLEX F(3)
C
      REAL GNT(0:N_GAUNT)
      INTEGER DQ,DM
      COMPLEX IC
      COMPLEX F1,F2,F3,F4,F5,F6
C
      DATA PI,SQR2 /3.141593,1.414214/
C
C
      IC=(0.,1.)
C
      DO I=1,3
        F(I)=(0.0,0.0)
      ENDDO
C
      DQ=DQ0
      IF(DQ.NE.1.AND.DQ.NE.2) THEN
        WRITE(*,2000)
 2000   FORMAT(//,2X,'WARNING!, only dipole or quadrupole',
     1            'approximation can be treated!')
        IF(ABS(DQ).LT.1) THEN
          DQ=1
          WRITE(*,2100)
 2100     FORMAT(/,2X,'Change to Dipole approximation calculation')
        ELSE
          DQ=2
          WRITE(*,2200)
 2200     FORMAT(/,2X,'Change to Quadrupole approximation calculation')
        ENDIF
      ENDIF
C
      CALL GAUNT(L2,M2,L1,M1,GNT)
C
      IF(DQ.EQ.1) THEN
        SQR2PIS3=SQRT(2.0*PI/3.0)
C
        DM=M1-M2
        IF(DM.EQ.1) THEN
          A1=GNT(1)
          A2=0.0
          A3=0.0
        ELSEIF(DM.EQ.0) THEN
          A1=0.0
          A2=GNT(1)
          A3=0.0
        ELSEIF(DM.EQ.-1) THEN
          A1=0.0
          A2=0.0
          A3=GNT(1)
        ELSE
          A1=0.0
          A2=0.0
          A3=0.0
        ENDIF
C
        F(1)=(A3-A1)
        F(2)=-IC*(A3+A1)
        F(3)=SQR2*A2
        DO I=1,3
          F(I)=SQR2PIS3*F(I)
        ENDDO
C
      ELSE
C
        DM=M1-M2
        IF(DM.EQ.2) THEN
          A1=GNT(2)
          A2=0.0
          A3=0.0
          A4=0.0
          A5=0.0
          B=0.0
        ELSEIF(DM.EQ.1) THEN
          A1=0.0
          A2=GNT(2)
          A3=0.0
          A4=0.0
          A5=0.0
          B=0.0
        ELSEIF(DM.EQ.0) THEN
          A1=0.0
          A2=0.0
          A3=GNT(2)
          A4=0.0
          A5=0.0
          B=GNT(0)
        ELSEIF(DM.EQ.-1) THEN
          A1=0.0
          A2=0.0
          A3=0.0
          A4=GNT(2)
          A5=0.0
          B=0.0
        ELSEIF(DM.EQ.-2) THEN
          A1=0.0
          A2=0.0
          A3=0.0
          A4=0.0
          A5=GNT(2)
          B=0.0
        ELSE
          A1=0.0
          A2=0.0
          A3=0.0
          A4=0.0
          A5=0.0
          B=0.0
        ENDIF
C
        C1=SQRT(2.0/15.0)
        C3=2.0/3.0
        C2=SQRT(C3)
        C4=1.0/3.0
C
        F1=C1*A1-C1*C2*A3+C1*A5+C3*B
        F1=F1*DIR(1)
C
        F2=C1*A1+C1*C2*A3+C1*A5-C3*B
        F2=F2*(-1.0)*DIR(2)
C
        F3=C4*B+C1*C2*A3
        F3=F3*2.0*DIR(3)
C
        F4=C1*(A1-A5)
C
        F5=C1*(A4-A2)
C
        F6=C1*(A2+A4)
C
        F(1)=F1+F4*DIR(2)*IC+F5*DIR(3)
        F(2)=F2+F4*DIR(1)*IC+F6*DIR(3)*(-1.0)*IC
        F(3)=F3+F5*DIR(1)+F6*DIR(2)*(-1.0)*IC
        DO I=1,3
          F(I)=SQRT(PI)*F(I)
        ENDDO
      ENDIF
C
      RETURN
C
      END
C
C
C======================================================================
C
      SUBROUTINE XYZ2THPH(A0,THETA,PHI)
C
C       This subrouitne transfers the Cartesian coordinates A0(3) to
C     polar angles (theta,phi)
C     X=A0(1),Y=A0(2),Z=A0(3)
C
C                                            last modified 2013.12.20
C
      REAL A0(3),THETA,PHI
C
      REAL B,A(3)
      REAL SMALL
C
      DATA PI /3.141593/
      DATA SMALL /0.0001/
C
      THETA=0.0
      PHI=0.0
C
      B=0.0
      DO I=1,3
        B=B+A0(I)*A0(I)
      ENDDO
C
      IF(B.LE.SMALL) THEN
        THETA=0.0
        PHI=0.0
        WRITE(*,1000) THETA,PHI
 1000   FORMAT(//,2X,'WARNING! THETA AND PHI SET TO BE (',F6.3,','
     1                                                   ,F6.3,')',//)
      ELSE
        DO I=1,3
          A(I)=A0(I)/B
        ENDDO
C
        THETA=ACOS(A(3))
        IF(A(1).LE.SMALL) THEN
          IF(A(2).GT.0.0) THEN
            PHI=PI/2.0
          ELSE
            PHI=-1.0*PI/2.0
          ENDIF
        ELSE
          B=A(2)/A(1)
          PHI=ATAN(B)
          IF(A(1).LT.0.0) THEN
            IF(A(2).GT.0.0) THEN
              PHI=PHI+PI
            ELSE
              PHI=PHI-PI
            ENDIF
          ENDIF            
        ENDIF
      ENDIF
C
      RETURN
C
      END
C
C
C======================================================================
C
      SUBROUTINE RESDIF_MI(NPLAN,VAL,ZEM,IPHA,NAT2,COORD,NATYP,RHOK,
     1                     NATCLU,NFICHLEC,JFICH,NP)
C
C
C
C     the subroutines is modifed based on PHD which comes from
C    file 'phd_note.f'
C
C
C   This subroutine computes the REXS formula in the spin-independent case
C        from a non spin-orbit resolved initial core state LI.
C
C   The calculation is performed using a matrix inversion for the
C           expression of the scattering path operator
C
C   The matrix inversion is performed using the LAPACK inversion
C          routines for a general complex matrix
C
C                                        Last modified : 15 Jan 2014
C
      INCLUDE 'spec.inc'
C
      REAL LUM(3),AXE(3),EPS(3),DIRLUM(3),E_PH(NE_M)
C
      COMPLEX IC,ONEC,ZEROC,VK,DLT,TL,COEF,PW(0:NDIF_M),DELTA
      COMPLEX TLT(0:NT_M,4,NATM,NE_M),RHO01,RHOMI
      COMPLEX TAU(LINMAX,LINFMAX,NATCLU_M)
      COMPLEX YLMR(0:NL_M,-NL_M:NL_M),MATRIX(3,2,2)
      COMPLEX YLME(0:NL_M,-NL_M:NL_M)
      COMPLEX R2,MLFLI(2,-LI_M:LI_M,3,2,3,2)
      COMPLEX SJDIR_1,SJDIR_2,SJDIF_1,SJDIF_2
      COMPLEX RHOK(NE_M,NATM,0:18,5,NSPIN2_M),RD
      COMPLEX SLJDIF,ATT_M,MLIL0(2,-LI_M:LI_M,6,2),SLF_1,SLF_2
      COMPLEX SDIRTHM_1,SDIRTHM_2
      COMPLEX SL0DIF,SMJDIF
C
      DIMENSION VAL(NATCLU_M),NATYP(NATM),DIRPOL(3,2,2)
      DIMENSION EMET(3),R_L(9),COORD(3,NATCLU_M)
      DIMENSION R(NDIF_M),XR(NDIF_M),YR(NDIF_M),ZR(NDIF_M)
      DIMENSION JPOS(NDIF_M,3),JPA(NDIF_M)
C
      CHARACTER*2 ALGO1,ALGO2,ALGO3,ALGO4
      CHARACTER*3 UNIT,S_O,SPECTRO,STEREO
      CHARACTER*7 INTERACT,STAT
      CHARACTER*7 INTERACT_I,INTERACT_O
      CHARACTER*13 OUTDATA1,OUTDATA2
      CHARACTER*24 INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6
      CHARACTER*24 INFILE7,INFILE8,INFILE9
      CHARACTER*24 OUTFILE,OUTFILE1,OUTFILE2,OUTFILE3,OUTFILE4
      CHARACTER*24 AMPFILE
C
      COMPLEX SMIDIR_1,SMIDIR_2,SMIDIF_1,SMIDIF_2
      REAL AA(3)
      REAL DIRQIN(3),DIRQOU(3)
C
      COMMON /ALGORITHM/ ALGO1,ALGO2,ALGO3,ALGO4
      COMMON /AMPLI/ I_AMP
      COMMON /APPROX/ NDIF,NO,ISPHER,IFWD,NTHOUT,RTHFWD(NATP_M),
     1                IBWD(NATP_M),RTHBWD(NATP_M),IPW,NCUT,PCTINT,IPP,
     2                ISPEED,IATTS,ILENGTH,RLENGTH
      COMMON /COOR/ NTCLU,N_PROT,NTP(NATM),NCHTYP(NATP_M),
     1              NCORR(NAT_EQ_M,NATP_M),INEW_AT(NATCLU_M),
     2              SYM_AT(3,NATCLU_M)
      COMMON /DEBWAL/ IDCM,IDWSPH,TD,QD,TEMP,RSJ,UJ2(NATM)
      COMMON /DIRECT/ DIRANA(3,49),ANADIR(3,49),RTHETA,
     1                RPHI,THETAR(49),PHIR(49)
      COMMON /EXTREM/ FMIN(0:NDIF_M),FMAX(0:NDIF_M),IREF
      COMMON /FIXSCAN/ N_FIXED,N_SCAN,IPH_1,FIX0,FIX1,SCAN0,SCAN1
      COMMON /INFILES/ INFILE1,INFILE2,INFILE3,INFILE4,INFILE5,INFILE6,
     1                 INFILE7,INFILE8,INFILE9
      COMMON /INUNITS/ IUI1,IUI2,IUI3,IUI4,IUI5,IUI6,IUI7,IUI8,IUI9
      COMMON /INIT_L/ LI,INITL,NNL,LF1,LF2,ISTEP_LF
      COMMON /INIT_J/ JF1,JF2,I_SO,S_O
      COMMON /LIMAMA/ NIV,COUPUR
      COMMON /MOYEN/ IMOY,NDIR,ACCEPT,ICHKDIR
      COMMON /OUTFILES/ OUTFILE1,OUTFILE2,OUTFILE3,OUTFILE4
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /PARCAL/ NPHI0,NE0,NTHETA0,NFTHET0,NEPS0
      COMMON /RXSGEN/ IE_RX,IPHI_RX,ITHETA_RX,NE,NPHI,NTHETA,
     1                E_INI_RX,E_FIN_RX,PHI_FIN_RX,THETA_FIN_RX
      COMMON /RXSINI/ THLUM_I,PHILUM_I,IPOL_I,NEPS_I,INTERACT_I
      COMMON /RXSFIN/ THLUM_O,PHILUM_O,IPOL_O,NEPS_O,INTERACT_O
C
      COMMON /RESEAU/ NCRIST,NCENTR,IBAS,NAT,A,BSURA,CSURA,UNIT
      COMMON /SPIN/ ISPIN,IDICHR,NSPIN,NSPIN2,ISFLIP,IR_DIA,NSTEP
      COMMON /TESTPB/ RHO01,TH01,PHI01
      COMMON /TESTS/ ITEST,IPRINT,ISORT1,NPATHP,ISOM
      COMMON /TRANS/ DLT(NE_M,NATM,0:18,5),TL(0:NT_M,4,NATM,NE_M),
     1               VK(NE_M),VK2(NE_M),IPOTC,ITL,LMAX(NATM,NE_M)
      COMMON /TYPCAL/ IPHI,IE,ITHETA,IFTHET,IMOD,IPOL,I_CP,I_EXT,I_TEST
      COMMON /TYPEM/ NEMET,IESURF,IEMET(NEMET_M)
      COMMON /TYPEXP/ SPECTRO,INTERACT,STEREO
      COMMON /VALIN/ PHI0,E0,THETA0,THLUM,PHLUM,ELUM,VINT,NONVOL(100)
      COMMON /VALIN_AV/ I_SET,TH_0(NTH_M),PH_0(NPH_M)
      COMMON /VALFIN/ PHI1,EFIN,THETA1
C
      DATA PI,PIS180,CONV /3.141593,0.017453,0.512314/
C
C       For light, when E=1eV, k=5.0677288*10^(-4)A^(-1)
C
      DATA CONVPH1,CONVPH2/5.0677288,10000.0/
      DATA FINSTRUC,CVECT,SMALL /0.007297,1.0,0.0001/
C
C      R0E is the classic radius for electron which is 2.81794*10^(-15)m 
C      I use anstrone unit here and R0E=R0E1*R0E2, the reason why I write
C     in this way is for the storage problem.
C
C
      REAL R0E1,R0E2
C
C      DATA R0E1 R0E2 /2.81794,0.00001/
C
      DATA R0E1,R0E2/1.0,1.0/
C
      ALGO1='MI'
      ALGO2='  '
      ALGO3='  '
      ALGO4='  '
C
      I_DIR=0
      NSET=1
      JEL=1
      OUTDATA1='CROSS-SECTION'
C
      IF(I_AMP.EQ.1) THEN
        I_MI=1
        OUTDATA2='MS AMPLITUDES'
      ELSE
        I_MI=0
      ENDIF
C
      IF(SPECTRO.EQ.'RES') THEN
        IOUT=IUO2
        OUTFILE=OUTFILE2
        STAT='UNKNOWN'
        IF(I_MI.EQ.1) THEN
          IOUT2=IUSCR2+1
          N_DOT=1
          DO J_CHAR=1,24
            IF(OUTFILE(J_CHAR:J_CHAR).EQ.'.') GOTO 888
            N_DOT=N_DOT+1
          ENDDO
  888     CONTINUE
          AMPFILE=OUTFILE(1:N_DOT)//'amp'
          OPEN(UNIT=IOUT2, FILE=AMPFILE, STATUS=STAT)
        ENDIF
      ENDIF
C
      RTHLUM=THLUM*PIS180
      RPHLUM=PHLUM*PIS180
      X_LUM_Z=SIN(RTHLUM)*COS(RPHLUM)
      Y_LUM_Z=SIN(RTHLUM)*SIN(RPHLUM)
      Z_LUM_Z=COS(RTHLUM)
C
      IF(IMOD.EQ.0) THEN
        DIRLUM(1)=X_LUM_Z
        DIRLUM(2)=Y_LUM_Z
        DIRLUM(3)=Z_LUM_Z
      ELSE
        IF(I_EXT.EQ.0) THEN
          RTH0=THETA0*PIS180
          RPH0=PHI0*PIS180
          RTH=RTH0
          RPH=RPH0
          R_L(1)=COS(RTH0)*COS(RPH0)*COS(RPH0)+SIN(RPH0)*SIN(RPH0)
          R_L(2)=COS(RTH0)*SIN(RPH0)*COS(RPH0)-SIN(RPH0)*COS(RPH0)
          R_L(3)=SIN(RTH0)*COS(RPH0)
          R_L(4)=COS(RTH0)*SIN(RPH0)*COS(RPH0)-SIN(RPH0)*COS(RPH0)
          R_L(5)=COS(RTH0)*SIN(RPH0)*SIN(RPH0)+COS(RPH0)*COS(RPH0)
          R_L(6)=SIN(RTH0)*SIN(RPH0)
          R_L(7)=-SIN(RTH0)*COS(RPH0)
          R_L(8)=-SIN(RTH0)*SIN(RPH0)
          R_L(9)=COS(RTH0)
          LUM(1)=X_LUM_Z*R_L(1)+Y_LUM_Z*R_L(2)+Z_LUM_Z*R_L(3)
          LUM(2)=X_LUM_Z*R_L(4)+Y_LUM_Z*R_L(5)+Z_LUM_Z*R_L(6)
          LUM(3)=X_LUM_Z*R_L(7)+Y_LUM_Z*R_L(8)+Z_LUM_Z*R_L(9)
C
        ENDIF
      ENDIF
C
      IC=(0.,1.)
      ONEC=(1.,0.)
      ZEROC=(0.,0.)
      NSCAT=NATCLU-1
      ATTSE=1.
      ATTSJ=1.
      ZSURF=VAL(1)
      IF((ISOM.EQ.0).OR.(JFICH.EQ.1)) THEN
        OPEN(UNIT=IOUT, FILE=OUTFILE, STATUS=STAT)
      ENDIF
C
      CALL HEADERS(IOUT)
C
      IF((ISOM.EQ.0).OR.((ISOM.GT.0).AND.(JFICH.EQ.1))) THEN
        WRITE(IOUT,12) SPECTRO,OUTDATA1
        WRITE(IOUT,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE,
     1                IPH_1,I_EXT
        IF(I_MI.EQ.1) THEN
          WRITE(IOUT2,12) SPECTRO,OUTDATA2
          WRITE(IOUT2,12) STEREO
          WRITE(IOUT2,19) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,
     1                    ITHETA,IE,IPH_1,I_EXT
          WRITE(IOUT2,20) PHI0,THETA0,PHI1,THETA1,NONVOL(1)
        ENDIF
      ENDIF
C
      IF(ISOM.EQ.0) THEN
        WRITE(IOUT,79) NPLAN,NEMET,NTHETA,NPHI,NE
        IF(I_MI.EQ.1) THEN
          WRITE(IOUT2,79) NPLAN,NEMET,NTHETA,NPHI,NE
        ENDIF
      ELSEIF((ISOM.NE.0).AND.(JFICH.EQ.1)) THEN
        WRITE(IOUT,11) NTHETA,NPHI,NE
        IF(I_MI.EQ.1) THEN
          WRITE(IOUT2,11) NTHETA,NPHI,NE
        ENDIF
      ENDIF
      IJK=0
C
C   Loop over the planes
C
      DO JPLAN=1,NPLAN
        Z=VAL(JPLAN)
        IF((IPHA.EQ.1).OR.(IPHA.EQ.2)) THEN
          DZZEM=ABS(Z-ZEM)
          IF(DZZEM.LT.SMALL) GOTO 10
          GOTO 1
        ENDIF
 10     CONTINUE
C
C   Loop over the different absorbers in a given plane
C
        DO JEMET=1,NEMET
          CALL EMETT(JEMET,IEMET,Z,SYM_AT,NATYP,EMET,NTYPEM,
     1             JNEM,*4)
          GO TO 2
C
   4      IF((ISORT1.EQ.0).AND.(IPRINT.GT.0)) THEN
            IF(I_TEST.NE.2) WRITE(IUO1,51) JPLAN,NTYPEM
          ENDIF
          GO TO 3
C
   2      IF((ABS(EMET(3)).GT.COUPUR).AND.(IBAS.EQ.1)) GOTO 5
          IF((ISORT1.EQ.0).AND.(IPRINT.GT.0)) THEN
            IF(I_TEST.NE.2) THEN
              WRITE(IUO1,52) JPLAN,EMET(1),EMET(2),EMET(3),NTYPEM
            ENDIF
          ENDIF
          IF(ISOM.EQ.1) NP=JPLAN
          ZSURFE=VAL(1)-EMET(3)
          JATLEM=JNEM
C
C   Loop over the energies
C
          DO JE=1,NE
            FMIN(0)=1.
            FMAX(0)=1.
            IF(NE.GT.1) THEN
              ECIN=E0+FLOAT(JE-1)*(EFIN-E0)/FLOAT(NE-1)
              E_PH(JE)=ELUM+FLOAT(JE-1)*(EFIN-E0)/FLOAT(NE-1)
            ELSEIF(NE.EQ.1) THEN
              ECIN=E0
              E_PH(JE)=ELUM
            ENDIF
C
            QIN=E_PH(JE)*CONVPH1*A/CONVPH2
            QOU=QIN
C
            IF(I_TEST.NE.1) THEN
              CFM=(-0.5)*E_PH(JE)*E_PH(JE)
              CFM=CFM/(-4.0*PI)
            ELSE
              CFM=1.
            ENDIF
C
            CALL LPM(ECIN,XLPM,*6)
            XLPM1=XLPM/A
C
            IF(IPRINT.GT.0) WRITE(IUO1,56) A,XLPM1
            IF((IPRINT.GT.0).AND.(IBAS.EQ.1)) THEN
              IF(I_TEST.NE.2) WRITE(IUO1,57) COUPUR
            ENDIF
C
            IF(ITL.EQ.0) THEN
              VK(JE)=SQRT(ECIN+VINT)*CONV*A*(1.,0.)
              VK2(JE)=CABS(VK(JE)*VK(JE))
            ENDIF
C
            GAMMA=1./(2.*XLPM1)
            IF(IPOTC.EQ.0) THEN
              VK(JE)=VK(JE)+IC*GAMMA
            ENDIF
C
            IF(I_TEST.NE.1) THEN
              VKR=REAL(VK(JE))
            ELSE
              VKR=1.
            ENDIF
C
            IF(I_MI.EQ.1) THEN
              WRITE(IOUT2,21) ECIN,VKR*CFM
            ENDIF
C
            IF((IDWSPH.EQ.1).AND.(ISPEED.EQ.1)) THEN
              IF(IDCM.GE.1) WRITE(IUO1,22)
              DO JAT=1,N_PROT
                IF(IDCM.EQ.0) THEN
                  XK2UJ2=VK2(JE)*UJ2(JAT)
                ELSE
                  XK2UJ2=VK2(JE)*UJ_SQ(JAT)
                  WRITE(IUO1,23) JAT,UJ_SQ(JAT)*A*A
                ENDIF
                CALL DWSPH(JAT,JE,XK2UJ2,TLT,ISPEED)
                DO LAT=0,LMAX(JAT,JE)
                  TL(LAT,1,JAT,JE)=TLT(LAT,1,JAT,JE)
                ENDDO
              ENDDO
            ENDIF
C
            IF(ABS(I_EXT).GE.1) THEN
              OPEN(UNIT=IUI6, FILE=INFILE6, STATUS='OLD')
              READ(IUI6,13) I_DIR,NSET,N_DUM1
              READ(IUI6,14) I_DUM1,N_DUM2,N_DUM3
            ENDIF
C
C   Initialization of TAU(INDJ,LINFMAX,JTYP)
C
            JATL=0
            DO JTYP=1,N_PROT
              NBTYP=NATYP(JTYP)
              LMJ=LMAX(JTYP,JE)
              DO JNUM=1,NBTYP
                JATL=JATL+1
                DO LF=LF1,LF2,ISTEP_LF
                  ILF=LF*LF+LF+1
                  DO MF=-LF,LF
                    INDF=ILF+MF
                    DO LJ=0,LMJ
                      ILJ=LJ*LJ+LJ+1
                      DO MJ=-LJ,LJ
                        INDJ=ILJ+MJ
                        TAU(INDJ,INDF,JATL)=ZEROC
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
C
C
C   Matrix inversion for the calculation of TAU
C
            IF(I_TEST.EQ.2) GOTO 666
C
            CALL INV_MAT_MS(JE,TAU)
C
 666        CONTINUE
C
C
C    Calculation of the REXS formula
C
C
C    Loop over the 'fixed' angle
C
  15        DO J_FIXED=1,N_FIXED
              IF(N_FIXED.GT.1) THEN
                IF(I_EXT.EQ.0) THEN
                  FIX_STEP=(FIX1-FIX0)/FLOAT(N_FIXED-1)
                  XINCRF=FLOAT(J_FIXED-1)*FIX_STEP
                ELSE
                  XINCRF=0.
                ENDIF
              ELSEIF(N_FIXED.EQ.1) THEN
                XINCRF=0.
              ENDIF
C
              IF(ABS(I_EXT).GE.1) THEN
                READ(IUI6,86) JSET,JLINE,THD,PHD
                IF(I_EXT.EQ.-1) BACKSPACE IUI6
                THETA0=THD
                PHI0=PHD
              ENDIF
C
              IF(IPH_1.EQ.1) THEN
                IF(I_EXT.EQ.0) THEN
                  DPHI=PHI0+XINCRF
                ELSE
                  DPHI=PHD
                ENDIF
                RPHI=DPHI*PIS180
                IF(IPRINT.GT.0) WRITE(IUO1,66) DPHI
              ELSE
                ISAUT=0
                IF(I_EXT.EQ.0) THEN
                  DTHETA=THETA0+XINCRF
                ELSE
                  DTHETA=THD
                ENDIF
                RTHETA=DTHETA*PIS180
                IF(ABS(DTHETA).GT.90.) ISAUT=ISAUT+1
                IF(I_EXT.GE.1) ISAUT=0
                IF(I_TEST.EQ.2)  ISAUT=0
                IF(ISAUT.GT.0) GOTO 8
                IF(IPRINT.GT.0) WRITE(IUO1,65) DTHETA
                IF((IPRINT.GT.0).AND.(I_TEST.NE.2)) WRITE(IUO1,59)
                IF((IPRINT.EQ.1).AND.(I_TEST.NE.2)) WRITE(IUO1,60)
                IF(STEREO.EQ.'YES') THEN
                  N_SCAN=INT((SCAN1-SCAN0)*SIN(RTHETA)/FIX_STEP+SMALL)+1
                ENDIF
              ENDIF
C
              IF((N_FIXED.GT.1).AND.(IMOD.EQ.1)) THEN
                IF(IPH_1.EQ.1) THEN
                  RTH=THETA0*PIS180
                  RPH=RPHI
                ELSE
                  RTH=RTHETA
                  RPH=PHI0*PIS180
                ENDIF
C
                R_L(1)=COS(RTH)*COS(RPH)
                R_L(2)=-SIN(RPH)
                R_L(3)=SIN(RTH)*COS(RPH)
                R_L(4)=COS(RTH)*SIN(RPH)
                R_L(5)=COS(RPH)
                R_L(6)=SIN(RTH)*SIN(RPH)
                R_L(7)=-SIN(RTH)
                R_L(8)=0.
                R_L(9)=COS(RTH)
C
                LUM(1)=X_LUM_Z*R_L(1)+Y_LUM_Z*R_L(2)+Z_LUM_Z*R_L(3)
                LUM(2)=X_LUM_Z*R_L(4)+Y_LUM_Z*R_L(5)+Z_LUM_Z*R_L(6)
                LUM(3)=X_LUM_Z*R_L(7)+Y_LUM_Z*R_L(8)+Z_LUM_Z*R_L(9)
              ENDIF
C
C    Loop over the scanned angle
C
              DO J_SCAN=1,N_SCAN
                IF(N_SCAN.GT.1) THEN
                  XINCRS=FLOAT(J_SCAN-1)*(SCAN1-SCAN0)/FLOAT(N_SCAN-1)
                ELSEIF(N_SCAN.EQ.1) THEN
                  XINCRS=0.
                ENDIF
                IF(I_EXT.EQ.-1) THEN
                  READ(IUI6,86) JSET,JLINE,THD,PHD
                  BACKSPACE IUI6
                ENDIF
C
                IF(IPH_1.EQ.1) THEN
                  ISAUT=0
                  IF(I_EXT.EQ.0) THEN
                    DTHETA=THETA0+XINCRS
                  ELSE
                    DTHETA=THD
                  ENDIF
                  RTHETA=DTHETA*PIS180
                  IF(ABS(DTHETA).GT.90.) ISAUT=ISAUT+1
                  IF(I_EXT.GE.1) ISAUT=0
                  IF(I_TEST.EQ.2)  ISAUT=0
                  IF(ISAUT.GT.0) GOTO 8
                  IF(IPRINT.GT.0) WRITE(IUO1,65) DTHETA
                  IF((IPRINT.GT.0).AND.(I_TEST.NE.2)) WRITE(IUO1,59)
                  IF((IPRINT.EQ.1).AND.(I_TEST.NE.2)) WRITE(IUO1,60)
                ELSE
                  IF(I_EXT.EQ.0) THEN
                    DPHI=PHI0+XINCRS
                  ELSE
                    DPHI=PHD
                  ENDIF
                  RPHI=DPHI*PIS180
                  IF(IPRINT.GT.0) WRITE(IUO1,66) DPHI
                ENDIF
C
C  Loop over the sets of directions to average over (for gaussian average)
C
C
                SSETDIR_1=0.
                SSETDIF_1=0.
                SSETDIR_2=0.
                SSETDIF_2=0.
C
                SSET2DIR_1=0.
                SSET2DIF_1=0.
                SSET2DIR_2=0.
                SSET2DIF_2=0.
C
                IF(I_EXT.EQ.-1) THEN
                  JREF=INT(NSET)/2+1
                ELSE
                  JREF=1
                ENDIF
C
                DO J_SET=1,NSET
                  IF(I_EXT.EQ.-1) THEN
                    READ(IUI6,86) JSET,JLINE,THD,PHD,W
                    DTHETA=THD
                    DPHI=PHD
                    RTHETA=DTHETA*PIS180
                    RPHI=DPHI*PIS180
C
                    RTH=TH_0(J_SET)*PIS180
                    RPH=PH_0(J_SET)*PIS180
C
                    IF(IMOD.EQ.1) THEN
                      R_L(1)=COS(RTH)*COS(RPH)
                      R_L(2)=-SIN(RPH)
                      R_L(3)=SIN(RTH)*COS(RPH)
                      R_L(4)=COS(RTH)*SIN(RPH)
                      R_L(5)=COS(RPH)
                      R_L(6)=SIN(RTH)*SIN(RPH)
                      R_L(7)=-SIN(RTH)
                      R_L(8)=0.
                      R_L(9)=COS(RTH)
C
                      LUM(1)=X_LUM_Z*R_L(1)+Y_LUM_Z*R_L(2)
     1                                     +Z_LUM_Z*R_L(3)
                      LUM(2)=X_LUM_Z*R_L(4)+Y_LUM_Z*R_L(5)
     1                                     +Z_LUM_Z*R_L(6)
                      LUM(3)=X_LUM_Z*R_L(7)+Y_LUM_Z*R_L(8)
     1                                     +Z_LUM_Z*R_L(9)
                    ENDIF
                  ELSE
                    W=1.
                  ENDIF
C
                  IF(I_EXT.EQ.-1) PRINT 89
C
                  CALL DIRAN(VINT,ECIN,JEL)
C
                  IF(J_SET.EQ.JREF) THEN
                    DTHETAP=DTHETA
                    DPHIP=DPHI
                  ENDIF
C
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO1,88) DTHETA,DPHI
                  ENDIF
C
C
C     ..........          Case IMOD=1 only          ..........
C
C  Calculation of the position of the light when the analyzer is at
C   (THETA,PHI). DIRLUM is the direction of the light and its initial
C   value (at (THETA0,PHI0)) is LUM. AXE is the direction of the theta
C   rotation axis and EPS is defined so that (AXE,DIRLUM,EPS) is a
C   direct orthonormal basis. The transform of a vector R by a rotation
C   of OMEGA about AXE is then given by
C
C     R' = R COS(OMEGA) + (AXE.R)(1-COS(OMEGA)) AXE + (AXE^R) SIN(OMEGA)
C
C   Here, DIRANA is the internal direction of the analyzer and ANADIR
C                      its external position
C
C   Note that when the initial position of the analyzer is (RTH,RPH)
C    which coincides with (RTH0,RPH0) only for the first fixed angle
C
                  IF(IMOD.EQ.1) THEN
                    IF(ITHETA.EQ.1) THEN
                      AXE(1)=-SIN(RPH)
                      AXE(2)=COS(RPH)
                      AXE(3)=0.
                      RANGLE=RTHETA-RTH
                    ELSEIF(IPHI.EQ.1) THEN
                      AXE(1)=0.
                      AXE(2)=0.
                      AXE(3)=1.
                      RANGLE=RPHI-RPH
                    ENDIF
                    CALL PRVECT(AXE,LUM,EPS,CVECT)
                    PRS=PRSCAL(AXE,LUM)
                    IF(J_SCAN.EQ.1) THEN
                      DIRLUM(1)=LUM(1)
                      DIRLUM(2)=LUM(2)
                      DIRLUM(3)=LUM(3)
                    ELSE
                      DIRLUM(1)=LUM(1)*COS(RANGLE)+PRS*(1.-COS(RANGLE))
     1                   *AXE(1)+SIN(RANGLE)*EPS(1)
                      DIRLUM(2)=LUM(2)*COS(RANGLE)+PRS*(1.-COS(RANGLE))
     1                   *AXE(2)+SIN(RANGLE)*EPS(2)
                      DIRLUM(3)=LUM(3)*COS(RANGLE)+PRS*(1.-COS(RANGLE))
     1                   *AXE(3)+SIN(RANGLE)*EPS(3)
                    ENDIF
                  ENDIF
                  IF(DIRLUM(3).GT.1.) DIRLUM(3)=1.
                  IF(DIRLUM(3).LT.-1.) DIRLUM(3)=-1.
                  THETALUM=ACOS(DIRLUM(3))
                  IF(I_TEST.EQ.2) THETALUM=-THETALUM
                  COEF=DIRLUM(1)+IC*DIRLUM(2)
C
                  CALL ARCSIN(COEF,DIRLUM(3),PHILUM)
C
                  ANALUM=ANADIR(1,1)*DIRLUM(1) +
     1            ANADIR(2,1)*DIRLUM(2) +
     2            ANADIR(3,1)*DIRLUM(3)
C
                  SEPSDIR_1=0.
                  SEPSDIF_1=0.
                  SEPSDIR_2=0.
                  SEPSDIF_2=0.
C
C
C    Loop over the directions of polarization for incoming photon
C
                  DO JEPS=1,NEPS_I
                    IF((JEPS.EQ.1).AND.(IPOL.GE.0)) THEN
                      DIRPOL(1,JEPS,1)=COS(THETALUM)*COS(PHILUM)
                      DIRPOL(2,JEPS,1)=COS(THETALUM)*SIN(PHILUM)
                      DIRPOL(3,JEPS,1)=-SIN(THETALUM)
                    ELSE
                      DIRPOL(1,JEPS,1)=-SIN(PHILUM)
                      DIRPOL(2,JEPS,1)=COS(PHILUM)
                      DIRPOL(3,JEPS,1)=0.
                    ENDIF
                    IF(ABS(IPOL).EQ.1) THEN
                      IF(IPRINT.GT.0) THEN
                        WRITE(IUO1,61) (DIRANA(J,1),J=1,3),
     1                         (DIRLUM(K),K=1,3),
     2                         (DIRPOL(K,1,1),K=1,3),
     3                          ANALUM
                      ENDIF
                    ELSE
                      IF((JEPS.EQ.1).AND.(IPRINT.GT.0)) THEN
                        WRITE(IUO1,63) (DIRANA(J,1),J=1,3),
     1                         (DIRLUM(K),K=1,3),ANALUM
                      ENDIF
                    ENDIF
                    IF((JEPS.EQ.1).AND.(I_EXT.EQ.-1)) PRINT 89
C
                    DO K=1,3
                      DIRQIN(K)=DIRLUM(K)
                      DIRQOU(K)=DIRANA(K,1)
C                      DIRPOLIN(K)=DIRPOL(K,1,1)
                    ENDDO
C
                    DO K=1,3
                      AA(K)=DIRANA(K,1)
                    ENDDO
                    CALL XYZ2THPH(AA,THETAOU,PHIOU)
C
C
C    Loop over the directions of polarization for outgoing photon
C
                    DO JEPS2=1,NEPS_O
                      IF((JEPS2.EQ.1).AND.(IPOL.GE.0)) THEN
                        DIRPOL(1,JEPS2,2)=COS(THETAOU)*COS(PHIOU)
                        DIRPOL(2,JEPS2,2)=COS(THETAOU)*SIN(PHIOU)
                        DIRPOL(3,JEPS2,2)=-SIN(THETAOU)
                      ELSE
                        DIRPOL(1,JEPS2,2)=-SIN(PHIOU)
                        DIRPOL(2,JEPS2,2)=COS(PHIOU)
                        DIRPOL(3,JEPS2,2)=0.
                      ENDIF
                      IF(ABS(IPOL).EQ.1) THEN
                        IF(IPRINT.GT.0) THEN
                          WRITE(IUO1,61) (DIRANA(J,1),J=1,3),
     1                         (DIRLUM(K),K=1,3),
     2                         (DIRPOL(K,1,2),K=1,3),
     3                          ANALUM
                        ENDIF
                      ELSE
                        IF((JEPS2.EQ.1).AND.(IPRINT.GT.0)) THEN
                          WRITE(IUO1,63) (DIRANA(J,1),J=1,3),
     1                         (DIRLUM(K),K=1,3),ANALUM
                        ENDIF
                      ENDIF
                      IF((JEPS2.EQ.1).AND.(I_EXT.EQ.-1)) PRINT 89
C
C
C   Calculation of the coupling matrix MLIL0
C
                      DO MI=-LI,LI
                        DO LF=LF1,LF2,ISTEP_LF
                          LR=3+LF-LI
                          LRR=5*(LR-1)
                          DELTA=DLT(JE,NTYPEM,NNL,LR)
                          RD=RHOK(JE,NTYPEM,NNL,LR,1)
C
                          DO MF=-LF,LF
                            MR=3+MF-MI
                            IF((MF.LT.MI-2).OR.(MF.GT.MI+2)) GOTO 777
                            IF((INITL.EQ.0).AND.(MF.NE.MI)) GOTO 777
                            LMR=LRR+MR
C
C
C   Storage of the coupling matrix elements MLFLI along the basis
C                       directions X,Y ET Z
C
C   These basis directions refer to the polarization if IDICHR = 0
C                 but to the light when IDICHR = 1
C
C                      JBASE = 1 : X
C                      JBASE = 2 : Y
C                      JBASE = 3 : Z
C
                            CALL COUMAT(ITL,MI,LF,MF,
     1                                  DIRQIN,DIRQOU,QIN,QOU,
     2                                  DELTA,RD,MATRIX)
C
                            DO JBASE=1,3
                              DO IO=1,2
                                MLFLI(1,MI,MR,LR,JBASE,IO)=
     1                                          MATRIX(JBASE,1,IO)
                              ENDDO
C
C       in the test, IDICHR is 0
C
                              IF(IDICHR.GE.1) THEN
                                DO IO=1,2
                                  MLFLI(2,MI,MR,LR,JBASE,IO)=
     1                                          MATRIX(JBASE,2,IO)
                                ENDDO
                              ENDIF
                            ENDDO
C
                            IF(IDICHR.EQ.0) THEN
                              IF(I_TEST.NE.1) THEN
                                DO IO=1,2
                                  MLIL0(1,MI,LMR,IO)=
     1                                   MLFLI(1,MI,MR,LR,1,IO)*
     1                                   DIRPOL(1,JEPS,IO) +
     2                                   MLFLI(1,MI,MR,LR,2,IO)*
     3                                   DIRPOL(2,JEPS,IO) +
     4                                   MLFLI(1,MI,MR,LR,3,IO)*
     5                                   DIRPOL(3,JEPS,IO)
                                ENDDO
                              ELSE
                                DO IO=1,2
                                  MLIL0(1,MI,LMR,IO)=ONEC
                                ENDDO
                              ENDIF
                            ELSEIF(IDICHR.GE.1) THEN
                              IF(I_TEST.NE.1) THEN
                                DO IO=1,2
                                  MLIL0(1,MI,LMR,IO)=
     1                                   MLFLI(1,MI,MR,LR,1,IO)*
     1                                   DIRLUM(1) +
     2                                   MLFLI(1,MI,MR,LR,2,IO)*
     3                                   DIRLUM(2) +
     4                                   MLFLI(1,MI,MR,LR,3,IO)*
     5                                   DIRLUM(3)
                                ENDDO
                                DO IO=1,2
                                  MLIL0(2,MI,LMR,IO)=
     1                                   MLFLI(2,MI,MR,LR,1,IO)*
     1                                   DIRLUM(1) +
     2                                   MLFLI(2,MI,MR,LR,2,IO)*
     3                                   DIRLUM(2) +
     4                                   MLFLI(2,MI,MR,LR,3,IO)*
     5                                   DIRLUM(3)
                                ENDDO
                              ELSE
                                DO IO=1,2
                                  MLIL0(1,MI,LMR,IO)=ONEC
                                ENDDO
                              ENDIF
                            ENDIF
  777                       CONTINUE
                          ENDDO
                        ENDDO
                      ENDDO
C
                      SRDIF_1=0.
                      SRDIR_1=0.
                      SRDIF_2=0.
                      SRDIR_2=0.
C
C
C    Loop over the different directions of the analyzer contained in a cone
C
C
                      DO JDIR=1,NDIR
                        IF(IATTS.EQ.1) THEN
                          ATTSE=EXP(-ZSURFE*GAMMA/DIRANA(3,JDIR))
                        ENDIF
C
                        SMIDIR_1=(0.0,0.0)
                        SMIDIF_1=(0.0,0.0)
                        SMIDIR_2=(0.0,0.0)
                        SMIDIF_2=(0.0,0.0)
C
C      SUM OF MI WITH DISTRIBUTION IN THOMAS SCATTERING
C
C
                        LME=LMAX(1,JE)
                        DO MI=-LI,LI
                          SJDIR_1=ZEROC
                          SJDIF_1=ZEROC
                          SJDIR_2=ZEROC
                          SJDIF_2=ZEROC
C
                          CALL THOMPSON(DIRPOL,QIN,DIRQIN,QOU,DIRQOU,3,
     5                                 SDIRTHM_1,SDIRTHM_2)
C
C
C    Calculation of the direct emission (used a a reference for the output)
C
                          DO L1=LF1,LF2,ISTEP_LF
                            IL1=L1*L1+L1+1
                            LR1=3+L1-LI
                            LRR1=5*(LR1-1)
                            IF(ISPEED.EQ.1) THEN
                              R2=TL(L1,1,1,JE)
                            ELSE
                              R2=TLT(L1,1,1,JE)
                            ENDIF
                            DO M1=-L1,L1
                              IND1=IL1+M1
                              MR1=3+M1-MI
                              LMR1=LRR1+MR1
                              INDF=ILF+MF
                              IF((M1.LT.MI-2).OR.(M1.GT.MI+2)) GOTO 444
                              IF((INITL.EQ.0).AND.(M1.NE.MI)) GOTO 444
C
                              SJDIR_1=SJDIR_1+MLIL0(1,MI,LMR1,2)*R2*
     1                            MLIL0(1,MI,LMR1,1)
                              IF(IDICHR.GE.1) THEN
                                SJDIR_2=SJDIR_2+MLIL0(2,MI,LMR1,2)*R2*
     1                              MLIL0(2,MI,LMR1,1)
                              ENDIF
C
C    Contribution of the absorber to TAU (initialization of SJDIF)
C
                              IF(I_TEST.EQ.2) GOTO 444
C
                              DO L2=LF1,LF2,ISTEP_LF
                                LR2=3+L2-LI
                                LRR2=5*(LR2-1)
                                IL2=L2*L2+L2+1
C
                                DO M2=-L2,L2
                                  MR2=3+M2-MI
                                  LMR2=LRR2+MR2
                                  IND2=IL2+M2
C
                                  IF((M2.LT.MI-2).OR.(M2.GT.MI+2))
     1                                                      GOTO 442
                                  IF((INITL.EQ.0).AND.(M2.NE.MI))
     1                                                      GOTO 442
C
                                  SJDIF_1=SJDIF_1+MLIL0(1,MI,LMR1,2)*
     1                                TAU(IND1,IND2,1)*
     2                                MLIL0(1,MI,LMR2,1)
                                  IF(IDICHR.GE.1) THEN
                                    SJDIF_2=SJDIF_2+MLIL0(2,MI,LMR1,2)*
     1                                  TAU(IND1,IND2,1)*
     2                                  MLIL0(2,MI,LMR2,1)
                                  ENDIF
 442                              CONTINUE
                                ENDDO
                              ENDDO
C
 444                          CONTINUE
                            ENDDO
                          ENDDO
                        
                          SJDIF_1=SJDIF_1*ATTSE
                          IF(IDICHR.GE.1) THEN
                            SJDIF_2=SJDIF_2*ATTSE
                          ENDIF
C
 111                      CONTINUE
C
                          SMIDIF_1=SMIDIF_1+SJDIF_1*VKR*CFM+SDIRTHM_1
                          SMIDIR_1=SMIDIR_1+SJDIR_1*VKR*CFM+SDIRTHM_1
     1                          +SDIRTHM_1
                          IF(IDICHR.GE.1) THEN
                            SMIDIF_2=SMIDIF_2+SJDIF_2*VKR*CFM+SDIRTHM_2
                            SMIDIR_2=SMIDIR_2+SJDIR_2*VKR*CFM+SDIRTHM_2
     1                          +SDIRTHM_2
                          ENDIF
C
C
C     End of the loop over MI
C
                        ENDDO
C
                        SRDIR_1=SRDIR_1+ABS(SMIDIR_1)*ABS(SMIDIR_1)
                        SRDIF_1=SRDIF_1+ABS(SMIDIF_1)*ABS(SMIDIF_1)
                        IF(IDICHR.GE.1) THEN
                          SRDIR_2=SRDIR_2+ABS(SMIDIR_2)*ABS(SMIDIR_2)
                          SRDIF_2=SRDIF_2+ABS(SMIDIF_2)*ABS(SMIDIF_2)
                        ENDIF
 220                    CONTINUE
C
C     End of the loop on the directions of the analyzer
C
                      ENDDO
C
                      SEPSDIF_1=SEPSDIF_1+SRDIF_1/NDIR
                      SEPSDIR_1=SEPSDIR_1+SRDIR_1/NDIR
                      IF(IDICHR.GE.1) THEN
                        SEPSDIF_2=SEPSDIF_2+SRDIF_2/NDIR
                        SEPSDIR_2=SEPSDIR_2+SRDIR_2/NDIR
                      ENDIF
C
C
C     End of the loop on the polarization for outgoing photon
C
                    ENDDO
C
C
C     End of the loop on the polarization for incident photon
C
                  ENDDO
C
                  SSETDIR_1=SSETDIR_1+SEPSDIR_1*W*R0E1*R0E2
                  SSETDIF_1=SSETDIF_1+SEPSDIF_1*W*R0E1*R0E2
                  IF(ICHKDIR.EQ.2) THEN
                    IF(JSET.EQ.JREF) THEN
                      SSET2DIR_1=SEPSDIR_1
                      SSET2DIF_1=SEPSDIF_1
                    ENDIF
                  ENDIF
                  IF(IDICHR.GE.1) THEN
                    SSETDIR_2=SSETDIR_2+SEPSDIR_2*W*R0E1*R0E2
                    SSETDIF_2=SSETDIF_2+SEPSDIF_2*W*R0E1*R0E2
                    IF(ICHKDIR.EQ.2) THEN
                      IF(JSET.EQ.JREF) THEN
                        SSET2DIR_2=SEPSDIR_2
                        SSET2DIF_2=SEPSDIF_2
                      ENDIF
                    ENDIF
                  ENDIF
C
C     End of the loop on the set averaging
C
                ENDDO
C
                IF(IDICHR.EQ.0) THEN
                  IF(ISOM.EQ.2) THEN
                    WRITE(IOUT,67) JPLAN,JFICH,DTHETAP,DPHIP,ECIN,
     1                       SSETDIR_1,SSETDIF_1
                    IF(ICHKDIR.EQ.2) THEN
                      WRITE(IOUT,67) JPLAN,JFICH,DTHETAP,DPHIP,ECIN,
     1                        SSET2DIR_1,SSET2DIF_1
                    ENDIF
                  ELSE
                    WRITE(IOUT,67) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,
     1                       SSETDIR_1,SSETDIF_1
                    IF(ICHKDIR.EQ.2) THEN
                      WRITE(IOUT,67) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,
     1                        SSET2DIR_1,SSET2DIF_1
                    ENDIF
                  ENDIF
                ELSE
                  IF(ISOM.EQ.2) THEN
                    WRITE(IOUT,72) JPLAN,JFICH,DTHETAP,DPHIP,ECIN,
     1                       SSETDIR_1,SSETDIF_1,
     2                       SSETDIR_2,SSETDIF_2
                    IF(ICHKDIR.EQ.2) THEN
                      WRITE(IOUT,72) JPLAN,JFICH,DTHETAP,DPHIP,ECIN,
     1                        SSET2DIR_1,SSET2DIF_1,
     2                        SSET2DIR_2,SSET2DIF_2
                    ENDIF
                  ELSE
                    WRITE(IOUT,72) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,
     1                       SSETDIR_1,SSETDIF_1,SSETDIR_2,SSETDIF_2
                    IF(ICHKDIR.EQ.2) THEN
                      WRITE(IOUT,72) JPLAN,JEMET,DTHETAP,DPHIP,ECIN,
     1                        SSET2DIR_1,SSET2DIF_1,
     2                        SSET2DIR_2,SSET2DIF_2
                    ENDIF
                  ENDIF
                ENDIF
 222            CONTINUE
C
C     End of the loop on the scanned angle
C
              ENDDO
C
   8          CONTINUE
C
C     End of the loop on the fixed angle
C
            ENDDO
C
C     End of the loop on the energy
C
            CLOSE(IUI6)
          ENDDO
C
   3      CONTINUE
C
C     End of the loop on the emitters
C
        ENDDO
C
        GO TO 1
   5    IPLAN=JPLAN-1
        IJK=IJK+1
        IF((IJK.EQ.1).AND.(IPRINT.GT.0)) THEN
          IF(I_TEST.NE.2) WRITE(IUO1,54) IPLAN
        ENDIF
   1    CONTINUE
C
C     End of the loop on the planes
C
      ENDDO
C
      IF(ABS(I_EXT).GE.1) CLOSE(IUI6)
      IF((ISOM.EQ.0).OR.(JFICH.EQ.NFICHLEC)) WRITE(IOUT,*)
      IF(SPECTRO.EQ.'APC') CLOSE(IOUT)
      IF(SPECTRO.EQ.'APC') GOTO 7
c      IF(((NEMET.GT.1).OR.(NPLAN.GT.1)).AND.(ISOM.EQ.0)) THEN
      IF(((NEMET.GT.1).OR.(NPLAN.GT.0)).AND.(ISOM.EQ.0)) THEN
        NP=0
        CALL TREAT_PHD(ISOM,NFICHLEC,JFICH,NP)
      ENDIF
      IF(I_EXT.EQ.2) THEN
        CALL WEIGHT_SUM(ISOM,I_EXT,0,1)
      ENDIF
      GOTO 7
   6  WRITE(IUO1,55)
C
   9  FORMAT(9(2X,I1),2X,I2)
  11  FORMAT(I4,2X,I4,2X,I4)
  12  FORMAT(2X,A3,11X,A13)
  13  FORMAT(6X,I1,1X,I3,2X,I4)
  14  FORMAT(6X,I1,1X,I3,3X,I3)
  19  FORMAT(2(2X,I1),1X,I2,6(2X,I1),2X,I2)
  20  FORMAT(2(5X,F6.2,2X,F6.2),2X,I1)
  21  FORMAT(10X,E12.6,3X,E12.6)
  22  FORMAT(16X,'INTERNAL CALCULATION OF MEAN SQUARE DISPLACEMENTS',/,
     1       25X,' BY DEBYE UNCORRELATED MODEL:',/)
  23  FORMAT(21X,'ATOM TYPE ',I5,'  MSD = ',F8.6,' ANG**2')
  51  FORMAT(/////,2X,'******* PLANE NUMBER ',I3,' DOES NOT CONTAIN ',
     *'ANY ABSORBER OF TYPE ',I2,' *******')
  52  FORMAT(/////,2X,'******* PLANE NUMBER ',I3,' POSITION OF ',
     1'THE ABSORBER : (',F6.3,',',F6.3,',',F6.3,') *******',/,2X,
     2'******* ',19X,'THIS ABSORBER IS OF TYPE ',I2,20X,' *******')
  53  FORMAT(//,2X,'ORDER ',I2,'  TOTAL NUMBER OF PATHS     : ',F15.1,
     1             /,10X,'  EFFECTIVE NUMBER OF PATHS : ',F15.1,
     2             /,10X,'  MINIMAL INTENSITY         : ',E12.6,
     3             2X,'No OF THE PATH : ',F15.1,
     4             /,10X,'  MAXIMAL INTENSITY         : ',E12.6,
     5             2X,'No OF THE PATH : ',F15.1)
  54  FORMAT(//,7X,'DUE TO THE SIZE OF THE CLUSTER, THE SUMMATION',
     *' HAS BEEN TRUNCATED TO THE ',I2,' TH PLANE')
  55  FORMAT(///,12X,' <<<<<<<<<<  THIS VALUE OF ILPM IS NOT',
     *'AVAILABLE  >>>>>>>>>>')
  56  FORMAT(4X,'LATTICE PARAMETER A = ',F6.3,' ANGSTROEMS',4X,
     *'MEAN FREE PATH = ',F6.3,' * A',//)
  57  FORMAT(25X,'CLUSTER RADIUS = ',F6.3,' *A')
  58  FORMAT(//,2X,'ORDER ',I2,'  TOTAL NUMBER OF PATHS     : ',I10,
     1             /,10X,'  EFFECTIVE NUMBER OF PATHS : ',I10,
     2             /,10X,'  MINIMAL INTENSITY         : ',E12.6,
     3             2X,'No OF THE PATH : ',I10,
     4             /,10X,'  MAXIMAL INTENSITY         : ',E12.6,
     5             2X,'No OF THE PATH : ',I10)
  59  FORMAT(//,15X,'THE SCATTERING DIRECTION IS GIVEN INSIDE ',
     *'THE CRYSTAL')
  60  FORMAT(7X,'THE POSITIONS OF THE ATOMS ARE GIVEN WITH RESPECT ',
     *'TO THE ABSORBER')
  61  FORMAT(///,4X,'..........  DIRECTION OF THE DETECTOR      : (',
     1      F6.3,',',F6.3,',',F6.3,
     2       ')  ..........',/,16X,'DIRECTION OF THE LIGHT    ',
     3       '     : (',F6.3,',',F6.3,',',F6.3,
     4      ')',/,16X,'DIRECTION OF THE POLARIZATION  : (',
     5      F6.3,',',F6.3,',',F6.3,')',/,16X,'ANALYZER.LIGHT       ',
     6       '          :        ',F7.4)
  63  FORMAT(///,4X,'..........  DIRECTION OF THE DETECTOR      : (',
     1       F6.3,',',F6.3,',',F6.3,
     2       ')  ..........',/,16X,'DIRECTION OF THE LIGHT    ',
     3       '     : (',F6.3,',',F6.3,',',F6.3,')',/,16X,
     4       'ANALYZER.LIGHT               :        ',F7.4)
  65  FORMAT(////,3X,'++++++++++++++++++',9X,
     *'THETA = ',F6.2,' DEGREES',9X,'++++++++',
     *'++++++++++',///)
  66  FORMAT(////,3X,'++++++++++++++++++',9X,
     *'PHI = ',F6.2,' DEGREES',9X,'++++++++++',
     *'++++++++++',///)
  67  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,
     1       2X,E12.6)
  68  FORMAT(10X,'  CUT-OFF INTENSITY       : ',E12.6)
  69  FORMAT(9X,I2,2X,E12.6,7X,E12.6,1X,F6.3,1X,10(I3,2X))
  70  FORMAT(2X,I2,2X,I10,7X,E12.6,2X,F6.3,7X,I2,7X,10(I3,2X))
  71  FORMAT(//,1X,'JDIF',4X,'No OF THE PATH',2X,
     1       'INTENSITY',3X,'LENGTH',4X,'ABSORBER',2X,
     2       'ORDER OF THE SCATTERERS',/)
  72  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,
     1       2X,E12.6,2X,E12.6,2X,E12.6)
  74  FORMAT(10X,'<===== NUMBER OF PATHS TOO LARGE FOR PRINTING ',
     1       '=====>')
  76  FORMAT(2X,I2,2X,E12.6,7X,E12.6,2X,F6.3,7X,I2,7X,10(I3,2X))
  77  FORMAT('   ')
  79  FORMAT(2X,I3,2X,I2,2X,I4,2X,I4,2X,I4)
  80  FORMAT(///)
  81  FORMAT(//,1X,'RANK',1X,'ORDER',4X,'No PATH',3X,
     1       'INTENSITY',3X,'LENGTH',4X,'ABS',3X,
     2       'ORDER OF THE SCATTERERS',/)
  82  FORMAT(I3,4X,I2,1X,E12.6,3X,E12.6,2X,F6.3,4X,I2,4X,10(I3,1X))
  83  FORMAT(I3,4X,I2,1X,I10,3X,E12.6,2X,F6.3,4X,I2,4X,10(I3,1X))
  84  FORMAT(/////,18X,'THE ',I3,' MORE INTENSE PATHS BY DECREASING',
     1      ' ORDER :',/,24X,'(THE LENGTH IS GIVEN IN UNITS ',
     2      'OF A)')
  85  FORMAT(/////,25X,' PATHS USED IN THE CALCULATION : ',
     1       /,24X,'(THE LENGTH IS GIVEN IN UNITS OF A)')
  86  FORMAT(2X,I3,1X,I4,5X,F8.3,3X,F8.3,3X,E12.6)
  87  FORMAT(2X,I2,2X,I3,2X,I2,2X,I3,2X,I3,2X,I3,2X,I1,2X,I2,2X,I2,
     1       2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
  88  FORMAT(/,19X,'TILTED THETA =',F6.2,5X,'TILTED PHI =',
     1         F6.2)
  89  FORMAT(/,4X,'..........................................',
     1         '.....................................')
C
   7  RETURN
C
      END
C
C
C======================================================================
C
      SUBROUTINE THOMPSON(DIRPOL,
     1                  QIN,DIRQIN,
     2                  QOU,DIRQOU,
     4                  DQ0,
     5                  SDIRTHM_1,SDIRTHM_2)
C
C
C       This subroutine calculates Thomas scattering contribution
C    which is (1/2)(e_in*e_out)...
C
C                         last modifed 2013.12.20
C
C
      REAL DIRPOL(3,2,2)
      REAL QIN,DIRQIN(3)
      REAL QOU,DIRQOU(3)
C
      INTEGER DQ0
      COMPLEX SDIRTHM_1,SDIRTHM_2
C
      REAL Q,Q2,DIRQ(3)
      COMPLEX F,F1
      INTEGER DQ
C
C
      F=(0.0,0.0)
      DO K=1,3
        F=F+DIRPOL(K,1,1)*DIRPOL(K,1,2)
      ENDDO
      F=0.5*F
C
      DQ=DQ0
      IF(DQ.NE.1.AND.DQ.NE.2.AND.DQ.NE.3) THEN
        WRITE(*,2000)
 2000   FORMAT(//,2X,'WARNING!, only dipole, quadrupole or octpole',
     1            'approximation can be treated!')
        IF(ABS(DQ).LT.1) THEN
          DQ=1
          WRITE(*,2100)
 2100     FORMAT(/,2X,'Change to Dipole approximation calculation')
        ELSEIF(ABS(DQ).LT.2) THEN
          DQ=2
          WRITE(*,2200)
 2200     FORMAT(/,2X,'Change to Quadrupole approximation calculation')
        ELSE
          DQ=3
          WRITE(*,2300)
 2300     FORMAT(/,2X,'Change to Octpole approximation calculation')
        ENDIF
      ENDIF
C
      IF(DQ.EQ.1.OR.DQ.EQ.2) THEN
        SDIRTHM_1=F
        SDIRTHM_2=F
      ELSE
C
        F1=(0.0,0.0)
        CALL VECTDIFD(QIN,DIRQIN,QOU,DIRQOU,Q,DIRQ)
        Q2=Q*Q
        F1=(1.0/6.0)*Q2*F1
        F=F*(1.0-F1)
        SDIRTHM_1=F
        SDIRTHM_2=SDIRTHM_1
      ENDIF
C
      RETURN
C
      END
C
C
C======================================================================
C
C
      SUBROUTINE VECTDIFD(V1,DV1,V2,DV2,V,DV)
C
C       It calculate the vector V=V1-V2
C       Not finished yet.
C
C                                   last modified 2013.12.20
C
      REAL V1,DV1(3)
      REAL V2,DV2(3)
      REAL V,DV(3)
C
      V=V1-V2
      V=0.0
      DO K=1,3
        DV(K)=DV1(K)-DV2(K)
      ENDDO
      DV(1)=0.0
      DV(2)=0.0
      DV(3)=1.0
C
      RETURN
C
      END
C
C
C=======================================================================
C
      SUBROUTINE TREAT_PHD(ISOM,NFICHLEC,JFICH,NP)
C
C  This routine sums up the calculations corresponding to different
C        absorbers or different planes when this has to be done
C             (parameter ISOM in the input data file).
C
C                                       Last modified : 24 Jan 2013
C
      INCLUDE 'spec.inc'
C
      PARAMETER(N_HEAD=5000,N_FILES=1000)
C
      CHARACTER*3 SPECTRO,STEREO,DUMMY
      CHARACTER*7 INTERACT
      CHARACTER*13 OUTDATA
      CHARACTER*72 HEAD(N_HEAD,N_FILES)
C
      REAL TAB(NDIM_M,4)
      REAL ECIN(NE_M),DTHETA(NTH_M),DPHI(NPH_M)
C
      COMMON /OUTUNITS/ IUO1,IUO2,IUO3,IUO4,IUSCR,IUSCR2
      COMMON /TYPEXP/ DUMMY,INTERACT,STEREO
      COMMON /VALIN/ PHI0,E0,THETA0,THLUM,PHILUM,ELUM,VINT,NONVOL(100)
      COMMON /VALFIN/ PHI1,EFIN,THETA1
C
      DATA JVOL,JTOT/0,-1/
C
      REWIND IUO2
C
C  Reading and storing the headers:
C
      NHEAD=0
      DO JLINE=1,N_HEAD
        READ(IUO2,888) HEAD(JLINE,JFICH)
        NHEAD=NHEAD+1
        IF(HEAD(JLINE,JFICH)(1:6).EQ.'      ') GOTO 333
      ENDDO
C
 333  CONTINUE
C
      READ(IUO2,15) SPECTRO,OUTDATA
      READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE,
     1             IPH_1,I_EXT
C
      IF(I_EXT.EQ.2) THEN
        IPH_1=0
      ENDIF
C
      IF(ISOM.EQ.0) THEN
C
C........  ISOM = 0 : case of independent input files  .................
C
        READ(IUO2,1) NPLAN,NEMET,NTHETA,NPHI,NE
C
        IF(IPH_1.EQ.1) THEN
          N_FIXED=NPHI
          FIX0=PHI0
          FIX1=PHI1
          N_SCAN=NTHETA
        ELSE
          N_FIXED=NTHETA
          FIX0=THETA0
          FIX1=THETA1
          IF(STEREO.EQ.'YES') THEN
            NPHI=INT((PHI1-PHI0)*FLOAT(NTHETA-1)/
     1               (THETA1-THETA0)+0.0001)+1
            IF(NTHETA*NPHI.GT.NPH_M) GOTO 37
          ENDIF
          N_SCAN=NPHI
        ENDIF
C
        IF(I_EXT.EQ.-1) THEN
          N_SCAN=2*N_SCAN
        ENDIF
C
        IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
          NDP=NEMET*NTHETA*NPHI*NE
        ELSEIF(I_EXT.EQ.-1) THEN
          NDP=NEMET*NTHETA*NPHI*NE*2
        ELSEIF(I_EXT.EQ.2) THEN
          NDP=NEMET*NTHETA*NE
          N_FIXED=NTHETA
          N_SCAN=NPHI
          IF((N_FIXED.GT.NTH_M).OR.(N_FIXED.GT.NPH_M)) GOTO 35
        ENDIF
C
        NTT=NPLAN*NDP
        IF(NTT.GT.NDIM_M) GOTO 5
C
        DO JPLAN=1,NPLAN
          DO JEMET=1,NEMET
            DO JE=1,NE
C
              DO J_FIXED=1,N_FIXED
                IF(N_FIXED.GT.1) THEN
                  XINCRF=FLOAT(J_FIXED-1)*
     1               (FIX1-FIX0)/FLOAT(N_FIXED-1)
                ELSEIF(N_FIXED.EQ.1) THEN
                  XINCRF=0.
                ENDIF
                IF(IPH_1.EQ.1) THEN
                  JPHI=J_FIXED
                ELSE
                  THETA=THETA0+XINCRF
                  JTHETA=J_FIXED
                  IF((ABS(THETA).GT.90.).AND.(I_EXT.NE.2)) GOTO 11
                ENDIF
                IF(STEREO.EQ.' NO') THEN
                  N_SCAN_R=N_SCAN
                ELSE
                  RTHETA=THETA*0.017453
                  FIX_STEP=(FIX1-FIX0)/FLOAT(N_FIXED-1)
                  N_SCAN_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)
     1                     +1
                ENDIF
C
                DO J_SCAN=1,N_SCAN_R
                  IF(IPH_1.EQ.1) THEN
                    JTHETA=J_SCAN
                  ELSE
                    JPHI=J_SCAN
                  ENDIF
C
                  JLIN=(JPLAN-1)*NDP +
     1            (JEMET-1)*NE*N_FIXED*N_SCAN +
     2            (JE-1)*N_FIXED*N_SCAN +
     3            (JTHETA-1)*NPHI + JPHI
C
                  IF(I_EXT.LE.0) THEN
                    IF(STEREO.EQ.' NO') THEN
                      JPHI2=JPHI
                    ELSE
                      JPHI2=(JTHETA-1)*NPHI+JPHI
                    ENDIF
                  ELSE
                    JPHI2=JTHETA
                  ENDIF
C
                  READ(IUO2,2) JPL
                  IF(JPLAN.EQ.JPL) THEN
                    BACKSPACE IUO2
                    IF(IDICHR.EQ.0) THEN
                      READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),
     1                        ECIN(JE),TAB(JLIN,1),TAB(JLIN,2)
                      IF(I_EXT.EQ.-1) THEN
                        JLIN2=NTT+JLIN
                        READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)
                      ENDIF
                    ELSE
                      READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),
     1                         ECIN(JE),TAB(JLIN,1),TAB(JLIN,2),
     2                         TAB(JLIN,3),TAB(JLIN,4)
                      IF(I_EXT.EQ.-1) THEN
                        JLIN2=NTT+JLIN
                        READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),
     1                         DPHI(JPHI2),ECIN(JE),TAB(JLIN2,1),
     2                         TAB(JLIN2,2),TAB(JLIN2,3),TAB(JLIN2,4)
                      ENDIF
                    ENDIF
                  ELSE
                    BACKSPACE IUO2
                    DO JL=JLIN,JPLAN*NDP
                      TAB(JL,1)=0.0
                      TAB(JL,2)=0.0
                      TAB(JL,3)=0.0
                      TAB(JL,4)=0.0
                    ENDDO
                    GOTO 10
                  ENDIF
                ENDDO
              ENDDO
  11          CONTINUE
            ENDDO
          ENDDO
  10      CONTINUE
        ENDDO
C
        REWIND IUO2
C
C  Skipping the NHEAD lines of headers before rewriting:
C
        DO JLINE=1,NHEAD
          READ(IUO2,888) HEAD(JLINE,JFICH)
        ENDDO
C
        WRITE(IUO2,15) SPECTRO,OUTDATA
        WRITE(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        WRITE(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM
C
        DO JE=1,NE
          DO JTHETA=1,NTHETA
            IF(STEREO.EQ.' NO') THEN
              NPHI_R=NPHI
            ELSE
              RTHETA=DTHETA(JTHETA)*0.017453
              FIX_STEP=(THETA1-THETA0)/FLOAT(NTHETA-1)
              NPHI_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
              NPHI=INT((PHI1-PHI0)/FIX_STEP+0.0001)+1
            ENDIF
            DO JPHI=1,NPHI_R
              TOTDIF_1=0.
              TOTDIR_1=0.
              VOLDIF_1=0.
              VOLDIR_1=0.
              TOTDIF_2=0.
              TOTDIR_2=0.
              VOLDIF_2=0.
              VOLDIR_2=0.
              IF(I_EXT.EQ.-1) THEN
                TOTDIF2_1=0.
                TOTDIR2_1=0.
                VOLDIF2_1=0.
                VOLDIR2_1=0.
                TOTDIF2_2=0.
                TOTDIR2_2=0.
                VOLDIF2_2=0.
                VOLDIR2_2=0.
              ENDIF
C
              DO JPLAN=1,NPLAN
C
                SF_1=0.
                SR_1=0.
                SF_2=0.
                SR_2=0.
                IF(I_EXT.EQ.-1) THEN
                  SF2_1=0.
                  SR2_1=0.
                  SF2_2=0.
                  SR2_2=0.
                ENDIF
C
                DO JEMET=1,NEMET
                  JLIN=(JPLAN-1)*NDP +
     1            (JEMET-1)*NE*NTHETA*NPHI +
     2            (JE-1)*NTHETA*NPHI +
     3            (JTHETA-1)*NPHI + JPHI
                  SF_1=SF_1+TAB(JLIN,2)
                  SR_1=SR_1+TAB(JLIN,1)
                  IF(I_EXT.EQ.-1) THEN
                    JLIN2=NTT+JLIN
                    SF2_1=SF2_1+TAB(JLIN2,2)
                    SR2_1=SR2_1+TAB(JLIN2,1)
                  ENDIF
                  IF(IDICHR.GE.1) THEN
                    SF_2=SF_2+TAB(JLIN,4)
                    SR_2=SR_2+TAB(JLIN,3)
                    IF(I_EXT.EQ.-1) THEN
                      JLIN2=NTT+JLIN
                      SF2_2=SF2_2+TAB(JLIN2,4)
                      SR2_2=SR2_2+TAB(JLIN2,3)
                    ENDIF
                  ENDIF
                ENDDO
                IF(I_EXT.LE.0) THEN
                  IF(STEREO.EQ.' NO') THEN
                    JPHI2=JPHI
                  ELSE
                    JPHI2=(JTHETA-1)*NPHI+JPHI
                  ENDIF
                ELSE
                  JPHI2=JTHETA
                ENDIF
                IF(IDICHR.EQ.0) THEN
                  WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                      ECIN(JE),SR_1,SF_1
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                        ECIN(JE),SR2_1,SF2_1
                  ENDIF
                ELSE
                  WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                       ECIN(JE),SR_1,SF_1,SR_2,SF_2
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                         ECIN(JE),SR2_1,SF2_1,SR2_2,SF2_2
                  ENDIF
                ENDIF
                IF(JPLAN.GT.NONVOL(JFICH)) THEN
                  VOLDIF_1=VOLDIF_1+SF_1
                  VOLDIR_1=VOLDIR_1+SR_1
                  IF(I_EXT.EQ.-1) THEN
                    VOLDIF2_1=VOLDIF2_1+SF2_1
                    VOLDIR2_1=VOLDIR2_1+SR2_1
                  ENDIF
                  IF(IDICHR.GE.1) THEN
                    VOLDIF_2=VOLDIF_2+SF_2
                    VOLDIR_2=VOLDIR_1+SR_2
                    IF(I_EXT.EQ.-1) THEN
                      VOLDIF2_2=VOLDIF2_2+SF2_2
                      VOLDIR2_2=VOLDIR2_1+SR2_2
                    ENDIF
                  ENDIF
                ENDIF
                TOTDIF_1=TOTDIF_1+SF_1
                TOTDIR_1=TOTDIR_1+SR_1
                IF(I_EXT.EQ.-1) THEN
                  TOTDIF2_1=TOTDIF2_1+SF2_1
                  TOTDIR2_1=TOTDIR2_1+SR2_1
                ENDIF
                IF(IDICHR.GE.1) THEN
                  TOTDIF_2=TOTDIF_2+SF_2
                  TOTDIR_2=TOTDIR_2+SR_2
                  IF(I_EXT.EQ.-1) THEN
                    TOTDIF2_2=TOTDIF2_2+SF2_2
                    TOTDIR2_2=TOTDIR2_2+SR2_2
                  ENDIF
                ENDIF
              ENDDO
              IF(IDICHR.EQ.0) THEN
                WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),
     1                     VOLDIR_1,VOLDIF_1
                IF(I_EXT.EQ.-1) THEN
                  WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),
     1                       ECIN(JE),VOLDIR2_1,VOLDIF2_1
                ENDIF
                WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),
     1                     TOTDIR_1,TOTDIF_1
                IF(I_EXT.EQ.-1) THEN
                  WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),
     1                       ECIN(JE),TOTDIR2_1,TOTDIF2_1
                ENDIF
              ELSE
                WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),
     1                      VOLDIR_1,VOLDIF_1,VOLDIR_2,VOLDIF_2
                IF(I_EXT.EQ.-1) THEN
                  WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),
     1                  ECIN(JE),VOLDIR2_1,VOLDIF2_1,VOLDIR2_2,VOLDIF2_2
                ENDIF
                WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),ECIN(JE),
     1                      TOTDIR_1,TOTDIF_1,TOTDIR_2,TOTDIF_2
                IF(I_EXT.EQ.-1) THEN
                  WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),
     1                  ECIN(JE),TOTDIR2_1,TOTDIF2_1,TOTDIR2_2,TOTDIF2_2
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C
      ELSE
C
C........ ISOM not= 0 : multiple input files to be summed up  ..........
C
        READ(IUO2,7) NTHETA,NPHI,NE
C
        IF(IPH_1.EQ.1) THEN
          N_FIXED=NPHI
          FIX0=PHI0
          FIX1=PHI1
          N_SCAN=NTHETA
        ELSE
          N_FIXED=NTHETA
          FIX0=THETA0
          FIX1=THETA1
          IF(STEREO.EQ.'YES') THEN
            NPHI=INT((PHI1-PHI0)*FLOAT(NTHETA-1)/
     1               (THETA1-THETA0)+0.0001)+1
            IF(NTHETA*NPHI.GT.NPH_M) GOTO 37
          ENDIF
          N_SCAN=NPHI
        ENDIF
C
        IF(I_EXT.EQ.-1) THEN
          N_SCAN=2*N_SCAN
        ENDIF
C
        IF((I_EXT.EQ.0).OR.(I_EXT.EQ.1)) THEN
          NDP=NTHETA*NPHI*NE
        ELSEIF(I_EXT.EQ.-1) THEN
          NDP=NTHETA*NPHI*NE*2
        ELSEIF(I_EXT.EQ.2) THEN
          NDP=NTHETA*NE
          N_FIXED=NTHETA
          N_SCAN=NPHI
          IF((N_FIXED.GT.NTH_M).OR.(N_FIXED.GT.NPH_M)) GOTO 35
        ENDIF
C
        NTT=NFICHLEC*NDP
        IF(NTT.GT.NDIM_M) GOTO 5
C
        IF(ISOM.EQ.1) THEN
          NPLAN=NP
          NF=NP
        ELSEIF(ISOM.EQ.2) THEN
          NEMET=NFICHLEC
          NF=NFICHLEC
          NPLAN=1
        ENDIF
C
        DO JF=1,NF
C
C  Reading the headers for each file:
C
          IF(JF.GT.1) THEN
            DO JLINE=1,NHEAD
              READ(IUO2,888) HEAD(JLINE,JF)
            ENDDO
          ENDIF
C
          DO JE=1,NE
C
            DO J_FIXED=1,N_FIXED
              IF(N_FIXED.GT.1) THEN
                XINCRF=FLOAT(J_FIXED-1)*
     1               (FIX1-FIX0)/FLOAT(N_FIXED-1)
              ELSEIF(N_FIXED.EQ.1) THEN
                XINCRF=0.
              ENDIF
              IF(IPH_1.EQ.1) THEN
                JPHI=J_FIXED
              ELSE
                THETA=THETA0+XINCRF
                JTHETA=J_FIXED
                IF((ABS(THETA).GT.90.).AND.(I_EXT.NE.2)) GOTO 12
              ENDIF
              IF(STEREO.EQ.' NO') THEN
                N_SCAN_R=N_SCAN
              ELSE
                RTHETA=THETA*0.017453
                FIX_STEP=(FIX1-FIX0)/FLOAT(N_FIXED-1)
                N_SCAN_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
              ENDIF
C
              DO J_SCAN=1,N_SCAN_R
                IF(IPH_1.EQ.1) THEN
                  JTHETA=J_SCAN
                ELSE
                  JPHI=J_SCAN
                ENDIF
C
                JLIN=(JF-1)*NDP + (JE-1)*N_FIXED*N_SCAN +
     1            (JTHETA-1)*NPHI + JPHI

                IF(I_EXT.LE.0) THEN
                  IF(STEREO.EQ.' NO') THEN
                    JPHI2=JPHI
                  ELSE
                    JPHI2=(JTHETA-1)*NPHI+JPHI
                  ENDIF
                ELSE
                  JPHI2=JTHETA
                ENDIF
C
                IF(ISOM.EQ.1) THEN
                  READ(IUO2,2) JPL
                  IF(JF.EQ.JPL) THEN
                    BACKSPACE IUO2
                    IF(IDICHR.EQ.0) THEN
                      READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),
     1                          ECIN(JE),TAB(JLIN,1),TAB(JLIN,2)
                      IF(I_EXT.EQ.-1) THEN
                        JLIN2=NTT+JLIN
                        READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)
                      ENDIF
                    ELSE
                      READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),
     1                           ECIN(JE),TAB(JLIN,1),TAB(JLIN,2),
     2                           TAB(JLIN,3),TAB(JLIN,4)
                      IF(I_EXT.EQ.-1) THEN
                        JLIN2=NTT+JLIN
                        READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),
     1                             DPHI(JPHI2),ECIN(JE),
     2                             TAB(JLIN2,1),TAB(JLIN2,2),
     3                             TAB(JLIN2,3),TAB(JLIN2,4)
                      ENDIF
                    ENDIF
                  ELSE
                    BACKSPACE IUO2
                    DO JLINE=1,NHEAD
                      BACKSPACE IUO2
                    ENDDO
                    DO JL=JLIN,JF*NDP
                      TAB(JL,1)=0.0
                      TAB(JL,2)=0.0
                      TAB(JL,3)=0.0
                      TAB(JL,4)=0.0
                    ENDDO
                    GOTO 13
                  ENDIF
                ELSEIF(ISOM.EQ.2) THEN
                  IF(IDICHR.EQ.0) THEN
                    READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),
     1                        ECIN(JE),TAB(JLIN,1),TAB(JLIN,2)
                    IF(I_EXT.EQ.-1) THEN
                      JLIN2=NTT+JLIN
                      READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)
                    ENDIF
                  ELSE
                    READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),
     1                         ECIN(JE),TAB(JLIN,1),TAB(JLIN,2),
     2                         TAB(JLIN,3),TAB(JLIN,4)
                    IF(I_EXT.EQ.-1) THEN
                      JLIN2=NTT+JLIN
                      READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),
     1                           DPHI(JPHI2),ECIN(JE),
     2                           TAB(JLIN2,1),TAB(JLIN2,2),
     3                           TAB(JLIN2,3),TAB(JLIN2,4)
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
 12           CONTINUE
            ENDDO
          ENDDO
 13       CONTINUE
        ENDDO
C
        REWIND IUO2
C
C  Writing the headers:
C
        DO JLINE=1,2
          WRITE(IUO2,888) HEAD(JLINE,1)
        ENDDO
        DO JF=1,NFICHLEC
          DO JLINE=3,6
            WRITE(IUO2,888) HEAD(JLINE,JF)
          ENDDO
          WRITE(IUO2,888) HEAD(2,JF)
        ENDDO
        DO JLINE=7,NHEAD
          WRITE(IUO2,888) HEAD(JLINE,1)
        ENDDO
C
        WRITE(IUO2,15) SPECTRO,OUTDATA
        WRITE(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        WRITE(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM
C
        IF(ISOM.EQ.1) THEN
C
          DO JE=1,NE
C
            DO JTHETA=1,NTHETA
              IF(STEREO.EQ.' NO') THEN
                NPHI_R=NPHI
              ELSE
                RTHETA=DTHETA(JTHETA)*0.017453
                FIX_STEP=(THETA1-THETA0)/FLOAT(NTHETA-1)
                NPHI_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
                NPHI=INT((PHI1-PHI0)/FIX_STEP+0.0001)+1
              ENDIF
              DO JPHI=1,NPHI_R
C
                TOTDIF_1=0.
                TOTDIR_1=0.
                VOLDIF_1=0.
                VOLDIR_1=0.
                TOTDIF_2=0.
                TOTDIR_2=0.
                VOLDIF_2=0.
                VOLDIR_2=0.
                IF(I_EXT.EQ.-1) THEN
                  TOTDIF2_1=0.
                  TOTDIR2_1=0.
                  VOLDIF2_1=0.
                  VOLDIR2_1=0.
                  TOTDIF2_2=0.
                  TOTDIR2_2=0.
                  VOLDIF2_2=0.
                  VOLDIR2_2=0.
                ENDIF
C
                DO JPLAN=1,NPLAN
                  JF=JPLAN
C
                  JLIN=(JF-1)*NDP + (JE-1)*NTHETA*NPHI +
     1             (JTHETA-1)*NPHI + JPHI
C
                  SR_1=TAB(JLIN,1)
                  SF_1=TAB(JLIN,2)
                  IF(I_EXT.EQ.-1) THEN
                    JLIN2=NTT+JLIN
                    SF2_1=TAB(JLIN2,2)
                    SR2_1=TAB(JLIN2,1)
                  ENDIF
                  IF(I_EXT.LE.0) THEN
                    IF(STEREO.EQ.' NO') THEN
                      JPHI2=JPHI
                    ELSE
                      JPHI2=(JTHETA-1)*NPHI+JPHI
                    ENDIF
                  ELSE
                    JPHI2=JTHETA
                  ENDIF
                  IF(IDICHR.EQ.0) THEN
                    WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                        ECIN(JE),SR_1,SF_1
                    IF(I_EXT.EQ.-1) THEN
                      WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                          ECIN(JE),SR2_1,SF2_1
                    ENDIF
                  ELSE
                    SR_2=TAB(JLIN,3)
                    SF_2=TAB(JLIN,4)
                    IF(I_EXT.EQ.-1) THEN
                      JLIN2=NTT+JLIN
                      SF2_2=TAB(JLIN2,4)
                      SR2_2=TAB(JLIN2,3)
                    ENDIF
                    WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                         ECIN(JE),SR_1,SF_1,SR_2,SF_2
                    IF(I_EXT.EQ.-1) THEN
                      WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                           ECIN(JE),SR2_1,SF2_1,SR2_2,SF2_2
                    ENDIF
                  ENDIF
                  IF(NONVOL(JPLAN).EQ.0) THEN
                    VOLDIF_1=VOLDIF_1+SF_1
                    VOLDIR_1=VOLDIR_1+SR_1
                    IF(I_EXT.EQ.-1) THEN
                      VOLDIF2_1=VOLDIF2_1+SF2_1
                      VOLDIR2_1=VOLDIR2_1+SR2_1
                    ENDIF
                    IF(IDICHR.GE.1) THEN
                      VOLDIF_2=VOLDIF_2+SF_2
                      VOLDIR_2=VOLDIR_2+SR_2
                      IF(I_EXT.EQ.-1) THEN
                        VOLDIF2_2=VOLDIF2_2+SF2_2
                        VOLDIR2_2=VOLDIR2_1+SR2_2
                      ENDIF
                    ENDIF
                  ENDIF
                  TOTDIF_1=TOTDIF_1+SF_1
                  TOTDIR_1=TOTDIR_1+SR_1
                  IF(I_EXT.EQ.-1) THEN
                    TOTDIF2_1=TOTDIF2_1+SF2_1
                    TOTDIR2_1=TOTDIR2_1+SR2_1
                  ENDIF
                  IF(IDICHR.GE.1) THEN
                    TOTDIF_2=TOTDIF_2+SF_2
                    TOTDIR_2=TOTDIR_2+SR_2
                    IF(I_EXT.EQ.-1) THEN
                      TOTDIF2_2=TOTDIF2_2+SF2_2
                      TOTDIR2_2=TOTDIR2_2+SR2_2
                    ENDIF
                  ENDIF
                ENDDO
C
                IF(IDICHR.EQ.0) THEN
                  WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),
     1                       ECIN(JE),VOLDIR_1,VOLDIF_1
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),
     1                         ECIN(JE),VOLDIR2_1,VOLDIF2_1
                  ENDIF
                  WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),
     1                       ECIN(JE),TOTDIR_1,TOTDIF_1
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),
     1                         ECIN(JE),TOTDIR2_1,TOTDIF2_1
                  ENDIF
                ELSE
                  WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),
     1                    ECIN(JE),VOLDIR_1,VOLDIF_1,VOLDIR_2,VOLDIF_2
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),
     1                          ECIN(JE),VOLDIR2_1,VOLDIF2_1,
     3                          VOLDIR2_2,VOLDIF2_2
                  ENDIF
                  WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),
     1                    ECIN(JE),TOTDIR_1,TOTDIF_1,TOTDIR_2,TOTDIF_2
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),
     1                          ECIN(JE),TOTDIR2_1,TOTDIF2_1,
     3                          TOTDIR2_2,TOTDIF2_2
                  ENDIF
                ENDIF
C
              ENDDO
            ENDDO
          ENDDO
        ELSEIF(ISOM.EQ.2) THEN
          DO JE=1,NE
C
            DO JTHETA=1,NTHETA
              IF(STEREO.EQ.' NO') THEN
                NPHI_R=NPHI
              ELSE
                RTHETA=DTHETA(JTHETA)*0.017453
                FIX_STEP=(THETA1-THETA0)/FLOAT(NTHETA-1)
                NPHI_R=INT((PHI1-PHI0)*SIN(RTHETA)/FIX_STEP+0.0001)+1
                NPHI=INT((PHI1-PHI0)/FIX_STEP+0.0001)+1
              ENDIF
              DO JPHI=1,NPHI_R
C
                SF_1=0.
                SR_1=0.
                SF_2=0.
                SR_2=0.
                IF(I_EXT.EQ.-1) THEN
                  SF2_1=0.
                  SR2_1=0.
                  SF2_2=0.
                  SR2_2=0.
                ENDIF
C
                DO JEMET=1,NEMET
                  JF=JEMET
C
                  JLIN=(JF-1)*NDP + (JE-1)*NTHETA*NPHI +
     1              (JTHETA-1)*NPHI + JPHI
C
                  SF_1=SF_1+TAB(JLIN,2)
                  SR_1=SR_1+TAB(JLIN,1)
                  IF(I_EXT.EQ.-1) THEN
                    JLIN2=NTT+JLIN
                    SF2_1=SF2_1+TAB(JLIN2,2)
                    SR2_1=SR2_1+TAB(JLIN2,1)
                  ENDIF
                  IF(IDICHR.GE.1) THEN
                    SF_2=SF_2+TAB(JLIN,4)
                    SR_2=SR_2+TAB(JLIN,3)
                    IF(I_EXT.EQ.-1) THEN
                      JLIN2=NTT+JLIN
                      SF2_2=SF2_2+TAB(JLIN2,4)
                      SR2_2=SR2_2+TAB(JLIN2,3)
                    ENDIF
                  ENDIF
                ENDDO
                IF(I_EXT.LE.0) THEN
                  IF(STEREO.EQ.' NO') THEN
                    JPHI2=JPHI
                  ELSE
                    JPHI2=(JTHETA-1)*NPHI+JPHI
                  ENDIF
                ELSE
                  JPHI2=JTHETA
                ENDIF
                IF(IDICHR.EQ.0) THEN
                  WRITE(IUO2,3) JPL,DTHETA(JTHETA),DPHI(JPHI2),
     1                        ECIN(JE),SR_1,SF_1
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                          ECIN(JE),SR2_1,SF2_1
                  ENDIF
                ELSE
                  WRITE(IUO2,23) JPL,DTHETA(JTHETA),DPHI(JPHI2),
     1                        ECIN(JE),SR_1,SF_1,SR_2,SF_2
                  IF(I_EXT.EQ.-1) THEN
                    WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),
     1                           ECIN(JE),SR2_1,SF2_1,SR2_2,SF2_2
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
      GOTO 6
C
   5  WRITE(IUO1,4)
      STOP
  35  WRITE(IUO1,36) N_FIXED
      STOP
  37  WRITE(IUO1,38) NTHETA*NPHI
      STOP
C
   1  FORMAT(2X,I3,2X,I2,2X,I4,2X,I4,2X,I4)
   2  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,
     1       2X,E12.6)
   3  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
   4  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL ',
     1       'IN THE TREAT_PHD SUBROUTINE - INCREASE NDIM_M ',
     2       '>>>>>>>>>>')
   7  FORMAT(I4,2X,I4,2X,I4)
   8  FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
   9  FORMAT(9(2X,I1),2X,I2)
  15  FORMAT(2X,A3,11X,A13)
  22  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,
     1       2X,E12.6,2X,E12.6,2X,E12.6)
  23  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,
     1       2X,E12.6,2X,E12.6)
  25  FORMAT(37X,E12.6,2X,E12.6)
  36  FORMAT(//,4X,'<<<<<<<<<<  DIMENSION OF NTH_M OR NPH_M TOO SMALL ',
     1       'IN THE INCLUDE FILE >>>>>>>>>>',/,4X,
     2        '<<<<<<<<<<                 SHOULD BE AT LEAST ',I6,
     3        '                 >>>>>>>>>>')
  38  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF NPH_M TOO SMALL ',
     1       'IN THE INCLUDE FILE >>>>>>>>>>',/,8X,
     2        '<<<<<<<<<<             SHOULD BE AT LEAST ',I6,
     3        '             >>>>>>>>>>')
 888  FORMAT(A72)
C
   6  RETURN
C
      END


