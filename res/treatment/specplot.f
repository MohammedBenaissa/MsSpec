C
C
C  ************************************************************************
C  *                                                                      *
C  *        TREATMENT OF THE ELECTRON SPECTROSCOPY CALCULATED DATA        *
C  *                                                                      *
C  ************************************************************************
C
C
C    This program reads the result file of the PhD/LEED/XAS/AED/APECS code
C      and selects in this file what is going to be plotted by the user
C      The resulting file is a two column file corresponding to the
C      selected scan. Note that this code cannot accomodate three column
C      output files yet.
C
C    Author : D. Sebilleau
C
C                                            Last modified :  2 Jun 2021
C
C
      PROGRAM MAIN
C
      PARAMETER (LINE_M=300000,N_HEAD=100)
C
      CHARACTER*1 REP1,REP2,REP3,REP4,REP5,REP6,REP7,EP
      CHARACTER*3 SPECTRO,DUMMY
      CHARACTER*6 HEAD
      CHARACTER*24 INFILE,OUTFILE
C
      COMMON /TYPE_EXP/ IPHIR,ITHETAR,IER
      COMMON /NUMBERS/ NTHETA,NPHI,NTHETA_A,NPHI_A
C
      DATA ZERO /1.E-35/
C
      I_JUMP=0
      I_AUTO=0
      NS_E=0
C
      WRITE(6,114)
C
 14   WRITE(6,*) '  '
      WRITE(6,*) '  '
      WRITE(6,*) 'NAME OF THE INPUT FILE : '
      WRITE(6,*) '  '
      READ(5,10) INFILE
C
      OPEN(UNIT=1,FILE=INFILE,STATUS='OLD')
C
C  Skipping the headers: NHEAD2 number of lines before
C          the first angular/energy record
C
      NHEAD=0
      DO JLINE=1,N_HEAD
        READ(1,888) HEAD
        NHEAD=NHEAD+1
        IF(HEAD.EQ.'      ') GOTO 333
      ENDDO
C
 333  NHEAD2=NHEAD+3
C
C  Checking if all header lines have been read
C
      READ(1,2) EP
      IF(EP.EQ.'!') THEN
        PRINT 6
        STOP
      ELSE
        BACKSPACE 1
      ENDIF
C
      READ(1,16) SPECTRO
C
      IF(SPECTRO(3:3).EQ.'D') THEN
        READ(1,33) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        READ(1,44) NPHI,NTHETA,NE,NPLAN,ISOM
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        READ(1,33) ISPIN,IDICHR,I_SO,ISFLIP
        READ(1,45) NE,NPLAN,ISOM
        ICHKDIR=0
        IPHI=0
        ITHETA=0
        IE=1
        THETA=1.0
        PHI=1.0
        IPHIR=IPHI
        ITHETAR=ITHETA
        IER=IE
      ELSEIF(SPECTRO.EQ.'APC') THEN
        READ(1,33) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        READ(1,33) ISPIN_A,IDICHR_A,I_SO_A,ISFLIP_A,ICHKDIR_A,
     1             IPHI_A,ITHETA_A,IE_A
        READ(1,44) NPHI,NTHETA,NE,NPLAN,ISOM
        READ(1,44) NPHI_A,NTHETA_A
        NE_A=1
      ENDIF
C
      ITEST=0
      IS_SP=0
      IS_DI=0
      IS_SO=0
      IF((ISPIN.EQ.1).AND.(IDICHR.EQ.1)) THEN
        I_SUM=1
      ELSE
        I_SUM=0
      ENDIF
      I_PO=0
 13   NSP=0
C
      IF(ISOM.LT.2) THEN
        IF(SPECTRO(3:3).EQ.'D') THEN
          NLINE=NPLAN+2
        ELSEIF(SPECTRO.EQ.'XAS') THEN
          NLINE=NPLAN+1
        ELSEIF(SPECTRO.EQ.'APC') THEN
          NLINE=NPLAN+2
        ENDIF
      ELSE
        NLINE=1
      ENDIF
      IF(ISPIN.EQ.1) NLINE=NLINE*2
      IF(I_SO.EQ.2) NLINE=NLINE*2
      IF(ICHKDIR.EQ.2) NLINE=NLINE*2
      IF(SPECTRO(3:3).EQ.'D')THEN
        LINE=NE*NTHETA*NPHI*NLINE
      ELSEIF(SPECTRO.EQ.'XAS') THEN
        LINE=NE*NLINE
      ELSEIF(SPECTRO.EQ.'APC') THEN
        LINE=NE*NTHETA*NPHI*NTHETA_A*NPHI_A*NLINE
      ENDIF
      IF(LINE.GT.LINE_M) THEN
        WRITE(6,40) LINE
        STOP
      ENDIF
      IF(ICHKDIR.EQ.1) THEN
        IPHIR=IPHI
        ITHETAR=ITHETA
        IER=IE
        IF(IPHI.EQ.1) THEN
          LINE2=(NPHI*IPHI+NTHETA*ITHETA)*NE*NLINE
        ELSEIF(ITHETA.EQ.1) THEN
          LINE2=NPLAN+2
        ENDIF
        GOTO 23
      ENDIF
      IF(SPECTRO.EQ.'XAS') GOTO 23
      IF(SPECTRO.EQ.'APC') GOTO 22
C
C    PhD/LEED/AED case :
C
      IF(ISPIN.EQ.0) THEN
        IF(I_SO.EQ.0) THEN
          IF(IDICHR.EQ.0) THEN
            READ(1,7) JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1
          ELSE
            READ(1,17) JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1,SR1_2,SF1_2
          ENDIF
        ELSE
          IF(IDICHR.EQ.0) THEN
            READ(1,27) JSO,JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1
          ELSE
            READ(1,37) JSO,JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1,
     2                 SR1_2,SF1_2
          ENDIF
        ENDIF
      ELSE
        IF(IDICHR.EQ.0) THEN
          READ(1,3) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,SR1_1,SF1_1
        ELSE
          READ(1,43) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,SR1_1,SF1_1,
     2               SR1_2,SF1_2
        ENDIF
      ENDIF
C
      E_C=ECIN
      THETA=DTHETA
      PHI=DPHI
      IPHIR=IPHI
      ITHETAR=ITHETA
      IER=IE
      BACKSPACE 1
C
C   Check of the different E,THETA,PHI present in the result file
C            when those variables are not varied
C                   (PhD, LEED and AED case)
C
      CALL FIXED_VAR(IE,NE,ITHETA,NTHETA,IPHI,NPHI,NS_E,
     1               ISPIN,I_SO,NLINE,IEC,ITH,IPH,NHEAD2)
C
      I_PO=IEC+ITH+IPH
      IF(I_PO.EQ.0) GOTO 23
C
 19   CONTINUE
      CALL SELEC(NS_E,NE,NTHETA,NPHI,IEC,ITH,IPH,E_C,THETA,PHI)
      GOTO 23
C
C    APECS case :  P ---> photoelectron, A ---> Auger electron
C
  22  WRITE(6,*) '  '
      WRITE(6,*) 'WHICH ELECTRON DO YOU WISH TO SCAN ? '
      WRITE(6,*) '  '
      WRITE(6,*) '     *  1  : PHOTOELECTRON'
      WRITE(6,*) '     *  2  : AUGER ELECTRON'
      WRITE(6,*) '  '
      READ(5,11) NS_E
      IF((NS_E.LT.1).OR.(NS_E.GT.2)) THEN
        WRITE(6,*) '     <<<<<  WRONG CHOICE  >>>>>'
        GOTO 22
      ENDIF
      IF(NS_E.EQ.1) THEN
        IEC=ABS(IE-1)
        ITH=ABS(ITHETA-1)
        IPH=ABS(IPHI-1)
        IF((NTHETA_A.GT.1).OR.(NPHI_A.GT.1)) THEN
          WRITE(6,*) '  '
          WRITE(6,31) NTHETA_A*NPHI_A
C          IF((ICHKDIR.EQ.0).AND.(ICHKDIR_A.EQ.0)) THEN
C            IF(NTHETA_A*NPHI_A.LE.91) I_AUTO=1
C          ENDIF
          IF(I_AUTO.EQ.1) THEN
            WRITE(6,*) 'DO YOU WANT TO : '
            WRITE(6,*) '  '
            WRITE(6,*) '     *  1  : PLOT ALL OF THEM'
            WRITE(6,*) '     *  2  : PLOT SELECTED ONES'
            WRITE(6,*) '  '
            READ(5,11) N_PL
          ELSE
            N_PL=2
          ENDIF
          IF(N_PL.EQ.2) THEN
C
C  Selecting one direction for the fixed detector
C
            WRITE(6,*) '  '
            WRITE(6,*) '   --------------------------------------------'
            WRITE(6,*) '                 FIXED DETECTOR : '
            CALL FIXED_VAR(0,NE_A,0,NTHETA_A,0,NPHI_A,
     1                     NS_E,ISPIN_A,I_SO_A,NLINE,I_EC,I_TH,I_PH,
     2                     NHEAD2)
            CALL SELEC(NS_E,NE_A,NTHETA_A,NPHI_A,I_EC,I_TH,I_PH,
     1                 E_C,THETAF,PHIF)
            IF(NS_E.EQ.1) THEN
              WRITE(6,35) THETAF,PHIF
            ELSEIF(NS_E.EQ.2) THEN
              WRITE(6,34) E_C,THETAF,PHIF
            ENDIF
            WRITE(6,*) '   --------------------------------------------'
            WRITE(6,*) '  '
C
C  Selecting the energy or the complementary angle for
C     the scanned detector
C
            WRITE(6,*) '  '
            WRITE(6,*) '   --------------------------------------------'
            WRITE(6,*) '               SCANNED DETECTOR : '
            CALL FIXED_VAR(IE,NE,ITHETA,NTHETA,IPHI,NPHI,
     1                     2,ISPIN,I_SO,NLINE,IEC,ITH,IPH,NHEAD2)
            CALL SELEC(2,NE,NTHETA,NPHI,IEC,ITH,IPH,
     1                 E_C,THETA,PHI)
            WRITE(6,*) '   --------------------------------------------'
            WRITE(6,*) '  '
          ENDIF
        ELSE
          IF(ISPIN.EQ.0) THEN
            IF(I_SO.EQ.0) THEN
              IF(IDICHR.EQ.0) THEN
                READ(1,8) JPLAN,DTHETA,DPHI,ECIN,DTHETAA,DPHIA,
     1                    SR1_1,SF1_1
              ELSE
                READ(1,18) JPLAN,DTHETA,DPHI,ECIN,DTHETAA,DPHIA,
     1                     SR1_1,SF1_1,SR1_2,SF1_2
              ENDIF
            ELSE
              IF(IDICHR.EQ.0) THEN
                READ(1,28) JSO,JPLAN,DTHETA,DPHI,ECIN,DTHETAA,
     1                     DPHIA,SR1_1,SF1_1
              ELSE
                READ(1,38) JSO,JPLAN,DTHETA,DPHI,ECIN,DTHETAA,
     1                     DPHIA,SR1_1,SF1_1,SR1_2,SF1_2
              ENDIF
            ENDIF
          ELSE
            IF(IDICHR.EQ.0) THEN
              READ(1,1) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,DTHETAA,
     1                  DPHIA,SR1_1,SF1_1
            ELSE
              READ(1,41) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,DTHETAA,
     1                   DPHIA,SR1_1,SF1_1,SR1_2,SF1_2
            ENDIF
          ENDIF
C
          E_C=ECIN
          THETA=DTHETA
          PHI=DPHI
          THETAF=DTHETAA
          PHIF=DPHIA
          IPHIR=IPHI
          ITHETAR=ITHETA
          IER=IE
          BACKSPACE 1
        ENDIF
      ELSEIF(NS_E.EQ.2) THEN
        IEC=ABS(IE_A-1)
        ITH=ABS(ITHETA_A-1)
        IPH=ABS(IPHI_A-1)
        IF((NTHETA.GT.1).OR.(NPHI.GT.1)) THEN
          WRITE(6,*) '  '
          WRITE(6,32) NTHETA*NPHI
          IF((ICHKDIR.EQ.0).AND.(ICHKDIR_A.EQ.0)) THEN
            IF(NTHETA*NPHI.LE.91) I_AUTO=1
          ENDIF
          IF(I_AUTO.EQ.1) THEN
            WRITE(6,*) 'DO YOU WANT TO : '
            WRITE(6,*) '  '
            WRITE(6,*) '     *  1  : PLOT ALL OF THEM'
            WRITE(6,*) '     *  2  : PLOT SELECTED ONES'
            WRITE(6,*) '  '
            READ(5,11) N_PL
          ELSE
            N_PL=2
          ENDIF
          IF(N_PL.EQ.2) THEN
C
C  Selecting one direction for the fixed detector
C
            WRITE(6,*) '  '
            WRITE(6,*) '   --------------------------------------------'
            WRITE(6,*) '                 FIXED DETECTOR : '
            CALL FIXED_VAR(0,NE,0,NTHETA,0,NPHI,NS_E,
     1                     ISPIN,I_SO,NLINE,I_EC,I_TH,I_PH,NHEAD2)
            CALL SELEC(NS_E,NE,NTHETA,NPHI,I_EC,I_TH,I_PH,
     1                 E_C,THETAF,PHIF)
            WRITE(6,34) E_C,THETAF,PHIF
            WRITE(6,*) '   --------------------------------------------'
            WRITE(6,*) '  '
C
C  Selecting the energy or the complementary angle for
C     the scanned detector
C
            WRITE(6,*) '  '
            WRITE(6,*) '   --------------------------------------------'
            WRITE(6,*) '               SCANNED DETECTOR : '
            CALL FIXED_VAR(IE_A,NE_A,ITHETA_A,NTHETA_A,IPHI_A,
     1                     NPHI_A,1,ISPIN_A,I_SO_A,NLINE,IEC,ITH,IPH,
     2                     NHEAD2)
            CALL SELEC(1,NE_A,NTHETA_A,NPHI_A,IEC,ITH,IPH,
     1                 E_C,THETA,PHI)
            WRITE(6,*) '   --------------------------------------------'
            WRITE(6,*) '  '
          ENDIF
        ELSE
          IF(ISPIN.EQ.0) THEN
            IF(I_SO.EQ.0) THEN
              IF(IDICHR.EQ.0) THEN
                READ(1,8) JPLAN,DTHETA,DPHI,ECIN,DTHETAA,DPHIA,
     1                    SR1_1,SF1_1
              ELSE
                READ(1,18) JPLAN,DTHETA,DPHI,ECIN,DTHETAA,DPHIA,
     1                     SR1_1,SF1_1,SR1_2,SF1_2
              ENDIF
            ELSE
              IF(IDICHR.EQ.0) THEN
                READ(1,28) JSO,JPLAN,DTHETA,DPHI,ECIN,DTHETAA,
     1                     DPHIA,SR1_1,SF1_1
              ELSE
                READ(1,38) JSO,JPLAN,DTHETA,DPHI,ECIN,DTHETAA,
     1                     DPHIA,SR1_1,SF1_1,SR1_2,SF1_2
              ENDIF
            ENDIF
          ELSE
            IF(IDICHR.EQ.0) THEN
              READ(1,1) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,DTHETAA,
     1                  DPHIA,SR1_1,SF1_1
            ELSE
              READ(1,41) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,DTHETAA,
     1                   DPHIA,SR1_1,SF1_1,SR1_2,SF1_2
            ENDIF
          ENDIF
C
          E_C=ECIN
          THETA=DTHETAA
          PHI=DPHIA
          THETAF=DTHETA
          PHIF=DPHI
          IPHIR=IPHI_A
          ITHETAR=ITHETA_A
          IER=IE_A
          BACKSPACE 1
        ENDIF
      ENDIF
C
 23   IF(ICHKDIR.EQ.0) THEN
        WRITE(6,*) '  '
        WRITE(6,*) 'NAME OF THE OUTPUT FILE : '
        WRITE(6,*) '  '
        READ(5,10) OUTFILE
        OPEN(UNIT=2,FILE='plot/'//OUTFILE,STATUS='UNKNOWN')
      ELSEIF(ICHKDIR.EQ.1) THEN
        WRITE(6,*) '  '
        WRITE(6,*) 'NAME OF THE OUTPUT FILE : '
        WRITE(6,*) '  '
        READ(5,10) OUTFILE
        OPEN(UNIT=11,FILE='plot/plot_-4_'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=12,FILE='plot/plot_-3_'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=13,FILE='plot/plot_-2_'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=14,FILE='plot/plot_-1_'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=15,FILE='plot/plot_+0_'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=16,FILE='plot/plot_+1_'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=17,FILE='plot/plot_+2_'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=18,FILE='plot/plot_+3_'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=19,FILE='plot/plot_+4_'//OUTFILE,STATUS='UNKNOWN')
      ELSEIF(ICHKDIR.EQ.2) THEN
        WRITE(6,*) '  '
        WRITE(6,*) 'NAME OF THE OUTPUT FILE : '
        WRITE(6,*) '  '
        READ(5,10) OUTFILE
        OPEN(UNIT=2,FILE='plot/'//OUTFILE,STATUS='UNKNOWN')
        OPEN(UNIT=3,FILE='plot/ref_'//OUTFILE,STATUS='UNKNOWN')
      ENDIF
      IF(NSP.EQ.1) GOTO 15
      NSP=1
C
 15   IF((ISPIN.EQ.1).AND.(IS_SP.LE.3)) THEN
        IF(I_SUM.EQ.0) THEN
          WRITE(6,*) '  '
          WRITE(6,*) 'WHICH SPIN CONFIGURATION WOULD YOU LIKE ',
     1               'TO PLOT ? '
          WRITE(6,*) '     *  1  : (+)'
          WRITE(6,*) '     *  2  : (-)'
          WRITE(6,*) '     *  0  : SUM '
          WRITE(6,*) '  '
          READ(5,11) NS2
          IS_SP=IS_SP+1
        ELSE
          NS2=0
        ENDIF
      ENDIF
C
      NS1=1
      IF((IDICHR.GE.1).AND.(IS_DI.LE.3)) THEN
        WRITE(6,*) '  '
        WRITE(6,*) 'WHICH POLARIZATION COMPONENT WOULD YOU LIKE ',
     1             'TO PLOT ? '
        WRITE(6,*) '     *  1  : +/X'
        WRITE(6,*) '     *  2  : -/Y'
        WRITE(6,*) '     *  0  : DICHROIC SIGNAL'
        WRITE(6,*) '  '
        READ(5,11) NS1
        IS_DI=IS_DI+1
      ENDIF
C
      NSO=JSO
      IF(ISPIN.EQ.0) THEN
        IF((I_SO.EQ.2).AND.(IS_SO.LE.3)) THEN
          WRITE(6,*) '  '
          WRITE(6,*) 'WHICH SPIN-ORBIT COMPONENT WOULD YOU LIKE ',
     1               'TO PLOT ? '
          WRITE(6,*) '     *  1  : LI + 1/2'
          WRITE(6,*) '     *  2  : LI - 1/2'
          WRITE(6,*) '     *  0  : SUM'
          WRITE(6,*) '  '
          READ(5,11) NSO
          IS_SO=IS_SO+1
         ENDIF
      ELSE
        IF((I_SO.EQ.2).AND.(IS_SO.LE.2)) THEN
          WRITE(6,*) '  '
          WRITE(6,*) 'WHICH SPIN-ORBIT COMPONENT WOULD YOU LIKE ',
     1               'TO PLOT ? '
          WRITE(6,*) '     *  1  : LI + 1/2'
          WRITE(6,*) '     *  2  : LI - 1/2'
          WRITE(6,*) '  '
          READ(5,11) NSO
          IS_SO=IS_SO+1
        ENDIF
      ENDIF
C
      WRITE(6,*) '  '
      WRITE(6,*) 'WHICH PLANE WOULD YOU LIKE TO PLOT ? '
      WRITE(6,*) '     *  N  : PLANE NUMBER N'
      IF((SPECTRO.EQ.'PED').OR.(SPECTRO.EQ.'AED')) THEN
        WRITE(6,*) '     *  0  : BULK'
      ENDIF
      WRITE(6,*) '     * -1  : TOTAL'
      WRITE(6,*) '  '
      READ(5,11) NUMPLAN
C
      ITEST=ITEST+1
      IF(ITEST.GT.1) GOTO 56
C
      IF((SPECTRO.EQ.'PED').AND.(IE.EQ.0)) THEN
        WRITE(6,*) '  '
        WRITE(6,*) 'WOULD YOU LIKE TO SUBSTRACT THE DIRECT SIGNAL',
     1           ' ? (Y/N) '
        WRITE(6,*) '  '
        READ(5,2) REP1
        IF((REP1.EQ.'y').OR.(REP1.EQ.'Y')) THEN
          ISUB=1
        ELSE
          ISUB=0
        ENDIF
      ELSEIF((SPECTRO.EQ.'AED').AND.(IE.EQ.0)) THEN
        WRITE(6,*) '  '
        WRITE(6,*) 'WOULD YOU LIKE TO SUBSTRACT THE DIRECT SIGNAL',
     1           ' ? (Y/N) '
        WRITE(6,*) '  '
        READ(5,2) REP1
        IF((REP1.EQ.'y').OR.(REP1.EQ.'Y')) THEN
          ISUB=1
        ELSE
          ISUB=0
        ENDIF
      ELSEIF((SPECTRO.EQ.'LED').AND.(IE.EQ.0)) THEN
        ISUB=0
      ELSEIF((SPECTRO.EQ.'XAS').OR.(IE.EQ.1)) THEN
        WRITE(6,*) '  '
        WRITE(6,*) 'WOULD YOU LIKE TO PLOT (I - Io)/Io',
     1           ' ? (Y/N) '
        WRITE(6,*) '  '
        READ(5,2) REP1
        IF((REP1.EQ.'y').OR.(REP1.EQ.'Y')) THEN
          ISUB=2
        ELSE
          ISUB=0
        ENDIF
      ENDIF
C
      IF(SPECTRO.EQ.'LED') THEN
        WRITE(6,*) '  '
        WRITE(6,*) 'WOULD YOU LIKE TO PLOT SEPARATELY THE ',
     1             'KINEMATIC TERM ? (Y/N) '
        WRITE(6,*) '  '
        READ(5,2) REP7
        IF((REP7.EQ.'y').OR.(REP7.EQ.'Y')) THEN
          IDIR=1
          IF(ICHKDIR.EQ.0) THEN
            IUNIT2=4
            OPEN(UNIT=4, FILE='plot/kin_term.dat', STATUS='UNKNOWN')
            WRITE(6,*) '  '
            WRITE(6,*) ' -  --> THE KINEMATIC TERM IS STORED IN FILE ',
     1                 'plot/kin_sterm.dat'
            WRITE(6,*) '  '
          ELSEIF(ICHKDIR.EQ.1) THEN
            OPEN(UNIT=21,FILE='plot/plot_-4_kin.dat',STATUS='UNKNOWN')
            OPEN(UNIT=22,FILE='plot/plot_-3_kin.dat',STATUS='UNKNOWN')
            OPEN(UNIT=23,FILE='plot/plot_-2_kin.dat',STATUS='UNKNOWN')
            OPEN(UNIT=24,FILE='plot/plot_-1_kin.dat',STATUS='UNKNOWN')
            OPEN(UNIT=25,FILE='plot/plot_+0_kin.dat',STATUS='UNKNOWN')
            OPEN(UNIT=26,FILE='plot/plot_+1_kin.dat',STATUS='UNKNOWN')
            OPEN(UNIT=27,FILE='plot/plot_+2_kin.dat',STATUS='UNKNOWN')
            OPEN(UNIT=28,FILE='plot/plot_+3_kin.dat',STATUS='UNKNOWN')
            OPEN(UNIT=29,FILE='plot/plot_+4_kin.dat',STATUS='UNKNOWN')
            WRITE(6,*) '  '
            WRITE(6,*) ' -  --> THE KINEMATIC TERM IS STORED IN THE',
     1                 ' FILES plot/plot_x_kin.dat'
            WRITE(6,*) '  '
          ELSEIF(ICHKDIR.EQ.2) THEN
            OPEN(UNIT=12,FILE='plot/kin_'//OUTFILE,STATUS='UNKNOWN')
            OPEN(UNIT=13,FILE='plot/kin_ref_'//OUTFILE,STATUS='UNKNOWN')
            WRITE(6,*) '  '
            WRITE(6,*) ' -  --> THE KINEMATIC IS STORED IN THE',
     1                 ' FILES plot/kin_xxx.dat'
            WRITE(6,*) '  '
          ENDIF
        ELSE
          IDIR=0
        ENDIF
      ELSE
        WRITE(6,*) '  '
        WRITE(6,*) 'WOULD YOU LIKE TO PLOT SEPARATELY THE DIRECT ',
     1             'SIGNAL ? (Y/N) '
        WRITE(6,*) '  '
        READ(5,2) REP7
        IF((REP7.EQ.'y').OR.(REP7.EQ.'Y')) THEN
          IDIR=1
          IF(ICHKDIR.EQ.0) THEN
            IUNIT2=4
            OPEN(UNIT=4, FILE='plot/dir_signal.dat', STATUS='UNKNOWN')
            WRITE(6,*) '  '
            WRITE(6,*) ' -  --> THE DIRECT SIGNAL IS STORED IN FILE ',
     1                 'plot/dir_signal.dat'
            WRITE(6,*) '  '
          ELSEIF(ICHKDIR.EQ.1) THEN
            OPEN(UNIT=21,FILE='plot/plot_-4_dir.dat',STATUS='UNKNOWN')
            OPEN(UNIT=22,FILE='plot/plot_-3_dir.dat',STATUS='UNKNOWN')
            OPEN(UNIT=23,FILE='plot/plot_-2_dir.dat',STATUS='UNKNOWN')
            OPEN(UNIT=24,FILE='plot/plot_-1_dir.dat',STATUS='UNKNOWN')
            OPEN(UNIT=25,FILE='plot/plot_+0_dir.dat',STATUS='UNKNOWN')
            OPEN(UNIT=26,FILE='plot/plot_+1_dir.dat',STATUS='UNKNOWN')
            OPEN(UNIT=27,FILE='plot/plot_+2_dir.dat',STATUS='UNKNOWN')
            OPEN(UNIT=28,FILE='plot/plot_+3_dir.dat',STATUS='UNKNOWN')
            OPEN(UNIT=29,FILE='plot/plot_+4_dir.dat',STATUS='UNKNOWN')
            WRITE(6,*) '  '
            WRITE(6,*) ' -  --> THE DIRECT SIGNAL IS STORED IN THE',
     1                 ' FILES plot/plot_x_dir.dat'
            WRITE(6,*) '  '
          ELSEIF(ICHKDIR.EQ.2) THEN
            OPEN(UNIT=12,FILE='plot/dir_'//OUTFILE,STATUS='UNKNOWN')
            OPEN(UNIT=13,FILE='plot/dir_ref_'//OUTFILE,STATUS='UNKNOWN')
            WRITE(6,*) '  '
            WRITE(6,*) ' -  --> THE DIRECT SIGNAL IS STORED IN THE',
     1                 ' FILES plot/dir_xxx.dat'
            WRITE(6,*) '  '
          ENDIF
        ELSE
          IDIR=0
        ENDIF
      ENDIF
C
 56   IF(ISPIN.EQ.1) THEN
C
C   Spin-polarized case
C
        DO K=1,LINE_M
          IF(ICHKDIR.EQ.0) THEN
C
C     One direction calculation
C
            IUNIT=2
            IF((I_SUM.EQ.0).AND.(NS2.NE.0)) THEN
C
C       Spin detection
C
              IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                IF(SPECTRO(3:3).EQ.'D') THEN
                  READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                             SR,SF
                ELSEIF(SPECTRO.EQ.'XAS') THEN
                  READ(1,30,END=55) JSO,JPLAN,JS2,ECIN,
     1                              SR,SF
                  DTHETA=1.
                  DPHI=1.
                  E_C=ECIN
                ELSEIF(SPECTRO.EQ.'APC') THEN
                  READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,
     1                             ECIN,DTHETAA,DPHIA,SR,SF
                  IF(NS_E.EQ.1) THEN
                    DTHETA=DTHETAP
                    DPHI=DPHIP
                    DTHETAR=DTHETAA
                    DPHIR=DPHIA
                  ELSEIF(NS_E.EQ.2) THEN
                    DTHETA=DTHETAA
                    DPHI=DPHIA
                    DTHETAR=DTHETAP
                    DPHIR=DPHIP
                  ENDIF
                  IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                    I_JUMP=0
                  ELSE
                    I_JUMP=1
                  ENDIF
                ENDIF
                IF(I_JUMP.EQ.1) GOTO 80
                IF(JS2.EQ.NS2) THEN
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(JSO.EQ.NSO) THEN
                      IF(ISUB.EQ.1) THEN
                        SF=SF-SR
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR).GT.ZERO) THEN
                          SF=(SF-SR)/SR
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      CALL WR(IUNIT,DPHI,DTHETA,ECIN,PHI,THETA,E_C,SF)
                      IF(IDIR.EQ.1) THEN
                        CALL WR(IUNIT2,DPHI,DTHETA,ECIN,PHI,THETA,
     1                          E_C,SR)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
  80            CONTINUE
              ELSE
C
C          Dichroic calculation
C
                IF(SPECTRO(3:3).EQ.'D') THEN
                  READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                              SR1,SF1,SR2,SF2
                ELSEIF(SPECTRO.EQ.'XAS') THEN
                  READ(1,46,END=55) JSO,JPLAN,JS2,ECIN,
     1                              SR1,SF1,SR2,SF2
                  DTHETA=1.
                  DPHI=1.
                  E_C=ECIN
                ELSEIF(SPECTRO.EQ.'APC') THEN
                  READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                              DTHETAA,DPHIA,SR1,SF1,SR2,SF2
                  IF(NS_E.EQ.1) THEN
                    DTHETA=DTHETAP
                    DPHI=DPHIP
                    DTHETAR=DTHETAA
                    DPHIR=DPHIA
                  ELSEIF(NS_E.EQ.2) THEN
                    DTHETA=DTHETAA
                    DPHI=DPHIA
                    DTHETAR=DTHETAP
                    DPHIR=DPHIP
                  ENDIF
                  IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                    I_JUMP=0
                  ELSE
                    I_JUMP=1
                  ENDIF
                ENDIF
                IF(I_JUMP.EQ.1) GOTO 81
                IF(JS2.EQ.NS2) THEN
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(JSO.EQ.NSO) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1=SF1-SR1
                        SF2=SF2-SR2
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1).GT.ZERO) THEN
                          SF1=(SF1-SR1)/SR1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2).GT.ZERO) THEN
                          SF2=(SF2-SR2)/SR2
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      CALL WR_DI(IUNIT,NS1,DPHI,DTHETA,ECIN,PHI,THETA,
     1                           E_C,SF1,SF2)
                      IF(IDIR.EQ.1) THEN
                        CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,PHI,
     1                             THETA,E_C,SR1,SR2)
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
  81            CONTINUE
              ENDIF
            ELSE
C
C        No spin detection ---> sum over the spins
C
              IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                IF(SPECTRO(3:3).EQ.'D') THEN
                  READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                             SR1_1,SF1_1
                  READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                             SR2_1,SF2_1
                ELSEIF(SPECTRO.EQ.'XAS') THEN
                  READ(1,30,END=55) JSO,JPLAN,JS2,ECIN,
     1                             SR1_1,SF1_1
                  READ(1,30,END=55) JSO,JPLAN,JS2,ECIN,
     1                             SR2_1,SF2_1
                  DTHETA=1.
                  DPHI=1.
                  E_C=ECIN
                ELSEIF(SPECTRO.EQ.'APC') THEN
                  READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                             DTHETAA,DPHIA,SR1_1,SF1_1
                  READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                             DTHETAA,DPHIA,SR2_1,SF2_1
                  IF(NS_E.EQ.1) THEN
                    DTHETA=DTHETAP
                    DPHI=DPHIP
                    DTHETAR=DTHETAA
                    DPHIR=DPHIA
                  ELSEIF(NS_E.EQ.2) THEN
                    DTHETA=DTHETAA
                    DPHI=DPHIA
                    DTHETAR=DTHETAP
                    DPHIR=DPHIP
                  ENDIF
                  IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                    I_JUMP=0
                  ELSE
                    I_JUMP=1
                  ENDIF
                ENDIF
                IF(I_JUMP.EQ.1) GOTO 82
                IF(JPLAN.EQ.NUMPLAN) THEN
                  IF(JSO.EQ.NSO) THEN
                    IF(ISUB.EQ.1) THEN
                      SF1_1=SF1_1-SR1_1
                      SF2_1=SF2_1-SR2_1
                    ELSEIF(ISUB.EQ.2) THEN
                      IF(ABS(SR1_1).GT.ZERO) THEN
                        SF1_1=(SF1_1-SR1_1)/SR1_1
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR2_1).GT.ZERO) THEN
                        SF2_1=(SF2_1-SR2_1)/SR2_1
                      ELSE
                        PRINT 9
                      ENDIF
                    ENDIF
                    CALL WR(IUNIT,DPHI,DTHETA,ECIN,PHI,THETA,E_C,
     1                      SF1_1+SF2_1)
                    IF(IDIR.EQ.1) THEN
                      CALL WR(IUNIT2,DPHI,DTHETA,ECIN,PHI,THETA,E_C,
     1                        SR1_1+SR2_1)
                    ENDIF
                  ENDIF
                ENDIF
  82            CONTINUE
              ELSE
C
C          Dichroic calculation
C
                IF(SPECTRO(3:3).EQ.'D') THEN
                  READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                              SR1_1,SF1_1,SR1_2,SF1_2
                  READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                              SR2_1,SF2_1,SR2_2,SF2_2
                ELSEIF(SPECTRO.EQ.'XAS') THEN
                  READ(1,46,END=55) JSO,JPLAN,JS2,ECIN,
     1                              SR1_1,SF1_1,SR1_2,SF1_2
                  READ(1,46,END=55) JSO,JPLAN,JS2,ECIN,
     1                              SR2_1,SF2_1,SR2_2,SF2_2
                  DTHETA=1.
                  DPHI=1.
                  E_C=ECIN
                ELSEIF(SPECTRO.EQ.'APC') THEN
                  READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                              ECIN,DTHETAA,DPHIA,SR1_1,SF1_1,
     2                              SR1_2,SF1_2
                  READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                              ECIN,DTHETAA,DPHIA,SR2_1,SF2_1,
     2                              SR2_2,SF2_2
                  IF(NS_E.EQ.1) THEN
                    DTHETA=DTHETAP
                    DPHI=DPHIP
                    DTHETAR=DTHETAA
                    DPHIR=DPHIA
                  ELSEIF(NS_E.EQ.2) THEN
                    DTHETA=DTHETAA
                    DPHI=DPHIA
                    DTHETAR=DTHETAP
                    DPHIR=DPHIP
                  ENDIF
                  IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                    I_JUMP=0
                  ELSE
                    I_JUMP=1
                  ENDIF
                ENDIF
                IF(I_JUMP.EQ.1) GOTO 83
                IF(JPLAN.EQ.NUMPLAN) THEN
                  IF(JSO.EQ.NSO) THEN
                    IF(ISUB.EQ.1) THEN
                      SF1_1=SF1_1-SR1_1
                      SF2_1=SF2_1-SR2_1
                      SF1_2=SF1_2-SR1_2
                      SF2_2=SF2_2-SR2_2
                    ELSEIF(ISUB.EQ.2) THEN
                      IF(ABS(SR1_1).GT.ZERO) THEN
                        SF1_1=(SF1_1-SR1_1)/SR1_1
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR2_1).GT.ZERO) THEN
                        SF2_1=(SF2_1-SR2_1)/SR2_1
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR1_2).GT.ZERO) THEN
                        SF1_2=(SF1_2-SR1_2)/SR1_2
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR2_2).GT.ZERO) THEN
                        SF2_2=(SF2_2-SR2_2)/SR2_2
                      ELSE
                        PRINT 9
                      ENDIF
                    ENDIF
                    CALL WR_SO_SP(IUNIT,NS1,DPHI,DTHETA,ECIN,PHI,
     1                            THETA,E_C,SF1_1,SF1_2,SF2_1,SF2_2)
                    IF(IDIR.EQ.1) THEN
                      CALL WR_SO_SP(IUNIT2,NS1,DPHI,DTHETA,ECIN,PHI,
     1                              THETA,E_C,SR1_1,SR1_2,SR2_1,SR2_2)
                    ENDIF
                  ENDIF
                ENDIF
  83            CONTINUE
              ENDIF
            ENDIF
          ELSEIF(ICHKDIR.EQ.1) THEN
C
C     Nine directions calculation
C
            DO KANG=11,19
             DO JLN=1,LINE2
              IF((I_SUM.EQ.0).AND.(NS2.NE.0)) THEN
C
C       Spin detection
C
                IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                               SR,SF
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                               DTHETAA,DPHIA,SR,SF
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 84
                  IF(JS2.EQ.NS2) THEN
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(JSO.EQ.NSO) THEN
                        IF(ISUB.EQ.1) THEN
                          SF=SF-SR
                        ELSEIF(ISUB.EQ.2) THEN
                          IF(ABS(SR).GT.ZERO) THEN
                            SF=(SF-SR)/SR
                          ELSE
                            PRINT 9
                          ENDIF
                        ENDIF
                        CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                          ECIN,SF)
                        IF(IDIR.EQ.1) THEN
                          IUNIT2=KANG+10
                          CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                            ECIN,SR)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
  84              CONTINUE
                ELSE
C
C         Dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                                SR1,SF1,SR2,SF2
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,
     1                                ECIN,DTHETAA,DPHIA,SR1,SF1,
     2                                SR2,SF2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 85
                  IF(JS2.EQ.NS2) THEN
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(JSO.EQ.NSO) THEN
                        IF(ISUB.EQ.1) THEN
                          SF1=SF1-SR1
                          SF2=SF2-SR2
                        ELSEIF(ISUB.EQ.2) THEN
                          IF(ABS(SR1).GT.ZERO) THEN
                            SF1=(SF1-SR1)/SR1
                          ELSE
                            PRINT 9
                          ENDIF
                          IF(ABS(SR2).GT.ZERO) THEN
                            SF2=(SF2-SR2)/SR2
                          ELSE
                            PRINT 9
                          ENDIF
                        ENDIF
                        CALL WR_DI(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                             DTHETA,ECIN,SF1,SF2)
                        IF(IDIR.EQ.1) THEN
                          IUNIT2=KANG+10
                          CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                               DTHETA,ECIN,SR1,SR2)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
  85              CONTINUE
                ENDIF
              ELSE
C
C        No spin detection ---> sum over the spins
C
                IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                               SR1_1,SF1_1
                    READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                               SR2_1,SF2_1
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                               DTHETAA,DPHIA,SR1_1,SF1_1
                    READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                               DTHETAA,DPHIA,SR2_1,SF2_1
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 86
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(JSO.EQ.NSO) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1_1=SF1_1-SR1_1
                        SF2_1=SF2_1-SR2_1
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1_1).GT.ZERO) THEN
                          SF1_1=(SF1_1-SR1_1)/SR1_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_1).GT.ZERO) THEN
                          SF2_1=(SF2_1-SR2_1)/SR2_1
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      SF=SF1_1+SF2_1
                      CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,ECIN,
     1                        SF)
                      IF(IDIR.EQ.1) THEN
                        IUNIT2=KANG+10
                        CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                          ECIN,SR)
                      ENDIF
                    ENDIF
                  ENDIF
  86              CONTINUE
                ELSE
C
C         Dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                                SR1_1,SF1_1,SR1_2,SF1_2
                    READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                                SR2_1,SF2_1,SR2_2,SF2_2
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR1_1,SF1_1,
     2                                SR1_2,SF1_2
                    READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR2_1,SF2_1,
     2                                SR2_2,SF2_2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 87
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(JSO.EQ.NSO) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1_1=SF1_1-SR1_1
                        SF2_1=SF2_1-SR2_1
                        SF1_2=SF1_2-SR1_2
                        SF2_2=SF2_2-SR2_2
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1_1).GT.ZERO) THEN
                          SF1_1=(SF1_1-SR1_1)/SR1_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_1).GT.ZERO) THEN
                          SF2_1=(SF2_1-SR2_1)/SR2_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR1_2).GT.ZERO) THEN
                          SF1_2=(SF1_2-SR1_2)/SR1_2
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_2).GT.ZERO) THEN
                          SF2_2=(SF2_2-SR2_2)/SR2_2
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      CALL WR_SO_SP(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                              DTHETA,ECIN,SF1_1,SF1_2,SF2_1,
     2                              SF2_2)
                      IF(IDIR.EQ.1) THEN
                        IUNIT2=KANG+10
                        CALL WR_SO_SP(IUNIT2,NS1,DPHI,DTHETA,ECIN,
     1                                DPHI,DTHETA,ECIN,SR1_1,SR1_2,
     2                                SR2_1,SR2_2)
                      ENDIF
                    ENDIF
                  ENDIF
  87              CONTINUE
                ENDIF
              ENDIF
             ENDDO
            ENDDO
          ELSEIF(ICHKDIR.EQ.2) THEN
C
C     Gaussian averaging calculation ---> average file + normal file
C
            DO KANG=2,3
              IF((I_SUM.EQ.0).AND.(NS2.NE.0)) THEN
C
C       Spin detection
C
                IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                               SR,SF
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,ECIN,
     1                               DTHETAA,DPHIA,SR,SF
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                   ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 88
                  IF(JS2.EQ.NS2) THEN
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(JSO.EQ.NSO) THEN
                        IF(ISUB.EQ.1) THEN
                          SF=SF-SR
                        ELSEIF(ISUB.EQ.2) THEN
                          IF(ABS(SR).GT.ZERO) THEN
                            SF=(SF-SR)/SR
                          ELSE
                            PRINT 9
                          ENDIF
                        ENDIF
                        CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                          ECIN,SF)
                        IF(IDIR.EQ.1) THEN
                          IUNIT2=KANG+10
                          CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,
     1                            DTHETA,ECIN,SR)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
  88              CONTINUE
                ELSE
C
C         Dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                                SR1,SF1,SR2,SF2
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,
     1                                ECIN,DTHETAA,DPHIA,SR1,SF1,
     2                                SR2,SF2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 89
                  IF(JS2.EQ.NS2) THEN
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(JSO.EQ.NSO) THEN
                        IF(ISUB.EQ.1) THEN
                          SF1=SF1-SR1
                          SF2=SF2-SR2
                        ELSEIF(ISUB.EQ.2) THEN
                          IF(ABS(SR1).GT.ZERO) THEN
                            SF1=(SF1-SR1)/SR1
                          ELSE
                            PRINT 9
                          ENDIF
                          IF(ABS(SR2).GT.ZERO) THEN
                            SF2=(SF2-SR2)/SR2
                          ELSE
                            PRINT 9
                          ENDIF
                        ENDIF
                        CALL WR_DI(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                             DTHETA,ECIN,SF1,SF2)
                        IF(IDIR.EQ.1) THEN
                          IUNIT2=KANG+10
                          CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,
     1                               DPHI,DTHETA,ECIN,SR1,SR2)
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
  89              CONTINUE
                ENDIF
              ELSE
C
C        No spin detection ---> sum over the spins
C
                IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                               SR1_1,SF1_1
                    READ(1,3,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                               SR2_1,SF2_1
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,
     1                               ECIN,DTHETAA,DPHIA,SR1_1,SF1_1
                    READ(1,1,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,
     1                               ECIN,DTHETAA,DPHIA,SR2_1,SF2_1
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                     I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 90
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(JSO.EQ.NSO) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1_1=SF1_1-SR1_1
                        SF2_1=SF2_1-SR2_1
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1_1).GT.ZERO) THEN
                          SF1_1=(SF1_1-SR1_1)/SR1_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_1).GT.ZERO) THEN
                          SF2_1=(SF2_1-SR2_1)/SR2_1
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      SF=SF1_1+SF2_1
                      CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,ECIN,
     1                        SF)
                      IF(IDIR.EQ.1) THEN
                        IUNIT2=KANG+10
                        CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                          ECIN,SR)
                      ENDIF
                    ENDIF
                  ENDIF
  90              CONTINUE
                ELSE
C
C         Dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                                SR1_1,SF1_1,SR1_2,SF1_2
                    READ(1,43,END=55) JSO,JPLAN,JS2,DTHETA,DPHI,ECIN,
     1                                SR2_1,SF2_1,SR2_2,SF2_2
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,
     1                                ECIN,DTHETAA,DPHIA,SR1_1,SF1_1,
     2                                SR1_2,SF1_2
                    READ(1,41,END=55) JSO,JPLAN,JS2,DTHETAP,DPHIP,
     1                                ECIN,DTHETAA,DPHIA,SR2_1,SF2_1,
     2                                SR2_2,SF2_2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 91
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(JSO.EQ.NSO) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1_1=SF1_1-SR1_1
                        SF2_1=SF2_1-SR2_1
                        SF1_2=SF1_2-SR1_2
                        SF2_2=SF2_2-SR2_2
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1_1).GT.ZERO) THEN
                          SF1_1=(SF1_1-SR1_1)/SR1_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_1).GT.ZERO) THEN
                          SF2_1=(SF2_1-SR2_1)/SR2_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR1_2).GT.ZERO) THEN
                          SF1_2=(SF1_2-SR1_2)/SR1_2
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_2).GT.ZERO) THEN
                          SF2_2=(SF2_2-SR2_2)/SR2_2
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      CALL WR_SO_SP(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                              DTHETA,ECIN,SF1_1,SF1_2,SF2_1,
     2                              SF2_2)
                      IF(IDIR.EQ.1) THEN
                        IUNIT2=KANG+10
                        CALL WR_SO_SP(IUNIT2,NS1,DPHI,DTHETA,ECIN,
     1                                DPHI,DTHETA,ECIN,SR1_1,SR1_2,
     2                                SR2_1,SR2_2)
                      ENDIF
                    ENDIF
                  ENDIF
  91              CONTINUE
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ELSE
C
C   Spin-independent case
C
        DO K=1,LINE_M
          IF(ICHKDIR.EQ.0) THEN
C
C     One direction calculation
C
            IUNIT=2
            IF(I_SO.EQ.0) THEN
C
C       No spin-orbit resolved initial core state
C
              IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                IF(SPECTRO(3:3).EQ.'D') THEN
                  READ(1,7,END=55) JPLAN,DTHETA,DPHI,ECIN,SR,SF
                ELSEIF(SPECTRO.EQ.'XAS') THEN
                  READ(1,47,END=55) JPLAN,ECIN,SR,SF
                  DTHETA=1.
                  DPHI=1.
                  E_C=ECIN
                ELSEIF(SPECTRO.EQ.'APC') THEN
                  READ(1,8,END=55) JPLAN,DTHETAP,DPHIP,ECIN,
     1                             DTHETAA,DPHIA,SR,SF
                  IF(NS_E.EQ.1) THEN
                    DTHETA=DTHETAP
                    DPHI=DPHIP
                    DTHETAR=DTHETAA
                    DPHIR=DPHIA
                  ELSEIF(NS_E.EQ.2) THEN
                    DTHETA=DTHETAA
                    DPHI=DPHIA
                    DTHETAR=DTHETAP
                    DPHIR=DPHIP
                  ENDIF
                  IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                    I_JUMP=0
                  ELSE
                    I_JUMP=1
                  ENDIF
                ENDIF
                IF(I_JUMP.EQ.1) GOTO 92
                IF(JPLAN.EQ.NUMPLAN) THEN
                  IF(ISUB.EQ.1) THEN
                    SF=SF-SR
                  ELSEIF(ISUB.EQ.2) THEN
                    IF(ABS(SR).GT.ZERO) THEN
                      SF=(SF-SR)/SR
                    ELSE
                      PRINT 9
                    ENDIF
                  ENDIF
                  CALL WR(IUNIT,DPHI,DTHETA,ECIN,PHI,THETA,E_C,SF)
                  IF(IDIR.EQ.1) THEN
                    CALL WR(IUNIT2,DPHI,DTHETA,ECIN,PHI,THETA,E_C,SR)
                  ENDIF
                ENDIF
  92            CONTINUE
              ELSE
C
C         Dichroic calculation
C
                IF(SPECTRO(3:3).EQ.'D') THEN
                  READ(1,17,END=55) JPLAN,DTHETA,DPHI,ECIN,SR1,SF1,
     1                              SR2,SF2
                ELSEIF(SPECTRO.EQ.'XAS') THEN
                  READ(1,57,END=55) JPLAN,ECIN,SR1,SF1,
     1                              SR2,SF2
                  DTHETA=1.
                  DPHI=1.
                  E_C=ECIN
                ELSEIF(SPECTRO.EQ.'APC') THEN
                  READ(1,18,END=55) JPLAN,DTHETAP,DPHIP,ECIN,
     1                              DTHETAA,DPHIA,SR1,SF1,SR2,SF2
                  IF(NS_E.EQ.1) THEN
                    DTHETA=DTHETAP
                    DPHI=DPHIP
                    DTHETAR=DTHETAA
                    DPHIR=DPHIA
                  ELSEIF(NS_E.EQ.2) THEN
                    DTHETA=DTHETAA
                    DPHI=DPHIA
                    DTHETAR=DTHETAP
                    DPHIR=DPHIP
                  ENDIF
                  IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                    I_JUMP=0
                  ELSE
                    I_JUMP=1
                  ENDIF
                ENDIF
                IF(I_JUMP.EQ.1) GOTO 93
                IF(JPLAN.EQ.NUMPLAN) THEN
                  IF(ISUB.EQ.1) THEN
                    SF1=SF1-SR1
                    SF2=SF2-SR2
                  ELSEIF(ISUB.EQ.2) THEN
                    IF(ABS(SR1).GT.ZERO) THEN
                      SF1=(SF1-SR1)/SR1
                    ELSE
                      PRINT 9
                    ENDIF
                    IF(ABS(SR2).GT.ZERO) THEN
                      SF2=(SF2-SR2)/SR2
                    ELSE
                      PRINT 9
                    ENDIF
                  ENDIF
                  CALL WR_DI(IUNIT,NS1,DPHI,DTHETA,ECIN,PHI,THETA,E_C,
     1                     SF1,SF2)
                  IF(IDIR.EQ.1) THEN
                    CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,PHI,THETA,
     1                         E_C,SR1,SR2)
                  ENDIF
                ENDIF
  93            CONTINUE
              ENDIF
            ELSE
C
C       Spin-orbit resolved initial core state
C
              IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                IF(NSO.NE.0) THEN
C
C           One spin-orbit component
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                SR,SF
                  ELSEIF(SPECTRO.EQ.'XAS') THEN
                    READ(1,67,END=55) JSO,JPLAN,ECIN,SR,SF
                    DTHETA=1.
                    DPHI=1.
                    E_C=ECIN
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR,SF
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 94
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(JSO.EQ.NSO) THEN
                      IF(ISUB.EQ.1) THEN
                        SF=SF-SR
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR).GT.ZERO) THEN
                          SF=(SF-SR)/SR
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      CALL WR(IUNIT,DPHI,DTHETA,ECIN,PHI,THETA,E_C,
     1                        SF)
                      IF(IDIR.EQ.1) THEN
                        CALL WR(IUNIT2,DPHI,DTHETA,ECIN,PHI,THETA,
     1                          E_C,SR)
                      ENDIF
                    ENDIF
                  ENDIF
  94              CONTINUE
                ELSE
C
C           Two spin-orbit components ---> sum over them
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                SR1,SF1
                    READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                SR2,SF2
                  ELSEIF(SPECTRO.EQ.'XAS') THEN
                    READ(1,67,END=55) JSO,JPLAN,ECIN,SR1,SF1
                    READ(1,67,END=55) JSO,JPLAN,ECIN,SR2,SF2
                    DTHETA=1.
                    DPHI=1.
                    E_C=ECIN
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR1,SF1
                    READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR2,SF2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 95
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(ISUB.EQ.1) THEN
                      SF1=SF1-SR1
                      SF2=SF2-SR2
                    ELSEIF(ISUB.EQ.2) THEN
                      IF(ABS(SR1).GT.ZERO) THEN
                        SF1=(SF1-SR1)/SR1
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR2).GT.ZERO) THEN
                        SF2=(SF2-SR2)/SR2
                      ELSE
                        PRINT 9
                      ENDIF
                    ENDIF
                    SF=SF1+SF2
                    CALL WR(IUNIT,DPHI,DTHETA,ECIN,PHI,THETA,E_C,SF)
                    IF(IDIR.EQ.1) THEN
                      SR=SR1+SR2
                      CALL WR(IUNIT2,DPHI,DTHETA,ECIN,PHI,THETA,
     1                        E_C,SR)
                    ENDIF
                  ENDIF
  95              CONTINUE
                ENDIF
              ELSE
C
C         Dichroic calculation
C
                IF(NSO.NE.0) THEN
C
C           One spin-orbit component
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                SR1,SF1,SR2,SF2
                  ELSEIF(SPECTRO.EQ.'XAS') THEN
                    READ(1,77,END=55) JSO,JPLAN,ECIN,
     1                                SR1,SF1,SR2,SF2
                    DTHETA=1.
                    DPHI=1.
                    E_C=ECIN
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR1,SF1,SR2,SF2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 96
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(JSO.EQ.NSO) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1=SF1-SR1
                        SF2=SF2-SR2
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1).GT.ZERO) THEN
                          SF1=(SF1-SR1)/SR1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2).GT.ZERO) THEN
                          SF2=(SF2-SR2)/SR2
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      CALL WR_DI(IUNIT,NS1,DPHI,DTHETA,ECIN,PHI,THETA,
     1                           E_C,SF1,SF2)
                      IF(IDIR.EQ.1) THEN
                        CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,PHI,
     1                             THETA,E_C,SR1,SR2)
                      ENDIF
                    ENDIF
                  ENDIF
  96              CONTINUE
                ELSE
C
C           Two spin-orbit components ---> sum over them
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                SR1_1,SF1_1,SR1_2,SF1_2
                    READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                SR2_1,SF2_1,SR2_2,SF2_2
                  ELSEIF(SPECTRO.EQ.'XAS') THEN
                    READ(1,77,END=55) JSO,JPLAN,ECIN,
     1                                SR1_1,SF1_1,SR1_2,SF1_2
                    READ(1,77,END=55) JSO,JPLAN,ECIN,
     1                                SR2_1,SF2_1,SR2_2,SF2_2
                    DTHETA=1.
                    DPHI=1.
                    E_C=ECIN
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR1_1,SF1_1,
     2                                SR1_2,SF1_2
                    READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR2_1,SF2_1,
     2                                SR2_2,SF2_2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 97
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(ISUB.EQ.1) THEN
                      SF1_1=SF1_1-SR1_1
                      SF1_2=SF1_2-SR1_2
                      SF2_1=SF2_1-SR2_1
                      SF2_2=SF2_2-SR2_2
                    ELSEIF(ISUB.EQ.2) THEN
                      IF(ABS(SR1_1).GT.ZERO) THEN
                        SF1_1=(SF1_1-SR1_1)/SR1_1
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR2_1).GT.ZERO) THEN
                        SF2_1=(SF2_1-SR2_1)/SR2_1
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR1_2).GT.ZERO) THEN
                        SF1_2=(SF1_2-SR1_2)/SR1_2
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR2_2).GT.ZERO) THEN
                        SF2_2=(SF2_2-SR2_2)/SR2_2
                      ELSE
                        PRINT 9
                      ENDIF
                    ENDIF
                    SF_1=SF1_1+SF2_1
                    SF_2=SF1_2+SF2_2
                    CALL WR_DI(IUNIT,NS1,DPHI,DTHETA,ECIN,PHI,THETA,
     1                         E_C,SF_1,SF_2)
                    IF(IDIR.EQ.1) THEN
                      SR_1=SR1_1+SR2_1
                      SR_2=SR1_2+SR2_2
                      CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,PHI,
     1                           THETA,E_C,SR_1,SR_2)
                    ENDIF
                  ENDIF
  97              CONTINUE
C
                ENDIF
              ENDIF
            ENDIF
          ELSEIF(ICHKDIR.EQ.1) THEN
C
C     Nine directions calculation
C
           DO KANG=11,19
            DO JLN=1,LINE2
              IF(I_SO.EQ.0) THEN
C
C       No spin-orbit resolved initial core state
C
                IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,7,END=55) JPLAN,DTHETA,DPHI,ECIN,SR,SF
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,8,END=55) JPLAN,DTHETAP,DPHIP,ECIN,
     1                               DTHETAA,DPHIA,SR,SF
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 98
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(ISUB.EQ.1) THEN
                      SF=SF-SR
                    ELSEIF(ISUB.EQ.2) THEN
                      IF(ABS(SR).GT.ZERO) THEN
                        SF=(SF-SR)/SR
                      ELSE
                        PRINT 9
                      ENDIF
                    ENDIF
                    CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,ECIN,SF)
                    IF(IDIR.EQ.1) THEN
                      IUNIT2=KANG+10
                      CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                        ECIN,SR)
                    ENDIF
                  ENDIF
  98              CONTINUE
                ELSE
C
C         Dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,17,END=55) JPLAN,DTHETA,DPHI,ECIN,SR1,SF1,
     1                                SR2,SF2
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,18,END=55) JPLAN,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR1,SF1,SR2,SF2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                    IF(ISUB.EQ.1) THEN
                      SF1=SF1-SR1
                      SF2=SF2-SR2
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 99
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(ISUB.EQ.1) THEN
                      SF1=SF1-SR1
                      SF2=SF2-SR2
                    ELSEIF(ISUB.EQ.2) THEN
                      IF(ABS(SR1).GT.ZERO) THEN
                        SF1=(SF1-SR1)/SR1
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR2).GT.ZERO) THEN
                        SF2=(SF2-SR2)/SR2
                      ELSE
                        PRINT 9
                      ENDIF
                    ENDIF
                    CALL WR_DI(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                         ECIN,SF1,SF2)
                    IF(IDIR.EQ.1) THEN
                      IUNIT2=KANG+10
                      CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                           DTHETA,ECIN,SR1,SR2)
                    ENDIF
                  ENDIF
  99              CONTINUE
                ENDIF
              ELSE
C
C       Spin-orbit resolved initial core state
C
                IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                  IF(NSO.NE.0) THEN
C
C           One spin-orbit component
C
                    IF(SPECTRO(3:3).EQ.'D') THEN
                      READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR,SF
                    ELSEIF(SPECTRO.EQ.'APC') THEN
                      READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR,SF
                      IF(NS_E.EQ.1) THEN
                        DTHETA=DTHETAP
                        DPHI=DPHIP
                        DTHETAR=DTHETAA
                        DPHIR=DPHIA
                      ELSEIF(NS_E.EQ.2) THEN
                        DTHETA=DTHETAA
                        DPHI=DPHIA
                        DTHETAR=DTHETAP
                        DPHIR=DPHIP
                      ENDIF
                      IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                        I_JUMP=0
                      ELSE
                        I_JUMP=1
                      ENDIF
                    ENDIF
                    IF(I_JUMP.EQ.1) GOTO 100
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(JSO.EQ.NSO) THEN
                        IF(ISUB.EQ.1) THEN
                          SF=SF-SR
                        ELSEIF(ISUB.EQ.2) THEN
                          IF(ABS(SR).GT.ZERO) THEN
                            SF=(SF-SR)/SR
                          ELSE
                            PRINT 9
                          ENDIF
                        ENDIF
                        CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                          ECIN,SF)
                        IF(IDIR.EQ.1) THEN
                          IUNIT2=KANG+10
                          CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,
     1                            DTHETA,ECIN,SR)
                        ENDIF
                      ENDIF
                    ENDIF
  100               CONTINUE
                  ELSE
C
C           Two spin-orbit components ---> sum over them
C
                    IF(SPECTRO(3:3).EQ.'D') THEN
                      READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR1,SF1
                      READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR2,SF2
                    ELSEIF(SPECTRO.EQ.'APC') THEN
                      READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR1,SF1
                      READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR2,SF2
                      IF(NS_E.EQ.1) THEN
                        DTHETA=DTHETAP
                        DPHI=DPHIP
                        DTHETAR=DTHETAA
                        DPHIR=DPHIA
                      ELSEIF(NS_E.EQ.2) THEN
                        DTHETA=DTHETAA
                        DPHI=DPHIA
                        DTHETAR=DTHETAP
                        DPHIR=DPHIP
                      ENDIF
                      IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                        I_JUMP=0
                      ELSE
                        I_JUMP=1
                      ENDIF
                    ENDIF
                    IF(ISUB.EQ.1) THEN
                      SF1=SF1-SR1
                      SF2=SF2-SR2
                    ENDIF
                    IF(I_JUMP.EQ.1) GOTO 101
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1=SF1-SR1
                        SF2=SF2-SR2
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1).GT.ZERO) THEN
                          SF1=(SF1-SR1)/SR1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2).GT.ZERO) THEN
                          SF2=(SF2-SR2)/SR2
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      SF=SF1+SF2
                      CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                        ECIN,SF)
                      IF(IDIR.EQ.1) THEN
                        IUNIT2=KANG+10
                        SR=SR1+SR2
                        CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                          ECIN,SR)
                      ENDIF
                    ENDIF
  101               CONTINUE
                  ENDIF
                ELSE
C
C         Dichroic calculation
C
                  IF(NSO.NE.0) THEN
C
C           One spin-orbit component
C
                    IF(SPECTRO(3:3).EQ.'D') THEN
                      READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR1,SF1,SR2,SF2
                    ELSEIF(SPECTRO.EQ.'APC') THEN
                      READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR1,SF1,SR2,SF2
                      IF(NS_E.EQ.1) THEN
                        DTHETA=DTHETAP
                        DPHI=DPHIP
                        DTHETAR=DTHETAA
                        DPHIR=DPHIA
                      ELSEIF(NS_E.EQ.2) THEN
                        DTHETA=DTHETAA
                        DPHI=DPHIA
                        DTHETAR=DTHETAP
                        DPHIR=DPHIP
                      ENDIF
                      IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                        I_JUMP=0
                      ELSE
                        I_JUMP=1
                      ENDIF
                    ENDIF
                    IF(I_JUMP.EQ.1) GOTO 102
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(JSO.EQ.NSO) THEN
                        IF(ISUB.EQ.1) THEN
                          SF1=SF1-SR1
                          SF2=SF2-SR2
                        ELSEIF(ISUB.EQ.2) THEN
                          IF(ABS(SR1).GT.ZERO) THEN
                            SF1=(SF1-SR1)/SR1
                          ELSE
                            PRINT 9
                          ENDIF
                          IF(ABS(SR2).GT.ZERO) THEN
                            SF2=(SF2-SR2)/SR2
                          ELSE
                            PRINT 9
                          ENDIF
                        ENDIF
                        CALL WR_DI(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                             DTHETA,ECIN,SF1,SF2)
                        IF(IDIR.EQ.1) THEN
                          IUNIT2=KANG+10
                          CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,
     1                               DPHI,DTHETA,ECIN,SR1,SR2)
                        ENDIF
                      ENDIF
                    ENDIF
  102               CONTINUE
                  ELSE
C
C           Two spin-orbit components ---> sum over them
C
                    IF(SPECTRO(3:3).EQ.'D') THEN
                      READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR1_1,SF1_1,SR1_2,SF1_2
                      READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR2_1,SF2_1,SR2_2,SF2_2
                    ELSEIF(SPECTRO.EQ.'APC') THEN
                      READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR1_1,SF1_1,
     2                                  SR1_2,SF1_2
                      READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR2_1,SF2_1,
     2                                  SR2_2,SF2_2
                      IF(NS_E.EQ.1) THEN
                        DTHETA=DTHETAP
                        DPHI=DPHIP
                        DTHETAR=DTHETAA
                        DPHIR=DPHIA
                      ELSEIF(NS_E.EQ.2) THEN
                        DTHETA=DTHETAA
                        DPHI=DPHIA
                        DTHETAR=DTHETAP
                        DPHIR=DPHIP
                      ENDIF
                      IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                        I_JUMP=0
                      ELSE
                        I_JUMP=1
                      ENDIF
                    ENDIF
                    IF(I_JUMP.EQ.1) GOTO 103
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1_1=SF1_1-SR1_1
                        SF1_2=SF1_2-SR1_2
                        SF2_1=SF2_1-SR2_1
                        SF2_2=SF2_2-SR2_2
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1_1).GT.ZERO) THEN
                          SF1_1=(SF1_1-SR1_1)/SR1_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_1).GT.ZERO) THEN
                          SF2_1=(SF2_1-SR2_1)/SR2_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR1_2).GT.ZERO) THEN
                          SF1_2=(SF1_2-SR1_2)/SR1_2
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_2).GT.ZERO) THEN
                          SF2_2=(SF2_2-SR2_2)/SR2_2
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      SF1=SF1_1+SF2_1
                      SF2=SF1_2+SF2_2
                      CALL WR_DI(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                           DTHETA,ECIN,SF1,SF2)
                      IF(IDIR.EQ.1) THEN
                        IUNIT2=KANG+10
                        SR1=SR1_1+SR2_1
                        SR2=SR1_2+SR2_2
                        CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                             DTHETA,ECIN,SR1,SR2)
                      ENDIF
                    ENDIF
  103               CONTINUE
                  ENDIF
                ENDIF
              ENDIF
             ENDDO
            ENDDO
          ELSEIF(ICHKDIR.EQ.2) THEN
C
C     Gaussian averaging calculation ---> average file + normal file
C
            DO KANG=2,3
              IF(I_SO.EQ.0) THEN
C
C       No spin-orbit resolved initial core state
C
                IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,7,END=55) JPLAN,DTHETA,DPHI,ECIN,SR,SF
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,8,END=55) JPLAN,DTHETAP,DPHIP,ECIN,
     1                               DTHETAA,DPHIA,SR,SF
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 104
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(ISUB.EQ.1) THEN
                      SF=SF-SR
                    ELSEIF(ISUB.EQ.2) THEN
                      IF(ABS(SR).GT.ZERO) THEN
                        SF=(SF-SR)/SR
                      ELSE
                        PRINT 9
                      ENDIF
                    ENDIF
                    CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,ECIN,SF)
                    IF(IDIR.EQ.1) THEN
                      IUNIT2=KANG+10
                      CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                        ECIN,SR)
                    ENDIF
                  ENDIF
  104             CONTINUE
                ELSE
C
C         Dichroic calculation
C
                  IF(SPECTRO(3:3).EQ.'D') THEN
                    READ(1,17,END=55) JPLAN,DTHETA,DPHI,ECIN,SR1,SF1,
     1                                SR2,SF2
                  ELSEIF(SPECTRO.EQ.'APC') THEN
                    READ(1,18,END=55) JPLAN,DTHETAP,DPHIP,ECIN,
     1                                DTHETAA,DPHIA,SR1,SF1,SR2,SF2
                    IF(NS_E.EQ.1) THEN
                      DTHETA=DTHETAP
                      DPHI=DPHIP
                      DTHETAR=DTHETAA
                      DPHIR=DPHIA
                    ELSEIF(NS_E.EQ.2) THEN
                      DTHETA=DTHETAA
                      DPHI=DPHIA
                      DTHETAR=DTHETAP
                      DPHIR=DPHIP
                    ENDIF
                    IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                      I_JUMP=0
                    ELSE
                      I_JUMP=1
                    ENDIF
                  ENDIF
                  IF(ISUB.EQ.1) THEN
                    SF1=SF1-SR1
                    SF2=SF2-SR2
                  ENDIF
                  IF(I_JUMP.EQ.1) GOTO 105
                  IF(JPLAN.EQ.NUMPLAN) THEN
                    IF(ISUB.EQ.1) THEN
                      SF1=SF1-SR1
                      SF2=SF2-SR2
                    ELSEIF(ISUB.EQ.2) THEN
                      IF(ABS(SR1).GT.ZERO) THEN
                        SF1=(SF1-SR1)/SR1
                      ELSE
                        PRINT 9
                      ENDIF
                      IF(ABS(SR2).GT.ZERO) THEN
                        SF2=(SF2-SR2)/SR2
                      ELSE
                        PRINT 9
                      ENDIF
                    ENDIF
                    CALL WR_DI(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                         ECIN,SF1,SF2)
                    IF(IDIR.EQ.1) THEN
                      IUNIT2=KANG+10
                      CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                           DTHETA,ECIN,SR1,SR2)
                    ENDIF
                  ENDIF
  105             CONTINUE
                ENDIF
              ELSE
C
C       Spin-orbit resolved initial core state
C
                IF(IDICHR.EQ.0) THEN
C
C         No dichroic calculation
C
                  IF(NSO.NE.0) THEN
C
C           One spin-orbit component
C
                    IF(SPECTRO(3:3).EQ.'D') THEN
                      READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR,SF
                    ELSEIF(SPECTRO.EQ.'APC') THEN
                      READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR,SF
                      IF(NS_E.EQ.1) THEN
                        DTHETA=DTHETAP
                        DPHI=DPHIP
                        DTHETAR=DTHETAA
                        DPHIR=DPHIA
                      ELSEIF(NS_E.EQ.2) THEN
                        DTHETA=DTHETAA
                        DPHI=DPHIA
                        DTHETAR=DTHETAP
                        DPHIR=DPHIP
                      ENDIF
                      IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                        I_JUMP=0
                      ELSE
                        I_JUMP=1
                      ENDIF
                    ENDIF
                    IF(I_JUMP.EQ.1) GOTO 106
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(JSO.EQ.NSO) THEN
                        IF(ISUB.EQ.1) THEN
                          SF=SF-SR
                        ELSEIF(ISUB.EQ.2) THEN
                          IF(ABS(SR).GT.ZERO) THEN
                            SF=(SF-SR)/SR
                          ELSE
                            PRINT 9
                          ENDIF
                        ENDIF
                        CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                          ECIN,SF)
                        IF(IDIR.EQ.1) THEN
                          IUNIT2=KANG+10
                          CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,
     1                            DTHETA,ECIN,SR)
                        ENDIF
                      ENDIF
                    ENDIF
  106               CONTINUE
                  ELSE
C
C           Two spin-orbit components ---> sum over them
C
                    IF(SPECTRO(3:3).EQ.'D') THEN
                      READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR1,SF1
                      READ(1,27,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR2,SF2
                    ELSEIF(SPECTRO.EQ.'APC') THEN
                      READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR1,SF1
                      READ(1,28,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR2,SF2
                      IF(NS_E.EQ.1) THEN
                        DTHETA=DTHETAP
                        DPHI=DPHIP
                        DTHETAR=DTHETAA
                        DPHIR=DPHIA
                      ELSEIF(NS_E.EQ.2) THEN
                        DTHETA=DTHETAA
                        DPHI=DPHIA
                        DTHETAR=DTHETAP
                        DPHIR=DPHIP
                      ENDIF
                      IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                        I_JUMP=0
                      ELSE
                        I_JUMP=1
                      ENDIF
                    ENDIF
                    IF(ISUB.EQ.1) THEN
                      SF1=SF1-SR1
                      SF2=SF2-SR2
                    ENDIF
                    IF(I_JUMP.EQ.1) GOTO 107
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1=SF1-SR1
                        SF2=SF2-SR2
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1).GT.ZERO) THEN
                          SF1=(SF1-SR1)/SR1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2).GT.ZERO) THEN
                          SF2=(SF2-SR2)/SR2
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      SF=SF1+SF2
                      CALL WR(KANG,DPHI,DTHETA,ECIN,DPHI,DTHETA,
     1                        ECIN,SF)
                      IF(IDIR.EQ.1) THEN
                        IUNIT2=KANG+10
                        CALL WR(IUNIT2,DPHI,DTHETA,ECIN,DPHI,
     1                          DTHETA,ECIN,SR)
                      ENDIF
                    ENDIF
  107               CONTINUE
                  ENDIF
                ELSE
C
C         Dichroic calculation
C
                  IF(NSO.NE.0) THEN
C
C           One spin-orbit component
C
                    IF(SPECTRO(3:3).EQ.'D') THEN
                      READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR1,SF1,SR2,SF2
                    ELSEIF(SPECTRO.EQ.'APC') THEN
                      READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR1,SF1,SR2,SF2
                      IF(NS_E.EQ.1) THEN
                        DTHETA=DTHETAP
                        DPHI=DPHIP
                        DTHETAR=DTHETAA
                        DPHIR=DPHIA
                      ELSEIF(NS_E.EQ.2) THEN
                        DTHETA=DTHETAA
                        DPHI=DPHIA
                        DTHETAR=DTHETAP
                        DPHIR=DPHIP
                      ENDIF
                      IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                        I_JUMP=0
                      ELSE
                        I_JUMP=1
                      ENDIF
                    ENDIF
                    IF(I_JUMP.EQ.1) GOTO 108
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(JSO.EQ.NSO) THEN
                        IF(ISUB.EQ.1) THEN
                          SF1=SF1-SR1
                          SF2=SF2-SR2
                        ELSEIF(ISUB.EQ.2) THEN
                          IF(ABS(SR1).GT.ZERO) THEN
                            SF1=(SF1-SR1)/SR1
                          ELSE
                            PRINT 9
                          ENDIF
                          IF(ABS(SR2).GT.ZERO) THEN
                            SF2=(SF2-SR2)/SR2
                          ELSE
                            PRINT 9
                          ENDIF
                        ENDIF
                        CALL WR_DI(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                             DTHETA,ECIN,SF1,SF2)
                        IF(IDIR.EQ.1) THEN
                          IUNIT2=KANG+10
                          CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,
     1                               DPHI,THETA,ECIN,SR1,SR2)
                        ENDIF
                      ENDIF
                    ENDIF
  108               CONTINUE
                  ELSE
C
C           Two spin-orbit components ---> sum over them
C
                    IF(SPECTRO(3:3).EQ.'D') THEN
                      READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR1_1,SF1_1,SR1_2,SF1_2
                      READ(1,37,END=55) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                                  SR2_1,SF2_1,SR2_2,SF2_2
                    ELSEIF(SPECTRO.EQ.'APC') THEN
                      READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR1_1,SF1_1,
     2                                  SR1_2,SF1_2
                      READ(1,38,END=55) JSO,JPLAN,DTHETAP,DPHIP,ECIN,
     1                                  DTHETAA,DPHIA,SR2_1,SF2_1,
     2                                  SR2_2,SF2_2
                      IF(NS_E.EQ.1) THEN
                        DTHETA=DTHETAP
                        DPHI=DPHIP
                        DTHETAR=DTHETAA
                        DPHIR=DPHIA
                      ELSEIF(NS_E.EQ.2) THEN
                        DTHETA=DTHETAA
                        DPHI=DPHIA
                        DTHETAR=DTHETAP
                        DPHIR=DPHIP
                      ENDIF
                      IF((DTHETAR.EQ.THETAF).AND.(DPHIR.EQ.PHIF)) THEN
                        I_JUMP=0
                      ELSE
                        I_JUMP=1
                      ENDIF
                    ENDIF
                    IF(I_JUMP.EQ.1) GOTO 109
                    IF(JPLAN.EQ.NUMPLAN) THEN
                      IF(ISUB.EQ.1) THEN
                        SF1_1=SF1_1-SR1_1
                        SF1_2=SF1_2-SR1_2
                        SF2_1=SF2_1-SR2_1
                        SF2_2=SF2_2-SR2_2
                      ELSEIF(ISUB.EQ.2) THEN
                        IF(ABS(SR1_1).GT.ZERO) THEN
                          SF1_1=(SF1_1-SR1_1)/SR1_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_1).GT.ZERO) THEN
                          SF2_1=(SF2_1-SR2_1)/SR2_1
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR1_2).GT.ZERO) THEN
                          SF1_2=(SF1_2-SR1_2)/SR1_2
                        ELSE
                          PRINT 9
                        ENDIF
                        IF(ABS(SR2_2).GT.ZERO) THEN
                          SF2_2=(SF2_2-SR2_2)/SR2_2
                        ELSE
                          PRINT 9
                        ENDIF
                      ENDIF
                      SF1=SF1_1+SF2_1
                      SF2=SF1_2+SF2_2
                      CALL WR_DI(KANG,NS1,DPHI,DTHETA,ECIN,DPHI,
     1                           DTHETA,ECIN,SF1,SF2)
                      IF(IDIR.EQ.1) THEN
                        IUNIT2=KANG+10
                        SR1=SR1_1+SR2_1
                        SR2=SR1_2+SR2_2
                        CALL WR_DI(IUNIT2,NS1,DPHI,DTHETA,ECIN,
     1                             DPHI,DTHETA,ECIN,SR1,SR2)
                      ENDIF
                    ENDIF
  109               CONTINUE
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
 55   CLOSE(2)
      WRITE(6,5) OUTFILE
      WRITE(6,*) '   '
C
      IF(ICHKDIR.EQ.2) THEN
        PRINT *,'               '
        PRINT *,'               '
        PRINT *,'        Gaussian averaged output file : ',
     1          OUTFILE
        PRINT *,'        Non averaged output file      : ',
     1          'ref_'//OUTFILE
        PRINT *,'               '
        PRINT *,'               '
      ENDIF
C
      IF(ISPIN.EQ.0) GOTO 20
      IF(I_SUM.EQ.1) GOTO 12
C
      WRITE(6,*) '   '
      WRITE(6,*) 'WOULD YOU LIKE TO PLOT ANOTHER SPIN CONFIGURATION',
     1           ' ? (Y/N) '
      WRITE(6,*) '   '
C
      READ(5,2) REP4
      IF((REP4.EQ.'y').OR.(REP4.EQ.'Y')) THEN
        REWIND 1
        DO JLINE=1,NHEAD2
          READ(1,*) HEAD
        ENDDO
        GOTO 23
      ENDIF
C
 12   IF(I_SO.EQ.2) THEN
        WRITE(6,*) '   '
        WRITE(6,*) 'WOULD YOU LIKE TO PLOT ANOTHER SPIN-ORBIT ',
     1             'COMPONENT ? (Y/N) '
        WRITE(6,*) '   '
C
        READ(5,2) REP6
        IF((REP6.EQ.'y').OR.(REP6.EQ.'Y')) THEN
          REWIND 1
          DO JLINE=1,NHEAD2
            READ(1,*) HEAD
          ENDDO
          GOTO 23
        ENDIF
      ENDIF
C
 20   IF(I_PO.GT.0) THEN
        WRITE(6,*) '   '
        WRITE(6,*) 'WOULD YOU LIKE TO PLOT ANOTHER EC/THETA/PHI',
     1             ' CONFIGURATION ? (Y/N) '
        WRITE(6,*) '   '
C
        READ(5,2) REP5
        IF((REP5.EQ.'y').OR.(REP5.EQ.'Y')) THEN
          REWIND 1
          DO JLINE=1,NHEAD2
            READ(1,*) HEAD
          ENDDO
          GOTO 19
        ENDIF
      ENDIF
C
      WRITE(6,*) '   '
      WRITE(6,*) 'WOULD YOU LIKE TO PLOT ANOTHER PLANE ? (Y/N) '
      WRITE(6,*) '   '
C
      READ(5,2) REP2
      IF((REP2.EQ.'y').OR.(REP2.EQ.'Y')) THEN
        REWIND 1
        DO JLINE=1,NHEAD2
          READ(1,*) HEAD
        ENDDO
        GOTO 13
      ENDIF
      CLOSE(1)
C
      WRITE(6,*) '   '
      WRITE(6,*) 'WOULD YOU LIKE TO TREAT ANOTHER FILE ? (Y/N) '
      WRITE(6,*) '   '
C
      READ(5,2) REP3
      IF((REP3.EQ.'y').OR.(REP3.EQ.'Y')) GOTO 14
C
  1   FORMAT(2X,I3,2X,I2,2X,I1,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,
     1       F6.2,2X,E12.6,2X,E12.6)
  2   FORMAT(A1)
  3   FORMAT(2X,I3,2X,I2,2X,I1,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,
     1       E12.6)
  5   FORMAT(/,'     ---> THE RESULT HAS BEEN WRITTEN IN FILE ',A24)
  6   FORMAT(//,'!!!  DIMENSIONING ERROR : INCREASE N_HEAD !!!',//)
  7   FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
  8   FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,2X,E12.6,
     1       2X,E12.6)
  9   FORMAT(//,'!!!  ERROR : THE OUTPUT IS NULL !!!',//)
 10   FORMAT(A24)
 11   FORMAT(I3)
 16   FORMAT(2X,A3)
 17   FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,
     1       2X,E12.6,2X,E12.6)
 18   FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,2X,
     1       E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
 27   FORMAT(2X,I1,2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
 28   FORMAT(2X,I1,2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,2X,
     1       E12.6,2X,E12.6)
 30   FORMAT(2X,I3,2X,I2,2X,I1,2X,F8.2,2X,E12.6,2X,E12.6)
 31   FORMAT(' THERE ARE ',I4,' AUGER ELECTRON DETECTOR DIRECTIONS.')
 32   FORMAT(' THERE ARE ',I4,' PHOTOELECTRON DETECTOR DIRECTIONS.')
 33   FORMAT(8(2X,I1))
 34   FORMAT(/,' THE DIRECTION YOU HAVE CHOSEN FOR THE FIXED DETECTOR',
     1       ' IS :',//,'     Ec = ',F8.2,'   THETA = ',F6.2,
     2       '   PHI = ',F6.2)
 35   FORMAT(/,' THE DIRECTION YOU HAVE CHOSEN FOR THE FIXED DETECTOR',
     1       ' IS :',//,'   THETA = ',F6.2,'   PHI = ',F6.2)
 37   FORMAT(2X,I1,2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,
     1       E12.6,2X,E12.6,2X,E12.6)
 38   FORMAT(2X,I1,2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,2X,
     1       E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
 40   FORMAT(///,'<<<<<<<<<<  LINE_M SHOULD BE AT LEAST : ',I8,
     1       '  >>>>>>>>>>')
 41   FORMAT(2X,I1,2X,I3,2X,I1,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,
     1       2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
 43   FORMAT(2X,I1,2X,I3,2X,I1,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,
     1       E12.6,2X,E12.6,2X,E12.6)
 44   FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
 45   FORMAT(I4,2X,I3,2X,I1)
 46   FORMAT(2X,I1,2X,I3,2X,I1,2X,2X,F8.2,2X,E12.6,2X,
     1       E12.6,2X,E12.6,2X,E12.6)
 47   FORMAT(2X,I3,2X,F8.2,2X,E12.6,2X,E12.6)
 57   FORMAT(2X,I3,2X,F8.2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
 67   FORMAT(2X,I1,2X,I3,2X,F8.2,2X,E12.6,2X,E12.6)
 77   FORMAT(2X,I1,2X,I3,2X,F8.2,2X,E12.6,2X,
     1       E12.6,2X,E12.6,2X,E12.6)
 114  FORMAT(' ',/,
     1       '     ===============================================',/,
     2       '                                                    ',/,
     3       '             PROCESSING A CROSS-SECTION FILE        ',/,
     4       '                     BEFORE PLOTTING                ',/,
     5       '                                                    ',/,
     6       '     ===============================================',/,
     7       ' ',//)
 888  FORMAT(A6)
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE WR_SO(IUNIT,NS1,VAR,RES1,RES2)
C
      DATA ZERO /1.E-35/
C
      IF(NS1.EQ.1) THEN
        WRITE(IUNIT,8) VAR,RES1
      ELSEIF(NS1.EQ.2) THEN
        WRITE(IUNIT,8) VAR,RES2
      ELSE
        IF(ABS(RES1+RES2).LT.ZERO) THEN
          PRINT 7
          STOP
        ELSE
          RES=(RES1-RES2)/(RES1+RES2)
          WRITE(IUNIT,8) VAR,RES
        ENDIF
      ENDIF
C
  7   FORMAT('!!!  ERROR : THE OUTPUT IS NULL !!!')
  8   FORMAT(3X,F6.2,3X,E12.6)
C
      RETURN
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE WR(IUNIT,VAR1,VAR2,VAR3,VAL1,VAL2,VAL3,RES)
C
      COMMON /TYPE_EXP/ IPHI,ITHETA,IE
C
      IF(IPHI.EQ.1) THEN
        IF((VAR3.EQ.VAL3).AND.(VAR2.EQ.VAL2)) THEN
          WRITE(IUNIT,8) VAR1,RES
        ENDIF
      ELSEIF(ITHETA.EQ.1) THEN
        IF((VAR3.EQ.VAL3).AND.(VAR1.EQ.VAL1)) THEN
          WRITE(IUNIT,8) VAR2,RES
        ENDIF
      ELSEIF(IE.EQ.1) THEN
        IF((VAR2.EQ.VAL2).AND.(VAR1.EQ.VAL1)) THEN
          WRITE(IUNIT,8) VAR3,RES
        ENDIF
      ENDIF
C
  8   FORMAT(3X,F6.2,3X,E12.6)
C
      RETURN
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE WR_DI(IUNIT,NS1,VAR1,VAR2,VAR3,VAL1,VAL2,VAL3,
     1                 RES1,RES2)
C
      COMMON /TYPE_EXP/ IPHI,ITHETA,IE
C
      IF(IPHI.EQ.1) THEN
        IF((VAR3.EQ.VAL3).AND.(VAL2.EQ.VAR2)) THEN
          CALL WR_SO(IUNIT,NS1,VAR1,RES1,RES2)
        ENDIF
      ELSEIF(ITHETA.EQ.1) THEN
        IF((VAR3.EQ.VAL3).AND.(VAR1.EQ.VAL1)) THEN
          CALL WR_SO(IUNIT,NS1,VAR2,RES1,RES2)
        ENDIF
      ELSEIF(IE.EQ.1) THEN
        IF((VAR2.EQ.VAL2).AND.(VAR1.EQ.VAL1)) THEN
          CALL WR_SO(IUNIT,NS1,VAR3,RES1,RES2)
        ENDIF
      ENDIF
C
      RETURN
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE WR_SO_SP(IUNIT,NS1,VAR1,VAR2,VAR3,VAL1,VAL2,VAL3,
     1                    RES1_1,RES1_2,RES2_1,RES2_2)
C
C   This subroutines writes the output in the dichroic case
C   but with no spin detection. Hence, it sums over the spins.
C
C                  NS1 = 1 : +/x polarization plotted
C                  NS1 = 2 : -/y polarization plotted
C                  NS1 = 0 : dichroic signal plotted
C
C   The input results are RESn_m whre n corresponds to the spin
C               and m to the polarization
C
C   VARn is the scanned variable read in the output file from
C   the calculation code and VALn is the value selected by the
C                 user when running this program.
C
C
      COMMON /TYPE_EXP/ IPHI,ITHETA,IE
C
      DATA ZERO /1.E-35/
C
      IF(NS1.EQ.1) THEN
        RES=RES1_1+RES2_1
        IF(IPHI.EQ.1) THEN
          IF((VAR3.EQ.VAL3).AND.(VAR2.EQ.VAL2)) THEN
            WRITE(IUNIT,8) VAR1,RES
          ENDIF
        ELSEIF(ITHETA.EQ.1) THEN
          IF((VAR3.EQ.VAL3).AND.(VAR1.EQ.VAL1)) THEN
            WRITE(IUNIT,8) VAR2,RES
          ENDIF
        ELSEIF(IE.EQ.1) THEN
          IF((VAR2.EQ.VAL2).AND.(VAR1.EQ.VAL1)) THEN
            WRITE(IUNIT,8) VAR3,RES
          ENDIF
        ENDIF
      ELSEIF(NS1.EQ.2) THEN
        RES=RES1_2+RES2_2
        IF(IPHI.EQ.1) THEN
          IF((VAR3.EQ.VAL3).AND.(VAR2.EQ.VAL2)) THEN
            WRITE(IUNIT,8) VAR1,RES
          ENDIF
        ELSEIF(ITHETA.EQ.1) THEN
          IF((VAR3.EQ.VAL3).AND.(VAR1.EQ.VAL1)) THEN
            WRITE(IUNIT,8) VAR2,RES
          ENDIF
        ELSEIF(IE.EQ.1) THEN
          IF((VAR2.EQ.VAL2).AND.(VAR1.EQ.VAL1)) THEN
            WRITE(IUNIT,8) VAR3,RES
          ENDIF
        ENDIF
      ELSE
        RES1=RES1_1+RES2_1
        RES2=RES1_2+RES2_2
        IF(ABS(RES1+RES2).LT.ZERO) THEN
          PRINT 7
          STOP
        ELSE
          RES=(RES1-RES2)/(RES1+RES2)
        ENDIF
        IF(IPHI.EQ.1) THEN
          IF((VAR3.EQ.VAL3).AND.(VAR2.EQ.VAL2)) THEN
            WRITE(IUNIT,8) VAR1,RES
          ENDIF
        ELSEIF(ITHETA.EQ.1) THEN
          IF((VAR3.EQ.VAL3).AND.(VAR1.EQ.VAL1)) THEN
            WRITE(IUNIT,8) VAR2,RES
          ENDIF
        ELSEIF(IE.EQ.1) THEN
          IF((VAR2.EQ.VAL2).AND.(VAR1.EQ.VAL1)) THEN
            WRITE(IUNIT,8) VAR3,RES
          ENDIF
        ENDIF
      ENDIF
C
  7   FORMAT('!!!  ERROR : THE OUTPUT IS NULL !!!')
  8   FORMAT(3X,F6.2,3X,E12.6)
C
      RETURN
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE FIXED_VAR(IE,NE,ITHETA,NTHETA,IPHI,NPHI,NS_E,
     1                     ISPIN,I_SO,NLINE,I_EC,I_TH,I_PH,NHEAD2)
C
C  This subroutine checks the different values available for
C     the fixed variables in the spec result file.
C
C     Electron corresponding to the fixed variables :
C
C            NS_E = 0 ---> PhD/LEED/AED case
C            NS_E = 1 ---> photoelectron in the APECS case
C            NS_E = 2 ---> Auger electron in the APECS case
C
C
      PARAMETER (NPOINTS=1000)
C
      CHARACTER*1 NDUM
C
      COMMON /NUMBERS/ NTHETAP,NPHIP,NTHETAA,NPHIA
      COMMON /VALUES/ EC(NPOINTS),TH(NPOINTS),PH(NPOINTS)
C
      I_STOP=0
C
      I_EC=0
      I_TH=0
      I_PH=0
C
      IF((NPHIP.GT.NPOINTS).OR.(NPHIA.GT.NPOINTS)) I_STOP=1
      IF((NTHETAP.GT.NPOINTS).OR.(NTHETAA.GT.NPOINTS)) I_STOP=1
      IF(NE.GT.NPOINTS) I_STOP=1
      IF(I_STOP.EQ.1) THEN
         WRITE(*,10) MAX(NPHIP,NPHIA,NTHETAP,NTHETAA,NE)
         STOP
      ENDIF
C
      DO JP=1,NPOINTS
        EC(JP)=0.
        TH(JP)=0.
        PH(JP)=0.
      ENDDO
C
      IF((IE.EQ.0).AND.(NE.GT.1)) THEN
        REWIND 1
        DO JDUM=1,NHEAD2
          READ(1,30) NDUM
        ENDDO
        IF(NS_E.GT.0) THEN
          READ(1,30) NDUM
          READ(1,30) NDUM
        ENDIF
        I_EC=1
        IF(NS_E.EQ.0) THEN
          NSTEP=NTHETA*NPHI*NLINE
        ELSE
          NSTEP=NTHETAP*NPHIP*NTHETAA*NPHIA*NLINE
        ENDIF
        DO JE=1,NE
          IF(ISPIN.EQ.0) THEN
            IF(I_SO.EQ.0) THEN
              IF(NS_E.EQ.0) THEN
                READ(1,7) JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1
              ELSE
                READ(1,8) JPLAN,DTHETA,DPHI,ECIN,
     1                    DTHETAA,DPHIA,SR1_1,SF1_1
              ENDIF
            ELSE
              IF(NS_E.EQ.0) THEN
                READ(1,27) JSO,JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1
              ELSE
                READ(1,28) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                     DTHETAA,DPHIA,SR1_1,SF1_1
              ENDIF
            ENDIF
          ELSE
            IF(NS_E.EQ.0) THEN
              READ(1,3) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,SR1_1,SF1_1
            ELSE
              READ(1,1) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,
     1                  DTHETAA,DPHIA,SR1_1,SF1_1
            ENDIF
          ENDIF
          DO JSTEP=1,NSTEP-1
            READ(1,*) JPLAN
          ENDDO
          EC(JE)=ECIN
        ENDDO
        TH(1)=DTHETA
        PH(1)=DPHI
      ENDIF
      IF((ITHETA.EQ.0).AND.(NTHETA.GT.1))THEN
        REWIND 1
        DO JDUM=1,NHEAD2
          READ(1,30) NDUM
        ENDDO
        IF(NS_E.GT.0) THEN
          READ(1,30) NDUM
          READ(1,30) NDUM
        ENDIF
        I_TH=1
        IF(NS_E.EQ.0) THEN
          NSTEP=NPHI*NLINE
        ELSEIF(NS_E.EQ.1) THEN
          NSTEP=NPHIA*NLINE
        ELSEIF(NS_E.EQ.2) THEN
          NSTEP=NPHIP*NTHETAA*NPHIA*NLINE
        ENDIF
        DO JTHETA=1,NTHETA
          IF(ISPIN.EQ.0) THEN
            IF(I_SO.EQ.0) THEN
              IF(NS_E.EQ.0) THEN
                READ(1,7) JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1
              ELSE
                READ(1,8) JPLAN,DTHETA,DPHI,ECIN,
     1                    DTHETAA,DPHIA,SR1_1,SF1_1
              ENDIF
            ELSE
              IF(NS_E.EQ.0) THEN
                READ(1,27) JSO,JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1
              ELSE
                READ(1,28) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                     DTHETAA,DPHIA,SR1_1,SF1_1
              ENDIF
            ENDIF
          ELSE
            IF(NS_E.EQ.0) THEN
              READ(1,3) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,SR1_1,SF1_1
            ELSE
              READ(1,1) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,
     1                  DTHETAA,DPHIA,SR1_1,SF1_1
            ENDIF
          ENDIF
          DO JSTEP=1,NSTEP-1
            READ(1,*) JPLAN
          ENDDO
          IF(NS_E.EQ.0) THEN
            TH(JTHETA)=DTHETA
          ELSEIF(NS_E.EQ.1) THEN
            TH(JTHETA)=DTHETAA
          ELSEIF(NS_E.EQ.2) THEN
            TH(JTHETA)=DTHETA
          ENDIF
        ENDDO
        EC(1)=ECIN
        PH(1)=DPHI
      ENDIF
      IF((IPHI.EQ.0).AND.(NPHI.GT.1)) THEN
        REWIND 1
        DO JDUM=1,NHEAD2
          READ(1,30) NDUM
        ENDDO
        IF(NS_E.GT.0) THEN
          READ(1,30) NDUM
          READ(1,30) NDUM
        ENDIF
        I_PH=1
        IF(NS_E.EQ.0) THEN
          NSTEP=NLINE
        ELSEIF(NS_E.EQ.1) THEN
          NSTEP=NLINE
        ELSEIF(NS_E.EQ.2) THEN
          NSTEP=NTHETAA*NPHIA*NLINE
        ENDIF
        DO JPHI=1,NPHI
          IF(ISPIN.EQ.0) THEN
            IF(I_SO.EQ.0) THEN
              IF(NS_E.EQ.0) THEN
                READ(1,7) JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1
              ELSE
                READ(1,8) JPLAN,DTHETA,DPHI,ECIN,
     1                    DTHETAA,DPHIA,SR1_1,SF1_1
              ENDIF
            ELSE
              IF(NS_E.EQ.0) THEN
                READ(1,27) JSO,JPLAN,DTHETA,DPHI,ECIN,SR1_1,SF1_1
              ELSE
                READ(1,28) JSO,JPLAN,DTHETA,DPHI,ECIN,
     1                     DTHETAA,DPHIA,SR1_1,SF1_1
              ENDIF
            ENDIF
          ELSE
            IF(NS_E.EQ.0) THEN
              READ(1,3) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,SR1_1,SF1_1
            ELSE
              READ(1,1) JSO,JPLAN,JS,DTHETA,DPHI,ECIN,
     1                  DTHETAA,DPHIA,SR1_1,SF1_1
            ENDIF
          ENDIF
          DO JSTEP=1,NSTEP-1
            READ(1,*) JPLAN
          ENDDO
          IF(NS_E.EQ.0) THEN
            PH(JPHI)=DPHI
          ELSEIF(NS_E.EQ.1) THEN
            PH(JPHI)=DPHIA
          ELSEIF(NS_E.EQ.2) THEN
            PH(JPHI)=DPHI
          ENDIF
        ENDDO
        EC(1)=ECIN
        IF(I_TH.EQ.0) TH(1)=DTHETA
      ENDIF
C
      REWIND 1
      DO JDUM=1,NHEAD2
        READ(1,30) NDUM
      ENDDO
      IF(NS_E.GT.0) THEN
        READ(1,30) NDUM
        READ(1,30) NDUM
      ENDIF
C
  1   FORMAT(2X,I3,2X,I2,2X,I1,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,
     1       F6.2,2X,E12.6,2X,E12.6)
  3   FORMAT(2X,I3,2X,I2,2X,I1,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,
     1       E12.6)
  7   FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
  8   FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,2X,E12.6,
     1       2X,E12.6)
 10   FORMAT(/,'<<<<<  NPOINTS IN SUBROUTINES FIXED_VAR AND ',
     1      'SELEC SHOULD BE AT LEAST ',I7,'  >>>>>',/,'  ')
 27   FORMAT(2X,I1,2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
 28   FORMAT(2X,I1,2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,F6.2,2X,
     1       E12.6,2X,E12.6)
 30   FORMAT(A1)
C
      RETURN
C
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE SELEC(NS_E,NE,NTHETA,NPHI,I_EC,I_TH,I_PH,
     1                 E_C,THETA,PHI)
C
      PARAMETER (NPOINTS=1000)
C
      COMMON /VALUES/ EC(NPOINTS),TH(NPOINTS),PH(NPOINTS)
C
      E_C=EC(1)
      THETA=TH(1)
      PHI=PH(1)
C
      IF(I_EC.EQ.1) THEN
        WRITE(6,*) '  '
        WRITE(6,16) NE
        WRITE(6,*) '  '
        DO JE=1,NE
          WRITE(6,4) EC(JE)
        ENDDO
        WRITE(6,*) '  '
        WRITE(6,*) '  '
        WRITE(6,*) 'WHICH ONE DO YOU WISH TO SELECT ?  '
        WRITE(6,*) '  '
        READ(5,24) E_C
      ENDIF
      IF(I_TH.EQ.1) THEN
        WRITE(6,*) '  '
        IF(NS_E.EQ.0) THEN
          WRITE(6,13) NTHETA
        ELSEIF(NS_E.EQ.1) THEN
          WRITE(6,14) NTHETA
        ELSEIF(NS_E.EQ.2) THEN
          WRITE(6,15) NTHETA
        ENDIF
        LIN=NTHETA/5
        NREM=MOD(NTHETA,5)
        WRITE(6,*) '  '
        DO JL=1,LIN
          JTHETA=(JL-1)*5+1
          WRITE(6,5) TH(JTHETA),TH(JTHETA+1),TH(JTHETA+2),
     1               TH(JTHETA+3),TH(JTHETA+4)
        ENDDO
        IF(NREM.GT.0) THEN
          JL=5*LIN+1
          WRITE(6,5) (TH(JTHETA), JTHETA=JL,NTHETA)
        ENDIF
        WRITE(6,*) '  '
        WRITE(6,*) '  '
        WRITE(6,*) 'WHICH ONE DO YOU WISH TO SELECT ?  '
        WRITE(6,*) '  '
        READ(5,25) THETA
      ENDIF
      IF(I_PH.EQ.1) THEN
        WRITE(6,*) '  '
        IF(NS_E.EQ.0) THEN
          WRITE(6,10) NPHI
        ELSEIF(NS_E.EQ.1) THEN
          WRITE(6,11) NPHI
        ELSEIF(NS_E.EQ.2) THEN
          WRITE(6,12) NPHI
        ENDIF
        LIN=NPHI/5
        NREM=MOD(NPHI,5)
        WRITE(6,*) '  '
        DO JL=1,LIN
          JPHI=(JL-1)*5+1
          WRITE(6,5) PH(JPHI),PH(JPHI+1),PH(JPHI+2),
     1               PH(JPHI+3),PH(JPHI+4)
        ENDDO
        IF(NREM.GT.0) THEN
          JL=5*LIN+1
          WRITE(6,5) (PH(JPHI), JPHI=JL,NPHI)
        ENDIF
        WRITE(6,*) '  '
        WRITE(6,*) '  '
        WRITE(6,*) 'WHICH ONE DO YOU WISH TO SELECT ?  '
        WRITE(6,*) '  '
        READ(5,25) PHI
      ENDIF
C
  4   FORMAT(4X,'*  ',F8.2)
  5   FORMAT(5('    *  ',F6.2),'   *')
 10   FORMAT(' THE ',I5,' POSSIBLE PHI ANGLES ARE : ')
 11   FORMAT(' THE ',I5,' POSSIBLE AUGER ELECTRON PHI ANGLES ',
     1       'ARE : ')
 12   FORMAT(' THE ',I5,' POSSIBLE PHOTOELECTRON PHI ANGLES ',
     1       'ARE : ')
 13   FORMAT(' THE ',I5,' POSSIBLE THETA ANGLES ARE : ')
 14   FORMAT(' THE ',I5,' POSSIBLE AUGER ELECTRON THETA ',
     1       'ANGLES ARE : ')
 15   FORMAT(' THE ',I5,' POSSIBLE PHOTOELECTRON THETA ',
     1       'ANGLES ARE : ')
 16   FORMAT(' THE ',I5,' POSSIBLE ENERGIES ARE : ')
 24   FORMAT(F8.2)
 25   FORMAT(F6.2)
C
      RETURN
C
      END

