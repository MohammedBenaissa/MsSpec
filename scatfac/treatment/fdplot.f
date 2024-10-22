C
C
C  ************************************************************************
C  *                                                                      *
C  *    TREATMENT OF THE SCATTERING FACTOR CALCULATED DATA FOR PLOTTING   *
C  *                                                                      *
C  ************************************************************************
C
C
C
C     This program reads the result file of the scattering factor 
C      computed by the MsSpec code and allows to select
C      interactively in this file what is going to be plotted by the user.
C
C  Author : D. Sebilleau
C
C                                              Last modified : 13 Nov 2007
C
      PROGRAM MAIN
C
      PARAMETER(NDIM=100000)
C
      CHARACTER*18 TITRE1,TITRE2
      CHARACTER*1 REP1,REP2,REP3
C
      REAL ANGLE(2,NDIM),RAD(2,NDIM),E(NDIM),FD(4,NDIM)
C
      INTEGER MAT(4,NDIM)
C
      COMPLEX FSPH
C
      DATA PI/ 3.141593/
C
      WRITE(6,20)
C
      WRITE(6,*) '  '
  13  WRITE(6,*) '  '
      PRINT *, '   NAME OF THE INPUT FILE (18 CHARACTERS ',
     1         'MAXIMUM) :'
      WRITE(6,*) '  '
      READ(5,1) TITRE1
      OPEN (UNIT=1, FILE=TITRE1, STATUS='OLD')
C
  21  WRITE(6,*) '  '
      PRINT *, '   NAME OF THE OUTPUT FILE (18 CHARACTERS ',
     1         'MAXIMUM) :'
      WRITE(6,*) '  '
      READ(5,1) TITRE2
      WRITE(6,*) '  '
      OPEN (UNIT=7, FILE=TITRE2, STATUS='UNKNOWN')
C
      READ(1,2) ISPHER,NL,NAT,L,NTHT,NPHI,NE,E0,EFIN
      IF(ISPHER.LE.1) THEN
        LMAX=0
      ELSE
        LMAX=L
      ENDIF
 
      PRINT *, '     '
      PRINT *, '   GIVE THE NUMBER CORRESPONDING TO WHAT ',
     1         'YOU WANT TO PLOT :'
      PRINT *, '         '
      PRINT *, '       1 = MODULUS OF THE SCATTERING FACTOR '
      PRINT *, '       2 = PHASE OF THE SCATTERING FACTOR'
      PRINT *, '       3 = REAL PART OF THE SCATTERING FACTOR '
      PRINT *, '       4 = IMAGINARY PART OF THE SCATTERING FACTOR '
      WRITE(6,*) '  '
      READ(5,*) NR1
C
      PRINT *, '   '
      PRINT *, '   WITH RESPECT TO WHICH PARAMETER(S) ?'
      PRINT *, '                    '
      PRINT *, '       1 = THETA AND PHI '
      PRINT *, '       2 = THETA '
      PRINT *, '       3 = PHI '
      PRINT *, '       4 = ENERGY '
      WRITE(6,*) '  '
      READ(5,*) NR2
      IF(NR2.EQ.4) GOTO 16
C
      PRINT *, '     '
      PRINT *, '   HOW WILL YOU PLOT THE FILE ?'
      PRINT *, '      '
      PRINT *, '       1 = POLAR REPRESENTATION (R,THETA,PHI)'
      PRINT *, '       2=  CARTESIAN REPRESENTATION (X,Y,Z)'
      WRITE(6,*) '  '
      READ(5,*) NR3
      GOTO 17
  16  NR3=1
C
  17  PRINT *, '     '
      PRINT *, '   NUMBER OF THE ATOM : '
      WRITE(6,*) '  '
      READ(5,*) JAT
      IF(JAT.GT.NAT) GOTO 3
C
      IF(NR2.EQ.4) GOTO 4
C
   6  IF(NE.GT.1) THEN
        WRITE(6,*) '  '
        PRINT *,'   NUMBER OF THE ENERGY POINT : '
        WRITE(6,*) '  '
        READ(5,*) JE
      ELSE
        JE=1
      ENDIF
      IF(JE.GT.NE) GO TO 3
      IF(NE.EQ.1) GO TO 4
      ECIN=E0+FLOAT(JE-1)*(EFIN-E0)/FLOAT(NE)
      WRITE(6,12) ECIN
C
      PRINT *, '   IS IT OK ? (Y/N) '
      WRITE(6,*) '  '
      READ(5,5) REP1
      IF((REP1.NE.'Y').AND.(REP1.NE.'y')) GO TO 6
C
   4  IF(ISPHER.LE.1) GOTO 15
C
      WRITE(6,*) '  '
      PRINT *, '   VALUE OF M (AZIMUTHAL QUANTUM NUMBER) :'
      WRITE(6,*) '  '
      READ(5,*) M
      IF(ABS(M).GT.L) GO TO 3
C
 15   NPOINT=NE*NAT*(2*LMAX+1)*NTHT*NPHI
      IF(NPOINT.GT.NDIM) THEN
        WRITE(6,9) TITRE1
        STOP
      ENDIF
C
      DO 7 J=1,NPOINT
        READ(1,8) (MAT(JCOL,J),JCOL=1,4),(FD(JCOL,J),JCOL=3,4),
     1            (ANGLE(JCOL,J),JCOL=1,2),E(J)
        FD(1,J)=SQRT(FD(3,J)**2+FD(4,J)**2)
        FSPH=FD(3,J)+(0.,1.)*FD(4,J)
        CALL ARCSIN(FSPH,PHASE)
        FD(2,J)=PHASE
        RAD(1,J)=ANGLE(1,J)*PI/180.
        RAD(2,J)=ANGLE(2,J)*PI/180.
        IF((MAT(1,J).EQ.JE).OR.(NR2.EQ.4)) THEN
          IF(MAT(2,J).EQ.JAT) THEN
            IF((MAT(3,J).EQ.L).OR.(ISPHER.EQ.0)) THEN
              IF((MAT(4,J).EQ.M).OR.(ISPHER.EQ.0)) THEN
                IF(NR2.EQ.1) THEN
                  IF(NR3.EQ.1) THEN
                    X=FD(NR1,J)
                    Y=ANGLE(1,J)
                    Z=ANGLE(2,J)
                  ELSE
                    X=FD(NR1,J)*SIN(RAD(1,J))*COS(RAD(2,J))
                    Y=FD(NR1,J)*SIN(RAD(1,J))*SIN(RAD(2,J))
                    Z=FD(NR1,J)*COS(RAD(1,J))
                  ENDIF
                ELSE
                  IF(NR3.EQ.1) THEN
                    Y=FD(NR1,J)
                    IF(NR2.EQ.4) THEN
                      X=E(J)
                    ELSE
                      X=ANGLE(NR2-1,J)
                    ENDIF
                  ELSE
                    X=FD(NR1,J)*COS(RAD(NR2-1,J))
                    Y=FD(NR1,J)*SIN(RAD(NR2-1,J))
                    Z=0.
                  ENDIF
                ENDIF
                WRITE(7,10) X,Y,Z
              ENDIF
            ENDIF
          ENDIF
        ENDIF
   7  CONTINUE
      CLOSE(7)
C
      WRITE(6,*) '  '
      PRINT *, '   WOULD YOU LIKE TO USE THIS INPUT FILE '
      PRINT *, '   FOR ANOTHER PLOT ? (Y/N)'
      WRITE(6,*) '  '
      READ(5,5) REP3
      IF((REP3.EQ.'y').OR.(REP3.EQ.'Y')) THEN
        REWIND 1
        GOTO 21
      ENDIF
C
      WRITE(6,*) '  '
      PRINT *, '   WOULD YOU LIKE TO PLOT ANOTHER FILE ? (Y/N)'
      WRITE(6,*) '  '
      READ(5,5) REP2
      IF((REP2.EQ.'y').OR.(REP2.EQ.'Y')) THEN
        NREP=1
        REWIND 1
      ELSE
        NREP=2
      ENDIF
      GO TO (13,11) NREP
C
  11  WRITE(6,*) '  '
      PRINT *, '   END OF THE TREATMENT'
      WRITE(6,*) '  '
      CLOSE(1)
      GO TO 14
C
   3  WRITE(6,*) '  '
      PRINT *,'      --->  ERROR IN THE DATA'
      WRITE(6,*) '  '
C
   1  FORMAT(A18)
   2  FORMAT(5X,I1,2X,I2,2X,I4,2X,I2,2X,I3,2X,I3,2X,I3,2X,F8.2,
     1       2X,F8.2)
   5  FORMAT(A1)      
   8  FORMAT(1X,I3,1X,I4,1X,I2,1X,I3,1X,F6.3,1X,F6.3,1X,
     1         F6.2,1X,F6.2,1X,F8.2)
   9    FORMAT(//,15X,'   <<<<< THERE ARE TOO MANY LINES ',
     1         'IN THE FILE ',A15,'  >>>>>',//)
  10  FORMAT(3(5X,E12.6))
  12  FORMAT('   THE ENERGY YOU ARE REQUESTING IS ',F8.2,' eV')
  20  FORMAT(' ',/,
     1       '    ===============================================',/,
     2       '                                                   ',/,
     3       '              PLOT OF A SCATTERING FACTOR          ',/,
     5       '                                                   ',/,
     6       '    ===============================================',/,
     7       ' ')
C
  14  END
C
C
C
      SUBROUTINE ARCSIN(U,RANGLE)
C
      COMPLEX U,CANGLE
      IF(CABS(U).LT.0.0001) THEN
        RANGLE=0.
      ELSE
        CANGLE=(0.,-1.)*CLOG(U/CABS(U))
        RANGLE=REAL(CANGLE)
      ENDIF
C
      RETURN
C
      END
