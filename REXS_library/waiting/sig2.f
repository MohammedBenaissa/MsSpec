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

