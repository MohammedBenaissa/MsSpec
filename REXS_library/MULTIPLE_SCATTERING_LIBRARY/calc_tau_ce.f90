!
!======================================================================
!
MODULE CALC_TAU0J_CE
!
!  This module provides the subroutine to compute
!    Tau^{0 j} or Tau^{j 0} in the correlation expansion case
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE CALC_TAU_CE(TAU_TYPE,IATL,JE,TAU)
!
!
!   This subroutine calculates the scattering path operator by
!    the correlation expansion method.
!
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
!
!   The scattering path operator matrix of each small atom group
!    is obtained by using the LU decomposition method.
!
!   The running time of matrix inversion subroutine used in this program
!    scales with N^3, the memory occupied scales with N^2. We advise user to
!    use full MS method to get the scattering path operator, i.e. directly
!    with matrix inversion method if NGR is larger than 3. If NGR is less
!    than 4 (i.e <=3) this subroutine will gain time.
!
!   This subroutine never gains memory compared to the subrourine INV_MAT_MS
!    as I use three large matrices stored in common, each matrix is larger or
!    as large as the matrix used in INV_MAT_MS.
!
!   As I don't find a good way to solve the group problem, where all the contribution
!    of group IGR<=NGR are collected and each small contribution has to be stored
!    for the further larger-atom-group contribution, it's better that users change the
!    parameter NGR_M which is set in included file 'spec.inc' to be NGR or NGR+1
!    where NGR is the cut-off.user insterested. this subrouitne works for NGR is less
!    than 6(<=5), if users want to calculate larger NGR, they should modify the code here
!    to make them workable, the code is marked by 'C' in each lines (about 300 lines
!    below here), users just release them until to the desired cut-off, the maximum is
!    9, however, users can enlarge it if they want to. Warning ! NGR_M set in
!    included file should be larger than NGR and the figure listed below, don't forget
!    to compile the code after modification.
!
!   Users can modify the code to make it less memory-occupied, however, no matter they
!    do, the memories that used are more than full MS method used, so the only advantage
!    that this code has is to gain time when NGR<=3, with command 'common' used here,
!    the code will run faster.
!
!
!   Author :  H.-F. Zhao : 2007
!
!
!                                                      Last modified (DS) : 10 Jun 2021
!
!
      USE COMPLEX_NUMBERS,           ONLY : ZEROC,ONEC
!
      USE CLUSTER_COORD
      USE CLUSTER_LIMITS,            ONLY : VAL
      USE CURRENT_BEAM
      USE CURRENT_CALC
      USE CURRENT_T_MATRIX
      USE INIT_L
!
      USE ANGULAR_MOMENTUM,          ONLY : STORE_GAUNT
!
      USE CORR_EXP
      USE Q_ARRAY
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)      ::  TAU_TYPE
!
      INTEGER, INTENT(INOUT)    ::  IATL
      INTEGER, INTENT(IN)       ::  JE
!
      INTEGER                   ::  NTYPEM,JNEM
      INTEGER                   ::  LOGF
      INTEGER                   ::  NLM(NGR_M)
      INTEGER                   ::  ITYP(NGR_M)
      INTEGER                   ::  IGS(NGR_M)
      INTEGER                   ::  LM0,NRHS,NGR_MAX,NGR
      INTEGER                   ::  LJ,MJ,ILJ,INDJ
      INTEGER                   ::  INDK,NBTYPK
      INTEGER                   ::  IGR,LMJ,INDJM
      INTEGER                   ::  I,J
      INTEGER                   ::  KTYP,KNUM,KATL
      INTEGER                   ::  LMK,INDKM
!
      INTEGER                   ::  MIN
!
      COMPLEX (WP), INTENT(OUT) ::  TAU(LINMAX,LINFMAX,NATCLU_M)
!
      COMPLEX (WP)              ::  TLJ
      COMPLEX (WP)              ::  TAU1(LINMAX,LINFMAX,NATCLU_M)
!
      NTYPEM = 1                                                    ! number of absorbers
      JNEM   = 1                                                    ! type of current absorber
!
      LOGF = 6                                                      ! log file unit
!
      LM0  = LMAX(1,JE)                                             !
      LM0  = MIN(LM0,LF2)                                           !
      NRHS = (LM0 + 1) * (LM0 + 1)                                  !
!
      NGR_MAX = NGR_M                                               !
      NGR     = N_SCAT                                              !
!
!..........    Consistency checks    ..........
!
      IF(NGR_M > NATCLU) THEN                                       !
        WRITE(LOGF,10) NATCLU                                       !
        NGR_MAX = NATCLU                                            !
      END IF                                                        !
!
      IF(NGR < 1) THEN                                              !
        WRITE(LOGF,20)                                              !
        STOP                                                        !
      ELSE                                                          !
        IF(NGR > NGR_MAX) THEN                                      !
          WRITE(LOGF,30) NGR_MAX                                    !
          NGR = NGR_MAX                                             !
        END IF                                                      !
      END IF                                                        !
!
!..........    Storage of Gaunt coefficients    ..........
!
      CALL STORE_GAUNT                                              !
!
!  Case NGR = 1
!
      IF(NGR == 1) THEN                                             !
        DO LJ = 0, LM0                                              !
          ILJ = LJ * LJ + LJ + 1                                    !
          TLJ = TL(LJ,1,1,JE)                                       !
          DO MJ = - LJ, LJ                                          !
            INDJ = ILJ + MJ                                         !
            TAU(INDJ,INDJ,1) = TLJ                                  !
          END DO                                                    !
        END DO                                                      !
!
        GO TO 100                                                   !
      END IF                                                        !
!
!  NGR >=2 case
!
!
      DO INDJ = 1, NRHS                                             !
        TAU1(INDJ,INDJ,1) = Q(1) * ONEC                             !
      END DO                                                        !
!
!  Constructs the group matrix and inverts it
!
      IGR       = 1                                                 !
      LMJ       = LMAX(1,JE)                                        !
      NLM(IGR)  = LMJ                                               !
      INDJM     = (LMJ + 1) * (LMJ + 1)                             !
      ITYP(IGR) = 1                                                 !
      IGS(IGR)  = 1                                                 !
!
!  Initialization of matrix A to unit matrix
!
      DO I = 1, INDJM                                               !
        DO J = 1, INDJM                                             !
          IF (J == I) THEN                                          !
            A(J,I) = ONEC                                           !
          ELSE                                                      !
            A(J,I) = ZEROC                                          !
          END IF                                                    !
        END DO                                                      !
      END DO                                                        !
!
!  Call of recursive subroutine COREXP_SAVM to compute (I - G_o T)^{-1}
!
      IGR = IGR + 1                                                 !
      CALL COREXP_SAVM(JE,IGR,NGR,NLM,ITYP,IGS,TAU1)                !
      IGR = IGR - 1                                                 !
!
!     TAU = TAU * tj
!
      DO KTYP = 1, N_PROT                                           !
        NBTYPK = NATYP(KTYP)                                        !
        LMK    = LMAX(KTYP,JE)                                      !
        INDKM  = (LMK + 1) * (LMK + 1)                              !
        DO KNUM = 1, NBTYPK                                         !
          KATL = NCORR(KNUM,KTYP)                                   !
!
          DO LJ = 0, LM0                                            !
            ILJ = LJ * LJ + LJ + 1                                  !
            TLJ = TL(LJ,1,1,JE)                                     !
            DO MJ = -LJ, LJ                                         !
              INDJ = ILJ + MJ                                       !
!
              DO INDK = 1, INDKM                                    !
                TAU(INDK,INDJ,KATL) = TAU1(INDK,INDJ,KATL) * TLJ    !
              END DO                                                !
!
            END DO                                                  !
          END DO                                                    !
!
        END DO                                                      !
      END DO                                                        !
!
 100  CONTINUE                                                      !
!
!  Formats:
!
  10  FORMAT(//,5X,'     --->  NGR_M should be smaller than NATCLU',/, &
                5X,'     --->  it is reduced to NATCLU = ',I6)
  20  FORMAT(//,5X,'     --->  NGR < 1, no expansion is done')
  30  FORMAT(//,5X,'     --->  NGR is too large, reduce to NGR_M = ',I6)
!
      END SUBROUTINE CALC_TAU_CE
!
!======================================================================
!
      RECURSIVE SUBROUTINE COREXP_SAVM(JE,IGR,NGR,NLM,ITYPE,IGS,TAU)
!
!  This subroutine call the correlation matrices calculations
!    for a given order IGR
!
!
!  Input variables:
!
!     JE       : energy point
!     IGR      : index of the correlation matrix
!     NGR      : maximum size of correlation matrix (?)
!
!  Output variable:
!
!     NLM     :
!     ITYP    :
!     IGS     :
!     TAU     : scattering path operator
!
!
!   Author :  H.-F. Zhao : 2007
!
!
!                                    Last modified (DS) :  4 Jun 2021
!
!
      USE CLUSTER_COORD
      USE CURRENT_T_MATRIX
!
      USE Q_ARRAY
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  JE,NGR
      INTEGER, INTENT(OUT)      ::  NLM(NGR_M)
      INTEGER, INTENT(OUT)      ::  ITYPE(NGR_M)
      INTEGER, INTENT(OUT)      ::  IGS(NGR_M)
      INTEGER, INTENT(INOUT)    ::  IGR
!
      INTEGER                   ::  ITYP,NBTYP,NUM
!
      REAL (WP)                 ::  QI
!
      COMPLEX (WP), INTENT(OUT) ::  TAU(LINMAX,LINFMAX,NATCLU_M)
!
      DO ITYP = 1, N_PROT                                           !
!
        NBTYP      = NATYP(ITYP)                                    !
        NLM(IGR)   = LMAX(ITYP,JE)                                  !
        ITYPE(IGR) = ITYP                                           !
!
        DO NUM = 1, NBTYP                                           !
!
          IGS(IGR) = NCORR(NUM,ITYP)                                !
!
          IF(IGS(IGR) > IGS(IGR-1)) THEN                            !
            QI = Q(IGR)                                             !
!
            CALL MPIS(IGR,NLM,ITYPE,IGS,JE,QI,TAU)                  !
!
            IGR = IGR + 1                                           !
!
!  Recursive call of the subroutine
!
            IF(IGR <= NGR) THEN                                     !
              CALL COREXP_SAVM(JE,IGR,NGR,NLM,ITYPE,IGS,TAU)        !
            END IF                                                  !
!
            IGR = IGR - 1                                           !
!
          END IF                                                    !
!
        END DO                                                      !
!
      END DO                                                        !
!
      END SUBROUTINE COREXP_SAVM
!
!======================================================================
!
      SUBROUTINE MPIS(N,NLM,ITYP,IGS,JE,QI,TAU)
!
!
!   This subroutine constructs the correlation matrices and uses
!     the LU decomposition method to perform the matrix inversion.
!
!   The inverse matrix which is the contribution of a small group of atoms
!     is kept for further use.
!
!
!  Input variables:
!
!     N       : index of the correlation matrix
!     NLM     :
!     ITYP    :
!     IGS     :
!     JE      : energy point
!     QI      :
!
!  Output variable:
!
!     TAU     : scattering path operator
!
!
!   Author :  H.-F. Zhao : 2007
!
!
!                                    Last modified (DS) : 16 Jun 2021
!
!
      USE COMPLEX_NUMBERS
!
      USE CLUSTER_COORD
      USE CURRENT_T_MATRIX
      USE INIT_L
!
      USE SPHERICAL_BESSEL
      USE SPHERICAL_HARMONICS
!
      USE CORR_EXP
      USE GAUNT_C
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  N,JE
      INTEGER, INTENT(IN)       ::  NLM(NGR_M)
      INTEGER, INTENT(IN)       ::  ITYP(NGR_M)
      INTEGER, INTENT(IN)       ::  IGS(NGR_M)
!
      INTEGER                   ::  LOGF
      INTEGER                   ::  ONE_L
      INTEGER                   ::  IBESS
      INTEGER                   ::  LM0,NRHS,INDJ,NM
      INTEGER                   ::  I,J
      INTEGER                   ::  L,LNMAX,NM1,NML,NTYP
      INTEGER                   ::  JATL,JTYP,LJMAX,J1
      INTEGER                   ::  II,JJ,IN1
      INTEGER                   ::  LN,ILN,MLN,JJ0,JJ1,JJ1N
      INTEGER                   ::  INDN
      INTEGER                   ::  LJ,ILJ,L_MIN,L_MAX
      INTEGER                   ::  MLJ,INDJN,MA,MB
      INTEGER                   ::  INFO,INFO1
      INTEGER                   ::  IPIV(NLMM)
      INTEGER                   ::  KLIN,K,KATL,LMK
      INTEGER                   ::  INDK,INDKM
!
      INTEGER                   ::  MIN,ABS
!
      REAL (WP), INTENT(IN)     ::  QI
!
      REAL (WP)                 ::  PI4
      REAL (WP)                 ::  XJN,YJN,ZJN,RJN
      REAL (WP)                 ::  KRJN,ZDJN
      REAL (WP)                 ::  IM_VK,RE_VK
      REAL (WP)                 ::  XN,YN,ZN
!
      REAL (WP)                 ::  AIMAG,REAL,SQRT,EXP
!
      COMPLEX (WP), INTENT(OUT) ::  TAU(LINMAX,LINFMAX,NATCLU_M)
!
      COMPLEX (WP)              ::  ATTL(0:NT_M,NATM)
      COMPLEX (WP)              ::  EXPJN,ATTJN
      COMPLEX (WP)              ::  CN1,IC_REF
      COMPLEX (WP)              ::  PI4_IC,IC_L
      COMPLEX (WP)              ::  YLM(0:NLTWO,-NLTWO:NLTWO)
      COMPLEX (WP)              ::  HL1(0:NLTWO)
      COMPLEX (WP)              ::  SUM_L,SUM_L2
      COMPLEX (WP)              ::  TEMP,TEMP1,TEMP2
      COMPLEX (WP)              ::  AINV(NLMM,NLMM)
      COMPLEX (WP)              ::  IN(NLMM,LINFMAX)
      COMPLEX (WP)              ::  SUM_L_A,SUM_L2_A
      COMPLEX (WP)              ::  SUM_L_B,SUM_L2_B
!
      DATA PI4/ 12.5663706143591729538505735331180E0_WP /           ! 4 pi
!
      LOGF = 6                                                      ! log file unit
!
      IBESS  = 3                                                    !
      PI4_IC = - IC * PI4                                           !
!
      LM0  = LMAX(1,JE)                                             !
      LM0  = MIN(LM0,LF2)                                           !
      NRHS = (LM0 + 1) * (LM0 + 1)                                  !
      INDJ = 0               !                                      !
!
      NM = 0                                                        !
      DO I = 1, N-1                                                 !
        J  = NLM(I) + 1                                             !
        NM = NM + J * J                                             !
      END DO                                                        !
!
      L     = NLM(N)                                                !
      LNMAX = L                                                     !
      L     = (L + 1) * (L + 1)                                     !
      NM1   = NM + 1                                                !
      NML   = NM + L                                                !
      NTYP  = ITYP(N)                                               !
!
      DO L = 0, LNMAX                                               !
        ATTL(L,N) = TL(L,1,NTYP,JE)                                 !
      END DO                                                        !
!
      IM_VK = - AIMAG(VK(JE))                                       !
      RE_VK =   REAL(VK(JE),KIND=WP)                                !
!
!  Setting up matrix blocks C((N-1)*1) and D(1*(N-1))
!
      I  = IGS(N)                                                   !
      XN = SYM_AT(1,I)                                              !
      YN = SYM_AT(2,I)                                              !
      ZN = SYM_AT(3,I)                                              !
!
      DO J = 1, N-1                                                 !
        JATL  = IGS(J)                                              !
        LJMAX = NLM(J)                                              !
        JTYP  = ITYP(J)                                             !
        J1    = J - 1                                               !
!
        XJN = (SYM_AT(1,JATL) - XN)                                 !
        YJN = (SYM_AT(2,JATL) - YN)                                 !
        ZJN = (SYM_AT(3,JATL) - ZN)                                 !
        RJN = SQRT(XJN * XJN + YJN * YJN + ZJN * ZJN)               !
!
        KRJN  = RE_VK * RJN                                         !
        ATTJN = PI4_IC * EXP(IM_VK * RJN)                           !
        EXPJN = (XJN + IC * YJN) / RJN                              !
        ZDJN  = ZJN / RJN                                           !
!
        CALL SPH_HAR(2*NL_M,ZDJN,EXPJN,YLM,LNMAX+LJMAX)             !
        CALL BESPHE(LNMAX+LJMAX,IBESS,KRJN,HL1)                     !
!
        DO L = 0, LJMAX                                             !
          ATTL(L,J) = ATTJN * TL(L,1,JTYP,JE)                       !
        END DO                                                      !
!
        II  = NM                                                    !
        IN1 = - 1                                                   !
        CN1 = IC                                                    !
        JJ  = 0                                                     !
!
        DO LN = 0, LNMAX                                            !
          ILN = LN * LN + LN + 1                                    !
          IN1 = - IN1                                               ! (-1)^LN
          CN1 = - CN1 * IC                                          ! (-i)^LN
!
          DO MLN = -LN, LN                                          !
            INDN   = ILN + MLN                                      !
            II     = II + 1                                         !
            JJ0    = J1 * INDJ                                      !
            ONE_L  = - IN1                                          !
            IC_REF = - CN1 * IC                                     !
!
            DO LJ = 0, LJMAX                                        !
              ILJ    = LJ * LJ + LJ + 1                             !
              L_MIN  = ABS(LJ - LN)                                 !
              L_MAX  = LJ + LN                                      !
              ONE_L  = - ONE_L                                      !
              IC_REF = IC_REF * IC                                  !
!
!  Case MLJ equal to zero
!
              JJ1 = JJ0 + ILJ                                       !
              IF(LJ >= LN) THEN                                     !
                IC_L = - IC_REF                                     !
              ELSE                                                  !
                IC_L = - ONEC / IC_REF                              !
              END IF                                                !
!
              SUM_L  = ZEROC                                        !
              SUM_L2 = ZEROC                                        !
!
              DO L = L_MIN, L_MAX, 2                                !
                IC_L = - IC_L                                       !
                IF(ABS(MLN) <= L) THEN                              !
                  TEMP   = IC_L * HL1(L) * GNT(L,ILJ,INDN)          !
                  SUM_L  = SUM_L + YLM(L,MLN) * TEMP                !
                  SUM_L2 = SUM_L2 + CONJG(YLM(L,MLN)) * TEMP        !
                END IF                                              !
              END DO                                                !
!
              IF(ONE_L == -1) SUM_L2 = - SUM_L2                     !
!
              A(JJ1,II) = ATTL(LJ,J) * SUM_L                        !
              A(II,JJ1) = ATTJN * ATTL(LN,N) * SUM_L2               !
!
!  Case MLJ not equal to zero
!
              DO MLJ = 1, LJ                                        !
                INDJ  = ILJ + MLJ                                   !
                INDJN = ILJ - MLJ                                   !
                JJ1   = JJ0 + INDJ                                  !
                JJ1N  = JJ0 + INDJN                                 !
                MA    = MLN - MLJ                                   !
                MB    = MLN + MLJ                                   !
                IF(LJ >= LN) THEN                                   !
                  IC_L = - IC_REF                                   !
                ELSE                                                !
                  IC_L = - ONEC / IC_REF                            !
                END IF                                              !
!
                SUM_L_A  = ZEROC                                    !
                SUM_L2_A = ZEROC                                    !
                SUM_L_B  = ZEROC                                    !
                SUM_L2_B = ZEROC                                    !
!
                DO L = L_MIN ,L_MAX,2                               !
                  IC_L = - IC_L                                     !
                  IF(ABS(MA) <= L) THEN                             !
                    TEMP1    = IC_L * HL1(L) * GNT(L,INDJ,INDN)     !
                    SUM_L_A  = SUM_L_A + YLM(L,MA) * TEMP1          !
                    SUM_L2_A = SUM_L2_A + CONJG(YLM(L,MA)) * TEMP1  !
                  END IF                                            !
                  IF(ABS(MB) <= L) THEN                             !
                    TEMP2    = IC_L * HL1(L) * GNT(L,INDJN,INDN)    !
                    SUM_L_B  = SUM_L_B + YLM(L,MB) * TEMP2          !
                    SUM_L2_B = SUM_L2_B + CONJG(YLM(L,MB)) * TEMP2  !
                  END IF                                            !
                END DO                                              !
!
                IF(ONE_L == -1) THEN                                !
                   SUM_L2_A = - SUM_L2_A                            !
                   SUM_L2_B = - SUM_L2_B                            !
                END IF                                              !
!
                A(JJ1,II)  = ATTL(LJ,J) * SUM_L_A                   !
                A(II,JJ1)  = ATTJN * ATTL(LN,N) * SUM_L2_A          !
                A(JJ1N,II) = ATTL(LJ,J) * SUM_L_B                   !
                A(II,JJ1N) = ATTJN * ATTL(LN,N) * SUM_L2_B          !
              END DO                                                ! end of MLJ loop
!
            END DO                                                  ! end of LJ loop
!
            JJ = JJ0 + INDJ                                         !
!
          END DO                                                    ! end of MLN loop
!
        END DO                                                      ! end of LN loop
!
        JJ = JJ - INDN                                              !
!
      END DO                                                        ! end of J loop
!
!  Adding B to A
!
      DO I = NM1, NML                                               !
        DO J = NM1, NML                                             !
          IF(J == I) THEN                                           !
            A(J,I) = ONEC                                           !
          ELSE                                                      !
            A(J,I) = ZEROC                                          !
          END IF                                                    !
        END DO                                                      !
      END DO                                                        !
!
!  Constructing AINV matrix
!
      DO I = 1, NML                                                 !
        DO J = 1, NML                                               !
          AINV(J,I) = A(J,I)                                        !
        END DO                                                      !
      END DO                                                        !
!
!  Matrix inversion using A * X = B
!
!    1) Decomposition into L * U form:
!
      CALL ZGETRF(NML,NML,AINV,NLMM,IPIV,INFO1)                     !
!
      IF(INFO1 /= 0) THEN                                           !
        WRITE(LOGF,10) INFO1                                        ! error in the LU factorization
      ELSE                                                          !
!
        DO I = 1 ,NRHS                                              !
          DO J = 1, NML                                             !
            IF(J == I) THEN                                         !
              IN(J,I) = ONEC                                        !
            ELSE                                                    !
              IN(J,I) = ZEROC                                       !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
!
!    2) Solving A * X = B with A in L*U form
!
        CALL ZGETRS('N',NML,NRHS,AINV,NLMM,IPIV,IN,NLMM,INFO)       !
        IF(INFO /= 0) THEN                                          !
          WRITE(LOGF,20) INFO                                       ! error in the inversion
        END IF                                                      !
      END IF                                                        !
!
!  Contribution of X to the scattering path operator TAU
!
      KLIN = 0                                                      !
      DO K = 1, N                                                   !
        KATL  = IGS(K)                                              !
        LMK   = NLM(K)                                              !
        INDKM = (LMK + 1) * (LMK + 1)                               !
!
        DO INDJ = 1 ,NRHS                                           !
!
          DO INDK = 1, INDKM                                        !
            KLIN = KLIN + 1                                         !
!
            TAU(INDK,INDJ,KATL) = TAU(INDK,INDJ,KATL) +           & !
                                  QI * IN(KLIN,INDJ)                !
!
          END DO                                                    !
!
          KLIN = KLIN - INDKM                                       !
!
        END DO                                                      !
!
        KLIN = KLIN + INDKM                                         !
!
      END DO
!
!  Formats:
!
  10  FORMAT(//,5X,'    --->  INFO  = ',I2)
  20  FORMAT(//,5X,'    --->  INFO1 = ',I2)
!
      END SUBROUTINE MPIS
!
END MODULE CALC_TAU0J_CE
