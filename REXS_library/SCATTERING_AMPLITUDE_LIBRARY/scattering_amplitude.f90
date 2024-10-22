!
!=======================================================================
!
MODULE SCATTERING_AMPLITUDE
!
!  This module provide the routine to compute the scattering amplitude
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE CALC_SF(J_EL)
!
!  This routine calculates the scattering factor and prepares the output
!      for the plot
!
!
!  Input variables :
!
!                       J_EL      :  electron index
!                                       J_EL = 1  : incoming electron
!                                       J_EL = 2  : excited  electron
!                                       J_EL = 3  : outgoing electron 1
!                                       J_EL = 4  : outgoing electron 2
!                                       J_EL = 5  : outgoing electron 3
!
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 28 May 2021
!
!
      USE REAL_NUMBERS,              ONLY : ZERO,ONE
      USE PI_ETC,                    ONLY : PI
!
      USE CLUSTER
      USE CURRENT_BEAM
      USE CURRENT_CALC
      USE CURRENT_INIT_VAL
      USE CURRENT_FINA_VAL
      USE CURRENT_F_TH
      USE CURRENT_T_MATRIX
      USE INIT_L
      USE OUTFILES
      USE OUTUNITS
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 27)  ::  SCATFILE
      CHARACTER (LEN =  3)  ::  EXT(5)
!
      INTEGER, INTENT(IN)   ::  J_EL
!
      INTEGER, PARAMETER    ::  MAX_LENGTH = 50
!
      INTEGER               ::  N_CHAR,N_DOT,J_CHAR
      INTEGER               ::  I_BASIS_OLD,NFTHET_OLD
      INTEGER               ::  ITHETA_OLD,NTHETA_OLD
      INTEGER               ::  IPHI_OLD,NPHI_OLD
      INTEGER               ::  IE_OLD,NE_OLD
      INTEGER               ::  L_MAX,L_OLD
      INTEGER               ::  JT,JP,JE,JAT,M
      INTEGER               ::  NPHIM,NTHT
!
      REAL (WP)             ::  SSMALL
      REAL (WP)             ::  THETA0_OLD,THETA1_OLD
      REAL (WP)             ::  PHI0_OLD,PHI1_OLD
      REAL (WP)             ::  E0_OLD,E1_OLD
      REAL (WP)             ::  PHITOT,THTOT
      REAL (WP)             ::  DTHETA,RTHETA,TEST
      REAL (WP)             ::  DPHI,RPHI
      REAL (WP)             ::  POZ,EPS,BETA,GAMMA
      REAL (WP)             ::  EKIN
      REAL (WP)             ::  REFTH,XIMFTH
!
      REAL (WP)             ::  FLOAT,SIN,ABS,MAX,REAL,AIMAG
!
      COMPLEX (WP)          ::  FSPH
!
      DATA  SSMALL  / 0.0001E0_WP /
      DATA  EXT     / '_in', '_ex', '_o1', '_o2', 'o_3' /
!
!  Constructing the name of the output file
!
!  1) length of the filename : N_CHAR
!
      N_CHAR = 0                                                    !
      DO J_CHAR = 1, MAX_LENGTH                                     !
        IF(OUTFILE3(J_CHAR:J_CHAR) == ' ') GO TO 400                !
        N_CHAR = N_CHAR + 1                                         !
      END DO                                                        !
 400  CONTINUE                                                      !
!
!  2) finding the position of the dot : N_DOT + 1
!
      N_DOT = 0                                                     !
      DO J_CHAR = 1, MAX_LENGTH                                     !
        IF(OUTFILE3(J_CHAR:J_CHAR) == '.') GO TO 500                !
        N_DOT = N_DOT + 1                                           !
      END DO                                                        !
 500  CONTINUE                                                      !
!
      SCATFILE = OUTFILE3(1:N_DOT)//EXT(J_EL)//OUTFILE3(N_DOT+1:N_CHAR)
!
!  Opening the scattering factor file
!
      OPEN(UNIT=IUO3, FILE=SCATFILE, STATUS='UNKNOWN')              !
!
!  Setting up the standard values (IFTHET = 0 case)
!
      IF(IFTHET == 0) THEN                                          !
!
        I_BASIS_OLD = I_BASIS                                       !
        NFTHET_OLD  = NFTHET                                        !
!
        I_BASIS = 0                                                 !
!
        ITHETA_OLD = ITHETA                                         !
        NTHETA_OLD = NTHETA                                         !
        IPHI_OLD   = IPHI                                           !
        NPHI_OLD   = NPHI                                           !
        IE_OLD     = IE                                             !
        NE_OLD     = NE                                             !
        THETA0_OLD = THETA0                                         !
        THETA1_OLD = THETA1                                         !
        PHI0_OLD   = PHI0                                           !
        PHI1_OLD   = PHI1                                           !
        E0_OLD     = E0                                             !
        E1_OLD     = E1                                             !
        L_OLD      = LI                                             !
!
        IF(IE == 0) THEN                                            !
          ITHETA = 1                                                !
          NTHETA = 1                                                !
          IPHI   = 0                                                !
        ELSE                                                        !
          ITHETA = 0                                                !
          NTHETA = 1                                                !
          IPHI   = 0                                                !
        END IF                                                      !
!
        THETA1 = ZERO                                               !
        PHI1   = ZERO                                               !
!
      END IF                                                        !
!
      IF(I_BASIS == 0) THEN                                         !
         LI   = 0                                                   !
         L_MAX = 0                                                  !
      ELSE                                                          !
        L_MAX  = LI                                                 !
      END IF                                                        !
!
      PHITOT = 360.0E0_WP                                           !
      THTOT  = 360.0E0_WP * ITHETA * (1 - IPHI) +                 & !
               180.0E0_WP * ITHETA * IPHI                           !
      NPHI   = (NFTHET + 1) * IPHI + (1 - IPHI)                     !
      NTHT   = (NFTHET + 1) * ITHETA * (1 - IPHI) +               & !
               (NFTHET / 2 + 1) * ITHETA * IPHI   +               & !
               (1 - ITHETA)                                         !
      NE     = NFTHET * IE + (1 - IE)                               !
!
      WRITE(IUO3,10) I_BASIS,L_MAX,NAT,LI,NTHT,NPHI,NE,E0,E1        !
!
!  Loop over the theta angle
!
      DO JT = 1 , NTHT                                              !
!
        DTHETA = THETA1 + FLOAT(JT - 1) * THTOT /                 & !
                          FLOAT(MAX(NTHT-1,1))                      !
        RTHETA = DTHETA * PI / 180.0E0_WP                           !
        TEST   = SIN(RTHETA)                                        !
!
        IF(TEST >= ZERO) THEN                                       !
          POZ = PI                                                  !
          EPS = ONE                                                 !
        ELSE                                                        !
          POZ = ZERO                                                !
          EPS = - ONE                                               !
        END IF                                                      !
!
        BETA = RTHETA * EPS                                         !
        IF(ABS(TEST) < SSMALL) THEN                                 !
          NPHIM = 1                                                 !
        ELSE                                                        !
          NPHIM = NPHI                                              !
        END IF                                                      !
!
!  Loop over the phi angle
!
        DO JP = 1, NPHIM                                            !
!
          DPHI  = PHI1 + FLOAT(JP - 1) * PHITOT /                 & !
                         FLOAT(MAX(NPHI-1,1))                       !
          RPHI  = DPHI * PI / 180.0E0_WP                            !
          GAMMA = POZ - RPHI                                        !
!
!  Loop over the energies
!
          DO JE = 1, NE                                             !
!
            IF(NE == 1) THEN                                        !
              EKIN = E0                                             !
            ELSE                                                    !
              EKIN = E0 + FLOAT(JE - 1) * (E1 - E0) / FLOAT(NE - 1) !
            END IF                                                  !
!
!  Loop over the prototypical atoms
!
            DO JAT = 1, NAT                                         !
!
              IF(LI > LMAX(JAT,JE)) GO TO 200                       !
!
              DO M = - L_MAX, L_MAX                                 !
!
                CALL SW_SF(R0,R1,THETA0,PHI0,BETA,GAMMA,LI,M,     & !
                           FSPH,JAT,JE,*201)                        !
!
                REFTH  = REAL(FSPH,KIND=WP)                         !
                XIMFTH = AIMAG(FSPH)                                !
                WRITE(IUO3,20) JE,JAT,LI,M,REFTH,XIMFTH,          & !
                               DTHETA,DPHI,EKIN                     !
              END DO                                                !
!
!  End of loop on prototypical atoms
!
            END DO                                                  !
!
!  End of loop over energies
!
          END DO                                                    !
!
!  End of loop on phi angle
!
        END DO                                                      !
!
!  End of loop on theta angle
!
      END DO                                                        !
!
      CLOSE(IUO3)                                                   !
!
      WRITE(IUO1,300) EXT(J_EL)                                     !
      WRITE(IUO1,301) SCATFILE                                      !
!
!  Recovering the original parameters
!
      IF(IFTHET == 0) THEN                                          !
        I_BASIS = I_BASIS_OLD                                       !
        NFTHET  = NFTHET_OLD                                        !
        ITHETA  = ITHETA_OLD                                        !
        NTHETA  = NTHETA_OLD                                        !
        IPHI    = IPHI_OLD                                          !
        NPHI    = NPHI_OLD                                          !
        IE      = IE_OLD                                            !
        NE      = NE_OLD                                            !
        THETA0  = THETA0_OLD                                        !
        THETA1  = THETA1_OLD                                        !
        PHI0    = PHI0_OLD                                          !
        PHI1    = PHI1_OLD                                          !
        E0      = E0_OLD                                            !
        E1      = E1_OLD                                            !
        LI      = L_OLD                                             !
!
      END IF                                                        !
!
      GO TO 888                                                     !
!
!  Stops:
!
 200  WRITE(IUO1,100) JAT                                           !
      STOP                                                          !
 201  WRITE(IUO1,101)                                               !
      STOP                                                          !
!
!  Formats
!
  10  FORMAT(5X,I1,2X,I2,2X,I4,2X,I2,2X,I3,2X,I3,2X,I3,2X,F8.2,2X,F8.2)
  20  FORMAT(1X,I3,1X,I4,1X,I2,1X,I3,1X,F6.3,1X,F6.3,1X,F6.2,       &
             1X,F6.2,1X,F8.2)
!
 100  FORMAT(15X,'<<<<<  THE VALUE OF L EST IS TOO LARGE FOR ATOM', &
                 ' : ',I2,'  >>>>>')
 101  FORMAT(15X,'<<<<<  WRONG VALUE OF THETA0 : THE DENOMINATOR ', &
                 'IS ZERO  >>>>>')
!
 300  FORMAT(///,17X,'---> CALCULATION OF THE SCATTERING FACTOR ',  &
                     'FOR BEAM ',A3,' DONE')
 301  FORMAT(17X,'---> STORED IN FILE ',A24)
!
 888  RETURN                                                        !
!
      END SUBROUTINE CALC_SF
!
!=======================================================================
!
      SUBROUTINE SW_SF(RJ,RJK,THRJ,PHIRJ,BETA,GAMMA,L,M,           &
                       FSPH,JAT,JE,*)
!
!  This routine computes a spherical wave scattering factor
!
!
!  Input variables:
!
!     RJ       :
!     RJK      :
!     THRJ     :
!     PHIRJ    :
!     BETA     :
!     GAMMA    :
!     L        :
!     M        :
!     JAT      :
!     JE       :
!
!  Output variable:
!
!     FSPH     : spherical-wave scattering amplitude
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 19 May 2021
!
!
      USE REAL_NUMBERS,              ONLY : ZERO,ONE,FOUR,HALF
      USE COMPLEX_NUMBERS
      USE PI_ETC,                    ONLY : PI
!
      USE CALC_TYPE
      USE CURRENT_BEAM
      USE CURRENT_F_TH
      USE CURRENT_T_MATRIX
      USE STORE_COEF,                ONLY : EXPF
!
      USE HANKEL_POLYNOMIALS
      USE LEGENDRE_FUNCTIONS,        ONLY : PLM
      USE WIGNER_ROTATIONS
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  L,M
      INTEGER, INTENT(IN)       ::  JAT,JE
!
      INTEGER                   ::  INTER,IEM
      INTEGER                   ::  NO1,NDUM
      INTEGER                   ::  MUMAX,NUMAX
      INTEGER                   ::  MU,NU,NU1,LP
      INTEGER                   ::  LPMIN
!
      INTEGER                   ::  MOD,INT
!
      REAL (WP), INTENT(IN)     ::  RJ,RJK,THRJ,PHIRJ
      REAL (WP), INTENT(IN)     ::  BETA,GAMMA
!
      REAL (WP)                 ::  SSMALL,CSTH
      REAL (WP)                 ::  PLMM(0:150,0:150)
      REAL (WP)                 ::  D(1-NL_M:NL_M-1,1-NL_M:NL_M-1,0:NL_M-1)
      REAL (WP)                 ::  A,B,C,A1
!
      REAL (WP)                 ::  COS,SIN,FLOAT,SQRT
!
      COMPLEX (WP), INTENT(OUT) ::  FSPH
!
      COMPLEX (WP)              ::  VKE,RHOJ,RHOJK
      COMPLEX (WP)              ::  HLM1,HLM2,HLM3,HLM4
      COMPLEX (WP)              ::  ALMU,BLMU
      COMPLEX (WP)              ::  HLM(0:NO_ST_M,0:NL_M-1)
      COMPLEX (WP)              ::  HLN(0:NO_ST_M,0:NL_M-1)
      COMPLEX (WP)              ::  SLP,SNU,SMU
!
      COMPLEX (WP)              ::  EXP,CMPLX
!
      DATA  SSMALL  / 0.0001E0_WP /
!
      IF(ITL == 1) VKE = VK(JE)                                     !
!
      A     = ONE                                                   !
      INTER = 0                                                     !
      RHOJ  = VKE * RJ                                              !
      RHOJK = VKE * RJK                                             !
      HLM1  = ONEC                                                  !
      HLM2  = ONEC                                                  !
      HLM3  = ONEC                                                  !
      HLM4  = ONEC                                                  !
      IEM   = 1                                                     !
      CSTH  = COS(BETA)                                             !
!
      IF((IFTHET == 0) .OR. (THRJ < SSMALL)) THEN                   !
        INTER = 1                                                   !
        BLMU  = SQRT( FOUR * PI / FLOAT(L + L + 1) ) *            & !
                EXP( - IC * M * (PHIRJ - PI) )                      !
      END IF                                                        !
!
      CALL PLM(LMAX(JAT,JE),CSTH,PLMM)                              !
!
      IF(I_BASIS == 0) NO1 = 0                                      !
      IF(I_BASIS == 1) THEN                                         !
        IF(NO == 8) THEN                                            !
          NO1 = LMAX(JAT,JE) + 1                                    !
        ELSE                                                        !
          NO1 = NO                                                  !
        END IF                                                      !
        CALL POLHAN(I_BASIS,NO1,LMAX(JAT,JE),RHOJ,HLM)              !
         IF(IEM == 0) THEN                                          !
          HLM4 = HLM(0,L)                                           !
        END IF                                                      !
        IF(RJK > SSMALL) THEN                                       !
          NDUM = 0                                                  !
          CALL POLHAN(I_BASIS,NDUM,LMAX(JAT,JE),RHOJK,HLN)          !
        END IF                                                      !
        CALL DJMN2(THRJ,D,L,0)                                      !
        A1 = ABS(D(0,M,L))                                          !
        IF( ((A1 < SSMALL) .AND. (IFTHET == 1))                   & !
                           .AND. (INTER  == 0) ) RETURN 1           !
      END IF                                                        !
!
      MUMAX = MIN(L,NO1)                                            !
      SMU   = ZEROC                                                 !
!
      DO MU = 0, MUMAX                                              !
        IF(MOD(MU,2) == 0) THEN                                     !
          B =   ONE                                                 !
        ELSE                                                        !
          B = - ONE                                                 !
          IF(SIN(BETA) < ZERO) THEN                                 !
            A = - ONE                                               !
          END IF                                                    !
        END IF                                                      !
        IF(I_BASIS <= 1) THEN                                       !
          ALMU = ONEC                                               !
          C    = ONE                                                !
        END IF                                                      !
        IF(I_BASIS == 0) GO TO 40                                   !
        IF(INTER   == 0) BLMU = CMPLX(D(M,0,L))                     !
        IF(MU > 0) THEN                                             !
          C    = B * FLOAT(L + L + 1) / EXPF(MU,L)                  !
          ALMU = (D(M,MU,L) * EXP( - IC * MU * GAMMA ) +          & !
                 B * EXP( IC * MU * GAMMA) * D(M,-MU,L) ) / BLMU    !
        ELSE                                                        !
          C    = ONE                                                !
          ALMU = CMPLX(D(M,0,L)) / BLMU                             !
        END IF                                                      !
!
  40    SNU   = ZEROC                                               !
        NU1   = INT(HALF * (NO1 - MU) + SSMALL)                     !
        NUMAX = MIN(NU1,L-MU)                                       !
        DO NU = 0, NUMAX                                            !
          SLP   = ZEROC                                             !
          LPMIN = MAX(MU,NU)                                        !
          DO LP = LPMIN, LMAX(JAT,JE)                               !
            IF(I_BASIS == 1) THEN                                   !
              HLM1 = HLM(NU,LP)                                     !
              IF(RJK > SSMALL) HLM3 = HLN(0,LP)                     !
            END IF                                                  !
            SLP = SLP + FLOAT(LP + LP + 1) * TL(LP,1,JAT,JE) *    & !
                        HLM1 * PLMM(LP,MU) * HLM3                   !
          END DO                                                    !
          IF(I_BASIS == 1) THEN                                     !
            HLM2 = HLM(MU+NU,L)                                     !
          END IF                                                    !
          SNU = SNU + SLP * HLM2                                    !
        END DO                                                      !
        SMU = SMU + SNU * C * ALMU * A * B                          !
      END DO                                                        !
      FSPH = SMU / (VKE * HLM4)                                     !
!
      END SUBROUTINE SW_SF
!
END MODULE SCATTERING_AMPLITUDE
