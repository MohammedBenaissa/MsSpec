!
!=======================================================================
!
MODULE TREATMENT_PED
!
!  This module provides subroutine to treat the PED data
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NDIM_M,NE_M,N_TH_M,N_PH_M
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE TREAT_PED(ISOM,N_INFILES,J_FILE,NP)
!
!  This routine sums up the calculations corresponding to different
!        absorbers or different planes when this has to be done
!             (parameter ISOM in the input data file).
!
!
!  Input variables :
!
!                       ISOM      :
!                       N_INFILES :  number of input files
!                       J_FILE    :  current input file index
!                       NP        :
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 28 May 2021
!
!
      USE REAL_NUMBERS,              ONLY : ZERO
!
      USE CURRENT_FINA_VAL
      USE CURRENT_INIT_VAL
      USE EXP_TYPE
      USE NON_BULK
      USE OUTUNITS
!
      IMPLICIT NONE
!
      INTEGER, PARAMETER         ::  N_HEAD  = 5000
      INTEGER, PARAMETER         ::  N_FILES = 1000
!
      CHARACTER (LEN = 72)       ::  HEAD(N_HEAD,N_FILES)
      CHARACTER (LEN = 13)       ::  OUTDATA
!
      INTEGER, INTENT(IN)        ::  ISOM,N_INFILES,J_FILE,NP
!
      INTEGER                    ::  JVOL,JTOT
      INTEGER                    ::  NHEAD
      INTEGER                    ::  JLINE
      INTEGER                    ::  NDP,NTT
      INTEGER                    ::  JPLAN,JEMET,JE,JEM,JPL
      INTEGER                    ::  J_FIXED,J_SCAN
      INTEGER                    ::  JPHI,JPHI2,JTHETA
      INTEGER                    ::  N_FIXED,N_SCAN,N_SCAN_R,NPHI_R
      INTEGER                    ::  JLIN,JLIN2,JL
      INTEGER                    ::  JF,NF
      INTEGER                    ::  I_EXT,I_SO,ICHKDIR,IDICHR
      INTEGER                    ::  ISPIN,IE,IPH_1,IPHI,ISFLIP,ITHETA
      INTEGER                    ::  NPLAN,NEMET,NTHETA,NPHI,NE
!
      INTEGER                    ::  INT
!
      REAL (WP)                  ::  TAB(NDIM_M,4)
      REAL (WP)                  ::  ECIN(NE_M),DTHETA(N_TH_M),DPHI(N_PH_M)
      REAL (WP)                  ::  SSMALL
      REAL (WP)                  ::  XINCRF,FIX_STEP
      REAL (WP)                  ::  THETA,RTHETA
      REAL (WP)                  ::  TOTDIF_1,TOTDIR_1
      REAL (WP)                  ::  TOTDIF_2,TOTDIR_2
      REAL (WP)                  ::  VOLDIF_1,VOLDIR_1
      REAL (WP)                  ::  VOLDIF_2,VOLDIR_2
      REAL (WP)                  ::  TOTDIF2_1,TOTDIR2_1
      REAL (WP)                  ::  TOTDIF2_2,TOTDIR2_2
      REAL (WP)                  ::  VOLDIF2_1,VOLDIR2_1
      REAL (WP)                  ::  VOLDIF2_2,VOLDIR2_2
      REAL (WP)                  ::  SF_1,SR_1,SF_2,SR_2
      REAL (WP)                  ::  SF2_1,SR2_1,SF2_2,SR2_2
      REAL (WP)                  ::  FIX0,FIX1
!
      DATA JVOL,JTOT  / 0, -1       /
      DATA SSMALL     / 0.0001E0_WP /
!
      REWIND IUO2                                                   !
!
!  Reading and storing the headers:
!
      NHEAD = 0                                                     !
      DO JLINE = 1 ,N_HEAD                                          !
        READ(IUO2,888) HEAD(JLINE,J_FILE)                           !
        NHEAD = NHEAD + 1                                           !
        IF(HEAD(JLINE,J_FILE)(1:6) == '      ') GO TO 333           !
      END DO                                                        !
!
 333  CONTINUE                                                      !
!
      READ(IUO2,15) SPECTRO,OUTDATA                                 !
      READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,  & !
                   IE,IPH_1,I_EXT                                   !
!
      IF(I_EXT == 2) THEN                                           !
         IPH_1 = 0                                                  !
      END IF                                                        !
!
      IF(ISOM == 0) THEN                                            !
!
!........  ISOM = 0 : case of independent input files  .................
!
        READ(IUO2,1) NPLAN,NEMET,NTHETA,NPHI,NE                     !
!
        IF(IPH_1 == 1) THEN                                         !
          N_FIXED = NPHI                                            !
          FIX0    = PHI0                                            !
          FIX1    = PHI1                                            !
          N_SCAN  = NTHETA                                          !
        ELSE
          N_FIXED = NTHETA                                          !
          FIX0    = THETA0                                          !
          FIX1    = THETA1                                          !
          IF(STEREO == 'YES') THEN                                  !
            NPHI = INT( (PHI1 - PHI0) * FLOAT(NTHETA - 1) /       & !
                        (THETA1 - THETA0) + SSMALL ) + 1            !
            IF(NTHETA * NPHI > N_PH_M) GO TO 37                     !
          END IF                                                    !
          N_SCAN = NPHI                                             !
        END IF                                                      !
!
        IF(I_EXT == -1) THEN                                        !
          N_SCAN = 2 * N_SCAN                                       !
        END IF                                                      !
!
        IF((I_EXT == 0) .OR. (I_EXT == 1)) THEN                     !
          NDP     = NEMET * NTHETA * NPHI * NE                      !
        ELSE IF(I_EXT == -1) THEN                                   !
          NDP     = NEMET * NTHETA * NPHI * NE * 2                  !
        ELSE IF(I_EXT == 2) THEN                                    !
          NDP     = NEMET * NTHETA * NE                             !
          N_FIXED = NTHETA                                          !
          N_SCAN  = NPHI                                            !
          IF((N_FIXED > N_TH_M) .OR. (N_FIXED > N_PH_M)) GO TO 35    !
        END IF                                                      !
!
        NTT = NPLAN * NDP                                           !
        IF(NTT > NDIM_M) GO TO 5                                    !
!
        DO JPLAN = 1, NPLAN                                         !
         DO JEMET = 1, NEMET                                        !
          DO JE = 1, NE                                             !
!
           DO J_FIXED = 1, N_FIXED                                  !
            IF(N_FIXED > 1) THEN                                    !
              XINCRF = FLOAT(J_FIXED - 1) * (FIX1 - FIX0) /       & !
                       FLOAT(N_FIXED - 1)                           !
            ELSE IF(N_FIXED == 1) THEN                              !
              XINCRF = ZERO                                         !
            END IF                                                  !
            IF(IPH_1 == 1) THEN                                     !
              JPHI = J_FIXED                                        !
            ELSE                                                    !
              THETA  = THETA0 + XINCRF                              !
              JTHETA = J_FIXED                                      !
              IF((ABS(THETA) > 90.E0_WP) .AND. (I_EXT /= 2)) GO TO 11
            END IF                                                  !
            IF(STEREO == ' NO') THEN                                !
              N_SCAN_R = N_SCAN                                     !
            ELSE                                                    !
              RTHETA   = THETA * 0.017453E0_WP                      !
              FIX_STEP = (FIX1 - FIX0) / FLOAT(N_FIXED - 1)         !
              N_SCAN_R = INT( (PHI1 - PHI0) * SIN(RTHETA) /       & !
                               FIX_STEP + SSMALL ) + 1              !
            END IF                                                  !
!
            DO J_SCAN = 1, N_SCAN_R                                 !
             IF(IPH_1 == 1) THEN                                    !
              JTHETA = J_SCAN                                       !
             ELSE                                                   !
              JPHI   = J_SCAN                                       !
             END IF                                                 !
!
             JLIN = (JPLAN - 1) * NDP +                           & !
                    (JEMET - 1) * NE * N_FIXED * N_SCAN +         & !
                    (JE - 1) * N_FIXED * N_SCAN +                 & !
                    (JTHETA - 1) * NPHI + JPHI                      !
!
             IF(I_EXT <= 0) THEN                                    !
               IF(STEREO == ' NO') THEN                             !
                 JPHI2 = JPHI                                       !
               ELSE                                                 !
                 JPHI2 = (JTHETA - 1) * NPHI + JPHI                 !
               END IF                                               !
             ELSE                                                   !
               JPHI2 = JTHETA                                       !
             END IF                                                 !
!
             READ(IUO2,2) JPL                                       !
             IF(JPLAN == JPL) THEN                                  !
               BACKSPACE IUO2                                       !
               IF(IDICHR == 0) THEN                                 !
                 READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2), & !
                              ECIN(JE),TAB(JLIN,1),TAB(JLIN,2)      !
                 IF(I_EXT == -1) THEN                               !
                   JLIN2 = NTT + JLIN                               !
                   READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)          !
                 END IF                                             !
               ELSE                                                 !
                 READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),& !
                               ECIN(JE),TAB(JLIN,1),TAB(JLIN,2),  & !
                               TAB(JLIN,3),TAB(JLIN,4)              !
                 IF(I_EXT == -1) THEN                               !
                   JLIN2 = NTT + JLIN                               !
                   READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),&
                                 ECIN(JE),TAB(JLIN2,1),TAB(JLIN2,2),&
                                 TAB(JLIN2,3),TAB(JLIN2,4)          !
                 END IF                                             !
               END IF                                               !
             ELSE                                                   !
               BACKSPACE IUO2                                       !
               DO JL = JLIN, JPLAN * NDP                            !
                 TAB(JL,1) = ZERO                                   !
                 TAB(JL,2) = ZERO                                   !
                 TAB(JL,3) = ZERO                                   !
                 TAB(JL,4) = ZERO                                   !
               END DO                                               !
               GO TO 10                                             !
             END IF                                                 !
            END DO                                                  !
           END DO                                                   !
  11       CONTINUE                                                 !
          END DO                                                    !
         END DO                                                     !
  10     CONTINUE                                                   !
        END DO                                                      !
!
        REWIND IUO2                                                 !
!
!  Skipping the NHEAD lines of headers before rewriting:
!
        DO JLINE = 1, NHEAD                                         !
          READ(IUO2,888) HEAD(JLINE,J_FILE)                         !
        END DO                                                      !
!
        WRITE(IUO2,15) SPECTRO,OUTDATA                              !
        WRITE(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        WRITE(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM                     !
!
        DO JE = 1, NE                                               !
         DO JTHETA = 1, NTHETA                                      !
          IF(STEREO == ' NO') THEN                                  !
            NPHI_R = NPHI                                           !
          ELSE                                                      !
            RTHETA   = DTHETA(JTHETA) * 0.017453E0_WP               !
            FIX_STEP = (THETA1 - THETA0) / FLOAT(NTHETA - 1)        !
            NPHI_R   = INT( (PHI1 - PHI0) * SIN(RTHETA) /         & !
                             FIX_STEP + SSMALL ) + 1                !
            NPHI     = INT( (PHI1 - PHI0) / FIX_STEP + SSMALL ) + 1 !
          END IF                                                    !
          DO JPHI = 1, NPHI_R                                       !
           TOTDIF_1 = ZERO                                          !
           TOTDIR_1 = ZERO                                          !
           VOLDIF_1 = ZERO                                          !
           VOLDIR_1 = ZERO                                          !
           TOTDIF_2 = ZERO                                          !
           TOTDIR_2 = ZERO                                          !
           VOLDIF_2 = ZERO                                          !
           VOLDIR_2 = ZERO                                          !
           IF(I_EXT == -1) THEN                                     !
             TOTDIF2_1 = ZERO                                       !
             TOTDIR2_1 = ZERO                                       !
             VOLDIF2_1 = ZERO                                       !
             VOLDIR2_1 = ZERO                                       !
             TOTDIF2_2 = ZERO                                       !
             TOTDIR2_2 = ZERO                                       !
             VOLDIF2_2 = ZERO                                       !
             VOLDIR2_2 = ZERO                                       !
           END IF                                                   !
!
           DO JPLAN = 1, NPLAN                                      !
!
            SF_1 = ZERO                                             !
            SR_1 = ZERO                                             !
            SF_2 = ZERO                                             !
            SR_2 = ZERO                                             !
            IF(I_EXT == -1) THEN                                    !
              SF2_1 = ZERO                                          !
              SR2_1 = ZERO                                          !
              SF2_2 = ZERO                                          !
              SR2_2 = ZERO                                          !
            END IF                                                  !
!
            DO JEMET = 1, NEMET                                     !
             JLIN = (JPLAN - 1) * NDP +                           & !
                    (JEMET - 1) * NE * NTHETA * NPHI +            & !
                    (JE - 1) * NTHETA * NPHI +                    & !
                    (JTHETA - 1) * NPHI + JPHI                      !
             SF_1 = SF_1 + TAB(JLIN,2)                              !
             SR_1 = SR_1 + TAB(JLIN,1)                              !
             IF(I_EXT == -1) THEN                                   !
               JLIN2 = NTT + JLIN                                   !
               SF2_1 = SF2_1 + TAB(JLIN2,2)                         !
               SR2_1 = SR2_1 + TAB(JLIN2,1)                         !
             END IF                                                 !
             IF(IDICHR >= 1) THEN                                   !
               SF_2 = SF_2 + TAB(JLIN,4)                            !
               SR_2 = SR_2 + TAB(JLIN,3)                            !
               IF(I_EXT == -1) THEN                                 !
                 JLIN2 = NTT + JLIN                                 !
                 SF2_2 = SF2_2 + TAB(JLIN2,4)                       !
                 SR2_2 = SR2_2 + TAB(JLIN2,3)                       !
               END IF                                               !
             END IF                                                 !
            END DO                                                  !
            IF(I_EXT <= 0) THEN                                     !
               IF(STEREO == ' NO') THEN                             !
                 JPHI2 = JPHI                                       !
               ELSE                                                 !
                 JPHI2 = (JTHETA - 1) * NPHI + JPHI                 !
               END IF                                               !
             ELSE                                                   !
              JPHI2 = JTHETA                                        !
            END IF                                                  !
            IF(IDICHR == 0) THEN                                    !
              WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),     & !
                            ECIN(JE),SR_1,SF_1                      !
              IF(I_EXT == -1) THEN                                  !
                WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),   & !
                              ECIN(JE),SR2_1,SF2_1                  !
              END IF                                                !
            ELSE                                                    !
              WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),    & !
                             ECIN(JE),SR_1,SF_1,SR_2,SF_2           !
              IF(I_EXT == -1) THEN                                  !
                WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),  & !
                               ECIN(JE),SR2_1,SF2_1,SR2_2,SF2_2     !
              END IF                                                !
            END IF                                                  !
            IF(JPLAN > NONVOL(J_FILE)) THEN                         !
              VOLDIF_1 = VOLDIF_1 + SF_1                            !
              VOLDIR_1 = VOLDIR_1 + SR_1                            !
              IF(I_EXT == -1) THEN                                  !
                VOLDIF2_1 = VOLDIF2_1 + SF2_1                       !
                VOLDIR2_1 = VOLDIR2_1 + SR2_1                       !
              END IF                                                !
              IF(IDICHR >= 1) THEN                                  !
                VOLDIF_2 = VOLDIF_2 + SF_2                          !
                VOLDIR_2 = VOLDIR_1 + SR_2                          !
                IF(I_EXT == -1) THEN                                !
                  VOLDIF2_2 = VOLDIF2_2 + SF2_2                     !
                  VOLDIR2_2 = VOLDIR2_1 + SR2_2                     !
                END IF                                              !
              END IF                                                !
            END IF                                                  !
            TOTDIF_1 = TOTDIF_1 + SF_1                              !
            TOTDIR_1 = TOTDIR_1 + SR_1                              !
            IF(I_EXT == -1) THEN                                    !
              TOTDIF2_1 = TOTDIF2_1 + SF2_1                         !
              TOTDIR2_1 = TOTDIR2_1 + SR2_1                         !
            END IF                                                  !
            IF(IDICHR >= 1) THEN                                    !
              TOTDIF_2 = TOTDIF_2 + SF_2                            !
              TOTDIR_2 = TOTDIR_2 + SR_2                            !
              IF(I_EXT == -1) THEN                                  !
                TOTDIF2_2 = TOTDIF2_2 + SF2_2                       !
                TOTDIR2_2 = TOTDIR2_2 + SR2_2                       !
              END IF                                                !
            END IF                                                  !
           END DO                                                   !
           IF(IDICHR == 0) THEN                                     !
             WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),       & !
                           ECIN(JE),VOLDIR_1,VOLDIF_1               !
             IF(I_EXT == -1) THEN                                   !
               WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),     & !
                             ECIN(JE),VOLDIR2_1,VOLDIF2_1           !
             END IF                                                 !
             WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),       & !
                           ECIN(JE),TOTDIR_1,TOTDIF_1               !
             IF(I_EXT == -1) THEN                                   !
               WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),     & !
                             ECIN(JE),TOTDIR2_1,TOTDIF2_1           !
             END IF                                                 !
           ELSE                                                     !
             WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),      & !
                            ECIN(JE),VOLDIR_1,VOLDIF_1,           & !
                            VOLDIR_2,VOLDIF_2                       !
            IF(I_EXT == -1) THEN                                    !
               WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),    & !
                              ECIN(JE),VOLDIR2_1,VOLDIF2_1,       & !
                              VOLDIR2_2,VOLDIF2_2                   !
             END IF                                                 !
             WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),      & !
                            ECIN(JE),TOTDIR_1,TOTDIF_1,           & !
                            TOTDIR_2,TOTDIF_2                       !
             IF(I_EXT == -1) THEN                                   !
               WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),    & !
                              ECIN(JE),TOTDIR2_1,TOTDIF2_1,       & !
                              TOTDIR2_2,TOTDIF2_2                   !
             END IF                                                 !
           END IF                                                   !
          END DO                                                    !
         END DO                                                     !
        END DO                                                      !
!
      ELSE                                                          !
!
!........ ISOM not= 0 : multiple input files to be summed up  ..........
!
        READ(IUO2,7) NTHETA,NPHI,NE                                 !
!
        IF(IPH_1 == 1) THEN                                         !
          N_FIXED = NPHI                                            !
          FIX0    = PHI0                                            !
          FIX1    = PHI1                                            !
          N_SCAN  = NTHETA                                          !
        ELSE                                                        !
          N_FIXED = NTHETA                                          !
          FIX0    = THETA0                                          !
          FIX1    = THETA1                                          !
          IF(STEREO == 'YES') THEN                                  !
            NPHI = INT( (PHI1 - PHI0) * FLOAT(NTHETA - 1) /       & !
                        (THETA1 - THETA0) + SSMALL ) + 1            !
            IF(NTHETA*NPHI > N_PH_M) GO TO 37                       !
          END IF                                                    !
          N_SCAN = NPHI                                             !
        END IF                                                      !
!
        IF(I_EXT == -1) THEN                                        !
          N_SCAN = 2 * N_SCAN                                       !
        END IF                                                      !
!
        IF((I_EXT == 0) .OR. (I_EXT == 1)) THEN                     !
          NDP     = NTHETA * NPHI * NE                              !
        ELSE IF(I_EXT == -1) THEN                                   !
          NDP     = NTHETA * NPHI * NE * 2                          !
        ELSE IF(I_EXT == 2) THEN                                    !
          NDP     = NTHETA * NE                                     !
          N_FIXED = NTHETA                                          !
          N_SCAN  = NPHI                                            !
          IF((N_FIXED > N_TH_M) .OR. (N_FIXED > N_PH_M)) GO TO 35   !
        END IF                                                      !
!
        NTT = N_INFILES * NDP                                       !
        IF(NTT > NDIM_M) GO TO 5                                    !
!
        IF(ISOM == 1) THEN                                          !
          NPLAN = NP                                                !
          NF    = NP                                                !
        ELSE IF(ISOM == 2) THEN                                     !
          NEMET = N_INFILES                                         !
          NF    = N_INFILES                                         !
          NPLAN = 1                                                 !
        END IF                                                      !
!
        DO JF = 1, NF                                               !
!
!  Reading the headers for each file:
!
          IF(JF > 1) THEN                                           !
            DO JLINE = 1, NHEAD                                     !
              READ(IUO2,888) HEAD(JLINE,JF)                         !
            END DO                                                  !
          END IF                                                    !
!
          DO JE = 1,NE                                              !
!
           DO J_FIXED = 1, N_FIXED                                  !
            IF(N_FIXED > 1) THEN                                    !
              XINCRF = FLOAT(J_FIXED - 1) * (FIX1 - FIX0) /       & !
                       FLOAT(N_FIXED - 1)                           !
            ELSE IF(N_FIXED == 1) THEN                              !
              XINCRF = ZERO                                         !
            END IF                                                  !
            IF(IPH_1 == 1) THEN                                     !
              JPHI = J_FIXED                                        !
            ELSE                                                    !
              THETA  = THETA0 + XINCRF                              !
              JTHETA = J_FIXED                                      !
              IF((ABS(THETA) > 90.E0_WP) .AND. (I_EXT /= 2)) GO TO 12
            END IF                                                  !
            IF(STEREO == ' NO') THEN                                !
              N_SCAN_R = N_SCAN                                     !
            ELSE                                                    !
              RTHETA   = THETA * 0.017453E0_WP                      !
              FIX_STEP = (FIX1 - FIX0) / FLOAT(N_FIXED - 1)         !
              N_SCAN_R = INT( (PHI1 - PHI0) * SIN(RTHETA )/       & !
                               FIX_STEP + SSMALL ) + 1              !
            END IF                                                  !
!
            DO J_SCAN = 1, N_SCAN_R                                 !
             IF(IPH_1 == 1) THEN                                    !
              JTHETA = J_SCAN                                       !
             ELSE                                                   !
              JPHI   = J_SCAN                                       !
             END IF                                                 !
!
             JLIN = (JF - 1) * NDP + (JE - 1) * N_FIXED * N_SCAN +& !
                    (JTHETA - 1) * NPHI + JPHI                      !

             IF(I_EXT <= 0) THEN                                    !
               IF(STEREO == ' NO') THEN                             !
                 JPHI2 = JPHI                                       !
               ELSE
                 JPHI2 = (JTHETA - 1) * NPHI + JPHI                 !
               END IF                                               !
             ELSE                                                   !
               JPHI2 = JTHETA                                       !
             END IF                                                 !
!
             IF(ISOM == 1) THEN                                     !
               READ(IUO2,2) JPL                                     !
               IF(JF == JPL) THEN                                   !
                 BACKSPACE IUO2                                     !
                 IF(IDICHR == 0) THEN                               !
                   READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2), &
                                ECIN(JE),TAB(JLIN,1),TAB(JLIN,2)    !
                   IF(I_EXT == -1) THEN                             !
                     JLIN2 = NTT + JLIN                             !
                     READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)        !
                   END IF                                           !
                 ELSE                                               !
                   READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),&
                                 ECIN(JE),TAB(JLIN,1),TAB(JLIN,2),& !
                                 TAB(JLIN,3),TAB(JLIN,4)            !
                   IF(I_EXT == -1) THEN                             !
                     JLIN2 = NTT + JLIN                             !
                     READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),        & !
                                   DPHI(JPHI2),ECIN(JE),          & !
                                   TAB(JLIN2,1),TAB(JLIN2,2),     & !
                                   TAB(JLIN2,3),TAB(JLIN2,4)        !
                   END IF                                           !
                 END IF                                             !
               ELSE                                                 !
                 BACKSPACE IUO2                                     !
                 DO JLINE = 1, NHEAD                                !
                   BACKSPACE IUO2                                   !
                 END DO                                             !
                 DO JL = JLIN, JF * NDP                             !
                   TAB(JL,1) = ZERO                                 !
                   TAB(JL,2) = ZERO                                 !
                   TAB(JL,3) = ZERO                                 !
                   TAB(JL,4) = ZERO                                 !
                 END DO                                             !
                 GO TO 13                                           !
               END IF                                               !
             ELSE IF(ISOM == 2) THEN                                !
               IF(IDICHR == 0) THEN                                 !
                 READ(IUO2,2) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2), & !
                              ECIN(JE),TAB(JLIN,1),TAB(JLIN,2)      !
                 IF(I_EXT == -1) THEN                               !
                   JLIN2 = NTT + JLIN                               !
                   READ(IUO2,25) TAB(JLIN2,1),TAB(JLIN2,2)          !
                 END IF                                             !
               ELSE                                                 !
                 READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),DPHI(JPHI2),& !
                               ECIN(JE),TAB(JLIN,1),TAB(JLIN,2),  & !
                               TAB(JLIN,3),TAB(JLIN,4)              !
                 IF(I_EXT == -1) THEN                               !
                   JLIN2 = NTT + JLIN                               !
                   READ(IUO2,22) JPL,JEM,DTHETA(JTHETA),          & !
                                 DPHI(JPHI2),ECIN(JE),            & !
                                 TAB(JLIN2,1),TAB(JLIN2,2),       & !
                                 TAB(JLIN2,3),TAB(JLIN2,4)          !
                   END IF                                           !
               END IF                                               !
             END IF                                                 !
            END DO                                                  !
 12         CONTINUE                                                !
           END DO                                                   !
          END DO                                                    !
 13       CONTINUE                                                  !
        END DO                                                      !
!
        REWIND IUO2                                                 !
!
!  Writing the headers:
!
        DO JLINE = 1, 2                                             !
          WRITE(IUO2,888) HEAD(JLINE,1)                             !
        END DO                                                      !
        DO JF = 1, N_INFILES                                        !
          DO JLINE = 3, 6                                           !
            WRITE(IUO2,888) HEAD(JLINE,JF)                          !
          END DO                                                    !
          WRITE(IUO2,888) HEAD(2,JF)                                !
        END DO                                                      !
        DO JLINE = 7, NHEAD                                         !
          WRITE(IUO2,888) HEAD(JLINE,1)                             !
        END DO                                                      !
!
        WRITE(IUO2,15) SPECTRO,OUTDATA                              !
        WRITE(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE
        WRITE(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM                     !
!
        IF(ISOM == 1) THEN                                          !
!
          DO JE = 1, NE                                             !
!
           DO JTHETA = 1, NTHETA                                    !
            IF(STEREO == ' NO') THEN                                !
              NPHI_R   = NPHI                                       !
            ELSE                                                    !
              RTHETA   = DTHETA(JTHETA) * 0.017453E0_WP             !
              FIX_STEP = (THETA1 - THETA0) / FLOAT(NTHETA - 1)      !
              NPHI_R   = INT( (PHI1 - PHI0) * SIN(RTHETA) /       & !
                               FIX_STEP + SSMALL ) + 1              !
              NPHI     = INT( (PHI1 - PHI0) / FIX_STEP + SSMALL ) + 1
            END IF                                                  !
            DO JPHI = 1, NPHI_R                                     !
!
             TOTDIF_1 = ZERO                                        !
             TOTDIR_1 = ZERO                                        !
             VOLDIF_1 = ZERO                                        !
             VOLDIR_1 = ZERO                                        !
             TOTDIF_2 = ZERO                                        !
             TOTDIR_2 = ZERO                                        !
             VOLDIF_2 = ZERO                                        !
             VOLDIR_2 = ZERO                                        !
             IF(I_EXT == -1) THEN                                   !
               TOTDIF2_1 = ZERO                                     !
               TOTDIR2_1 = ZERO                                     !
               VOLDIF2_1 = ZERO                                     !
               VOLDIR2_1 = ZERO                                     !
               TOTDIF2_2 = ZERO                                     !
               TOTDIR2_2 = ZERO                                     !
               VOLDIF2_2 = ZERO                                     !
               VOLDIR2_2 = ZERO                                     !
             END IF                                                 !
!
             DO JPLAN = 1, NPLAN                                    !
              JF = JPLAN                                            !
!
              JLIN = (JF - 1) * NDP + (JE - 1) * NTHETA * NPHI +  & !
                     (JTHETA - 1) * NPHI + JPHI                     !
!
              SR_1 = TAB(JLIN,1)                                    !
              SF_1 = TAB(JLIN,2)                                    !
              IF(I_EXT == -1) THEN                                  !
                 JLIN2 = NTT + JLIN                                 !
                 SF2_1 = TAB(JLIN2,2)                               !
                 SR2_1 = TAB(JLIN2,1)                               !
              END IF                                                !
              IF(I_EXT <= 0) THEN                                   !
                IF(STEREO == ' NO') THEN                            !
                  JPHI2 = JPHI                                      !
                ELSE                                                !
                  JPHI2 = (JTHETA - 1) * NPHI + JPHI                !
                END IF                                              !
              ELSE                                                  !
                JPHI2 = JTHETA                                      !
              END IF                                                !
              IF(IDICHR == 0) THEN                                  !
                WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),   & !
                              ECIN(JE),SR_1,SF_1                    !
                IF(I_EXT == -1) THEN                                !
                  WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2), & !
                                ECIN(JE),SR2_1,SF2_1                !
                END IF                                              !
              ELSE                                                  !
                SR_2 = TAB(JLIN,3)                                  !
                SF_2 = TAB(JLIN,4)                                  !
                IF(I_EXT == -1) THEN                                !
                  JLIN2 = NTT + JLIN                                !
                  SF2_2 = TAB(JLIN2,4)                              !
                  SR2_2 = TAB(JLIN2,3)                              !
                END IF                                              !
                WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),  & !
                               ECIN(JE),SR_1,SF_1,SR_2,SF_2         !
                IF(I_EXT == -1) THEN                                !
                  WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),& !
                                 ECIN(JE),SR2_1,SF2_1,SR2_2,SF2_2   !
                END IF                                              !
              END IF                                                !
              IF(NONVOL(JPLAN) == 0) THEN                           !
                VOLDIF_1 = VOLDIF_1 + SF_1                          !
                VOLDIR_1 = VOLDIR_1 + SR_1                          !
                IF(I_EXT == -1) THEN                                !
                  VOLDIF2_1 = VOLDIF2_1 + SF2_1                     !
                  VOLDIR2_1 = VOLDIR2_1 + SR2_1                     !
                END IF                                              !
                IF(IDICHR >= 1) THEN                                !
                  VOLDIF_2 = VOLDIF_2 + SF_2                        !
                  VOLDIR_2 = VOLDIR_2 + SR_2                        !
                  IF(I_EXT == -1) THEN                              !
                    VOLDIF2_2 = VOLDIF2_2 + SF2_2                   !
                    VOLDIR2_2 = VOLDIR2_1 + SR2_2                   !
                  END IF                                            !
                END IF                                              !
              END IF                                                !
              TOTDIF_1 = TOTDIF_1 + SF_1                            !
              TOTDIR_1 = TOTDIR_1 + SR_1                            !
              IF(I_EXT == -1) THEN                                  !
                TOTDIF2_1 = TOTDIF2_1 + SF2_1                       !
                TOTDIR2_1 = TOTDIR2_1 + SR2_1                       !
              END IF                                                !
              IF(IDICHR >= 1) THEN                                  !
                TOTDIF_2 = TOTDIF_2 + SF_2                          !
                TOTDIR_2 = TOTDIR_2 + SR_2                          !
                IF(I_EXT == -1) THEN                                !
                  TOTDIF2_2 = TOTDIF2_2 + SF2_2                     !
                  TOTDIR2_2 = TOTDIR2_2 + SR2_2                     !
                END IF                                              !
              END IF                                                !
             END DO                                                 !
!
             IF(IDICHR == 0) THEN                                   !
               WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),     & !
                             ECIN(JE),VOLDIR_1,VOLDIF_1             !
               IF(I_EXT == -1) THEN                                 !
                 WRITE(IUO2,3) JVOL,DTHETA(JTHETA),DPHI(JPHI2),   & !
                               ECIN(JE),VOLDIR2_1,VOLDIF2_1         !
               END IF                                               !
               WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),     & !
                             ECIN(JE),TOTDIR_1,TOTDIF_1             !
               IF(I_EXT == -1) THEN                                 !
                 WRITE(IUO2,3) JTOT,DTHETA(JTHETA),DPHI(JPHI2),   & !
                               ECIN(JE),TOTDIR2_1,TOTDIF2_1         !
               END IF                                               !
             ELSE                                                   !
               WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),    & !
                              ECIN(JE),VOLDIR_1,VOLDIF_1,         & !
                              VOLDIR_2,VOLDIF_2                     !
               IF(I_EXT == -1) THEN                                 !
                 WRITE(IUO2,23) JVOL,DTHETA(JTHETA),DPHI(JPHI2),  & !
                                ECIN(JE),VOLDIR2_1,VOLDIF2_1,     & !
                                VOLDIR2_2,VOLDIF2_2                 !
               END IF                                               !
               WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),    & !
                              ECIN(JE),TOTDIR_1,TOTDIF_1,         & !
                              TOTDIR_2,TOTDIF_2                     !
               IF(I_EXT == -1) THEN                                 !
                 WRITE(IUO2,23) JTOT,DTHETA(JTHETA),DPHI(JPHI2),  & !
                                ECIN(JE),TOTDIR2_1,TOTDIF2_1,     & !
                                TOTDIR2_2,TOTDIF2_2                 !
               END IF                                               !
             END IF                                                 !
!
            END DO                                                  !
           END DO                                                   !
          END DO                                                    !
        ELSE IF(ISOM == 2) THEN                                     !
          DO JE = 1, NE                                             !
!
           DO JTHETA = 1, NTHETA                                    !
            IF(STEREO == ' NO') THEN                                !
              NPHI_R   = NPHI                                       !
            ELSE                                                    !
              RTHETA   = DTHETA(JTHETA) * 0.017453E0_WP             !
              FIX_STEP = (THETA1 - THETA0) / FLOAT(NTHETA - 1)      !
              NPHI_R   = INT( (PHI1 - PHI0) * SIN(RTHETA) /       & !
                               FIX_STEP + SSMALL ) + 1              !
              NPHI      = INT( (PHI1 - PHI0) / FIX_STEP + SSMALL ) + 1
            END IF                                                  !
            DO JPHI = 1, NPHI_R                                     !
!
              SF_1 = ZERO                                           !
              SR_1 = ZERO                                           !
              SF_2 = ZERO                                           !
              SR_2 = ZERO                                           !
              IF(I_EXT == -1) THEN                                  !
                SF2_1 = ZERO                                        !
                SR2_1 = ZERO                                        !
                SF2_2 = ZERO                                        !
                SR2_2 = ZERO                                        !
              END IF                                                !
!
              DO JEMET = 1, NEMET                                   !
               JF = JEMET                                           !
!
               JLIN = (JF - 1) * NDP + (JE - 1) * NTHETA * NPHI + & !
                      (JTHETA - 1) * NPHI + JPHI                    !
!
               SF_1 = SF_1 + TAB(JLIN,2)                            !
               SR_1 = SR_1 + TAB(JLIN,1)                            !
               IF(I_EXT == -1) THEN                                 !
                 JLIN2 = NTT + JLIN                                 !
                 SF2_1 = SF2_1 + TAB(JLIN2,2)                       !
                 SR2_1 = SR2_1 + TAB(JLIN2,1)                       !
               END IF                                               !
               IF(IDICHR >= 1) THEN                                 !
                 SF_2 = SF_2 + TAB(JLIN,4)                          !
                 SR_2 = SR_2 + TAB(JLIN,3)                          !
                 IF(I_EXT == -1) THEN                               !
                   JLIN2 = NTT + JLIN                               !
                   SF2_2 = SF2_2 + TAB(JLIN2,4)                     !
                   SR2_2 = SR2_2 + TAB(JLIN2,3)                     !
                 END IF                                             !
               END IF                                               !
              END DO                                                !
              IF(I_EXT <= 0) THEN                                   !
                IF(STEREO == ' NO') THEN                            !
                  JPHI2 = JPHI                                      !
                ELSE                                                !
                  JPHI2 = (JTHETA - 1) * NPHI + JPHI                !
                END IF                                              !
              ELSE                                                  !
                JPHI2 = JTHETA                                      !
              END IF                                                !
              IF(IDICHR == 0) THEN                                  !
                WRITE(IUO2,3) JPL,DTHETA(JTHETA),DPHI(JPHI2),     & !
                              ECIN(JE),SR_1,SF_1                    !
                IF(I_EXT == -1) THEN                                !
                  WRITE(IUO2,3) JPLAN,DTHETA(JTHETA),DPHI(JPHI2), & !
                                ECIN(JE),SR2_1,SF2_1                !
                END IF                                              !
              ELSE                                                  !
                WRITE(IUO2,23) JPL,DTHETA(JTHETA),DPHI(JPHI2),    & !
                              ECIN(JE),SR_1,SF_1,SR_2,SF_2          !
                IF(I_EXT == -1) THEN                                !
                  WRITE(IUO2,23) JPLAN,DTHETA(JTHETA),DPHI(JPHI2),& !
                                 ECIN(JE),SR2_1,SF2_1,SR2_2,SF2_2   !
                END IF                                              !
              END IF                                                !
            END DO                                                  !
           END DO                                                   !
          END DO                                                    !
        END IF                                                      !
      END IF                                                        !
!
      GO TO 6                                                       !
!
   5  WRITE(IUO1,4)                                                 !
      STOP                                                          !
  35  WRITE(IUO1,36) N_FIXED                                        !
      STOP                                                          !
  37  WRITE(IUO1,38) NTHETA * NPHI                                  !
      STOP                                                          !
!
   1  FORMAT(2X,I3,2X,I2,2X,I4,2X,I4,2X,I4)
   2  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,          &
             2X,E12.6)
   3  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
   4  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL ',&
             'IN THE TREAT_PED SUBROUTINE - INCREASE NDIM_M ',      &
             '>>>>>>>>>>')
   7  FORMAT(I4,2X,I4,2X,I4)
   8  FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
   9  FORMAT(9(2X,I1),2X,I2)
  15  FORMAT(2X,A3,11X,A13)
  22  FORMAT(2X,I3,2X,I2,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,          &
             2X,E12.6,2X,E12.6,2X,E12.6)
  23  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,       &
             2X,E12.6,2X,E12.6)
  25  FORMAT(37X,E12.6,2X,E12.6)
  36  FORMAT(//,4X,'<<<<<<<<<<  DIMENSION OF NTH_M OR NPH_M TOO SMALL ', &
             'IN THE INCLUDE FILE >>>>>>>>>>',/,4X,                      &
              '<<<<<<<<<<                 SHOULD BE AT LEAST ',I6,       &
              '                 >>>>>>>>>>')
  38  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF NPH_M TOO SMALL ',     &
             'IN THE INCLUDE FILE >>>>>>>>>>',/,8X,                 &
              '<<<<<<<<<<             SHOULD BE AT LEAST ',I6,      &
              '             >>>>>>>>>>')
 888  FORMAT(A72)
!
   6  RETURN                                                        !
!
      END SUBROUTINE TREAT_PED
!
!=======================================================================
!
      SUBROUTINE WEIGHT_SUM(I_EXT,I_EXT_A,JEL)
!
!   This subroutine performs a weighted sum of the results
!     corresponding to different directions of the detector.
!     The directions and weights are read from an external input file
!
!   JEL is the electron undetected (i.e. for which the outgoing
!     directions are integrated over the unit sphere). It is always
!     1 for one electron spectroscopies (PED). For APECS, It can be
!     1 (photoelectron) or 2 (Auger electron) or even 0 (no electron
!     detected)
!
!  Input variables :
!
!                       ISOM      :
!                       I_EXT     :
!                       I_EXT_A   :
!                       JEL       :
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 19 May 2021
!
!
      USE REAL_NUMBERS,              ONLY : ZERO
!
      USE INFILES
      USE INUNITS
      USE OUTUNITS
!
      IMPLICIT NONE
!
      CHARACTER (LEN = 13)       ::  OUTDATA
      CHARACTER (LEN =  5)       ::  LIKE
      CHARACTER (LEN =  3)       ::  SPECTRO,SPECTRO2
!
      INTEGER, INTENT(IN)        ::  I_EXT,I_EXT_A,JEL
!
      INTEGER, PARAMETER         ::  N_MAX = 5810
      INTEGER, PARAMETER         ::  NPM   =   20
!
      INTEGER                    ::  JVOL,JTOT
      INTEGER                    ::  I_DIM
      INTEGER                    ::  I_SO,ICHKDIR,IDICHR
      INTEGER                    ::  IE,IPHI,ITHETA,ISOM
      INTEGER                    ::  I_SO_A,ICHKDIR_A,IDICHR_A
      INTEGER                    ::  IE_A,IPHI_A,ITHETA_A
      INTEGER                    ::  ISFLIP,ISPIN
      INTEGER                    ::  ISFLIP_A,ISPIN_A
      INTEGER                    ::  JDUM,NDUM,N_DUM1,N_DUM2
      INTEGER                    ::  NTHETA0,NPHI0
      INTEGER                    ::  NANGLE
      INTEGER                    ::  JE,JANGLE,JANGLE_A,JPLAN
      INTEGER                    ::  N_POINTS,N_POINTS_A
      INTEGER                    ::  NPHI,NTHETA,NE,NPLAN
      INTEGER                    ::  NPHI_A,NTHETA_A,NE_A
!
      REAL (WP)                  ::  SUMR_1(NPM,NE_M,N_MAX)
      REAL (WP)                  ::  SUMR_2(NPM,NE_M,N_MAX)
      REAL (WP)                  ::  SUMF_1(NPM,NE_M,N_MAX)
      REAL (WP)                  ::  SUMF_2(NPM,NE_M,N_MAX)
      REAL (WP)                  ::  SR_1,SF_1,SR_2,SF_2
      REAL (WP)                  ::  DTHETA(N_MAX),DPHI(N_MAX)
      REAL (WP)                  ::  DTHETAA(N_MAX),DPHIA(N_MAX)
      REAL (WP)                  ::  W(N_MAX),W_A(N_MAX),ECIN(NE_M)
      REAL (WP)                  ::  THETA,PHI,TH,PH,THA,PHA
!
      DATA JVOL,JTOT    / 0, -1 /
      DATA LIKE         / '-like' /
!
      REWIND IUO2                                                   !
!
      READ(IUO2,15) SPECTRO,OUTDATA                                 !
      IF(SPECTRO /= 'APC') THEN                                     !
        READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE!
        READ(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM                      !
        SPECTRO2 = 'XAS'                                            !
      ELSE                                                          !
        READ(IUO2,9) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,ITHETA,IE!
        READ(IUO2,9) ISPIN_A,IDICHR_A,I_SO_A,ISFLIP_A,ICHKDIR_A,  & !
                     IPHI_A,ITHETA_A,IE_A                           !
        READ(IUO2,8) NPHI,NTHETA,NE,NPLAN,ISOM                      !
        READ(IUO2,8) NPHI_A,NTHETA_A                                !
        IF(JEL == 1) THEN                                           !
          SPECTRO2 = 'AED'                                          !
        ELSE IF(JEL == 2) THEN                                      !
          SPECTRO2 = 'PED'                                          !
        ELSE IF(JEL == 0) THEN                                      !
          SPECTRO2 = 'XAS'                                          !
        END IF                                                      !
      END IF                                                        !
!
      IF(NPLAN > NPM) THEN                                          !
        WRITE(IUO1,4) NPLAN + 2                                     !
        STOP                                                        !
      END IF                                                        !
!
!  Reading the number of angular points
!
      IF(SPECTRO /= 'APC') THEN                                     !
        OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')               !
        READ(IUI11,1) N_POINTS                                      !
        READ(IUI11,5) I_DIM,N_DUM1,N_DUM2                           !
        N_POINTS_A = 1                                              !
      ELSE                                                          !
        IF(JEL == 1) THEN                                           !
          OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')             !
          READ(IUI11,1) N_POINTS                                    !
          READ(IUI11,5) I_DIM,N_DUM1,N_DUM2                         !
          IF(I_EXT_A == 0) THEN                                     !
            N_POINTS_A = NTHETA_A * NPHI_A                          !
          ELSE                                                      !
            OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')             !
            READ(IUI9,1) N_POINTS_A                                 !
            READ(IUI9,5) I_DIM,N_DUM1,N_DUM2                        !
          END IF                                                    !
          NTHETA0 = NTHETA_A                                        !
          NPHI0   = NPHI_A                                          !
        ELSE IF(JEL == 2) THEN                                      !
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')               !
          READ(IUI9,1) N_POINTS_A                                   !
          READ(IUI9,5) I_DIM,N_DUM1,N_DUM2                          !
          IF(I_EXT == 0) THEN                                       !
            N_POINTS = NTHETA * NPHI                                !
          ELSE                                                      !
            OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')           !
            READ(IUI11,1) N_POINTS                                  !
            READ(IUI11,5) I_DIM,N_DUM1,N_DUM2                       !
          END IF                                                    !
          NTHETA0 = NTHETA                                          !
          NPHI0   = NPHI                                            !
        ELSE IF(JEL == 0) THEN                                      !
          OPEN(UNIT=IUI11, FILE=INFILE11, STATUS='OLD')             !
          OPEN(UNIT=IUI9, FILE=INFILE9, STATUS='OLD')               !
          READ(IUI11,1) N_POINTS                                    !
          READ(IUI9,1) N_POINTS_A                                   !
          READ(IUI11,5) I_DIM,N_DUM1,N_DUM2                         !
          READ(IUI9,5) I_DIM,N_DUM1,N_DUM2                          !
        END IF                                                      !
      END IF                                                        !
!
      IF(SPECTRO /= 'APC') THEN                                     !
        NANGLE = 1                                                  !
      ELSE                                                          !
        IF(JEL == 1) THEN                                           !
          NANGLE = N_POINTS_A                                       !
        ELSE IF(JEL == 2) THEN                                      !
          NANGLE = N_POINTS                                         !
        ELSE IF(JEL == 0) THEN                                      !
          NANGLE = 1                                                !
        END IF                                                      !
      END IF                                                        !
!
!  Initialization of the arrays
!
      DO JE = 1, NE                                                 !
        DO JANGLE = 1, NANGLE                                       !
          DO JPLAN = 1, NPLAN + 2                                   !
            SUMR_1(JPLAN,JE,JANGLE) = ZERO                          !
            SUMF_1(JPLAN,JE,JANGLE) = ZERO                          !
            IF(IDICHR > 0) THEN                                     !
              SUMR_2(JPLAN,JE,JANGLE) = ZERO                        !
              SUMF_2(JPLAN,JE,JANGLE) = ZERO                        !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
!
!  Reading of the data to be angle integrated
!
      DO JE = 1, NE                                                 !
!
        DO JANGLE = 1, N_POINTS                                     !
          IF(I_EXT /= 0) READ(IUI6,2) TH,PH,W(JANGLE)               !
          DO JANGLE_A = 1, N_POINTS_A                               !
            IF((I_EXT_A /= 0) .AND. (JANGLE == 1)) THEN             !
               READ(IUI9,2) THA,PHA,W_A(JANGLE_A)                   !
            END IF                                                  !
!
            DO JPLAN = 1, NPLAN + 2                                 !
!
              IF(IDICHR == 0) THEN                                  !
                IF(SPECTRO /= 'APC') THEN                           !
                  READ(IUO2,3) JDUM,DTHETA(JANGLE),DPHI(JANGLE),  & !
                               ECIN(JE),SR_1,SF_1                   !
                ELSE                                                !
                  READ(IUO2,13) JDUM,DTHETA(JANGLE),DPHI(JANGLE), & !
                                ECIN(JE),DTHETAA(JANGLE_A),       & !
                                DPHIA(JANGLE_A),SR_1,SF_1           !
                END IF                                              !
              ELSE                                                  !
                IF(SPECTRO /= 'APC') THEN                           !
                  READ(IUO2,23) JDUM,DTHETA(JANGLE),DPHI(JANGLE), & !
                                ECIN(JE),SR_1,SF_1,SR_2,SF_2        !
                ELSE                                                !
                  READ(IUO2,24) JDUM,DTHETA(JANGLE),DPHI(JANGLE), & !
                                ECIN(JE),DTHETAA(JANGLE_A),       & !
                                DPHIA(JANGLE_A),SR_1,SF_1,SR_2,SF_2 !
                END IF                                              !
              END IF                                                !
!
              IF(JEL == 1) THEN                                     !
                SUMR_1(JPLAN,JE,JANGLE_A) = SUMR_1(JPLAN,JE,JANGLE_A) + &
                                            SR_1 * W(JANGLE)        !
                SUMF_1(JPLAN,JE,JANGLE_A) = SUMF_1(JPLAN,JE,JANGLE_A) + &
                                            SF_1 * W(JANGLE)        !
              ELSE IF(JEL == 2) THEN                                !
                SUMR_1(JPLAN,JE,JANGLE) = SUMR_1(JPLAN,JE,JANGLE) +&!
                                          SR_1 * W_A(JANGLE_A)      !
                SUMF_1(JPLAN,JE,JANGLE) = SUMF_1(JPLAN,JE,JANGLE) +&!
                                          SF_1 * W_A(JANGLE_A)      !
              ELSE IF(JEL == 0) THEN                                !
                SUMR_1(JPLAN,JE,1) = SUMR_1(JPLAN,JE,1) +         & !
                                     SR_1 * W(JANGLE) * W_A(JANGLE_A)
                SUMF_1(JPLAN,JE,1) = SUMF_1(JPLAN,JE,1) +         & !
                                     SF_1 * W(JANGLE) * W_A(JANGLE_A)
              END IF                                                !
              IF(IDICHR > 0) THEN                                   !
                IF(JEL == 1) THEN                                   !
                  SUMR_2(JPLAN,JE,JANGLE_A) = SUMR_2(JPLAN,JE,JANGLE_A) + &
                                              SR_2 * W(JANGLE)      !
                  SUMF_2(JPLAN,JE,JANGLE_A) = SUMF_2(JPLAN,JE,JANGLE_A) + &
                                              SF_2 * W(JANGLE)      !
                ELSE IF(JEL == 2) THEN                              !
                  SUMR_2(JPLAN,JE,JANGLE) = SUMR_2(JPLAN,JE,JANGLE) + &
                                            SR_2 * W_A(JANGLE_A)    !
                  SUMF_2(JPLAN,JE,JANGLE) = SUMF_2(JPLAN,JE,JANGLE) + &
                                            SF_2 * W_A(JANGLE_A)    !
                ELSE IF(JEL == 0) THEN                              !
                  SUMR_2(JPLAN,JE,1) = SUMR_2(JPLAN,JE,1) +       & !
                                       SR_2 * W(JANGLE) * W_A(JANGLE_A)
                  SUMF_2(JPLAN,JE,1) = SUMF_2(JPLAN,JE,1) +       & !
                                       SF_2 * W(JANGLE) * W_A(JANGLE_A)
                END IF                                              !
              END IF                                                !
!
            END DO                                                  !
!
          END DO                                                    !
          IF(I_EXT_A /= 0) THEN                                     !
            REWIND IUI9                                             !
            READ(IUI9,1) NDUM                                       !
            READ(IUI9,1) NDUM                                       !
          END IF                                                    !
        END DO                                                      !
!
        IF(I_EXT /= 0) THEN                                         !
          REWIND IUI11                                              !
          READ(IUI11,1) NDUM                                        !
          READ(IUI11,1) NDUM                                        !
        END IF                                                      !
      END DO                                                        !
!
      CLOSE(IUI11)                                                  !
      CLOSE(IUI9)                                                   !
      REWIND IUO2                                                   !
!
      WRITE(IUO2,16) SPECTRO2,LIKE,SPECTRO,OUTDATA                  !
      IF((SPECTRO /= 'APC') .OR. (JEL == 0)) THEN                   !
        WRITE(IUO2,19) ISPIN,IDICHR,I_SO,ISFLIP                     !
        WRITE(IUO2,18) NE,NPLAN,ISOM                                !
      ELSE IF(JEL == 1) THEN                                        !
        WRITE(IUO2,20) ISPIN_A,IDICHR_A,I_SO_A,ISFLIP_A,ICHKDIR_A,& !
                       IPHI_A,ITHETA_A,IE_A                         !
        WRITE(IUO2,21) NPHI0,NTHETA0,NE,NPLAN,ISOM                  !
      ELSE IF(JEL == 2) THEN                                        !
        WRITE(IUO2,20) ISPIN,IDICHR,I_SO,ISFLIP,ICHKDIR,IPHI,     & !
                       ITHETA,IE                                    !
        WRITE(IUO2,21) NPHI0,NTHETA0,NE,NPLAN,ISOM                  !
      END IF                                                        !
!
      DO JE = 1, NE                                                 !
        DO JANGLE = 1, NANGLE                                       !
          IF(SPECTRO == 'APC') THEN                                 !
            IF(JEL == 1) THEN                                       !
              THETA = DTHETAA(JANGLE)                               !
              PHI   = DPHIA(JANGLE)                                 !
            ELSE IF(JEL == 2) THEN                                  !
              THETA = DTHETA(JANGLE)                                !
              PHI   = DPHI(JANGLE)                                  !
            END IF                                                  !
          END IF                                                    !
!
          DO JPLAN = 1, NPLAN                                       !
            IF(IDICHR == 0) THEN                                    !
              IF((SPECTRO /= 'APC') .OR. (JEL == 0)) THEN           !
                WRITE(IUO2,33) JPLAN,ECIN(JE),                    & !
                               SUMR_1(JPLAN,JE,JANGLE),           & !
                               SUMF_1(JPLAN,JE,JANGLE)              !
              ELSE                                                  !
                WRITE(IUO2,34) JPLAN,THETA,PHI,ECIN(JE),          & !
                               SUMR_1(JPLAN,JE,JANGLE),           & !
                               SUMF_1(JPLAN,JE,JANGLE)              !
              END IF                                                !
            ELSE                                                    !
              IF((SPECTRO /= 'APC') .OR. (JEL == 0)) THEN           !
                WRITE(IUO2,43) JPLAN,ECIN(JE),                    & !
                               SUMR_1(JPLAN,JE,JANGLE),           & !
                               SUMF_1(JPLAN,JE,JANGLE),           & !
                               SUMR_2(JPLAN,JE,JANGLE),           & !
                               SUMF_2(JPLAN,JE,JANGLE)              !
              ELSE                                                  !
                WRITE(IUO2,44) JPLAN,THETA,PHI,ECIN(JE),          & !
                               SUMR_1(JPLAN,JE,JANGLE),           & !
                               SUMF_1(JPLAN,JE,JANGLE),           & !
                               SUMR_2(JPLAN,JE,JANGLE),           & !
                               SUMF_2(JPLAN,JE,JANGLE)              !
              END IF                                                !
            END IF                                                  !
          END DO                                                    !
!
          IF(IDICHR == 0) THEN                                      !
            IF((SPECTRO /= 'APC').OR.(JEL == 0)) THEN               !
              WRITE(IUO2,33) JVOL,ECIN(JE),                       & !
                             SUMR_1(NPLAN+1,JE,JANGLE),           & !
                             SUMF_1(NPLAN+1,JE,JANGLE)              !
              WRITE(IUO2,33) JTOT,ECIN(JE),                       & !
                             SUMR_1(NPLAN+2,JE,JANGLE),           & !
                             SUMF_1(NPLAN+2,JE,JANGLE)              !
            ELSE                                                    !
              WRITE(IUO2,34) JVOL,THETA,PHI,ECIN(JE),             & !
                             SUMR_1(NPLAN+1,JE,JANGLE),           & !
                             SUMF_1(NPLAN+1,JE,JANGLE)              !
              WRITE(IUO2,34) JTOT,THETA,PHI,ECIN(JE),             & !
                             SUMR_1(NPLAN+2,JE,JANGLE),           & !
                             SUMF_1(NPLAN+2,JE,JANGLE)              !
            END IF                                                  !
          ELSE                                                      !
            IF((SPECTRO /= 'APC') .OR. (JEL == 0)) THEN             !
              WRITE(IUO2,43) JVOL,ECIN(JE),                       & !
                             SUMR_1(NPLAN+1,JE,JANGLE),           & !
                             SUMF_1(NPLAN+1,JE,JANGLE),           & !
                             SUMR_2(NPLAN+1,JE,JANGLE),           & !
                             SUMF_2(NPLAN+1,JE,JANGLE)              !
              WRITE(IUO2,43) JTOT,ECIN(JE),                       & !
                             SUMR_1(NPLAN+2,JE,JANGLE),           & !
                             SUMF_1(NPLAN+2,JE,JANGLE),           & !
                             SUMR_2(NPLAN+2,JE,JANGLE),           & !
                             SUMF_2(NPLAN+2,JE,JANGLE)              !
            ELSE                                                    !
              WRITE(IUO2,44) JVOL,THETA,PHI,ECIN(JE),             & !
                             SUMR_1(NPLAN+1,JE,JANGLE),           & !
                             SUMF_1(NPLAN+1,JE,JANGLE),           & !
                             SUMR_2(NPLAN+1,JE,JANGLE),           & !
                             SUMF_2(NPLAN+1,JE,JANGLE)              !
              WRITE(IUO2,44) JTOT,THETA,PHI,ECIN(JE),             & !
                             SUMR_1(NPLAN+2,JE,JANGLE),           & !
                             SUMF_1(NPLAN+2,JE,JANGLE),           & !
                             SUMR_2(NPLAN+2,JE,JANGLE),           & !
                             SUMF_2(NPLAN+2,JE,JANGLE)              !
            END IF                                                  !
          END IF                                                    !
!
        END DO                                                      !
      END DO                                                        !
!
!  Formats:
!
   1  FORMAT(13X,I4)
   2  FORMAT(15X,F8.3,3X,F8.3,3X,E12.6)
   3  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
   4  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL ',  &
             'IN THE WEIGHT_SUM SUBROUTINE - INCREASE NPM TO ',I3,    &
             '>>>>>>>>>>')
   5  FORMAT(6X,I1,1X,I3,3X,I3)
   8  FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
   9  FORMAT(9(2X,I1),2X,I2)
  13  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,                &
             F6.2,2X,E12.6,2X,E12.6)
  15  FORMAT(2X,A3,11X,A13)
  16  FORMAT(2X,A3,A5,1X,A3,2X,A13)
  18  FORMAT(I4,2X,I3,2X,I1)
  19  FORMAT(4(2X,I1))
  20   FORMAT(8(2X,I1))
  21   FORMAT(I4,2X,I4,2X,I4,2X,I3,2X,I1)
  23  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,         &
             2X,E12.6,2X,E12.6)
  24  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,F6.2,2X,                &
             F6.2,2X,E12.6,2X,E12.6,2X,E12.6,2X,E12.6)
  33  FORMAT(2X,I3,2X,F8.2,2X,E12.6,2X,E12.6)
  34  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6)
  43  FORMAT(2X,I3,2X,F8.2,2X,E12.6,2X,E12.6,                         &
             2X,E12.6,2X,E12.6)
  44  FORMAT(2X,I3,2X,F6.2,2X,F6.2,2X,F8.2,2X,E12.6,2X,E12.6,         &
             2X,E12.6,2X,E12.6)
!
      END SUBROUTINE WEIGHT_SUM
!
!=======================================================================
!
      SUBROUTINE STOP_TREAT(N_INFILES,NPLAN,NEMET,NE,NTHETA,NTHETA_A, &
                            NPHI,NPHI_A,ISOM,I_EXT,I_EXT_A,SPECTRO)
!
!   This subroutine stops the code before the long MS calculations
!     when the dimensioning NDIM_M of the treatment routines
!     (treat_aed,treat_apc,treat_phd,treat_xas) is insufficient.
!
!  Input variables :
!
!                       N_INFILES :  number of input data files
!                       NPLAN     :  number of planes
!                       NEMET     :  number of absorbers
!                       NE        :  number of energy points
!                       NTHETA    :  number of theta points
!                       NTHETA_A  :  number of theta points for second electron
!                       NPHI      :  number of phi points
!                       NPHI_A    :  number of phi points for second electron
!                       ISOM      :  summation index
!                       I_EXT     :  external directions index
!                       I_EXT_1   :  external directions index (2nd electron)
!                       SPECTRO    type of spectroscopy
!
!
!
!
!   Author :  D. Sébilleau
!
!
!                                           Last modified : 18 May 2016
!
!
      USE OUTUNITS
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)       ::  SPECTRO
!
      INTEGER, INTENT(IN)        ::  N_INFILES,NPLAN,NEMET,NE
      INTEGER, INTENT(IN)        ::  NTHETA,NTHETA_A,NPHI,NPHI_A
      INTEGER, INTENT(IN)        ::  ISOM,I_EXT,I_EXT_A
!
      INTEGER                    ::  NDP,NTT
!
      IF(ISOM == 0) THEN                                            !
!1
!   Photoelectron diffraction case
!
        IF(SPECTRO == 'PED') THEN                                   !
          IF((I_EXT == 0) .OR. (I_EXT == 1)) THEN                   !
            NDP = NEMET * NTHETA * NPHI * NE                        !
          ELSE IF(I_EXT == -1) THEN                                 !
            NDP = NEMET * NTHETA * NPHI * NE * 2                    !
          ELSE IF(I_EXT == 2) THEN                                  !
            NDP = NEMET * NTHETA * NE                               !
          END IF                                                    !
          NTT = NPLAN * NDP                                         !
          IF(NTT > NDIM_M) GO TO 10                                 !
          IF((NTHETA > N_TH_M) .OR. (NPHI > N_PH_M)) GO TO 50       !
!
!   Auger electron diffraction case
!
        ELSE IF(SPECTRO == 'AED') THEN                              !
          IF((I_EXT_A == 0) .OR. (I_EXT_A == 1)) THEN               !
            NDP = NEMET * NTHETA_A * NPHI_A * NE                    !
          ELSE IF(I_EXT_A == -1) THEN                               !
            NDP = NEMET * NTHETA_A * NPHI_A * NE * 2                !
          ELSE IF(I_EXT_A == 2) THEN                                !
            NDP = NEMET * NTHETA_A * NE                             !
          END IF                                                    !
          NTT = NPLAN * NDP                                         !
          IF(NTT > NDIM_M) GO TO 20                                 !
          IF((NTHETA_A > N_TH_M) .OR. (NPHI_A > N_PH_M)) GO TO 50   !
!
!   X-ray absorption case
!
        ELSE IF(SPECTRO == 'XAS') THEN                              !
          NDP = NEMET * NE                                          !
          NTT = NPLAN * NDP                                         !
          IF(NTT > NDIM_M) GO TO 30                                 !
!
!   Auger Photoelectron coincidence spectroscopy case
!
        ELSE IF(SPECTRO == 'APC') THEN                              !
          IF((I_EXT == 0) .OR. (I_EXT == 1)) THEN                   !
            IF((I_EXT_A == 0) .OR. (I_EXT_A == 1)) THEN             !
              NDP = NEMET * NTHETA * NPHI * NE * NTHETA_A * NPHI_A  !
            ELSE IF(I_EXT_A == -1) THEN                             !
              NDP = NEMET * NTHETA * NPHI * NE * NTHETA_A * NPHI_A * 2
            ELSE IF(I_EXT_A == 2) THEN                              !
              NDP = NEMET * NTHETA * NPHI * NE * NTHETA_A           !
            END IF                                                  !
          ELSE IF(I_EXT == -1) THEN
            IF((I_EXT_A == 0) .OR. (I_EXT_A == 1)) THEN             !
              NDP = NEMET * NTHETA * NPHI * NE * 2 * NTHETA_A * NPHI_A
            ELSE IF(I_EXT_A == -1) THEN                             !
              NDP = NEMET * NTHETA * NPHI * NE * 2 * NTHETA_A * NPHI_A * 2
            ELSE IF(I_EXT_A == 2) THEN                              !
              NDP = NEMET * NTHETA * NPHI * NE * 2 * NTHETA_A       !
            END IF                                                  !
          ELSE IF(I_EXT == 2) THEN                                  !
            IF((I_EXT_A == 0) .OR. (I_EXT_A == 1)) THEN             !
              NDP = NEMET * NTHETA * NE * NTHETA_A * NPHI_A         !
            ELSE IF(I_EXT_A == -1) THEN                             !
              NDP = NEMET * NTHETA * NE * NTHETA_A * NPHI_A * 2     !
            ELSE IF(I_EXT_A == 2) THEN                              !
              NDP = NEMET * NTHETA * NE * NTHETA_A                  !
            END IF                                                  !
          END IF                                                    !
          NTT = NPLAN * NDP                                         !
          IF(NTT > NDIM_M) GO TO 40                                 !
          IF((NTHETA > N_TH_M) .OR. (NPHI > N_PH_M))     GO TO 50   !
          IF((NTHETA_A > N_TH_M) .OR. (NPHI_A > N_PH_M)) GO TO 50   !
          END IF                                                    !
!
      ELSE                                                          !
!
!   Photoelectron diffraction case
!
        IF(SPECTRO == 'PED') THEN                                   !
          IF((I_EXT == 0) .OR. (I_EXT == 1)) THEN                   !
            NDP = NTHETA * NPHI * NE                                !
          ELSE IF(I_EXT == -1) THEN                                 !
            NDP = NTHETA * NPHI * NE * 2                            !
          ELSE IF(I_EXT == 2) THEN                                  !
            NDP = NTHETA * NE                                       !
          END IF                                                    !
          NTT = N_INFILES * NDP                                     !
          IF(NTT > NDIM_M) GO TO 10                                 !
          IF((NTHETA > N_TH_M) .OR. (NPHI > N_PH_M)) GO TO 50       !
!
!   Auger electron diffraction case
!
        ELSE IF(SPECTRO == 'AED') THEN                              !
          IF((I_EXT_A == 0) .OR. (I_EXT_A == 1)) THEN               !
            NDP = NTHETA_A * NPHI_A * NE                            !
          ELSE IF(I_EXT_A  == -1) THEN                              !
            NDP = NTHETA_A * NPHI_A * NE * 2                        !
          ELSE IF(I_EXT_A == 2) THEN                                !
            NDP = NTHETA_A * NE                                     !
          END IF                                                    !
          NTT = N_INFILES * NDP                                     !
          IF(NTT > NDIM_M) GO TO 20                                 !
          IF((NTHETA_A > N_TH_M) .OR. (NPHI_A > N_PH_M)) GO TO 50   !
!
!   X-ray absorption case
!
        ELSE IF(SPECTRO == 'XAS') THEN                              !
          NDP = NE                                                  !
          NTT = N_INFILES * NDP                                     !
          IF(NTT > NDIM_M) GO TO 30                                 !
!
!   Auger Photoelectron coincidence spectroscopy case
!
        ELSE IF(SPECTRO == 'APC') THEN                              !
          IF((I_EXT == 0) .OR. (I_EXT == 1)) THEN                   !
            IF((I_EXT_A ==  0) .OR. (I_EXT_A == 1)) THEN            !
              NDP = NTHETA * NPHI * NE * NTHETA_A * NPHI_A          !
            ELSE IF(I_EXT_A == -1) THEN                             !
              NDP = NTHETA * NPHI * NE * NTHETA_A * NPHI_A * 2      !
            ELSE IF(I_EXT_A == 2) THEN                              !
              NDP = NTHETA * NPHI * NE * NTHETA_A                   !
            END IF                                                  !
          ELSE IF(I_EXT == -1) THEN                                 !
            IF((I_EXT_A == 0) .OR. (I_EXT_A == 1)) THEN             !
              NDP = NTHETA * NPHI * NE * 2 * NTHETA_A * NPHI_A      !
            ELSE IF(I_EXT_A == -1) THEN                             !
              NDP = NTHETA * NPHI * NE * 2 * NTHETA_A * NPHI_A * 2  !
            ELSE IF(I_EXT_A == 2) THEN                              !
              NDP = NTHETA * NPHI * NE * 2 * NTHETA_A               !
            END IF                                                  !
          ELSE IF(I_EXT == 2) THEN                                  !
            IF((I_EXT_A == 0) .OR. (I_EXT_A == 1)) THEN             !
              NDP = NTHETA * NE * NTHETA_A * NPHI_A                 !
            ELSE IF(I_EXT_A == -1) THEN                             !
              NDP = NTHETA * NE * NTHETA_A * NPHI_A * 2             !
            ELSE IF(I_EXT_A == 2) THEN                              !
              NDP = NTHETA * NE * NTHETA_A                          !
            END IF                                                  !
          END IF                                                    !
          NTT = N_INFILES * NDP                                     !
          IF(NTT > NDIM_M) GO TO 40                                 !
          IF((NTHETA > N_TH_M) .OR. (NPHI > N_PH_M))     GO TO 50   !
          IF((NTHETA_A > N_TH_M) .OR. (NPHI_A > N_PH_M)) GO TO 50   !
        END IF                                                      !
      END IF                                                        !
!
      GO TO 5                                                       !
!
!  Stops:
!
  10  WRITE(IUO1,11) NTT                                            !
      STOP                                                          !
  20  WRITE(IUO1,21) NTT                                            !
      STOP                                                          !
  30  WRITE(IUO1,31) NTT                                            !
      STOP                                                          !
  40  WRITE(IUO1,41) NTT                                            !
      STOP                                                          !
  50  WRITE(IUO1,51) MAX(NTHETA,NPHI,NTHETA_A,NPHI_A)               !
      STOP                                                          !
!
!  Formats:
!
  11  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',    &
             ' IN THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ', &
             'FOR THE TREAT_PED SUBROUTINE   >>>>>>>>>>',/,8X,         &
             '<<<<<<<<<<          NDIM_M BE AT LEAST ',I8,             &
             '         >>>>>>>>>>')
  21  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',    &
             ' IN THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ', &
             'FOR THE TREAT_AED SUBROUTINE   >>>>>>>>>>',/,8X,         &
             '<<<<<<<<<<          NDIM_M BE AT LEAST ',I8,             &
             '         >>>>>>>>>>')
  31  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',    &
             ' IN THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ', &
             'FOR THE TREAT_XAS SUBROUTINE   >>>>>>>>>>',/,8X,         &
             '<<<<<<<<<<          NDIM_M BE AT LEAST ',I8,             &
             '         >>>>>>>>>>')
  41  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF THE ARRAYS TOO SMALL',    &
             ' IN THE    >>>>>>>>>>',/,8X,'<<<<<<<<<<  INCLUDE FILE ', &
             'FOR THE TREAT_APC SUBROUTINE   >>>>>>>>>>',/,8X,         &
             '<<<<<<<<<<          NDIM_M BE AT LEAST ',I8,             &
             '         >>>>>>>>>>')
  51  FORMAT(//,8X,'<<<<<<<<<<  DIMENSION OF NTH_M OR NPH_M TOO SMALL',&
             'IN THE INCLUDE FILE - SHOULD BE AT LEAST ',I6,           &
             '  >>>>>>>>>>')
!
   5  RETURN                                                        !
!
      END SUBROUTINE STOP_TREAT
!
END MODULE TREATMENT_PED
