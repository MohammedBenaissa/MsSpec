!
!=======================================================================
!
MODULE ELECTRON_PHOTON
!
!  This module provides routines to compute the photon-electron interaction
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : LI_M,LINFMAX,N_GAUNT
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE CALC_MLFLI(JE,E_LIGHT,JEL,Q,MLFLI_E1,MLFLI_E2,MLFLI_E3, &
                            MLFLI_M1,MLFLI_M2)
!
!  This subroutine stores the photon-electron excitation matrix elements
!       into the MLFLI arrays for further use. These matrix elements are
!       computed for the x, y and z directions of the polarization and
!       of the direction of the light (JBASE1,JBASE2). They are stored
!       as a function of MI,LF and MF.
!
!
!  Input variables :
!
!                       JE        :  energy point
!                       E_LIGHT   :  energy of the photon beam (in eV)
!                       JEL       :  electron number --> 1 : IN  incoming
!                                                        2 : EX  excited
!                                                        3 : O1  1st outgoing
!                                                        4 : O2  2nd outgoing
!                                                        5 : O3  3rd outgoing
!                       Q         :  wave number of the electromagnetic field
!                                      in Rydberg atomic units
!
!  Output variables :
!
!                       MLFLI_E1  :  electric dipole matrix elements
!                       MLFLI_E2  :  electric quadrupole matrix elements
!                       MLFLI_E3  :  electric octupole matrix elements
!                       MLFLI_M1  :  magnetic dipole matrix elements
!                       MLFLI_M2  :  magnetic quadrupole + polar toroidal dipole matrix elements
!
!
!  Types of matrix elements calculated :
!
!              * I_C1 = 1 : single channel
!              * I_E1 = 1 : electric dipole E1
!              * I_E2 = 1 : electric quadrupole E2
!              * I_E3 = 1 : electric octupole E3
!              * I_M1 = 1 : magnetic dipole M1
!              * I_M2 = 1 : magnetic quadrupole M2 + polar toroidal dipole T1
!
!
!  Conventions for the storage index LR :
!
!              * LR = 1 : electric dipole channel     LF = LI-1
!              * LR = 2 : electric dipole channel     LF = LI+1
!
!              * LR = 3 : electric quadrupole channel LF = LI-2
!              * LR = 4 : electric quadrupole channel LF = LI
!              * LR = 5 : electric quadrupole channel LF = LI+2
!
!              * LR = 6 : magnetic dipole channel     LF = LI
!
!              * LR = 7 : electric quadrupole channel LF = LI-3
!              * LR = 8 : electric quadrupole channel LF = LI-1
!              * LR = 9 : electric quadrupole channel LF = LI+1
!              * LR = 10: electric quadrupole channel LF = LI+3
!
!              * LR = 11: magnetic quadrupole channel LF = LI-1
!              * LR = 12: magnetic quadrupole channel LF = LI+1
!
!
!
!
!  Author : D. Sébilleau
!
!                                     Last modified : 19 May 2021
!
!
      USE COMPLEX_NUMBERS,      ONLY : ZEROC
!
      USE DICHROISM
      USE EL_PH_CORRECTION
      USE EXP_TYPE
      USE INIT_L
      USE ONE_CHANNEL
      USE OPTICAL_ELTS
      USE RADINT_R
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  JE,JEL
!
      INTEGER                   ::  L_STORE(-3:3)
      INTEGER                   ::  MI,LF,MF,MD,LD
      INTEGER                   ::  JBASE1,JBASE2
      INTEGER                   ::  LA,MA,LR
!
      REAL (WP), INTENT(IN)     ::  E_LIGHT,Q
!
      COMPLEX (WP), INTENT(OUT) ::  MLFLI_E1(2,-LI_M:LI_M,-LI_M-1:LI_M+1,0:LI_M+1,3)
      COMPLEX (WP), INTENT(OUT) ::  MLFLI_E2(2,-LI_M:LI_M,-LI_M-2:LI_M+2,0:LI_M+2,3,3)
      COMPLEX (WP), INTENT(OUT) ::  MLFLI_E3(2,-LI_M:LI_M,-LI_M-3:LI_M+3,0:LI_M+3,3,3)
      COMPLEX (WP), INTENT(OUT) ::  MLFLI_M1(2,-LI_M:LI_M,-LI_M-1:LI_M+1,0:LI_M,3)
      COMPLEX (WP), INTENT(OUT) ::  MLFLI_M2(2,-LI_M:LI_M,-LI_M-2:LI_M+2,0:LI_M+1,3,3)
!
      COMPLEX (WP)              ::  RD
      COMPLEX (WP)              ::  MATRIX_E1(3,2),MATRIX_E2(3,3,2)
      COMPLEX (WP)              ::  MATRIX_E3(3,3,2),MATRIX_M1(3,2)
      COMPLEX (WP)              ::  MATRIX_M2(3,3,2)
      COMPLEX (WP)              ::  MAT_E1_E1(3),MAT_E2_E2(3,3)
      COMPLEX (WP)              ::  MAT_E3_E3(3,3),MAT_M1_M1(3)
      COMPLEX (WP)              ::  MAT_M2_M2(3,3),MAT_E1_E2(3,3)
      COMPLEX (WP)              ::  MAT_M1_E1(3)
!
!   Initialization of (LF-LI) storage index L_STORE
!     (order in which the 10 radial integrals are stored in phagen_scf.f)
!
      L_STORE(-3) = 7                                               !
      L_STORE(-2) = 3                                               !
      L_STORE(-1) = 1                                               !
      L_STORE(0)  = 4                                               !
      L_STORE(1)  = 2                                               !
      L_STORE(2)  = 5                                               !
      L_STORE(3)  = 10                                              !
!
!  Third order correction for E1
!
      IF(I_E1 == 2) THEN                                            !
        I_CORR = 1                                                  !
      ELSE                                                          !
        I_CORR = 0                                                  !
      END IF                                                        !
!
!  Initialization of the matrix elements MLFLI_Xn
!
      DO MI = - LI, LI                                              !
        DO LF = LF1, LF2, ISTEP_LF                                  !
          LD = ABS(LF - LI)                                         !
          DO MF = -LF, LF                                           !
            MD = ABS(MF - MI)                                       !
            DO JBASE1 = 1, 3                                        !
!
              IF((LD <= 1) .AND. (MD <= 1)) THEN                    !
                MLFLI_E1(1,MI,MF,LF,JBASE1) = ZEROC                 !
                MLFLI_E1(2,MI,MF,LF,JBASE1) = ZEROC                 !
                MLFLI_M1(1,MI,MF,LF,JBASE1) = ZEROC                 !
              END IF                                                !
!
              DO JBASE2 = 1, 3
!
                IF((LD <= 1) .AND. (MD <= 2)) THEN                  !
                  MLFLI_M2(1,MI,MF,LF,JBASE1,JBASE2) = ZEROC        !
                END IF                                              !
                IF((LD <= 2) .AND. (MD <= 2)) THEN                  !
                  MLFLI_E2(1,MI,MF,LF,JBASE1,JBASE2) = ZEROC        !
                END IF                                              !
                IF((LD <= 3) .AND. (MD <= 3)) THEN                  !
                  MLFLI_E3(1,MI,MF,LF,JBASE1,JBASE2) = ZEROC        !
                END IF                                              !
!
              END DO                                                !
            END DO                                                  !
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
!
!  Electric dipole case E1
!
      IF(I_E1 >= 1) THEN                                            !
        DO MI = - LI, LI                                            !
          DO LF = LF1, LF2, ISTEP_LF                                !
            LA = ABS(LF - LI)                                       !
            IF(LA == 1) THEN                                        ! E1 L-selection rule
              LR = L_STORE(LF-LI)                                   !
              RD = RHOR(JE,LR,1,JEL)                                !
              DO MF = - LF, LF                                      !
                MA = ABS(MF - MI)                                   !
                IF(MA > 1) GO TO 100                                ! E1 M-selection rule
                IF((SELRULE == '00000') .AND. (MF /= MI)) GO TO 100 !
!
                CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,       & !
                                  MATRIX_E1,MATRIX_E2,MATRIX_E3,  & !
                                  MATRIX_M1,MATRIX_M2)              !
                IF((SPECTRO == 'XAS') .OR. (SPECTRO == 'RES')) THEN !
                  CALL COUMAT_IRREG(MI,LF,MF,RD,LR,MAT_E1_E1,     & !
                                    MAT_E2_E2,MAT_E3_E3,MAT_M1_M1,& !
                                    MAT_M2_M2,MAT_E1_E2,MAT_M1_E1)  !
                END IF                                              !
!
                DO JBASE1 = 1, 3                                    !
                  MLFLI_E1(1,MI,MF,LF,JBASE1) = MATRIX_E1(JBASE1,1) !
                  IF(I_DICHR >= 1) THEN                             !
                   MLFLI_E1(2,MI,MF,LF,JBASE1) = MATRIX_E1(JBASE1,2)!
                  END IF                                            !
                END DO                                              !
!
 100            CONTINUE                                            !
              END DO                                                !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
!  Electric quadrupole case E2
!
      IF(I_E2 >= 1) THEN                                            !
        DO MI = - LI, LI                                            !
          DO LF = LF1, LF2, ISTEP_LF                                !
            LA= ABS(LF - LI)                                        !
            IF((LA == 0) .OR. (LA == 2)) THEN                       ! E2 L-selection rule
              LR = L_STORE(LF-LI)                                   !
              RD = RHOR(JE,LR,1,JEL)                                !
              DO MF = - LF, LF                                      !
                MA = ABS(MF - MI)                                   !
                IF(MA > 2) GO TO 200                                ! E2 M-selection rule
                IF((SELRULE == '00000').AND.(MF /=  MI)) GO TO 200   !
!
                CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,       & !
                                  MATRIX_E1,MATRIX_E2,MATRIX_E3,  & !
                                  MATRIX_M1,MATRIX_M2)              !
                IF((SPECTRO == 'XAS').OR.(SPECTRO == 'RES')) THEN   !
                  CALL COUMAT_IRREG(MI,LF,MF,RD,LR,MAT_E1_E1,     & !
                                    MAT_E2_E2,MAT_E3_E3,MAT_M1_M1,& !
                                    MAT_M2_M2,MAT_E1_E2,MAT_M1_E1)  !
                END IF                                              !
!
                DO JBASE1 = 1, 3                                    !
                  DO JBASE2 = 1, 3                                  !
                    MLFLI_E2(1,MI,MF,LF,JBASE1,JBASE2) =          & !
                                        MATRIX_E2(JBASE1,JBASE2,1)  !
                  END DO                                            !
                END DO                                              !
!
 200            CONTINUE                                            !
              END DO                                                !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
!  Electric octupole case E3
!
      IF(I_E3 == 1) THEN                                            !
        DO MI = - LI, LI                                            !
          DO LF = LF1, LF2, ISTEP_LF                                !
            LA = ABS(LF - LI)                                       !
            IF((LA == 1) .OR. (LA == 3)) THEN                       ! E3 L-selection rule
              LR = L_STORE(LF-LI)                                   !
              IF(LR == 1) LR = 8                                    ! LF=LI-1 for E3
              IF(LR == 2) LR = 9                                    ! LF=LI+1 for E3
              RD = RHOR(JE,LR,1,JEL)                                !
              DO MF = - LF, LF                                      !
                MA = ABS(MF - MI)                                   !
                IF(MA > 3) GO TO 300                                ! E3 M-selection rule
!
                CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,       & !
                                  MATRIX_E1,MATRIX_E2,MATRIX_E3,  & !
                                  MATRIX_M1,MATRIX_M2)              !
                IF((SPECTRO == 'XAS') .OR. (SPECTRO == 'RES')) THEN !
                  CALL COUMAT_IRREG(MI,LF,MF,RD,LR,MAT_E1_E1,     & !
                                    MAT_E2_E2,MAT_E3_E3,MAT_M1_M1,& !
                                    MAT_M2_M2,MAT_E1_E2,MAT_M1_E1)  !
                END IF                                              !
!
                DO JBASE1 = 1, 3                                    !
                  DO JBASE2 = 1 ,3                                  !
                    MLFLI_E3(1,MI,MF,LF,JBASE1,JBASE2) =          & !
                                        MATRIX_E3(JBASE1,JBASE2,1)  !
                  END DO                                            !
                END DO                                              !
!
 300            CONTINUE                                            !
              END DO                                                !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
!  Magnetic dipole case M1
!
      IF(I_M1 == 1) THEN                                            !
        DO MI = - LI, LI                                            !
          LF = LI                                                   ! M1 L-selection rule
          LR = 6                                                    !
          RD = RHOR(JE,LR,1,JEL)                                    !
          DO MF = - LF, LF                                          !
            MA = ABS(MF - MI)                                       !
            IF(MA > 1) GO TO 400                                    ! M1 M-selection rule
!
            CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,MATRIX_E1, & !
                              MATRIX_E2,MATRIX_E3,MATRIX_M1,      & !
                              MATRIX_M2)                            !
            IF((SPECTRO == 'XAS') .OR. (SPECTRO == 'RES')) THEN     !
              CALL COUMAT_IRREG(MI,LF,MF,RD,LR,MAT_E1_E1,         & !
                                MAT_E2_E2,MAT_E3_E3,MAT_M1_M1,    & !
                                MAT_M2_M2,MAT_E1_E2,MAT_M1_E1)      !
            END IF                                                  !
!
            DO JBASE1 = 1, 3                                        !
              MLFLI_M1(1,MI,MF,LF,JBASE1) = MATRIX_M1(JBASE1,1)     !
            END DO                                                  !
!
 400        CONTINUE                                                !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
!  Magnetic quadrupole case and polar toroidal case M2 + T1
!
      IF(I_M2 == 1) THEN                                            !
        DO MI = - LI, LI                                            !
          DO LF = LF1, LF2, ISTEP_LF                                !
            LA = ABS(LF - LI)                                       !
            IF(LA == 1) THEN                                        ! M2 L-selection rule
              LR = L_STORE(LF-LI)                                   !
              RD = RHOR(JE,LR,1,JEL)                                ! same radial integral as E1
              IF(LR == 1) LR = 11                                   ! LF=LI-1 for M2
              IF(LR == 2) LR = 12                                   ! LF=LI+1 for M2
              DO MF = - LF, LF                                      !
                MA = ABS(MF - MI)                                   !
                IF(MA > 2) GO TO 500                                ! M2 M-selection rule
!
                CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,       & !
                                  MATRIX_E1,MATRIX_E2,MATRIX_E3,  & !
                                  MATRIX_M1,MATRIX_M2)              !
                IF((SPECTRO == 'XAS') .OR. (SPECTRO == 'RES')) THEN !
                  CALL COUMAT_IRREG(MI,LF,MF,RD,LR,MAT_E1_E1,     & !
                                    MAT_E2_E2,MAT_E3_E3,MAT_M1_M1,& !
                                    MAT_M2_M2,MAT_E1_E2,MAT_M1_E1)  !
                END IF                                              !
!
                DO JBASE1 = 1 ,3                                    !
                  DO JBASE2 = 1, 3                                  !
                    MLFLI_M2(1,MI,MF,LF,JBASE1,JBASE2) =          & !
                                        MATRIX_M2(JBASE1,JBASE2,1)  !
                  END DO                                            !
                END DO                                              !
!
 500            CONTINUE                                            !
              END DO                                                !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
      END IF                                                        !
!
!  Single channel cases
!                                                                   !
      IF(I_C1 == 1) THEN
!
        IF(TYPCHAN == 'E1') THEN                                    !
!
!.........Electric dipole single channel
!
          IF(CHANNEL == '-1') THEN                                  !
            LF = LI - 1                                             !
            LR = 1                                                  !
          ELSE IF(CHANNEL == '+1') THEN                             !
            LF = LI + 1                                             !
            LR = 2                                                  !
          END IF                                                    !
          RD = RHOR(JE,LR,1,JEL)                                    !
          DO MI = - LI, LI                                          !
            DO MF = - LF, LF                                        !
              MA = ABS(MF - MI)                                     !
              IF(MA > 1) GO TO 600                                  !
!
              CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,         & !
                                MATRIX_E1,MATRIX_E2,MATRIX_E3,    & !
                                MATRIX_M1,MATRIX_M2)                !

              DO JBASE1 = 1, 3                                      !
                MLFLI_E1(1,MI,MF,LF,JBASE1) = MATRIX_E1(JBASE1,1)   !
              END DO                                                !
!
 600          CONTINUE                                              !
            END DO                                                  !
          END DO                                                    !
!
        ELSE IF(TYPCHAN == 'E2') THEN                               !
!
!.........Electric quadrupole single channel
!
          IF(CHANNEL == '-2') THEN                                  !
            LF = LI - 2                                             !
            LR = 3                                                  !
          ELSE IF(CHANNEL == '+0') THEN                             !
            LF = LI                                                 !
            LR = 4                                                  !
          ELSE IF(CHANNEL == '+2') THEN                             !
            LF = LI + 1                                             !
            LR = 5                                                  !
          END IF                                                    !
          RD = RHOR(JE,LR,1,JEL)                                    !
          DO MI = - LI, LI                                          !
            DO MF = - LF, LF                                        !
              MA = ABS(MF - MI)                                     !
              IF(MA > 2) GO TO 700                                  !
!
              CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,         & !
                                MATRIX_E1,MATRIX_E2,MATRIX_E3,    & !
                                MATRIX_M1,MATRIX_M2)                !
              DO JBASE1 = 1, 3                                      !
                DO JBASE2 = 1, 3                                    !
                  MLFLI_E2(1,MI,MF,LF,JBASE1,JBASE2) =            & !
                                      MATRIX_E2(JBASE1,JBASE2,1)    !
                END DO                                              !
              END DO                                                !
!
 700          CONTINUE                                              !
            END DO                                                  !
          END DO                                                    !
!
        ELSE IF(TYPCHAN == 'E3') THEN                               !
!
!.........Electric octupole single channel
!
          IF(CHANNEL == '-3') THEN                                  !
            LF = LI - 3                                             !
            LR = 7                                                  !
          ELSE IF(CHANNEL == '-1') THEN                             !
            LF = LI - 1                                             !
            LR = 8                                                  !
          ELSE IF(CHANNEL == '+1') THEN                             !
            LF = LI + 1                                             !
            LR = 9                                                  !
          ELSE IF(CHANNEL == '+3') THEN                             !
            LF = LI + 3                                             !
            LR = 10                                                 !
          END IF                                                    !
          RD = RHOR(JE,LR,1,JEL)                                    !
          DO MI = - LI ,LI                                          !
            DO MF = - LF, LF                                        !
              MA = ABS(MF - MI)                                     !
              IF(MA > 3) GO TO 800                                  !
!
              CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,         & !
                                MATRIX_E1,MATRIX_E2,MATRIX_E3,    & !
                                MATRIX_M1,MATRIX_M2)                !
              DO JBASE1 = 1, 3                                      !
                DO JBASE2 = 1, 3                                    !
                  MLFLI_E3(1,MI,MF,LF,JBASE1,JBASE2) =            & !
                                      MATRIX_E3(JBASE1,JBASE2,1)    !
                END DO                                              !
              END DO                                                !
!
 800          CONTINUE                                              !
            END DO                                                  !
          END DO                                                    !
!
        ELSE IF(TYPCHAN == 'M1') THEN                               !
!
!.........Magnetic dipole single channel
!
          IF(CHANNEL == '+0') THEN                                  !
            LF = LI                                                 !
            LR = 6                                                  !
          END IF                                                    !
          RD = RHOR(JE,LR,1,JEL)                                    !
          DO MI = - LI, LI                                          !
            DO MF = - LF, LF                                        !
              MA = ABS(MF - MI)                                     !
              IF(MA > 1) GO TO 900                                  !
!
              CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,         & !
                                MATRIX_E1,MATRIX_E2,MATRIX_E3,    & !
                                MATRIX_M1,MATRIX_M2)                !
              DO JBASE1 = 1, 3                                      !
                MLFLI_M1(1,MI,MF,LF,JBASE1) = MATRIX_M1(JBASE1,1)   !
              END DO                                                !
!
 900          CONTINUE                                              !
            END DO                                                  !
          END DO                                                    !
!
        ELSE IF(TYPCHAN == 'M2') THEN                               !
!
!.........Magnetic quadrupole and polar toroidal single channel
!
          IF(CHANNEL == '-1') THEN                                  !
            LF = LI - 1                                             !
            LR = 1                                                  !
          ELSE IF(CHANNEL == '+1') THEN                             !
            LF = LI + 1                                             !
            LR = 2                                                  !
          END IF                                                    !
          RD = RHOR(JE,LR,1,JEL)                                    !
          DO MI = - LI, LI                                          !
            DO MF = - LF, LF                                        !
              MA = ABS(MF - MI)                                     !
              IF(MA > 2) GO TO 1000                                 !
!
              CALL COUMAT_REGUL(MI,LF,MF,RD,LR,Q,E_LIGHT,         & !
                                MATRIX_E1,MATRIX_E2,MATRIX_E3,    & !
                                MATRIX_M1,MATRIX_M2)                !

              DO JBASE1 = 1, 3                                      !
                DO JBASE2 = 1, 3                                    !
                  MLFLI_M2(1,MI,MF,LF,JBASE1,JBASE2) =            & !
                                MATRIX_M2(JBASE1,JBASE2,1)          !
                END DO                                              !
              END DO                                                !
!
 1000         CONTINUE                                              !
            END DO                                                  !
          END DO                                                    !
        END IF                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE CALC_MLFLI
!
!=======================================================================
!
      SUBROUTINE COUMAT_REGUL(MI,LF,MF,RADIAL,LR,Q,E_LIGHT,MATRIX_E1, &
                              MATRIX_E2,MATRIX_E3,MATRIX_M1,MATRIX_M2)
!
!  This routine calculates the spin-independent regular optical matrix
!   elements for dipolar, quadrupolar and octupolar  excitations. They
!   are calculated separately for x, y and z directions of the polarization
!   (J_PO) and of the propagation of the photon field (J_LI), for a given
!   type of polarization (JPOL).
!
!  They are stored in MATRIX_XX(J_PO,J_LI,JPOL)
!
!
!  Input variables :
!
!                       MI        :  azimuthal quantum number of the core state
!                       LF        :  quantum number of the excited state
!                       MF        :  azimuthal quantum number of the excited state
!                       RADIAL    :  regular radial integral of the corresponding channel
!                       LR        :  storage index of the radial integral
!                       Q         :  photon wave number in Rydberg atomic units
!                       E_LIGHT   :  photon energy (in eV)
!
!  Output variables :
!
!                       MATRIX_E1 :  electric dipole matrix element
!                       MATRIX_E2 :  electric quadrupole matrix element
!                       MATRIX_E3 :  electric octupole matrix element
!                       MATRIX_M1 :  magnetic dipole matrix element
!                       MATRIX_M2 :  magnetic quadrupole matrix element +
!                                      polar toroidal dipole matrix element
!
!
!
!  Here, the conventions are :
!
!             IPOL=1  : linearly polarized light
!             IPOL=2  : circularly polarized light
!
!             JPOL=1  : +/x polarization for circular/linear light
!             JPOL=2  : -/y polarization for circular/linear light
!
!  When I_DICHR=0, J_PO = 1,2 and 3 correspond respectively to the x,y
!                         and z directions for the linear polarization
!                         (except for magnetic dipole where it corresponds
!                          to the x,y and z directions of the cluster)
!                 J_LI = 1,2 and 3 correspond respectively to the x,y
!                         and z directions for the position of the light
!
!  When I_DICHR=1, J_PO= 1,2 and 3 correspond respectively to the x,y
!                         and z directions for the position of the light
!
!
!  Types of matrix elements calculated:
!
!              * I_C1 = 1 : single channel
!              * I_E1 = 1 : electric dipole E1
!              * I_E2 = 1 : electric quadrupole E2
!              * I_E3 = 1 : electric octupole E3
!              * I_M1 = 1 : magnetic dipole M1
!              * I_M2 = 1 : magnetic quadrupole M2 + polar toroidal dipole T1
!
!
!  Conventions for the storage index LR :
!
!              * LR = 1 : electric dipole channel     LF = LI-1
!              * LR = 2 : electric dipole channel     LF = LI+1
!
!              * LR = 3 : electric quadrupole channel LF = LI-2
!              * LR = 4 : electric quadrupole channel LF = LI
!              * LR = 5 : electric quadrupole channel LF = LI+2
!
!              * LR = 6 : magnetic dipole channel     LF = LI
!
!              * LR = 7 : electric octupole channel LF = LI-3
!              * LR = 8 : electric octupole channel LF = LI-1
!              * LR = 9 : electric octupole channel LF = LI+1
!              * LR = 10: electric octupole channel LF = LI+3
!
!              * LR = 11: magnetic quadrupole channel LF = LI-1
!              * LR = 12: magnetic quadrupole channel LF = LI+1
!
!
!  Note: the LIGHT x POLARIZATION cross product is calculated outside this subroutine.
!        The matrix elements here do not take it into account in M1 and M2
!
!
!  All output regular matrix elements MATRIX_Xn are written as
!
!                  MATRIX_Xn = COEF(Xn) * RADIAL * ANG(Xn)
!
!  COEF(Xn) and ANG(Xn) are stored in the common block /COEFMAT/ for further use
!    in the calculation of the irregular matrix elements
!
!
!
!  Author : D. Sébilleau
!
!                                     Last modified : 20 May 2021
!
!
      USE REAL_NUMBERS,         ONLY : ZERO,ONE,THREE,SIX,HALF
      USE COMPLEX_NUMBERS
      USE CONSTANTS_P2,         ONLY : ALPHA
      USE ENE_CHANGE,           ONLY : RYD
!
      USE EL_PH_CORRECTION
      USE ANGULAR_MOMENTUM,     ONLY : GAUNT
      USE INIT_L
      USE MATRIX_COEF
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  MI,LF,MF,LR
!
      INTEGER                   ::  INDI,INDF
      INTEGER                   ::  J_PO,JPOL,J_LI,J_TE
      INTEGER                   ::  I,J
      INTEGER                   ::  L_MIN,L_MAX
      INTEGER                   ::  LB_MIN,LB_MAX,LB_LIN
      INTEGER                   ::  L1_MIN,L1_MAX
      INTEGER                   ::  MP,MPP,M,L1,M1,LB,MB,MD
!
      REAL (WP), INTENT(IN)     ::  E_LIGHT,Q
!
      REAL (WP)                 ::  FPI3,KOE1,KOE2
      REAL (WP)                 ::  GG,CORR
      REAL (WP)                 ::  GNT1(0:N_GAUNT),GNT2(0:N_GAUNT),GNT3(0:N_GAUNT)
      REAL (WP)                 ::  A_MINUS,A_ZERO,A_PLUS
!
      REAL (WP)                 ::  FLOAT,SQRT
!
      COMPLEX (WP), INTENT(IN)  ::  RADIAL
      COMPLEX (WP), INTENT(OUT) ::  MATRIX_E1(3,2)
      COMPLEX (WP), INTENT(OUT) ::  MATRIX_E2(3,3,2)
      COMPLEX (WP), INTENT(OUT) ::  MATRIX_E3(3,3,2)
      COMPLEX (WP), INTENT(OUT) ::  MATRIX_M1(3,2)
      COMPLEX (WP), INTENT(OUT) ::  MATRIX_M2(3,3,2)
!
      COMPLEX (WP)              ::  PROD
      COMPLEX (WP)              ::  SUM_MM,SUM_M0,SUM_MP,SUM_M
      COMPLEX (WP)              ::  YLM1(3,-1:1),YLM2(3,-1:1)
!
      COMPLEX (WP)              ::  CMPLX,CONJG
!
      DATA FPI3  / 4.1887902047863909846168578443726705122628925325E0_WP /  ! 4 * pi / 3
      DATA KOE1  / 0.3454941494713354792652446460318896831393773703E0_WP /  ! 1/2 * sqrt( 3 / 2 pi)
      DATA KOE2  / 0.4886025119029199215863846228383470045758856082E0_WP /  ! 1/2 * sqrt( 3 / pi)
!
!  Linear indices for (LI,MI) and (LF,MF)
!
      INDI = LI * LI + LI + MI + 1                                  !
      INDF = LF * LF + LF + MF + 1                                  !
!
!  Third order correction term for dipole E1 term
!
      CORR = ONE                                                    !
!
      IF(I_CORR == 1) THEN                                          !
        CORR = ONE + E_LIGHT * ALPHA * ALPHA / (12.0E0_WP * RYD)    ! Dimensionless correction term
      END IF                                                        !
!
!  Initialisation of all the matrix elements                        ! J_PO : directions of polarization ! 1 = x
!                                                                   ! J_LI : directions of light        ! 2 = y
      DO J_PO = 1, 3                                                                                    ! 3 = z
        DO JPOL = 1 ,2                                              ! JPOL : polarization type
          MATRIX_E1(J_PO,JPOL) = ZEROC                              !
          MATRIX_M1(J_PO,JPOL) = ZEROC                              ! in this case, just x, y and z
          DO J_LI = 1 ,3                                            !
            MATRIX_E2(J_PO,J_LI,JPOL) = ZEROC                       !
            MATRIX_E3(J_PO,J_LI,JPOL) = ZEROC                       !
            MATRIX_M2(J_PO,J_LI,JPOL) = ZEROC                       !
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
!
      DO J_PO = 1, 3                                                !
        DO J_LI = 1 ,3                                              !
         DO J_TE = 1, 5                                             ! 1 : E1, 2 : E2, 3 : E3, 4 : M1, 5 : M2+T1
           DO I = 1, LINIMAX                                        ! (LI,MI) storage
             DO J = 1, LINFMAX                                      ! (LF,MF) storage
               ANG(J_TE,J_PO,J_LI,I,J) = ZEROC                      ! angular component of matrix element
             END DO                                                 ! stored for re-use in irregular terms
            END DO                                                  ! in module MATRIX_COEF
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
!
!  Spherical harmonics Y_1^m along the basis directions: YLMn(JDIR,M)
!
!  YLM1 : spherical harmonics for the (x,y,z) directions J_PO
!         of the polarization
!
!  YLM2 : spherical harmonics for the (x,y,z) directions J_LI
!         of the direction of the light
!
      YLM1(1,-1) =   CMPLX(KOE1,KIND=WP)                            !
      YLM1(1,0)  =   ZEROC                                          ! x direction
      YLM1(1,1)  = - CMPLX(KOE1,KIND=WP)                            !
!
      YLM1(2,-1) = - IC * CMPLX(KOE1,KIND=WP)                       !
      YLM1(2,0)  =   ZEROC                                          ! y direction
      YLM1(2,1)  = - IC * CMPLX(KOE1,KIND=WP)                       !
!
      YLM1(3,-1) =   ZEROC                                          !
      YLM1(3,0)  =   CMPLX(KOE2,KIND=WP)                            ! z direction
      YLM1(3,1)  =   ZEROC                                          !
!
      YLM2(1,-1) =   CMPLX(KOE1,KIND=WP)                            !
      YLM2(1,0)  =   ZEROC                                          ! x direction
      YLM2(1,1)  = - CMPLX(KOE1,KIND=WP)                            !
!
      YLM2(2,-1) = - IC * CMPLX(KOE1,KIND=WP)                       !
      YLM2(2,0)  =   ZEROC                                          ! y direction
      YLM2(2,1)  = - IC * CMPLX(KOE1,KIND=WP)                       !
!
      YLM2(3,-1) =   ZEROC                                          !
      YLM2(3,0)  =   CMPLX(KOE2,KIND=WP)                            ! z direction
      YLM2(3,1)  =   ZEROC                                          !
!
!
!  Storage of the coefficients of the different matrix elements COEF(JJ)
!          with:
!                  JJ = 1 --> E1
!                  JJ = 2 --> E2
!                  JJ = 3 --> E3
!                  JJ = 4 --> M1
!                  JJ = 5 --> M2 + T1
!
!  They are stored in the module MATRIX_COEF for further use
!    in the calculation of the irregular matrix elements
!    by subroutine COUMAT_IRREG
!
      COEF(1) =   CMPLX(FPI3 * CORR,KIND=WP)                        ! E1 (with/without correction term)
      COEF(2) =   CMPLX(FPI3 * FPI3 * Q * HALF,KIND=WP) * IC        ! E2
      COEF(3) = - CMPLX(FPI3 * FPI3 * FPI3 * Q * Q / SIX,KIND=WP)   ! E3
      COEF(4) =   CMPLX(ALPHA * HALF,KIND=WP)                       ! M1
      COEF(5) =   CMPLX(FPI3 * Q * ALPHA / THREE,KIND=WP) * IC      ! M2 + T1
!
      MD = MF - MI                                                  ! for the m-selection rule
!
!    Electric dipole case : LR = 1-2
!
      IF(LR <= 2) THEN                                              !
!
        IF(ABS(MD) <= 1) THEN                                       ! selection rule on MF
!
          CALL GAUNT(LI,MI,LF,MF,GNT1)                              !
          PROD = COEF(1) * RADIAL * GNT1(1)                         !
!
          DO J_PO = 1, 3                                            !
            MATRIX_E1(J_PO,1)       = PROD * CONJG(YLM1(J_PO,MD))   !
            ANG(1,J_PO,1,INDI,INDF) = CONJG(YLM1(J_PO,MD)) * GNT1(1)!
          END DO                                                    !
!
        END IF                                                      !
!
!     Electric quadrupole case : LR = 3-5
!
      ELSE IF((LR > 2) .AND. (LR <= 5)) THEN                        !
!
        IF(ABS(MD) <= 2) THEN                                       ! selection rule on MF
!
          L_MIN = MAX(0,ABS(LF-LI))                                 !
          L_MAX = MIN(2,LF+LI)                                      !
!
          CALL GAUNT(LI,MI,LF,MF,GNT1)                              !
          PROD = COEF(2) * RADIAL                                   !
!
          DO J_PO = 1, 3                                            !
            DO J_LI = 1, 3                                          !
!
              SUM_M = ZEROC                                         !
              DO MP = -1, 1                                         !
                DO M = -1, 1                                        !
!
                  MB = M - MP                                       !
                  IF(MB /= MD) GO TO 10                             !
                  CALL GAUNT(1,MP,1,M,GNT2)                         !
!
                  GG = ZERO                                         !
                  DO LB = L_MIN, L_MAX, 2                           !
                    GG = GG + GNT1(LB) * GNT2(LB)                   !
                  END DO
!
                  SUM_M = SUM_M + GG * CONJG(YLM1(J_PO,M)) *      & !
                                       YLM2(J_LI,MP)                !
!
  10              CONTINUE                                          !
!
                END DO                                              !
              END DO                                                !
!
              MATRIX_E2(J_PO,J_LI,1)     = PROD * SUM_M             !
              ANG(2,J_PO,J_LI,INDI,INDF) = SUM_M                    !
!
            END DO                                                  !
          END DO                                                    !
!
        END IF                                                      !
!
!     Magnetic dipole case : LR = 6                                 ! does not take into account
!                                                                   ! the LIGHT x POLARIZATION cross product
      ELSE IF(LR == 6) THEN                                         !
!
        IF(ABS(MD) <= 1) THEN                                       ! selection rule on MF
!
          IF(LF == LI) THEN                                         ! selection rule on LF
!
            A_PLUS  = SQRT(FLOAT((LI-MI)*(LI+MI+1))) * HALF         !
            A_ZERO  = FLOAT(MI)                                     !
            A_MINUS = SQRT(FLOAT((LI+MI)*(LI-MI+1))) * HALF         !
!
            PROD    = COEF(4) * RADIAL                              !
!
            IF(MF == (MI - 1)) THEN                                 !
              MATRIX_M1(1,1)       =   PROD * CMPLX(A_MINUS)        !
              MATRIX_M1(2,1)       = - PROD * CMPLX(A_MINUS) * IC   !
              ANG(4,1,1,INDI,INDF) =   CMPLX(A_MINUS)               !
              ANG(4,2,1,INDI,INDF) = - CMPLX(A_MINUS) * IC          !
            ELSE IF(MF == MI) THEN                                  !
              MATRIX_M1(3,1)       =   PROD * CMPLX(A_ZERO)         !
              ANG(4,3,1,INDI,INDF) =   CMPLX(A_ZERO)                !
            ELSE IF(MF == (MI + 1)) THEN                            !
              MATRIX_M1(1,1)       =   PROD * CMPLX(A_PLUS)         !
              MATRIX_M1(2,1)       =   PROD * CMPLX(A_PLUS) * IC    !
              ANG(4,1,1,INDI,INDF) =   CMPLX(A_PLUS)                !
              ANG(4,2,1,INDI,INDF) =   CMPLX(A_PLUS) * IC           !
            END IF                                                  !
!
          END IF                                                    !
!
        END IF                                                      !
!
!     Electric octupole case : LR = 7-10
!
      ELSE IF((LR > 6) .AND. (LR <= 10)) THEN                       !
!
        IF(ABS(MD) <= 3) THEN                                       ! selection rule on MF
!
          CALL GAUNT(LI,MI,LF,MF,GNT1)                              !
          PROD = COEF(3) * RADIAL                                   !
!
          MB     = MF - MI                                          !
          LB_MIN = ABS(LF - LI)                                     !
          LB_MAX = LF + LI                                          !
!
          DO J_PO = 1, 3                                            !
            DO J_LI = 1 ,3                                          !
!
              SUM_M = ZEROC                                         !
              DO M = - 1, 1                                         !
                DO MP = -1, 1                                       !
                  M1  = M - MP                                      !
                  MPP = MB - M1                                     !
                  IF(ABS(MPP) > 1) GO TO 20                         !
                  CALL GAUNT(1,MP,1,M,GNT2)                         !
!
                  GG = ZERO                                         !
                  DO LB = LB_MIN, LB_MAX, 2                         !
                    CALL GAUNT(1,MPP,LB,MB,GNT3)                    !
                    L1_MIN = MAX0(0,ABS(LB-1))                      !
                    L1_MAX = MIN0(2,LB+1)                           !
                    DO L1 = L1_MIN, L1_MAX, 2                       !
!
                      GG = GG + GNT1(LB) * GNT2(L1) * GNT3(L1)      !
!
                    END DO                                          !
                  END DO                                            !
!
                  SUM_M = SUM_M + GG * CONJG(YLM1(J_PO,M)) *      & !
                          YLM2(J_LI,MP) * CONJG(YLM2(J_LI,MPP))     !
!
  20              CONTINUE                                          !
!
                END DO                                              !
              END DO                                                !
!
              MATRIX_E3(J_PO,J_LI,1)     = PROD * SUM_M             !
              ANG(3,J_PO,J_LI,INDI,INDF) = SUM_M                    !
!
            END DO                                                  !
          END DO                                                    !
!
        END IF                                                      !
!
!     Magnetic quadrupole + polar toroidal case : LR = 11-12  ! does not take into account
!                                                             ! the LIGHT x POLARIZATION cross product
!
      ELSE IF((LR > 10) .AND. (LR <= 12)) THEN                      !
!
        IF(ABS(MD) <= 2) THEN                         ! selection rule on MF
!
          IF(ABS(LF-LI) == 1) THEN                    ! selection rule on LF
!
            A_PLUS  = SQRT(FLOAT((LI - MI) * (LI + MI + 1))) * HALF !
            A_ZERO  = FLOAT(MI)                                     !
            A_MINUS = SQRT(FLOAT((LI + MI) * (LI - MI + 1))) * HALF !
!
            PROD = COEF(5) * RADIAL                                 !
!
            CALL GAUNT(LI,MI-1,LF,MF,GNT1)                          !
            CALL GAUNT(LI,MI,LF,MF,GNT2)                            !
            CALL GAUNT(LI,MI+1,LF,MF,GNT3)                          !
!
            DO J_LI = 1, 3                                          !
!
              SUM_MM = ZEROC                                        !
              SUM_M0 = ZEROC                                        !
              SUM_MP = ZEROC                                        !
!
              DO M = -1, 1                                          !
!
                IF(M == (MD+1)) THEN                                ! -1 spherical tensor term
                  SUM_MM = SUM_MM + CONJG(YLM1(J_LI,M)) * GNT1(1)   !
                ELSE IF(M == MD) THEN                               !  0 spherical tensor term
                  SUM_M0 = SUM_M0 + CONJG(YLM1(J_LI,M)) * GNT2(1)   !
                ELSE IF(M == (MD-1)) THEN                           ! +1 spherical tensor term
                  SUM_MP = SUM_MP + CONJG(YLM1(J_LI,M)) * GNT3(1)   !
                END IF                                              !
!
              END DO                                                !
!
              MATRIX_M2(1,J_LI,1) = PROD * (                      & !
                                    CMPLX(A_PLUS) * SUM_MP +      & !
                                    CMPLX(A_MINUS)  *SUM_MM       & !
                                           )                        !
              MATRIX_M2(2,J_LI,1) = PROD * (                      & !
                                    CMPLX(A_PLUS) * SUM_MP -      & !
                                    CMPLX(A_MINUS) * SUM_MM       & !
                                           ) * IC                   !
              MATRIX_M2(3,J_LI,1) = PROD * CMPLX(A_ZERO) * SUM_M0   !
!
              ANG(5,1,J_LI,INDI,INDF) = CMPLX(A_PLUS)  * SUM_MP + & !
                                        CMPLX(A_MINUS) * SUM_MM     !
              ANG(5,2,J_LI,INDI,INDF) = ( CMPLX(A_PLUS)* SUM_MP - & !
                                          CMPLX(A_MINUS)*SUM_MM   & !
                                        ) * IC                      !
              ANG(5,3,J_LI,INDI,INDF) = CMPLX(A_ZERO) * SUM_M0      !
!
            END DO                                                  !
!
          END IF                                                    !
!
        END IF                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE COUMAT_REGUL
!
!=======================================================================
!
      SUBROUTINE COUMAT_IRREG(MI,LF,MF,RADIAL,LR,MAT_E1_E1,MAT_E2_E2,  &
                              MAT_E3_E3,MAT_M1_M1,MAT_M2_M2,MAT_E1_E2, &
                              MAT_M1_E1)
!
!  This routine calculates the spin-independent irregular optical matrix
!   elements for dipolar and quadrupolar excitations. They are calculated
!   separately for x, y and z directions of the polarization (J_PO) and
!   of the propagation of the photon field (J_LI), for a given type of
!   polarization (JPOL).
!
!  They are stored in MAT_XX_YY(J_PO,J_LI,JPOL)
!
!
!  Input variables :
!
!                       MI        :  azimuthal quantum number of the core state
!                       LF        :  quantum number of the excited state
!                       MF        :  azimuthal quantum number of the excited state
!                       RADIAL    :  irregular radial integral of the corresponding channel
!                       LR        :  storage index of the radial integral
!
!  Output variables :
!                    1) diagonal terms:
!
!                       MAT_E1_E1 :  diagonal electric dipole matrix element
!                       MAT_E2_E2 :  diagonal electric quadrupole matrix element
!                       MAT_E3_E3 :  diagonal electric octupole matrix element
!                       MAT_M1_M1 :  diagonal magnetic dipole matrix element
!                       MAT_M2_M2 :  diagonal magnetic quadrupole matrix element +
!                                      polar toroidal dipole matrix element
!
!                    2) non diagonal terms:
!
!                       MAT_E1_E2 :  electric dipole - electric quadrupole matrix element
!                       MAT_M1_E1 :  magnetic dipole - electric dipole matrix element
!
!
!
!  Here, the conventions are :
!
!             IPOL=1  : linearly polarized light
!             IPOL=2  : circularly polarized light
!
!             JPOL=1  : +/x polarization for circular/linear light
!             JPOL=2  : -/y polarization for circular/linear light
!
!  When I_DICHR=0, J_PO = 1,2 and 3 correspond respectively to the x,y
!                         and z directions for the linear polarization
!                         (except for magnetic dipole where it corresponds
!                          to the x,y and z directions of the cluster)
!                 J_LI = 1,2 and 3 correspond respectively to the x,y
!                         and z directions for the position of the light
!
!  When I_DICHR=1, J_PO = 1,2 and 3 correspond respectively to the x,y
!                         and z directions for the position of the light
!
!
!  Types of matrix elements calculated:
!
!              * I_C1 = 1 : single channel
!              * I_E1 = 1 : electric dipole E1
!              * I_E2 = 1 : electric quadrupole E2
!              * I_E3 = 1 : electric octupole E3
!              * I_M1 = 1 : magnetic dipole M1
!              * I_M2 = 1 : magnetic quadrupole M2 + polar toroidal dipole T1
!
!
!  Conventions for the storage index LR :
!
!              * LR = 1 : electric dipole channel     LF = LI-1
!              * LR = 2 : electric dipole channel     LF = LI+1
!
!              * LR = 3 : electric quadrupole channel LF = LI-2
!              * LR = 4 : electric quadrupole channel LF = LI
!              * LR = 5 : electric quadrupole channel LF = LI+2
!
!              * LR = 6 : magnetic dipole channel     LF = LI
!
!              * LR = 7 : electric octupole channel LF = LI-3
!              * LR = 8 : electric octupole channel LF = LI-1
!              * LR = 9 : electric octupole channel LF = LI+1
!              * LR = 10: electric octupole channel LF = LI+3
!
!              * LR = 11: magnetic quadrupole channel LF = LI-1
!              * LR = 12: magnetic quadrupole channel LF = LI+1
!
!
!  Note: the LIGHT x POLARIZATION cross product is calculated outside this subroutine.
!        The matrix elements here do not take it into account
!
!
!
!
!  Author : D. Sébilleau
!
!                                     Last modified : 19 May 2021
!
!
      USE COMPLEX_NUMBERS,      ONLY : ZEROC
!
      USE INIT_L
      USE MATRIX_COEF
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  MI,LF,MF,LR
!
      INTEGER                   ::  INDI,INDF
      INTEGER                   ::  J_PO,J_LI,JDIR
      INTEGER                   ::  MD
!
      COMPLEX (WP), INTENT(IN)  ::  RADIAL
      COMPLEX (WP), INTENT(OUT) ::  MAT_E1_E1(3),MAT_E2_E2(3,3)
      COMPLEX (WP), INTENT(OUT) ::  MAT_E3_E3(3,3)
      COMPLEX (WP), INTENT(OUT) ::  MAT_M1_M1(3),MAT_M2_M2(3,3)
      COMPLEX (WP), INTENT(OUT) ::  MAT_E1_E2(3,3),MAT_M1_E1(3)
!
      COMPLEX (WP)              ::  ANG_E1(3),ANG_E2(3,3)
      COMPLEX (WP)              ::  ANG_E3(3,3),ANG_M1(3)
      COMPLEX (WP)              ::  ANG_M2(3,3)
!
      COMPLEX (WP)              ::  CONJG
!
!  Linear indices for (LI,MI) and (LF,MF)
!
      INDI = LI * LI + LI + MI + 1                                  !
      INDF = LF * LF + LF + MF + 1                                  !
!
!  Initialisation of angular the matrix elements
!
      DO J_PO = 1, 3                                                ! x, y, z directions of the polarization
        ANG_E1(J_PO) = ZEROC                                        !
        ANG_M1(J_PO) = ZEROC                                        !
        DO J_LI = 1, 3                                              ! x, y, z  directions of the light
          ANG_E2(J_PO,J_LI) = ZEROC                                 !
          ANG_E3(J_PO,J_LI) = ZEROC                                 !
          ANG_M2(J_PO,J_LI) = ZEROC                                 !
        END DO                                                      !
      END DO                                                        !
!
!  Calculation of the angular parts of the E1,E2,E3,M1 and M2-T1 matrix elements
!
      MD = MF - MI                                                  ! for the m-selection rule
!
!  Electric dipole case : LR = 1-2
!
      IF(LR <= 2) THEN                                              !
!
        IF(ABS(MD) <= 1) THEN                                       ! selection rule on MF
!
          DO J_PO = 1, 3                                            !
            ANG_E1(J_PO) = ANG(1,J_PO,1,INDI,INDF)                  !
          END DO                                                    !
!
        END IF                                                      !
!
!     Electric quadrupole case : LR = 3-5
!
      ELSE IF((LR > 2) .AND. (LR <= 5)) THEN                        !
!
        IF(ABS(MD) <= 2) THEN                                       ! selection rule on MF
!
          DO J_PO = 1, 3                                            !
            DO J_LI = 1, 3                                          !
              ANG_E2(J_PO,J_LI) = ANG(2,J_PO,J_LI,INDI,INDF)        !
            END DO                                                  !
          END DO                                                    !
!
        END IF                                                      !
!
!     Magnetic dipole case : LR = 6
!
      ELSE IF(LR == 6) THEN                                         !
!
        IF(ABS(MD) <= 1) THEN                                       ! selection rule on MF
!
          IF(LF == LI) THEN                                         ! selection rule on LF
            DO JDIR = 1, 3                                          !
              ANG_M1(JDIR) = ANG(4,JDIR,1,INDI,INDF)                !
            END DO                                                  !
          END IF                                                    !
!
        END IF                                                      !
!
!     Electric octupole case : LR = 7-10
!
      ELSE IF((LR > 6) .AND. (LR <= 10)) THEN                       !
!
        IF(ABS(MD) <= 3) THEN                                       ! selection rule on MF
!
          DO J_PO = 1, 3                                            !
            DO J_LI = 1 ,3                                          !
!
              ANG_E3(J_PO,J_LI) = ANG(3,J_PO,J_LI,INDI,INDF)        !
!
            END DO                                                  !
          END DO                                                    !
!
        END IF                                                      !
!
!     Magnetic quadrupole + polar toroidal dipole case : LR = 11-12
!
      ELSE IF((LR > 10) .AND. (LR <= 12)) THEN                      !
!
        IF(ABS(MD) <= 2) THEN                                       ! selection rule on MF
!
          IF(ABS(LF-LI) == 1) THEN                                  ! selection rule on LF
!
            DO JDIR = 1, 3                                          !
              DO J_LI = 1, 3                                        !
                ANG_M2(JDIR,J_LI) = ANG(5,JDIR,J_LI,INDI,INDF)      !
              END DO                                                !
            END DO                                                  !
!
          END IF                                                    !
!
        END IF                                                      !
!
      END IF                                                        !
!
!  Calculation of the irregular matrix elements along the (POL,LIGHT) directions:
!
!
!  1) Diagonal terms e.g. Xn-Xn:
!
!
!        Electric dipole case : LR = 1-2
!
      IF(LR <= 2) THEN                                              !
!
        IF(ABS(MD) <= 1) THEN                                       ! selection rule on MF
          DO J_PO = 1, 3                                            !
            MAT_E1_E1(J_PO) = CONJG(COEF(1)) * COEF(1) * RADIAL * & !
                              CONJG(ANG_E1(J_PO)) * ANG_E1(J_PO)    !
          END DO                                                    !
        END IF                                                      !
!
!        Electric quadrupole case : LR = 3-5
!
      ELSE IF((LR > 2) .AND. (LR <= 5)) THEN                        !
!
        IF(ABS(MD) <= 2) THEN                                       ! selection rule on MF
          DO J_PO = 1, 3                                            !
            DO J_LI = 1, 3                                          !
              MAT_E2_E2(J_PO,J_LI) = CONJG(COEF(2)) * COEF(2) *   & !
                                     RADIAL *                     & !
                                     CONJG(ANG_E2(J_PO,J_LI)) *   & !
                                     ANG_E2(J_PO,J_LI)              !
            END DO                                                  !
          END DO                                                    !
        END IF                                                      !
!
!        Magnetic dipole case : LR = 6
!
      ELSE IF(LR == 6) THEN                                         !
!
        IF(ABS(MD) <= 1) THEN                                       ! selection rule on MF
!
          IF(LF == LI) THEN                                         ! selection rule on LF
            DO J_PO=1,3
              MAT_M1_M1(J_PO) = CONJG(COEF(4)) * COEF(4) *        & !
                                RADIAL *                          & !
                                CONJG(ANG_M1(J_PO)) * ANG_M1(J_PO)  !
            END DO                                                  !
          END IF                                                    !
!
        END IF                                                      !
!
!        Electric octupole case : LR = 7-10
!
      ELSE IF((LR > 6) .AND. (LR <= 10)) THEN                       !
!
        IF(ABS(MD) <= 3) THEN                                       ! selection rule on MF
          DO J_PO = 1, 3                                            !
            DO J_LI = 1, 3                                          !
              MAT_E3_E3(J_PO,J_LI) = CONJG(COEF(3)) * COEF(3) *   & !
                                     RADIAL *                     & !
                                     CONJG(ANG_E3(J_PO,J_LI)) *   & !
                                     ANG_E3(J_PO,J_LI)              !
            END DO                                                  !
          END DO                                                    !
        END IF                                                      !
!
!        Magnetic quadrupole + polar toroidal dipole case : LR = 11-12
!
      ELSE IF((LR > 10) .AND. (LR <= 12)) THEN                      !
!
        IF(ABS(MD) <= 2) THEN                                       ! selection rule on MF
!
          IF(ABS(LF-LI) == 1) THEN                                  ! selection rule on LF
            DO J_PO = 1, 3
              DO J_LI = 1, 3
                MAT_M2_M2(J_PO,J_LI) = CONJG(COEF(5)) * COEF(5) * & !
                                       RADIAL *                   & !
                                       CONJG(ANG_M2(J_PO,J_LI)) * & !
                                       ANG_M2(J_PO,J_LI)            !
              END DO                                                !
            END DO                                                  !
          END IF                                                    !
!
        END IF                                                      !
!
      END IF                                                        !
!
!  1) Non diagonal terms e.g. Xn-Ym
!
!
      RETURN                                                        !
!
      END SUBROUTINE COUMAT_IRREG
!
!=======================================================================
!
      SUBROUTINE CALC_MLIL0(JEPS,I_TEST,DIRPOL,DIRLUM,MLFLI_E1,       &
                            MLFLI_E2,MLFLI_E3,MLFLI_M1,MLFLI_M2,MLIL0)
!
!  This subroutine computes the photon-electron excitation matrix elements
!    for a given direction DIRPOL of the polarization of the light and
!    a given direction DIRLUM of the incoming photons
!
!
!
!  Input variables :
!
!                       JEPS      :  polarization index of the light
!                       I_TEST    :  switch used for testing purposes
!                       DIRPOL    :  direction of the polarization
!                       DIRLUM    :  direction of the incident light
!                       MLFLI_E1  :  electric dipole matrix elements
!                       MLFLI_E2  :  electric quadrupole matrix elements
!                       MLFLI_E3  :  electric octopole matrix elements
!                       MLFLI_M1  :  magnetic dipole matrix elements
!                       MLFLI_M2  :  magnetic quadrupole and polar toroidal
!                                             matrix elements
!
!  Output variables :
!
!                       MLIL0     :  excitation matrix elements
!
!  Types of matrix elements calculated :
!
!              * I_C1 = 1 : single channel
!              * I_E1 = 1 : electric dipole E1
!              * I_E2 = 1 : electric quadrupole E2
!              * I_E3 = 1 : electric octupole E3
!              * I_M1 = 1 : magnetic dipole M1
!              * I_M2 = 1 : magnetic quadrupole M2 + polar toroidal dipole T1
!
!
!
!  Author : D. Sébilleau
!
!                                     Last modified : 19 May 2021
!
      USE REAL_NUMBERS,         ONLY : ONE
      USE COMPLEX_NUMBERS,      ONLY : ZEROC,ONEC
!
      USE DICHROISM
      USE INIT_L
      USE ONE_CHANNEL
      USE OPTICAL_ELTS
      USE VECTOR
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)       ::  JEPS,I_TEST
!
      INTEGER                   ::  JJ,MI,LD,MD
      INTEGER                   ::  LF,ILF,MF,INDF
!
      REAL (WP), INTENT(IN)     ::  DIRPOL(3,2),DIRLUM(3)
!
      REAL (WP)                 ::   POL(3),QXE(3)
!
      COMPLEX (WP), INTENT(IN)  ::  MLFLI_E1(2,-LI_M:LI_M,-LI_M-1:LI_M+1,0:LI_M+1,3)
      COMPLEX (WP), INTENT(IN)  ::  MLFLI_E2(2,-LI_M:LI_M,-LI_M-2:LI_M+2,0:LI_M+2,3,3)
      COMPLEX (WP), INTENT(IN)  ::  MLFLI_E3(2,-LI_M:LI_M,-LI_M-3:LI_M+3,0:LI_M+3,3,3)
      COMPLEX (WP), INTENT(IN)  ::  MLFLI_M1(2,-LI_M:LI_M,-LI_M-1:LI_M+1,0:LI_M,3)
      COMPLEX (WP), INTENT(IN)  ::  MLFLI_M2(2,-LI_M:LI_M,-LI_M-2:LI_M+2,0:LI_M+1,3,3)
      COMPLEX (WP), INTENT(OUT) ::  MLIL0(2,-LI_M:LI_M,LINFMAX)
!
      COMPLEX (WP)              ::  E_DIPOLE,E_QUADRU,E_OCTUPO
      COMPLEX (WP)              ::  M_DIPOLE,M_QUADRU
!
!  Vector product LIGHT x POLARIZATION = QXE for magnetic dipole and quadrupole terms
!
      DO JJ = 1, 3                                                  !
        POL(JJ)=DIRPOL(JJ,JEPS)                                     !
      END DO                                                        !
!
      CALL PRVECT(DIRLUM,POL,QXE,ONE)                               !
!
!   Calculation of the coupling matrix MLIL0 (E1,E2,E3,M1 and M2+T1 terms)
!
      DO MI = - LI, LI                                              !
        DO LF = LF1, LF2, ISTEP_LF                                  !
          ILF = LF * LF + LF + 1                                    !
!
!  Absolute difference term for selection rule on L
!
          LD = ABS(LF - LI)                                         !
!
          DO MF = - LF, LF                                          !
            INDF = ILF + MF                                         !
            IF(I_DICHR == 0) THEN                                   !
              IF(I_TEST /=  1) THEN                                 !
!
!  Absolute difference term for selection rule on M
!
              MD = ABS(MF - MI)                                     !
!
!  Single channel case
!
                IF(I_C1 == 1) THEN                                  !
                  IF(TYPCHAN == 'E1') THEN                          !
                    I_E1 = 1                                        !
                    LD   = 1                                        !
                  ELSE IF(TYPCHAN == 'E2') THEN                     !
                    I_E2 = 1                                        !
                    LD   = 2                                        !
                  ELSE IF(TYPCHAN == 'E3') THEN                     !
                    I_E3 = 1                                        !
                    LD   = 3                                        !
                  ELSE IF(TYPCHAN == 'M1') THEN                     !
                    I_M1 = 1                                        !
                    LD   = 0                                        !
                  ELSE IF(TYPCHAN == 'M2') THEN                     !
                    I_M2 = 1                                        !
                    LD   = 1                                        !
                  ELSE IF(TYPCHAN == '  ') THEN                     !
                    LD   = 0                                        !
                  END IF                                            !
                END IF                                              !
!
!  Electric dipole term E1
!
                IF(I_E1 == 1) THEN                                  !
                  IF(LD == 1) THEN                                  ! E1 L-selection rule
                    IF(MD <= 1) THEN                                ! E1 M-selection rule
                      E_DIPOLE = MLFLI_E1(1,MI,MF,LF,1) *         & !
                                 DIRPOL(1,JEPS) +                 & !
                                 MLFLI_E1(1,MI,MF,LF,2) *         & !
                                 DIRPOL(2,JEPS) +                 & !
                                 MLFLI_E1(1,MI,MF,LF,3) *         & !
                                 DIRPOL(3,JEPS)                     !
                    ELSE                                            !
                      E_DIPOLE = ZEROC                              !
                    END IF                                          !
                  ELSE                                              !
                    E_DIPOLE = ZEROC                                !
                  END IF                                            !
                ELSE                                                !
                  E_DIPOLE = ZEROC                                  !
                END IF                                              !
!
!  Electric quadrupole term E2
!
                IF(I_E2 == 1) THEN                                  !
                  IF((LD == 2) .OR. (LD == 0)) THEN                 ! E2 L-selection rule
                    IF(MD <= 2) THEN                                ! E2 M-selection rule
                      E_QUADRU = MLFLI_E2(1,MI,MF,LF,1,1) *       & !
                                 DIRPOL(1,JEPS) * DIRLUM(1) +     & !
                                 MLFLI_E2(1,MI,MF,LF,1,2) *       & !
                                 DIRPOL(1,JEPS) * DIRLUM(2) +     & !
                                 MLFLI_E2(1,MI,MF,LF,1,3) *       & !
                                 DIRPOL(1,JEPS) * DIRLUM(3) +     & !
                                 MLFLI_E2(1,MI,MF,LF,2,1) *       & !
                                 DIRPOL(2,JEPS) * DIRLUM(1) +     & !
                                 MLFLI_E2(1,MI,MF,LF,2,2) *       & !
                                 DIRPOL(2,JEPS) * DIRLUM(2) +     & !
                                 MLFLI_E2(1,MI,MF,LF,2,3) *       & !
                                 DIRPOL(2,JEPS) * DIRLUM(3) +     & !
                                 MLFLI_E2(1,MI,MF,LF,3,1) *       & !
                                 DIRPOL(3,JEPS) * DIRLUM(1) +     & !
                                 MLFLI_E2(1,MI,MF,LF,3,2) *       & !
                                 DIRPOL(3,JEPS) * DIRLUM(2) +     & !
                                 MLFLI_E2(1,MI,MF,LF,3,3) *       & !
                                DIRPOL(3,JEPS) * DIRLUM(3)          !
                    ELSE                                            !
                      E_QUADRU = ZEROC                              !
                    END IF                                          !
                  ELSE                                              !
                    E_QUADRU = ZEROC                                !
                  END IF                                            !
                ELSE                                                !
                  E_QUADRU = ZEROC                                  !
                END IF                                              !
!
!  Electric octupole term E3
!
                IF(I_E3 == 1) THEN                                  !
                  IF((LD == 3) .OR. (LD == 1)) THEN                 ! E3 L-selection rule
                    IF(MD <= 3) THEN                                ! E3 M-selection rule
                      E_OCTUPO = MLFLI_E3(1,MI,MF,LF,1,1) *       & !
                                 DIRPOL(1,JEPS) * DIRLUM(1) +     & !
                                 MLFLI_E3(1,MI,MF,LF,1,2) *       & !
                                 DIRPOL(1,JEPS) * DIRLUM(2) +     & !
                                 MLFLI_E3(1,MI,MF,LF,1,3) *       & !
                                 DIRPOL(1,JEPS) * DIRLUM(3) +     & !
                                 MLFLI_E3(1,MI,MF,LF,2,1) *       & !
                                 DIRPOL(2,JEPS) * DIRLUM(1) +     & !
                                 MLFLI_E3(1,MI,MF,LF,2,2) *       & !
                                 DIRPOL(2,JEPS) * DIRLUM(2) +     & !
                                 MLFLI_E3(1,MI,MF,LF,2,3) *       & !
                                 DIRPOL(2,JEPS) * DIRLUM(3) +     & !
                                 MLFLI_E3(1,MI,MF,LF,3,1) *       & !
                                 DIRPOL(3,JEPS) * DIRLUM(1) +     & !
                                 MLFLI_E3(1,MI,MF,LF,3,2) *       & !
                                 DIRPOL(3,JEPS) * DIRLUM(2) +     & !
                                 MLFLI_E3(1,MI,MF,LF,3,3) *       & !
                                 DIRPOL(3,JEPS) * DIRLUM(3)         !
                    ELSE                                            !
                      E_OCTUPO = ZEROC                              !
                    END IF                                          !
                  ELSE                                              !
                    E_OCTUPO = ZEROC                                !
                  END IF                                            !
                ELSE                                                !
                  E_OCTUPO = ZEROC                                  !
                END IF                                              !
!
!  Magnetic dipole term M1
!
                IF(I_M1 == 1) THEN                                  !
                  IF(LD == 0) THEN                                  ! M1 L-selection rule
                    IF(MD <= 1) THEN                                ! M1 M-selection rule
                      M_DIPOLE = MLFLI_M1(1,MI,MF,LF,1) * QXE(1) +& !
                                 MLFLI_M1(1,MI,MF,LF,2) * QXE(2) +& !
                                 MLFLI_M1(1,MI,MF,LF,3) * QXE(3)    !
                    ELSE                                            !
                      M_DIPOLE = ZEROC                              !
                    END IF                                          !
                  ELSE                                              !
                    M_DIPOLE = ZEROC                                !
                  END IF                                            !
                ELSE                                                !
                  M_DIPOLE = ZEROC                                  !
                END IF                                              !
!
!  Magnetic quadrupole term M2 and polar toroidal dipole term T1
!
                IF(I_M1 == 1) THEN                                  !
                  IF(LD == 1) THEN                                  ! M2 L-selection rule
                    IF(MD <= 2) THEN                                ! M2 M-selection rule
                      M_QUADRU =                                  & !
                  MLFLI_M2(1,MI,MF,LF,1,1) * DIRLUM(1) * QXE(1) + & !
                  MLFLI_M2(1,MI,MF,LF,1,2) * DIRLUM(1) * QXE(2) + & !
                  MLFLI_M2(1,MI,MF,LF,1,3) * DIRLUM(1) * QXE(3) + & !
                  MLFLI_M2(1,MI,MF,LF,2,1) * DIRLUM(2) * QXE(1) + & !
                  MLFLI_M2(1,MI,MF,LF,2,2) * DIRLUM(2) * QXE(2) + & !
                  MLFLI_M2(1,MI,MF,LF,2,3) * DIRLUM(2) * QXE(3) + & !
                  MLFLI_M2(1,MI,MF,LF,3,1) * DIRLUM(3) * QXE(1) + & !
                  MLFLI_M2(1,MI,MF,LF,3,2) * DIRLUM(3) * QXE(2) + & !
                  MLFLI_M2(1,MI,MF,LF,3,3) * DIRLUM(3) * QXE(3)     !
                    ELSE                                            !
                      M_QUADRU = ZEROC                              !
                    END IF                                          !
                  ELSE                                              !
                    M_QUADRU = ZEROC                                !
                  END IF                                            !
                ELSE                                                !
                  M_QUADRU = ZEROC                                  !
                END IF                                              !
!
!  Total photon-electron coupling matrix element
!
                MLIL0(1,MI,INDF) = E_DIPOLE + E_QUADRU + E_OCTUPO +&!
                                   M_DIPOLE + M_QUADRU              !
!
              ELSE                                                  !
                MLIL0(1,MI,INDF) = ONEC                             !
              END IF                                                !
            ELSE IF(I_DICHR >= 1) THEN                              !
              IF(I_TEST /=  1) THEN                                 !
                MLIL0(1,MI,INDF) = MLFLI_E1(1,MI,MF,LF,1) *       & !
                                   DIRLUM(1) +                    & !
                                   MLFLI_E1(1,MI,MF,LF,2) *       & !
                                   DIRLUM(2) +                    & !
                                   MLFLI_E1(1,MI,MF,LF,3) *       & !
                                   DIRLUM(3)                        !
                MLIL0(2,MI,INDF) = MLFLI_E1(2,MI,MF,LF,1) *       & !
                                   DIRLUM(1) +                    & !
                                   MLFLI_E1(2,MI,MF,LF,2) *       & !
                                   DIRLUM(2) +                    & !
                                   MLFLI_E1(2,MI,MF,LF,3) *       & !
                                   DIRLUM(3)                        !
              ELSE                                                  !
                MLIL0(1,MI,INDF) = ONEC                             !
              END IF                                                !
            END IF                                                  !
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
!
      END SUBROUTINE CALC_MLIL0
!
END MODULE ELECTRON_PHOTON
