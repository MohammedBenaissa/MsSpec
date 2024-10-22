!
!=======================================================================
!
MODULE DATA_TREATMENT
!
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE TREAT_DATA(N_INFILES,NP)
!
!  This routine calles the treatment routine that transforms the
!  calculated data to prepare it for the code that will process
!  it before plotting
!
!
!  Input variables :
!
!                       N_INFILES :  number of input data files
!
!
!
!   Author :  D. Sébilleau
!
!
!                                           Last modified : 19 May 2021
!
!
      USE EXP_TYPE
      USE CURRENT_SPIN
      USE TESTS
!
      USE TREATMENT_PED
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)        ::  N_INFILES,NP
!
      INTEGER                    ::  JFF
!
      IF(ISOM /= 0) THEN                                            !
!
        JFF = 1                                                     !
        IF(ISPIN == 0) THEN                                         !
          IF(SPECTRO /= 'XAS') THEN                                 !
            CALL TREAT_PED(ISOM,N_INFILES,JFF,NP)                   !
          ELSE                                                      !
!            CALL TREAT_XAS(ISOM,N_INFILES,NP)                     !
          END IF                                                    !
        ELSE IF(ISPIN == 1) THEN                                    !
!          IF((SPECTRO == 'PED') .OR. (SPECTRO == 'AED')) THEN     !
!            CALL TREAT_PED_SP(ISOM,N_INFILES,JFF,NP)              !
!          ELSE IF(SPECTRO == 'XAS') THEN                          !
!            CALL TREAT_XAS_SP(ISOM,N_INFILES,NP)                  !
!          END IF                                                  !
          CONTINUE                                                  !
        END IF                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE TREAT_DATA
!
END MODULE DATA_TREATMENT
