!
!=======================================================================
!
MODULE AUGER_MULTIPLET
!
!  This module computes the Auger multiplets
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE,       ONLY : NL_M
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE AUGER_MULT
!
!  This subroutine computes all the possible multiplets that are
!     contained in a given Auger transition line. It assumes that
!     the atom has closed shells only.
!
!
!   Author :  D. Sébilleau
!
!                                        Last modified : 19 May 2021
!
!
      USE INIT_A
      USE OUTUNITS
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)  ::  MULTIPLET(112)
      CHARACTER (LEN =  1)  ::  SC(0:1),LC(0:6),JC(0:7)
!
      INTEGER               ::  N_MULT,NS,L,J
!
      DATA SC /'1','3'/
      DATA LC /'S','P','D','F','G','H','I'/
      DATA JC /'0','1','2','3','4','5','6','7'/
!
      WRITE(IUO1,10)                                                !
!
      N_MULT = 0                                                    !
      DO NS = 0, 1                                                  !
        DO L = ABS(LI_A - LI_I), LI_A + LI_I                        !
          DO J = ABS(L - NS), L + NS                                !
            N_MULT            = N_MULT + 1                          !
            MULTIPLET(N_MULT) = SC(NS)//LC(L)//JC(J)                !
            WRITE(IUO1,20) MULTIPLET(N_MULT)                        !
          END DO                                                    !
        END DO                                                      !
      END DO                                                        !
!
!  Formats:
!
  10  FORMAT(///,26X,'THE POSSIBLE MULTIPLETS ARE :',/,'  ')
  20  FORMAT(58X,A3)
!
      END SUBROUTINE AUGER_MULT
!
END MODULE AUGER_MULTIPLET
