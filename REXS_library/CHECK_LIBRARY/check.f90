MODULE CHECK
!
!  This module provides tools for checking different properties
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE CHECK_VIB(NAT2)
!
!  This subroutines checks the geometrical environment of each atom
!     to identify those which can move "freely" in one direction, in
!     order to see whether the mean square displacement in this
!     direction is of bulk type or surface type
!
!  An atom is considered to move freely in one direction if no other
!     atom is present in the tetragonal cell of height ALENGTH * A
!     and base edge 2 * A, whose base is centered on the atom considered
!
!  Only prototypical atoms are considered as all equivalent atoms are
!     in the same geometrical environment
!
!  Surface-like atoms are then identified as having I_FREE = 1
!
!
!
!   Author :  D. Sébilleau
!
!                                     Last modified : 19 May 2021
!
!
      USE REAL_NUMBERS,           ONLY : ONE,FOUR,SMALL
!
      USE CLUSTER_COORD
      USE OUTUNITS
      USE VIBR_TYPE
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  NAT2
!
      INTEGER               ::  JTYP,JAT,J,N_SUR,JAT0
      INTEGER               ::  I_ACC
      INTEGER               ::  NSUR(NATP_M)
!
      REAL (WP)             ::  ALENGTH
      REAL (WP)             ::  XA,YA,ZA
      REAL (WP)             ::  X,Y,Z
      REAL (WP)             ::  D_LAT,D_VER
!
      ALENGTH = FOUR                                                !
!
!.................... Checking the z direction ....................
!
      WRITE(IUO1,11)                                                !
      N_SUR = 0                                                     !
!
!  Loop on the prototypical atoms
!
      DO JTYP = 1, N_PROT                                           !
!
        I_FREE(JTYP) = 0                                            !
        JAT0         = NCORR(1,JTYP)                                !
        XA           = SYM_AT(1,JAT0)                               !
        YA           = SYM_AT(2,JAT0)                               !
        ZA           = SYM_AT(3,JAT0)                               !
!
!  Loop on the surrounding atoms
!
        I_ACC = 0                                                   !
!
        DO JAT = 1, NAT2                                            !
!
          IF(JAT == JAT0) GO TO 10                                  !
!
          X = SYM_AT(1,JAT)                                         !
          Y = SYM_AT(2,JAT)                                         !
          Z = SYM_AT(3,JAT)                                         !
!
!  Considering only atoms with Z > ZA
!
          IF(Z < (ZA + SMALL)) GO TO 10                             !
!
!  Lateral and vertical distances between the two atoms
!
          D_LAT = (X - XA) * (X - XA) + (Y - YA) * (Y - YA)          !
          D_VER = (Z - ZA) * (Z - ZA)                                !
!
          IF(D_VER < (ALENGTH + SMALL)) THEN                         !
            IF(D_LAT < (ONE + SMALL)) THEN                           !
              I_ACC = I_ACC + 1                                      !
            END IF                                                   !
          END IF                                                     !
!
          IF(I_ACC >= 1) GO TO 10                                    !
!
  10      CONTINUE                                                   !
!
        END DO                                                       !
!
        IF(I_ACC == 0) THEN                                          !
          I_FREE(JTYP) = 1                                           !
          N_SUR        = N_SUR + 1                                   !
          NSUR(N_SUR)  = JTYP                                        !
         END IF                                                      !
!
      END DO                                                         !
!
      WRITE(IUO1,12) (NSUR(J),J=1,N_SUR)                             !
!
!  Formats:
!
  11  FORMAT(//,18X,'SURFACE-LIKE ATOMS FOR MSD CALCULATIONS: ',/)
  12  FORMAT(20X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5,2X,I5)
!
      END SUBROUTINE CHECK_VIB
!
!=======================================================================
!
      SUBROUTINE CHECK_EMPTY_SPHERES(N_L)

!  This subroutine checks for the presence of empty spheres.
!
!    If the top layers are only composed of empty spheres,
!    N_ESL is set to the number of top layers containing
!    exclusively empty spheres (layers of empty spheres
!    that might be inside the cluster are not counted)
!
!    For empty spheres intercalated within the cluster as well
!    as for other empty spheres, this subroutine associates
!    to each atom a switch I_ES that indicates whether it is
!    an empty sphere or a real atom
!
!
!  Input parameters:
!
!       * N_L      : number of layers in the cluster
!
!
!  Intermediate parameters:
!
!       * VALZ     : z position of the different layers
!       * SYM_AT   : coordinates of the atoms
!                    ATOM(J,JAT): coordinate J (J=1,2,3) of atom number JAT
!       * NZAT     : atomic number of atoms
!       * NATP_M   : number of prototypical atoms
!       * NATYP    : number of equivalent atoms
!
!
!  Output parameters:
!
!       * N_ESL    : number of top layers composed only of empty spheres
!       * I_ES     : switch to indicate whether the atom is
!                       -  a real atom                  (=0)
!                       -  a top layers empty sphere    (=1)
!                       -  an intercalated empty sphere (=2)
!
!                    These parameters are stored in the common block /EMPTY_S/
!
!
!   Author :  D. Sébilleau
!
!                                          Last modified :  19 May 2021
!
!
      USE ATOMS
      USE EMPTY_SPHERES
      USE CLUSTER_COORD
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)   ::  N_L
!
      INTEGER               ::  NAL,NES
      INTEGER               ::  JL,JTYP,JAT,JEQ
      INTEGER               ::  NEQ,NEQJ
      INTEGER               ::  NRN,KTYP,KAT,NEQK,KEQ
!
      REAL (WP)             ::  SSMALL
      REAL (WP)             ::  Z_LAY,Z,ZJ,ZK
!
      REAL (WP)             ::  ABS
!
      DATA  SSMALL     / 0.001E0_WP /
!
!  1) Checking for the presence of overlayers of empty spheres
!
      N_ESL = 0                                                     !
!
      NAL   = 0                                                     !
      NES   = 0                                                     !
!
      DO JL = 1 ,N_L                                                !
!
        Z_LAY = VALZ(JL)                                            !
        JAT   = 0                                                   !
!
        DO JTYP = 1, NATP_M                                         !
!
          NEQ = NATYP(JTYP)                                         !
!
          DO JEQ = 1, NEQ                                           !
!
            JAT = JAT + 1                                           !
!
!  Computing the number of atoms in layer JL (NAL)
!    and the number of empty spheres in this layer (NES)
!
            Z = SYM_AT(3,JAT)                                       !
            IF(ABS(Z - Z_LAY) < SSMALL) THEN                        !
              NAL = NAL + 1                                         !
              IF(NZAT(JTYP) == 0) THEN                              !
                NES = NES + 1                                       !
              END IF                                                !
            END IF                                                  !
!
          END DO                                                    !
!
        END DO                                                      !
!
!  If layer of empty spheres found, increment the number N_ESL
!       of top layers of empty spheres
!
        IF(NAL == NES) THEN                                         !
           N_ESL = N_ESL + 1                                        !
        ELSE                                                        !
           GO TO 10                                                 !
        END IF                                                      !
!
      END DO                                                        !
!
  10  CONTINUE                                                      !
!
!  2) Identifying atoms either as real atoms or as empty spheres
!       and storage of the outcome in the I_ES array
!
      JAT = 0                                                       !
      DO JTYP = 1, NATP_M                                           !
!
        NEQJ = NATYP(JTYP)                                          !
!
        DO JEQ = 1, NEQJ                                            !
!
          JAT = JAT + 1                                             !
          ZJ  = SYM_AT(3,JAT)                                       !
!
          IF(NZAT(JTYP) > 0) THEN                                   !
!
!.............  Real atom  .............
!
            I_ES(JAT) = 0                                           !
!
          ELSE                                                      !
!
!.............  Empty sphere  .............
!
!    Checking if in top layer or intercalated
!
!     (intercalated: either real atom with higher z
!      or real atom with same z)
!     NRN: number of real atoms with z >= current empty sphere
!
!
            NRN = 0                                                 !
            KAT = 0                                                 !
            DO KTYP = 1, NATP_M                                     !
!
              NEQK = NATYP(KTYP)                                    !
!
              DO KEQ = 1, NEQK                                      !
!
                KAT = KAT + 1                                       !
                IF(KAT == JAT) GO TO 20                             !
                ZK = SYM_AT(3,KAT)                                  !
                IF(NZAT(KTYP) > 0) THEN                             !
                  IF( (ABS(ZK - ZJ) < SSMALL) .OR.                & !
                      (ZK > ZJ + SSMALL) ) THEN                     !
                    NRN = NRN + 1                                   !
                  END IF                                            !
                END IF                                              !
  20            CONTINUE                                            !
!
              END DO                                                !
!
            END DO                                                  !
!
            IF(NRN == 0) THEN                                       !
!
!............. Top layers empty sphere  .............
!
              I_ES(JAT) = 1                                         !
!
            ELSE                                                    !
!
!............. Intercalated empty sphere  .............
!
              I_ES(JAT) = 2                                         !
!
            END IF                                                  !
!
          END IF                                                    !
!
        END DO                                                      !
!
      END DO                                                        !
!
      END SUBROUTINE CHECK_EMPTY_SPHERES
!
END MODULE CHECK
