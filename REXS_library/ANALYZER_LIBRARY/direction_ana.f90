!
!=======================================================================
!
MODULE DIRECTION_ANA
!
!  This module computes the direction of the analyzer
!
      USE ACCURACY_REAL
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE INIT_INC_BEAM(IS,TH_IN,PH_IN,TH_0,PH_0,IMOD,I_EXT, &
                               I_PH1,RTH,RPH,DIR_BEAM,BEAM)
!
!  This subroutine initializes the position of the incoming beam
!
!  Three cases are considered in this subroutine:
!
!      1) first initialization (from position of detector along 0z)
!      2) new initialization for each value of the fixed angle (N_FIXED > 1)
!      3) new initialization for each new set of scans (NSET > 1)
!
!  Input variables:
!
!     IS        : switch to select the type of initialization
!                    = 1 --> case 1) (see above)
!                    = 2 --> case 2) (see above)
!                    = 3 --> case 3) (see above)
!     TH_IN     : theta angle of incoming beam when detector is along z: case 1)
!                 actual theta angle of the detector: case 2)
!     PH_IN     : phi angle of incoming beam when detector is along z: case 1)
!                 actual phi angle of the detector: case 2)
!     TH_0      : starting value of theta for the detector (case 3): for the actual set)
!     PH_0      : starting value of phi for the detector (case 3): for the actual set)
!     IMOD      : rotation mode
!                    = 0 --> detector is rotated
!                    = 1 --> sample is rotated
!     I_EXT     : switch for external reading of the detector positions
!                    = 0 --> detector angles generated internally
!                    = 1 --> detector angles read externally
!                    =-1 --> detector angles read externally + averaging
!     I_PH1     : switch to indicate whether phi is the first angle looped
!                   on or not (the angle in the outer angle loop is called
!                   the 'fixed angle')
!                    = 0 --> theta angle is the 'fixed angle'
!                    = 1 --> phi angle is the 'fixed angle'
!
!
!  Output variable:
!
!     RTH       : starting value of theta (radians) for the detector
!                      (kept for further use)
!     RPH       : starting value of phi (radians) for the detector
!                      (kept for further use)
!     DIR_BEAM  : position of the incoming beam when the detector is along 0z
!     BEAM      : position of the incoming beam when the detector is along (TH_0,PH_0)
!
!
!  Notes: the denomination 'fixed' angle, as opposed to 'scanned angle'
!           refers to that of theta and phi that has the less values
!
!         angles starting by 'R' are in radians
!
!
!  Author : D. Sébilleau
!
!                                           Last modified:  7 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)    ::  IS,IMOD,I_EXT,I_PH1
!
      REAL (WP), INTENT(IN)  ::  TH_IN,PH_IN,TH_0,PH_0
      REAL (WP), INTENT(OUT) ::  RTH,RPH,DIR_BEAM(3),BEAM(3)
!
      REAL (WP)              ::  PIS180
      REAL (WP)              ::  R_L(9)
      REAL (WP)              ::  RTH_IN,RPH_IN
      REAL (WP)              ::  X_BEAM_Z,Y_BEAM_Z,Z_BEAM_Z
      REAL (WP)              ::  RTH_0,RPH_0
!
      REAL (WP)              ::  SIN,COS
!
      DATA PIS180     / 0.01745329251994329576923690768488612713E0_WP /
!
      IF(IS == 1) THEN                                              !
!
!  Case 1) : first initialization (detector along 0z ---> detector at (TH_0,PH_0)
!
       RTH_IN = TH_IN * PIS180                                      !
       RPH_IN = PH_IN * PIS180                                      !
!
       X_BEAM_Z = SIN(RTH_IN) * COS(RPH_IN)                         !
       Y_BEAM_Z = SIN(RTH_IN) * SIN(RPH_IN)                         !
       Z_BEAM_Z = COS(RTH_IN)                                       !
!
       IF(IMOD == 0) THEN                                           !
!
!  The analyzer is rotated:
!
!  Position of the incoming beam when the analyzer is along the z axis :
!                  (X_LUM_Z,Y_LUM_Z,Z_LUM_Z)
!
         DIR_BEAM(1) = X_BEAM_Z                                     !
         DIR_BEAM(2) = Y_BEAM_Z                                     !
         DIR_BEAM(3) = Z_BEAM_Z                                     !
!
       ELSE                                                         !
!
!  The sample is rotated ---> light and analyzer rotated
!
         IF(I_EXT == 0) THEN                                        !
!
           RTH_0 = TH_0 * PIS180                                    !
           RPH_0 = PH_0 * PIS180                                    !
           RTH   = RTH_0                                            !
           RPH   = RPH_0                                            !
!
!  R_L is the rotation matrix from 0z to (TH_0,PH_0) expressed as
!    a function of the Euler angles ALPHA=PH_0, BETA=TH_0, GAMMA=-PH_0
!    It is stored as (1 2 3)
!                    (4 5 6)
!                    (7 8 9)
!
           R_L(1) =   COS(RTH_0) * COS(RPH_0) * COS(RPH_0) +      & !
                                   SIN(RPH_0) * SIN(RPH_0)          !
           R_L(2) =   COS(RTH_0) * SIN(RPH_0) * COS(RPH_0) -      & !
                                   SIN(RPH_0) * COS(RPH_0)          !
           R_L(3) =   SIN(RTH_0) * COS(RPH_0)                       !
           R_L(4) =   COS(RTH_0) * SIN(RPH_0) * COS(RPH_0) -      & !
                                   SIN(RPH_0) * COS(RPH_0)          !
           R_L(5) =   COS(RTH_0) * SIN(RPH_0) * SIN(RPH_0) +      & !
                                   COS(RPH_0) * COS(RPH_0)          !
           R_L(6) =   SIN(RTH_0) * SIN(RPH_0)                       !
           R_L(7) = - SIN(RTH_0) * COS(RPH_0)                       !
           R_L(8) = - SIN(RTH_0) * SIN(RPH_0)                       !
           R_L(9) =   COS(RTH_0)                                    !
!
!  Position of the incoming beam when the detector is along (TH_0,PH_0) : BEAM(3)
!
           BEAM(1) = X_BEAM_Z * R_L(1) +                          & !
                     Y_BEAM_Z * R_L(2) +                          & !
                     Z_BEAM_Z * R_L(3)                              !
           BEAM(2) = X_BEAM_Z * R_L(4) +                          & !
                     Y_BEAM_Z * R_L(5) +                          & !
                     Z_BEAM_Z * R_L(6)                              !
           BEAM(3) = X_BEAM_Z * R_L(7) +                          & !
                     Y_BEAM_Z * R_L(8) +                          & !
                     Z_BEAM_Z * R_L(9)                              !
!
         END IF                                                     !
!
       END IF                                                       !
!
      ELSE IF(IS == 2) THEN                                         !
!
!  Case 2) second initialization (initial position BEAM
!            of the incoming beam is recalculated
!            for each initial position (TH,PH) of the analyzer
!
        IF(I_PH1 == 1) THEN                                         !
          RTH = TH_0  * PIS180                                      !
          RPH = PH_IN * PIS180                                      !
        ELSE                                                        !
          RTH = TH_IN * PIS180                                      !
          RPH = PH_0  * PIS180                                      !
        END IF                                                      !
!
        R_L(1) =   COS(RTH) * COS(RPH)                              !
        R_L(2) = - SIN(RPH)                                         !
        R_L(3) =   SIN(RTH) * COS(RPH)                              !
        R_L(4) =   COS(RTH) * SIN(RPH)                              !
        R_L(5) =   COS(RPH)                                         !
        R_L(6) =   SIN(RTH) * SIN(RPH)                              !
        R_L(7) = - SIN(RTH)                                         !
        R_L(8) =   ZERO                                             !
        R_L(9) =   COS(RTH)                                         !
!
        BEAM(1) = X_BEAM_Z * R_L(1) +                             & !
                  Y_BEAM_Z * R_L(2) +                             & !
                  Z_BEAM_Z * R_L(3)                                 !
        BEAM(2) = X_BEAM_Z * R_L(4) +                             & !
                  Y_BEAM_Z * R_L(5) +                             & !
                  Z_BEAM_Z * R_L(6)                                 !
        BEAM(3) = X_BEAM_Z * R_L(7) +                             & !
                  Y_BEAM_Z * R_L(8) +                             & !
                  Z_BEAM_Z * R_L(9)                                 !
!
      ELSE IF(IS == 3) THEN
!
!  Case 3) third initialization (initial position BEAM recalculated
!            each initial position of the analyzer)
!
        RTH = TH_0 * PIS180                                         !
        RPH = PH_0 * PIS180                                         !
!
        IF(IMOD == 1) THEN                                          !
!
          R_L(1) =   COS(RTH) * COS(RPH)                            !
          R_L(2) = - SIN(RPH)                                       !
          R_L(3) =   SIN(RTH) * COS(RPH)                            !
          R_L(4) =   COS(RTH) * SIN(RPH)                            !
          R_L(5) =   COS(RPH)                                       !
          R_L(6) =   SIN(RTH) * SIN(RPH)                            !
          R_L(7) = - SIN(RTH)                                       !
          R_L(8) =   ZERO                                           !
          R_L(9) =   COS(RTH)                                       !
!
          BEAM(1) = X_BEAM_Z * R_L(1) +                           & !
                    Y_BEAM_Z * R_L(2) +                           & !
                    Z_BEAM_Z * R_L(3)                               !
          BEAM(2) = X_BEAM_Z * R_L(4) +                           & !
                    Y_BEAM_Z * R_L(5) +                           & !
                    Z_BEAM_Z * R_L(6)                               !
          BEAM(3) = X_BEAM_Z * R_L(7) +                           & !
                    Y_BEAM_Z * R_L(8) +                           & !
                    Z_BEAM_Z * R_L(9)                               !
!
        END IF                                                      !
!
      END IF                                                        !
!
      END SUBROUTINE INIT_INC_BEAM
!
!=======================================================================
!
      SUBROUTINE DIR_ANA(V_INT,E_K,RTH_EXT,RPH_EXT)
!
!  This subroutine calculates the direction(s) of the analyzer with
!             or without angular averaging.
!
!
!  Input parameters:
!
!       * V_INT    : inner potential (difference between vacuum
!                      level and muffin-tin zero (in eV))
!       * E_K      : kinetic energy of the electron (in eV)
!       * RTH_EXT  : internal theta value
!       * RPH_EXT  : internal phi value
!
!
!
!  Output parameters:
!
!       * DIRANA    : internal directions of the analyzer
!       * ANADIR    : external directions of the analyzer
!       * RTH_IN    : internal theta values for angular averaging over
!                       the entrance slit of the analyzer
!       * RPH_IN    : internal phi values for angular averaging over
!                       the entrance slit of the analyzer
!
!                  stored in module ANA_DIR
!
!   Author :  D. Sébilleau
!
!                                          Last modified : 26 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : HALF,SMALL
      USE COMPLEX_NUMBERS,   ONLY : IC
      USE PI_ETC,            ONLY : PI
!
      USE CURRENT_AVER
      USE CALC_TYPE
      USE EMPTY_SPHERES
      USE OUTUNITS
      USE TESTS
!
      USE ANA_DIR
      USE ARC_SIN
!
      IMPLICIT NONE
!
      INTEGER                  ::  I_PRINT
      INTEGER                  ::  I_PHI,I_THETA
      INTEGER                  ::  N,J,K1,K2
!
      REAL (WP), INTENT(IN)    ::  E_K
      REAL (WP), INTENT(INOUT) ::  V_INT
      REAL (WP), INTENT(INOUT) ::  RTH_EXT,RPH_EXT
!
      REAL (WP)                ::  PIS2,PIS180
      REAL (WP)                ::  S,RN,RK1,RK2,C,D
      REAL (WP)                ::  DTH_INT,DTH_EXT
      REAL (WP)                ::  RTH_INT,RPH_INT
      REAL (WP)                ::  THETA_R
      REAL (WP)                ::  PHINT_R,THINT_R
!
      REAL (WP)                ::  SIN,COS,FLOAT,SQRT,ACOS
!
      COMPLEX (WP)             ::  COEF
!
      DATA PIS180     / 0.01745329251994329576923690768488612713E0_WP /
!
      PIS2    = HALF * PI                                           !
!
      I_PRINT = IPRINT                                              !
      I_PHI   = IPHI                                                !
      I_THETA = ITHETA                                              !
!
!  Computing and storing the internal/external direction of the analyzer
!
      ANADIR(1,1) = SIN(RTH_EXT) * COS(RPH_EXT)                     !
      ANADIR(2,1) = SIN(RTH_EXT) * SIN(RPH_EXT)                     !
      ANADIR(3,1) = COS(RTH_EXT)                                    !
!
      IF((ABS(I_EXT) <= 1) .AND. (I_TEST /= 2)) THEN                !
        IF(N_ESL == 0) THEN                                         !
          CALL REFRAC(V_INT,E_K,RTH_EXT,RTH_INT)                    !
        ELSE                                                        !
          RTH_INT = RTH_EXT                                         !
        END IF                                                      !
        RPH_INT = RPH_EXT                                           !
      ELSE                                                          !
        RTH_INT = RTH_EXT                                           !
        RPH_INT = RPH_EXT                                           !
      END IF                                                        !
      IF((I_PRINT > 0) .AND. (I_EXT /= 2)) THEN                     !
        DTH_EXT = RTH_EXT / PIS180                                  !
        DTH_INT = RTH_INT / PIS180                                  !
        IF(I_TEST /= 2) WRITE(IUO1,20) DTH_EXT,DTH_INT              !
      END IF                                                        !
!
      DIRANA(1,1) = SIN(RTH_INT) * COS(RPH_INT)                     !
      DIRANA(2,1) = SIN(RTH_INT) * SIN(RPH_INT)                     !
      DIRANA(3,1) = COS(RTH_INT)                                    !
!
      RTH_IN(1) = RTH_INT                                           !
      RPH_IN(1) = RPH_EXT                                           !
!
      IF(I_THETA == 1) THEN                                         !
        IF(RPH_EXT > PIS2) THEN                                     !
          RTH_EXT = - RTH_EXT                                       !
          RPH_EXT =   RPH_EXT - PI                                  !
        ELSE IF(RPH_EXT < - PIS2) THEN                              !
          RTH_EXT = - RTH_EXT                                       !
          RPH_EXT =   RPH_EXT + PI                                  !
        END IF                                                      !
      END IF                                                        !
!
      IF(I_AVER >= 1) THEN                                          !
!
!  Computing the other internal/external directions for angular averaging
!          (to account for the acceptance of the analyzer)
!
        N  = 2**(I_AVER - 1)                                        !
        S  = SIN(ACCEPT * PI / 180.0E0_WP)                          !
        RN = FLOAT(N)                                               !
        J  = 1                                                      !
        DO K1 = - N, N                                              !
          RK1 = FLOAT(K1)                                           !
          DO K2 = - N, N                                            !
            RK2 = FLOAT(K2)                                         !
            D   = SQRT(RK1 * RK1 + RK2 * RK2)                       !
            IF((D - RN) > SMALL) GO TO 10                           !
            IF((K1 == 0) .AND. (K2 == 0)) GO TO 10                  !
            C = SQRT( RN * RN - (RK1 * RK1 + RK2 * RK2) * S*  S )   !
            J = J + 1                                               !
!
            ANADIR(1,J) = ( RK1 * S * COS(RTH_EXT) * COS(RPH_EXT) & !
                          - RK2 * S * SIN(RPH_EXT) +              & !
                            C * ANADIR(1,1)                       & !
                          ) / RN                                    !
            ANADIR(2,J) = ( RK1 * S * COS(RTH_EXT) * SIN(RPH_EXT) & !
                          + RK2 * S * COS(RPH_EXT) +              & !
                            C * ANADIR(2,1)                       & !
                          ) / RN                                    !
            ANADIR(3,J) = ( - RK1 * S * SIN(RTH_EXT) +            & !
                            C * ANADIR(3,1)                       & !
                          ) / RN                                    !
!
            THETA_R = ACOS(ANADIR(3,J))                             !
            COEF    = ANADIR(1,J) + IC * ANADIR(2,J)                !
            CALL ARCSIN(COEF,ANADIR(3,J),PHINT_R)                   !
!
            IF((ABS(I_EXT) <= 1) .AND. (I_TEST /= 2)) THEN          !
              IF(N_ESL == 0) THEN                                   !
                CALL REFRAC(V_INT,E_K,THETA_R,THINT_R)              !
              ELSE                                                  !
                THINT_R = THETA_R                                   !
              END IF                                                !
            ELSE                                                    !
              THINT_R = THETA_R                                     !
            END IF                                                  !
!
            DIRANA(1,J) = SIN(THINT_R) * COS(PHINT_R)               !
            DIRANA(2,J) = SIN(THINT_R) * SIN(PHINT_R)               !
            DIRANA(3,J) = COS(THINT_R)                              !
!
            RTH_IN(J) = THINT_R                                     !
            RPH_IN(J) = PHINT_R                                     !
!
  10        CONTINUE                                                !
          END DO                                                    !
        END DO                                                      !
!
      END IF                                                        !
!
!  Formats:
!
  20  FORMAT(/,10X,'ELECTRON EXTERNAL THETA  =',F7.2,5X,       &
                   'INTERNAL THETA =', F7.2)
!
      END SUBROUTINE DIR_ANA
!
!=======================================================================
!
      SUBROUTINE REFRAC(VINT,EKIN,RTHETA,RTHINT)
!
!  This routine calculates the refraction of a plane wave beam induced
!     by the surface potential barrier VINT. EKIN is the kinetic energy
!     outside the crystal.
!
!                                          Last modified :  7 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE
!
      IMPLICIT NONE
!
      REAL (WP), INTENT(IN)    ::  EKIN,RTHETA
      REAL (WP), INTENT(OUT)   ::  RTHINT
      REAL (WP), INTENT(INOUT) ::  VINT
!
      REAL (WP)                ::  PIS180,SSMALL
      REAL (WP)                ::  U,DTHETA,REFRA
!
      REAL (WP)                ::  ABS,SIN,ASIN,SQRT
!
!
      DATA PIS180     / 0.01745329251994329576923690768488612713E0_WP /
      DATA SSMALL     / 0.001E0_WP /
!
      IF(VINT < ZERO) VINT = ABS(VINT)                              !
!
      IF(ABS(VINT) < SSMALL) THEN                                   !
        RTHINT = RTHETA                                             !
      ELSE                                                          !
        U      = VINT / (EKIN + VINT)                               !
        DTHETA = RTHETA / PIS180                                    !
        REFRA  = SIN(RTHETA) * SIN(RTHETA) * (ONE - U)              !
        RTHINT = ASIN(SQRT(REFRA))                                  !
        IF(DTHETA < ZERO) THEN                                      !
          RTHINT = - RTHINT                                         !
        END IF                                                      !
      END IF                                                        !
!
      END SUBROUTINE REFRAC
!
!======================================================================
!
      SUBROUTINE ROTATE_INC_BEAM(RTHETA,RPHI,RTH,RPH,ITHETA,IPHI,  &
                                 J_SCAN,BEAM,DIR_BEAM)
!
!  This subroutine computes the position of the incoming beam
!    when the detector is at position (TH,PH) and it is the sample
!    that is rotated.
!
!  Input variables:
!
!     RTHETA    : actual theta angle of the detector
!     RPHI      : actual phi angle of the detector
!     RTH       : initial theta angle of the detector for the current fixed angle
!     RPH       : initial phi angle of the detector for the current fixed angle
!     ITHETA    : switch for polar scan
!     IPHI      : switch for azimutal scan
!     J_SCAN    : position of direction in the scan
!     BEAM      : initial position of the incoming beam when the
!                   detector is at its initial position (TH_0,PH_0)
!
!
!  Output variable:
!
!     DIR_BEAM  : position of the incoming beam corresponding
!                  to the detector at (TH,PH)
!
!
!  Note: angles starting by 'R' are in radians
!
!  Author : D. Sébilleau
!
!                                           Last modified: 19 May 2021
!
!
      USE REAL_NUMBERS,      ONLY : ZERO,ONE
      USE VECTOR
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)      ::  ITHETA,IPHI,J_SCAN
!
      REAL (WP), INTENT(IN)    ::  RTHETA,RPHI,RTH,RPH,BEAM(3)
      REAL (WP), INTENT(OUT)   ::  DIR_BEAM(3)
!
      REAL (WP)                ::  AXIS(3),EPS(3)
      REAL (WP)                ::  RANGLE,CVECT,PRS
!
      REAL (WP)                ::  SIN,COS
!
!  AXIS is the direction of the theta rotation axis and
!     EPS is defined so that (AXIS,DIR_BEAM,EPS) is a
!     direct orthonormal basis. The transform of a vector R
!     by a rotation of OMEGA about AXIS is then given by
!
!     R' = R COS(OMEGA) + (AXIS.R)(1-COS(OMEGA)) AXIS + (AXIS^R) SIN(OMEGA)
!
!  Note that the initial position of the detector is (RTH,RPH)
!     which coincides with (RTH0,RPH0) only for the first fixed angle
!
      CVECT = ONE                                                   !
!
      IF(ITHETA == 1) THEN                                          !
        AXIS(1) = - SIN(RPH)                                        !
        AXIS(2) =   COS(RPH)                                        !
        AXIS(3) =   ZERO                                            !
        RANGLE  = RTHETA - RTH                                      !
      ELSE IF(IPHI == 1) THEN                                       !
        AXIS(1) = ZERO                                              !
        AXIS(2) = ZERO                                              !
        AXIS(3) = ONE                                               !
        RANGLE  = RPHI - RPH                                        !
      END IF
!
      CALL PRVECT(AXIS,BEAM,EPS,CVECT)                              !
!
      PRS = PRSCAL(AXIS,BEAM)                                       !
!
      IF(J_SCAN == 1) THEN                                          !
        DIR_BEAM(1) = BEAM(1)                                       !
        DIR_BEAM(2) = BEAM(2)                                       !
        DIR_BEAM(3) = BEAM(3)                                       !
      ELSE                                                          !
        DIR_BEAM(1) = BEAM(1) * COS(RANGLE) +                     & !
                      PRS * (ONE - COS(RANGLE)) * AXIS(1) +       & !
                      SIN(RANGLE) * EPS(1)                          !
        DIR_BEAM(2) = BEAM(2) * COS(RANGLE) +                     & !
                      PRS * (ONE - COS(RANGLE)) * AXIS(2) +       & !
                      SIN(RANGLE) * EPS(2)                          !
        DIR_BEAM(3) = BEAM(3) * COS(RANGLE) +                     & !
                      PRS * (ONE - COS(RANGLE)) * AXIS(3) +       & !
                      SIN(RANGLE) * EPS(3)                          !
      END IF                                                        !
!
      END SUBROUTINE ROTATE_INC_BEAM
!
END MODULE DIRECTION_ANA
