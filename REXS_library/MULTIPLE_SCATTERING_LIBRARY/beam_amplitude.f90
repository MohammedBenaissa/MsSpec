!
!======================================================================
!
MODULE BEAM_AMPLITUDE
!
!  This mudule contains subroutine to compute the MS beam amplitude B_L
!
      USE ACCURACY_REAL
      USE DIMENSION_CODE
!
CONTAINS
!
!======================================================================
!
      SUBROUTINE BL_AMP(TAU_TYPE,ALGO,JE,VALS,NA_TYP,NTYPEM,JBL,JATBL, &
                        JDIR,ZSURF,GAMMA,TAU,AL,BL)
!
!  This subroutine computes the multiple scattering amplitude
!    for a beam of electron outgoing from or incoming on a given atom
!
!  Input variables:
!
!     TAU_TYPE : type of scattering path operator
!                    'ABS' --> absorption case
!                    'INC' --> incoming beam
!                    'OUT' --> outgoing beam
!     ALGO     : type of algorithm used
!                    'SE'  --> series expansion
!                    'MI'  --> matrix inversion
!                    'CE'  --> correlation expansion
!                    'RE'  --> renormalized expansion
!     JE       : energy point for which BL is computed
!     VALS     : z position of surface
!     NA_TYP   : prototypical atoms indices array
!     NTYPEM   : index of prototypical atom corresponding to the absorber
!     JBL      : index of prototypical atom for which BL is computed
!     JATBL    : index number of the atom for which BL is computed
!     ZSURF    : z position of the atom for which BL is computed
!     JDIR     : index of the incoming/outgoing direction
!                    is equal to 1 for incoming beam
!                    varies from 1 to 49 for outgoing beam (averaging in a cone)
!     GAMMA    : imaginary part ok VK when the potential is real
!     TAU      : scattering path operator
!
!
!     Note: Tau^{j 0} and Tau^{0 j} are stored in the same way
!                          to avoid having 2 TAU arrays
!
!  Output variables:
!
!     AL       : single atom amplitude
!     BL       : multiple scattering amplitude corresponding to atom JATBL
!
!
!  Authors : Xu Junqing and D. Sébilleau
!
!
!                                      Last modified (DS) : 31 May 2021
!
!
      USE REAL_NUMBERS,         ONLY : ZERO,ONE,TWO
      USE COMPLEX_NUMBERS,      ONLY : ZEROC,IC
!
      USE ANA_DIR
      USE CLUSTER_COORD
      USE CURRENT_BEAM
      USE DEB_WAL_CLU
      USE INIT_L
      USE SPHERICAL_HARMONICS
      USE CURRENT_T_MATRIX
!
      IMPLICIT NONE
!
      CHARACTER (LEN =  3)      ::  TAU_TYPE
      CHARACTER (LEN =  2)      ::  ALGO
!
      INTEGER, INTENT(IN)       ::  JE,NA_TYP(NATM),NTYPEM
      INTEGER, INTENT(IN)       ::  JBL,JATBL,JDIR
!
      INTEGER                   ::  L,L1,L2,ISTEP_L
      INTEGER                   ::  LME,IL,IND,JATL,JTYP,JNUM
      INTEGER                   ::  LMJ,NBTYP,M
      INTEGER                   ::  LJ,MJ,ILJ,INDJ
!
      REAL (WP), INTENT(IN)     ::  VALS,ZSURF,GAMMA
!
      REAL (WP)                 ::  ATTSE,ATTSJ
      REAL (WP)                 ::  XOJ,YOJ,ZOJ,ROJ
      REAL (WP)                 ::  ZSURFJ
      REAL (WP)                 ::  DWTER,CSTHJR,COS_TH
      REAL (WP)                 ::  SSMALL,CSKZ2J,UJJ
!
      REAL (WP)                 ::  SQRT,ABS
!
      COMPLEX (WP), INTENT(IN)  ::  TAU(LINMAX,LINFMAX,NATCLU_M)
      COMPLEX (WP), INTENT(OUT) ::  AL(LINMAX),BL(LINMAX)
!
      COMPLEX (WP)              ::  YLM(0:NL_M,-NL_M:NL_M)
      COMPLEX (WP)              ::  YLME(0:NL_M,-NL_M:NL_M)
      COMPLEX (WP)              ::  ATT_M
      COMPLEX (WP)              ::  SLJ,SUM_AT,SUM_J
      COMPLEX (WP)              ::  TL_AT
!
      COMPLEX (WP)              ::  CONJG
!
      DATA SSMALL / 0.0001E0_WP /
!
!  Initialization of damping values and distance to surface
!
      ATTSE  = ONE                                                  !
      ATTSJ  = ONE                                                  !
      ZSURFJ = ZERO                                                 !
!
!  Initialization of AL and BL arrays
!
      DO L = 1, LINMAX                                              !
        AL(L) = ZEROC                                               !
        BL(L) = ZEROC                                               !
      END DO                                                        !
!
!  Damping of the direct electron wave between atom JATBL
!    and the surface plane
!
      IF(IATTS == 1) THEN                                           !
        ATTSE = EXP(- ZSURF * GAMMA / DIRANA(3,JDIR))               !
      END IF                                                        !
!
!  Variables dependent on the atom where BL is computed
!
      IF(JBL == NTYPEM) THEN                                        !
        L1      = LF1                                               !
        L2      = LF2                                               !
        ISTEP_L = ISTEP_LF                                          !
      ELSE                                                          !
        L1      = 0                                                 !
        L2      = LMAX(JBL,JE)                                      !
        ISTEP_L = 1                                                 !
      END IF                                                        !
!
!  The atomic contribution (not contained in series
!   expansion calculation of Tau)
!
      LME = LMAX(JBL,JE)                                            !
      CALL HARSPH(NL_M,RTH_IN(JDIR),RPH_IN(JDIR),YLME,LME)          !
!
!  Loop on angular momentum on atom where BL is computed
!
      DO L = L1, L2, ISTEP_L                                        !
        IL    = L * L + L + 1                                       !
        TL_AT = TL(L,1,JBL,JE)                                      !
        DO M = - L, L                                               !
          IND = IL + M                                              !
!
!  Inner loop on the atoms J of the cluster
!
          JATL = 0                                                  !
          SUM_J = ZEROC                                             !
          DO JTYP = 1, N_PROT                                       !
            NBTYP = NA_TYP(JTYP)                                    !
            LMJ   = LMAX(JTYP,JE)                                   !
            DO JNUM = 1, NBTYP                                      !
              JATL = JATL + 1                                       ! atom number in the loop
!
              IF(JATL == JATBL) GO TO 10                            !
!
              XOJ    = SYM_AT(1,JATL) - SYM_AT(1,JATBL)             !
              YOJ    = SYM_AT(2,JATL) - SYM_AT(2,JATBL)             !
              ZOJ    = SYM_AT(3,JATL) - SYM_AT(3,JATBL)             !
              ROJ    = SQRT(XOJ * XOJ + YOJ * YOJ + ZOJ * ZOJ)      !
!
              ZSURFJ = VALS - SYM_AT(3,JATL)                        !
!
!  Damping of the direct electron wave between atom JATL
!    and the surface plane
!
              IF(IATTS == 1) THEN                                   !
                ATTSJ = EXP(- ZSURFJ * GAMMA / DIRANA(3,JDIR))      !
              END IF                                                !
!
  10          CONTINUE                                              !
!
              CALL HARSPH(NL_M,RTH_IN(JDIR),RPH_IN(JDIR),YLM,LMJ)   !
!
              IF(JATL == JATBL) THEN                                !
                DWTER = ONE                                         !
                GO TO 30                                            !
              ELSE IF(JATL /= JATBL) THEN                           !
                CSTHJR = ( XOJ * DIRANA(1,JDIR) +                 & !
                           YOJ * DIRANA(2,JDIR) +                 & !
                           ZOJ * DIRANA(3,JDIR)                   & !
                         ) / ROJ                                    !
                IF(IDWSPH == 1) GO TO 20                            !
                COS_TH = ZOJ / ROJ                                  !
                IF(COS_TH > ONE)  THEN                              !
                  COS_TH =   ONE                                    !
                ELSE IF(COS_TH < - ONE) THEN                        !
                  COS_TH = - ONE                                    !
                END IF                                              !
                IF(ABS(ZSURFJ) <= SSMALL) THEN                      !
                  IF(ABS(CSTHJR - ONE) > SSMALL) THEN               !
                    CSKZ2J = ( DIRANA(3,JDIR) - COS_TH ) *        & !
                             ( DIRANA(3,JDIR) - COS_TH ) /        & !
                             (TWO - TWO *CSTHJR)                    !
                  ELSE                                              !
                    CSKZ2J = ONE                                    !
                  END IF                                            !
                  UJJ = UJ2(JTYP) * ( ONE + CSKZ2J * (RSJ - ONE) )  !
                ELSE                                                !
                  UJJ = UJ2(JTYP)                                   !
                END IF                                              !
              END IF                                                !
!
  20          IF(IDWSPH == 1) THEN                                  !
                DWTER = ONE                                         !
              ELSE                                                  !
                DWTER = EXP( - VK2(JE) * UJJ * (ONE - CSTHJR) )     !
              END IF                                                !
!
  30          IF(TAU_TYPE == 'OUT') THEN                            !
                IF(JATL == JATBL) THEN                              !
                  ATT_M = ATTSE * DWTER                             !
                ELSE                                                !
                  ATT_M = ATTSJ * DWTER *                         & !
                          EXP( - IC * VK(JE) * ROJ * CSTHJR )       !
                END IF                                              !
              ELSE IF(TAU_TYPE == 'INC') THEN                       !
                IF(JATL == JATBL) THEN                              !
                  ATT_M = ATTSE * DWTER                             !
                ELSE                                                !
                  ATT_M = ATTSJ * DWTER *                         & !
                            EXP( IC * VK(JE) * ROJ * CSTHJR )       !
                END IF                                              !
              END IF                                                !
!
!  Inner loop on the angular momentum of atoms J of the cluster
!
              SLJ = ZEROC                                           !
              DO LJ = 0, LMJ                                        !
                ILJ = LJ * LJ + LJ + 1                              !
                DO MJ = - LJ, LJ                                    !
                  INDJ = ILJ + MJ                                   !
!
                  IF(TAU_TYPE == 'OUT') THEN                        !
                    SLJ = SLJ + TAU(INDJ,IND,JATL) * YLM(LJ,MJ)     !
                  ELSE IF(TAU_TYPE == 'INC') THEN                   !
                    SLJ = SLJ + TAU(INDJ,IND,JATL) *              & !
                                CONJG(YLM(LJ,MJ))                   !
                  END IF                                            !
!
!  End of inner loop on the angular momentum of atoms J of the cluster
!
                END DO                                              !
              END DO                                                !
!
              SUM_J = SUM_J + SLJ * ATT_M                           !
!
!  End of inner loop on the atoms of the cluster
!
            END DO                                                  !
          END DO                                                    !
!
!  Computation of atomic (direct) contribution: AL  <-- 0th order scattering
!
          IF(JBL == NTYPEM) THEN                                    !
            SUM_AT  = YLME(L,M) * ATTSE * TL_AT                     !
            AL(IND) = SUM_AT                                        !
          END IF                                                    !
!
!  Computation of the MS amplitude: BL
!
          IF(ALGO == 'SE') THEN                                     !
            BL(IND) = SUM_J + SUM_AT                                !
          ELSE IF(ALGO == 'RE') THEN                                !
            BL(IND) = SUM_J + SUM_AT                                !
          ELSE IF(ALGO == 'MI') THEN                                !
            BL(IND) = SUM_J                                         !
          ELSE IF(ALGO == 'CE') THEN                                !
            BL(IND) = SUM_J                                         !
          END IF
!
!  End of loop on angular momentum on atom where BL is computed
!
        END DO                                                      !
      END DO                                                        !
!
      END SUBROUTINE BL_AMP
!
END MODULE BEAM_AMPLITUDE
