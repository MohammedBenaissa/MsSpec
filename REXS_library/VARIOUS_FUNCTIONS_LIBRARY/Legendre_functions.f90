!
!=======================================================================
!
MODULE LEGENDRE_FUNCTIONS 
!
!  This module provides Legendre polynomials and functions
!
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE POLLEG(NC,X,PL)                                        
!
!  This routine computes the Legendre polynomials up to order NC
!    using the standard Bonnet recurrence:
!
!       (n+1) P_(n+1)(x) = (2n+1)x P_(n)(x) - n P_(n-1)(x)
!
!   starting from P_0(x) = 1
!                 P_1(x) = x
!
!
!   Author :  D. Sébilleau
!
!                                           Last modified : 14 Aug 2020
!
      USE REAL_NUMBERS,        ONLY : ONE
!
      IMPLICIT NONE
!
      REAL (WP)             ::  X
      REAL (WP)             ::  XL,XL1,XL3
      REAL (WP)             ::  PL(0:150) 
!
      INTEGER               ::  NC,L
      INTEGER               ::  L1,L2
      INTEGER               ::  LOGF
!
      LOGF = 6                                                      !
!
      IF(NC > 150) THEN                                             !
        WRITE(LOGF,10)                                              !
        STOP                                                        !
      END IF                                                        !
!
      PL(0) = ONE                                                   !
      PL(1) = X                                                     !                                                        
!
      DO L=2,NC                                                     !
!
        L1    = L - 1                                               !
        L2    = L - 2                                               !
        XL    = DFLOAT(L)                                           ! L
        XL1   = DFLOAT(L1)                                          ! L+1
        XL3   = XL+XL+ONE                                           ! 2L+1                                                        
        PL(L) = (X*XL3*PL(L1)-XL1*PL(L2))/XL                        !
!
!  Format:
!
  10  FORMAT(5X,'<<<<<  DIMENSION ERROR IN POLLEG  >>>>>',/,      &
             5X,'<<<<<  L > 150. RE-DIMENSION PL   >>>>>',//)
!
      END DO                                                        !
!
      END SUBROUTINE POLLEG                                                              
!
!=======================================================================
!
      SUBROUTINE PLM(NC,X,PLMM)
!
!  This routine computes the associated Legendre functions 
!    of the first kind. It is a modified version of that written by
!
!    W.H. Press, B.P. Flannery, S.A. Teukolsky and W.T. Vetterling
!       in "Numerical Recipes : The Art of Scientific Computing" 
!               (Cambridge University Press 1992).
!
!  It computes all values of P_l^m(x) up to l = NC 
!       and stores them as PLMM(L,M)
!
!
!  Input variables :
!
!       * NC       : upper value of l
!       * X        : argument of P_l^m
!
!  Output variables :
!
!       * PLMM     : P_l^m(x) for l = 0 to l = NC
!
!
!                                           Last modified : 14 Aug 2020
!
!
      USE REAL_NUMBERS,        ONLY : ONE,TWO
!
      IMPLICIT NONE
!
      REAL (WP)             ::  X
      REAL (WP)             ::  PLMM(0:150,0:150)
      REAL (WP)             ::  PMM,FACT,SOMX2,PMMP1,PLL
!
      INTEGER               ::  NC,L,I,M
      INTEGER               ::  LOGF
!
      LOGF = 6                                                      !
!
!  Initialization with Legendre polynomials PLMM(L,0)
!                (recurrence on L)
!
      PLMM(0,0) = ONE                                               !
      PLMM(1,0) = X                                                 ! 
!
      DO L=2,NC                                                     !
        PLMM(L,0)=( X * DFLOAT(L+L-1) * PLMM(L-1,0) -             & !
                        DFLOAT(L-1)   * PLMM(L-2,0)               & !
                  )   / DFLOAT(L)                                   !
      END DO                                                        !
!
      DO M=1,NC                                                     !
!
        PMM   = ONE                                                 !                                                       
        FACT  = ONE                                                 !                                                     
        SOMX2 = DSQRT(ONE - X*X)                                    !                                     
        FACT  = ONE                                                 !
!
        DO I=1,M                                                    !                                               
          PMM  = -PMM * FACT * SOMX2                                !                                         
          FACT =  FACT+TWO                                          !                                     
        END DO                                                      ! 
!
        PMMP1       = X* FACT * PMM                                 !
        PLMM(M,M)   = PMM                                           !
        PLMM(M+1,M) = PMMP1                                         !
!
        IF(M < NC-1) THEN                                           !
!
          DO L=M+2,NC                                               !                                            
            PLL=( X*DFLOAT(L+L-1) * PMMP1 -                       & ! 
                    DFLOAT(L+M-1) * PMM                           & ! 
                ) / DFLOAT(L-M)                                     !            
            PMM       = PMMP1                                       !                                             
            PMMP1     = PLL                                         !                                       
            PLMM(L,M) = PLL                                         ! 
!
          END DO                                                    !  
!
        END IF                                                      ! 
!
      END DO                                                        ! 
!
      END SUBROUTINE PLM
!
!=======================================================================
!
      SUBROUTINE LQMN(MM,M,N,X,QM,QD)
!
!
!       ==========================================================
!       Purpose: Compute the associated Legendre functions of the
!                second kind, Qmn(x) and Qmn'(x)
!       Input :  x  --- Argument of Qmn(x) 
!                m  --- Order of Qmn(x)  ( m = 0,1,2,...)
!                n  --- Degree of Qmn(x) ( n = 0,1,2,...)
!                mm --- Physical dimension of QM and QD
!       Output:  QM(m,n) --- Qmn(x)
!                QD(m,n) --- Qmn'(x)
!       ==========================================================
!
!  From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
!
!  The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."
!
!
!                                     Last modified (DS) : 14 Aug 2020
!
!
      USE REAL_NUMBERS,        ONLY : ZERO,ONE,TWO,THREE,HALF,INF
!
      IMPLICIT NONE
!  
      INTEGER               ::  MM,M,N
      INTEGER               ::  LS,I,J,K,KM
!
      REAL (WP)             ::  X
      REAL (WP)             ::  QM(0:MM,0:N),QD(0:MM,0:N)
      REAL (WP)             ::  XS,XQ,Q0,Q1,Q10,QF
      REAL (WP)             ::  XI,XJ,XK
      REAL (WP)             ::  QF0,QF1,QF2
!
!  Trivial cas X = 1:
!
      IF (DABS(X) == ONE) THEN                                      !
!
        DO I=0,M                                                    !
          DO J=0,N
            QM(I,J) = INF                                           !
            QD(I,J) = INF                                           !
          END DO                                                    !
        END DO                                                      !
!
        RETURN                                                      !
!
      END IF                                                        !
!
      LS = 1                                                        !
      IF(DABS(X) > ONE) LS = -1                                     !
      XS = LS * (ONE - X*X)                                         !
      XQ = DSQRT(XS)                                                !
      Q0 = HALF * DLOG(DABS((X+ONE) / (X-ONE)))                     !
!
      IF(DABS(X) < 1.0001E0_WP) THEN                                !
        QM(0,0) = Q0                                                !
        QM(0,1) = X * Q0 - ONE                                      !
        QM(1,0) = -ONE / XQ                                         !
        QM(1,1) = -XQ * (Q0 + X / (ONE - X*X))                      !
!
        DO I=0,1                                                    !
          XI = DFLOAT(I)                                            !
          DO J=2,N                                                  !
            XJ      = DFLOAT(J)                                     !
            QM(I,J )= ( (TWO*XJ-ONE) * X * QM(I,J-1)              & !
                       -(XJ+XI-ONE)*QM(I,J-2)                     & !
                      ) / (XJ-XI)                                   !
          END DO                                                    !
        END DO                                                      !
!
        DO J=0,N
          XJ = DFLOAT(J)                                            !
          DO I=2,M
            XI      = DFLOAT(I)                                     !
            QM(I,J) = -TWO*(XI-ONE) * X / XQ * QM(I-1,J) - LS *   & !
                       (XJ+XI-ONE) * (XJ-XI+TWO) * QM(I-2,J)        !
          END DO                                                    !
        END DO                                                      !
!
      ELSE                                                          !
!
        IF(DABS(X) > 1.1E0_WP) THEN                                 !
          KM = 40 + M + N
        ELSE
          KM = (40 + M + N) * INT(- ONE - 1.8E0_WP * DLOG(X-ONE))   !
        END IF                                                      !
!
        QF2 = ZERO                                                  !
        QF1 = ONE                                                   !
!
        DO K=KM,0,-1                                                !
          XK  = DFLOAT(K)                                           !
          QF0 = ( (XK + XK + THREE)*X*QF1-(XK+TWO)*QF2 ) / (XK+ONE) !
          IF(K <= N) QM(0,K) = QF0                                  !
          QF2 = QF1                                                 !
          QF1 = QF0                                                 !
        END DO                                                      !
!
        DO K=0,N                                                    !
          QM(0,K) = Q0 * QM(0,K) / QF0                              !
        END DO                                                      !
!
        QF2 = ZERO                                                  !
        QF1 = ONE                                                   !
!
        DO K=KM,0,-1
          XK  = DFLOAT(K)                                           !
          QF0 = ( (XK + XK + THREE)*X*QF1-(XK+ONE)*QF2 ) / (XK+TWO) !
          IF(K <= N) QM(1,K) = QF0                                  !
          QF2 = QF1                                                 !
          QF1 = QF0                                                 !
        END DO                                                      !
!
        Q10 = -ONE / XQ                                             !
!
        DO K=0,N                                                    !
          QM(1,K) = Q10 * QM(1,K) / QF0                             !
        END DO                                                      !
!
        DO J=0,N                                                    !
          XJ = DFLOAT(J)                                            !
          Q0 = QM(0,J)                                              !
          Q1 = QM(1,J)                                              !
          DO I=0,M-2                                                !
            XI        = DFLOAT(I)                                   !
            QF        = -TWO*(XI+1)*X/XQ*Q1+(XJ-XI)*(XJ+XI+ONE)*Q0  !
            QM(I+2,J) = QF                                          !
            Q0        = Q1                                          !
            Q1        = QF                                          !
          END DO                                                    !  
        END DO                                                      !
!
      END IF                                                        !
!
      QD(0,0) = DFLOAT(LS) / XS                                     !
      DO J=1,N                                                      !
        QD(0,J) = LS * J * ( QM(0,J-1) - X * QM(0,J) ) / XS         !
      END DO                                                        !
      DO J=0,N                                                      !
        XJ = DFLOAT(J)                                              !
        DO I=1,M                                                    !
          XI      = DFLOAT(I)                                       !
          QD(I,J) = LS*XI*X/XS*QM(I,J) + (XI+XJ)*(XJ-XI+ONE) /    & !
                    XQ*QM(I-1,J)                                    !
        END DO                                                      !
      END DO                                                        !
!
      END SUBROUTINE LQMN
!
END MODULE LEGENDRE_FUNCTIONS
