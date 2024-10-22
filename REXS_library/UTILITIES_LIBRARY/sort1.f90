!
!=======================================================================
!
MODULE SORT1
!
!  This module provides several routines for sorting REAL arrays
!
!
!   1) SUBROUTINE SORT_INC(NINI,VALIN,NFIN,VALFIN)
!
!   2) SUBROUTINE SORT_DEC(NINI,VALIN,NFIN,VALFIN)
!
!   3) SUBROUTINE SORT_1(ARRV,N)
!
!   4) SUBROUTINE SSORT(X,Y,N,KFLAG)
!
!   5) SUBROUTINE QUICK_SORT(LIST,ORDER)
!
!
      USE ACCURACY_REAL
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE SORT_INC(NINI,VALIN,NFIN,VALFIN)
!
!   Given a set of **real** numbers VALINI, this routine orders them and
!       suppresses the values appearing more than once. The remaining
!       values are stored in VALFIN.
!
!       VALINI(K+1) > VALINI(K) : decreasing order
!       VALINI(K+1) < VALINI(K) : increasing order
!
!
!                                     Last modified (DS) :  1 Sep 2020
!
!
      USE REAL_NUMBERS,     ONLY : SMALL
!
      IMPLICIT NONE
!
      INTEGER, INTENT (IN)   ::  NINI
      INTEGER, INTENT (OUT)  ::  NFIN
!
      INTEGER                ::  I,J,K
      INTEGER                ::  JFIN
!
      REAL (WP), INTENT(IN)  ::  VALIN(NINI)
      REAL (WP), INTENT(OUT) ::  VALFIN(NINI)
!
      REAL (WP)              ::  VALINI(NINI)
      REAL (WP)              ::  R1
!
      LOGICAL                ::  BUBBLE
!
!...Store input array
!
      DO I = 1,NINI                                                 !
        VALINI(I) = VALIN(I)                                        !
      END DO                                                        !
!
      DO J = 1,NINI-1                                               !
        K      = J                                                  !
        BUBBLE = .TRUE.                                             !
!
  150   IF(K >= 1 .AND. BUBBLE) THEN                                !
          IF(VALINI(K+1) < VALINI(K)) THEN                          !
            R1          = VALINI(K)                                 !
            VALINI(K)   = VALINI(K+1)                               !
            VALINI(K+1) = R1                                        !
          ELSE                                                      !
            BUBBLE = .FALSE.                                        !
          END IF                                                    !
          K = K - 1                                                 !
          GO TO 150                                                 !
        END IF                                                      !
      END DO                                                        !
!
      JFIN      = 1                                                 !
      VALFIN(1) = VALINI(1)                                         !
      DO J = 1,NINI-1                                               !
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)) > SMALL) THEN              !
          JFIN         = JFIN+1                                     !
          VALFIN(JFIN) = VALINI(J+1)                                !
        END IF                                                      !
      END DO                                                        !
      NFIN = JFIN                                                   !
!
      RETURN                                                        !
!
      END SUBROUTINE SORT_INC
!
!=======================================================================
!
      SUBROUTINE SORT_DEC(NINI,VALIN,NFIN,VALFIN)
!
!   Given a set of **real** numbers VALINI, this routine orders them and
!       suppresses the values appearing more than once. The remaining
!       values are stored in VALFIN.
!
!       VALINI(K+1) > VALINI(K) : decreasing order
!       VALINI(K+1) < VALINI(K) : increasing order
!
!
!                                     Last modified (DS) :  1 Sep 2020
!
!
      USE REAL_NUMBERS,     ONLY : SMALL
!
      IMPLICIT NONE
!
      INTEGER, INTENT (IN)   ::  NINI
      INTEGER, INTENT (OUT)  ::  NFIN
!
      INTEGER                ::  I,J,K
      INTEGER                ::  JFIN
!
      REAL (WP), INTENT(IN)  ::  VALIN(NINI)
      REAL (WP), INTENT(OUT) ::  VALFIN(NINI)
!
      REAL (WP)              ::  VALINI(NINI)
      REAL (WP)              ::  R1
!
      LOGICAL                ::  BUBBLE
!
!...Store input array
!
      DO I = 1,NINI                                                 !
        VALINI(I) = VALIN(I)                                        !
      END DO                                                        !
!
      DO J = 1,NINI-1                                               !
        K      = J                                                  !
        BUBBLE = .TRUE.                                             !
!
  150   IF(K >= 1 .AND. BUBBLE) THEN                                !
          IF(VALINI(K+1) > VALINI(K)) THEN                          !
            R1          = VALINI(K)                                 !
            VALINI(K)   = VALINI(K+1)                               !
            VALINI(K+1) = R1                                        !
          ELSE                                                      !
            BUBBLE = .FALSE.                                        !
          END IF                                                    !
          K = K - 1                                                 !
          GO TO 150                                                 !
        END IF                                                      !
      END DO                                                        !
!
      JFIN      = 1                                                 !
      VALFIN(1) = VALINI(1)                                         !
      DO J = 1,NINI-1                                               !
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)) > SMALL) THEN              !
          JFIN         = JFIN+1                                     !
          VALFIN(JFIN) = VALINI(J+1)                                !
        END IF                                                      !
      END DO                                                        !
      NFIN = JFIN                                                   !
!
      RETURN                                                        !
!
      END SUBROUTINE SORT_DEC
!
!=======================================================================
!
      SUBROUTINE SORT_1(ARRV,N)
!
!     Sorting elements of array ARRV in ascending order. ARRV is overwritten by
!     reordered elements.
!
!                                     Last modified (DS) :  1 Sep 2020
!
!
      USE REAL_NUMBERS,     ONLY : SMALL
!
      IMPLICIT NONE
!
      INTEGER, INTENT (IN)   ::  N
!
      INTEGER                ::  I,J,IPTR
!
      REAL (WP)              ::  ARRV(N)
      REAL (WP)              ::  TEMP
!
      DO I = 1, N-1                                                 !
!
!   Find the maximum value in ARR(I) through ARR(N)
!
         IPTR = I                                                   !
         DO J = I+1, N                                              !
            IF(ARRV(J) < ARRV(IPTR)) THEN                           !
               IPTR = J                                             !
            END IF                                                  !
         END DO                                                     !
!
!   IPTR now points to the minimum value, so swap ARR(IPTR) with ARR(I)
!
         IF(I /= IPTR) THEN                                         !
            TEMP       = ARRV(I)                                    !
            ARRV(I)    = ARRV(IPTR)                                 !
            ARRV(IPTR) = TEMP                                       !
         END IF                                                     !
      END DO                                                        !
!
      END SUBROUTINE SORT_1
!
!=======================================================================
!
      SUBROUTINE SSORT(X,Y,N,KFLAG)
!
!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters:
!
!      X - array of values to be sorted   (usually abscissas)
!      Y - array to be (optionally) carried along
!      N - number of values in array X to be sorted
!      KFLAG - control parameter
!            =  2  means sort X in increasing order and carry Y along.
!            =  1  means sort X in increasing order (ignoring Y)
!            = -1  means sort X in decreasing order (ignoring Y)
!            = -2  means sort X in decreasing order and carry Y along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  SSORT
!
!                                     Last modified (DS) :  1 Sep 2020
!
!
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)        ::  N
      INTEGER, INTENT(IN)        ::  KFLAG
!
      INTEGER                    ::  NN,KK
      INTEGER                    ::  I,J,K,L,M,IJ
!
      REAL (WP)                  ::  X(*),Y(*)
      REAL (WP)                  ::  R,T,TT,TY,TTY
      REAL (WP)                  ::  IL(21), IU(21)
!
!***FIRST EXECUTABLE STATEMENT SSORT
!
      NN = N                                                        !
      IF(NN < 1) THEN                                               !
         WRITE(6,10)                                                !
         RETURN                                                     !
      END IF                                                        !
!
      KK = ABS(KFLAG)                                               !
      IF(KK /= 1 .AND. KK /= 2) THEN                                !
         WRITE(6,11)
         RETURN                                                     !
      END IF                                                        !
!
!     Alter array X to get decreasing order if needed
!
      IF(KFLAG <= -1) THEN                                          !
         DO I = 1,NN                                                !
            X(I) = - X(I)                                           !
         END DO                                                     !
      END IF                                                        !
!
      IF(KK == 2) GO TO 100                                         !
!
!     Sort X only
!
      M = 1                                                         !
      I = 1                                                         !
      J = NN                                                        !
      R = 0.375E0_WP                                                !
!
   20 IF(I == J) GO TO 60                                           !
      IF(R <= 0.5898437E0_WP) THEN                                  !
         R = R + 3.90625E-2_WP                                      !
      ELSE                                                          !
         R = R - 0.21875E0_WP                                       !
      END IF                                                        !
!
   30 K = I                                                         !
!
!     Select a central element of the array and save it in location T
!
      IJ = I + INT((J - I) * R)                                     !
      T  = X(IJ)                                                    !
!
!     If first element of array is greater than T, interchange with T
!
      IF(X(I) > T) THEN                                             !
         X(IJ) = X(I)                                               !
         X(I)  = T                                                  !
         T     = X(IJ)                                              !
      END IF                                                        !
      L = J                                                         !
!
!     If last element of array is less than than T, interchange with T
!
      IF(X(J) < T) THEN                                             !
         X(IJ) = X(J)                                               !
         X(J)  = T                                                  !
         T     = X(IJ)                                              !
!
!        If first element of array is greater than T, interchange with T
!
         IF(X(I) > T) THEN                                          !
            X(IJ) = X(I)                                            !
            X(I)  = T                                               !
            T     = X(IJ)                                           !
         END IF                                                     !
      END IF                                                        !
!
!     Find an element in the second half of the array which is smaller
!     than T
!
   40 L = L - 1                                                     !
      IF(X(L) > T) GO TO 40                                         !
!
!     Find an element in the first half of the array which is greater
!     than T
!
   50 K = K + 1                                                     !
      IF(X(K) < T) GO TO 50                                         !
!
!     Interchange these elements
!
      IF(K <= L) THEN                                               !
         TT   = X(L)                                                !
         X(L) = X(K)                                                !
         X(K) = TT                                                  !
         GO TO 40                                                   !
      END IF                                                        !
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF(L-I > J-K) THEN                                            !
         IL(M) = I                                                  !
         IU(M) = L                                                  !
         I     = K                                                  !
         M     = M + 1                                              !
      ELSE                                                          !
         IL(M) = K                                                  !
         IU(M) = J                                                  !
         J     = L                                                  !
         M     = M + 1                                              !
      END IF                                                        !
      GO TO 70                                                      !
!
!     Begin again on another portion of the unsorted array
!
   60 M = M - 1                                                     !
      IF(M == 0) GO TO 190                                          !
      I = IL(M)                                                     !
      J = IU(M)                                                     !
!
   70 IF(J-I >= 1) GO TO 30
      IF(I == 1) GO TO 20
      I = I-1
!
   80 I = I + 1                                                     !
      IF(I == J) GO TO 60                                           !
      T = X(I+1)                                                    !
      IF(X(I) <= T) GO TO 80                                        !
      K = I                                                         !
!
   90 X(K+1) = X(K)                                                 !
      K      = K - 1                                                !
      IF(T < X(K)) GO TO 90                                         !
      X(K+1) = T                                                    !
      GO TO 80                                                      !
!
!     Sort X and carry Y along
!
  100 M = 1                                                         !
      I = 1                                                         !
      J = NN                                                        !
      R = 0.375E0_WP                                                !
!
  110 IF(I == J) GO TO 150                                          !
      IF(R <= 0.5898437E0_WP) THEN                                  !
         R = R + 3.90625E-2_WP                                      !
      ELSE                                                          !
         R = R - 0.21875E0_WP                                       !
      END IF                                                        !
!
  120 K = I                                                         !
!
!     Select a central element of the array and save it in location T
!
      IJ = I + INT((J - I) * R)                                     !
      T  = X(IJ)                                                    !
      TY = Y(IJ)                                                    !
!
!     If first element of array is greater than T, interchange with T
!
      IF(X(I) > T) THEN                                             !
         X(IJ) = X(I)                                               !
         X(I)  = T                                                  !
         T     = X(IJ)                                              !
         Y(IJ) = Y(I)                                               !
         Y(I)  = TY                                                 !
         TY    = Y(IJ)                                              !
      END IF                                                        !
      L = J                                                         !
!
!     If last element of array is less than T, interchange with T
!
      IF(X(J) < T) THEN                                             !
         X(IJ) = X(J)                                               !
         X(J)  = T                                                  !
         T     = X(IJ)                                              !
         Y(IJ) = Y(J)                                               !
         Y(J)  = TY                                                 !
         TY    = Y(IJ)                                              !
!
!        If first element of array is greater than T, interchange with T
!
         IF(X(I) > T) THEN                                          !
            X(IJ) = X(I)                                            !
            X(I)  = T                                               !
            T     = X(IJ)                                           !
            Y(IJ) = Y(I)                                            !
            Y(I)  = TY                                              !
            TY    = Y(IJ)                                           !
         END IF                                                     !
      END IF                                                        !
!
!     Find an element in the second half of the array which is smaller
!     than T
!
  130 L = L - 1                                                     !
      IF(X(L) > T) GO TO 130                                        !
!
!     Find an element in the first half of the array which is greater
!     than T
!
  140 K = K + 1                                                     !
      IF(X(K) < T) GO TO 140                                        !
!
!     Interchange these elements
!
      IF(K <= L) THEN                                               !
         TT   = X(L)                                                !
         X(L) = X(K)                                                !
         X(K) = TT                                                  !
         TTY  = Y(L)                                                !
         Y(L) = Y(K)                                                !
         Y(K) = TTY                                                 !
         GO TO 130                                                  !
      END IF                                                        !
!
!     Save upper and lower subscripts of the array yet to be sorted
!
      IF(L-I > J-K) THEN                                            !
         IL(M) = I                                                  !
         IU(M) = L                                                  !
         I     = K                                                  !
         M     = M + 1                                              !
      ELSE                                                          !
         IL(M) = K                                                  !
         IU(M) = J                                                  !
         J     = L                                                  !
         M     = M + 1                                              !
      END IF
      GO TO 160                                                     !
!
!     Begin again on another portion of the unsorted array
!
  150 M = M - 1                                                     !
      IF(M == 0) GO TO 190                                          !
      I = IL(M)                                                     !
      J = IU(M)                                                     !
!
  160 IF(J-I >= 1) GO TO 120                                        !
      IF(I == 1) GO TO 110                                          !
      I = I - 1                                                     !
!
  170 I = I + 1                                                     !
      IF(I == J) GO TO 150                                          !
      T  = X(I+1)                                                   !
      TY = Y(I+1)                                                   !
      IF(X(I) <= T) GO TO 170                                       !
      K = I                                                         !
!
  180 X(K+1) = X(K)                                                 !
      Y(K+1) = Y(K)                                                 !
      K = K - 1                                                     !
      IF(T < X(K)) GO TO 180                                        !
      X(K+1) = T                                                    !
      Y(K+1) = TY                                                   !
      GO TO 170                                                     !
!
!     Clean up
!
  190 IF(KFLAG <= -1) THEN                                          !
         DO I = 1,NN                                                !
            X(I) = - X(I)                                           !
        END DO                                                      !
      END IF                                                        !
!
      RETURN                                                        !
!
!  Format:
!
  10  FORMAT('The number of values to be sorted is not positive.')
  11  FORMAT('The sort control parameter, K, is not 2, 1, -1, or -2.')
!
      END SUBROUTINE SSORT
!
!=======================================================================
!
      RECURSIVE SUBROUTINE QUICK_SORT(LIST,ORDER)
!
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
!
!                                     Last modified (DS) :  1 Sep 2020
!
!
      IMPLICIT NONE
!
      REAL, DIMENSION (:), INTENT(IN OUT)  :: LIST
!
      INTEGER, DIMENSION (:), INTENT(OUT)  :: ORDER
!
      INTEGER                              :: I
!
      DO I = 1, SIZE(LIST)                                          !
        ORDER(I) = I                                                !
      END DO                                                        !
!
      CALL QUICK_SORT_1(1,SIZE(LIST))                               !
!
CONTAINS
!
!-----------------------------------------------------------------------
!
      RECURSIVE SUBROUTINE QUICK_SORT_1(LEFT_END,RIGHT_END)
!
      INTEGER, INTENT(IN) :: LEFT_END,RIGHT_END
!
!     Local variables
!
      INTEGER             :: I,J,ITEMP
      REAL                :: REFERENCE,TEMP
      INTEGER, PARAMETER  :: MAX_SIMPLE_SORT_SIZE = 6
!
      IF(RIGHT_END < LEFT_END + MAX_SIMPLE_SORT_SIZE) THEN          !
!
! Use interchange sort for small lists
!
        CALL INTERCHANGE_SORT(LEFT_END, RIGHT_END)                  !
!
      ELSE
!
! Use partition ("quick") sort
!
        REFERENCE = LIST((LEFT_END + RIGHT_END)/2)                  !
        I         = LEFT_END  - 1                                   !
        J         = RIGHT_END + 1                                   !
!
        DO
!
! Scan list from left end until element >= reference is found
!
          DO                                                        !
            I = I + 1                                               !
            IF(LIST(I) >= REFERENCE) GO TO 10                       !
          END DO                                                    !
  10      CONTINUE                                                  !
!
! Scan list from right end until element <= reference is found
!
          DO                                                        !
            J = J - 1                                               !
            IF(LIST(J) <= REFERENCE) GO TO 20                       !
          END DO                                                    !
  20      CONTINUE                                                  !
!
          IF (I < J) THEN                                           !
!
! Swap two out-of-order elements
!
            TEMP     = LIST(I)                                      !
            LIST(I)  = LIST(J)                                      !
            LIST(J)  = TEMP                                         !
            ITEMP    = ORDER(I)                                     !
            ORDER(I) = ORDER(J)                                     !
            ORDER(J) = ITEMP                                        !
          ELSE IF (I == J) THEN                                     !
            I = I + 1                                               !
            GO TO 30                                                !
          ELSE                                                      !
            GO TO 30                                                !
          END IF                                                    !
        END DO                                                      !
  30    CONTINUE                                                    !
!
        IF(LEFT_END < J)  CALL QUICK_SORT_1(LEFT_END,J)             !
        IF(I < RIGHT_END) CALL QUICK_SORT_1(I,RIGHT_END)            !
      END IF                                                        !
!
      END SUBROUTINE QUICK_SORT_1
!
!-----------------------------------------------------------------------
!
      SUBROUTINE INTERCHANGE_SORT(LEFT_END, RIGHT_END)
!
      INTEGER, INTENT(IN) :: LEFT_END, RIGHT_END
!
!     Local variables
!
      INTEGER             :: I, J, ITEMP
      REAL                :: TEMP
!
      DO I = LEFT_END, RIGHT_END - 1                                !
        DO J = I+1, RIGHT_END                                       !
          IF (LIST(I) > LIST(J)) THEN                               !
            TEMP     = LIST(I)                                      !
            LIST(I)  = LIST(J)                                      !
            LIST(J)  = TEMP                                         !
            ITEMP    = ORDER(I)                                     !
            ORDER(I) = ORDER(J)                                     !
            ORDER(J) = ITEMP                                        !
          END IF                                                    !
        END DO                                                      !
      END DO                                                        !
!
      END SUBROUTINE INTERCHANGE_SORT
!
      END SUBROUTINE QUICK_SORT
!
END MODULE SORT1
