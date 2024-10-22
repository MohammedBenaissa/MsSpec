!
!=======================================================================
!
MODULE ACCURACY_REAL
!
      INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,  37)   ! single precision
      INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)  ! double precision
      INTEGER, PARAMETER :: QP = SELECTED_REAL_KIND(33, 4931) ! quadruple precision
!
      INTEGER, PARAMETER :: WP = DP                           ! selected value for code
!
END MODULE ACCURACY_REAL
!
!=======================================================================
!
MODULE ACCURACY_INTEGER
!
      INTEGER, PARAMETER :: I1 = SELECTED_INT_KIND(8)  ! default precision
      INTEGER, PARAMETER :: I2 = SELECTED_INT_KIND(16) ! single precision
      INTEGER, PARAMETER :: I4 = SELECTED_INT_KIND(32) ! double precision
      INTEGER, PARAMETER :: I8 = SELECTED_INT_KIND(64) ! quadruple precision
!
      INTEGER, PARAMETER :: IW = I1                    ! selected value for code
!
END MODULE ACCURACY_INTEGER
!
!=======================================================================
!
MODULE MINMAX_VALUES
!
      USE ACCURACY_REAL
      USE ACCURACY_INTEGER
!
      INTEGER (IW)            :: III
!
      REAL (WP)               :: XXX
!
      INTEGER (IW), PARAMETER ::  INT_MAX = HUGE(III)              ! maximal value of integer        
!
      REAL (WP), PARAMETER    ::  LN2     = 0.6931471805599453094172321214581765681D0  ! ln(2)              
 
      REAL (WP), PARAMETER    ::  MAX_2XP = MAXEXPONENT(XXX)       ! max value of y so that 2^x  is defined
      REAL (WP), PARAMETER    ::  MIN_2XP = MINEXPONENT(XXX)       ! max value of y so that 2^-x is defined
      REAL (WP), PARAMETER    ::  REL_MIN = TINY(XXX)              ! minimum value of real number
      REAL (WP), PARAMETER    ::  REL_MAX = HUGE(XXX)              ! maximum value of real number
      REAL (WP), PARAMETER    ::  EPS_MIN = EPSILON(XXX)           ! smallest value real such that x + epsilon /= x and x = 1
      REAL (WP), PARAMETER    ::  DGT_SIG = DIGITS(XXX)            ! number of significant digits
!
CONTAINS
!
!=======================================================================
!
      SUBROUTINE MINMAX_EXP(MAX_EXP,MIN_EXP)
!
!  This module computes the maximal and minimal exponent 
!    so that e^x is defined
!
      IMPLICIT NONE 
!
      REAL (WP), INTENT(OUT)  ::  MAX_EXP,MIN_EXP
!
      REAL (WP), PARAMETER    ::  LN2     = 0.6931471805599453094172321214581765681D0  ! ln(2)              
!
      MAX_EXP = INT(MAXEXPONENT(XXX) * LN2)         ! max value of y so that e^x  is defined
      MIN_EXP = INT(MINEXPONENT(XXX) * LN2)         ! max value of y so that e^-x is defined
!
      END SUBROUTINE MINMAX_EXP
!
END MODULE MINMAX_VALUES
!
!=======================================================================
!
MODULE MACHINE_ACCURACY 
!
!  This module provides the AMOS legacy routines for machine accuracy:
!
!                 * FUNCTION D1MACH(I)  --> double precision reals
!
!                 * FUNCTION I1MACH(I)  --> integers
!
!                 * FUNCTION R1MACH(I)  --> single precision reals
!
!
      USE ACCURACY_REAL
      USE ACCURACY_INTEGER
!
CONTAINS
!
!=======================================================================
!
      FUNCTION D1MACH(I)
!
      IMPLICIT NONE
!
      INTEGER (IW)          ::  I
!
      REAL (WP)             ::  D1MACH
      REAL (WP)             ::  B,X
!
!***BEGIN PROLOGUE  D1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (D1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH(3) = B**(-T), the smallest relative spacing.
!   D1MACH(4) = B**(1-T), the largest relative spacing.
!   D1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EHG)      
!***END PROLOGUE  D1MACH
!      
      X = 1.0E0_WP
      B = RADIX(X)
!
      SELECT CASE (I)
        CASE (1)
          D1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)               ! the largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
        CASE (4)
          D1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          WRITE (*,10)
          STOP
      END SELECT
!
!  Formats:
!
  10  FORMAT ('1ERROR    1 in D1MACH - I out of bounds')
!
      END FUNCTION D1MACH
!
!=======================================================================
!
      FUNCTION I1MACH(I)
!
      IMPLICIT NONE
!
      INTEGER               ::  I,I1MACH
!
      REAL (SP)             ::  X
!
      REAL (WP)             ::  XX
!
!***BEGIN PROLOGUE  I1MACH
!***PURPOSE  Return integer machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      INTEGER (I1MACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!        K = I1MACH(I)
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit.
!     I1MACH( 4) = the standard error message unit.
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0 .LE. X(I) .LT. B for I=1,...,T,
!                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   960411  Modified for Fortran 90 (BE after suggestions by EHG).   
!   980727  Modified value of I1MACH(6) (BE after suggestion by EHG).   
!***END PROLOGUE  I1MACH
!
      X  = 1.0      
      XX = 1.0E0_WP

      SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          I1MACH = 0 ! Punch unit is no longer used
        CASE (4)
          I1MACH = 0 ! Error message unit
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          I1MACH = 4            ! Characters per integer is hopefully no
                                ! longer used. 
                                ! If it is used it has to be set manually.
                                ! The value 4 is correct on IEEE-machines.
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X)
        CASE (11)
          I1MACH = DIGITS(X)
        CASE (12)
          I1MACH = MINEXPONENT(X)
        CASE (13)
          I1MACH = MAXEXPONENT(X)
        CASE (14)
          I1MACH = DIGITS(XX)
        CASE (15)
          I1MACH = MINEXPONENT(XX)
        CASE (16)
          I1MACH = MAXEXPONENT(XX) 
        CASE DEFAULT
          WRITE (*,10)
          STOP
        END SELECT
!
!  Formats
!
  10  FORMAT ('Fatal in I1MACH - I out of bounds')
!
      END FUNCTION I1MACH
!
!=======================================================================
!
      FUNCTION R1MACH (I)
!
      IMPLICIT NONE
!
      INTEGER               ::  I
!
      REAL (SP)             ::  B,X,R1MACH
!
!***BEGIN PROLOGUE  R1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   R1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = R1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   R1MACH(3) = B**(-T), the smallest relative spacing.
!   R1MACH(4) = B**(1-T), the largest relative spacing.
!   R1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EG)      
!***END PROLOGUE  R1MACH
!      
      X = 1.0
      B = RADIX(X)
!
      SELECT CASE (I)
        CASE (1)
          R1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
        CASE (2)
          R1MACH = HUGE(X)               ! the largest magnitude.
        CASE (3)
          R1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
        CASE (4)
          R1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
        CASE (5)
          R1MACH = LOG10(B)
        CASE DEFAULT
          WRITE (*,10)
          STOP
      END SELECT
!
!  Formats:
!
  10  FORMAT ('1ERROR    1 IN R1MACH - I out of bounds')
!
      END FUNCTION R1MACH
!
END MODULE MACHINE_ACCURACY
