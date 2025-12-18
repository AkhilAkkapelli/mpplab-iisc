PROGRAM test_multi_precision_float
  USE multi_precision_float_mod
  USE multi_precision_integer_mod, ONLY: mpi_t, new_mpi_from_integer, new_mpi_from_coeffs, new_mpi_from_string, &
                                      MULTI_PRECISION_BASE, mpi_to_string, mpi_is_zero, mpi_shift_bits_left, &
                                      OPERATOR(==), mpi_power, OPERATOR(+), OPERATOR(*), new_mpi_from_string
  IMPLICIT NONE

  WRITE(*, '(/A)') "Running tests for MultiPrecisionFloat module..."
  WRITE(*, '(/A)') "=============================================="

  CALL test_constructors()
  CALL test_arithmetic()
  CALL test_large_arithmetic()

  WRITE(*, '(/A)') "=============================================="
  WRITE(*, '(A/)') "All MultiPrecisionFloat tests completed."

CONTAINS

  ! A simple assertion utility to make tests cleaner.
  SUBROUTINE assert(condition, message)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(LEN=*), INTENT(IN) :: message
    IF (condition) THEN
      WRITE(*, '(A, A)') "  [PASS] ", message
    ELSE
      WRITE(*, '(A, A)') "  [FAIL] ", message
      STOP "Test assertion failed."
    END IF
  END SUBROUTINE assert

  ! Compares the numerical value of two HPFs by aligning their exponents.
  ! This is different from the built-in `==` which does a structural comparison.
  LOGICAL FUNCTION mpf_value_equal(a, b)
    TYPE(mpf), INTENT(IN) :: a, b
    INTEGER :: common_exp, diff_a, diff_b
    TYPE(mpi_t) :: scaled_a, scaled_b

    IF (mpi_is_zero(a%mantissa) .AND. mpi_is_zero(b%mantissa)) THEN
        mpf_value_equal = .TRUE.
        RETURN
    END IF

    ! To correctly compare values, we must align their exponents first.
    common_exp = MIN(a%exponent, b%exponent)
    diff_a = a%exponent - common_exp
    diff_b = b%exponent - common_exp

    scaled_a = mpi_shift_bits_left(a%mantissa, diff_a)
    scaled_b = mpi_shift_bits_left(b%mantissa, diff_b)

    mpf_value_equal = (scaled_a == scaled_b)
  END FUNCTION mpf_value_equal

  ! Tests constructors and the normalization process.
  SUBROUTINE test_constructors()
    TYPE(mpf) :: a, b, c, mpf1, mpf2, mpf_zero, mpf_pos, mpf_neg, mpf_copy
    TYPE(mpi_t)                :: five, mpi_unnormalized, mpi_val, large_mpi
    CHARACTER(LEN=100)         :: large_pos_str, large_neg_str, str_out

    WRITE(*, '(A)') "Testing constructors and I/O..."

    ! Test 1: Create -5 * 2^-100 and convert back to string
    CALL new_mpi_from_integer(-5_8, five)
    a = new_mpf_from_mpi_exp(five, -100)
    print*, "  Created MPF for -5 * 2^-100", a
    str_out = mpf_to_string(a)
    print*, "  -5 * 2^-100 as string: ", TRIM(str_out)
    CALL assert(TRIM(str_out) == &
    "-3.944304526105059027058642826413931148366e-30", "mpf_to_string: Small negative number")

    ! Test 2: Create from string "0.0048828125" (which is 5 * 2^-10)
    b = new_mpf_from_string("0.0048828125") ! This is 5 / 1024 = 5 * 2^-10
    print*, "  Created MPF from string '0.0048828125'", b
    CALL new_mpi_from_integer(5_8, five)
    CALL assert(b%mantissa == five .AND. b%exponent == -10, "new_mpf_from_string: Exact power of 2 fraction")

    ! Test construction from a large positive float string
    large_pos_str = "12345678901234567890.0"
    mpf1 = new_mpf_from_string(large_pos_str)
    CALL new_mpi_from_string("12345678901234567890", large_mpi)
    CALL assert(mpf_value_equal(mpf1, new_mpf_from_mpi_exp(large_mpi, 0)), &
                "new_mpf_from_string: Large integer as float string")

    ! Test round-trip for a large number
    str_out = mpf_to_string(mpf1)
    print*, "  Large integer as string: ", TRIM(str_out)
    CALL assert(TRIM(str_out) == "1.234567890123456789000000000000000000000e+19", "mpf_to_string: Large integer")
    
    ! Test zero variations
    mpf_zero = new_mpf_from_string("0.0")
    CALL assert(mpi_is_zero(mpf_zero%mantissa) .AND. mpf_zero%exponent == 0, "new_mpf_from_string: '0.0'")
    mpf_zero = new_mpf_from_string("-0")
    CALL assert(mpi_is_zero(mpf_zero%mantissa) .AND. mpf_zero%exponent == 0, "new_mpf_from_string: '-0'")

    ! Test simple negative number
    mpf_neg = new_mpf_from_string("-1.5")
    CALL new_mpi_from_integer(3_8, mpi_val)
    b = new_mpf_from_mpi_exp(mpi_val, -1)
    CALL assert(mpf_value_equal(mpf_neg, -b), "new_mpf_from_string: Negative float '-1.5'")

    ! Test scientific notation (positive exponent)
    a = new_mpf_from_string("1.25e3") ! 1250
    CALL new_mpi_from_integer(1250_8, mpi_val)
    b = new_mpf_from_mpi_exp(mpi_val, 0)
    CALL assert(mpf_value_equal(a, b), "new_mpf_from_string: Scientific notation '1.25e3'")

    ! Test scientific notation (negative exponent)
    a = new_mpf_from_string("-1.25e-2") ! -0.0125
    CALL new_mpi_from_integer(125_8, mpi_val)
    b = new_mpf_from_mpi_exp(mpi_val, -10) ! -125 * 2^-10 = -0.0122... No, -1.25e-2 = -125/10000
    c = new_mpf_from_string("-0.0125")
    CALL assert(mpf_value_equal(a, c), "new_mpf_from_string: Scientific notation '-1.25e-2'")

    ! Test round-trip for a number with many decimal places
    str_out = "3.1415926535897932384626433832795"
    a = new_mpf_from_string(TRIM(str_out))
    b = new_mpf_from_string(mpf_to_string(a))
    CALL assert(mpf_value_equal(a, b), "I/O round-trip: Multi-precision Pi")

  END SUBROUTINE test_constructors

  ! SUBROUTINE test_string_symmetry()
  !   TYPE(multi_precision_float) :: mpf_pos, mpf_neg

  !   WRITE(*, '(A)') "Testing construction symmetry..."

  !   ! Create HPF for "1.5" and "-1.5"
  !   mpf_pos = new_mpf_from_string("1.5")
  !   mpf_neg = new_mpf_from_string("-1.5")

  !   print*, "  MPF for 1.5: ", mpf_pos
  !   print*, "  MPF for -1.5: ", mpf_neg
  !   ! Assert that the mpf from "-1.5" is numerically equal to the negation of the mpf from "1.5"
  !   CALL assert(mpf_value_equal(mpf_neg, -mpf_pos), "Construction symmetry for 1.5 and -1.5")

  ! END SUBROUTINE test_string_symmetry

  SUBROUTINE test_arithmetic()
    TYPE(mpf) :: a, b, c, res
    TYPE(mpi_t) :: mpi_val

    WRITE(*, '(A)') "Testing arithmetic operators..."

    ! --- Addition Tests ---
    a = new_mpf_from_string("1.5")   ! 3 * 2^-1
    b = new_mpf_from_string("2.75")  ! 11 * 2^-2
    c = new_mpf_from_string("4.25")  ! 17 * 2^-2
    res = a + b
    CALL assert(mpf_value_equal(res, c), "Addition: 1.5 + 2.75 = 4.25")

    a = new_mpf_from_string("10.5")
    b = new_mpf_from_string("-2.5")
    c = new_mpf_from_string("8.0")
    res = a + b
    print*, "Addition: 10.5 + (-2.5) = ", res
    PRINT*, c
    CALL assert(mpf_value_equal(res, c), "Addition: 10.5 + (-2.5) = 8.0")

    a = new_mpf_from_string("1.23")
    b = new_mpf_from_string("-1.23")
    res = a + b
    CALL assert(mpi_is_zero(res%mantissa), "Addition: 1.23 + (-1.23) = 0.0")

    ! --- Subtraction Tests ---
    a = new_mpf_from_string("4.25")
    b = new_mpf_from_string("1.5")
    c = new_mpf_from_string("2.75")
    res = a - b
    CALL assert(mpf_value_equal(res, c), "Subtraction: 4.25 - 1.5 = 2.75")

    res = b - a
    CALL assert(mpf_value_equal(res, -c), "Subtraction: 1.5 - 4.25 = -2.75")

    ! --- Multiplication Tests ---
    a = new_mpf_from_string("1.5")   ! 3 * 2^-1
    b = new_mpf_from_string("2.0")   ! 1 * 2^1
    c = new_mpf_from_string("3.0")
    res = a * b
    CALL assert(mpf_value_equal(res, c), "Multiplication: 1.5 * 2.0 = 3.0")

    a = new_mpf_from_string("-1.25") ! -5 * 2^-2
    b = new_mpf_from_string("1.25")  !  5 * 2^-2
    c = new_mpf_from_string("-1.5625")! -25 * 2^-4
    res = a * b
    CALL assert(mpf_value_equal(res, c), "Multiplication: -1.25 * 1.25 = -1.5625")

    ! --- Division Tests ---
    a = new_mpf_from_string("3.0")
    b = new_mpf_from_string("2.0")
    c = new_mpf_from_string("1.5")
    res = a / b
    print*, "res = ", res
    print*, "a: ", a
    print*, "b: ", b
    print*, "c: ", c
    CALL assert(mpf_value_equal(res, c), "Division: 3.0 / 2.0 = 1.5")

    a = new_mpf_from_string("-1.5625")
    b = new_mpf_from_string("1.25")
    c = new_mpf_from_string("-1.25")
    res = a / b
    CALL assert(mpf_value_equal(res, c), "Division: -1.5625 / 1.25 = -1.25")

  END SUBROUTINE test_arithmetic

  SUBROUTINE test_large_arithmetic()
    TYPE(mpf) :: a, b, c, res

    WRITE(*, '(A)') "Testing arithmetic with large numbers (multi-coefficient)..."

    ! --- Addition Test ---
    ! Use numbers large enough to require > 1 coefficient (10^10 > 2^32)
    a = new_mpf_from_string("1.0e10")
    b = new_mpf_from_string("2.0e10")
    c = new_mpf_from_string("3.0e10")
    res = a + b
    CALL assert(mpf_value_equal(res, c), "Large Addition: 1e10 + 2e10 = 3e10")

    ! --- Subtraction Test ---
    a = new_mpf_from_string("5.0e12")
    b = new_mpf_from_string("2.0e12")
    c = new_mpf_from_string("3.0e12")
    res = a - b
    CALL assert(mpf_value_equal(res, c), "Large Subtraction: 5e12 - 2e12 = 3e12")

    ! --- Multiplication Test ---
    ! Product will require multiple coefficients
    a = new_mpf_from_string("1234567.0")
    b = new_mpf_from_string("7654321.0")
    c = new_mpf_from_string("9449772114007.0")
    res = a * b
    print*, "a: ", a
    print*, "b: ", b
    print*, "c: ", c
    print*, "Multiplication Result: ", res
    CALL assert(mpf_value_equal(res, c), "Large Multiplication: 1234567 * 7654321")

    ! --- Division Test ---
    ! Use large numbers for division
    a = new_mpf_from_string("9.449772114007e18") ! A very large number
    b = new_mpf_from_string("1234567.0")
    c = new_mpf_from_string("7.654321e12")
    res = a / b
    CALL assert(mpf_value_equal(res, c), "Large Division: 9.44e18 / 1.23e6")

  END SUBROUTINE test_large_arithmetic

END PROGRAM test_multi_precision_float