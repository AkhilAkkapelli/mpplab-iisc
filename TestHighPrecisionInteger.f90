PROGRAM test_multi_precision_integer
  USE multi_precision_integer_mod
  IMPLICIT NONE

  INTEGER(KIND=8), PARAMETER :: MASK32 = INT(Z'FFFFFFFF', KIND=8)

  CALL test_normalization()
  CALL test_new_mpi_from_coeffs()
  CALL test_new_mpi_from_integer()
  CALL test_mpi_to_integer()
  CALL test_mpi_scalar_arithmetic()
  ! CALL test_string_conversions_edge_cases()
  ! CALL test_arithmetic_operators()
  ! CALL test_large_string_conversion()
CONTAINS


  SUBROUTINE check(description, condition)
    CHARACTER(LEN=*), INTENT(IN) :: description
    LOGICAL, INTENT(IN)          :: condition
    IF (condition) THEN
      WRITE(*, '(A, A)') "[PASS] ", description
    ELSE
      WRITE(*, '(A, A)') "[FAIL] ", description
    END IF
  END SUBROUTINE check

  SUBROUTINE check_string(description, actual, expected)
    CHARACTER(LEN=*), INTENT(IN) :: description, actual, expected
    CALL check(description, TRIM(actual) == TRIM(expected))
    IF (TRIM(actual) /= TRIM(expected)) THEN
        WRITE(*, '(2A)') "       Expected: ", TRIM(expected)
        WRITE(*, '(2A)') "       Actual:   ", TRIM(actual)
    END IF
  END SUBROUTINE check_string

  SUBROUTINE check_mpi(description, mpi_val, expected_coeffs)
    CHARACTER(LEN=*), INTENT(IN)             :: description
    TYPE(mpi), INTENT(IN)                  :: mpi_val
    INTEGER(KIND=8), DIMENSION(:), INTENT(IN) :: expected_coeffs
    LOGICAL                                  :: is_ok
    CHARACTER(LEN=:), ALLOCATABLE            :: actual_str, expected_str
    CHARACTER(LEN=22)                        :: temp_str
    INTEGER                                  :: i

    is_ok = (ALL(mpi_val%coeffs(1:SIZE(expected_coeffs)) == expected_coeffs)) .AND. &
            (ALL(mpi_val%coeffs(SIZE(expected_coeffs)+1:) == 0_8))

    CALL check(description, is_ok)

    IF (.NOT. is_ok) THEN
        ! Build the 'actual' string
        actual_str = "["
        DO i = 1, COEFFS_LIMIT
            WRITE(temp_str, '(I0)') mpi_val%coeffs(i)
            actual_str = TRIM(actual_str) // TRIM(ADJUSTL(temp_str))
            IF (i < COEFFS_LIMIT) actual_str = TRIM(actual_str) // ", "
        END DO
        actual_str = TRIM(actual_str) // "]"
        WRITE(*, '(A, A)') "       Actual coeffs:   ", TRIM(actual_str)

        ! Build the 'expected' string
        expected_str = mpi_to_string(new_mpi_from_coeffs(expected_coeffs))
        WRITE(*, '(A, A, A, A)') "       Expected value:  ", TRIM(expected_str), " (coeffs: ", TRIM(actual_str), ")"
    END IF
  END SUBROUTINE check_mpi

SUBROUTINE test_normalization()
    TYPE(mpi) :: mpi_val

    PRINT *, ""
    PRINT *, "--- Testing Normalization ---"

    ! Test 1: All zeros
    mpi_val = new_mpi_from_integer(0_8)
    CALL normalize_mpi(mpi_val)
    CALL check_mpi("Normalization: All zeros", mpi_val, [0_8, 0_8, 0_8, 0_8])

    ! Test 2: No carry, already normalized
    mpi_val = new_mpi_from_coeffs([1_8, 2_8, 3_8, 0_8])
    CALL check_mpi("Normalization: No carry", mpi_val, [1_8, 2_8, 3_8, 0_8])

    ! Test 3: Simple carry
    mpi_val%coeffs = [MULTI_PRECISION_BASE + 5_8, 1_8, 0_8, 0_8]
    CALL normalize_mpi(mpi_val)
    CALL check_mpi("Normalization: Simple positive carry", mpi_val, [5_8, 2_8, 0_8, 0_8])

    ! Test 4: Multiple carries 
    mpi_val%coeffs = [MULTI_PRECISION_BASE + 1_8, MULTI_PRECISION_BASE - 1_8, 0_8, 0_8]
    CALL normalize_mpi(mpi_val)
    CALL check_mpi("Normalization: Cascading positive carry", mpi_val, [1_8, 0_8, 1_8, 0_8])

    ! Test 5: Carry out of the last coefficient
    mpi_val%coeffs = [0_8, 0_8, 0_8, MULTI_PRECISION_BASE]
    CALL normalize_mpi(mpi_val)
    CALL check_mpi("Normalization: Carry out of last coeff (overflow)", mpi_val, [0_8, 0_8, 0_8, 0_8])

    ! Test 6: Simple negative number
    mpi_val = new_mpi_from_coeffs([-1_8, -2_8, 0_8, 0_8])
    CALL check_mpi("Normalization: Simple negative, no borrow", mpi_val, [-1_8, -2_8, 0_8, 0_8])

    ! Test 7: Normalizing a negative value
    mpi_val%coeffs = [-(MULTI_PRECISION_BASE + 5_8), -(MULTI_PRECISION_BASE + 2_8), 0_8, 0_8]
    CALL normalize_mpi(mpi_val)
    CALL check_mpi("Normalization: Large negative, multiple borrows", mpi_val, [-5_8, -3_8, -1_8, 0_8])

    PRINT *, ""

END SUBROUTINE test_normalization

SUBROUTINE test_new_mpi_from_coeffs()
    TYPE(mpi) :: mpi_val

    PRINT *, ""
    PRINT *, "--- Testing new_mpi_from_coeffs ---"

    ! Test 1: Simple positive coefficients
    mpi_val = new_mpi_from_coeffs([123_8, 456_8])
    CALL check_mpi("From Coeffs: Simple positive", mpi_val, [123_8, 456_8, 0_8, 0_8])

    ! Test 2: Simple negative coefficients (note: normalize_mpi handles borrows, not sign uniformity)
    mpi_val = new_mpi_from_coeffs([-123_8, -456_8])
    CALL check_mpi("From Coeffs: Simple negative", mpi_val, [-123_8, -456_8, 0_8, 0_8])

    ! Test 3: Coefficients requiring normalization (carry)
    mpi_val = new_mpi_from_coeffs([MULTI_PRECISION_BASE + 10_8, 1_8])
    CALL check_mpi("From Coeffs: With positive normalization (carry)", mpi_val, [10_8, 2_8, 0_8, 0_8])

    ! Test 4: Empty coefficient array
    mpi_val = new_mpi_from_coeffs(coeffs=[INTEGER(KIND=8) ::])
    CALL check_mpi("From Coeffs: Empty input", mpi_val, [0_8, 0_8, 0_8, 0_8])

    ! Test 5: Single coefficient
    mpi_val = new_mpi_from_coeffs([999_8])
    CALL check_mpi("From Coeffs: Single coefficient", mpi_val, [999_8, 0_8, 0_8, 0_8])

    ! Test 6: Full coefficient array
    mpi_val = new_mpi_from_coeffs([1_8, 2_8, 3_8, 4_8])
    CALL check_mpi("From Coeffs: Full array", mpi_val, [1_8, 2_8, 3_8, 4_8])

    ! Test 7: Coefficients that normalize to zero
    mpi_val = new_mpi_from_coeffs([-MULTI_PRECISION_BASE, 1_8])
    CALL check_mpi("From Coeffs: Normalization to zero", mpi_val, [0_8, 0_8, 0_8, 0_8])

    PRINT *, ""

END SUBROUTINE test_new_mpi_from_coeffs

SUBROUTINE test_new_mpi_from_integer()
    TYPE(mpi) :: mpi_val
    INTEGER(KIND=8) :: int_val

    PRINT *, ""
    PRINT *, "--- Testing new_mpi_from_integer ---"

    ! Test 1: Zero
    mpi_val = new_mpi_from_integer(0_8)
    CALL check_mpi("From Integer: 0", mpi_val, [0_8, 0_8, 0_8, 0_8])

    ! Test 2: Small positive
    mpi_val = new_mpi_from_integer(12345_8)
    CALL check_mpi("From Integer: Small positive", mpi_val, [12345_8, 0_8, 0_8, 0_8])

    ! Test 3: Small negative
    mpi_val = new_mpi_from_integer(-54321_8)
    CALL check_mpi("From Integer: Small negative", mpi_val, [-54321_8, 0_8, 0_8, 0_8])

    ! Test 4: Positive, crossing one coefficient boundary
    int_val = MULTI_PRECISION_BASE + 100_8
    mpi_val = new_mpi_from_integer(int_val)
    CALL check_mpi("From Integer: Large positive", mpi_val, [100_8, 1_8, 0_8, 0_8])

    ! Test 5: Negative, crossing one coefficient boundary
    int_val = -(MULTI_PRECISION_BASE + 100_8)
    mpi_val = new_mpi_from_integer(int_val)
    CALL check_mpi("From Integer: Large negative", mpi_val, [-100_8, -1_8, 0_8, 0_8])

    ! Test 6: Maximum 64-bit integer
    int_val = HUGE(0_8)
    mpi_val = new_mpi_from_integer(int_val)
    CALL check_mpi("From Integer: HUGE(0_8)", mpi_val, [MASK32, ISHFT(int_val, -32), 0_8, 0_8])

    ! Test 7: Minimum 64-bit integer
    int_val = -HUGE(0_8) - 1_8
    mpi_val = new_mpi_from_integer(int_val)
    CALL check_mpi("From Integer: Most negative", mpi_val, [0_8, -ISHFT(MULTI_PRECISION_BASE, -1), 0_8, 0_8])

    ! Test 8: Round trip check for a large value
    int_val = 123456789012345_8
    mpi_val = new_mpi_from_integer(int_val)
    CALL check("From/To Integer: Round trip", mpi_to_integer(mpi_val) == int_val)

    PRINT *, ""

END SUBROUTINE test_new_mpi_from_integer

SUBROUTINE test_mpi_to_integer()
    TYPE(mpi) :: mpi_val
    INTEGER(KIND=8) :: int_val, expected_val

    PRINT *, ""
    PRINT *, "--- Testing mpi_to_integer ---"

    ! Test 1: Zero
    mpi_val = new_mpi_from_coeffs([0_8, 0_8, 0_8, 0_8])
    int_val = mpi_to_integer(mpi_val)
    CALL check("To Integer: 0", int_val == 0_8)

    ! Test 2: Small positive
    mpi_val = new_mpi_from_coeffs([12345_8, 0_8, 0_8, 0_8])
    int_val = mpi_to_integer(mpi_val)
    CALL check("To Integer: Small positive", int_val == 12345_8)

    ! Test 3: Small negative
    mpi_val = new_mpi_from_coeffs([-54321_8, 0_8, 0_8, 0_8])
    int_val = mpi_to_integer(mpi_val)
    CALL check("To Integer: Small negative", int_val == -54321_8)

    ! Test 4: Large positive (two coefficients)
    expected_val = MULTI_PRECISION_BASE * 2_8 + 123_8
    mpi_val = new_mpi_from_coeffs([123_8, 2_8, 0_8, 0_8])
    int_val = mpi_to_integer(mpi_val)
    CALL check("To Integer: Large positive", int_val == expected_val)

    ! Test 5: Large negative (two coefficients)
    expected_val = -(MULTI_PRECISION_BASE * 3_8 + 456_8)
    mpi_val = new_mpi_from_coeffs([-456_8, -3_8, 0_8, 0_8])
    int_val = mpi_to_integer(mpi_val)
    CALL check("To Integer: Large negative", int_val == expected_val)

    ! Test 6: Maximum 64-bit integer
    mpi_val = new_mpi_from_integer(HUGE(0_8))
    int_val = mpi_to_integer(mpi_val)
    CALL check("To Integer: HUGE(0_8)", int_val == HUGE(0_8))

    ! Test 7: Minimum 64-bit integer
    mpi_val = new_mpi_from_integer(-HUGE(0_8) - 1_8)
    int_val = mpi_to_integer(mpi_val)
    CALL check("To Integer: Most negative", int_val == -HUGE(0_8)-1_8)

END SUBROUTINE test_mpi_to_integer

SUBROUTINE test_mpi_scalar_arithmetic()
    TYPE(mpi) :: mpi_val
    INTEGER(KIND=8) :: scalar, remainder

    PRINT *, ""
    PRINT *, "--- Testing MPI Scalar Arithmetic (Multiply/Divide) ---"

    ! --- Multiplication Tests ---
    PRINT *, "  Testing Multiplication..."
    mpi_val = new_mpi_from_integer(12345_8)
    CALL mpi_multiply_by_scalar(mpi_val, 0_8)
    CALL check_mpi("Multiply by 0", mpi_val, [0_8, 0_8, 0_8, 0_8])

    mpi_val = new_mpi_from_integer(12345_8)
    CALL mpi_multiply_by_scalar(mpi_val, 1_8)
    CALL check_mpi("Multiply by 1", mpi_val, [12345_8, 0_8, 0_8, 0_8])

    mpi_val = new_mpi_from_integer(100_8)
    CALL mpi_multiply_by_scalar(mpi_val, 5_8)
    CALL check_mpi("Multiply positive by positive", mpi_val, [500_8, 0_8, 0_8, 0_8])

    mpi_val = new_mpi_from_coeffs([MULTI_PRECISION_BASE - 1_8, 1_8]) ! Represents 2*B - 1
    CALL mpi_multiply_by_scalar(mpi_val, 2_8)
    CALL check_mpi("Multiply with carry", mpi_val, [MULTI_PRECISION_BASE - 2_8, 3_8, 0_8, 0_8])

    mpi_val = new_mpi_from_integer(-100_8)
    CALL mpi_multiply_by_scalar(mpi_val, 5_8)
    CALL check_mpi("Multiply negative by positive", mpi_val, [-500_8, 0_8, 0_8, 0_8])

    ! --- Division Tests ---
    PRINT *, "  Testing Division..."
    mpi_val = new_mpi_from_integer(100_8)
    remainder = mpi_div_by_scalar(mpi_val, 10_8)
    CALL check_mpi("Div by scalar: Quotient for 100/10", mpi_val, [10_8, 0_8, 0_8, 0_8])
    CALL check("Div by scalar: Remainder for 100/10", remainder == 0_8)

    mpi_val = new_mpi_from_integer(105_8)
    remainder = mpi_div_by_scalar(mpi_val, 10_8)
    CALL check_mpi("Div by scalar: Quotient for 105/10", mpi_val, [10_8, 0_8, 0_8, 0_8])
    CALL check("Div by scalar: Remainder for 105/10", remainder == 5_8)

    mpi_val = new_mpi_from_coeffs([10_8, 5_8]) ! Represents 5*B + 10
    remainder = mpi_div_by_scalar(mpi_val, 3_8)
    CALL check("Div by scalar: Remainder for (5B+10)/3", remainder == 1_8)

    mpi_val = new_mpi_from_integer(-105_8)
    remainder = mpi_div_by_scalar(mpi_val, 10_8)
    CALL check_mpi("Div by scalar: Quotient for -105/10", mpi_val, [-10_8, 0_8, 0_8, 0_8])
    CALL check("Div by scalar: Remainder for -105/10", remainder == -5_8)

    ! --- Round-trip Test ---
    PRINT *, "  Testing Multiply/Divide Round-trip..."
    mpi_val = new_mpi_from_integer(1234567_8)
    scalar = 987_8
    CALL mpi_multiply_by_scalar(mpi_val, scalar)
    remainder = mpi_div_by_scalar(mpi_val, scalar)
    CALL check("Round-trip: Remainder is zero", remainder == 0_8)
    CALL check("Round-trip: Original value is recovered", mpi_val == new_mpi_from_integer(1234567_8))

END SUBROUTINE test_mpi_scalar_arithmetic

SUBROUTINE test_string_conversions_edge_cases()
    TYPE(mpi) :: mpi
    CHARACTER(LEN=:), ALLOCATABLE :: str_val

    PRINT *, ""
    PRINT *, "--- Testing String Conversion Edge Cases ---"

    ! Test leading/trailing spaces
    CALL new_mpi_from_string("  123  ", mpi)
    str_val = mpi_to_string(mpi)
    CALL check_string("String Edge Case: Leading/trailing spaces", str_val, "123")

    ! Test explicit positive sign
    CALL new_mpi_from_string("+456", mpi)
    str_val = mpi_to_string(mpi)
    CALL check_string("String Edge Case: Explicit positive sign", str_val, "456")

    ! Test leading zeros
    CALL new_mpi_from_string("000789", mpi)
    str_val = mpi_to_string(mpi)
    CALL check_string("String Edge Case: Leading zeros", str_val, "789")

    ! Test negative with leading zeros
    CALL new_mpi_from_string("-00789", mpi)
    str_val = mpi_to_string(mpi)
    CALL check_string("String Edge Case: Negative with leading zeros", str_val, "-789")

    ! Test string containing only zeros
    CALL new_mpi_from_string("000", mpi)
    str_val = mpi_to_string(mpi)
    CALL check_string("String Edge Case: Only zeros", str_val, "0")

    ! Test empty string
    CALL new_mpi_from_string("", mpi)
    str_val = mpi_to_string(mpi)
    CALL check_string("String Edge Case: Empty string", str_val, "0")

END SUBROUTINE test_string_conversions_edge_cases

SUBROUTINE test_arithmetic_operators()
    TYPE(mpi) :: a, b, c, expected_mpi

    PRINT *, ""
    PRINT *, "--- Testing Arithmetic Operators ---"

    ! --- Addition ---
    CALL new_mpi_from_string("100", a)
    CALL new_mpi_from_string("200", b)
    c = a + b
    CALL new_mpi_from_string("300", expected_mpi)
    CALL check_mpi("Addition: 100 + 200", c, expected_mpi%coeffs)

    CALL new_mpi_from_string("-100", a)
    CALL new_mpi_from_string("50", b)
    c = a + b
    CALL new_mpi_from_string("-50", expected_mpi)
    CALL check_mpi("Addition: -100 + 50", c, expected_mpi%coeffs)

    CALL new_mpi_from_string("10633823966279326983230400083533428776", a) 
    CALL new_mpi_from_string("2658455991569831745809959798659454619", b)
    c = a + b
    CALL new_mpi_from_string("13292279957849158729040359882192883395", expected_mpi)
    CALL check_mpi("Addition: Large with carry", c, expected_mpi%coeffs)

    ! --- Subtraction ---
    CALL new_mpi_from_string("300", a)
    CALL new_mpi_from_string("100", b)
    c = a - b
    CALL new_mpi_from_string("200", expected_mpi)
    CALL check_mpi("Subtraction: 300 - 100", c, expected_mpi%coeffs)

    CALL new_mpi_from_string("50", a)
    CALL new_mpi_from_string("100", b)
    c = a - b
    CALL new_mpi_from_string("-50", expected_mpi)
    CALL check_mpi("Subtraction: 50 - 100", c, expected_mpi%coeffs)

    ! --- Multiplication ---
    CALL new_mpi_from_string("-2305842971231290661", a)
    CALL new_mpi_from_string("4614031696526153371", b)
    c = a * b
    print*, "mult: ", mpi_to_string(c)
    CALL new_mpi_from_string("-10639232556473218309152790657965968231", expected_mpi)
    CALL check_mpi("Multiplication: 1000 * 2000", c, expected_mpi%coeffs)

END SUBROUTINE test_arithmetic_operators

SUBROUTINE test_large_string_conversion()
    TYPE(mpi) :: mpi
    CHARACTER(LEN=:), ALLOCATABLE :: str_val
    CHARACTER(LEN=100) :: expected_str

    PRINT *, ""
    PRINT *, "--- Testing Large String Conversion ---"

    expected_str = "-1329227995784915872903806277077091289" 

    CALL new_mpi_from_string(expected_str, mpi)
    CALL check_mpi("From String: ", mpi, [-2775761881_8, -4294967113_8, -4294967295_8, -16777215_8])
    str_val = mpi_to_string(mpi)
    CALL check_string("To String: ", str_val, expected_str)
END SUBROUTINE test_large_string_conversion

END PROGRAM test_multi_precision_integer