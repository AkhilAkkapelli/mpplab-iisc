PROGRAM test_mpi_subroutines
  USE multi_precision_integer_mod
  IMPLICIT NONE
  
  INTEGER :: test_count, pass_count, fail_count
  
  test_count = 0
  pass_count = 0
  fail_count = 0
  
  PRINT *, "============================================"
  PRINT *, "MPI Multi-Precision Integer Test Suite - Part 2"
  PRINT *, "Testing with COEFFS_LIMIT = 4"
  PRINT *, "BASE = 2^32 = ", MULTI_PRECISION_BASE
  PRINT *, "============================================"
  PRINT *, ""
  
  CALL test_mpi_to_integer(test_count, pass_count, fail_count)
  CALL test_mpi_multiply_by_scalar(test_count, pass_count, fail_count)
  CALL test_mpi_add(test_count, pass_count, fail_count)
  CALL test_mpi_subtract(test_count, pass_count, fail_count)
  CALL test_mpi_multiply(test_count, pass_count, fail_count)
  CALL test_mpi_div_rem(test_count, pass_count, fail_count)
  
  PRINT *, ""
  PRINT *, "============================================"
  PRINT *, "TEST SUMMARY"
  PRINT *, "============================================"
  PRINT *, "Total tests:  ", test_count
  PRINT *, "Passed:       ", pass_count
  PRINT *, "Failed:       ", fail_count
  PRINT *, "============================================"
  
CONTAINS

  SUBROUTINE test_mpi_to_integer(test_count, pass_count, fail_count)
    INTEGER, INTENT(INOUT) :: test_count, pass_count, fail_count
    TYPE(mpi) :: test_mpi
    INTEGER(KIND=8) :: result
    
    PRINT *, "Testing mpi_to_integer..."
    PRINT *, "----------------------------------------"
    
    ! Test 1: Zero
    test_count = test_count + 1
    test_mpi%coeffs = [0_8, 0_8, 0_8, 0_8]
    result = mpi_to_integer(test_mpi)
    IF (result == 0_8) THEN
      PRINT *, "Test 1 PASSED: Zero conversion"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 1 FAILED: Expected 0, got ", result
      fail_count = fail_count + 1
    END IF
    
    ! Test 2: Small positive value (fits in first coefficient)
    test_count = test_count + 1
    test_mpi%coeffs = [42_8, 0_8, 0_8, 0_8]
    result = mpi_to_integer(test_mpi)
    IF (result == 42_8) THEN
      PRINT *, "Test 2 PASSED: Small positive value"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 2 FAILED: Expected 42, got ", result
      fail_count = fail_count + 1
    END IF
    
    ! Test 3: Value using both first and second coefficients
    test_count = test_count + 1
    test_mpi%coeffs = [MULTI_PRECISION_BASE - 1_8, 1_8, 0_8, 0_8]
    result = mpi_to_integer(test_mpi)
    IF (result == MULTI_PRECISION_BASE + (MULTI_PRECISION_BASE - 1_8)) THEN
      PRINT *, "Test 3 PASSED: Two-coefficient value"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 3 FAILED: Two-coefficient value"
      fail_count = fail_count + 1
    END IF
    
    ! Test 4: Maximum 32-bit value in first coefficient
    test_count = test_count + 1
    test_mpi%coeffs = [INT(Z'FFFFFFFF', KIND=8), 0_8, 0_8, 0_8]
    result = mpi_to_integer(test_mpi)
    IF (result == INT(Z'FFFFFFFF', KIND=8)) THEN
      PRINT *, "Test 4 PASSED: Maximum 32-bit value"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 4 FAILED: Maximum 32-bit value"
      fail_count = fail_count + 1
    END IF
    
    ! Test 5: Power of 2 in second coefficient
    test_count = test_count + 1
    test_mpi%coeffs = [0_8, 1_8, 0_8, 0_8]
    result = mpi_to_integer(test_mpi)
    IF (result == MULTI_PRECISION_BASE) THEN
      PRINT *, "Test 5 PASSED: Power of 2 (BASE)"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 5 FAILED: Expected ", MULTI_PRECISION_BASE, ", got ", result
      fail_count = fail_count + 1
    END IF
    
    ! Test 6: Combined value from coefficients
    test_count = test_count + 1
    test_mpi%coeffs = [123456789_8, 987654321_8, 0_8, 0_8]
    result = mpi_to_integer(test_mpi)
    IF (result == 123456789_8 + ISHFT(987654321_8, 32)) THEN
      PRINT *, "Test 6 PASSED: Combined coefficient value"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 6 FAILED: Combined coefficient value"
      fail_count = fail_count + 1
    END IF
    
    ! Test 7: Second coefficient at maximum
    test_count = test_count + 1
    test_mpi%coeffs = [0_8, MULTI_PRECISION_BASE - 1_8, 0_8, 0_8]
    result = mpi_to_integer(test_mpi)
    IF (result == ISHFT(MULTI_PRECISION_BASE - 1_8, 32)) THEN
      PRINT *, "Test 7 PASSED: Second coefficient at maximum"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 7 FAILED: Second coefficient at maximum"
      fail_count = fail_count + 1
    END IF
    
    ! Test 8: Value with pattern 0x0000000100000000
    test_count = test_count + 1
    test_mpi%coeffs = [0_8, 1_8, 0_8, 0_8]
    result = mpi_to_integer(test_mpi)
    IF (result == 4294967296_8) THEN
      PRINT *, "Test 8 PASSED: 2^32 value"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 8 FAILED: Expected 4294967296, got ", result
      fail_count = fail_count + 1
    END IF
    
    PRINT *, ""
  END SUBROUTINE test_mpi_to_integer

  SUBROUTINE test_mpi_multiply_by_scalar(test_count, pass_count, fail_count)
    INTEGER, INTENT(INOUT) :: test_count, pass_count, fail_count
    TYPE(mpi) :: test_mpi
    
    PRINT *, "Testing mpi_multiply_by_scalar..."
    PRINT *, "----------------------------------------"
    
    ! Test 1: Multiply by zero
    test_count = test_count + 1
    test_mpi%coeffs = [100_8, 200_8, 300_8, 400_8]
    CALL mpi_multiply_by_scalar(test_mpi, 0_8)
    IF (ALL(test_mpi%coeffs == 0_8)) THEN
      PRINT *, "Test 1 PASSED: Multiply by zero"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 1 FAILED: Multiply by zero"
      fail_count = fail_count + 1
    END IF
    
    ! Test 2: Multiply by one (identity)
    test_count = test_count + 1
    test_mpi%coeffs = [123_8, 456_8, 789_8, 0_8]
    CALL mpi_multiply_by_scalar(test_mpi, 1_8)
    IF (test_mpi%coeffs(1) == 123_8 .AND. test_mpi%coeffs(2) == 456_8 &
        .AND. test_mpi%coeffs(3) == 789_8) THEN
      PRINT *, "Test 2 PASSED: Multiply by one"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 2 FAILED: Multiply by one"
      fail_count = fail_count + 1
    END IF
    
    ! Test 3: Simple multiplication without carry
    test_count = test_count + 1
    test_mpi%coeffs = [100_8, 0_8, 0_8, 0_8]
    CALL mpi_multiply_by_scalar(test_mpi, 2_8)
    IF (test_mpi%coeffs(1) == 200_8 .AND. test_mpi%coeffs(2) == 0_8) THEN
      PRINT *, "Test 3 PASSED: Simple multiplication (100 * 2)"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 3 FAILED: Expected [200, 0, 0, 0], got [", &
               test_mpi%coeffs(1), test_mpi%coeffs(2), test_mpi%coeffs(3), test_mpi%coeffs(4), "]"
      fail_count = fail_count + 1
    END IF
    
    ! Test 4: Multiplication causing carry to next position
    test_count = test_count + 1
    test_mpi%coeffs = [MULTI_PRECISION_BASE / 2_8, 0_8, 0_8, 0_8]
    CALL mpi_multiply_by_scalar(test_mpi, 3_8)
    IF (test_mpi%coeffs(1) == MULTI_PRECISION_BASE / 2_8 .AND. test_mpi%coeffs(2) == 1_8) THEN
      PRINT *, "Test 4 PASSED: Multiplication with carry"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 4 FAILED: Multiplication with carry"
      fail_count = fail_count + 1
    END IF
    
    ! Test 5: Multiply maximum single coefficient by 2
    test_count = test_count + 1
    test_mpi%coeffs = [MULTI_PRECISION_BASE - 1_8, 0_8, 0_8, 0_8]
    CALL mpi_multiply_by_scalar(test_mpi, 2_8)
    IF (test_mpi%coeffs(1) == MULTI_PRECISION_BASE - 2_8 .AND. test_mpi%coeffs(2) == 1_8) THEN
      PRINT *, "Test 5 PASSED: Maximum coefficient * 2"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 5 FAILED: Expected [", MULTI_PRECISION_BASE - 2_8, ", 1, 0, 0], got [", &
               test_mpi%coeffs(1), test_mpi%coeffs(2), test_mpi%coeffs(3), test_mpi%coeffs(4), "]"
      fail_count = fail_count + 1
    END IF
    
    ! Test 6: Multiply by large scalar
    test_count = test_count + 1
    test_mpi%coeffs = [1000_8, 0_8, 0_8, 0_8]
    CALL mpi_multiply_by_scalar(test_mpi, 1000000_8)
    IF (test_mpi%coeffs(1) == 1000000000_8 .AND. test_mpi%coeffs(2) == 0_8) THEN
      PRINT *, "Test 6 PASSED: Large scalar multiplication"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 6 FAILED: Large scalar multiplication"
      fail_count = fail_count + 1
    END IF
    
    ! Test 7: Cascading carries through multiple positions
    test_count = test_count + 1
    test_mpi%coeffs = [MULTI_PRECISION_BASE - 1_8, MULTI_PRECISION_BASE - 1_8, 0_8, 0_8]
    CALL mpi_multiply_by_scalar(test_mpi, 2_8)
    IF (test_mpi%coeffs(1) == MULTI_PRECISION_BASE - 2_8 .AND. &
        test_mpi%coeffs(2) == MULTI_PRECISION_BASE - 1_8 .AND. &
        test_mpi%coeffs(3) == 1_8) THEN
      PRINT *, "Test 7 PASSED: Cascading carries"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 7 FAILED: Cascading carries"
      fail_count = fail_count + 1
    END IF
    
    ! Test 8: Multiply across all four coefficients
    test_count = test_count + 1
    test_mpi%coeffs = [1_8, 2_8, 3_8, 4_8]
    CALL mpi_multiply_by_scalar(test_mpi, 10_8)
    IF (test_mpi%coeffs(1) == 10_8 .AND. test_mpi%coeffs(2) == 20_8 &
        .AND. test_mpi%coeffs(3) == 30_8 .AND. test_mpi%coeffs(4) == 40_8) THEN
      PRINT *, "Test 8 PASSED: Multiply all four coefficients"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 8 FAILED: Expected [10, 20, 30, 40], got [", &
               test_mpi%coeffs(1), test_mpi%coeffs(2), test_mpi%coeffs(3), test_mpi%coeffs(4), "]"
      fail_count = fail_count + 1
    END IF
    
    ! Test 9: Zero mpi multiplied by non-zero scalar
    test_count = test_count + 1
    test_mpi%coeffs = [0_8, 0_8, 0_8, 0_8]
    CALL mpi_multiply_by_scalar(test_mpi, 999999_8)
    IF (ALL(test_mpi%coeffs == 0_8)) THEN
      PRINT *, "Test 9 PASSED: Zero * large scalar"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 9 FAILED: Zero * large scalar"
      fail_count = fail_count + 1
    END IF
    
    ! Test 10: Multiply by power of 2 (bit shift equivalent)
    test_count = test_count + 1
    test_mpi%coeffs = [1_8, 0_8, 0_8, 0_8]
    CALL mpi_multiply_by_scalar(test_mpi, 256_8)
    IF (test_mpi%coeffs(1) == 256_8 .AND. test_mpi%coeffs(2) == 0_8) THEN
      PRINT *, "Test 10 PASSED: Multiply by 256 (2^8)"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 10 FAILED: Multiply by 256"
      fail_count = fail_count + 1
    END IF
    
    PRINT *, ""
  END SUBROUTINE test_mpi_multiply_by_scalar

  SUBROUTINE test_mpi_add(test_count, pass_count, fail_count)
    INTEGER, INTENT(INOUT) :: test_count, pass_count, fail_count
    TYPE(mpi) :: a, b, result_mpi
    
    PRINT *, "Testing mpi_add..."
    PRINT *, "----------------------------------------"
    
    ! Test 1: Add two zeros
    test_count = test_count + 1
    a = new_mpi_from_integer(0_8)
    b = new_mpi_from_integer(0_8)
    result_mpi = a + b
    IF (mpi_is_zero(result_mpi)) THEN
      PRINT *, "Test 1 PASSED: 0 + 0 = 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 1 FAILED: 0 + 0"
      fail_count = fail_count + 1
    END IF
    
    ! Test 2: Add zero to non-zero (identity)
    test_count = test_count + 1
    a = new_mpi_from_integer(123_8)
    b = new_mpi_from_integer(0_8)
    result_mpi = a + b
    IF (mpi_to_integer(result_mpi) == 123_8) THEN
      PRINT *, "Test 2 PASSED: 123 + 0 = 123"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 2 FAILED: Expected 123, got ", mpi_to_integer(result_mpi)
      fail_count = fail_count + 1
    END IF
    
    ! Test 3: Simple addition without carry
    test_count = test_count + 1
    a = new_mpi_from_integer(100_8)
    b = new_mpi_from_integer(50_8)
    result_mpi = a + b
    IF (mpi_to_integer(result_mpi) == 150_8) THEN
      PRINT *, "Test 3 PASSED: 100 + 50 = 150"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 3 FAILED: 100 + 50"
      fail_count = fail_count + 1
    END IF
    
    ! Test 4: Addition causing carry to next coefficient
    test_count = test_count + 1
    a%coeffs = [MULTI_PRECISION_BASE - 10_8, 0_8, 0_8, 0_8]
    b%coeffs = [20_8, 0_8, 0_8, 0_8]
    result_mpi = a + b
    IF (result_mpi%coeffs(1) == 10_8 .AND. result_mpi%coeffs(2) == 1_8) THEN
      PRINT *, "Test 4 PASSED: Addition with carry"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 4 FAILED: Expected [10, 1, 0, 0], got [", &
               result_mpi%coeffs(1), result_mpi%coeffs(2), result_mpi%coeffs(3), result_mpi%coeffs(4), "]"
      fail_count = fail_count + 1
    END IF
    
    ! Test 5: Cascading carries through multiple positions
    test_count = test_count + 1
    a%coeffs = [MULTI_PRECISION_BASE - 1_8, MULTI_PRECISION_BASE - 1_8, 0_8, 0_8]
    b%coeffs = [1_8, 0_8, 0_8, 0_8]
    result_mpi = a + b
    IF (result_mpi%coeffs(1) == 0_8 .AND. result_mpi%coeffs(2) == 0_8 &
        .AND. result_mpi%coeffs(3) == 1_8) THEN
      PRINT *, "Test 5 PASSED: Cascading carries"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 5 FAILED: Cascading carries"
      fail_count = fail_count + 1
    END IF
    
    ! Test 6: Add positive and negative (should subtract)
    test_count = test_count + 1
    a = new_mpi_from_integer(100_8)
    b = new_mpi_from_integer(-30_8)
    result_mpi = a + b
    IF (mpi_to_integer(result_mpi) == 70_8) THEN
      PRINT *, "Test 6 PASSED: 100 + (-30) = 70"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 6 FAILED: 100 + (-30)"
      fail_count = fail_count + 1
    END IF
    
    ! Test 7: Add two negative numbers
    test_count = test_count + 1
    a = new_mpi_from_integer(-50_8)
    b = new_mpi_from_integer(-75_8)
    result_mpi = a + b
    IF (mpi_to_integer(result_mpi) == -125_8) THEN
      PRINT *, "Test 7 PASSED: (-50) + (-75) = -125"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 7 FAILED: (-50) + (-75)"
      fail_count = fail_count + 1
    END IF
    
    ! Test 8: Add to result in zero (inverse)
    test_count = test_count + 1
    a = new_mpi_from_integer(999_8)
    b = new_mpi_from_integer(-999_8)
    result_mpi = a + b
    IF (mpi_is_zero(result_mpi)) THEN
      PRINT *, "Test 8 PASSED: 999 + (-999) = 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 8 FAILED: 999 + (-999) should be 0"
      fail_count = fail_count + 1
    END IF
    
    ! Test 9: Add large values spanning all coefficients
    test_count = test_count + 1
    a%coeffs = [100_8, 200_8, 300_8, 400_8]
    b%coeffs = [50_8, 100_8, 150_8, 200_8]
    result_mpi = a + b
    IF (result_mpi%coeffs(1) == 150_8 .AND. result_mpi%coeffs(2) == 300_8 &
        .AND. result_mpi%coeffs(3) == 450_8 .AND. result_mpi%coeffs(4) == 600_8) THEN
      PRINT *, "Test 9 PASSED: Add all four coefficients"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 9 FAILED: Add all four coefficients"
      fail_count = fail_count + 1
    END IF
    
    ! Test 10: Commutativity (a + b = b + a)
    test_count = test_count + 1
    a = new_mpi_from_integer(12345_8)
    b = new_mpi_from_integer(67890_8)
    IF ((a + b) == (b + a)) THEN
      PRINT *, "Test 10 PASSED: Addition commutativity"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 10 FAILED: Addition not commutative"
      fail_count = fail_count + 1
    END IF
    
    PRINT *, ""
  END SUBROUTINE test_mpi_add

  SUBROUTINE test_mpi_subtract(test_count, pass_count, fail_count)
    INTEGER, INTENT(INOUT) :: test_count, pass_count, fail_count
    TYPE(mpi) :: a, b, result_mpi
    
    PRINT *, "Testing mpi_subtract..."
    PRINT *, "----------------------------------------"
    
    ! Test 1: Subtract zero from zero
    test_count = test_count + 1
    a = new_mpi_from_integer(0_8)
    b = new_mpi_from_integer(0_8)
    result_mpi = a - b
    IF (mpi_is_zero(result_mpi)) THEN
      PRINT *, "Test 1 PASSED: 0 - 0 = 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 1 FAILED: 0 - 0"
      fail_count = fail_count + 1
    END IF
    
    ! Test 2: Subtract zero from non-zero
    test_count = test_count + 1
    a = new_mpi_from_integer(100_8)
    b = new_mpi_from_integer(0_8)
    result_mpi = a - b
    IF (mpi_to_integer(result_mpi) == 100_8) THEN
      PRINT *, "Test 2 PASSED: 100 - 0 = 100"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 2 FAILED: 100 - 0"
      fail_count = fail_count + 1
    END IF
    
    ! Test 3: Simple subtraction without borrow
    test_count = test_count + 1
    a = new_mpi_from_integer(150_8)
    b = new_mpi_from_integer(50_8)
    result_mpi = a - b
    IF (mpi_to_integer(result_mpi) == 100_8) THEN
      PRINT *, "Test 3 PASSED: 150 - 50 = 100"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 3 FAILED: 150 - 50"
      fail_count = fail_count + 1
    END IF
    
    ! Test 4: Subtract to get negative result
    test_count = test_count + 1
    a = new_mpi_from_integer(50_8)
    b = new_mpi_from_integer(100_8)
    result_mpi = a - b
    IF (mpi_to_integer(result_mpi) == -50_8) THEN
      PRINT *, "Test 4 PASSED: 50 - 100 = -50"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 4 FAILED: 50 - 100"
      fail_count = fail_count + 1
    END IF
    
    ! Test 5: Subtract to get zero (self-subtraction)
    test_count = test_count + 1
    a = new_mpi_from_integer(999_8)
    b = new_mpi_from_integer(999_8)
    result_mpi = a - b
    IF (mpi_is_zero(result_mpi)) THEN
      PRINT *, "Test 5 PASSED: 999 - 999 = 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 5 FAILED: 999 - 999"
      fail_count = fail_count + 1
    END IF
    
    ! Test 6: Subtract negative number (double negative becomes addition)
    test_count = test_count + 1
    a = new_mpi_from_integer(100_8)
    b = new_mpi_from_integer(-50_8)
    result_mpi = a - b
    IF (mpi_to_integer(result_mpi) == 150_8) THEN
      PRINT *, "Test 6 PASSED: 100 - (-50) = 150"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 6 FAILED: 100 - (-50)"
      fail_count = fail_count + 1
    END IF
    
    ! Test 7: Subtract two negative numbers
    test_count = test_count + 1
    a = new_mpi_from_integer(-100_8)
    b = new_mpi_from_integer(-50_8)
    result_mpi = a - b
    IF (mpi_to_integer(result_mpi) == -50_8) THEN
      PRINT *, "Test 7 PASSED: (-100) - (-50) = -50"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 7 FAILED: (-100) - (-50)"
      fail_count = fail_count + 1
    END IF
    
    ! Test 8: Subtraction with borrow across coefficients
    test_count = test_count + 1
    a%coeffs = [10_8, 1_8, 0_8, 0_8]
    b%coeffs = [20_8, 0_8, 0_8, 0_8]
    result_mpi = a - b
    IF (result_mpi%coeffs(1) == MULTI_PRECISION_BASE - 10_8 .AND. result_mpi%coeffs(2) == 0_8) THEN
      PRINT *, "Test 8 PASSED: Subtraction with borrow"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 8 FAILED: Subtraction with borrow"
      fail_count = fail_count + 1
    END IF
    
    ! Test 9: Large subtraction spanning multiple coefficients
    test_count = test_count + 1
    a%coeffs = [200_8, 300_8, 400_8, 500_8]
    b%coeffs = [100_8, 150_8, 200_8, 250_8]
    result_mpi = a - b
    IF (result_mpi%coeffs(1) == 100_8 .AND. result_mpi%coeffs(2) == 150_8 &
        .AND. result_mpi%coeffs(3) == 200_8 .AND. result_mpi%coeffs(4) == 250_8) THEN
      PRINT *, "Test 9 PASSED: Multi-coefficient subtraction"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 9 FAILED: Multi-coefficient subtraction"
      fail_count = fail_count + 1
    END IF
    
    ! Test 10: Verify a - b = -(b - a)
    test_count = test_count + 1
    a = new_mpi_from_integer(1234_8)
    b = new_mpi_from_integer(5678_8)
    IF ((a - b) == (-(b - a))) THEN
      PRINT *, "Test 10 PASSED: a - b = -(b - a)"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 10 FAILED: Subtraction anti-commutativity"
      fail_count = fail_count + 1
    END IF
    
    PRINT *, ""
  END SUBROUTINE test_mpi_subtract

  SUBROUTINE test_mpi_multiply(test_count, pass_count, fail_count)
    INTEGER, INTENT(INOUT) :: test_count, pass_count, fail_count
    TYPE(mpi) :: a, b, result_mpi
    
    PRINT *, "Testing mpi_multiply..."
    PRINT *, "----------------------------------------"
    
    ! Test 1: Multiply by zero
    test_count = test_count + 1
    a = new_mpi_from_integer(12345_8)
    b = new_mpi_from_integer(0_8)
    result_mpi = a * b
    IF (mpi_is_zero(result_mpi)) THEN
      PRINT *, "Test 1 PASSED: 12345 * 0 = 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 1 FAILED: Multiply by zero"
      fail_count = fail_count + 1
    END IF
    
    ! Test 2: Multiply by one (identity)
    test_count = test_count + 1
    a = new_mpi_from_integer(999_8)
    b = new_mpi_from_integer(1_8)
    result_mpi = a * b
    IF (mpi_to_integer(result_mpi) == 999_8) THEN
      PRINT *, "Test 2 PASSED: 999 * 1 = 999"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 2 FAILED: Multiply by one"
      fail_count = fail_count + 1
    END IF
    
    ! Test 3: Simple multiplication
    test_count = test_count + 1
    a = new_mpi_from_integer(123_8)
    b = new_mpi_from_integer(456_8)
    result_mpi = a * b
    IF (mpi_to_integer(result_mpi) == 56088_8) THEN
      PRINT *, "Test 3 PASSED: 123 * 456 = 56088"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 3 FAILED: Expected 56088, got ", mpi_to_integer(result_mpi)
      fail_count = fail_count + 1
    END IF
    
    ! Test 4: Multiply positive by negative
    test_count = test_count + 1
    a = new_mpi_from_integer(100_8)
    b = new_mpi_from_integer(-50_8)
    result_mpi = a * b
    IF (mpi_to_integer(result_mpi) == -5000_8) THEN
      PRINT *, "Test 4 PASSED: 100 * (-50) = -5000"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 4 FAILED: 100 * (-50)"
      fail_count = fail_count + 1
    END IF
    
    ! Test 5: Multiply two negative numbers
    test_count = test_count + 1
    a = new_mpi_from_integer(-25_8)
    b = new_mpi_from_integer(-40_8)
    result_mpi = a * b
    IF (mpi_to_integer(result_mpi) == 1000_8) THEN
      PRINT *, "Test 5 PASSED: (-25) * (-40) = 1000"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 5 FAILED: (-25) * (-40)"
      fail_count = fail_count + 1
    END IF
    
    ! Test 6: Multiply by power of 2
    test_count = test_count + 1
    a = new_mpi_from_integer(1000_8)
    b = new_mpi_from_integer(1024_8)
    result_mpi = a * b
    IF (mpi_to_integer(result_mpi) == 1024000_8) THEN
      PRINT *, "Test 6 PASSED: 1000 * 1024 = 1024000"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 6 FAILED: Power of 2 multiplication"
      fail_count = fail_count + 1
    END IF
    
    ! Test 7: Multiplication causing overflow to second coefficient
    test_count = test_count + 1
    a = new_mpi_from_integer(1000000_8)
    b = new_mpi_from_integer(10000_8)
    result_mpi = a * b
    IF (mpi_to_integer(result_mpi) == 10000000000_8) THEN
      PRINT *, "Test 7 PASSED: 1000000 * 10000 = 10000000000"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 7 FAILED: Expected 10000000000, got ", mpi_to_integer(result_mpi)
      fail_count = fail_count + 1
    END IF
    
    ! Test 8: Large multiplication spanning multiple coefficients
    test_count = test_count + 1
    a%coeffs = [1000_8, 2000_8, 0_8, 0_8]
    CALL normalize_mpi(a)
    b = new_mpi_from_integer(3_8)
    result_mpi = a * b
    IF (result_mpi%coeffs(1) == 3000_8 .AND. result_mpi%coeffs(2) == 6000_8) THEN
      PRINT *, "Test 8 PASSED: Multi-coefficient multiplication"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 8 FAILED: Multi-coefficient multiplication"
      fail_count = fail_count + 1
    END IF
    
    ! Test 9: Commutativity (a * b = b * a)
    test_count = test_count + 1
    a = new_mpi_from_integer(789_8)
    b = new_mpi_from_integer(321_8)
    IF ((a * b) == (b * a)) THEN
      PRINT *, "Test 9 PASSED: Multiplication commutativity"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 9 FAILED: Multiplication not commutative"
      fail_count = fail_count + 1
    END IF
    
    ! Test 10: Multiply large numbers (BASE-1 * 2)
    test_count = test_count + 1
    a%coeffs = [MULTI_PRECISION_BASE - 1_8, 0_8, 0_8, 0_8]
    b = new_mpi_from_integer(2_8)
    result_mpi = a * b
    IF (result_mpi%coeffs(1) == MULTI_PRECISION_BASE - 2_8 .AND. result_mpi%coeffs(2) == 1_8) THEN
      PRINT *, "Test 10 PASSED: (BASE-1) * 2"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 10 FAILED: (BASE-1) * 2"
      fail_count = fail_count + 1
    END IF
    
    ! Test 11: Square of a number
    test_count = test_count + 1
    a = new_mpi_from_integer(111_8)
    result_mpi = a * a
    IF (mpi_to_integer(result_mpi) == 12321_8) THEN
      PRINT *, "Test 11 PASSED: 111^2 = 12321"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 11 FAILED: 111^2"
      fail_count = fail_count + 1
    END IF
    
    PRINT *, ""
  END SUBROUTINE test_mpi_multiply

  SUBROUTINE test_mpi_div_rem(test_count, pass_count, fail_count)
    INTEGER, INTENT(INOUT) :: test_count, pass_count, fail_count
    TYPE(mpi) :: numerator, denominator, quotient, remainder
    
    PRINT *, "Testing mpi_div_rem..."
    PRINT *, "----------------------------------------"
    
    ! Test 1: Divide zero by non-zero
    test_count = test_count + 1
    numerator = new_mpi_from_integer(0_8)
    denominator = new_mpi_from_integer(5_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_is_zero(quotient) .AND. mpi_is_zero(remainder)) THEN
      PRINT *, "Test 1 PASSED: 0 / 5 = 0 remainder 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 1 FAILED: 0 / 5"
      fail_count = fail_count + 1
    END IF
    
    ! Test 2: Divide by one (identity)
    test_count = test_count + 1
    numerator = new_mpi_from_integer(999_8)
    denominator = new_mpi_from_integer(1_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_to_integer(quotient) == 999_8 .AND. mpi_is_zero(remainder)) THEN
      PRINT *, "Test 2 PASSED: 999 / 1 = 999 remainder 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 2 FAILED: 999 / 1"
      fail_count = fail_count + 1
    END IF
    
    ! Test 3: Simple division with no remainder
    test_count = test_count + 1
    numerator = new_mpi_from_integer(100_8)
    denominator = new_mpi_from_integer(10_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_to_integer(quotient) == 10_8 .AND. mpi_is_zero(remainder)) THEN
      PRINT *, "Test 3 PASSED: 100 / 10 = 10 remainder 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 3 FAILED: 100 / 10"
      fail_count = fail_count + 1
    END IF
    
    ! Test 4: Division with remainder
    test_count = test_count + 1
    numerator = new_mpi_from_integer(100_8)
    denominator = new_mpi_from_integer(7_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_to_integer(quotient) == 14_8 .AND. mpi_to_integer(remainder) == 2_8) THEN
      PRINT *, "Test 4 PASSED: 100 / 7 = 14 remainder 2"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 4 FAILED: 100 / 7, got quotient=", mpi_to_integer(quotient), &
               " remainder=", mpi_to_integer(remainder)
      fail_count = fail_count + 1
    END IF
    
    ! Test 5: Numerator less than denominator
    test_count = test_count + 1
    numerator = new_mpi_from_integer(5_8)
    denominator = new_mpi_from_integer(10_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_is_zero(quotient) .AND. mpi_to_integer(remainder) == 5_8) THEN
      PRINT *, "Test 5 PASSED: 5 / 10 = 0 remainder 5"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 5 FAILED: 5 / 10"
      fail_count = fail_count + 1
    END IF
    
    ! Test 6: Self-division (a / a = 1)
    test_count = test_count + 1
    numerator = new_mpi_from_integer(777_8)
    denominator = new_mpi_from_integer(777_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_to_integer(quotient) == 1_8 .AND. mpi_is_zero(remainder)) THEN
      PRINT *, "Test 6 PASSED: 777 / 777 = 1 remainder 0"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 6 FAILED: 777 / 777"
      fail_count = fail_count + 1
    END IF
    
    ! Test 7: Positive divided by negative
    test_count = test_count + 1
    numerator = new_mpi_from_integer(100_8)
    denominator = new_mpi_from_integer(-7_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_to_integer(quotient) == -14_8 .AND. mpi_to_integer(remainder) == 2_8) THEN
      PRINT *, "Test 7 PASSED: 100 / -7 = -14 remainder 2"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 7 FAILED: 100 / -7"
      fail_count = fail_count + 1
    END IF
    
    ! Test 8: Negative divided by positive
    test_count = test_count + 1
    numerator = new_mpi_from_integer(-100_8)
    denominator = new_mpi_from_integer(7_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_to_integer(quotient) == -14_8 .AND. mpi_to_integer(remainder) == -2_8) THEN
      PRINT *, "Test 8 PASSED: -100 / 7 = -14 remainder -2"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 8 FAILED: -100 / 7, got quotient=", mpi_to_integer(quotient), &
               " remainder=", mpi_to_integer(remainder)
      fail_count = fail_count + 1
    END IF
    
    ! Test 9: Both negative
    test_count = test_count + 1
    numerator = new_mpi_from_integer(-100_8)
    denominator = new_mpi_from_integer(-7_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_to_integer(quotient) == 14_8 .AND. mpi_to_integer(remainder) == -2_8) THEN
      PRINT *, "Test 9 PASSED: -100 / -7 = 14 remainder -2"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 9 FAILED: -100 / -7"
      fail_count = fail_count + 1
    END IF
    
    ! Test 10: Large number division
    test_count = test_count + 1
    numerator = new_mpi_from_integer(1000000_8)
    denominator = new_mpi_from_integer(333_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (mpi_to_integer(quotient) == 3003_8 .AND. mpi_to_integer(remainder) == 1_8) THEN
      PRINT *, "Test 10 PASSED: 1000000 / 333 = 3003 remainder 1"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 10 FAILED: 1000000 / 333"
      fail_count = fail_count + 1
    END IF
    
    ! Test 11: Division with multi-coefficient numerator
    test_count = test_count + 1
    numerator%coeffs = [MULTI_PRECISION_BASE, 0_8, 0_8, 0_8]
    CALL normalize_mpi(numerator)
    denominator = new_mpi_from_integer(2_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (quotient%coeffs(1) == MULTI_PRECISION_BASE / 2_8 .AND. quotient%coeffs(2) == 0_8 &
        .AND. mpi_is_zero(remainder)) THEN
      PRINT *, "Test 11 PASSED: BASE / 2"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 11 FAILED: BASE / 2"
      fail_count = fail_count + 1
    END IF
    
    ! Test 12: Verify division property: numerator = quotient * denominator + remainder
    test_count = test_count + 1
    numerator = new_mpi_from_integer(12345_8)
    denominator = new_mpi_from_integer(678_8)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (numerator == (quotient * denominator + remainder)) THEN
      PRINT *, "Test 12 PASSED: Division property verified"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 12 FAILED: Division property not satisfied"
      fail_count = fail_count + 1
    END IF
    
    ! Test 13: Division by large denominator (spans 2 coefficients)
    test_count = test_count + 1
    numerator%coeffs = [100_8, 200_8, 0_8, 0_8]
    CALL normalize_mpi(numerator)
    denominator%coeffs = [50_8, 100_8, 0_8, 0_8]
    CALL normalize_mpi(denominator)
    CALL mpi_div_rem(numerator, denominator, quotient, remainder)
    IF (numerator == (quotient * denominator + remainder)) THEN
      PRINT *, "Test 13 PASSED: Multi-coefficient division property"
      pass_count = pass_count + 1
    ELSE
      PRINT *, "Test 13 FAILED: Multi-coefficient division"
      fail_count = fail_count + 1
    END IF
    
    PRINT *, ""
  END SUBROUTINE test_mpi_div_rem


END PROGRAM test_mpi_subroutines