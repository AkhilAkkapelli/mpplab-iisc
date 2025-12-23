MODULE multi_precision_float_mod
  USE multi_precision_integer_mod, ONLY : mpi, COEFFS_LIMIT, MULTI_PRECISION_BASE, &
                                     normalize_mpi, new_mpi_from_coeffs, &
                                     mpi_abs, new_mpi_from_integer, mpi_to_string, &
                                     OPERATOR(==), OPERATOR(<), OPERATOR(+), OPERATOR(-), &
                                     OPERATOR(*), new_mpi_from_string, mpi_is_zero, &
                                     mpi_div_rem, mpi_sign, mpi_shift_bits_coeffs, mpi_size
  IMPLICIT NONE

  TYPE mpf
    TYPE(mpi)              :: mantissa
    INTEGER                :: exponent
  END TYPE mpf

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE mpf_equal
  END INTERFACE

  INTERFACE OPERATOR(<)
    MODULE PROCEDURE mpf_less
  END INTERFACE

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE mpf_add
  END INTERFACE

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE mpf_subtract
    MODULE PROCEDURE mpf_unary_negate
  END INTERFACE

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE mpf_multiply
  END INTERFACE

  INTERFACE OPERATOR(/)
    MODULE PROCEDURE mpf_divide
  END INTERFACE

  PRIVATE
  PUBLIC :: mpf, OPERATOR(==), OPERATOR(<), OPERATOR(+), &
            OPERATOR(-), OPERATOR(*), OPERATOR(/), new_mpf_from_mpi_exp, mpf_from_real16, &
            new_mpf_from_integer, mpf_is_zero, nearest_mpi_to_mpf, truncate_mpf_to_mpi, mpf_to_real16, &
            normalize_mpf, mpf_to_string, mpf_value_equal, mpf_abs, new_mpf_from_string

CONTAINS

  FUNCTION new_mpf_from_mpi_exp(mantissa_in, exponent_in) RESULT(mpf_out)
    TYPE(mpi), INTENT(IN)              :: mantissa_in
    INTEGER, OPTIONAL, INTENT(IN)        :: exponent_in
    TYPE(mpf)           :: mpf_out

    mpf_out%mantissa = mantissa_in
    IF (PRESENT(exponent_in)) THEN
      mpf_out%exponent = exponent_in
    ELSE
      mpf_out%exponent = 0
    END IF
    CALL normalize_mpf(mpf_out)
  END FUNCTION new_mpf_from_mpi_exp

  FUNCTION mpf_from_real16(r16_in) RESULT(mpf_out)
    REAL(16), INTENT(IN) :: r16_in
    TYPE(mpf) :: mpf_out

    INTEGER(KIND=8) :: bits(2)
    INTEGER(KIND=8) :: hi48
    INTEGER :: exp_raw, exp_unbiased
    TYPE(mpi) :: mantissa_mpi

    IF(r16_in == 0.0_16 .OR. COEFFS_LIMIT < 4 ) THEN
      mpf_out = new_mpf_from_integer(0_8)
      RETURN
    END IF

    bits = TRANSFER(r16_in, [0_8, 0_8])
    exp_raw = IAND(ISHFT(bits(2), -48), Z'7FFF')
    hi48 = IAND(bits(2), Z'0000FFFFFFFFFFFF')

    IF (exp_raw /= 0) hi48 = IOR(hi48, Z'0001000000000000')

    mantissa_mpi = new_mpi_from_integer(0_8)

    mantissa_mpi%coeffs(1) = IAND(bits(1),        Z'FFFFFFFF')
    mantissa_mpi%coeffs(2) = IAND(ISHFT(bits(1), -32), Z'FFFFFFFF')
    mantissa_mpi%coeffs(3) = IAND(hi48,           Z'FFFFFFFF')
    mantissa_mpi%coeffs(4) = IAND(ISHFT(hi48, -32), Z'FFFFFFFF')

    IF (r16_in < 0.0_16) mantissa_mpi = -mantissa_mpi

    IF (exp_raw == 0) THEN
      exp_unbiased = 1 - 16383
    ELSE
      exp_unbiased = exp_raw - 16383
    END IF

    mpf_out%mantissa = mantissa_mpi
    mpf_out%exponent = exp_unbiased - 112

    CALL normalize_mpf(mpf_out)

  END FUNCTION mpf_from_real16

  FUNCTION mpf_to_real16(mpf_in) RESULT(r16_out)
    TYPE(mpf), INTENT(IN) :: mpf_in
    REAL(16) :: r16_out

    REAL(16) :: mantissa_r16
    INTEGER :: i
    TYPE(mpi) :: abs_mantissa_mpi
    LOGICAL :: is_negative

    IF (mpf_is_zero(mpf_in)) THEN
      r16_out = 0.0_16
      RETURN
    END IF

    is_negative = mpi_sign(mpf_in%mantissa)
    abs_mantissa_mpi = mpi_abs(mpf_in%mantissa)

    mantissa_r16 = 0.0_16
    DO i = 1, MIN(COEFFS_LIMIT,4)
      mantissa_r16 = mantissa_r16 + REAL(abs_mantissa_mpi%coeffs(i), 16) * (REAL(MULTI_PRECISION_BASE, 16)**(i-1))
    END DO

    r16_out = SCALE(mantissa_r16, mpf_in%exponent)
    IF (is_negative) r16_out = -r16_out

  END FUNCTION mpf_to_real16

  FUNCTION nearest_mpi_to_mpf(mpf_in) RESULT(mpi_out)
    TYPE(mpf), INTENT(IN) :: mpf_in
    TYPE(mpi)             :: mpi_out

    INTEGER               :: k, word_idx, bit_idx, i
    INTEGER(KIND=8)       :: val_at_k, rounding_bit
    LOGICAL               :: sticky

    IF (mpf_in%exponent >= 0) THEN
      mpi_out = mpf_in%mantissa
      CALL mpi_shift_bits_coeffs(mpi_out%coeffs, mpf_in%exponent)
      RETURN
    END IF

    k = -mpf_in%exponent
    
    word_idx = (k - 1) / 32 + 1
    bit_idx  = MOD(k - 1, 32)

    val_at_k = ABS(mpf_in%mantissa%coeffs(word_idx))
    rounding_bit = IAND(ISHFT(val_at_k, -bit_idx), 1_8)

    sticky = IAND(val_at_k, (ISHFT(1_8, bit_idx) - 1_8)) /= 0
    
    IF (.NOT. sticky .AND. word_idx > 1) THEN
      DO i = 1, word_idx - 1
        IF (mpf_in%mantissa%coeffs(i) /= 0_8) THEN
          sticky = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    mpi_out = mpf_in%mantissa
    CALL mpi_shift_bits_coeffs(mpi_out%coeffs, -k)

    IF ((rounding_bit == 1_8) .AND. (sticky .OR. IAND(ABS(mpi_out%coeffs(1)), 1_8) == 1_8)) THEN
      IF (mpi_sign(mpf_in%mantissa)) THEN
        mpi_out = mpi_out - new_mpi_from_integer(1_8)
      ELSE
        mpi_out = mpi_out + new_mpi_from_integer(1_8)
      END IF
    ELSE
      CALL normalize_mpi(mpi_out)
    END IF

  END FUNCTION nearest_mpi_to_mpf

  FUNCTION truncate_mpf_to_mpi(mpf_in) RESULT(mpi_out)
    TYPE(mpf), INTENT(IN) :: mpf_in
    TYPE(mpi) :: mpi_out

    mpi_out = mpf_in%mantissa
    call mpi_shift_bits_coeffs(mpi_out%coeffs, mpf_in%exponent)

  END FUNCTION truncate_mpf_to_mpi
  
  SUBROUTINE normalize_mpf(mpf_val)
    TYPE(mpf), INTENT(INOUT) :: mpf_val

    INTEGER(KIND=8) :: carry, current_coeff_val, sign
    INTEGER :: i, msb_idx, trailing_zeros, shift_bits
    INTEGER(KIND=8) :: combined_coeffs(COEFFS_LIMIT + 1)


    IF (mpi_is_zero(mpf_val%mantissa)) THEN
      mpf_val%exponent = 0
      RETURN
    END IF

    combined_coeffs(1:COEFFS_LIMIT) = mpf_val%mantissa%coeffs(1:COEFFS_LIMIT)
    combined_coeffs(COEFFS_LIMIT + 1) = 0_8
    
    carry = 0_8
    DO i = 1, COEFFS_LIMIT + 1
      current_coeff_val = combined_coeffs(i) + carry
      carry = current_coeff_val / MULTI_PRECISION_BASE
      combined_coeffs(i) = current_coeff_val - carry * MULTI_PRECISION_BASE
      IF (combined_coeffs(i) /= 0_8 .AND. carry == 0_8) msb_idx = i
    END DO

    sign = KISIGN(1_8, combined_coeffs(msb_idx))
    combined_coeffs = sign * combined_coeffs

    DO i = 1, msb_idx - 1
      IF (combined_coeffs(i) < 0_8) THEN
        combined_coeffs(i+1) = combined_coeffs(i+1) - sign
        combined_coeffs(i) = combined_coeffs(i) + sign * MULTI_PRECISION_BASE
      END IF
    END DO

    trailing_zeros = 0
    DO i = 1, COEFFS_LIMIT + 1
      IF (combined_coeffs(i) /= 0_8) THEN
        trailing_zeros = TRAILZ(combined_coeffs(i)) + 32 * (i - 1)
        EXIT
      END IF
    END DO

    IF (trailing_zeros > 0) THEN
      call mpi_shift_bits_coeffs(combined_coeffs, -trailing_zeros)
      mpf_val%exponent = mpf_val%exponent + trailing_zeros
    END IF

    IF (combined_coeffs(COEFFS_LIMIT + 1) /= 0_8) THEN
      shift_bits = 64 - LEADZ(combined_coeffs(COEFFS_LIMIT + 1))
      call mpi_shift_bits_coeffs(combined_coeffs, shift_bits)
      mpf_val%exponent = mpf_val%exponent + shift_bits
    END IF

    mpf_val%mantissa%coeffs(1:COEFFS_LIMIT) = sign*combined_coeffs(1:COEFFS_LIMIT)
    
  END SUBROUTINE normalize_mpf

  FUNCTION new_mpf_from_integer(x_in) RESULT(mpf_out)
    INTEGER(KIND=8), INTENT(IN) :: x_in
    TYPE(mpf)                   :: mpf_out

    mpf_out = new_mpf_from_mpi_exp(new_mpi_from_integer(x_in))
  END FUNCTION new_mpf_from_integer

  FUNCTION mpf_abs(mpf_in) RESULT(mpf_out)
    TYPE(mpf), INTENT(IN) :: mpf_in
    TYPE(mpf)             :: mpf_out
    mpf_out%exponent = mpf_in%exponent
    mpf_out%mantissa = mpi_abs(mpf_in%mantissa)
  END FUNCTION mpf_abs

  FUNCTION mpf_equal(a, b) RESULT(res)
    TYPE(mpf), INTENT(IN) :: a, b
    LOGICAL               :: res
    res = (a%mantissa == b%mantissa) .AND. (a%exponent == b%exponent)
  END FUNCTION mpf_equal

  FUNCTION mpf_less(a, b) RESULT(res)
    TYPE(mpf), INTENT(IN) :: a, b
    LOGICAL               :: res

    INTEGER               :: common_exponent
    TYPE(mpi)             :: scaled_a, scaled_b

    IF(mpi_sign(a%mantissa) .NEQV. mpi_sign(b%mantissa)) THEN
      res = mpi_sign(a%mantissa)
      RETURN
    ELSE IF(mpf_is_zero(a)) THEN
      res = .FALSE.
      RETURN
    END IF

    common_exponent = MIN(a%exponent, b%exponent)

    scaled_a = a%mantissa
    scaled_b = b%mantissa

    call mpi_shift_bits_coeffs(scaled_a%coeffs, a%exponent - common_exponent)
    call mpi_shift_bits_coeffs(scaled_b%coeffs, b%exponent - common_exponent)

    res = (scaled_a < scaled_b)
  END FUNCTION mpf_less

  LOGICAL FUNCTION mpf_value_equal(a, b)
    TYPE(mpf), INTENT(IN) :: a, b

    INTEGER :: common_exp

    TYPE(mpi) :: scaled_a, scaled_b

    IF (mpf_is_zero(a) .AND. mpf_is_zero(b)) THEN
        mpf_value_equal = .TRUE.
        RETURN
    END IF

    common_exp = MIN(a%exponent, b%exponent)

    scaled_a = a%mantissa
    scaled_b = b%mantissa

    call mpi_shift_bits_coeffs(scaled_a%coeffs, a%exponent - common_exp)
    call mpi_shift_bits_coeffs(scaled_b%coeffs, b%exponent - common_exp)

    mpf_value_equal = (scaled_a == scaled_b)
  END FUNCTION mpf_value_equal

  FUNCTION mpf_add(a, b) RESULT(mpf_sum)
    TYPE(mpf), INTENT(IN) :: a, b
    TYPE(mpf)             :: mpf_sum

    INTEGER               :: common_exponent
    TYPE(mpi)             :: scaled_a, scaled_b

    IF (mpf_is_zero(a)) THEN
      mpf_sum = b
      RETURN
    END IF
    IF (mpf_is_zero(b)) THEN
      mpf_sum = a
      RETURN
    END IF

    common_exponent = MIN(a%exponent, b%exponent)

    scaled_a = a%mantissa
    scaled_b = b%mantissa

    call mpi_shift_bits_coeffs(scaled_a%coeffs, a%exponent - common_exponent)
    call mpi_shift_bits_coeffs(scaled_b%coeffs, b%exponent - common_exponent)

    mpf_sum = new_mpf_from_mpi_exp(scaled_a + scaled_b, common_exponent)
  END FUNCTION mpf_add

  FUNCTION mpf_unary_negate(mpf_in) RESULT(mpf_out)
    TYPE(mpf), INTENT(IN) :: mpf_in
    TYPE(mpf)             :: mpf_out

    mpf_out%exponent = mpf_in%exponent
    mpf_out%mantissa = -mpf_in%mantissa
  END FUNCTION mpf_unary_negate

  FUNCTION mpf_subtract(a, b) RESULT(mpf_diff)
    TYPE(mpf), INTENT(IN) :: a, b
    TYPE(mpf)             :: mpf_diff
    mpf_diff = a + -b
  END FUNCTION mpf_subtract

  FUNCTION mpf_multiply(a, b) RESULT(mpf_prod)
    TYPE(mpf), INTENT(IN) :: a, b
    TYPE(mpf)             :: mpf_prod
    
    mpf_prod%mantissa = a%mantissa * b%mantissa
    mpf_prod%exponent = a%exponent + b%exponent
  END FUNCTION mpf_multiply

  FUNCTION mpf_divide(a, b) RESULT(mpf_quotient)
    TYPE(mpf), INTENT(IN) :: a, b
    TYPE(mpf)             :: mpf_quotient

    TYPE(mpi)             :: scaled_a, result_mantissa, remainder_ignored
    INTEGER               :: scale_bits
    INTEGER               :: i, bit_pos
    INTEGER(KIND=8)       :: temp_val

    IF (mpf_is_zero(b)) THEN
      STOP "Error: Division by zero MultiPrecisionFloat"
    END IF
    IF (mpf_is_zero(a)) THEN
      mpf_quotient = new_mpf_from_integer(0_8)
      RETURN
    END IF

    temp_val = 0
    bit_pos = -1
    DO i = COEFFS_LIMIT, 1, -1
      temp_val = IAND(ABS(a%mantissa%coeffs(i)), INT(Z'FFFFFFFF', KIND=8))
      IF (temp_val /= 0_8) THEN
        DO bit_pos = 31, 0, -1
          IF (BTEST(temp_val, bit_pos)) EXIT
        END DO
        bit_pos = (i - 1) * 32 + bit_pos
        EXIT
      END IF
    END DO

    scale_bits = MAX(0, (COEFFS_LIMIT * 32 - 8) - bit_pos)

    scaled_a = a%mantissa
    call mpi_shift_bits_coeffs(scaled_a%coeffs, scale_bits)

    CALL mpi_div_rem(scaled_a, b%mantissa, result_mantissa, remainder_ignored)

    mpf_quotient = new_mpf_from_mpi_exp(result_mantissa, a%exponent - b%exponent - scale_bits)

  END FUNCTION mpf_divide

  ELEMENTAL FUNCTION mpf_is_zero(val) RESULT(is_zero)
    TYPE(mpf), INTENT(IN) :: val
    LOGICAL :: is_zero
    is_zero = mpi_is_zero(val%mantissa) .AND. (val%exponent == 0)
  END FUNCTION mpf_is_zero

  FUNCTION new_mpf_from_string(s_in) RESULT(mpf_out)
    CHARACTER(LEN=*), INTENT(IN) :: s_in
    TYPE(mpf)   :: mpf_out
    CHARACTER(LEN=:), ALLOCATABLE :: cmd
    CHARACTER(LEN=2048) :: buffer
    CHARACTER(LEN=512) :: temp_file
    CHARACTER(LEN=10) :: coeffs_limit_str
    INTEGER :: stat, unit_num, io_stat
    REAL :: r
    INTEGER(KIND=8) :: c(4)
    INTEGER :: e, sign_flag

    CALL RANDOM_NUMBER(r)
    WRITE(temp_file, '(A,I0,A)') 'mpf_temp_', INT(r*1000000), '.txt'
    WRITE(coeffs_limit_str, '(I0)') COEFFS_LIMIT
    cmd = './mpf_from_string.sh "' // TRIM(s_in) // '" "' // TRIM(temp_file) // '" ' // TRIM(coeffs_limit_str)
    
    CALL EXECUTE_COMMAND_LINE(TRIM(cmd), WAIT=.TRUE., EXITSTAT=stat)
    
    IF (stat /= 0) THEN
      PRINT *, "ERROR: Julia command failed with exit status:", stat
      mpf_out%mantissa = new_mpi_from_coeffs([INT(0,8), INT(0,8), &
                                              INT(0,8), INT(0,8)])
      mpf_out%exponent = 0
      RETURN
    END IF
    
    OPEN(NEWUNIT=unit_num, FILE=TRIM(temp_file), STATUS='OLD', ACTION='READ', IOSTAT=io_stat)
    
    IF (io_stat /= 0) THEN
      PRINT *, "ERROR: Could not open temporary file:", TRIM(temp_file)
      mpf_out%mantissa = new_mpi_from_coeffs([INT(0,8), INT(0,8), INT(0,8), INT(0,8)])
      mpf_out%exponent = 0
      RETURN
    END IF
    
    READ(unit_num, '(A)', IOSTAT=io_stat) buffer
    CLOSE(unit_num, STATUS='DELETE')
    
    IF (io_stat /= 0) THEN
      PRINT *, "ERROR: Could not read from temporary file"
      mpf_out%mantissa = new_mpi_from_coeffs([INT(0,8), INT(0,8), &
                                              INT(0,8), INT(0,8)])
      mpf_out%exponent = 0
      RETURN
    END IF
    
    READ(buffer, *, IOSTAT=io_stat) c(1), c(2), c(3), c(4), e, sign_flag
    
    IF (io_stat /= 0) THEN
      PRINT *, "ERROR: Could not parse output:", TRIM(buffer)
      mpf_out%mantissa = new_mpi_from_coeffs([INT(0,8), INT(0,8), &
                                              INT(0,8), INT(0,8)])
      mpf_out%exponent = 0
      RETURN
    END IF
    
    mpf_out%mantissa = new_mpi_from_coeffs(c)
    mpf_out%exponent = e
    
    IF (sign_flag == 1) THEN
      mpf_out%mantissa = -mpf_out%mantissa 
    END IF

    CALL normalize_mpf(mpf_out)

  END FUNCTION new_mpf_from_string

  FUNCTION mpf_to_string(mpf_val) RESULT(str)
    TYPE(mpf), INTENT(IN) :: mpf_val
    CHARACTER(LEN=:), ALLOCATABLE :: str, cmd
    CHARACTER(LEN=2048) :: buffer
    CHARACTER(LEN=512) :: temp_file
    CHARACTER(LEN=32) :: coeff_str(4), exp_str, sign_str
    INTEGER :: stat, unit_num, io_stat, i
    REAL :: r
    INTEGER(KIND=8) :: c(4)
    INTEGER :: sign_flag
    TYPE(mpf) :: mpf_work
    
    mpf_work = mpf_val
    sign_flag = 0
    
    IF (mpi_sign(mpf_val%mantissa)) THEN
      sign_flag = 1
      mpf_work%mantissa = mpi_abs(mpf_val%mantissa)
    END IF
    
    IF (mpf_is_zero(mpf_work)) THEN
        str = "0.0"
        RETURN
    END IF

    c = mpf_work%mantissa%coeffs(1:4)
    
    DO i = 1, 4
      WRITE(coeff_str(i), '(I0)') c(i)
    END DO
    WRITE(exp_str, '(I0)') mpf_work%exponent
    WRITE(sign_str, '(I0)') sign_flag
    
    CALL RANDOM_NUMBER(r)
    WRITE(temp_file, '(A,I0,A)') 'mpf_temp_', INT(r*1000000), '.txt'

    cmd = './mpf_to_string.sh ' // &
        TRIM(coeff_str(1)) // ' ' // TRIM(coeff_str(2)) // ' ' // &
        TRIM(coeff_str(3)) // ' ' // TRIM(coeff_str(4)) // ' ' // &
        TRIM(exp_str)      // ' ' // TRIM(sign_str)     // ' "' // TRIM(temp_file) // '"'
    
    CALL EXECUTE_COMMAND_LINE(TRIM(cmd), WAIT=.TRUE., EXITSTAT=stat)  
    
    IF (stat /= 0) THEN
      str = "ERROR: Julia conversion failed"
      RETURN
    END IF
    
    OPEN(NEWUNIT=unit_num, FILE=TRIM(temp_file), STATUS='OLD', ACTION='READ', IOSTAT=io_stat)
    
    IF (io_stat /= 0) THEN
      str = "ERROR: Could not open output file"
      RETURN
    END IF
    
    READ(unit_num, '(A)', IOSTAT=io_stat) buffer
    CLOSE(unit_num, STATUS='DELETE') 
    
    IF (io_stat /= 0) THEN
      str = "ERROR: Could not read output"
      RETURN
    END IF
    
    str = TRIM(buffer)
    
  END FUNCTION mpf_to_string

END MODULE multi_precision_float_mod
