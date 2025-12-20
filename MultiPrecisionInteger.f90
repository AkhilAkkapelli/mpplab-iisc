MODULE multi_precision_integer_mod
  IMPLICIT NONE
  PRIVATE
  
  INTEGER,         PARAMETER :: COEFFS_LIMIT         = 4
  INTEGER(KIND=8), PARAMETER :: MULTI_PRECISION_BASE = 2_8**32

  TYPE mpi
    INTEGER(KIND=8) :: coeffs(COEFFS_LIMIT)
  END TYPE mpi



  INTERFACE OPERATOR(+)
    MODULE PROCEDURE mpi_add
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE mpi_unary_negate
    MODULE PROCEDURE mpi_subtract
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE mpi_multiply
  END INTERFACE
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE mpi_equal
  END INTERFACE
  INTERFACE OPERATOR(<)
    MODULE PROCEDURE mpi_less
  END INTERFACE

  PUBLIC :: COEFFS_LIMIT, MULTI_PRECISION_BASE, mpi, normalize_mpi, new_mpi_from_coeffs, new_mpi_from_integer, &
            mpi_to_integer, mpi_to_string, new_mpi_from_string, mpi_div_rem, &
            mpi_shift_bits_right, mpi_is_zero, mpi_sign, mpi_size, mpi_max_value, mpi_abs, mpi_shift_bits_left, &
            mpi_multiply_by_scalar, mpi_div_by_scalar, &
            OPERATOR(+), OPERATOR(-), OPERATOR(*), OPERATOR(==), OPERATOR(<)
CONTAINS

! meant to conevert mpi_val coeffs to be proper
SUBROUTINE normalize_mpi(mpi_val)
  TYPE(mpi), INTENT(INOUT) :: mpi_val

  INTEGER(KIND=8) :: carry, current_coeff_val
  INTEGER :: i, msb_idx, sign

  IF(mpi_is_zero(mpi_val)) RETURN

  carry = 0_8
  DO i = 1, COEFFS_LIMIT
    current_coeff_val = mpi_val%coeffs(i) + carry
    carry = current_coeff_val / MULTI_PRECISION_BASE
    mpi_val%coeffs(i) = current_coeff_val - carry * MULTI_PRECISION_BASE
    IF (mpi_val%coeffs(i) /= 0_8) msb_idx = i
  END DO

  IF(carry /= 0_8) print*, "overflow"

  IF (mpi_val%coeffs(msb_idx) > 0_8) THEN
    sign = 1_8
  ELSE
    sign = -1_8
  END IF

  DO i = 1, msb_idx - 1
    IF (sign * mpi_val%coeffs(i) < 0_8) THEN
      mpi_val%coeffs(i+1) = mpi_val%coeffs(i+1) - sign
      mpi_val%coeffs(i) = mpi_val%coeffs(i) + sign * MULTI_PRECISION_BASE
    END IF
  END DO

END SUBROUTINE normalize_mpi

! convert a set of coeffs to mpi
FUNCTION new_mpi_from_coeffs(coeffs) RESULT(mpi_out)
  INTEGER(KIND=8), INTENT(IN) :: coeffs(:)
  TYPE(mpi) :: mpi_out
  INTEGER :: n

  n = size(coeffs)
  mpi_out%coeffs = 0_8
  IF (n == 0) RETURN
  mpi_out%coeffs(1:n) = coeffs

  CALL normalize_mpi(mpi_out)

END FUNCTION new_mpi_from_coeffs

! convert an integer(8) to an mpi
FUNCTION new_mpi_from_integer(x_in) RESULT(mpi_out)
  INTEGER(KIND=8), INTENT(IN)   :: x_in
  TYPE(mpi)                     :: mpi_out

  INTEGER(KIND=8), PARAMETER :: MASK32 = INT(Z'FFFFFFFF', KIND=8)
  INTEGER(KIND=8), PARAMETER :: MAX_NEGATIVE_I8 = -HUGE(0_8) - 1_8
  INTEGER(KIND=8) :: mag_x

  mpi_out%coeffs = 0_8

  IF (x_in == 0_8) THEN
    RETURN
  ELSE IF (x_in == MAX_NEGATIVE_I8) THEN
    mpi_out%coeffs(1) = 0_8
    mpi_out%coeffs(2) = -ISHFT(MULTI_PRECISION_BASE, -1)
  ELSE
    mag_x = ABS(x_in)
    mpi_out%coeffs(1) = IAND(mag_x, MASK32)
    IF (mag_x >= MULTI_PRECISION_BASE) mpi_out%coeffs(2) = ISHFT(mag_x, -32)
    IF (x_in < 0_8) mpi_out%coeffs = -mpi_out%coeffs
    CALL normalize_mpi(mpi_out)
  END IF
  
END FUNCTION new_mpi_from_integer

! convert an mpi to an integer(8)
FUNCTION mpi_to_integer(mpi_in) RESULT(x)
  TYPE(mpi), INTENT(IN)   :: mpi_in
  integer(kind=8)         :: coeffs(COEFFS_LIMIT)
  INTEGER(KIND=8)         :: x, sign

  ! INTEGER :: i

  ! DO i = 3, COEFFS_LIMIT
  !   IF (mpi_in%coeffs(i) /= 0_8) PRINT*, "ERROR in mpi_to_integer: Value is too large to fit."
  !   EXIT
  ! END DO
  sign= merge(1,-1, all(mpi_in%coeffs >= 0_8))
  coeffs= abs(mpi_in%coeffs)
  x = coeffs(1) + ISHFT(coeffs(2), 32)
  x = x*sign
END FUNCTION mpi_to_integer

! multiply by a scalar
SUBROUTINE mpi_multiply_by_scalar(mpi_val, scalar)
  TYPE(mpi), INTENT(INOUT)      :: mpi_val
  INTEGER(KIND=8), INTENT(IN)   :: scalar
  
  INTEGER(KIND=8) :: carry, prod
  INTEGER :: i

  ! sign reslution variables
  integer(kind=8)               :: sign_info_pos, sign_info_neg
  integer(kind=8)               :: sign_mpi
  integer(kind=8)               :: sign_scl, abs_scl

  ! resolve sign of mpi_val
  sign_info_pos= merge(1, 0, all(mpi_val%coeffs >= 0_8)) 
  sign_info_neg= merge(1, 0, all(mpi_val%coeffs <=  0_8)) 
  if (.not.(sign_info_pos == 1 .or. sign_info_neg == 1))then
    print *, "Provide valid MPI Value!"
    stop
  end if
  sign_mpi = sign_info_pos - sign_info_neg
  mpi_val%coeffs= abs(mpi_val%coeffs)

  ! resolve sign of scalar
  sign_scl= merge(1,-1,scalar >= 0_8)
  abs_scl= abs(scalar)

  carry = 0_8
  DO i = 1, COEFFS_LIMIT
    prod = mpi_val%coeffs(i) * abs_scl + carry
    ! mpi_val%coeffs(i)= modulo(prod, MULTI_PRECISION_BASE)
    mpi_val%coeffs(i) = IAND(prod, MULTI_PRECISION_BASE - 1)
    ! carry= prod/MULTI_PRECISION_BASE
    carry = ISHFT(prod, -32)
  END DO

  ! ! resolve final sign
  mpi_val%coeffs= mpi_val%coeffs*sign_mpi*sign_scl
END SUBROUTINE mpi_multiply_by_scalar

! divide by a scalar
FUNCTION mpi_div_by_scalar(mpi_val, divisor) RESULT(remainder)
  TYPE(mpi), INTENT(INOUT)      :: mpi_val
  INTEGER(KIND=8), INTENT(IN)   :: divisor
  INTEGER(KIND=8)               :: remainder

  INTEGER(KIND=8) :: current_val
  INTEGER :: i

  remainder = 0_8
  DO i = COEFFS_LIMIT, 1, -1
      current_val = mpi_val%coeffs(i) + ISHFT(remainder, 32)
      mpi_val%coeffs(i) = current_val / divisor
      remainder = MODULO(current_val, divisor)
  END DO
END FUNCTION mpi_div_by_scalar

FUNCTION mpi_shift_bits_left(mpi_in, num_bits) RESULT(mpi_out)
  TYPE(mpi), INTENT(IN) :: mpi_in
  INTEGER, INTENT(IN)   :: num_bits
  TYPE(mpi)             :: mpi_out

  TYPE(mpi) :: mpi_abs_in
  LOGICAL :: is_negative
  INTEGER :: word_shift, bit_shift, i
  INTEGER(KIND=8) :: low_part, high_part

  IF (num_bits <= 0) THEN
    mpi_out = mpi_in
    RETURN
  END IF

  IF (mpi_is_zero(mpi_in)) THEN
      mpi_out%coeffs = 0_8
      RETURN
  END IF

  is_negative = mpi_sign(mpi_in)
  mpi_abs_in = mpi_abs(mpi_in)

  word_shift = SHIFTR(num_bits, 5)
  bit_shift = IAND(num_bits, 31)
  
  IF (word_shift > 0) THEN
    mpi_out%coeffs(1:min(word_shift, COEFFS_LIMIT)) = 0_8
  END IF

  IF (word_shift >= COEFFS_LIMIT) RETURN
  
  high_part = 0_8
  DO i = 1, COEFFS_LIMIT - word_shift
      low_part = IAND(ISHFT(mpi_abs_in%coeffs(i), bit_shift), INT(Z'FFFFFFFF', KIND=8))
      IF (i + word_shift <= COEFFS_LIMIT) mpi_out%coeffs(i + word_shift) = IOR(low_part, high_part)
      high_part = ISHFT(ISHFT(mpi_abs_in%coeffs(i), bit_shift), -32)
  END DO

  IF (is_negative) mpi_out = -mpi_out
  
END FUNCTION mpi_shift_bits_left

FUNCTION mpi_shift_bits_right(mpi_in, num_bits) RESULT(mpi_out)
  TYPE(mpi), INTENT(IN) :: mpi_in
  INTEGER, INTENT(IN)   :: num_bits
  TYPE(mpi)             :: mpi_out

  TYPE(mpi) :: mpi_abs_in
  LOGICAL :: is_negative
  INTEGER :: word_shift, bit_shift, i
  INTEGER(KIND=8) :: low_part, high_part

  IF (num_bits <= 0) THEN
    mpi_out = mpi_in
    RETURN
  END IF

  IF (mpi_is_zero(mpi_in)) THEN
    mpi_out%coeffs = 0_8
    RETURN
  END IF

  is_negative = mpi_sign(mpi_in)
  mpi_abs_in = mpi_abs(mpi_in)

  word_shift = SHIFTR(num_bits, 5)
  bit_shift = IAND(num_bits, 31)

  mpi_out%coeffs = 0_8
  IF (word_shift >= COEFFS_LIMIT) RETURN

  low_part = 0_8
  DO i = COEFFS_LIMIT, word_shift + 1, -1
    high_part = IAND(ISHFT(mpi_abs_in%coeffs(i), -bit_shift), INT(Z'FFFFFFFF', KIND=8))
    mpi_out%coeffs(i - word_shift) = IOR(high_part, low_part)
    low_part = ISHFT(IAND(mpi_abs_in%coeffs(i), ISHFT(1_8, bit_shift) - 1_8), 32 - bit_shift)
  END DO

  IF (is_negative) mpi_out = -mpi_out
    
END FUNCTION mpi_shift_bits_right

SUBROUTINE mpi_div_rem(numerator, denominator, quotient, remainder)
  TYPE(mpi), INTENT(IN)  :: numerator, denominator
  TYPE(mpi), INTENT(OUT) :: quotient, remainder

  TYPE(mpi)     :: num_abs, den_abs, q_abs, r_abs
  LOGICAL :: num_is_neg, den_is_neg
  LOGICAL :: q_is_neg, r_is_neg

  IF (mpi_is_zero(denominator)) STOP "mpi_div_rem: Division by zero."

  IF (mpi_is_zero(numerator)) THEN
    quotient = new_mpi_from_integer(0_8)
    remainder = new_mpi_from_integer(0_8)
    RETURN
  END IF

  num_abs = mpi_abs(numerator)
  den_abs = mpi_abs(denominator)

  IF (num_abs < den_abs) THEN
    quotient = new_mpi_from_integer(0_8)
    remainder = numerator
    RETURN
  END IF

  CALL div_rem_magnitude(num_abs, den_abs, q_abs, r_abs)

  num_is_neg = mpi_sign(numerator)
  den_is_neg = mpi_sign(denominator)

  q_is_neg = (num_is_neg .NEQV. den_is_neg)
  r_is_neg = num_is_neg

  quotient = q_abs
  IF (q_is_neg) quotient = -quotient

  remainder = r_abs
  IF (r_is_neg) remainder = -remainder

CONTAINS
  SUBROUTINE div_rem_magnitude(num, den, q, r)
    TYPE(mpi), INTENT(IN)  :: num, den
    TYPE(mpi), INTENT(OUT) :: q, r

    TYPE(mpi) :: current_rem, shifted_den, test_prod
    INTEGER :: i, n, m, k
    INTEGER(KIND=8) :: q_digit, low, high, mid

    n = COEFFS_LIMIT; DO WHILE (n > 1 .AND. num%coeffs(n) == 0_8); n = n - 1; END DO
    m = COEFFS_LIMIT; DO WHILE (m > 1 .AND. den%coeffs(m) == 0_8); m = m - 1; END DO

    q%coeffs = 0_8
    current_rem = num

    DO k = n - m, 0, -1
      shifted_den = mpi_scale_up_by_base_power(den, k)

      IF (current_rem < shifted_den) THEN
          q_digit = 0_8
      ELSE
          low = 1_8
          high = MULTI_PRECISION_BASE - 1
          q_digit = 1_8

          DO WHILE (low <= high)
              mid = low + ISHFT(high - low, -1)
              test_prod = new_mpi_from_integer(mid)
              test_prod = shifted_den * test_prod
              IF (current_rem < test_prod) THEN
                  high = mid - 1
              ELSE
                  q_digit = mid
                  low = mid + 1
              END IF
          END DO
      END IF

      IF (q_digit > 0) THEN
          q%coeffs(k+1) = q_digit
          test_prod = new_mpi_from_integer(q_digit)
          current_rem = current_rem - (shifted_den * test_prod)
      END IF
    END DO

    CALL normalize_mpi(q)
    r = current_rem
  END SUBROUTINE div_rem_magnitude

  FUNCTION mpi_scale_up_by_base_power(mpi_in, power) RESULT(mpi_out)
    TYPE(mpi), INTENT(IN) :: mpi_in
    INTEGER, INTENT(IN)     :: power
    TYPE(mpi)             :: mpi_out

    IF (power < 0 .OR. power >= COEFFS_LIMIT) THEN
      STOP "mpi_scale_up_by_base_power: Invalid power."
    END IF
    IF (power == 0) THEN
      mpi_out = mpi_in
      RETURN
    END IF

    mpi_out%coeffs = 0_8
    mpi_out%coeffs(power+1:COEFFS_LIMIT) = mpi_in%coeffs(1:COEFFS_LIMIT-power)
  END FUNCTION mpi_scale_up_by_base_power
END SUBROUTINE mpi_div_rem

PURE FUNCTION mpi_is_zero(mpi_val) RESULT(is_zero)
  TYPE(mpi), INTENT(IN) :: mpi_val
  LOGICAL :: is_zero

  is_zero = ALL(mpi_val%coeffs == 0_8)
END FUNCTION mpi_is_zero

FUNCTION mpi_sign(mpi_val) RESULT(is_negative)
  TYPE(mpi), INTENT(IN) :: mpi_val
  LOGICAL :: is_negative
  INTEGER :: i
  is_negative = .FALSE.
  DO i = COEFFS_LIMIT, 1, -1
    IF (mpi_val%coeffs(i) /= 0_8) THEN
      is_negative = (mpi_val%coeffs(i) < 0)
      EXIT
    END IF
  END DO
END FUNCTION mpi_sign

FUNCTION mpi_size(mpi_in) RESULT(size_val)
  TYPE(mpi), INTENT(IN) :: mpi_in
  INTEGER :: size_val
  INTEGER :: i

  size_val = 0
  DO i = COEFFS_LIMIT, 1, -1
    IF (mpi_in%coeffs(i) /= 0_8) THEN
      size_val = i
      EXIT
    END IF
  END DO
END FUNCTION mpi_size

FUNCTION mpi_max_value() RESULT(mpi_out)
  TYPE(mpi) :: mpi_out
  INTEGER :: i
  DO i = 1, COEFFS_LIMIT
    mpi_out%coeffs(i) = MULTI_PRECISION_BASE - 1
  END DO
END FUNCTION mpi_max_value

FUNCTION mpi_equal(a, b) RESULT(res)
  TYPE(mpi), INTENT(IN) :: a, b
  LOGICAL :: res

  res = mpi_is_zero(a - b)
END FUNCTION mpi_equal

FUNCTION mpi_less(a, b) RESULT(res)
  TYPE(mpi), INTENT(IN) :: a, b
  LOGICAL :: res
  TYPE(mpi) :: diff
  diff = b - a 
  res = (.NOT. mpi_is_zero(diff)) .AND. (.NOT. mpi_sign(diff))
END FUNCTION mpi_less

FUNCTION mpi_unary_negate(a) RESULT(res)
  TYPE(mpi), INTENT(IN) :: a
  TYPE(mpi) :: res
  res%coeffs = -a%coeffs
  CALL normalize_mpi(res)
END FUNCTION mpi_unary_negate

FUNCTION mpi_add(a, b) RESULT(res)
  TYPE(mpi), INTENT(IN) :: a, b
  TYPE(mpi) :: res
  res%coeffs = a%coeffs + b%coeffs
  CALL normalize_mpi(res)
END FUNCTION mpi_add

FUNCTION mpi_subtract(a, b) RESULT(res)
  TYPE(mpi), INTENT(IN) :: a, b
  TYPE(mpi) :: res
  res%coeffs = a%coeffs - b%coeffs
  CALL normalize_mpi(res)
END FUNCTION mpi_subtract

FUNCTION mpi_abs(mpi_in) RESULT(mpi_out)
    TYPE(mpi), INTENT(IN) :: mpi_in
    TYPE(mpi) :: mpi_out
    IF (.NOT. mpi_sign(mpi_in)) THEN
        mpi_out = mpi_in
    ELSE
        mpi_out%coeffs = -mpi_in%coeffs
        CALL normalize_mpi(mpi_out)
    END IF
END FUNCTION mpi_abs

FUNCTION mpi_multiply(a, b) RESULT(prod_mpi)
  TYPE(mpi), INTENT(IN) :: a, b
  TYPE(mpi) :: prod_mpi

  TYPE(mpi) :: abs_a, abs_b
  TYPE(mpi) :: temp_prod
  LOGICAL :: is_neg_a, is_neg_b, final_is_neg
  INTEGER :: i, j
  INTEGER(KIND=8) :: carry, prod_chunk

  IF (mpi_is_zero(a) .OR. mpi_is_zero(b)) THEN
    prod_mpi%coeffs = 0_8
    RETURN
  END IF

  is_neg_a = mpi_sign(a)
  is_neg_b = mpi_sign(b)
  final_is_neg = (is_neg_a .NEQV. is_neg_b)
  abs_a = mpi_abs(a)
  abs_b = mpi_abs(b)

  prod_mpi%coeffs = 0_8

  DO i = 1, COEFFS_LIMIT
      IF (abs_a%coeffs(i) == 0_8) CYCLE

      temp_prod%coeffs = 0_8
      carry = 0_8
      DO j = 1, COEFFS_LIMIT - i + 1
          prod_chunk = abs_a%coeffs(i) * abs_b%coeffs(j) + carry
          temp_prod%coeffs(j) = IAND(prod_chunk, MULTI_PRECISION_BASE - 1)
          carry = ISHFT(prod_chunk, -32)
      END DO

      IF (i > 1) THEN
          temp_prod%coeffs(i:COEFFS_LIMIT) = temp_prod%coeffs(1:COEFFS_LIMIT - i + 1)
          temp_prod%coeffs(1:i-1) = 0_8
      END IF

      prod_mpi = prod_mpi + temp_prod
      CALL normalize_mpi(prod_mpi)
  END DO

  IF (final_is_neg) THEN
    prod_mpi%coeffs = -prod_mpi%coeffs
    CALL normalize_mpi(prod_mpi)
  END IF

END FUNCTION mpi_multiply

FUNCTION mpi_to_string(mpi_in) RESULT(str_out)
  TYPE(mpi), INTENT(IN)       :: mpi_in
  CHARACTER(LEN=:), ALLOCATABLE :: str_out

  TYPE(mpi) :: temp_mpi
  INTEGER, PARAMETER :: NUM_DECIMAL_DIGITS = 9
  INTEGER(KIND=8), PARAMETER :: DECIMAL_BASE = 10_8**NUM_DECIMAL_DIGITS
  CHARACTER(LEN=NUM_DECIMAL_DIGITS) :: chunk_str_formatted
  INTEGER(KIND=8) :: digit_chunks(COEFFS_LIMIT * 10), remainder, current_val
  INTEGER :: i, chunk_count, final_str_len, start_pos, current_pos
  LOGICAL :: is_zero
  LOGICAL :: is_negative
  CHARACTER(LEN=:), ALLOCATABLE :: first_chunk_trimmed

  is_negative = mpi_sign(mpi_in)
  IF (mpi_is_zero(mpi_in)) THEN
    str_out = "0"
    RETURN
  END IF

  IF (is_negative) THEN
      temp_mpi = -mpi_in ! Correctly get magnitude of negative number
  ELSE
      temp_mpi = mpi_in
  END IF

  chunk_count = 0
  DO
      IF (mpi_is_zero(temp_mpi)) EXIT
      chunk_count = chunk_count + 1
      digit_chunks(chunk_count) = mpi_div_by_scalar(temp_mpi, DECIMAL_BASE)
  END DO

  WRITE(chunk_str_formatted, '(I0)') digit_chunks(chunk_count)
  first_chunk_trimmed = TRIM(ADJUSTL(chunk_str_formatted))
  final_str_len = LEN(first_chunk_trimmed) + (chunk_count - 1) * NUM_DECIMAL_DIGITS
  start_pos = 1
  IF (is_negative) THEN
    final_str_len = final_str_len + 1
    start_pos = 2
  END IF
  ALLOCATE(CHARACTER(LEN=final_str_len) :: str_out)

  current_pos = final_str_len
  DO i = 1, chunk_count - 1
    WRITE(chunk_str_formatted, '(I9.9)') digit_chunks(i)
    str_out(current_pos - NUM_DECIMAL_DIGITS + 1 : current_pos) = chunk_str_formatted
    current_pos = current_pos - NUM_DECIMAL_DIGITS
  END DO
  str_out(start_pos:current_pos) = first_chunk_trimmed

  IF (is_negative) str_out(1:1) = "-"

END FUNCTION mpi_to_string

SUBROUTINE new_mpi_from_string(str_in, mpi_out)
  CHARACTER(LEN=*), INTENT(IN) :: str_in
  TYPE(mpi), INTENT(OUT)     :: mpi_out
  
  CHARACTER(LEN=:), ALLOCATABLE :: num_str
  INTEGER :: i, j, len_str, numeric_len, first_chunk_len
  INTEGER(KIND=8) :: chunk_val
  INTEGER, PARAMETER :: NUM_DECIMAL_DIGITS = 9
  INTEGER(KIND=8), PARAMETER :: DECIMAL_CHUNK_BASE = 10_8**NUM_DECIMAL_DIGITS
  LOGICAL :: is_negative

  mpi_out%coeffs = 0_8
  num_str = TRIM(ADJUSTL(str_in))
  len_str = LEN(num_str)
  IF (len_str == 0) THEN; RETURN; END IF

  is_negative = .FALSE.
  i = 1
  IF (num_str(1:1) == '-') THEN; is_negative = .TRUE.; i = 2; END IF
  IF (num_str(1:1) == '+') THEN; i = 2; END IF

  DO WHILE (i <= len_str .AND. num_str(i:i) == '0'); i = i + 1; END DO
  IF (i > len_str) THEN; RETURN; END IF

  numeric_len = len_str - i + 1
  first_chunk_len = MOD(numeric_len, NUM_DECIMAL_DIGITS)
  IF (first_chunk_len == 0) THEN
      first_chunk_len = NUM_DECIMAL_DIGITS
  END IF

  chunk_val = 0_8
  DO j = i, i + first_chunk_len - 1
      chunk_val = chunk_val * 10_8 + (ICHAR(num_str(j:j)) - ICHAR('0'))
  END DO
  mpi_out%coeffs(1) = chunk_val
  i = i + first_chunk_len

  DO WHILE (i <= len_str)
      CALL mpi_multiply_by_scalar(mpi_out, DECIMAL_CHUNK_BASE)

      chunk_val = 0_8
      DO j = i, i + NUM_DECIMAL_DIGITS - 1
          chunk_val = chunk_val * 10_8 + (ICHAR(num_str(j:j)) - ICHAR('0'))
      END DO
      mpi_out%coeffs(1) = mpi_out%coeffs(1) + chunk_val
      CALL normalize_mpi(mpi_out)
      i = i + NUM_DECIMAL_DIGITS
  END DO
  IF (is_negative) mpi_out = -mpi_out

END SUBROUTINE new_mpi_from_string

END MODULE multi_precision_integer_mod