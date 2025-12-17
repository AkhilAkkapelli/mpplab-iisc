PROGRAM test_multi_precision_linear_algebra
  USE multi_precision_integer_mod, ONLY: COEFFS_LIMIT, MULTI_PRECISION_BASE, mpi, new_mpi_from_integer, new_mpi_from_coeffs, &
      mpi_is_zero, mpi_sign, mpi_size, mpi_shift_bits_left, mpi_to_string, new_mpi_from_string, mpi_div_rem, &
      mpi_max_value, OPERATOR(==)
  USE multi_precision_float_mod, ONLY: mpf, OPERATOR(+), OPERATOR(/), new_mpf_from_integer, &
      OPERATOR(*), new_mpf_from_mpi_exp, new_mpf_from_string, mpf_value_equal, mpf_abs, mpf_to_string, nearest_mpi_to_mpf, &
      truncate_mpf_to_mpi, OPERATOR(<), mpf_from_real16, mpf_to_real16
  USE multi_precision_linear_algebra_mod
  IMPLICIT NONE
  
  INTEGER :: i, seed
  TYPE(mpf) :: range, mpf_val
  REAL(KIND=16) :: r16

  WRITE(*, '(/A)') "Running tests for MultiPrecisionLinearAlgebra module..."
  WRITE(*, '(/A)') "===================================================="

  ! seed = 4567889
  ! range = new_mpf_from_mpi_exp(new_mpi_from_coeffs((/ MULTI_PRECISION_BASE-1, MULTI_PRECISION_BASE-1, 0_8, 0_8 /)), 0_8)
  ! ! range = new_mpf_from_mpi_exp(new_mpi_from_coeffs((/ MULTI_PRECISION_BASE-1, MULTI_PRECISION_BASE-1, MULTI_PRECISION_BASE-1, MULTI_PRECISION_BASE-1/)), 0_1)
  ! ! print*, range
  ! mpf_val = create_random_mpf(range, seed)
  ! DO i = 1, 10
  !   print('(A,I0,5I20)'), "mpf", i, mpf_val
  !   IF(.NOT. mpf_val < range) print*, "mpf_val out of range"
  ! END DO
  ! print*, "one", mpf_from_real16(1._16)
  ! print*, "two", mpf_from_real16(2._16)
  ! print*, "sixteen", mpf_from_real16(16._16)
  ! print*, "one half", 1.5_16, mpf_from_real16(1.5_16)

  ! print*, "random real"
  ! DO i = 1, 5
  !   print*, "i", i
  !   CALL random_number(r16)
  !   r16 = r16 * MULTI_PRECISION_BASE
  !   print*, "r16", r16
  !   mpf_val = mpf_from_real16(r16)
  !   print*, "mpf_val", mpf_to_string(mpf_val)
  !   r16 = mpf_to_real16(mpf_val)
  !   print*, "r16", r16
  !   print*, ""
  ! END DO
  call test_dot_product()
  ! call test_random_mpf_distribution()
  ! call test_new_mpf_vector_from_mpfs()

  WRITE(*, '(/A)') "===================================================="
  WRITE(*, '(A/)') "All MultiPrecisionLinearAlgebra tests completed"

CONTAINS

  FUNCTION RANDOM(seed) RESULT(u)
    !Uniform Distribution

    INTEGER, INTENT(INOUT) :: seed
    REAL(KIND=8)           :: u

    INTEGER, PARAMETER     :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    REAL(KIND=8),   SAVE          :: am
    INTEGER, SAVE          :: ix=-1,iy=-1,k

    !!!Initialise
    IF (seed <= 0 .OR. iy < 0) THEN
    am=nearest(1.0,-1.0)/IM
    iy=ior(ieor(888889999,abs(seed)),1)
    ix=ieor(777755555,abs(seed))
    seed=abs(seed)+1
    END IF

    !!! Marsaglia shift sequence with period 2^32 -1
    ix=ieor(ix,ishft(ix,13))
    ix=ieor(ix,ishft(ix,-17))
    ix=ieor(ix,ishft(ix,5))

    !!! Park-Miller sequence by Schrage’s method with period2^31−2
    k=iy/IQ
    iy=IA*(iy-k*IQ)-IR*k
    if (iy < 0) iy=iy+IM

    !!! Combine two Generators
    u=am*ior(iand(IM,ieor(ix,iy)),1)

  END FUNCTION RANDOM

  FUNCTION create_random_mpi(range, seed) RESULT(val)
    TYPE(mpi), INTENT(IN)  :: range
    INTEGER, INTENT(INOUT) :: seed
    TYPE(mpi) :: val

    INTEGER(KIND=8) :: sign
    INTEGER :: i

    sign = MERGE(-1_8, 1_8, mpi_sign(range))

    IF(mpi_is_zero(range)) THEN
      val = new_mpi_from_integer(0_8)
      RETURN
    END IF

    DO i = 1, COEFFS_LIMIT
      IF(i==mpi_size(range)) THEN
        val%coeffs(i) = sign * (RANDOM(seed) * range%coeffs(i))
        val%coeffs(i+1:) = 0_8
        EXIT
      ELSE
        val%coeffs(i) = sign * RANDOM(seed) * MULTI_PRECISION_BASE  
      END If
    END DO

  END FUNCTION create_random_mpi

  SUBROUTINE test_random_mpi_distribution()
    INTEGER, PARAMETER :: num_samples = 2**16
    INTEGER, PARAMETER :: num_buckets = 10

    TYPE(mpi) :: range_mpi, random_val, index_int
    TYPE(mpf) :: bucket_width, index
    INTEGER :: buckets(num_buckets), i, bucket_index
    INTEGER :: seed

    WRITE(*, '(A)') "Testing create_random_mpi for uniform distribution..."

    seed = 1234567
    range_mpi = new_mpi_from_coeffs( (/100_8, 0_8, 0_8, 0_8/) )
    bucket_width = new_mpf_from_mpi_exp(range_mpi, 0_8) / new_mpf_from_integer(num_buckets)
    print*, mpf_to_string(bucket_width)

    buckets = 0
    DO i = 1, num_samples
      random_val = create_random_mpi(range_mpi, seed)
      index = new_mpf_from_mpi_exp(random_val, 0_8) / bucket_width
      index_int = truncate_mpf_to_mpi(index)
      bucket_index = INT(index_int%coeffs(1)) + 1
      IF (bucket_index >= 1 .AND. bucket_index <= num_buckets) buckets(bucket_index) = buckets(bucket_index) + 1
    END DO

    DO i = 1, num_buckets
      WRITE(*, '(A,I2,A,I0)') "Bucket ", i, " count: ", buckets(i)
    END DO
    print*, "sum", sum(buckets), num_samples

  END SUBROUTINE test_random_mpi_distribution

  FUNCTION create_random_mpf(range, seed) RESULT(val)
    TYPE(mpf), INTENT(IN)  :: range    
    INTEGER, INTENT(INOUT) :: seed
    TYPE(mpf) :: val
    
    TYPE(mpi) :: random_mpi, max_mpi

    max_mpi = new_mpi_from_integer(MULTI_PRECISION_BASE)
    random_mpi = create_random_mpi(max_mpi, seed)

    val = new_mpf_from_mpi_exp(random_mpi, 0_8) / new_mpf_from_mpi_exp(max_mpi, 0_8)
    val = val * range

  END FUNCTION create_random_mpf

  SUBROUTINE test_random_mpf_distribution()
    INTEGER, PARAMETER :: num_samples = 2**16
    INTEGER, PARAMETER :: num_buckets = 10

    TYPE(mpi) :: index_int
    TYPE(mpf) :: random_val, range_mpf, bucket_width, index
    INTEGER :: buckets(num_buckets), i, bucket_index
    INTEGER :: seed

    WRITE(*, '(A)') "Testing create_random_mpf for uniform distribution..."

    seed = 1234567
    range_mpf = new_mpf_from_integer(1234567)
    bucket_width = range_mpf / new_mpf_from_integer(num_buckets)

    buckets = 0
    DO i = 1, num_samples
      random_val = create_random_mpf(range_mpf, seed)
      index = random_val / bucket_width
      index_int = truncate_mpf_to_mpi(index)
      bucket_index = INT(index_int%coeffs(1)) + 1
      IF (bucket_index >= 1 .AND. bucket_index <= num_buckets) buckets(bucket_index) = buckets(bucket_index) + 1
    END DO

    DO i = 1, num_buckets
      WRITE(*, '(A,I2,A,I0)') "Bucket ", i, " count: ", buckets(i)
    END DO
    print*, "sum", sum(buckets), num_samples

  END SUBROUTINE test_random_mpf_distribution

  SUBROUTINE test_nearest_mpi_to_mpf()
    TYPE(mpf) :: mpf_val
    TYPE(mpi) :: mpi_res, mpi_expected

    WRITE(*, '(A)') "Testing nearest_mpi_to_mpf..."

    mpf_val = new_mpf_from_string("-3438729169.465782")
    print*, "DEBUG", mpf_val
    mpi_res = nearest_mpi_to_mpf(mpf_val)
    mpi_expected = new_mpi_from_integer(-3438729169_8)
    print*, "DEBUG", mpi_to_string(mpi_res)
    CALL assert(mpi_res == mpi_expected, "Rounding 3.4 to 3")

  END SUBROUTINE test_nearest_mpi_to_mpf

  SUBROUTINE test_truncate_mpf_to_mpi()
    TYPE(mpf) :: mpf_val
    TYPE(mpi) :: mpi_res, mpi_expected

    WRITE(*, '(A)') "Testing truncate_mpf_to_mpi..."

    mpf_val = new_mpf_from_string("-3438729169.865782")
    print*, "DEBUG", mpf_val
    mpi_res = truncate_mpf_to_mpi(mpf_val)
    mpi_expected = new_mpi_from_integer(-3438729169_8)
    print*, "DEBUG", mpi_to_string(mpi_res)
    CALL assert(mpi_res == mpi_expected, "Truncating 3.7 to 3")

  END SUBROUTINE test_truncate_mpf_to_mpi

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

  ! SUBROUTINE test_mpi_shift_bits_left()
  !   TYPE(mpi) :: mpi_in, mpi_out
  !   TYPE(mpf) :: mpf_out, mpf_expected
  !   INTEGER :: num_bits
  !   INTEGER :: seed

  !   WRITE(*, '(A)') "Testing mpi_shift_bits_left..."

  !   seed = 1234567
  !   mpi_in = create_random_mpi(seed)
  !   num_bits = 0
  !   mpi_out = mpi_shift_bits_left(mpi_in, num_bits)
  !   CALL assert(mpi_out == mpi_in, "Shift by 0 bits")

  !   mpi_in = new_mpi_from_integer(7_8)
  !   num_bits = 10
  !   mpf_out = new_mpf_from_mpi_exp(mpi_shift_bits_left(mpi_in, num_bits), 0)
  !   mpf_expected = new_mpf_from_mpi_exp(mpi_in, num_bits)
  !   CALL assert(mpf_value_equal(mpf_out, mpf_expected), "Small positive shift")

  !   mpi_in = new_mpi_from_integer(13_8)
  !   num_bits = 70 
  !   mpf_out = new_mpf_from_mpi_exp(mpi_shift_bits_left(mpi_in, num_bits), 0)
  !   mpf_expected = new_mpf_from_mpi_exp(mpi_in, num_bits)
  !   CALL assert(mpf_value_equal(mpf_out, mpf_expected), "Large positive shift (cross-boundary)")

  !   mpi_in = new_mpi_from_integer(-19_8)
  !   num_bits = 40
  !   mpf_out = new_mpf_from_mpi_exp(mpi_shift_bits_left(mpi_in, num_bits), 0)
  !   mpf_expected = new_mpf_from_mpi_exp(mpi_in, num_bits)
  !   CALL assert(mpf_value_equal(mpf_out, mpf_expected), "Shift negative number")

  !   ! Test 5: Overflow shift (should result in zero)
  !   mpi_in = new_mpi_from_integer(1_8)
  !   num_bits = COEFFS_LIMIT * 32 ! A shift large enough to clear all bits
  !   mpi_out = mpi_shift_bits_left(mpi_in, num_bits)
  !   CALL assert(mpi_is_zero(mpi_out), "Overflow shift results in zero")

  !   ! Test 6: Negative shift (should do nothing)
  !   mpi_in = new_mpi_from_integer(98765_8)
  !   mpi_out = mpi_shift_bits_left(mpi_in, -10)
  !   CALL assert(mpi_out == mpi_in, "Shift by negative bits")

  ! END SUBROUTINE test_mpi_shift_bits_left

  SUBROUTINE test_new_mpf_vector_from_mpfs()
    TYPE(mpf), DIMENSION(:), ALLOCATABLE :: mpf_elements
    TYPE(mpf) :: complex_mpf(4)
    TYPE(mpf_vector) :: vec
    TYPE(mpi) :: mpi_val_1, mpi_val_2, mpi_val_3
    INTEGER :: seed

    TYPE(mpf) :: range_mpf
    
    WRITE(*, '(A)') "Testing new_mpf_vector_from_mpfs..."

    ALLOCATE(mpf_elements(0))
    vec = new_mpf_vector_from_mpfs(mpf_elements)
    CALL assert(vec%n == 0, "Vector Construction (Empty): Size is 0")
    CALL assert(.NOT. ALLOCATED(vec%mantissas), "Vector Construction (Empty): Mantissas not allocated")
    DEALLOCATE(mpf_elements)

    ALLOCATE(mpf_elements(3))
    mpf_elements(:) = new_mpf_from_string("0.0")
    vec = new_mpf_vector_from_mpfs(mpf_elements)
    CALL assert(vec%n == 3, "Vector Construction (Zeros): Size")
    CALL assert(vec%mantissas(1) == new_mpi_from_integer(0_8) .AND. &
                vec%mantissas(2) == new_mpi_from_integer(0_8) .AND. &
                vec%mantissas(3) == new_mpi_from_integer(0_8), &
                "Vector Construction (Zeros): All mantissas are zero")
    DEALLOCATE(mpf_elements)

    seed = 1234567
    range_mpf = new_mpf_from_mpi_exp(new_mpi_from_coeffs((/ MULTI_PRECISION_BASE-1, MULTI_PRECISION_BASE-1, MULTI_PRECISION_BASE-1, 0_8 /)), 0_8)

    complex_mpf(1) = create_random_mpf(range_mpf, seed)
    complex_mpf(2) = create_random_mpf(new_mpf_from_integer(MULTI_PRECISION_BASE), seed)
    complex_mpf(3) = create_random_mpf(new_mpf_from_integer(-MULTI_PRECISION_BASE), seed)
    complex_mpf(4) = create_random_mpf(-range_mpf, seed)

    vec = new_mpf_vector_from_mpfs(complex_mpf)

    print*, vec%exponent
    DO i = 1, vec%n
      print('(5I20)'), vec%mantissas(i)
    END DO

    CALL assert(vec%n == 4, "Vector Construction (Complex): Size")


    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec%mantissas(1), vec%exponent), complex_mpf(1)), &
                "Vector Construction (Complex): Value 1 preserved")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec%mantissas(2), vec%exponent), complex_mpf(2)), &
                "Vector Construction (Complex): Value 2 preserved")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec%mantissas(3), vec%exponent), complex_mpf(3)), &
                "Vector Construction (Complex): Value 3 preserved")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec%mantissas(4), vec%exponent), complex_mpf(4)), &
                "Vector Construction (Complex): Value 4 preserved")

    WRITE(*, '(A)') "Testing new_mpf_vector_from_mpfs completed successfully"

  END SUBROUTINE test_new_mpf_vector_from_mpfs

  SUBROUTINE test_vector_arithmetic()
    TYPE(mpf_vector) :: vec1, vec2, vec_res
    TYPE(mpf) :: dot_res, expected_res, scalar

    WRITE(*, '(A)') "Testing vector arithmetic..."

    vec1 = new_mpf_vector_from_mpfs([new_mpf_from_string("1.5"), new_mpf_from_string("2.0")])
    vec2 = new_mpf_vector_from_mpfs([new_mpf_from_string("3.0"), new_mpf_from_string("4.25")])

    PRINT *, "--- Inside test_vector_arithmetic ---"
    
    ! Dot Product
    dot_res = vec1 * vec2
    expected_res = new_mpf_from_string("13.0")
    CALL assert(mpf_value_equal(dot_res, expected_res), "Vector Dot Product")

    ! Addition
    vec_res = vec1 + vec2
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(1), vec_res%exponent), new_mpf_from_string("4.5")), "Vector Addition: Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(2), vec_res%exponent), new_mpf_from_string("6.25")), "Vector Addition: Element 2")

    ! Scalar Multiplication
    scalar = new_mpf_from_string("2.0")
    vec_res = vec1 * scalar
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(1), vec_res%exponent), new_mpf_from_string("3.0")), "Vector-Scalar Mul: Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(2), vec_res%exponent), new_mpf_from_string("4.0")), "Vector-Scalar Mul: Element 2")

  END SUBROUTINE test_vector_arithmetic

  SUBROUTINE test_matrix_arithmetic()
      TYPE(mpf_matrix) :: mat1, mat2, mat_res
      TYPE(mpf_vector) :: vec, vec_res

      WRITE(*, '(A)') "Testing matrix arithmetic..."

      mat1 = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("1.0"), new_mpf_from_string("3.0"), &
                                               new_mpf_from_string("2.0"), new_mpf_from_string("4.0")], [2,2]))
      vec = new_mpf_vector_from_mpfs([new_mpf_from_string("5.0"), new_mpf_from_string("6.0")])

      ! Matrix-Vector product
      vec_res = mat1 * vec
      ! Expected: [1*5+2*6, 3*5+4*6] = [17, 39]
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(1), vec_res%exponent), new_mpf_from_string("17.0")), "Matrix-Vector Mul: Element 1")
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec_res%mantissas(2), vec_res%exponent), new_mpf_from_string("39.0")), "Matrix-Vector Mul: Element 2")

      ! Matrix-Matrix product
      mat2 = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("1.0"), new_mpf_from_string("0.0"), &
                                               new_mpf_from_string("0.0"), new_mpf_from_string("1.0")], [2,2]))
      mat_res = mat1 * mat2
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_res%mantissas(1,1), mat_res%exponent), new_mpf_from_string("1.0")), "Matrix-Matrix Mul: Identity (1,1)")
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_res%mantissas(1,2), mat_res%exponent), new_mpf_from_string("2.0")), "Matrix-Matrix Mul: Identity (1,2)")

      ! Matrix Addition
      mat_res = mat1 + mat1
      CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_res%mantissas(1,1), mat_res%exponent), new_mpf_from_string("2.0")), "Matrix Addition")

  END SUBROUTINE test_matrix_arithmetic

  SUBROUTINE test_level1_blas()
    TYPE(mpf_vector) :: vec
    TYPE(mpf) :: asum_res, nrm2_res, expected_res
    INTEGER :: idx

    WRITE(*, '(A)') "Testing Level 1 BLAS routines..."

    vec = new_mpf_vector_from_mpfs([new_mpf_from_string("-3.0"), new_mpf_from_string("4.0")])

    ! ASUM
    asum_res = mpf_vector_asum(vec)
    expected_res = new_mpf_from_string("7.0")
    CALL assert(mpf_value_equal(asum_res, expected_res), "ASUM: |-3| + |4| = 7")
    
    ! IAMAX / IAMIN
    vec = new_mpf_vector_from_mpfs([new_mpf_from_string("1.0"), new_mpf_from_string("-5.0"), &
                                    new_mpf_from_string("2.0")])
    idx = mpf_vector_iamax(vec)
    CALL assert(idx == 2, "IAMAX: Index of max absolute value is 2")
    idx = mpf_vector_iamin(vec)
    CALL assert(idx == 1, "IAMIN: Index of min absolute value is 1")

  END SUBROUTINE test_level1_blas

  SUBROUTINE test_level2_blas()
    TYPE(mpf_matrix) :: mat
    TYPE(mpf_vector) :: vec, vec_res, b_vec

    WRITE(*, '(A)') "Testing Level 2 BLAS routines..."

    ! --- TRMV Test ---
    ! A = [1 2]  x = [3]  A*x = [1*3+2*4] = [11]
    !     [0 4]      [4]        [0*3+4*4]   [16]
    mat = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("1.0"), new_mpf_from_string("0.0"), &
                                            new_mpf_from_string("2.0"), new_mpf_from_string("4.0")], [2,2]))
    vec = new_mpf_vector_from_mpfs([new_mpf_from_string("3.0"), new_mpf_from_string("4.0")])
    CALL mpf_trmv(mat, vec, 'U')
    vec_res = new_mpf_vector_from_mpfs([new_mpf_from_string("11.0"), new_mpf_from_string("16.0")])
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec%mantissas(1), vec%exponent), new_mpf_from_string("11.0")), "TRMV (Upper): Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(vec%mantissas(2), vec%exponent), new_mpf_from_string("16.0")), "TRMV (Upper): Element 2")

    ! --- TRSV Test ---
    ! A*x = b.  x = inv(A)*b.  Let's solve A*x = [11, 16]
    ! We should get back x = [3, 4]
    b_vec = new_mpf_vector_from_mpfs([new_mpf_from_string("11.0"), new_mpf_from_string("16.0")])
    CALL mpf_trsv(mat, b_vec, 'U')
    vec_res = new_mpf_vector_from_mpfs([new_mpf_from_string("3.0"), new_mpf_from_string("4.0")])
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(b_vec%mantissas(1), b_vec%exponent), new_mpf_from_string("3.0")), "TRSV (Upper): Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(b_vec%mantissas(2), b_vec%exponent), new_mpf_from_string("4.0")), "TRSV (Upper): Element 2")

    ! Lower triangular
    mat = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("3.0"), new_mpf_from_string("1.0"), &
                                            new_mpf_from_string("0.0"), new_mpf_from_string("2.0")], [2,2]))
    b_vec = new_mpf_vector_from_mpfs([new_mpf_from_string("6.0"), new_mpf_from_string("8.0")])
    ! 3x1 = 6 -> x1=2.  1*x1 + 2*x2 = 8 -> 1*2 + 2*x2 = 8 -> 2*x2=6 -> x2=3
    CALL mpf_trsv(mat, b_vec, 'L')
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(b_vec%mantissas(1), b_vec%exponent), new_mpf_from_string("2.0")), "TRSV (Lower): Element 1")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(b_vec%mantissas(2), b_vec%exponent), new_mpf_from_string("3.0")), "TRSV (Lower): Element 2")

  END SUBROUTINE test_level2_blas

  SUBROUTINE test_level3_blas()
    TYPE(mpf_matrix) :: mat_a, mat_b, mat_res

    WRITE(*, '(A)') "Testing Level 3 BLAS routines..."

    ! --- TRMM Test (Upper) ---
    ! A = [1 2], B = [3 4], A*B = [1*3+2*5, 1*4+2*6] = [13, 16]
    !     [0 4]      [5 6]         [0*3+4*5, 0*4+4*6]   [20, 24]
    mat_a = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("1.0"), new_mpf_from_string("0.0"), &
                                              new_mpf_from_string("2.0"), new_mpf_from_string("4.0")], [2,2]))
    mat_b = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("3.0"), new_mpf_from_string("5.0"), &
                                              new_mpf_from_string("4.0"), new_mpf_from_string("6.0")], [2,2]))
    CALL mpf_trmm(mat_a, mat_b, 'U')
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(1,1), mat_b%exponent), new_mpf_from_string("13.0")), "TRMM (Upper): C(1,1)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(2,1), mat_b%exponent), new_mpf_from_string("20.0")), "TRMM (Upper): C(2,1)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(1,2), mat_b%exponent), new_mpf_from_string("16.0")), "TRMM (Upper): C(1,2)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(2,2), mat_b%exponent), new_mpf_from_string("24.0")), "TRMM (Upper): C(2,2)")

    ! --- TRSM Test (Lower) ---
    ! A = [2 0], X = [1 2], B = A*X = [2  4]
    !     [1 3]      [3 4]           [10 14]
    ! We solve A*X=B and expect to get X back.
    mat_a = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("2.0"), new_mpf_from_string("1.0"), &
                                              new_mpf_from_string("0.0"), new_mpf_from_string("3.0")], [2,2]))
    mat_b = new_mpf_matrix_from_mpfs(RESHAPE([new_mpf_from_string("2.0"), new_mpf_from_string("10.0"), &
                                              new_mpf_from_string("4.0"), new_mpf_from_string("14.0")], [2,2]))
    CALL mpf_trsm(mat_a, mat_b, 'L')
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(1,1), mat_b%exponent), new_mpf_from_string("1.0")), "TRSM (Lower): X(1,1)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(2,1), mat_b%exponent), new_mpf_from_string("3.0")), "TRSM (Lower): X(2,1)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(1,2), mat_b%exponent), new_mpf_from_string("2.0")), "TRSM (Lower): X(1,2)")
    CALL assert(mpf_value_equal(new_mpf_from_mpi_exp(mat_b%mantissas(2,2), mat_b%exponent), new_mpf_from_string("4.0")), "TRSM (Lower): X(2,2)")

  END SUBROUTINE test_level3_blas

  SUBROUTINE test_dot_product()
    INTEGER, PARAMETER :: n = 10

    TYPE(mpf_vector) :: vec1_mpf, vec2_mpf
    TYPE(mpf) :: dot_res_mpf
    TYPE(mpf), ALLOCATABLE :: mpf_elements1(:), mpf_elements2(:)

    REAL(16), ALLOCATABLE :: vec1_q(:), vec2_q(:)
    REAL(16) :: dot_res_q, rel_error

    REAL :: rand_val
    CHARACTER(LEN=50) :: str_val
    INTEGER :: i

    WRITE(*, '(A, I0, A)') "Testing large dot product (n=", n, ")..."

    ALLOCATE(vec1_q(n), vec2_q(n))
    ALLOCATE(mpf_elements1(n), mpf_elements2(n))

    CALL RANDOM_SEED()

    ! 1. Generate random vectors for both mpf and REAL(16)
    DO i = 1, n
      CALL RANDOM_NUMBER(rand_val)
      vec1_q(i) = (REAL(rand_val, 16) - 0.5_16) * 1000.0_16
      WRITE(str_val, '(F0.40)') vec1_q(i)
      mpf_elements1(i) = new_mpf_from_string(TRIM(str_val))
      WRITE(*, '(A,I2,A,F45.35)') "DEBUG: vec1_q(", i, ") = ", vec1_q(i)
      WRITE(*, '(A,I2,A,A)')     "DEBUG: vec1_mpf(", i, ") = ", mpf_to_string(mpf_elements1(i))

      CALL RANDOM_NUMBER(rand_val)
      vec2_q(i) = (REAL(rand_val, 16) - 0.5_16) * 200.0_16
      WRITE(str_val, '(F0.40)') vec2_q(i)
      mpf_elements2(i) = new_mpf_from_string(TRIM(str_val))
      WRITE(*, '(A,I2,A,F45.35)') "DEBUG: vec2_q(", i, ") = ", vec2_q(i)
      WRITE(*, '(A,I2,A,A)')     "DEBUG: vec2_mpf(", i, ") = ", mpf_to_string(mpf_elements2(i))
    END DO

    ! 2. Compute dot products
    vec1_mpf = new_mpf_vector_from_mpfs(mpf_elements1)
    vec2_mpf = new_mpf_vector_from_mpfs(mpf_elements2)
    dot_res_mpf = vec1_mpf * vec2_mpf
    dot_res_q = DOT_PRODUCT(vec1_q, vec2_q)

    ! 3. Compare the results by converting the REAL(16) result to an mpf
    WRITE(str_val, '(F0.40)') dot_res_q

    WRITE(*, '(A, A)') "  MPF Result:    ", mpf_to_string(dot_res_mpf)
    WRITE(*, '(A, F36.15)') "  REAL(16) Result: ", dot_res_q

    CALL assert(mpf_value_equal(dot_res_mpf, new_mpf_from_string(TRIM(str_val))), "Large Dot Product: mpf vs REAL(16)")

  END SUBROUTINE test_dot_product

END PROGRAM test_multi_precision_linear_algebra