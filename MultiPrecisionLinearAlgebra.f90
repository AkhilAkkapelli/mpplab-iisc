MODULE multi_precision_linear_algebra_mod
  USE multi_precision_integer_mod, ONLY: mpi, COEFFS_LIMIT, mpi_shift_bits_right, mpi_shift_bits_left, mpi_abs, normalize_mpi, &
                                     mpi_is_zero, new_mpi_from_integer, OPERATOR(*), OPERATOR(+), OPERATOR(-), OPERATOR(<)
  USE multi_precision_float_mod, ONLY: mpf, new_mpf_from_integer, new_mpf_from_mpi_exp, &
                                     normalize_mpf_float, mpf_abs, mpf_is_zero, OPERATOR(+), OPERATOR(*), OPERATOR(/), OPERATOR(-)

  IMPLICIT NONE
  PRIVATE


  TYPE :: mpf_vector
    TYPE(mpi), ALLOCATABLE :: mantissas(:)
    INTEGER :: exponent = 0
    INTEGER :: n = 0
  END TYPE mpf_vector

  TYPE :: mpf_matrix
    TYPE(mpi), ALLOCATABLE :: mantissas(:,:)
    INTEGER :: exponent = 0
    INTEGER :: m = 0, n = 0
  END TYPE mpf_matrix

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE mpf_vector_add
    MODULE PROCEDURE mpf_matrix_add
  END INTERFACE

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE mpf_vector_unary_negate
    MODULE PROCEDURE mpf_matrix_unary_negate
    MODULE PROCEDURE mpf_vector_subtract
    MODULE PROCEDURE mpf_matrix_subtract
  END INTERFACE

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE mpf_vector_dot_product
    MODULE PROCEDURE mpf_vector_scalar_multiply
    MODULE PROCEDURE mpf_matrix_scalar_multiply
    MODULE PROCEDURE mpf_matrix_vector_multiply
    MODULE PROCEDURE mpf_matrix_multiply
  END INTERFACE

  PUBLIC :: mpf_vector, mpf_matrix, new_mpf_vector_from_mpfs, new_mpf_matrix_from_mpfs, &
            OPERATOR(+), OPERATOR(-), OPERATOR(*), normalize_mpf_vector, normalize_mpf_matrix, &
            mpf_vector_asum, &
            mpf_vector_iamax, mpf_vector_iamin, &
            mpf_trmv, mpf_trsv, &
            mpf_trmm, mpf_trsm

CONTAINS

  FUNCTION new_mpf_vector_from_mpfs(mpf_array) RESULT(vec)
    TYPE(mpf), INTENT(IN) :: mpf_array(:)
    TYPE(mpf_vector) :: vec

    INTEGER :: i
    INTEGER :: shift

    vec%n = SIZE(mpf_array)
    IF (vec%n == 0) RETURN

    ALLOCATE(vec%mantissas(vec%n))

    IF (ALL(mpf_is_zero(mpf_array))) THEN
      vec%mantissas = new_mpi_from_integer(0_8)
      vec%exponent = 0
      RETURN
    END IF

    vec%exponent = MINVAL(mpf_array%exponent, MASK=.NOT. mpf_is_zero(mpf_array))

    DO i = 1, vec%n
      shift = mpf_array(i)%exponent - vec%exponent
      vec%mantissas(i) = mpi_shift_bits_left(mpf_array(i)%mantissa, shift)
    END DO

  END FUNCTION new_mpf_vector_from_mpfs

  FUNCTION new_mpf_matrix_from_mpfs(mpf_array) RESULT(mat)
    TYPE(mpf), INTENT(IN) :: mpf_array(:,:)
    TYPE(mpf_matrix) :: mat

    INTEGER :: i, j

    mat%m = SIZE(mpf_array, 1)
    mat%n = SIZE(mpf_array, 2)
    IF (mat%m == 0 .OR. mat%n == 0) RETURN

    ALLOCATE(mat%mantissas(mat%m, mat%n))

    mat%exponent = MINVAL(mpf_array%exponent, MASK=.NOT. mpf_is_zero(mpf_array))

    DO j = 1, mat%n
      DO i = 1, mat%m
        mat%mantissas(i,j) = mpi_shift_bits_left(mpf_array(i,j)%mantissa, mpf_array(i,j)%exponent - mat%exponent)
      END DO
    END DO

  END FUNCTION new_mpf_matrix_from_mpfs

  FUNCTION mpf_vector_dot_product(vec1, vec2) RESULT(dot_prod)
    TYPE(mpf_vector), INTENT(IN) :: vec1, vec2
    TYPE(mpf) :: dot_prod

    TYPE(mpi) :: total_sum
    INTEGER :: n, i, j, k
    INTEGER(KIND=8), ALLOCATABLE :: slice1(:), slice2(:)

    n = vec1%n
    IF (vec2%n /= n) THEN
      STOP "mpf_vector_dot_product: Vector dimensions do not match."
    END IF

    IF (n == 0) THEN
      dot_prod = new_mpf_from_integer(0_8)
      RETURN
    END IF

    dot_prod%exponent = vec1%exponent + vec2%exponent

    DO i = 1, COEFFS_LIMIT
      dot_prod%mantissa%coeffs(i) = DOT_PRODUCT(vec1%mantissas(:)%coeffs(i), vec2%mantissas(:)%coeffs(i))
    END DO

    CALL normalize_mpf_float(dot_prod)

  END FUNCTION mpf_vector_dot_product

  FUNCTION mpf_vector_add(vec1, vec2) RESULT(res)
    TYPE(mpf_vector), INTENT(IN) :: vec1, vec2
    TYPE(mpf_vector) :: res
    INTEGER :: i, common_exp

    IF (vec1%n /= vec2%n) STOP "mpf_vector_add: Vector dimensions do not match."
    res%n = vec1%n
    IF (res%n == 0) RETURN

    ALLOCATE(res%mantissas(res%n))

    IF (vec1%exponent <= vec2%exponent) THEN
      res%exponent = vec1%exponent
      DO i = 1, res%n
        res%mantissas(i) = vec1%mantissas(i) + mpi_shift_bits_left(vec2%mantissas(i), vec2%exponent - res%exponent)
      END DO
    ELSE
      res%exponent = vec2%exponent
      DO i = 1, res%n
        res%mantissas(i) = vec2%mantissas(i) + mpi_shift_bits_left(vec1%mantissas(i), vec1%exponent - res%exponent)
      END DO
    END IF

    CALL normalize_mpf_vector(res)

  END FUNCTION mpf_vector_add

  FUNCTION mpf_vector_subtract(vec1, vec2) RESULT(res)
    TYPE(mpf_vector), INTENT(IN) :: vec1, vec2
    TYPE(mpf_vector) :: res
    res = vec1 + (-vec2)
  END FUNCTION mpf_vector_subtract

  FUNCTION mpf_vector_scalar_multiply(vec, scalar) RESULT(res)
    TYPE(mpf_vector), INTENT(IN) :: vec
    TYPE(mpf), INTENT(IN) :: scalar  
    TYPE(mpf_vector) :: res
    
    INTEGER :: i
    res%n = vec%n
    IF (vec%n == 0) RETURN

    ALLOCATE(res%mantissas(res%n))
    res%exponent = vec%exponent + scalar%exponent

    DO i = 1, COEFFS_LIMIT
      res%mantissas(:)%coeffs(i) = vec%mantissas(:)%coeffs(i) * scalar%mantissa%coeffs(i)
    END DO

    CALL normalize_mpf_vector(res)

  END FUNCTION mpf_vector_scalar_multiply

  FUNCTION mpf_matrix_scalar_multiply(mat, scalar) RESULT(res)
    TYPE(mpf_matrix), INTENT(IN) :: mat
    TYPE(mpf), INTENT(IN) :: scalar
    TYPE(mpf_matrix) :: res

    INTEGER :: i
    res%m = mat%m
    res%n = mat%n
    IF (mat%m == 0 .OR. mat%n == 0) RETURN

    ALLOCATE(res%mantissas(res%m, res%n))
    res%exponent = mat%exponent + scalar%exponent

    DO i = 1, COEFFS_LIMIT
      res%mantissas(:,:)%coeffs(i) = mat%mantissas(:,:)%coeffs(i) * scalar%mantissa%coeffs(i)
    END DO

    CALL normalize_mpf_matrix(res)

  END FUNCTION mpf_matrix_scalar_multiply

  FUNCTION mpf_matrix_add(mat1, mat2) RESULT(res)
    TYPE(mpf_matrix), INTENT(IN) :: mat1, mat2
    TYPE(mpf_matrix) :: res
    INTEGER :: i, j, common_exp

    IF (mat1%m /= mat2%m .OR. mat1%n /= mat2%n) STOP "mpf_matrix_add: Matrix dimensions do not match."
    res%m = mat1%m
    res%n = mat1%n
    IF (res%m == 0 .OR. res%n == 0) RETURN

    ALLOCATE(res%mantissas(res%m, res%n))

    IF (mat1%exponent <= mat2%exponent) THEN
      res%exponent = mat1%exponent
      DO j = 1, res%n
        DO i = 1, res%m
          res%mantissas(i,j) = mat1%mantissas(i,j) + mpi_shift_bits_left(mat2%mantissas(i,j), mat2%exponent - res%exponent)
        END DO
      END DO
    ELSE
      res%exponent = mat2%exponent
      DO j = 1, res%n
        DO i = 1, res%m
          res%mantissas(i,j) = mat2%mantissas(i,j) + mpi_shift_bits_left(mat1%mantissas(i,j), mat1%exponent - res%exponent)
        END DO
      END DO
    END IF

    CALL normalize_mpf_matrix(res)
    
  END FUNCTION mpf_matrix_add

  FUNCTION mpf_matrix_subtract(mat1, mat2) RESULT(res)
    TYPE(mpf_matrix), INTENT(IN) :: mat1, mat2
    TYPE(mpf_matrix) :: res
    res = mat1 + (-mat2)
  END FUNCTION mpf_matrix_subtract

  FUNCTION mpf_vector_unary_negate(vec) RESULT(res)
    TYPE(mpf_vector), INTENT(IN) :: vec
    TYPE(mpf_vector) :: res
    INTEGER :: i

    res = vec
    ALLOCATE(res%mantissas(vec%n))
    
    DO i = 1, res%n
      res%mantissas(i) = -vec%mantissas(i)
    END DO

  END FUNCTION mpf_vector_unary_negate

  FUNCTION mpf_matrix_unary_negate(mat) RESULT(res)
    TYPE(mpf_matrix), INTENT(IN) :: mat
    TYPE(mpf_matrix) :: res
    INTEGER :: i, j

    res = mat
    ALLOCATE(res%mantissas(mat%m, mat%n))
    DO j = 1, res%n
      DO i = 1, res%m
        res%mantissas(i,j) = -mat%mantissas(i,j)
      END DO
    END DO

  END FUNCTION mpf_matrix_unary_negate

  FUNCTION mpf_matrix_vector_multiply(mat, vec) RESULT(res)
    TYPE(mpf_matrix), INTENT(IN) :: mat
    TYPE(mpf_vector), INTENT(IN) :: vec
    TYPE(mpf_vector) :: res
    INTEGER :: i, j, k, l

    IF (mat%n /= vec%n) STOP "mpf_matrix_vector_multiply: Inner dimensions do not match."
    res%n = mat%m
    IF (res%n == 0) RETURN

    ALLOCATE(res%mantissas(res%n))
    res%mantissas = new_mpi_from_integer(0_8)
    res%exponent = mat%exponent + vec%exponent

    DO j = 1, COEFFS_LIMIT
      DO i = 1, mat%m
          res%mantissas(i)%coeffs(j) = res%mantissas(i)%coeffs(j) + DOT_PRODUCT(mat%mantissas(i,:)%coeffs(j), vec%mantissas(:)%coeffs(j))
        END DO
      END DO

    CALL normalize_mpf_vector(res)
    
  END FUNCTION mpf_matrix_vector_multiply

  FUNCTION mpf_matrix_multiply(mat1, mat2) RESULT(res_mat)
    TYPE(mpf_matrix), INTENT(IN) :: mat1, mat2
    TYPE(mpf_matrix) :: res_mat
    INTEGER :: i, j, k, l, p

    IF (mat1%n /= mat2%m) STOP "mpf_matrix_multiply: Inner dimensions do not match."
    res_mat%m = mat1%m
    res_mat%n = mat2%n
    IF (res_mat%m == 0 .OR. res_mat%n == 0) RETURN

    ALLOCATE(res_mat%mantissas(res_mat%m, res_mat%n))
    res_mat%mantissas = new_mpi_from_integer(0_8)
    res_mat%exponent = mat1%exponent + mat2%exponent

    DO k = 1, COEFFS_LIMIT
      DO l = 1, COEFFS_LIMIT
        IF (k + l - 1 > COEFFS_LIMIT) CYCLE
        DO j = 1, res_mat%n
          DO i = 1, res_mat%m
            res_mat%mantissas(i,j)%coeffs(k + l - 1) = res_mat%mantissas(i,j)%coeffs(k + l - 1) + &
                DOT_PRODUCT(mat1%mantissas(i,:)%coeffs(k), mat2%mantissas(:,j)%coeffs(l))
          END DO
        END DO
      END DO
    END DO

   CALL normalize_mpf_matrix(res_mat)
    
  END FUNCTION mpf_matrix_multiply

  FUNCTION mpf_vector_asum(vec) RESULT(res)
    TYPE(mpf_vector), INTENT(IN) :: vec
    TYPE(mpf) :: res
    TYPE(mpi) :: sum_mantissa, temp_abs
    INTEGER :: i

    sum_mantissa%coeffs = 0_8
    DO i = 1, vec%n
      temp_abs = mpi_abs(vec%mantissas(i))
      sum_mantissa%coeffs = sum_mantissa%coeffs + temp_abs%coeffs
    END DO
    CALL normalize_mpi(sum_mantissa)
    res = new_mpf_from_mpi_exp(sum_mantissa, vec%exponent)
    
  END FUNCTION mpf_vector_asum

  FUNCTION mpf_vector_iamax(vec) RESULT(idx)
    TYPE(mpf_vector), INTENT(IN) :: vec
    INTEGER :: idx
    INTEGER :: i, k
    LOGICAL :: current_is_larger

    IF (vec%n < 1) THEN; idx = 0; RETURN; END IF

    idx = 1
    DO i = 2, vec%n
      current_is_larger = .FALSE.
      DO k = COEFFS_LIMIT, 1, -1
        IF (ABS(vec%mantissas(i)%coeffs(k)) > ABS(vec%mantissas(idx)%coeffs(k))) THEN
          current_is_larger = .TRUE.
          EXIT
        ELSE IF (ABS(vec%mantissas(i)%coeffs(k)) < ABS(vec%mantissas(idx)%coeffs(k))) THEN
          EXIT
        END IF
      END DO

      IF (current_is_larger) THEN
        idx = i
      END IF
    END DO
  END FUNCTION mpf_vector_iamax

  FUNCTION mpf_vector_iamin(vec) RESULT(idx)
    TYPE(mpf_vector), INTENT(IN) :: vec
    INTEGER :: idx
    INTEGER :: i, k
    LOGICAL :: current_is_smaller

    IF (vec%n < 1) THEN; idx = 0; RETURN; END IF

    idx = 1
    DO i = 2, vec%n
      current_is_smaller = .FALSE.
      DO k = COEFFS_LIMIT, 1, -1
        IF (ABS(vec%mantissas(i)%coeffs(k)) < ABS(vec%mantissas(idx)%coeffs(k))) THEN
          current_is_smaller = .TRUE.
          EXIT
        ELSE IF (ABS(vec%mantissas(i)%coeffs(k)) > ABS(vec%mantissas(idx)%coeffs(k))) THEN
          EXIT
        END IF
      END DO

      IF (current_is_smaller) THEN
        idx = i 
      END IF
    END DO
  END FUNCTION mpf_vector_iamin

  SUBROUTINE mpf_trmv(mat, vec, uplo)
    TYPE(mpf_matrix), INTENT(IN) :: mat
    TYPE(mpf_vector), INTENT(INOUT) :: vec
    CHARACTER(LEN=1), INTENT(IN) :: uplo ! 'U' for Upper, 'L' for Lower
    TYPE(mpf_vector) :: res_vec
    TYPE(mpi) :: row_sum
    INTEGER :: i, j

    IF (mat%m /= mat%n .OR. mat%n /= vec%n) STOP "mpf_trmv: Matrix must be square and dims must match."
    res_vec%n = vec%n
    ALLOCATE(res_vec%mantissas(res_vec%n))
    res_vec%exponent = mat%exponent + vec%exponent

    SELECT CASE (uplo)
    CASE ('U', 'u')
      DO i = 1, mat%m
        row_sum = new_mpi_from_integer(0_8)
        DO j = i, mat%n
          row_sum = row_sum + (mat%mantissas(i,j) * vec%mantissas(j))
        END DO
        res_vec%mantissas(i) = row_sum
      END DO
    CASE ('L', 'l')
      DO i = 1, mat%m
        row_sum = new_mpi_from_integer(0_8)
        DO j = 1, i
          row_sum = row_sum + (mat%mantissas(i,j) * vec%mantissas(j))
        END DO
        res_vec%mantissas(i) = row_sum
      END DO
    CASE DEFAULT
      STOP "mpf_trmv: uplo must be 'U' or 'L'."
    END SELECT
    vec = res_vec
  END SUBROUTINE mpf_trmv

  !> @brief Solves a triangular system of equations (TRSV). x <- inv(A)*x
  SUBROUTINE mpf_trsv(mat, vec, uplo)
    TYPE(mpf_matrix), INTENT(IN) :: mat
    TYPE(mpf_vector), INTENT(INOUT) :: vec
    CHARACTER(LEN=1), INTENT(IN) :: uplo ! 'U' for Upper, 'L' for Lower
    TYPE(mpf), ALLOCATABLE :: x(:), b(:)
    TYPE(mpf) :: sum_val
    INTEGER :: i, j

    IF (mat%m /= mat%n .OR. mat%n /= vec%n) STOP "mpf_trsv: Matrix must be square and dims must match."
    ALLOCATE(x(vec%n), b(vec%n))

    ! Convert from block-float to standard hpf for calculations
    DO i = 1, vec%n
      b(i) = new_mpf_from_mpi_exp(vec%mantissas(i), vec%exponent)
    END DO

    SELECT CASE (uplo)
    CASE ('L', 'l') ! Forward substitution
      DO i = 1, mat%m
        sum_val = new_mpf_from_integer(0_8)
        DO j = 1, i - 1
          sum_val = sum_val + new_mpf_from_mpi_exp(mat%mantissas(i,j), mat%exponent) * x(j)
        END DO
        x(i) = (b(i) - sum_val) / new_mpf_from_mpi_exp(mat%mantissas(i,i), mat%exponent)
      END DO
    CASE ('U', 'u') ! Backward substitution
      DO i = mat%m, 1, -1
        sum_val = new_mpf_from_integer(0_8)
        DO j = i + 1, mat%n
          sum_val = sum_val + new_mpf_from_mpi_exp(mat%mantissas(i,j), mat%exponent) * x(j)
        END DO
        x(i) = (b(i) - sum_val) / new_mpf_from_mpi_exp(mat%mantissas(i,i), mat%exponent)
      END DO
    CASE DEFAULT
      STOP "mpf_trsv: uplo must be 'U' or 'L'."
    END SELECT

    ! Convert result back to a block-float vector
    vec = new_mpf_vector_from_mpfs(x)
  END SUBROUTINE mpf_trsv

  !> @brief Computes matrix-matrix product where one matrix is triangular (TRMM). B <- A*B
  SUBROUTINE mpf_trmm(mat_a, mat_b, uplo)
    TYPE(mpf_matrix), INTENT(IN) :: mat_a
    TYPE(mpf_matrix), INTENT(INOUT) :: mat_b
    CHARACTER(LEN=1), INTENT(IN) :: uplo ! 'U' for Upper, 'L' for Lower
    TYPE(mpf_matrix) :: res_mat
    TYPE(mpi) :: element_sum
    INTEGER :: i, j, k

    IF (mat_a%n /= mat_b%m) STOP "mpf_trmm: Inner dimensions do not match."
    res_mat%m = mat_a%m; res_mat%n = mat_b%n
    ALLOCATE(res_mat%mantissas(res_mat%m, res_mat%n))
    res_mat%exponent = mat_a%exponent + mat_b%exponent

    SELECT CASE (uplo)
    CASE ('U', 'u')
      DO j = 1, res_mat%n
        DO i = 1, res_mat%m
          element_sum = new_mpi_from_integer(0_8)
          DO k = i, mat_a%n
            element_sum = element_sum + (mat_a%mantissas(i,k) * mat_b%mantissas(k,j))
          END DO
          res_mat%mantissas(i,j) = element_sum
        END DO
      END DO
    CASE ('L', 'l')
      DO j = 1, res_mat%n
        DO i = 1, res_mat%m
          element_sum = new_mpi_from_integer(0_8)
          DO k = 1, i
            element_sum = element_sum + (mat_a%mantissas(i,k) * mat_b%mantissas(k,j))
          END DO
          res_mat%mantissas(i,j) = element_sum
        END DO
      END DO
    CASE DEFAULT
      STOP "mpf_trmm: uplo must be 'U' or 'L'."
    END SELECT
    mat_b = res_mat
  END SUBROUTINE mpf_trmm

  !> @brief Solves a triangular matrix equation (TRSM). B <- inv(A)*B
  SUBROUTINE mpf_trsm(mat_a, mat_b, uplo)
    TYPE(mpf_matrix), INTENT(IN) :: mat_a
    TYPE(mpf_matrix), INTENT(INOUT) :: mat_b
    CHARACTER(LEN=1), INTENT(IN) :: uplo ! 'U' for Upper, 'L' for Lower
    TYPE(mpf_vector) :: b_col
    TYPE(mpf), ALLOCATABLE :: solution_mpfs(:,:)
    INTEGER :: i, j

    IF (mat_a%m /= mat_a%n) STOP "mpf_trsm: Matrix A must be square."
    IF (mat_a%n /= mat_b%m) STOP "mpf_trsm: Inner dimensions do not match."

    ALLOCATE(solution_mpfs(mat_b%m, mat_b%n))

    ! Solve A*x=b for each column of B
    DO j = 1, mat_b%n
      b_col = new_mpf_vector_from_mpfs( &
          (/ (new_mpf_from_mpi_exp(mat_b%mantissas(i,j), mat_b%exponent), i=1,mat_b%m) /) )

      CALL mpf_trsv(mat_a, b_col, uplo)

      DO i = 1, b_col%n
          solution_mpfs(i,j) = new_mpf_from_mpi_exp(b_col%mantissas(i), b_col%exponent)
      END DO
    END DO

    mat_b = new_mpf_matrix_from_mpfs(solution_mpfs)
  END SUBROUTINE mpf_trsm

  SUBROUTINE normalize_mpf_vector(vec)
    TYPE(mpf_vector), INTENT(INOUT) :: vec
    TYPE(mpf) :: temp_mpf
    INTEGER :: i

    IF (vec%n == 0) RETURN

    DO i = 1, vec%n
      CALL normalize_mpi(vec%mantissas(i))
    END DO
  END SUBROUTINE normalize_mpf_vector

  SUBROUTINE normalize_mpf_matrix(mat)
    TYPE(mpf_matrix), INTENT(INOUT) :: mat
    TYPE(mpf) :: temp_mpf
    INTEGER :: i, j

    IF (mat%m == 0 .OR. mat%n == 0) RETURN

    DO j = 1, mat%n
      DO i = 1, mat%m
        CALL normalize_mpi(mat%mantissas(i,j))
      END DO
    END DO
  END SUBROUTINE normalize_mpf_matrix

END MODULE multi_precision_linear_algebra_mod