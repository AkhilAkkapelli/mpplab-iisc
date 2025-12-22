program main

    USE multi_precision_integer_mod
    implicit none

    type(mpi)   :: mpi_val
    integer(8)  :: scalar
    integer(8)  :: ntests, npass= 0
    integer     :: unit_num= 10, iostatus
    integer(8)  :: k
    integer(8)  :: c1(COEFFS_LIMIT), c2(COEFFS_LIMIT)

#if defined(__debug__)
    integer(8)  :: i
#endif

    open(unit= unit_num, file='./bin/mpi_scalar_mult.txt', status='old', action='read', iostat=iostatus)
    IF (iostatus /= 0) THEN
        PRINT *, "Error: Could not open './bin/mpi_scalar_mult.txt'. Run python generator first."
        STOP
    END IF

    read(unit_num, *)ntests
    print '(A,I0,A)', "Loading ", ntests, " Test Cases from file!"

#if defined(__debug__)
    ntests= 1
#endif

    do k= 1, ntests
        read(unit_num, *)c1, scalar, c2
        mpi_val%coeffs= c1
    
        ! debug statements
#if defined(__debug__)
        do i= 1, COEFFS_LIMIT
            write(*, '(I0, A)', advance="no") mpi_val%coeffs(i), " "
        end do
        print *, ""

        print '(I0)', scalar
#endif

        call mpi_multiply_by_scalar(mpi_val, scalar)

        ! debug statements
#if defined(__debug__)
        do i= 1, COEFFS_LIMIT
            write(*, '(I0, A)', advance="no") mpi_val%coeffs(i), " "
        end do
        print *,""

        do i= 1, COEFFS_LIMIT
            write(*, '(I0, A)', advance="no") c2(i), " "
        end do
        print *,""
#endif

        if(any(mpi_val%coeffs /= c2))then
            continue
        else
            npass= npass + 1
        end if

    end do
    close(unit_num)

    if (npass == ntests) then
        print '(A,I0,A,I0,A)', "[PASS] Scalar Multiply: Randomly Generated cases: ", npass, "/", ntests, " cases passed" 
    else
        print '(A,I0,A,I0,A)', "[FAIL] Scalar Multiply: Randomly Generated cases: ", npass, "/", ntests, " cases passed" 
    end if

end program main