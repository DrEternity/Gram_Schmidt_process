program main
    implicit none
    integer :: n, i, j, lda, ldwork, info
    double precision, dimension (:, :), allocatable :: A
    double precision, dimension (:), allocatable :: tau, sample, work
    real :: begin, end
    n = 2816
    read *,n
    lda = n
    ldwork = -1

    allocate (A(n, n), tau(n), sample(n))
    call random_number(A)

    call dgeqrf(n, n, A, lda, tau, sample, ldwork, info)
    ldwork = sample(1)
    allocate(work(ldwork))

    call cpu_time(begin)
    call dgeqrf(n, n, A, lda, tau, work, ldwork, info)
    call cpu_time(end)

    write (*,*) "N = ", n, "The time: ", end - begin

end program main
