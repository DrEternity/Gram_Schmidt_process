program main
    implicit none
    integer :: n = 512
    integer :: s = 128 
    integer :: i, j
    real :: begin, end
    double precision, dimension(:), allocatable :: A, C, Q, R

    read *,n
    allocate (A(0 : n * n - 1), R(0 : n * n - 1), Q(0 :  n * n - 1))
    call random_number(A)
    A = A * 100
    Q = A
    
    call cpu_time(begin)
    call QR(Q, R, n, s)
    call cpu_time(end)
    write (*,*) "N = ", n, "The time: ", end - begin
        

    !allocate (C(0: n * n - 1))
    !call matmul_checker(Q, R, C, n)
    !write(*,*) "Frob:", Frob(A, C, n)
    
contains


subroutine matmul_checker(A, B, C, n)
    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: C(0:)
    double precision, intent(in) :: A(0:), B(0:)

    call dgemm('N', 'N', n, n, n, 1.0d0, A, n, B, n, 0.0d0, C, n);
end subroutine matmul_checker


function Frob(A, B, n) result(error)
    implicit none
    real :: error
    integer, intent(in) :: n
    double precision, dimension(0:), intent(in) :: A, B
    integer :: i, j
    error = 0
    do i = 0, n - 1
        do j = 0, n - 1
            error = error + (A(i + j * n) - B(i + j * n)) * (A(i + j * n) - B(i + j * n))
        end do
    end do
    error = sqrt(error)
end function Frob


function norm(a, n) result (ans)
    implicit none
    real :: ans
    double precision, dimension(0:), intent(in) :: a
    integer, intent(in) :: n
    integer :: i
    ans = 0
    do i = 0, n - 1
        ans = ans + a(i) * a(i)
    end do
    ans = sqrt(ans)
end function norm


recursive subroutine intraortho(Q, R, n, s) 
    integer, intent(in) :: n, s
    double precision, intent(inout), dimension(0:) :: Q, R
    
    integer :: i
    double precision :: norm_q_1
    
    norm_q_1 = norm(Q, n)
    R(0) = norm_q_1
    Q(:n - 1) = Q(:n - 1) / norm_q_1
    
    do i = 1, s - 1
        call DGEMV('T', n, i, 1.0d0, Q, n, Q(n * i :), 1, 1.0d0, R(n * i :), 1)
        call DGEMV('N', n, i, -1.0d0, Q, n, R(n * i :), 1, 1.0d0, Q(n * i :), 1)
        call intraortho(Q(n * i:), R(n * i + i:), n, 1)
    end do
end subroutine intraortho


subroutine QR(Q, R, n, s)
    integer, intent(in) :: n, s
    double precision, intent(inout), dimension(0:) :: Q, R
    integer :: k
    
    call intraortho(Q, R, n, s)
    do k=1, (n/s - 1)
        call dgemm('T', 'N', k * s, s, n, 1.0d0, Q, n, Q(n * s * k:), n, 1.0d0, R(n * s * k:), n)
        call dgemm('N', 'N', n, s, k * s, -1.0d0, Q, n, R(n * s * k:), n, 1.0d0, Q(n * s * k:), n)        
        call intraortho(Q(n * s * k:), R(n * s * k + s * k:), n, s)
    end do

end subroutine QR

end program main
