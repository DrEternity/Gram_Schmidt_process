program main
    implicit none
    integer , parameter :: n = 1024
    integer :: i, j
    double precision, dimension(:), allocatable :: A, Q, R, C

    allocate (A(n * n), Q(n * n), R(n * n), C(n * n))
    
    !A = [1, 5, -4, -8, -10, 4, 5, 6, 1, 3, 5, 8, 7, 9, -4, 2]
    !B = [4, 5, -4, -8, -10, 4, 6, 6, 1, 3, 2, 8, 7, 10, -4, 2]

    call random_number(A)
    A = A * 1000
    !write(*, *) A
    call QR(A, Q, R, n)
    call matmul_checker(Q, R, C, n)
    !write(*, *) Q
    write(*,*) Frob(A, C, n)
    
contains


subroutine matmul_checker(A, B, C, n)
    implicit none
    integer, intent(in) :: n

    double precision, allocatable, intent(inout) :: C(:)
    double precision, allocatable, intent(in) :: A(:)
    double precision, allocatable, intent(in) :: B(:)
    integer :: k, j, i

    do j = 1, n
        do k = 1, n
            do i = 1, n
                C(i + (j - 1) * n) = C(i + (j - 1) * n) + A(i + (k - 1) * n) * B(k + (j - 1) * n)
            end do
        end do
    end do

end subroutine matmul_checker


function Frob(A, B, n) result(error)
    implicit none
    real :: error
    integer, intent(in) :: n
    double precision, allocatable, dimension(:), intent(in) :: A, B
    integer :: i, j
    error = 0
    do i = 1, n
        do j = 1, n
            error = error + (A(i + (j - 1) * n) - B(i + (j - 1) * n)) * (A(i + (j - 1) * n) - B(i + (j - 1) * n))
        end do
    end do
    error = sqrt(error)

end function Frob


function norm(a, n) result (ans)
    implicit none
    real :: ans
    double precision, allocatable, dimension(:), intent(in) :: a
    integer, intent(in) :: n
    integer :: i
    ans = 0
    do i = 1, n
        ans = ans + a(i) * a(i) 
    end do
    ans = sqrt(ans)
    
end function norm


subroutine mult_A_b(A, B, res, n_A, m_A, n_B, transpose, offset_B, offset_res)
    implicit none
    integer, intent(in) :: n_A, m_A, n_B, offset_B, offset_res 
    integer :: i, j

    double precision, allocatable, intent(inout) :: res(:)
    double precision, allocatable, intent(in) :: A(:)
    double precision, allocatable, intent(in) :: B(:)
    logical :: transpose
    real :: tmp, sum

    double precision, dimension(4096) :: tmp_ram
    tmp_ram = tmp_ram * 0.0
    

    if (transpose) then
        do i = 1, n_A
            tmp_ram(i) = B(i + offset_B)
        end do 
        
        do i = 1, m_A
            sum = 0
            do j = 1, n_A
                sum = sum + A((i - 1) * n_A + j) * tmp_ram(j)
            end do
            res(i + offset_res) = sum
        end do
    else
        do i = 1, m_A
            tmp = B(i + offset_B)
            do j = 1, n_A
                tmp_ram(j) = tmp_ram(j) + A(j + n_A * (i - 1)) * tmp
            end do
        end do
        
        do i = 1, n_A
            res(i + offset_res) = tmp_ram(i)
        end do
    end if
end subroutine mult_A_b

    
subroutine iteration(A, Q, R, n, num_iter)
    implicit none
    integer, intent(in) :: n, num_iter
    integer :: i
    real :: norm_q

    double precision, allocatable, intent(in) :: A(:)
    double precision, allocatable, intent(inout) :: Q(:)
    double precision, allocatable, intent(inout) :: R(:)

    double precision, dimension(:), allocatable :: tmp
    allocate (tmp(n))


    call mult_A_b(Q, A, R, n, num_iter - 1, n, .true., (num_iter - 1) * n, (num_iter - 1) * n)
    call mult_A_b(Q, R, tmp, n, num_iter - 1, num_iter - 1, .false., (num_iter - 1) * n, 0)


    do i = 1, n
        tmp(i) = -1 * tmp(i)
        tmp(i) = tmp(i) + A(i + (num_iter - 1) * n)
    end do

    norm_q = norm(tmp, n)
    do i = 1, n
        Q((num_iter - 1) * n + i) = tmp(i) / norm_q
    end do
    R((num_iter - 1) * n + (num_iter)) = norm_q

end subroutine iteration


subroutine QR(A, Q, R, n)
    double precision, allocatable, intent(in) :: A(:)
    double precision, allocatable, intent(inout) :: Q(:)
    double precision, allocatable, intent(inout) :: R(:)
    integer, intent(in) :: n
    integer :: i
    
    real :: norm_q_1

    norm_q_1 = norm(A, n)
    R(1) = norm_q_1
    do i = 1, n
        Q(i) = A(i) / norm_q_1
    end do

    do i = 2, n
        call iteration(A, Q, R, n, i)
    end do
end subroutine QR
end program main
    
