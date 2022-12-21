program main
    implicit none
    integer :: n = 1024
    integer :: size_block = 32
    integer :: s = 128
    
    integer :: i, j
    double precision, dimension(:), allocatable :: A, C, Q, R

    allocate (A(0 : n * n - 1), R(0 : n * n - 1), Q(0 :  n * n - 1))
    
    call random_number(A)
    A = A * 100
    Q = A

    call QR(Q, R, n, s, size_block)

    !allocate (C(0: n * n - 1))
    !call matmul_checker(Q, R, C, n)
    !write(*,*) "Frob:", Frob(A, C, n)
    
contains


subroutine matmul_checker(A, B, C, n)
    implicit none
    integer, intent(in) :: n

    double precision, intent(inout) :: C(0:)
    double precision, intent(in) :: A(0:)
    double precision, intent(in) :: B(0:)
    integer :: k, j, i

    do j = 0, n - 1
        do k = 0, n - 1
            do i = 0, n - 1
                C(i + j * n) = C(i + j * n) + A(i + k * n) * B(k + j * n)
            end do
        end do
    end do
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


subroutine mult_block_A_B(A, B, res, n_A, m_A, n_B, m_B, i_B, j_B, row_B)
    implicit none
    integer, intent(in) :: n_A, m_A, n_B, m_B, i_B, j_B, row_B
    double precision, intent(in), dimension(0: ) :: A
    double precision, intent(in), dimension(0: ) :: B
    double precision, intent(inout), dimension(0: ) :: res
    integer :: i, k, j
    
    do k = 0, m_B - 1
        do i = 0, m_A - 1
            do j = 0, n_A - 1
                res(k * n_A + j) = res(k * n_A + j) + A(i * n_A + j) * B((k + j_B) * row_B + (i + i_B));
            end do
        end do
    end do
end subroutine mult_block_A_B


subroutine copy_to_ram(A, destination, i, j, n_A, m_A, row_block, col_block, transpose)
    integer, intent(in) :: i, j, n_A, m_A, row_block, col_block
    logical :: transpose
    double precision, intent(in), dimension(0: ) :: A
    double precision, intent(inout), dimension(0: ) :: destination
    integer :: v, w

    if (transpose) then
        do v = 0, row_block - 1
            do w = 0, col_block - 1
                destination(w * row_block + v) = A((i + v) * n_A + (j + w))
            end do
        end do
    else
        do w = 0, col_block - 1
            do v = 0, row_block - 1
                destination(w * row_block + v) = A((j + w) * n_A + (i + v))
            end do
        end do
    end if
end subroutine copy_to_ram


subroutine return_from_ram(C, ram, row_C, row_block, col_block, i, j)
    integer, intent(in) :: row_C, row_block, col_block, i, j
    double precision, intent(inout), dimension(0: ) :: C
    double precision, intent(inout), dimension(0: ) :: ram
    integer :: w1, w2

    do w1 = 0, col_block - 1
        do w2 = 0, row_block - 1
            C((j + w1) * row_C + (i + w2)) = C((j + w1) * row_C + (i + w2)) +  ram(w1 * row_block + w2)
            ram(w1 * row_block + w2) = 0
        end do
    end do
end subroutine return_from_ram


subroutine swap(a, b)
  integer, intent(inout) :: a, b
  integer :: temp
  temp = a
  a = b
  b = temp
end subroutine swap


subroutine mult_A_B(A, row_A, col_A, n_A_, m_A_, B, row_B, col_B, n_B, m_B, C, row_C, col_C, n_C, m_C, transpose, size_block)
    implicit none
    integer, intent(in) :: row_A, col_A, row_B, col_B, n_B, m_B, row_C, col_C, n_C, m_C 
    integer, intent(in) :: size_block
    integer, intent(in) :: n_A_, m_A_
    logical :: transpose
    
    double precision, intent(in), dimension(0: ) :: A
    double precision, intent(in), dimension(0: ) :: B
    double precision, intent(inout), dimension(0: ) :: C
    
    integer i, j, k, n_A, m_A
    
    double precision, dimension(0:32 * 32 - 1) :: block_A, block_res ! static
    block_A = 0
    block_res = 0
    n_A = n_A_
    m_A = m_A_

    if (transpose) then
        call swap(n_A, m_A)
    end if


    do i = 0, n_A - MOD(n_A, size_block) - 1, size_block
        do k = 0, m_A - MOD(m_A, size_block) - 1, size_block
            call copy_to_ram(A, block_A, i, k, row_A, col_A, size_block, size_block, transpose)
            do j = 0, m_B - MOD(m_B, size_block) - 1, size_block
                call mult_block_A_B(block_A, B, block_res, size_block, size_block, size_block, size_block, k, j, row_B)
                call return_from_ram(C, block_res, row_C, size_block, size_block, i, j)
            end do
            if (MOD(m_B, size_block) > 0) then
                call mult_block_A_B(block_A, B, block_res, size_block, size_block, size_block, &
                MOD(m_B, size_block), k, m_B - MOD(m_B, size_block), row_B)
            end if
            call return_from_ram(C, block_res, row_C, size_block, MOD(m_B, size_block), i, m_B - MOD(m_B, size_block))
        end do
        if (MOD(m_A, size_block) > 0) then
            call copy_to_ram(A, block_A, i, m_A - MOD(m_A, size_block), row_A, col_A, size_block, MOD(m_A, size_block), transpose)
            do j = 0, m_B - MOD(m_B, size_block) - 1, size_block
                call mult_block_A_B(block_A, B, block_res, size_block, MOD(m_A, size_block), &
                MOD(m_A, size_block), size_block, m_A - MOD(m_A, size_block), j, row_B)
                call return_from_ram(C, block_res, row_C, size_block, size_block, i, j)
            end do
            if (MOD(m_B, size_block) > 0) then  
                call mult_block_A_B(block_A, B, block_res, size_block, MOD(m_A, size_block), &
                MOD(m_A, size_block), MOD(m_B, size_block), m_A - MOD(m_A, size_block), m_B - MOD(m_B, size_block), row_B)
            end if
            call return_from_ram(C, block_res, row_C, size_block, MOD(m_B, size_block), i, m_B - MOD(m_B, size_block))
        end if
    end do

    if (MOD(n_A, size_block) > 0) then
        do k = 0, m_A - MOD(m_A, size_block) - 1, size_block
            call copy_to_ram(A, block_A, n_A - MOD(n_A, size_block), k, row_A, col_A, MOD(n_A, size_block), size_block, transpose)
            do j = 0, m_B - MOD(m_B, size_block) - 1, size_block
                call mult_block_A_B(block_A, B, block_res, MOD(n_A, size_block), size_block, size_block, size_block, k, j, row_B)
                call return_from_ram(C, block_res, row_C, MOD(n_A, size_block), size_block, n_A - MOD(n_A, size_block), j)
            end do
            if (MOD(m_B, size_block) > 0) then
                call mult_block_A_B(block_A, B, block_res, MOD(n_A, size_block), size_block, size_block, &
                MOD(m_B, size_block), k, m_B - MOD(m_B, size_block), row_B)
            end if
            call return_from_ram(C, block_res, row_C, MOD(n_A, size_block), MOD(m_B, size_block), &
            n_A - MOD(n_A, size_block), m_B - MOD(m_B, size_block))
        end do
        if (MOD(m_A, size_block) > 0) then
            call copy_to_ram(A, block_A, n_A - MOD(n_A, size_block), m_A - MOD(m_A, size_block), &
            row_A, col_A, MOD(n_A, size_block), MOD(m_A, size_block), transpose)
            do j = 0, m_B - MOD(m_B, size_block) - 1, size_block
                call mult_block_A_B(block_A, B, block_res, MOD(n_A, size_block), MOD(m_A, size_block), &
                MOD(m_A, size_block), size_block, m_A - MOD(m_A, size_block), j, row_B)
                call return_from_ram(C, block_res, row_C, MOD(n_A, size_block), size_block, n_A - MOD(n_A, size_block), j)
            end do
            if (MOD(m_B, size_block) > 0) then 
                call mult_block_A_B(block_A, B, block_res, MOD(n_A, size_block), MOD(m_A, size_block), &
                MOD(m_A, size_block), MOD(m_B, size_block), m_A - MOD(m_A, size_block), m_B - MOD(m_B, size_block), row_B) 
            end if
            call return_from_ram(C, block_res, row_C, MOD(n_A, size_block), MOD(m_B, size_block), n_A - MOD(n_A, size_block), &
            m_B - MOD(m_B, size_block))
        end if
    end if
end subroutine mult_A_B


recursive subroutine intraortho(Q, R, n, s, size_block) 
    integer, intent(in) :: n, s, size_block
    double precision, intent(inout), dimension(0:) :: Q, R
    
    integer :: i
    double precision :: norm_q_1
    
    norm_q_1 = norm(Q, n)
    R(0) = norm_q_1
    Q(:n - 1) = Q(:n - 1) / norm_q_1
    
    do i = 1, s - 1
        call mult_A_B(Q, n, n, n, i, &
                Q(n * i :), n, n, n, 1, &
                R(n * i :), n, n, i, 1, &
                .true., size_block)
        Q(n * i: n * i + n - 1) = Q(n * i: n * i + n - 1) * (-1)
        call mult_A_B(Q, n, n, n, i, &
                R(n * i:), n, n, i, 1, &
                Q(n * i:), n, n, n, 1, &
                .false., size_block)
        Q(n * i: n * i + n - 1) = Q(n * i: n * i + n - 1) * (-1)
        call intraortho(Q(n * i:), R(n * i + i:), n, 1, size_block)
    end do
end subroutine intraortho


subroutine QR(Q, R, n, s, size_block)
    integer, intent(in) :: n, s, size_block
    double precision, intent(inout), dimension(0:) :: Q, R
    integer :: k
    
    call intraortho(Q, R, n, s, size_block)
    do k=1, (n/s - 1)
        call mult_A_B(Q, n, n, n, k * s, &
                Q(n * s * k:), n, n, n, s, &
                R(n * s * k:), n, n, k * s, s, &
                .true., size_block)
        Q(n * s * k: n * s * k + n * s - 1) = Q(n * s * k: n * s * k + n * s - 1) * (-1)
        call mult_A_B(Q, n, n, n, k * s, &
                R(n * s * k:), n, n, k * s, s, &
                Q(n * s * k:), n, n, n, s, &
                .false., size_block)
        Q(n * s * k: n * s * k + n * s - 1) = Q(n * s * k: n * s * k + n * s - 1) * (-1)
        call intraortho(Q(n * s * k:), R(n * s * k + s * k:), n, s, size_block)
    end do

end subroutine QR

end program main
    
