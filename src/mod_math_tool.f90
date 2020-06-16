module mod_math_tool

    implicit none

    private

    integer, parameter :: r8 = selected_real_kind(15,300)

    public :: Det, Inverse, Factorial, IdMat, SimpleSort

    interface Factorial

        module procedure FactorialInt, FactorialReal

    end interface Factorial

    contains

    !********************

    function Det (mat, dim)

        implicit none

        integer, intent(in) :: dim
        integer :: i, j, pivot
        real(kind=r8), dimension(dim,dim), intent(in) ::  mat
        real(kind=r8), dimension(dim,dim)             :: dmat
        real(kind=r8), dimension(dim)   :: temp
        real(kind=r8) :: assumpZero, Det, maxComp, sign_

        dmat = mat

        assumpZero = maxval(dmat) * 10e-8

        sign_ = 1.0_r8

        do i = 1, dim

            ! pivoting

            pivot   = maxloc(abs(dmat(i:dim,i)),1) + i - 1 ! pivot과 maximum value를 찾는데 내장 함수 사용
            maxComp = dmat(pivot,i)

            if (pivot/=i) then

                temp          = dmat(i,:)           ! 행 교환, singular를 찾아내기 위하여
                dmat(i,:)     = dmat(pivot,:)
                dmat(pivot,:) = temp
                sign_         = -sign_

            end if

            if (abs(maxComp) < assumpZero) then; Det=0.0_r8; return; end if

            do j = i+1, dim    ! 가우스 소거법 적용
                dmat(j,:) = dmat(j,:) - dmat(i,:) * dmat(j,i)/maxComp
            end do

        end do

        det = sign_

        do i = 1, dim
            det = det * dmat(i,i)
        end do

    end function Det

    !********************

    subroutine Inverse (mat, iMat, singularQ, dim)

        implicit none

        integer, intent(in):: DIM ! dimension of matrix
        integer            :: pivot, i, j
        logical :: singularQ
        real(kind=r8), dimension(DIM,DIM), intent(in)  :: mat
        real(kind=r8), dimension(DIM,DIM), intent(out) :: iMat
        real(kind=r8), dimension(DIM,DIM) :: dMat ! dummy matrix
        real(kind=r8), dimension(DIM)   :: temp, itemp
        real(kind=r8) :: assumpZero, maxComp, coeff ! maxComp : max component, coeff : coefficient of elimination

        ! 초기화
        singularQ = .false.
        dMat = mat
        iMat = idMat(dim) ! identity matrix
        assumpZero = maxval(dMat) * 10e-8

        do i = 1, DIM

            ! pivoting

            pivot   = maxloc(abs(dMat(i:DIM,i)),1) + i - 1 ! pivot과 maximum value를 찾는데 내장 함수 사용
            maxComp = dMat(pivot,i)

            if (pivot /= i) then

                 temp = dMat(i,:)           ! 행 교환
                itemp = iMat(i,:)
                dMat(i,:) = dMat(pivot,:)
                iMat(i,:) = iMat(pivot,:)
                dMat(pivot,:) =  temp
                iMat(pivot,:) = itemp

            end if

            if (abs(maxComp) < assumpZero) then;
                singularQ = .true.; iMat = 0.0_r8; return
            end if

            do j = i+1, DIM    ! 가우스 소거법 적용

                coeff     = dMat(j,i)
                dMat(j,:) = dMat(j,:) - dMat(i,:) * coeff/maxComp
                iMat(j,:) = iMat(j,:) - iMat(i,:) * coeff/maxComp

            end do

        end do

        do i = DIM, 1, -1      ! 후진 대입법으로 계산

           coeff     = dMat(i,i)

           dMat(i,:) = dMat(i,:) / coeff
           iMat(i,:) = iMat(i,:) / coeff

           do j = i-1, 1, -1

               coeff     = dMat(j,i)

               dMat(j,:) = dMat(j,:) - coeff*dMat(i,:)
               iMat(j,:) = iMat(j,:) - coeff*iMat(i,:)

            end do

        end do

    end subroutine Inverse

    !********************

    function IdMat (n) ! make n x n identity matrix

        integer :: i, j
        integer, intent(in) :: n
        real(kind=r8), dimension(n,n) :: IdMat

        IdMat = 0.0_r8

        forall(i=1:n ,j=1:n, i==j)
            idMat(i,j) = 1.0_r8
        end forall

    end function IdMat

    !********************

    subroutine SimpleSort (array, n, mode)

        implicit none

        integer :: i, d1, d2
        integer, intent(in) :: n
        character(len=1), optional, intent(in) :: mode
        logical :: compliteQ
        real(kind=r8), dimension(n), intent(out) :: array
        real(kind=r8) :: temp

        d1 = 0; d2 = 1;
        if (present(mode)) then
            select case (mode)
                case('u','U')
                    d1 = 1; d2 = 0;
                case('d','D')
                    d1 = 0; d2 = 1;
                case default
                    stop 'unclassified sort mode'
            end select
        end if

        compliteQ = .false.

        do while(.not.compliteQ)

            compliteQ = .true.

            do i = 1, n-1

                if (array(i+d1) > array(i+d2)) then

                    temp        = array(i+d1)
                    array(i+d1) = array(i+d2)
                    array(i+d2) = temp
                    compliteQ   = .false.

                end if
                !write(*,"(A, 6F6.2, L2)") 'sorting... ', array, compliteQ
            end do
            !write(*,"(A)") 'loop end'
        end do

    end subroutine SimpleSort

    !********************

    function FactorialInt(x)

        implicit none

        integer :: dNum    ! dummy number
        integer, intent(in) :: x
        integer :: FactorialInt

        dNum = x; FactorialInt = 1

        do while (dNum>1); FactorialInt=FactorialInt*dNum; dNum=dNum-1; end do

    end function FactorialInt

    function FactorialReal(x)

        implicit none

        real(kind=r8) :: dNum    ! dummy number
        real(kind=r8), intent(in) :: x
        real(kind=r8) :: FactorialReal

        dNum = x; FactorialReal = 1.0_r8

        do while (dNum>1.0_r8); FactorialReal=FactorialReal*dNum; dNum=dNum-1.0_r8; end do

    end function FactorialReal

    !********************

end module mod_math_tool

!============================================================
