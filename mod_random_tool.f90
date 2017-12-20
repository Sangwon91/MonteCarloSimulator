module mod_random_tool

    implicit none
    
    private

    integer, parameter :: r8 = selected_real_kind(15,300)

    logical :: MakeSeedQ = .false.

    integer :: iSeed, iSeed0

    public  :: PrintSeed, RandomSettingAll, RanMat, RanSample, Ranf

    interface RanMat

        module procedure RanTen, RanVec, Ranf
    
    end interface RanMat

    contains
    
    !********************

    subroutine RanfSetting

        real(kind=r8) :: x

        if (iSeed >= 0) stop 'not negative!'
    
        x = Ranf() ! <= WHY? 
        iSeed0 = iSeed       ! initial seed값을 저장한다.
        iSeed  = 1

    end subroutine RanfSetting

    !********************

    subroutine GiveRandomSeed(seed_)
    
        implicit none

        integer, intent(out) :: seed_
        integer              :: i
        integer, dimension(8):: dt
        character(len=10) :: time_, date_, zone_

        call date_and_time(date_, time_, zone_, dt)

        seed_ = 1

        do i = 1, 8
            seed_ = ( mod(mod(seed_*dt(i), 1000)*dt(8), 100000) + 1 )
        end do

        seed_ = -abs(seed_)

    end subroutine GiveRandomSeed

    !********************

    subroutine PrintSeed(unit_)

        integer, optional, intent(in) :: unit_ 

        if (.not.present(unit_)) write(*,"(A, I7)")     'initial seed = ', iseed0
        if (present(unit_))       write(unit_,"(A, I7)") 'initial seed = ', iseed0

    end subroutine PrintSeed 

    !********************    

    subroutine RandomSettingAll

        call GiveRandomSeed (iseed)
        call RanfSetting
        call PrintSeed 
        call PrintSeed (1)

    end subroutine RandomSettingAll

    !********************   

    real(kind=r8) function Ranf() result(ran)                ! ran2 in Numerical recipe

        integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        real(KIND=r8) ::  AM,EPS,RNMX
        parameter(IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
          &IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
          &IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        integer idum2,j,k,iv(NTAB),iy
        save iv,iy,idum,idum2  ! changed
        data idum2/123456789/, iv/NTAB*0/, iy/0/

        if ( iseed < 0 ) then  ! changed
          idum = iseed
        endif

        if ( idum.le.0 ) then  ! initialize
          idum = max(-idum,1)  ! make idum positive
          idum2=idum
          do j = NTAB+8, 1, -1
            k = idum / IQ1
            idum = IA1 * ( idum - k * IQ1 ) - k * IR1
            if ( idum.lt.0 ) idum = idum + IM1
            if ( j.le.NTAB ) iv(j) = idum
          end do
          iy = iv(1)
        end if
        k = idum / IQ1
        idum = IA1 * ( idum - k * IQ1 ) - k * IR1
        if ( idum.lt.0 ) idum = idum + IM1
        k = idum2 / IQ2
        idum2 = IA2 * ( idum2 - k * IQ2 ) - k * IR2
        if ( idum2.lt.0 ) idum2 = idum2 + IM2
        j = 1 + iy / NDIV
        iy = iv(j) - idum2
        iv(j) = idum
        if ( iy.lt.1 ) iy = iy + IMM1
        ran = min( AM*iy , RNMX )

    end function Ranf

    !********************

    function RanTen(m, n, matType) ! make m by n matrix

        implicit none

        integer :: i, j
        integer, intent(in) :: m, n
        character(len=2), optional, intent(in) :: matType
        real(kind=r8), dimension(m,n) :: RanTen
       
        ! matType : 1. u = upper, l = lower, d = diagonal matrix
        
        RanTen = 0.0_r8;

        if (present(matType).and.(m==n)) then
        
            select case(matType)
            
                case('u','U') ! upper triangle matrix
                
                    do i = 1,m; do j = 1,n
                        if (j>=i) RanTen(i,j) = Ranf()
                    end do; end do; return
                
                case('u-','U-','-u','-U')
                
                    do i = 1,m; do j = 1,n
                        if (j>=i) RanTen(i,j) = 2.0_r8 * Ranf() - 1.0_r8
                    end do; end do; return
    
                case('l','L') ! lower triangle matrix
                
                    do i = 1,m; do j = 1,n
                        if (i>=j) RanTen(i,j) = Ranf()
                    end do; end do; return
                    
                case('d','D') ! diagonal matrix
                
                    do i = 1,m; do j = 1,n
                        if (j==i) RanTen(i,j) = 2.0_r8 * Ranf() - 1.0_r8
                    end do; end do; return
       
                case default

                    do i = 1,m; do j = 1,n
                        RanTen(i,j) = Ranf()
                    end do; end do; return                
                
            end select

        end if        
        
        do i = 1,m; do j = 1,n
            RanTen(i,j) = Ranf()
        end do; end do;
       
    end function RanTen

    function RanVec(m) ! make m by n matrix

        implicit none

        integer :: i
        integer, intent(in) :: m
        real(kind=r8), dimension(m) :: RanVec

        do i = 1,m; RanVec(i) = Ranf(); end do           

    end function RanVec

    !********************

    subroutine RanSample(data_, nSmple, output, nChse) ! data array, number of sample & Choose number
 
        implicit none

        integer :: index_, maxS, maxC, temp  ! maximum Sample, maximum Choose
        integer, intent(in) :: nSmple, nChse
        real(kind=r8), dimension(nSmple), intent(in) :: data_
        real(kind=r8), dimension(nChse), intent(out) :: output
        real(kind=r8), dimension(nSmple) :: dData    ! dummy data for change
 
        maxS = nSmple; maxC = nChse; dData = data_
        
        if (maxC > maxS) stop 'choose > sample'
        do index_ = 1, maxC

            temp = int(Ranf()*real(maxS))+1
            output(index_) = dData(temp)
            dData(temp:maxS-1) = dData(temp+1:maxS)
            maxS = maxS-1             

        end do               

    end subroutine RanSample

    !********************  

end module mod_random_tool

!============================================================