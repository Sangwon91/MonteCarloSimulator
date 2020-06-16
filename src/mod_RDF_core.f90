module mod_RDF_core

    use mod_univ_const, only : CHARGE_CONST, PI   ! to use univertial constant
    use mod_math_tool,  only : Det           ! to use erfc function
    use mod_system_information,   only : molN => molNum, rCut => rCut
    use mod_molecule_information, only : nA => nAtoms

    use mod_config_information,   only : rX => rX, rY => rY, rZ => rZ, &
                                         hBox => hBox, ihBox => ihBox, vol => volume   ! to use configuration information
    implicit none

    private

    integer, parameter :: r8 = selected_real_kind(15,300)

    integer :: maxBin, maxC, mxCmb, nSample ! max category, max category combination
    integer, dimension(:), allocatable :: atomLabel, repetition, nAnB
    integer, dimension(:,:), allocatable :: cate
    real(kind=r8), dimension(:,:), allocatable :: hist, gr
    real(kind=r8) :: delR, maxR
    character(len=2), dimension(:), allocatable :: atomName
    character(len=4), dimension(:), allocatable :: molAmolB
    character(len=4), dimension(:,:), allocatable :: cateName

    public :: GetRDFInf, PrintRDFInf, RDFDeallocate, SampleRDF, MakeRDF, PrintRDF, ResetRDF

    contains

    !########################################

    subroutine GetRDFInf (filename)

        implicit none

        integer                       :: i, j, dummyN, dummyT, nCate, dummyA, dummyB ! A: alpha, B: beta
        integer, parameter            :: UNIT_RDF = 10
        character(len=31), intent(in) :: filename
        character(len=2)              :: dummyC ! C: Character, T: Type

        if (allocated(atomName))  deallocate(atomName)
        if (allocated(atomLabel)) deallocate(atomLabel)
        allocate(atomName(nA),atomLabel(nA))

        open(UNIT_RDF, file=filename)

            read(UNIT_RDF,*) delR   ! it is not tested

            do i = 1, nA

                read(UNIT_RDF, *) dummyN, dummyT, dummyC
                !if (dummyC /= '0')
                atomName(i)  = dummyC
                atomLabel(i) = dummyT

            end do

        close(UNIT_RDF)

        ! initialize


        maxR   = rCut
        maxBin = int(rCut / delR) + 1

        maxC = maxval(atomLabel)

        if (allocated(hist)) deallocate(hist)
        if (allocated(gr))   deallocate(gr)
        if (allocated(cate)) deallocate(cate)
        if (allocated(cateName)) deallocate(cateName)
        if (allocated(repetition)) deallocate(repetition)
        if (allocated(nAnB)) deallocate(nAnB)
        if (allocated(molAmolB)) deallocate(molAmolB)

        mxCmb = (maxC+1)*maxC/2 ! max combination

        allocate(hist(0:mxCmb,maxBin),gr(0:mxCmb,maxBin))
        allocate(cate(maxC,maxC), cateName(maxC,maxC))
        allocate(repetition(maxC))
        allocate(nAnB(mxCmb))
        allocate(molAmolB(mxCmb))

        nSample    = 0
        repetition = 0

        do i = 1, maxC
            do j = 1, nA

                if (atomLabel(j) == i) repetition(i) = repetition(i) + 1

            end do
        end do

        hist  = 0.0
        nCate = 1

        do i = 1, maxC
            do j = i, maxC

                nAnB(nCate) = repetition(i) * repetition(j)
                if (i /= j) nAnB(nCate) = nAnB(nCate) * 2
                cate(i,j) = nCate
                cate(j,i) = nCate
                nCate = nCate + 1

            end do
        end do

        do i = 1, nA
            do j = 1, nA

                dummyA = atomLabel(i)
                dummyB = atomLabel(j)

                cateName(dummyA,dummyB) = trim(atomName(i))//trim(atomName(j))

            end do
        end do

        nCate = 1

        do i = 1, maxC
            do j = i, maxC

                molAmolB(nCate) = trim(cateName(i,j))
                nCate = nCate + 1

            end do
        end do

    end subroutine GetRDFInf

    !########################################

    subroutine PrintRDFInf (unit_)

        implicit none

        integer :: i

        integer, intent(in), optional :: unit_

        if (present(unit_)) then

            write(unit_,"(A/)") 'RDF information ->'

            write(unit_,"('¥ÄR = ', F10.5)") delR

            write(unit_,"(3A5)") 'index', 'type', 'name'
            do i = 1, nA
                write(unit_,"(2I5, A5)") i, atomLabel(i), atomName(i)
            end do

            write(unit_,"('max combination: ', I5)") mxCmb

            write(unit_,"('cate matrix: ')")
            do i = 1, maxC
                write(unit_,"(10I3)") cate(i,:)
            end do
            write(unit_,"('cate name matrix: ')")
            do i = 1, maxC
                write(unit_,"(10A6)") cateName(i,:)
            end do

            write(unit_,"()")

        else

            write(*,"(A/)") 'RDF information ->'

            write(*,"('¥ÄR = ', F10.5)") delR

            write(*,"(3A5)") 'index', 'type', 'name'
            do i = 1, nA
                write(*,"(2I5, A5)") i, atomLabel(i), atomName(i)
            end do

            write(*,"('max combination: ', I5)") mxCmb

            write(*,"('cate matrix: ')")
            do i = 1, maxC
                write(*,"(10I3)") cate(i,:)
            end do
            write(*,"('cate name matrix: ')")
            do i = 1, maxC
                write(*,"(10A6)") cateName(i,:)
            end do

            write(*,"()")

        end if

    end subroutine PrintRDFInf

    !########################################

    subroutine SampleRDF

        implicit none

        integer :: i, j, ii, jj, bin, ap, bt ! alpha, beta

        real(kind=r8) :: rcutsq
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33
        real(kind=r8) :: rxij, ryij, rzij, rij, rijsq, dxij, dyij, dzij, srxij, sryij, srzij, &
                         rrxij, rryij, rrzij

        rcutsq = rCut**2

        b11 = hbox(1,1)
        b12 = hbox(1,2)
        b13 = hbox(1,3)
        b22 = hbox(2,2)
        b23 = hbox(2,3)
        b33 = hbox(3,3)

        ib11 = ihbox(1,1)
        ib12 = ihbox(1,2)
        ib13 = ihbox(1,3)
        ib22 = ihbox(2,2)
        ib23 = ihbox(2,3)
        ib33 = ihbox(3,3)

        nSample = nSample + 1

        do i = 1, molN - 1

            do j = i+1, molN

                rrxij = rX(i,0) - rX(j,0)
                rryij = rY(i,0) - rY(j,0)
                rrzij = rZ(i,0) - rZ(j,0)

                dxij = ib11*rrxij + ib12*rryij + ib13*rrzij
                dyij =            + ib22*rryij + ib23*rrzij
                dzij =                         + ib33*rrzij

                if (dxij >  0.5_r8) dxij = dxij - 1.0_r8
                if (dxij < -0.5_r8) dxij = dxij + 1.0_r8

                if (dyij >  0.5_r8) dyij = dyij - 1.0_r8
                if (dyij < -0.5_r8) dyij = dyij + 1.0_r8

                if (dzij >  0.5_r8) dzij = dzij - 1.0_r8
                if (dzij < -0.5_r8) dzij = dzij + 1.0_r8

                srxij = b11*dxij + b12*dyij + b13*dzij
                sryij =            b22*dyij + b23*dzij
                srzij =                       b33*dzij

                rijsq = srxij*srxij + sryij*sryij + srzij*srzij

                if (rijsq < rcutsq) then

                    rij   = rijsq**0.5_r8
                    bin   = int(rij/delR) + 1
                    hist(0,bin) = hist(0,bin) + 2.  ! center-center distribution

                end if

                do ii = 1, nA; do jj = 1, nA

                rxij = rX(i,ii) - rX(j,jj)
                ryij = rY(i,ii) - rY(j,jj)
                rzij = rZ(i,ii) - rZ(j,jj)

                dxij = ib11*rxij + ib12*ryij + ib13*rzij
                dyij =           + ib22*ryij + ib23*rzij
                dzij =                       + ib33*rzij

                if (dxij >  0.5_r8) dxij = dxij - 1.0_r8
                if (dxij < -0.5_r8) dxij = dxij + 1.0_r8

                if (dyij >  0.5_r8) dyij = dyij - 1.0_r8
                if (dyij < -0.5_r8) dyij = dyij + 1.0_r8

                if (dzij >  0.5_r8) dzij = dzij - 1.0_r8
                if (dzij < -0.5_r8) dzij = dzij + 1.0_r8

                srxij = b11*dxij + b12*dyij + b13*dzij
                sryij =            b22*dyij + b23*dzij
                srzij =                       b33*dzij

                rijsq = srxij*srxij + sryij*sryij + srzij*srzij

                if (rijsq < rcutsq) then

                    rij = rijsq**0.5_r8
                    bin = int(rij/delR) + 1

                    ap = atomLabel(ii)
                    bt = atomLabel(jj)

                    hist(cate(ap,bt),bin) = hist(cate(ap,bt),bin) + 2.

                end if

                end do; end do

            end do

        end do

    end subroutine SampleRDF

    !########################################

    subroutine MakeRDF

        implicit none

        integer :: i, bin

        real(kind=r8), parameter :: VCONST = 4.0_r8 * PI / 3.0_r8
        real(kind=r8) :: upperR, lowerR, delV, rstep

        rstep = dble(nSample)

        do i = 0, mxCmb

            do bin = 1, maxBin

                lowerR = dble(bin - 1) * delR
                upperR = dble(bin) * delR
                delV = VCONST * (upperR**3 - lowerR**3)

                gr(i,bin) = hist(i,bin) / rstep / dble(molN) / dble(molN) / delV * vol

            end do

            if (i /= 0) then
                !write(*,"(I5)") i
                gr(i,:) = gr(i,:) / dble(nAnB(i))
            end if

        end do

    end subroutine MakeRDF

    !########################################

    subroutine ResetRDF

        implicit none

        nSample = 0
        hist    = 0.0
        gr      = 0.0

    end subroutine ResetRDF

    !########################################

    subroutine PrintRDF

        implicit none

        integer, parameter :: UNIT_RDF = 10
        integer :: i, bin

        open(UNIT_RDF, file='RDF_CenCen.txt')

            do bin = 1, maxBin - 1

                write(UNIT_RDF,"(2F20.6)") (dble(bin - 1) + 0.5)*delR, gr(0,bin)

            end do

        close(UNIT_RDF)

        do i = 1, mxCmb

            open(UNIT_RDF, file='RDF_'//trim(molAmolB(i))//'.txt')

            do bin = 1, maxBin - 1

                write(UNIT_RDF,"(2F20.6)") (dble(bin - 1) + 0.5)*delR, gr(i,bin)

            end do

            close(UNIT_RDF)

        end do

    end subroutine PrintRDF

    !########################################

    subroutine RDFDeallocate ! destructor

        implicit none

        if (allocated(atomName))   deallocate(atomName)
        if (allocated(atomLabel))  deallocate(atomLabel)
        if (allocated(hist))       deallocate(hist)
        if (allocated(gr))         deallocate(gr)
        if (allocated(cate))       deallocate(cate)
        if (allocated(cateName))   deallocate(cateName)
        if (allocated(repetition)) deallocate(repetition)
        if (allocated(nAnB))       deallocate(nAnB)
        if (allocated(molAmolB))   deallocate(molAmolB)

        write(*,"('RDF array is deallocated')")
        write(1,"('RDF array is deallocated')")

    end subroutine RDFDeallocate

    !########################################

end module mod_RDF_core
