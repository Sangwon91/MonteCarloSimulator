!------------------------------------------------------------------------
!subroutine GetCoordinate(filename)
!subroutine GetMoleculeBody()
!subroutine GetLatticePoint()
!
!subroutine CalMaxRcut(maxR)
!subroutine CalculateGroupCenter(status_)
!subroutine AccVolume(reset_)
!subroutine UpdateGroupCenter(i)
!subroutine UpdatePosition(i)
!subroutine UpdateCoordinate(i)
!subroutine GroupOldNew(i)
!
!subroutine MakeConfigView(filename,viewType)
!subroutine MakeLatticePoint()
!subroutine MakeGlobalGroupVariable()
!
!subroutine PrintConfigInf(unit_)
!subroutine PrintDensity(unit_)
!subroutine PrintLatticeInf(mat,unit_)
!subroutine PrintGroupCoordinate()
!
!subroutine WriteCoordinate(filename, dhbox)
!
!subroutine ConfigDeallocate()
!------------------------------------------------------------------------

module mod_config_information

! mod_config_information module have configuration information
! 1. box shape matrix, inverse matrix
! 2. molecule center of mass configuration
! 3. molecule orientation
! 4. indivisual atoms configuration

    use mod_math_tool         ,   only : Inverse, Det  ! to use matrix inversion, volume calculation
    use mod_system_information,   only : molNum_c => molNum, &   ! to use molecule number and etc, and volume
                                         simulTypeENUM, EXP_SOLUTE ! for coordinate reading
    use mod_molecule_information, only : nAtoms_c => nAtoms, color_c => color, &
                                         fRx => molCrdX, fRy => molCrdy, fRz => molCrdZ, molMass, &
                                         myType, nLocalGroup, localGroup, mass
    use mod_filesystem
    ! to use number of atoms and etc

    implicit none

    private

    integer, parameter :: r8 = selected_real_kind(15,300)

    logical :: newViewType = .true.

    integer, parameter :: UNIT_CONFIG = 10, CL = 10 ! character length
    character(len=CL)  :: crystalType, coordinateType
    logical :: configReadQ   ! input configuration이 있으면 true of false
    real(kind=r8) :: volume, newVolume, densKgM3, nacVol, acvolA3, maxRcut
    real(kind=r8), dimension(3,3), target :: hBox, ihBox, newhBox, newihBox
    real(kind=r8), allocatable, dimension(:,:), target :: rx, ry, rz, rxNew, ryNew, rzNew
    real(kind=r8), allocatable, dimension(:), target   :: rxiOld, ryiOld, rziOld, &
                                                          rxiNew, ryiNew, rziNew, &
                                                          ltRx , ltRy, ltRz
    real(kind=r8), dimension(3,3) :: acBox, avBox
    real(kind=r8) :: nacBox


    integer :: nGlobalGroup

    real(kind=r8), dimension(:), allocatable, target :: grpX, grpY, grpZ, &
                                    grpXNew,  grpYNew,  grpZNew,  &
                                    grpiXOld, grpiYOld, grpiZOld, &
                                    grpiXNew, grpiYNew, grpiZNew

    integer, dimension(:), allocatable :: myGroupType, myMolIndex
    integer, dimension(:,:), allocatable :: globalGroupIndex

!   ---- expanded ensemble module for liquid ----

    integer :: soluteIndex

!   ---------------------------------------------

    public :: nGlobalGroup, grpX, grpY, grpZ, &
                grpXNew,  grpYNew,  grpZNew,  &
                grpiXOld, grpiYOld, grpiZOld, &
                grpiXNew, grpiYNew, grpiZNew, &
                myGroupType, globalGroupIndex, myMolIndex

    public :: rX, rY, rZ, hBox, ihBox, newhBox, newihBox, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew, &
              volume, rxNew, ryNew, rzNew, newVolume, configReadQ, ltRx , ltRy, ltRz, acBox, avBox, nacBox

    public :: soluteIndex

    public :: GetCoordinate, PrintConfigInf, WriteCoordinate, MakeConfigView, MakeLatticePoint, &
              GetMoleculeBody, GetLatticePoint, ConfigDeallocate, PrintLatticeInf, PrintDensity, AccVolume, &
              CalculateGroupCenter, UpdateGroupCenter, PrintGroupCoordinate, UpdateCoordinate, UpdatePosition, &
              GroupOldNew

    contains

    !********************

    subroutine GetCoordinate (filename)

        use mod_random_tool, only : ranf

        integer :: i, j, molMax
        logical :: singularQ = .false.
        character(len=31), intent(in) :: filename
        real(kind=r8) :: ranx, rany, ranz

        open(UNIT_CONFIG, file=filename)

            ! expanded ensemble initialize
            soluteIndex = molNum_c

            ! memory deallocation
            if (allocated(rX)   )  deallocate(rX,rY,rZ)
            if (allocated(rxiOld)) deallocate(rxiOld,ryiOld,rziOld)
            if (allocated(rxiNew)) deallocate(rxiNew,ryiNew,rziNew)

            ! memory allocation
            allocate(rx(molNum_c,0:nAtoms_c), ry(molNum_c,0:nAtoms_c), rz(molNum_c,0:nAtoms_c), &
                     rxNew(molNum_c,0:nAtoms_c), ryNew(molNum_c,0:nAtoms_c), rzNew(molNum_c,0:nAtoms_c))
            allocate(rxiOld(0:nAtoms_c),ryiOld(0:nAtoms_c),rziOld(0:nAtoms_c),&
                     rxiNew(0:nAtoms_c),ryiNew(0:nAtoms_c),rziNew(0:nAtoms_c) )

            ! read crystal type, cubic, fcc, bcc...
            read(UNIT_CONFIG,*) crystalType

            ! read lattice information
            hBox = 0.0_r8 ! box initialize
            do i = 1, 3

                read(UNIT_CONFIG,*) hBox(i,1:3)

            end do

            ! make inverse matric and volume of system in number density unit.
            call Inverse(hBox, ihBox, singularQ, 3) ! make inverse matrix
            volume = Det (hBox,3)
            if (singularQ) stop 'singular box shape'

            ! configuration read test.
            read(UNIT_CONFIG,*) configReadQ
            read(uNIT_CONFIG,*) coordinateType

            ! configuration reading...
            if (configReadQ) then

                select case(coordinateType)
                    ! never used...
                    case('center','CENTER','Center') ! only center of mass

                        do i=1, molNum_c
                            read(UNIT_CONFIG,*) rx(i,0), ry(i,0), rz(i,0)
                        end do

                    case('atom','ATOM','Atom')       ! center of mass and all atoms configuration

                        molMax = molNum_c

                        if (simulTypeENUM == EXP_SOLUTE) then

                            if (molMax == 1) stop 'molecule number is too small (mol number == 1)'
                            molMax = molMax - 1 ! neglect solute coordinate.

                        end if

                        do i = 1, molMax

                            do j = 0, nAtoms_c

                                read(UNIT_CONFIG,*) rx(i,j), ry(i,j), rz(i,j)

                            end do

                        end do

                        ! for solute, use index 1 configuration with random translation.
                        ! equilibiriation step is needed.
                        if (simulTypeENUM == EXP_SOLUTE) then

                            ranx = 2.0_r8 * ranf() - 1.0_r8 ! [-0.5,0.5]
                            rany = 2.0_r8 * ranf() - 1.0_r8 ! [-0.5,0.5]
                            ranz = 2.0_r8 * ranf() - 1.0_r8 ! [-0.5,0.5]

                            rx(soluteIndex,:) = rx(1,:) - rx(1,0) &
                                              + hBox(1,1) * ranx + hBox(1,2) * rany + hBox(1,3) * ranz
                            ry(soluteIndex,:) = ry(1,:) - ry(1,0) &
                                              + hBox(2,1) * ranx + hBox(2,2) * rany + hBox(2,3) * ranz
                            rz(soluteIndex,:) = rz(1,:) - rz(1,0) &
                                              + hBox(3,1) * ranx + hBox(3,2) * rany + hBox(3,3) * ranz

                            ! random position in the simulation box

                        end if

                    case default

                        stop 'unclassfied configration type'

                end select

            end if

            acBox = 0.0_r8; nacBox = 0.0_r8

        close(UNIT_CONFIG)

        call MakeGlobalGroupVariable

        call CalMaxRcut(maxRcut)

        call copy_file(filename, 'backup_'//filename)
        ! ----------------------------

    end subroutine GetCoordinate

    !********************

    subroutine PrintConfigInf (unit_)

        implicit none

        integer, optional, intent(in) :: unit_

        if (present(unit_)) then

            write(unit_,"(A/)") 'configuration information ->'

            write(unit_,"(A,A10)")   'crystal type    = ', crystalType
            write(unit_,"(A/)")      'lattice vector  = '
            write(unit_,"(3F10.4//,3F10.4//,3F10.4/)") transpose(hBox)
            write(unit_,"(A,F15.4)") 'maximum rcut    = ', maxRcut
            write(unit_,"(A,F15.4)") 'system volume   = ', volume
            write(unit_,"(A,F15.4)") 'system density  = ', dble(molNum_c) / volume
            write(unit_,"(A,L10)")   'config read?    = ', configReadQ
            write(unit_,"(A,A10)")   'coordinate type = ', coordinateType
            write(unit_,"()")

        else

            write(*,"(A/)") 'configuration information ->'

            write(*,"(A,A10)")       'crystal type    = ', crystalType
            write(*,"(A/)")          'lattice vector  = '
            write(*,"(3F10.4//,3F10.4//,3F10.4/)") transpose(hBox)
            write(*,"(A,F15.4)") 'maximum rcut    = ', maxRcut
            write(*,"(A,F15.4)") 'system volume   = ', volume
            write(*,"(A,F15.4)") 'system density  = ', dble(molNum_c) / volume
            write(*,"(A,L10)")       'config read?    = ', configReadQ
            write(*,"(A,A10)")       'coordinate type = ', coordinateType
            write(*,"()")

        end if

    end subroutine PrintConfigInf

    !********************
    subroutine AccVolume (reset_)

        implicit none

        logical, optional, intent(in) :: reset_
        logical, save                 :: ini_

        if (ini_) then
          acvolA3  = 0.0
          nacVol   = 0.0
          ini_     = .false.
        end if

        acvolA3  = acvolA3 + Det (hBox,3)
        nacVol   = nacVol  + 1.0

        if (present(reset_)) then

          if (reset_) then
            acvolA3  = 0.0
            nacVol   = 0.0
          end if

        end if

    end subroutine AccVolume
    !********************
    subroutine PrintDensity (unit_)

        implicit none

        integer, optional, intent(in) :: unit_

        real(kind=r8) :: volA3 ! volume 1/Angstrom^3
        real(kind=r8) :: densNA3, densgcm3

        real(kind=r8), parameter :: dTocm3 = 1.66054 ! density to mole per cm^3

        densNA3 = real(molNum_c)/(acvolA3 / nacVol)
        densgcm3 = densNA3 * dTocm3 * molMass

        if (present(unit_)) then

          write(unit_,"('dens = ', (F15.6, ' 1/A^3  ;'),(F15.6, ' g/cm^3'))") densNA3,  densgcm3

        else

          write(*,"('dens = ', (F15.6, ' 1/A^3   ;'), (F15.6, ' g/cm^3'))") densNA3,  densgcm3

        end if

    end subroutine PrintDensity
    !********************

    subroutine PrintLatticeInf (mat,unit_)

        use mod_univ_const, only : RAD_TO_DEG

        implicit none

        integer, optional, intent(in) :: unit_
        real(kind=r8), dimension(3,3), intent(in) :: mat
        real(kind=r8) norm1, norm2, norm3, ang12, ang23, ang31

        norm1 = sum(mat(:,1)*mat(:,1))**0.5_r8
        norm2 = sum(mat(:,2)*mat(:,2))**0.5_r8
        norm3 = sum(mat(:,3)*mat(:,3))**0.5_r8

        ang12 = acos(sum(mat(:,1)*mat(:,2))/norm1/norm2)*RAD_TO_DEG
        ang23 = acos(sum(mat(:,2)*mat(:,3))/norm2/norm3)*RAD_TO_DEG
        ang31 = acos(sum(mat(:,3)*mat(:,1))/norm3/norm1)*RAD_TO_DEG

        if (present(unit_)) then

            write(unit_,"(6A10)") 'norm1','norm2','norm3','ang12','ang23','ang31'
            write(unit_,"(6F10.4)") norm1,norm2,norm3,ang12,ang23,ang31
            write(unit_,"()")

        else

            write(*,"(6A10)") 'norm1','norm2','norm3','ang12','ang23','ang31'
            write(*,"(6F10.4)") norm1,norm2,norm3,ang12,ang23,ang31
            write(*,"()")

        end if

    end subroutine PrintLatticeInf

    !********************

    subroutine WriteCoordinate (filename, dhbox)

        !use mod_math_tool, only : Inverse

        implicit none

        logical :: singularQ
        character(len=31), intent(in) :: filename
        integer, parameter :: UNIT_COOR = 10
        integer :: i, j
        real(kind=r8)                 :: ux, uy, uz
        real(kind=r8), dimension(3,3), intent(inout) :: dhbox
        real(kind=r8), dimension(:,:), allocatable :: dRx, dRy, dRz

        allocate(dRx(molNum_c,0:nAtoms_c), dRy(molNum_c,0:nAtoms_c), dRz(molNum_c,0:nAtoms_c))

        !dhbox = hbox
        dhbox(2,1) = 0.0_r8
        dhbox(3,1) = 0.0_r8
        dhbox(3,2) = 0.0_r8

        !call Inverse (dhbox, dihbox, singularQ, 3)

        do i = 1, molNum_c

            ux = ihbox(1,1) * rx(i,0) + ihbox(1,2) * ry(i,0) + ihbox(1,3) * rz(i,0)
            uy =                        ihbox(2,2) * ry(i,0) + ihbox(2,3) * rz(i,0)
            uz =                                               ihbox(3,3) * rz(i,0)

            dRx(i,0) = dhbox(1,1) * ux + dhbox(1,2) * uy + dhbox(1,3) * uz
            dRy(i,0) =                   dhbox(2,2) * uy + dhbox(2,3) * uz
            dRz(i,0) =                                     dhbox(3,3) * uz

            dRx(i,:) = rx(i,:) - rx(i,0) + dRx(i,0)
            dRy(i,:) = ry(i,:) - ry(i,0) + dRy(i,0)
            dRz(i,:) = rz(i,:) - rz(i,0) + dRz(i,0)

        end do

        open(UNIT_COOR,file=filename)

            write(UNIT_COOR,"(A/)") crystalType
            do i=1,3; write(UNIT_COOR,"(3F10.5/)") dhBox(i,:); end do

            write(UNIT_COOR,"(L10)") .true.
            write(uNIT_COOR,"(A10)") coordinateType

            select case(coordinateType)

                case('center','CENTER','Center') ! only center of mass

                    do i=1, molNum_c
                        write(UNIT_COOR,"(3F20.5)") drx(i,0), dry(i,0), drz(i,0)
                    end do

                case('atom','ATOM','Atom')       ! center of mass and all atoms configuration

                    do i=1, molNum_c

                        do j=0, nAtoms_c

                            write(UNIT_COOR,"(3F20.5)") drx(i,j), dry(i,j), drz(i,j)

                        end do

                    end do

                case default

                    stop 'unclassfied configration type'

            end select

        close(UNIT_COOR)

        deallocate(dRx,dRy,dRz)

    end subroutine WriteCoordinate

    !********************

    subroutine MakeConfigView (filename, viewType)

        use mod_math_tool, only : Det

        implicit none

        integer, parameter :: UNIT_VIEW = 10
        integer            :: i, j, index_ = 0;
        character(len=31), intent(in) :: filename
        character(len=10), optional, intent(in) :: viewType
        character(len=10) :: dummyType

        if(.not.present(viewType)) then
            dummyType = 'atom'
        else
            dummyType = viewType
        end if

        open(UNIT_VIEW,file=filename)

            write(UNIT_VIEW,*) molNum_c * nAtoms_c

            if (newViewType) then

                do i = 1, 3
                    write(UNIT_VIEW,"(3F20.5)") hBox(i,1), hBox(i,2), hBox(i,3)
                end do

            end if

            select case(dummyType)

                case('center','CENTER','Center') ! only center of mass

                    do i=1, molNum_c
                        write(UNIT_VIEW,"(4F20.5)") 1, rx(i,0), ry(i,0), rz(i,0)
                    end do

                case('atom','ATOM','Atom')       ! center of mass and all atoms configuration

                    do i = 1, molNum_c - 1

                        do j = 1, nAtoms_c

                            index_ = index_ + 1
                            write(UNIT_VIEW,"(4F20.5)") color_c(myType(j)), rx(i,j), ry(i,j), rz(i,j)

                        end do

                            !index_ = index_ + 1
                            !write(UNIT_VIEW,"(4F20.5)") color_c(myType(1)), rx(i,5), ry(i,5), rz(i,5)

                    end do

                    if (simulTypeENUM == EXP_SOLUTE) then

                        do j = 1, nAtoms_c

                            index_ = index_ + 1
                            write(UNIT_VIEW,"(4F20.5)") 8, rx(i,j), ry(i,j), rz(i,j)

                        end do

                    else

                        do j = 1, nAtoms_c

                            index_ = index_ + 1
                            write(UNIT_VIEW,"(4F20.5)") color_c(myType(j)), rx(i,j), ry(i,j), rz(i,j)

                        end do

                    end if

                case default

                    stop 'unclassfied configration type'

            end select

            write(UNIT_VIEW,*) 0
            write(UNIT_VIEW,*) Det (hBox,3)**(1.0_r8/3.0_r8)

        close(UNIT_VIEW)

    end subroutine MakeConfigView

    !********************

    subroutine MakeLatticePoint

        implicit none

        integer :: i, j, ix, iy, iz, iref, m, cellNum, nc
        real(kind=r8) :: cell, cell2, di, realNC
        real(kind=r8), allocatable, dimension(:) :: dLtRx, dLtRy, dLtRz ! dummy lattice point

        select case (crystalType)

            case ('FCC','fcc')

                cellNum = 4

                nc      = molNum_c / cellNum
                realNC  = real(nc) ** (1.0_r8/3.0_r8)
                nc      = nint(realnc)

                cell    = 1.0_r8 / realNC
                cell2   = 0.5_r8 * cell

                if (allocated(ltRx)) deallocate(ltRx,ltRy,ltRz)

                allocate(ltRx(molNum_c),ltRy(molNum_c),ltRz(molNum_c))

                ltRx(1) = 0.0_r8
                ltRy(1) = 0.0_r8
                ltRz(1) = 0.0_r8

                ltRx(2) = cell2
                ltRy(2) = cell2
                ltRz(2) = 0.0_r8

                ltRx(3) = 0.0_r8
                ltRy(3) = cell2
                ltRz(3) = cell2

                ltRx(4) = cell2
                ltRy(4) = 0.0_r8
                ltRz(4) = cell2

            case ('BCC','bcc')

                cellNum = 2

                nc      = molNum_c / cellNum
                realNC  = real(nc) ** (1.0_r8/3.0_r8)
                nc      = nint(realnc)

                cell    = 1.0_r8 / realNC
                cell2   = 0.5_r8 * cell

                if (allocated(ltRx)) deallocate(ltRx,ltRy,ltRz)

                allocate(ltRx(molNum_c),ltRy(molNum_c),ltRz(molNum_c))

                ltRx(1) = 0.0_r8
                ltRy(1) = 0.0_r8
                ltRz(1) = 0.0_r8

                ltRx(2) = cell2
                ltRy(2) = cell2
                ltRz(2) = cell2

            case ('cubic','CUBIC')

                cellNum = 1

                nc      = molNum_c / cellNum
                realNC  = real(nc) ** (1.0_r8/3.0_r8)
                nc      = nint(realnc)
                !if (abs(real(nc) - molNum_c / cellNum) > 10.0e-3) then
                !  nc = nc + 1
                !  realNC = realNC + 1.0
                !end if

                cell    = 1.0_r8 / realNC
                cell2   = 0.5_r8 * cell

                if (allocated(ltRx)) deallocate(ltRx,ltRy,ltRz)

                allocate(ltRx(molNum_c),ltRy(molNum_c),ltRz(molNum_c))

                ltRx(1) = 0.0_r8
                ltRy(1) = 0.0_r8
                ltRz(1) = 0.0_r8

            case default

                stop 'unclassfied crystal type'

        end select

        m = 0

outer : do iz=1,nc; do iy=1,nc; do ix=1,nc

            do iref=1, cellNum

                ltRx(iref+m) = ltRx(iref) + cell * real(ix-1)
                ltRy(iref+m) = ltRy(iref) + cell * real(iy-1)
                ltRz(iref+m) = ltRz(iref) + cell * real(iz-1)

            end do

            m = m + cellNum
            if (m >= molNum_c) exit outer

        end do; end do; end do outer

        di = cell * real(nc-1) + cell2
        di = 0.5_r8 * di

        ltRx = ltRx - di
        ltRy = ltRy - di
        ltRz = ltRz - di

        if (allocated(dLtRx)) deallocate(dLtRx,dLtRy,dLtRz)

        allocate(dltRx(molNum_c),dltRy(molNum_c),dltRz(molNum_c))

        dltRx = ltRx
        dltRy = ltRy
        dltRz = ltRz

        ltRx = hBox(1,1) * dltRx + hBox(1,2) * dltRy + hBox(1,3) * dltRz
        ltRy = hBox(2,1) * dltRx + hBox(2,2) * dltRy + hBox(2,3) * dltRz
        ltRz = hBox(3,1) * dltRx + hBox(3,2) * dltRy + hBox(3,3) * dltRz

        deallocate(dLtRx,dLtRy,dLtRz)

    end subroutine MakeLatticePoint

    !********************

    subroutine GetMoleculeBody()

        implicit none

        integer :: i

        rX(:,0) = 0.0_r8
        rY(:,0) = 0.0_r8
        rZ(:,0) = 0.0_r8

        forall(i = 1 : nAtoms_c) ! add molecule frame to lattice points

            rX(:,i) = fRx(i)
            rY(:,i) = fRy(i)
            rZ(:,i) = fRz(i)

        end forall

    end subroutine GetMoleculeBody

    !********************

    subroutine GetLatticePoint()

        implicit none

        integer :: i

        do i = 1, molNum_c
            rX(i,:) = rX(i,:) - rX(i,0) + ltRx(i)
            rY(i,:) = rY(i,:) - rY(i,0) + ltRy(i)
            rZ(i,:) = rZ(i,:) - rZ(i,0) + ltRz(i)
        end do

    end subroutine GetLatticePoint

    !********************

    subroutine CalMaxRcut (maxR)

        implicit none

        real(kind=r8), intent(out) :: maxR
        real(kind=r8) :: dvol
        real(kind=r8):: ax, bx, by, cx, cy, cz
        real(kind=r8), dimension(3) :: hh

        dvol = Det (hBox,3)


        ax = hbox(1,1)

        bx = hbox(1,2)
        by = hbox(2,2)

        cx = hbox(1,3)
        cy = hbox(2,3)
        cz = hbox(3,3)

        hh(1) = dvol / abs(ax*by)
        hh(2) = dvol / sqrt(abs(ax*cy)**2 + abs(ax*cz)**2)
        hh(3) = dvol / sqrt(abs(-(by*cx) + bx*cy)**2 + abs(bx*cz)**2 + abs(by*cz)**2)

        maxR = 0.5 * minval(hh)

    end subroutine CalMaxRcut

!   -------------------------------------------------
    subroutine MakeGlobalGroupVariable()
!   -------------------------------------------------
        implicit none

        integer :: mol, ng, index
!   -------------------------------------------------
        nGlobalGroup = molNum_c * nLocalGroup

        allocate(globalGroupIndex(molNum_c,nLocalGroup))
        allocate(myGroupType(nGlobalGroup), myMolIndex(nGlobalGroup))
        allocate(grpX(nGlobalGroup), grpY(nGlobalGroup), grpZ(nGlobalGroup))
        allocate(grpXNew(nGlobalGroup), grpYNew(nGlobalGroup), grpZNew(nGlobalGroup))
        allocate(grpiXOld(nLocalGroup), grpiYOld(nLocalGroup), grpiZOld(nLocalGroup))
        allocate(grpiXNew(nLocalGroup), grpiYNew(nLocalGroup), grpiZNew(nLocalGroup))

        index = 0
        do mol = 1, molNum_c

            do ng = 1, nLocalGroup

                index = index + 1
                globalGroupIndex(mol,ng) = index
                myGroupType(index) = ng   ! is needed?
                myMolIndex(index)  = mol  ! is needed? yes, is needed.

            end do ! ng-loop

        end do ! mol-loop

    end subroutine MakeGlobalGroupVariable
!   -------------------------------------------------

!   -------------------------------------------------
    subroutine PrintGroupCoordinate()
!   -------------------------------------------------
        implicit none

        integer, parameter :: UNIT_ = 10
        integer :: mol, ng, ig
!   -------------------------------------------------
        open(unit=UNIT_, file='group_coordinate.txt')

        do mol = 1, molNum_c

            do ng = 1, nLocalGroup

                ig = globalGroupIndex(mol,ng)
                write(UNIT_,"('mol: ', I3, ', group: ', I3, ', r = <', 3F10.3, ' >')") mol, ng, grpx(ig), grpy(ig), grpz(ig)

            end do ! ng-loop

        end do ! mol-loop

        close(UNIT_)

    end subroutine PrintGroupCoordinate

!   -------------------------------------------------------------------------------------
    subroutine CalculateGroupCenter(status_) ! i is optional
                                             ! if i exgist, calculate grp of molecule i only.
!   -------------------------------------------------------------------------------------
        implicit none

        integer, intent(in) :: status_

        logical :: allQ = .false.
        integer :: im, mol, ng, ig, nMem, m, ia, iat ! index atom type
        integer, dimension(:), pointer :: mem

        real(kind=r8) :: gmass, rxm, rym, rzm, imass ! group mass, position*mass
!   ------------------------------------------------------------------------------------

        if (status_ == 0) allQ = .true.
        if (status_ /= 0) allQ = .false.

        im = status_

        do mol = 1, molNum_c
            ! neglect when non-overall calculation and mol != i case
            if(.not.allQ .and. mol /= im) cycle

            do ng = 1, nLocalGroup

                ig  =  globalGroupIndex(mol,ng)
                mem => localGroup(ng)%member

                nMem  = mem(0)
                gmass = 0.0_r8
                rxm   = 0.0_r8
                rym   = 0.0_r8
                rzm   = 0.0_r8

                do m = 1, nMem

                    ia  = mem(m)
                    iat = myType(ia)
                    imass = mass(iat)
                    gmass = gmass + imass
                    rxm = rxm + imass * rx(mol,ia)
                    rym = rym + imass * ry(mol,ia)
                    rzm = rzm + imass * rz(mol,ia)

                end do ! m-loop

                grpX(ig) = rxm / gmass
                grpY(ig) = rym / gmass
                grpZ(ig) = rzm / gmass
                !write(*,"('mol: ', I3, ', group: ', I3, ', r = <', 3F10.3, ' >')") mol, ng, grpx(ig), grpy(ig), grpz(ig)

            end do ! ng-loop

        end do ! mol-loop

    end subroutine CalculateGroupCenter

!   -----------------------------------------------------------------------------------
    subroutine UpdateGroupCenter(i) ! i is optional
                                    ! if i exgist, update grp of molecule i only.
!   -----------------------------------------------------------------------------------
        implicit none

        integer, intent(in) :: i

        logical :: allQ = .true.
        integer :: mol, ng, ig, nMem, m, ia, iat ! index atom type
        integer, dimension(:), pointer :: mem
!   -----------------------------------------------------------------------------------

        if (i == 0) allQ = .true.
        if (i /= 0) allQ = .false.

        if (allQ) then

            grpX = grpXNew ! parallel calculation
            grpY = grpYNew
            grpZ = grpZNew

        else

            do ng = 1, nLocalGroup

                ig  =  globalGroupIndex(i,ng)
                grpX(ig) = grpiXNew(ng)
                grpY(ig) = grpiYNew(ng)
                grpZ(ig) = grpiZNew(ng)

            end do ! ng-loop

        end if

    end subroutine UpdateGroupCenter

!   -----------------------------------------------------------------------------------
    subroutine UpdatePosition(i) ! i is optional
                                   ! if i exgist, update position of molecule i only.
!   -----------------------------------------------------------------------------------
        implicit none

        integer, intent(in) :: i

        logical :: allQ = .true.

!   -------------------------------------------------

        if (i == 0) allQ = .true.
        if (i /= 0) allQ = .false.

        if (allQ) then

            rx = rxNew
            ry = ryNew
            rz = rzNew

        else

            rx(i,:) = rxinew(:)
            ry(i,:) = ryinew(:)
            rz(i,:) = rzinew(:)

        end if

    end subroutine UpdatePosition

!   -----------------------------------------------------------------------------------
    subroutine UpdateCoordinate(i) ! i is optional
                                   ! if i exgist, update position of molecule i only.
!   -----------------------------------------------------------------------------------
        implicit none

        integer, intent(in) :: i

        logical :: allQ = .true.
!   -----------------------------------------------------------------------------------

        if (i == 0) allQ = .true.
        if (i /= 0) allQ = .false.

!        if (i == 50) &
!          write(*,*) "im 50..."

        if (allQ) then

            call UpdatePosition(0)
            call UpdateGroupCenter(0)

        else
!            write(*,*) 'update routine'
            call UpdatePosition(i)
            call UpdateGroupCenter(i)

        end if

    end subroutine UpdateCoordinate

!   -----------------------------------------------------------------------------------
    subroutine GroupOldNew(i) ! make Group old new coordinate given rx, ry, rz
!   -----------------------------------------------------------------------------------
        implicit none

        integer, intent(in) :: i

        logical :: allQ
        integer :: im, ng, ig, nMem, m, ia, iat ! index atom type
        integer, dimension(:), pointer :: mem

        real(kind=r8) :: gmass, rxm, rym, rzm, imass ! group mass, position*mass
        real(kind=r8) :: rxi0, ryi0, rzi0, rxiNew0, ryiNew0, rziNew0
!   -----------------------------------------------------------------------------------

        if (i == 0) allQ = .true.
        if (i /= 0) allQ = .false.

        if (allQ) then  ! change only translational order

            do im = 1, molNum_c

                rxi0    = rx(im,0)   ! save old center
                ryi0    = ry(im,0)
                rzi0    = rz(im,0)

                rxiNew0 = rxNew(im,0) ! save new center
                ryiNew0 = ryNew(im,0)
                rziNew0 = rzNew(im,0)

                do ng = 1, nLocalGroup

                    ig = globalGroupIndex(im,ng)

                    grpxNew(ig) = grpX(ig) - rxi0 + rxiNew0 ! put new position using rNew..
                    grpyNew(ig) = grpY(ig) - ryi0 + ryiNew0
                    grpzNew(ig) = grpZ(ig) - rzi0 + rziNew0

                end do ! ng-loop

            end do ! im-loop

        else

            do ng = 1, nLocalGroup

                ig = globalGroupIndex(i,ng)

                grpixOld(ng) = grpX(ig)
                grpiyOld(ng) = grpY(ig)
                grpizOld(ng) = grpZ(ig)

                mem => localGroup(ng)%member

                nMem  = mem(0)
                gmass = 0.0_r8
                rxm   = 0.0_r8
                rym   = 0.0_r8
                rzm   = 0.0_r8

                do m = 1, nMem

                    ia  = mem(m)
                    iat = myType(ia)
                    imass = mass(iat)
                    gmass = gmass + imass
                    rxm = rxm + imass * rxiNew(ia)
                    rym = rym + imass * ryiNew(ia)
                    rzm = rzm + imass * rziNew(ia)

                end do ! m-loop

                grpixNew(ng) = rxm / gmass
                grpiyNew(ng) = rym / gmass
                grpizNew(ng) = rzm / gmass

            end do ! ng-loop

        end if

    end subroutine GroupOldNew

!   -----------------------------------------------
    subroutine ConfigDeallocate()
!   -----------------------------------------------
        implicit none
!   -----------------------------------------------


        !write(*,*) 'ConfigDeallocate'

        if (allocated(rX)    ) deallocate(rX,rY,rZ)
        if (allocated(rxNew) ) deallocate(rxNew,ryNew,rzNew)
        if (allocated(rxiOld)) deallocate(rxiOld,ryiOld,rziOld)
        if (allocated(rxiNew)) deallocate(rxiNew,ryiNew,rziNew)
        if (allocated(ltRx))   deallocate(ltRx,ltRy,ltRz)
        if (allocated(globalGroupIndex)) deallocate(globalGroupIndex)
        if (allocated(myGroupType))      deallocate(myGroupType)
        if (allocated(myMolIndex))      deallocate(myMolIndex)
        if (allocated(grpX))             deallocate(grpX, grpY, grpZ)
        if (allocated(grpXNew))          deallocate(grpXNew, grpYNew, grpZNew)
        if (allocated(grpiXOld))   deallocate(grpiXOld, grpiYOld, grpiZOld)
        if (allocated(grpiXNew))   deallocate(grpiXNew, grpiYNew, grpiZNew)

        write(*,"('Config array is deallocated')")
        write(1,"('Config array is deallocated')")

    end subroutine ConfigDeallocate

    !********************

end module mod_config_information

!================================================================
