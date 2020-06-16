!------------------------------------------------------------------------------------
!subroutine GetEinTransInf(filename)
!subroutine GetEinRotInf(filename)
!subroutine GetEinRotEuler(filename)
!subroutine GetEinRotLinear(filename)
!
!subroutine EinEnergy(i, rxOld_, ryOld_, rzOld_, deEin)
!subroutine EinSumup(deEin)
!
!subroutine MolCenterToEinLattice()
!
!subroutine AccEin(i, reset)
!subroutine AcEinLattice(i, reset)
!subroutine AcEinFrame(i, reset_)
!subroutine AcEinEuler(i, reset)  ! 1: Euler angle
!subroutine AcEinLinear(i, reset) ! 1: Linear by Vega
!
!subroutine FindEulerAngle(i)
!subroutine FindPsi(i, rxCur, ryCur, rzCur)
!
!subroutine MakeInitialFrame()
!subroutine MakeBodyFrame(index_, rxCur, ryCur, rzCur)
!subroutine MakeAngleHist()
!subroutine WriteAngleProb()
!
!subroutine PrintEinInf(unit_)
!subroutine WriteEinInf(filename)
!subroutine WriteEinTransInf(filename)
!subroutine WriteEinRotInf(filename)
!subroutine WriteEinEuler(filename)
!subroutine WriteEinLinear(filename)
!
!subroutine EinDeallocate
!------------------------------------------------------------------------------------


module mod_einstein_core

    use mod_system_information,   only : molNum_i => molNum, latticeType => einLatticeType
    use mod_molecule_information, only : nAtoms_i => nAtoms, &
                                         fRx_m => molCrdX, fRy_m => molCrdy, fRz_m => molCrdZ, &
                                         fA1_i => frameAtom1, fA2_i => frameAtom2, linearQ => linearQ, &
                                         flexQ
    use mod_config_information,   only : rX_ein => rX, rY_ein => rY, rZ_ein => rZ
    !use mod_movement_information

    implicit none

    private

    integer, parameter :: r8 = selected_real_kind(15,300)

    logical :: check_prob_Q = .false.

    real(kind=r8), allocatable, dimension(:)   :: einX, einY, einZ, nAcframe, nAcLattice, &
                                                  acEinX, acEinY, acEinZ, histOmega, histPhi, histTheta, histPsi
    real(kind=r8), dimension(3)                :: frameXI, frameYI, frameZI
    real(kind=r8), allocatable, dimension(:,:) :: frameX, frameY, frameZ, &
                                                  frameX0, frameY0, frameZ0, &
                                                  acFrameX, acFrameY, acFrameZ, &
                                                  frameA, frameB, frameA0, frameB0, acFrameA, acFrameB

    real(kind=r8) :: einTrans, einTheta, einOmega, einPsi, einTransK, einThetaK, einOmegaK, einPsiK, theta, phi, psi, &
                     eEin, eEinTot, eEinOld, eEinNew, deEin, delOme, delPhi, delThe, delPsi, cosPsi

    public :: eEinTot, eEinOld, eEinNew, deEin !, cosPsi !eEin ! cosPsi 나중에 지워
    public :: GetEinTransInf, GetEinRotInf, WriteEinTransInf, WriteEinRotInf,&
              AcEinFrame, AcEinLattice, EinEnergy, EinSumup, MakeInitialFrame, &
              MakeBodyFrame, FindEulerAngle, EinDeallocate, PrintEinInf, FindPsi, &
              MakeAngleHist, WriteAngleProb, MolCenterToEinLattice, &
              AccEin, WriteEinInf, GetEinInf, &
              get_orien_factor, &
              WriteOrientationFactor, &
              TurnOnCheckProb, &
              InitCheckProb, &
              ResetAngleHist

    contains

    !********************

    subroutine GetEinInf()

        implicit none

        character(len=31) :: filename

        filename = 'EinTransInf.txt'
        call GetEinTransInf(filename)
        filename = 'EinRotInf.txt'
        call GetEinRotInf(filename)

    end subroutine GetEinInf

    subroutine GetEinTransInf (filename)

        use mod_univ_const, only : RGAS

        implicit none

        character(len=31) :: filename

        integer :: i, mLattice
        integer, parameter :: UNIT_EIN = 10

        open(UNIT_EIN, file=filename)

            read(UNIT_EIN,*) einTrans

            einTransK = einTrans / RGAS

            if (allocated(einX)) deallocate(einX, einY, einZ)

            if (latticeType == 'mol' .or. latticeType == 'MOL') then
                mLattice = molNum_i
            else if (latticeType == 'atom' .or. latticeType == 'ATOM') then  ! for all atom spring
                mLattice = molNum_i * nAtoms_i
            else
                write(*,*) 'lattice type error at mod_einstein_core'
                call system("pause")
            end if

            allocate(einX(mLattice), einY(mLattice), einZ(mLattice))

            do i = 1, mLattice
                read(UNIT_EIN,*) einX(i), einY(i), einZ(i)
            end do

        close(UNIT_EIN)

    end subroutine GetEinTransInf

    !********************

    subroutine GetEinRotInf (filename)

        implicit none

        character(len=31) :: filename

        if (flexQ) return             ! return for flexible molecule

        if (linearQ)      call GetEinRotLinear (filename)
        if (.not.linearQ) call GetEinRotEuler (filename)

    end subroutine GetEinRotInf

    subroutine GetEinRotEuler (filename)

        use mod_univ_const, only : RGAS

        implicit none

        character(len=31) :: filename

        integer :: i
        integer, parameter :: UNIT_EIN = 10

        open(UNIT_EIN, file=filename)

            read(UNIT_EIN,*) einTheta
            read(UNIT_EIN,*) einOmega

            if (check_prob_Q) then
              ! so as not to skip frame reading process.
              einTheta = 1.e3 ! meaningless value bigger than 1.e-6.
              einOmega = 1.e3 ! meaningless value bigger than 1.e-6.
            end if

            einThetaK = einTheta / RGAS
            einOmegaK = einOmega / RGAS

            if(allocated(frameX)) deallocate(frameX, frameY, frameZ, frameX0, frameY0, frameZ0 )
            allocate(frameX(molNum_i,3), frameY(molNum_i,3), frameZ(molNum_i,3), &
                     frameX0(molNum_i,3), frameY0(molNum_i,3), frameZ0(molNum_i,3))

            if (abs(einTheta) < 1.e-6 .and. abs(einOmega) < 1.e-6) then
                close(UNIT_EIN)
                return
            end if

            do i = 1, molNum_i

                read(UNIT_EIN,*) frameX0(i,1), frameX0(i,2), frameX0(i,3)
                read(UNIT_EIN,*) frameY0(i,1), frameY0(i,2), frameY0(i,3)
                read(UNIT_EIN,*) frameZ0(i,1), frameZ0(i,2), frameZ0(i,3)

            end do

        close(UNIT_EIN)

    end subroutine GetEinRotEuler

    subroutine GetEinRotLinear (filename)

        use mod_univ_const, only : RGAS

        implicit none

        character(len=31) :: filename

        integer :: i
        integer, parameter :: UNIT_EIN = 10

        open(UNIT_EIN, file=filename)

            read(UNIT_EIN,*) einPsi

            einPsiK = einPsi / RGAS

            if(allocated(frameA)) deallocate(frameA, frameA0)
            allocate(frameA(molNum_i,3), frameA0(molNum_i,3))

            if (abs(einPsi) < 1.e-7) then
                close(UNIT_EIN)
                return
            end if

            do i = 1, molNum_i

                read(UNIT_EIN,*) frameA0(i,1), frameA0(i,2), frameA0(i,3)

            end do

        close(UNIT_EIN)

        do i = 1, molNum_i

            frameA0(i,:) = frameA0(i,:) / sqrt(frameA0(i,1)**2 + frameA0(i,2)**2 + frameA0(i,3)**2)

        end do

    end subroutine GetEinRotLinear

    !********************

    subroutine PrintEinInf (unit_)

        implicit none

        integer, optional, intent(in) :: unit_

        if (present(unit_)) then

            write(unit_,"(A/)") 'einstein molecule information =>'
            write(unit_,"(4A25)") 'trans energy[kJ/mol/A^2]', 'rot energy theta[kJ/mol]', 'omega[kJ/mol]', 'psi[kJ/mol]'
            write(unit_,"(4F25.6)") einTrans, einTheta, einOmega, einPsi
            write(unit_,"(4A25)") 'trans energy[K/A^2]', 'rot energy theta[K]', 'omega[K]', 'psi[K]'
            write(unit_,"(4F25.6)") einTransK, einThetaK, einOmegaK, einPsiK
            write(unit_,"()")

        else

            write(*,"(A/)") 'einstein molecule information =>'
            write(*,"(4A25)") 'trans energy[kJ/mol/A^2]', 'rot energy theta[kJ/mol]', 'omega[kJ/mol]', 'psi[kJ/mol]'
            write(*,"(4F25.6)") einTrans, einTheta, einOmega, einPsi
            write(*,"(4A25)") 'trans energy[K/A^2]', 'rot energy theta[K]', 'omega[K]', 'psi[K]'
            write(*,"(4F25.6)") einTransK, einThetaK, einOmegaK, einPsiK
            write(*,"()")

        end if

    end subroutine PrintEinInf

    !********************

    subroutine AccEin(i, reset)

        implicit none

        integer, optional, intent(in) :: i
        logical, optional, intent(in) :: reset

        call AcEinLattice(i, reset)
        call AcEinFrame(i, reset)

    end subroutine AccEin

    !********************

    subroutine AcEinLattice (i, reset)

        implicit none

        integer, optional, intent(in) :: i
        logical, optional, intent(in) :: reset

        logical, save :: iniQ = .true.
        integer, save :: mLattice
        integer :: start_, end_, j, at, idx

        if (iniQ) then

            if (latticeType == 'mol' .or. latticeType == 'MOL') then
                mLattice = molNum_i
            else if (latticeType == 'atom' .or. latticeType == 'ATOM') then  ! for all atom spring
                mLattice = molNum_i * nAtoms_i
            else
                write(*,*) 'lattice type error at mod_einstein_core'
                call system("pause")
            end if

            iniQ = .false.

            if (.not.allocated(nAcLattice)) allocate(nAcLattice(mLattice))
            if (.not.allocated(acEinX))     allocate(acEinX(mLattice), acEinY(mLattice), acEinZ(mLattice))

        end if


        if (present(i) .and. (.not.present(reset))) then

            if (i /= 0) then       ! for specific molecule i
                start_ = i
                end_   = i
            else if (i == 0) then  ! for all molecules in system
                start_ = 1
                end_   = molNum_i
            end if


            if (molNum_i == mLattice) then   ! for mol base lattice

                do j = start_, end_

                    acEinX(j) = acEinX(j) + rX_ein(j,0)
                    acEinY(j) = acEinY(j) + rY_ein(j,0)
                    acEinZ(j) = acEinZ(j) + rZ_ein(j,0)

                    nAcLattice(j) = nAcLattice(j) + 1.0_r8
                end do

            end if

            if (molNum_i < mLattice) then   ! for atom base lattice

                do j = start_, end_

                    do at = 1, nAtoms_i

                    idx = (j - 1) * nAtoms_i + at   ! linear mapping

                    acEinX(idx) = acEinX(idx) + rX_ein(j,at)
                    acEinY(idx) = acEinY(idx) + rY_ein(j,at)
                    acEinZ(idx) = acEinZ(idx) + rZ_ein(j,at)

                    nAcLattice(idx) = nAcLattice(idx) + 1.0_r8

                    end do ! at-loop

                end do ! j-loop

            end if

        end if

        if (present(reset)) then

            if (reset) then

                acEinX = 0.0_r8
                acEinY = 0.0_r8
                acEinZ = 0.0_r8

                nAcLattice = 0.0_r8

            end if

        end if

    end subroutine AcEinLattice
    !********************

    subroutine AcEinFrame (i, reset_)

        implicit none

        integer, optional, intent(in) :: i
        logical, optional, intent(in) :: reset_

        if (flexQ) return    ! not valid for flexible body

        if(present(i)) then

            if (linearQ)       call AcEinLinear (i)
            if (.not.linearQ)  call AcEinEuler (i)

        end if

        if(present(reset_)) then

            if (linearQ)       call AcEinLinear (1, reset_)
            if (.not.linearQ)  call AcEinEuler (1, reset_)

        end if

    end subroutine AcEinFrame

    subroutine AcEinEuler (i, reset) ! 1: Euler angle

        implicit none

        integer, optional, intent(in) :: i
        logical, optional, intent(in) :: reset

        if (.not.allocated(acFrameX)) allocate(acFrameX(molNum_i,0:nAtoms_i), &
                acFrameY(molNum_i,0:nAtoms_i), acFrameZ(molNum_i,0:nAtoms_i))
        if (.not.allocated(nAcframe)) allocate(nAcframe(molNum_i))

        if (present(i)) then

            acFrameX(i,:)  = acFrameX(i,:) + (rX_ein(i,:) - rX_ein(i,0))
            acFrameY(i,:)  = acFrameY(i,:) + (rY_ein(i,:) - rY_ein(i,0))
            acFrameZ(i,:)  = acFrameZ(i,:) + (rZ_ein(i,:) - rZ_ein(i,0))
            nAcFrame(i)    = nAcFrame(i) + 1.0_r8

        end if

        if (present(reset)) then

            if ( reset ) then

                acFrameX = 0.0_r8
                acFrameY = 0.0_r8
                acFrameZ = 0.0_r8
                nAcFrame = 0.0_r8

            end if

        end if

    end subroutine AcEinEuler

    subroutine AcEinLinear (i, reset) ! 1: Linear by Vega

        implicit none

        integer, optional, intent(in) :: i
        logical, optional, intent(in) :: reset

        if (.not.allocated(acFrameA)) allocate(acFrameA(molNum_i,3))
        if (.not.allocated(nAcframe)) allocate(nAcframe(molNum_i))

        if (present(i)) then

            acFrameA(i,1) = acFrameA(i,1) + (rX_ein(i,2) - rX_ein(i,1))
            acFrameA(i,2) = acFrameA(i,2) + (rY_ein(i,2) - rY_ein(i,1))
            acFrameA(i,3) = acFrameA(i,3) + (rZ_ein(i,2) - rZ_ein(i,1))
            nAcFrame(i)   = nAcFrame(i) + 1.0_r8

        end if

        if (present(reset)) then

            if ( reset ) then

                acFrameA = 0.0_r8
                nAcFrame = 0.0_r8

            end if

        end if

    end subroutine AcEinLinear

    !********************

    subroutine WriteEinInf()

        character(len=31) :: filename

        filename = 'EinTransInf.txt'
        call WriteEinTransInf (filename)
        filename = 'EinRotInf.txt'
        call WriteEinRotInf (filename)

    end subroutine WriteEinInf

    subroutine WriteEinTransInf (filename)

        implicit none

        integer :: i, mLattice
        character(len=31), intent(in) :: filename
        integer, parameter :: UNIT_EIN = 10

        if (.not.allocated(acEinX)) return

        if (latticeType == 'mol' .or. latticeType == 'MOL') then
            mLattice = molNum_i
        else if (latticeType == 'atom' .or. latticeType == 'ATOM') then  ! for all atom spring
            mLattice = molNum_i * nAtoms_i
        else
            write(*,*) 'lattice type error at mod_einstein_core'
            call system("pause")
        end if

        open(UNIT_EIN, file=filename)

            write(UNIT_EIN,"(F30.6/)") einTrans

            do i = 1, mLattice
                write(UNIT_EIN, "(3F30.6)") acEinX(i) / nAcLattice(i), acEinY(i) / nAcLattice(i), &
                                            acEinZ(i) / nAcLattice(i)
            end do

        close(UNIT_EIN)

    end subroutine WriteEinTransInf

    !********************

    subroutine WriteEinRotInf (filename)

        implicit none

        character(len=31), intent(in) :: filename

        if (linearQ)      call WriteEinLinear (filename)
        if (.not.linearQ) call WriteEinEuler (filename)

    end subroutine WriteEinRotInf

    subroutine WriteEinEuler (filename)

        implicit none

        integer :: i
        character(len=31), intent(in) :: filename
        integer, parameter :: UNIT_EIN = 10
        real(kind=r8), allocatable, dimension(:) :: avX, avY, avZ

        if (.not.allocated(acFrameX)) return

        allocate(avX(0:nAtoms_i), avY(0:nAtoms_i), avZ(0:nAtoms_i))
        if (.not.allocated(frameX)) &
        allocate(frameX(molNum_i,3), frameY(molNum_i,3), frameZ(molNum_i,3))

        open(UNIT_EIN, file=filename)

            write(UNIT_EIN,"(F30.6)") einTheta
            write(UNIT_EIN,"(F30.6/)") einOmega

            if (.not.flexQ) then

                do i = 1, molNum_i

                    !write(*,*) 'test  set number = ', i

                    avX = acFrameX(i,:) / nAcFrame(i)
                    avY = acFrameY(i,:) / nAcFrame(i)
                    avZ = acFrameZ(i,:) / nAcFrame(i)

                    call MakeBodyFrame(i, avX, avY, avZ)

                    !write(*,*) 'write start'; call system("pause")

                    write(UNIT_EIN,"(3F30.6)") frameX(i,1), frameX(i,2), frameX(i,3)
                    write(UNIT_EIN,"(3F30.6)") frameY(i,1), frameY(i,2), frameY(i,3)
                    write(UNIT_EIN,"(3F30.6)") frameZ(i,1), frameZ(i,2), frameZ(i,3)

                end do

            end if ! flexQ

        close(UNIT_EIN)

        deallocate(avX, avY, avZ)

    end subroutine WriteEinEuler

    subroutine WriteEinLinear (filename)

        implicit none

        integer :: i
        character(len=31), intent(in) :: filename
        integer, parameter :: UNIT_EIN = 10
        real(kind=r8)      :: avAX, avAY, avAZ

        if (.not.allocated(acFrameA)) return

        open(UNIT_EIN, file=filename)

            write(UNIT_EIN,"(F30.6)") einPsi

            if (.not.flexQ) then

                do i = 1, molNum_i

                    avAX = acFrameA(i,1) / nAcFrame(i)
                    avAY = acFrameA(i,2) / nAcFrame(i)
                    avAZ = acFrameA(i,3) / nAcFrame(i)

                    write(UNIT_EIN,"(3F30.6)") avAX, avAY, avAZ

                end do

            end if ! flexQ

        close(UNIT_EIN)

    end subroutine WriteEinLinear

    !********************

    subroutine EinEnergy (i, rxOld_, ryOld_, rzOld_, deEin)

        implicit none

        integer, intent(in) :: i

        real(kind=r8), intent(out) :: deEin
        real(kind=r8), dimension(0:nAtoms_i), intent(in) :: rxOld_, ryOld_, rzOld_
        integer :: idx, at
        real(kind=r8) :: rsq_

        deEin = 0.0_r8

        if (.not.flexQ) then

            if (.not.linearQ) then  ! only for rigid body

                if (einThetaK > 1.e-5 .or. einOmegaK > 1.e-5) then

                    call MakeBodyFrame(i, rxOld_, ryOld_, rzOld_)
                    call FindEulerAngle(i)
                    deEin = deEin + einThetaK * (1.0_r8 - Cos(theta)) &
                                  + einOmegaK * (1.0_r8 - Cos(phi + psi))

                end if

            end if ! non-linear

            if (linearQ) then

                if (einPsiK > 1.e-5) then

                    call FindPsi (i, rxOld_, ryOld_, rzOld_)

                    deEin = deEin + einPsiK * (1.0_r8 - cos(2.0 * psi))

                end if

            end if ! linearQ

        end if ! flexQ

        if (latticeType == 'mol' .or. latticeType == 'MOL') then

            rsq_ = (rxOld_(0) - einX(i)) ** 2 &
                 + (ryOld_(0) - einY(i)) ** 2 &
                 + (rzOld_(0) - einZ(i)) ** 2

        else if (latticeType == 'atom' .or. latticeType == 'ATOM') then  ! for all atom spring

            rsq_ = 0.0

            do at = 1, nAtoms_i

                idx = (i - 1) * nAtoms_i + at   ! linear mapping

                rsq_ = rsq_ + (rxOld_(at) - einX(idx)) * (rxOld_(at) - einX(idx)) &
                            + (ryOld_(at) - einY(idx)) * (ryOld_(at) - einY(idx)) &
                            + (rzOld_(at) - einZ(idx)) * (rzOld_(at) - einZ(idx))

            end do

        else

            write(*,*) 'ein energy error at mod_einstein_core'
            call system("pause")

        end if

        deEin = deEin + einTransK * rsq_

        if (deEin < 0.0_r8) stop 'ein energy error'

    end subroutine EinEnergy

    !********************
    subroutine EinSumup (deEin)

        implicit none

        integer :: i

        real(kind=r8), intent(out) :: deEin
        real(kind=r8) :: rsq_, del

        deEin = 0.0_r8

        do i = 1, molNum_i

            call EinEnergy(i, rx_ein(i,:), ry_ein(i,:), rz_ein(i,:), del)
            deEin = deEin + del

        end do

        if (deEin < 0.0_r8) stop 'ein energy error'

    end subroutine EinSumup

    !********************

    subroutine MakeInitialFrame()

        implicit none

        real(kind=r8), dimension(3) :: dX, dY, dZ ! d : dummy
        real(kind=r8)               :: norm

        dX = (/fRx_m(fA1_i),fRy_m(fA1_i),fRz_m(fA1_i)/)
        dY = (/fRx_m(fA2_i),fRy_m(fA2_i),fRz_m(fA2_i)/)

        dZ(1) = dX(2)*dY(3) - dX(3)*dY(2)
        dZ(2) = dX(3)*dY(1) - dX(1)*dY(3)
        dZ(3) = dX(1)*dY(2) - dX(2)*dY(1)

        norm = (dZ(1)**2 + dZ(2)**2 + dZ(3)**2) ** 0.5_r8

        frameZI = dZ / norm

        dY(1) = dZ(2)*dX(3) - dZ(3)*dX(2)
        dY(2) = dZ(3)*dX(1) - dZ(1)*dX(3)
        dY(3) = dZ(1)*dX(2) - dZ(2)*dX(1)

        norm = (dY(1)**2 + dY(2)**2 + dY(3)**2) ** 0.5_r8

        frameYI = dY / norm

        norm = (dX(1)**2 + dX(2)**2 + dX(3)**2) ** 0.5_r8

        frameXI = dX / norm

    end subroutine MakeInitialFrame

    !********************

    subroutine MakeBodyFrame (index_, rxCur, ryCur, rzCur) ! Cur : current

        integer, intent(in) :: index_

        real(kind=r8), dimension(0:nAtoms_i), optional, intent(in) :: rxCur, ryCur, rzCur
        real(kind=r8), dimension(3) :: dX, dY, dZ ! d : dummy
        real(kind=r8)               :: norm, origin

        !write(*,*) "i'm, in"; call system("pause")

        if (present(rxCur).and.present(ryCur).and.present(rzCur)) then

            dX = (/rXCur(fA1_i)-rXCur(0),rYCur(fA1_i)-rYCur(0), rZCur(fA1_i)-rZCur(0)/)
            dY = (/rXCur(fA2_i)-rXCur(0),rYCur(fA2_i)-rYCur(0), rZCur(fA2_i)-rZCur(0)/)

        else

            dX = (/rX_ein(index_,fA1_i)-rX_ein(index_,0),rY_ein(index_,fA1_i)-rY_ein(index_,0), &
                   rZ_ein(index_,fA1_i)-rZ_ein(index_,0)/)
            dY = (/rX_ein(index_,fA2_i)-rX_ein(index_,0),rY_ein(index_,fA2_i)-rY_ein(index_,0), &
                   rZ_ein(index_,fA2_i)-rZ_ein(index_,0)/)

        end if

        dZ(1) = dX(2)*dY(3) - dX(3)*dY(2)
        dZ(2) = dX(3)*dY(1) - dX(1)*dY(3)
        dZ(3) = dX(1)*dY(2) - dX(2)*dY(1)

        norm = (dZ(1)**2 + dZ(2)**2 + dZ(3)**2) ** 0.5_r8

        frameZ(index_,:) = dZ / norm

        dY(1) = dZ(2)*dX(3) - dZ(3)*dX(2)
        dY(2) = dZ(3)*dX(1) - dZ(1)*dX(3)
        dY(3) = dZ(1)*dX(2) - dZ(2)*dX(1)

        norm = (dY(1)**2 + dY(2)**2 + dY(3)**2) ** 0.5_r8

        frameY(index_,:) = dY / norm

        norm = (dX(1)**2 + dX(2)**2 + dX(3)**2) ** 0.5_r8

        frameX(index_,:) = dX / norm

    end subroutine MakeBodyFrame

    !********************

    subroutine FindEulerAngle (i)

        use mod_univ_const, only : TWOPI

        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(3) :: dX, dY, dZ ! d : dummy
        real(kind=r8)               :: norm

        real(kind=r8), parameter :: one_ = 0.9999999999999_r8
        real(kind=r8), dimension(3) :: frameXnew
        real(kind=r8) :: z0_z, xnew_x0, xnew_y0, xnew_x, xnew_y ! '_' = dot product

        theta = 0.0_r8
        phi   = 0.0_r8
        psi   = 0.0_r8

        z0_z = frameZ0(i,1)*frameZ(i,1) + frameZ0(i,2)*frameZ(i,2) + frameZ0(i,3)*frameZ(i,3)

        if (z0_z >  1.0_r8) z0_z =  one_
        if (z0_z < -1.0_r8) z0_z = -one_

        theta = acos(z0_z)

        frameXnew(1) = frameZ0(i,2)*frameZ(i,3) - frameZ0(i,3)*frameZ(i,2)
        frameXnew(2) = frameZ0(i,3)*frameZ(i,1) - frameZ0(i,1)*frameZ(i,3)
        frameXnew(3) = frameZ0(i,1)*frameZ(i,2) - frameZ0(i,2)*frameZ(i,1)

        norm = (frameXnew(1)**2+frameXnew(2)**2+frameXnew(3)**2)**0.5_r8

        if      (norm > 10E-10)   then; frameXnew = frameXnew / norm
        else if (z0_z > 0.999999) then; frameXnew = frameX(i,:)
        else if (z0_z < -0.999999)then; frameXnew = frameX(i,:); end if

        xnew_x0 = frameXnew(1)*frameX0(i,1)+frameXnew(2)*frameX0(i,2)+frameXnew(3)*frameX0(i,3)

        if      (xnew_x0 >  1.0_r8) then; xnew_x0 =  one_
        else if (xnew_x0 < -1.0_r8) then; xnew_x0 = -one_; end if

        phi = acos(xnew_x0)

        xnew_y0 = frameXnew(1)*frameY0(i,1)+frameXnew(2)*frameY0(i,2)+frameXnew(3)*frameY0(i,3)

        if (xnew_y0 < 0.0_r8) phi = TWOPI - phi

        xnew_x = frameXnew(1)*frameX(i,1)+frameXnew(2)*frameX(i,2)+frameXnew(3)*frameX(i,3)

        if      (xnew_x >  1.0_r8) then; xnew_x =  one_
        else if (xnew_X < -1.0_r8) then; xnew_x = -one_; end if

        psi = acos(xnew_x)

        xnew_y = frameXnew(1)*frameY(i,1)+frameXnew(2)*frameY(i,2)+frameXnew(3)*frameY(i,3)

        if (xnew_y > 0.0_r8) psi = TWOPI - psi

    end subroutine FindEulerAngle

    !********************

    subroutine FindPsi (i, rxCur, ryCur, rzCur)

        implicit none

        integer, intent(in)         :: i

        logical, save               :: normQ = .false.
        real(kind=r8), parameter    :: ONE = 0.9999999999999_r8, MINUS_ONE = -ONE
        real(kind=r8), save         :: normA, normA0
        real(kind=r8)               :: AdotA0
        real(kind=r8), dimension(3) :: dframeA
        real(kind=r8), dimension(0:nAtoms_i), optional, intent(in) :: rxCur, ryCur, rzCur

        if (present(rxCur)) then

            dframeA(1) = (rxCur(2) - rxCur(1))
            dframeA(2) = (ryCur(2) - ryCur(1))
            dframeA(3) = (rzCur(2) - rzCur(1))

        else

            dframeA(1) = (rX_ein(i,2) - rX_ein(i,1))
            dframeA(2) = (rY_ein(i,2) - rY_ein(i,1))
            dframeA(3) = (rZ_ein(i,2) - rZ_ein(i,1))

        end if

        if (.not.normQ) then

            normA  = sqrt(dframeA(1)**2 + dframeA(2)**2 + dframeA(3)**2)
            normQ  = .true.

        end if

        AdotA0 = dframeA(1)*frameA0(i,1) + dframeA(2)*frameA0(i,2) + dframeA(3)*frameA0(i,3)
        cosPsi = AdotA0 / normA
        if (cosPsi > ONE)        cosPsi = ONE
        if (cosPsi < MINUS_ONE ) cosPsi = MINUS_ONE
        psi = acos(cosPsi)
    !~
    !~    write(*,"('A    = ', 3F10.6)") dframeA
    !~    write(*,"('A0   = ', 3F10.6)") frameA0(i,:)
    !~    write(*,"('A.A0,  cos(psi) = ', 2F10.6)") AdotA0, cosPsi
    !~    write(*,"('psi  = ',  F10.6)") psi
    !~    call system('pause')

    end subroutine FindPsi

    !********************

    function get_orien_factor()

      implicit none

      integer :: i
      real(kind=r8), dimension(3) :: get_orien_factor
      real(kind=r8), dimension(3) :: sum_ = 0.0

      sum_ = 0.0

      do i = 1, 1

        call FindEulerAngle(i)
        sum_ = sum_ + [phi, theta, psi]

      end do

      get_orien_factor = sum_ !/ dble(molNum_i)

    end function get_orien_factor

    subroutine WriteOrientationFactor(step)

      implicit none

      integer, intent(in) :: step
      integer, parameter :: UNIT_FACT = 31

      open(unit = UNIT_FACT, file = "fact.txt", position = "append")

        write(UNIT_FACT, "(I10, 3F10.2)") step, get_orien_factor()

      close(UNIT_FACT)

    end subroutine WriteOrientationFactor

    subroutine InitCheckProb()

      implicit none
      character(len=31) :: filename
      filename = "eq_axis.txt"

      if (.not.check_prob_Q) &
          return

      call GetEinRotInf(filename)

    end subroutine InitCheckProb

    subroutine TurnOnCheckProb()
      implicit none
      check_prob_Q = .true.
      write(*,*) "--check-prob is turned on."
    end subroutine TurnOnCheckProb

    subroutine ResetAngleHist()

      if (check_prob_Q) then
        call MakeAngleHist(.true.)
      end if

    end subroutine ResetAngleHist

    subroutine MakeAngleHist(reset_flag)

        use mod_univ_const, only : PI, TWOPI

        implicit none

        logical, intent(in), optional :: reset_flag
        integer       :: i
        logical, save :: alloQ = .false.
        real(kind=r8), parameter :: ZERO = 0.0_r8
        real(kind=r8), parameter :: ONE = 1.0_r8
        real(kind=r8), parameter :: MAXBIN = 500.0_r8
        ! reset histogram if reset_flag exist and reset_flag == .true.

        if (.not.check_prob_Q) &
            return ! nothing but bypass.

        if (present(reset_flag)) then

            if (reset_flag .and. alloQ) then

                if (linearQ) then
                    histPsi = ZERO
                else if (.not.linearQ) then
                    histOmega = ZERO
                    histPhi   = ZERO
                    histTheta = ZERO
                    histPsi   = ZERO
                else
                    stop '[Error] MakeAngleHist() in mod_einstein_core.f90'
                end if

                write(*,*) "histgram reset is done."

            end if

        end if

        if (.not.alloQ) then

            write(*,*) "array of angle histrogram is allocated."

            if (linearQ) then

                delPsi = 2.0_r8 / MAXBIN
                allocate(histPsi(floor(-1.0/delPsi)+1:floor(1.0/delPsi)+1))

            end if


            if (.not.linearQ) then

            delOme = 2.0_r8 * TWOPI / MAXBIN
            delPhi = TWOPI / MAXBIN
            delThe = 2.0_r8 / MAXBIN
            delPsi = TWOPI / MAXBIN
            allocate(histOmega(floor(2.0_r8 * TWOPI/delOme)+1), histPhi(floor(TWOPI/delPhi)+1), &
                     histTheta(floor(-1.0/delThe)+1:floor(1.0/delThe)+1), histPsi(floor(TWOPI/delPsi)+1) )

            end if

            alloQ = .true.

        end if

        do i = 1, molNum_i

            if (linearQ) then

                call FindPsi(i) ! output : psi & cosPsi !!
                histPsi(floor(cosPsi/delPsi)+1) = histPsi(floor(cosPsi/delPsi)+1) + 1.0

            end if

            if (.not.linearQ) then

                call MakeBodyFrame(i)
                call FindEulerAngle(i)

                histOmega(floor((phi+psi)/delOme)+1)  = histOmega(floor((phi+psi)/delOme)+1) + 1.0
                histTheta(floor(cos(theta)/delThe)+1) = histTheta(floor(cos(theta)/delThe)+1) + 1.0
                histPhi(floor(phi/delPhi)+1) = histPhi(floor(phi/delPhi)+1)  + 1.0
                histPsi(floor(psi/delPsi)+1) = histPsi(floor(psi/delPsi)+1)  + 1.0

            end if

        end do

    end subroutine MakeAngleHist

    subroutine WriteAngleProb()

        use mod_univ_const, only : PI, TWOPI
        !use mod_character_tool

        implicit none

        integer, parameter :: unit_= 10
        !integer, save      :: simulStep = 0

        integer             :: i
        real(kind=r8), dimension(:), allocatable :: probPsi, probOmega, probPhi, probTheta

        if (.not.check_prob_Q) &
            return ! nothing but bypass.

          if (linearQ) then

                allocate(probPsi(floor(-1.0/delPsi)+1:floor(1.0/delPsi)+1))

                probPsi = histPsi / sum(histPsi) / delPsi

                !open(unit=unit_, file='probPsi'//IntToCha(simulStep)//'.txt')
                open(unit=unit_, file='probPsi.txt')
                    do i = floor(-1.0/delPsi)+1, floor(1.0/delPsi)+1

                        write(unit_,"(2F30.6)") (i-1)*delPsi, probPsi(i)

                    end do

                close(unit_)

                !simulStep = simulStep + 1

                deallocate(probPsi)

            end if


            if (.not.linearQ) then

                allocate(probOmega(floor(2.0_r8 * TWOPI/delOme)+1), probPhi(floor(TWOPI/delPhi)+1), &
                         probTheta(floor(-1.0/delThe)+1:floor(1.0/delThe)+1), probPsi(floor(TWOPI/delPsi)+1) )

                probOmega = histOmega / sum(histOmega) / delOme
                probPhi   = histPhi   / sum(histPhi)   / delPhi
                probPsi   = histPsi   / sum(histPsi)   / delPsi
                probTheta = histTheta / sum(histTheta) / delThe

                open(unit=unit_, file='probOmega.txt')

                    do i = 1, floor(2.0_r8 * TWOPI/delOme)+1

                        write(unit_,"(2F30.6)") (i-1)*delOme, probOmega(i)

                    end do

                close(unit_)

                open(unit=unit_, file='probPhi.txt')

                    do i = 1, floor(TWOPI/delPhi)+1

                        write(unit_,"(2F30.6)") (i-1)*delPhi, probPhi(i)

                    end do

                close(unit_)

                open(unit=unit_, file='probPsi.txt')

                    do i = 1, floor(TWOPI/delPsi)+1

                        write(unit_,"(2F30.6)") (i-1)*delPsi, probPsi(i)

                    end do

                close(unit_)

                open(unit=unit_, file='probTheta.txt')

                    do i = floor(-1.0/delThe)+1, floor(1.0/delThe)+1

                        write(unit_,"(2F30.6)") (i-1)*delThe, probTheta(i)

                    end do

                close(unit_)

                deallocate(probOmega, probPhi, probTheta, probPsi)

            end if

    end subroutine WriteAngleProb

    !********************

    subroutine MolCenterToEinLattice() ! should be changed.. for all atom lattice model

        implicit none

        integer :: i, j, idx



        if (latticeType == 'mol' .or. latticeType == 'MOL') then

            do i = 1, molNum_i

                rx_ein(i,:) = rx_ein(i,:) - rx_ein(i,0) + einX(i)    ! go to equilibrium position
                ry_ein(i,:) = ry_ein(i,:) - ry_ein(i,0) + einY(i)
                rz_ein(i,:) = rz_ein(i,:) - rz_ein(i,0) + einZ(i)

            end do

            write(*,"('molecular center of mass to Ein Lattice')")
            write(1,"('molecular center of mass to Ein Lattice')")

        else if (latticeType == 'atom' .or. latticeType == 'ATOM') then  ! for all atom spring

            do i = 1, molNum_i

                do j = 1, nAtoms_i

                    idx = (i - 1) * nAtoms_i + j  ! linear mapping

                    rx_ein(i,j) = einX(idx)   ! valid only for flexmolecule
                    ry_ein(i,j) = einY(idx)
                    rz_ein(i,j) = einZ(idx)

                end do

            end do

            write(*,"('molecular atoms to be Ein Lattice')")
            write(1,"('molecular atoms to be Ein Lattice')")

        else

            write(*,*) 'lattice rearrange error at mod_einstein_core'
            call system("pause")

        end if

    end subroutine MolCenterToEinLattice

    !********************

    subroutine EinDeallocate

        implicit none

        if (allocated(einX))        deallocate(einX, einY, einZ)
        if (allocated(frameX))      deallocate(frameX, frameY, frameZ)
        if (allocated(frameX0))     deallocate(frameX0, frameY0, frameZ0 )
        if (allocated(acFrameX))    deallocate(acFrameX, acFrameY, acFrameZ)
        if (allocated(nAcframe))    deallocate(nAcframe)
        if (allocated(nAcLattice))  deallocate(nAcLattice)
        if (allocated(acEinX))      deallocate(acEinX, acEinY, acEinZ)
        if (allocated(frameA))      deallocate(frameA)
        if (allocated(frameA0))     deallocate(frameA0)
        if (allocated(acframeA))    deallocate(acframeA)
        if (allocated(histOmega))   deallocate(histOmega, histPhi, histTheta, histPsi)
        if (allocated(histPsi))     deallocate(histPsi)

        write(*,"('Ein array is deallocated')")
        write(1,"('Ein array is deallocated')")

    end subroutine EinDeallocate

    !********************

end module mod_einstein_core
