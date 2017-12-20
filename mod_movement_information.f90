!----------------------------------------------------------------------------------
!subroutine GetMoveInf(filename)
!
!subroutine ChangeToCombo(step)
!subroutine AdjustNVTMax()
!subroutine AdjustNPTMax()
!subroutine InterOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine TransOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine RotateOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine IntraOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine BondOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine BendOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine TorsOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine ChangeBox(xNew, yNew, zNew, nBox, niBox, nVol)
!subroutine MoveBond(xt, yt, zt)
!subroutine MoveBending(xt, yt, zt)
!subroutine MoveTorsion(xt, yt, zt)
!subroutine ReCenterize(xt, yt, zt)
!subroutine SysCOMToOrigin()
!subroutine ShrinkSolute(ratio_, rxi_, ryi_, rzi_)
!subroutine ShrinkOldNew(lambdaNew, xold_, yold_, zold_, xnew_, ynew_, znew_)
!
!subroutine PrintMoveInf(unit_)
!subroutine PrintAcceptanceRatio(unit_)
!
!subroutine WriteMoveInf()
!
! non using subroutines... (will be deleted...) >>
!
!subroutine TransPeriodicOldNew (i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine ComboOldNew (i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine ComboPeriodicOldNew (i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!subroutine ChangeCubicBox (xNew, yNew, zNew, nBox, niBox, nVol)
!----------------------------------------------------------------------------------

module mod_movement_information

    use mod_math_tool         ,   only : Inverse  ! to use matrix inversion
    use mod_system_information,   only : molNum_m => molNum, preStep, &     ! to use molecule number and etc, and volume
                                         periodicQ, cubicBoxQ, fixMol1Q, fixAtom1Q, randomTorsionQ, &
                                         simulTypeENUM, EXP_SOLUTE
    use mod_molecule_information, only : nAtoms_m => nAtoms, &
                                         fRx_m => molCrdX, fRy_m => molCrdy, fRz_m => molCrdZ, &
                                         fA1_m => frameAtom1, fA2_m => frameAtom2, &
                                         nBend, Index1InBend, Index2InBend, Index3InBend, &
                                         myAtomNum, molTree, molMass, mass, myType, &
                                         index1InTorsion,  index2InTorsion, index3InTorsion, index4InTorsion, &
                                         nTorsion, rotateQ, bondQ, bendQ, torsQ, &
                                         index1InBond, index2InBond, nBond
    use mod_config_information,   only : rX_m => rX, rY_m => rY, rZ_m => rZ, &
                                         xOld_m => rxiOld, yOld_m => ryiOld, zOld_m => rziOld, &
                                         xNew_m => rxiNew, yNew_m => ryiNew, zNew_m => rziNew, &
                                         box_m => hBox, iBox_m => ihBox, GroupOldNew, CalculateGroupCenter, &
                                         soluteIndex

    ! to use number of atoms and etc

    implicit none

    private

    integer, parameter :: r8 = selected_real_kind(15,300)

    integer, parameter :: fixIntra = 50, fixTrans = 50, fixRot = 50, &
                          fixBox = 50, fixBond = 50, fixBend = 50, fixTors = 50, fixInter = 50, &
                          sequenceStep = 500

    logical, dimension(20), parameter :: monoBondQ = (/.true.,  .false., .false., .false., .false., &
    &                                                  .false., .false., .false., .true. , .false., &
    &                                                  .false., .false., .false., .false., .false., &
    &                                                  .false., .true. , .false., .false., .false./)
                                                      ! H = 1, F = 9 and Cl = 17 are mono bond atoms

    logical :: adjustNVT, adjustNPT, sequenceQ = .true.
    real(kind=r8), parameter :: ACCEPT_RATIO = 0.5_r8
    real(kind=r8) :: drMax, dRotMax, dBoxMax, nTrans, nRotate, tTrans, tRotate, nBox, tBox, & ! n : number, t : try
                     dBendMax, dTorsMax, nIntra, tIntra, nInter, tInter, nBend_, tBend, nTors, tTors, &
                     dBondMax, nBond_, tBond
    public :: drMax, dRotMax, dBoxMax, nTrans, nRotate, tTrans, tRotate, nBox, tBox, &
              dBendMax, dTorsMax, nIntra, tIntra, nInter, tInter, nBend_, tBend, nTors, tTors, &
              sequenceQ, dBondMax, nBond_, tBond
    public :: GetMoveInf, PrintMoveInf, RotateOldNew, &
              TransOldNew, ComboOldNew, ChangeBox, SysCOMToOrigin, &
              ChangeCubicBox, TransPeriodicOldNew, ComboPeriodicOldNew, &
              AdjustNVTMax, AdjustNPTMax, WriteMoveInf, &
              IntraOldNew, ChangeToCombo, InterOldNew, BondOldNew, BendOldNew, TorsOldNew, &
              PrintAcceptanceRatio, ShrinkSolute, ShrinkOldNew, RandomTranslation, &
              SpecialOldNew

    contains

    !********************

    subroutine GetMoveInf (filename)

        use mod_math_tool, only : Det

        implicit none

        character(len=31), intent(in) :: filename
        integer, parameter :: UNIT_MOVE = 10
        integer :: i
        real(kind=r8) :: oot = 1.0_r8/3.0_r8 ! one over three

        open(UNIT_MOVE,file=filename)

            read(UNIT_MOVE,*) drMax
            read(UNIT_MOVE,*) dRotMax
            read(UNIT_MOVE,*) dBoxMax
            read(UNIT_MOVE,*) dBondMax
            read(UNIT_MOVE,*) dBendMax
            read(UNIT_MOVE,*) dTorsMax

        close(UNIT_MOVE)

        write(*,"('if no initial movement data, do auto adjustment'/)")

        if (drMax   < 1.e-10)  drMax    = 0.004_r8    ! 0.4 %  of box length
        if (dRotMax < 1.e-10)  dRotMax  = 0.314_r8    ! 20 degree
        if (dBoxMax < 1.e-10)  dBoxMax  = Det (box_m,3)**oot * 0.001_r8
        if (dBondMax < 1.e-10) dBondMax = 0.01_r8   ! 0.01 Å
        if (dBendMax < 1.e-10) dBendMax = 0.07_r8   ! 4 degree
        if (dTorsMax < 1.e-10) dTorsMax = 0.26_r8   ! 15 degree

        sequenceQ = .true.

        nTrans = 0.0_r8
        nRotate = 0.0_r8
        nBox   = 0.0_r8
        nInter = 0.0_r8
        nIntra = 0.0_r8
        nBond_ = 0.0_r8
        nBend_ = 0.0_r8
        nTors  = 0.0_r8

        tTrans = 0.0_r8
        tRotate = 0.0_r8
        tBox   = 0.0_r8
        tInter = 0.0_r8
        tIntra = 0.0_r8
        tBond  = 0.0_r8
        tBend  = 0.0_r8
        tTors  = 0.0_r8

        adjustNVT = .true.
        adjustNPT = .true.

    end subroutine GetMoveInf

    !********************

    subroutine PrintMoveInf (unit_)

        implicit none
        integer, intent(in), optional :: unit_

        if (present(unit_)) then

            write(unit_,"(A/)") 'movement information ->'
            write(unit_,"(A,F15.6)") 'dr max        = ', drMax
            write(unit_,"(A,F15.6)") 'drotation max = ', dRotMax
            write(unit_,"(A,F15.6)") 'dbox max      = ', dBoxMax
            write(unit_,"(A,F15.6)") 'dbond max     = ', dBondMax
            write(unit_,"(A,F15.6)") 'dbend max     = ', dBendMax
            write(unit_,"(A,F15.6)") 'dtorsion max  = ', dTorsMax
            write(unit_,"()")

        else

            write(*,"(A/)") 'movement information ->'
            write(*,"(A,F15.6)") 'dr max        = ', drMax
            write(*,"(A,F15.6)") 'drotation max = ', dRotMax
            write(*,"(A,F15.6)") 'dbox max      = ', dBoxMax
            write(*,"(A,F15.6)") 'dbond max     = ', dBondMax
            write(*,"(A,F15.6)") 'dbend max     = ', dBendMax
            write(*,"(A,F15.6)") 'dtorsion max  = ', dTorsMax
            write(*,"()")

        end if

    end subroutine PrintMoveInf

    !********************

    subroutine WriteMoveInf()

        implicit none

        integer, parameter :: UNIT_ = 10

        open(UNIT_, file='movement.txt')

            write(UNIT_,"(F15.6)") drMax
            write(UNIT_,"(F15.6)") dRotMax
            write(UNIT_,"(F15.6)") dBoxMax
            write(UNIT_,"(F15.6)") dBondMax
            write(UNIT_,"(F15.6)") dBendMax
            write(UNIT_,"(F15.6)") dTorsMax

        close(UNIT_)

    end subroutine WriteMoveInf

    !********************

    subroutine ChangeToCombo(step)

        implicit none

        integer :: step

        if (sequenceQ) then

            if (step > sequenceStep) then
                sequenceQ = .false. ! start combo calculation
                write(*,*) '! ------------ Combo stage start ---------------'
                write(1,*) '! ------------ Combo stage start ---------------'
            end if

        else

            return

        end if
    end subroutine ChangeToCombo


    subroutine PrintAcceptanceRatio (unit_)

        integer, intent(in), optional :: unit_

        character(len=10), parameter :: DASHED = '----------', BLANK = '          '
        real(kind=r8) :: ratioTrans, ratioRot, ratioIntra, ratioInter, &
                         ratioBond, ratioBend, ratioTors, ratioBox

        ratioTrans = nTrans  / tTrans
        ratioRot   = nRotate / tRotate
        ratioBond  = nBond_  / tBond
        ratioBend  = nBend_  / tBend
        ratioTors  = nTors   / tTors

        ratioInter = nInter  / tInter
        ratioIntra = nIntra  / tIntra

        ratioBox   = nBox    / tBox

        if (present(unit_)) then

            if (sequenceQ) then

                write(unit_,"(7(A10))") DASHED, DASHED, DASHED, DASHED, DASHED, DASHED, DASHED
                write(unit_,"(7(A10))") BLANK, 'trans', 'rotate', 'box', 'bond', 'bend', 'torsion'
                write(unit_,"(7(A10))") DASHED, DASHED, DASHED, DASHED, DASHED, DASHED, DASHED
                write(unit_,"(A10, 6(F10.6))") 'Δmax', drMax, dRotMax, dBoxMax, dBondMax, dBendMax, dTorsMax
                write(unit_,"(A10, 6(F10.6))") 'acceptance', ratioTrans, ratioRot, ratioBox, ratioBond, ratioBend, ratioTors

            else ! for combo movement

                write(unit_,"(7(A10))") DASHED, DASHED, DASHED, DASHED, DASHED, DASHED, DASHED
                write(unit_,"(7(A10))") BLANK, 'inter', BLANK, 'box', 'intra', BLANK, BLANK
                write(unit_,"(7(A10))") DASHED, DASHED, DASHED, DASHED, DASHED, DASHED, DASHED
                write(unit_,"(A10, 6(F10.6))") 'Δmax', drMax, dRotMax, dBoxMax, dBondMax, dBendMax, dTorsMax
                write(unit_,"(A10, F10.6, A10, F10.6, F10.6)") 'acceptance', ratioInter, BLANK, ratioBox, ratioIntra

            end if

        else

            if (sequenceQ) then

                write(*,"(7(A10))") DASHED, DASHED, DASHED, DASHED, DASHED, DASHED, DASHED
                write(*,"(7(A10))") BLANK, 'trans', 'rotate', 'box', 'bond', 'bend', 'torsion'
                write(*,"(7(A10))") DASHED, DASHED, DASHED, DASHED, DASHED, DASHED, DASHED
                write(*,"(A10, 6(F10.6))") 'Δmax', drMax, dRotMax, dBoxMax, dBondMax, dBendMax, dTorsMax
                write(*,"(A10, 6(F10.6))") 'acceptance', ratioTrans, ratioRot, ratioBox, ratioBond, ratioBend, ratioTors

            else ! for combo movement

                write(*,"(7(A10))") DASHED, DASHED, DASHED, DASHED, DASHED, DASHED, DASHED
                write(*,"(7(A10))") BLANK, 'inter', BLANK, 'box', 'intra', BLANK, BLANK
                write(*,"(7(A10))") DASHED, DASHED, DASHED, DASHED, DASHED, DASHED, DASHED
                write(*,"(A10, 6(F10.6))") 'Δmax', drMax, dRotMax, dBoxMax, dBondMax, dBendMax, dTorsMax
                write(*,"(A10, F10.6, A10, F10.6, F10.6)") 'acceptance', ratioInter, BLANK, ratioBox, ratioIntra

            end if

        end if


    end subroutine PrintAcceptanceRatio

    !********************

    subroutine AdjustNVTMax()

        implicit none

!        logical :: bendAdjustQ, torsAdjustQ   ! not tested
!        integer, parameter :: myAdjust = 1000 ! each 20 adjustment
        integer, parameter :: times = 10
        integer, save :: step = 0
        integer :: timesStep
        real(kind=r8) :: ratioTrans, ratioRot, ratioIntra, ratioInter, &
                         ratioBond, ratioBend, ratioTors

        step = step + 1
        timesStep = times * step

        ratioTrans = nTrans  / tTrans
        ratioRot   = nRotate / tRotate
        ratioBond  = nBond_  / tBond
        ratioBend  = nBend_  / tBend
        ratioTors  = nTors   / tTors

        ratioInter = nInter  / tInter
        ratioIntra = nIntra  / tIntra

        if (step > preStep) return

        if (sequenceQ) then ! Sequence adjust

            if(mod(timesStep,fixTrans)==0) then

                if (ratioTrans > ACCEPT_RATIO) drMax = drMax*1.05_r8
                if (ratioTrans < ACCEPT_RATIO) drMax = drMax*0.95_r8
                if (drMax > 0.5_r8) drMax = 0.5_r8

                nTrans = 0.0_r8; tTrans = 0.0_r8

            end if

            if(mod(timesStep,fixRot)==0) then

                if (ratioRot > ACCEPT_RATIO) dRotMax = dRotMax*1.05_r8
                if (ratioRot < ACCEPT_RATIO) dRotMax = dRotMax*0.95_r8
                if (dRotMax > 1.57_r8) dRotMax = 1.57_r8

                nRotate = 0.0_r8; tRotate = 0.0_r8

            end if

            if(mod(timesStep,fixBond)==0) then

                if (ratioBond > ACCEPT_RATIO) dBondMax = dBondMax*1.05_r8
                if (ratioBond < ACCEPT_RATIO) dBondMax = dBondMax*0.95_r8
                !if (dBondMax > 1.57_r8) dBondMax = 1.57_r8

                nBond_ = 0.0_r8; tBond = 0.0_r8

            end if

            if(mod(timesStep,fixBend)==0) then

                if (ratioBend > ACCEPT_RATIO) dBendMax = dBendMax*1.05_r8
                if (ratioBend < ACCEPT_RATIO) dBendMax = dBendMax*0.95_r8
                if (dBendMax > 1.57_r8) dBendMax = 1.57_r8

                nBend_ = 0.0_r8; tBend = 0.0_r8

            end if

            if(mod(timesStep,fixTors)==0) then

                if (ratioTors > ACCEPT_RATIO) dTorsMax = dTorsMax*1.05_r8
                if (ratioTors < ACCEPT_RATIO) dTorsMax = dTorsMax*0.95_r8
                if (dTorsMax > 3.14_r8) dTorsMax = 3.14_r8

                nTors = 0.0_r8; tTors = 0.0_r8

            end if

        else ! Combo adjust

            if(mod(step,fixInter)==0) then

                if (ratioInter > ACCEPT_RATIO) then

                    drMax   = drMax   * 1.01_r8
                    dRotMax = dRotMax * 1.01_r8

                end if

                if (ratioInter < ACCEPT_RATIO) then

                    drMax   = drMax   * 0.99_r8
                    dRotMax = dRotMax * 0.99_r8

                end if

                if (drMax   > 0.5_r8)  drMax   = 0.5_r8
                if (dRotMax > 1.57_r8) dRotMax = 1.57_r8

                nInter = 0.0_r8; tInter = 0.0_r8

            end if
        ! intra movement
            if(mod(step,fixIntra)==0) then

                if (ratioIntra > ACCEPT_RATIO) then

                    dBondMax = dBondMax * 1.01_r8
                    dBendMax = dBendMax * 1.01_r8
                    dTorsMax = dTorsMax * 1.01_r8

                end if

                if (ratioIntra < ACCEPT_RATIO) then

                    dBondMax = dBondMax * 0.99_r8
                    dBendMax = dBendMax * 0.99_r8
                    dTorsMax = dTorsMax * 0.99_r8

                end if

                if (dBendMax > 1.57_r8) dBendMax = 1.57_r8 ! never occurs
                if (dTorsMax > 3.14_r8) dTorsMax = 3.14_r8

                nIntra = 0.0_r8; tIntra = 0.0_r8

            end if

        end if

    endsubroutine AdjustNVTMax

    !********************

    subroutine AdjustNPTMax()

        implicit none

        integer, save :: step = 0
        real(kind=r8) :: ratioBox

        step = step + 1

        ratioBox = nBox / tBox

        if (step > preStep) return

        if(mod(step,fixBox)==0) then

            if (ratioBox > ACCEPT_RATIO) dBoxMax = dBoxMax * 1.01_r8
            if (ratioBox < ACCEPT_RATIO) dBoxMax = dBoxMax * 0.99_r8

            nBox = 0.0_r8; tBox = 0.0_r8

        end if

    endsubroutine AdjustNPTMax

    !********************

!   ------------------------------------------------------------
    subroutine InterOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!   ------------------------------------------------------------
        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
!   ------------------------------------------------------------

        if (sequenceQ) stop 'sequence error'

        call TransOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
        if (rotateQ) call RotateOldNew (i, xNew_, yNew_, zNew_, xNew_, yNew_, zNew_)

        if (i == 1) then

            if (fixMol1Q) then          ! fix COM of molecule 1

                xNew_ = xNew_ - xNew_(0) + xOld_(0)
                yNew_ = yNew_ - yNew_(0) + yOld_(0)
                zNew_ = zNew_ - zNew_(0) + zOld_(0)

            else if (fixAtom1Q) then    ! fix position of atom 1 of molecule 1

                xNew_ = xNew_ - xNew_(1) + xOld_(1)
                yNew_ = yNew_ - yNew_(1) + yOld_(1)
                zNew_ = zNew_ - zNew_(1) + zOld_(1)

            end if ! fix kind

        end if ! i == 1

        call GroupOldNew(i) ! calculate group position of molecule i

    end subroutine InterOldNew
!   ------------------------------------------------------------

    subroutine RotateOldNew (i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)

        use mod_random_tool, only : Ranf

        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
        integer :: xyz ! xyz : x축 : x, y축 : y, z축 : z
        real(kind=r8), dimension(0:nAtoms_m) :: dX, dY, dZ
        real(kind=r8) :: cos_, sin_, ang, oX, oY, oZ ! dRotMax

        if (sequenceQ) then
            xOld_ = rX_m(i,:)
            yOld_ = rY_m(i,:)
            zOld_ = rZ_m(i,:)
        end if
        ! else use input old value

        oX = xOld_(0)
        oY = yOld_(0)
        oZ = zOld_(0)

        dX = xOld_ - oX
        dY = yOld_ - oY
        dZ = zOld_ - oZ

        xyz  = int(3.0_r8*Ranf() + 1.0_r8)
        ang  = dRotMax * (2.0_r8*Ranf()-1.0_r8)
        cos_ = cos(ang)
        sin_ = sin(ang)

        if (xyz == 1) then

            xNew_ = dX                     + oX
            yNew_ =     cos_*dY - sin_*dZ  + oY
            zNew_ =     sin_*dY + cos_*dZ  + oZ

        else if (xyz == 2) then

            xNew_ = cos_*dX      -sin_*dZ  + oX
            yNew_ =          dY            + oY
            zNew_ = sin_*dX      +cos_*dZ  + oZ

        else if (xyz == 3) then

            xNew_ = cos_*dX - sin_*dY      + oX
            yNew_ = sin_*dX + cos_*dY      + oY
            zNew_ =                    dZ  + oZ

        else; stop 'critical error occur in rotation'; end if

        if (sequenceQ) call GroupOldNew(i) ! calculate group position of molecule i

    end subroutine RotateOldNew

    !********************

    subroutine TransOldNew (i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)

        use mod_random_tool, only : Ranf

        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
        real(kind=r8) :: noX, noY, noZ
        real(kind=r8) :: dX, dY, dZ
        real(kind=r8) :: halfBox, Box

        xOld_ = rX_m(i,:)
        yOld_ = rY_m(i,:)
        zOld_ = rZ_m(i,:)

        dX = drMax * (2.0_r8 * Ranf() - 1.0_r8)
        dY = drMax * (2.0_r8 * Ranf() - 1.0_r8)
        dZ = drMax * (2.0_r8 * Ranf() - 1.0_r8)

        noX = box_m(1,1) * dX + box_m(1,2) * dY + box_m(1,3) * dZ
        noY = box_m(2,1) * dX + box_m(2,2) * dY + box_m(2,3) * dZ
        noZ = box_m(3,1) * dX + box_m(3,2) * dY + box_m(3,3) * dZ

        xNew_ = xOld_ + noX
        yNew_ = yOld_ + noY
        zNew_ = zOld_ + noZ

        if (periodicQ) then

            if (.not.cubicBoxQ) then

                write(*,*) 'use only cubic box in periodic condition'
                call system('pause')

            end if

            halfBox = 0.5_r8 * box_m(1,1)
            Box     = box_m(1,1)

            if      (xNew_(0) >  halfBox) then; xNew_ = xNew_ - Box
            else if (xNew_(0) < -halfBox) then; xNew_ = xNew_ + Box; end if

            if      (yNew_(0) >  halfBox) then; yNew_ = yNew_ - Box
            else if (yNew_(0) < -halfBox) then; yNew_ = yNew_ + Box; end if

            if      (zNew_(0) >  halfBox) then; zNew_ = zNew_ - Box
            else if (zNew_(0) < -halfBox) then; zNew_ = zNew_ + Box; end if

        end if

        if (sequenceQ) call GroupOldNew(i) ! calculate group position of molecule i

    end subroutine TransOldNew


    subroutine RandomTranslation(xNew_, yNew_, zNew_)

        use mod_random_tool, only : Ranf

        real(kind=r8), dimension(0:nAtoms_m), intent(inout) :: xNew_, yNew_, zNew_

        real(kind=r8) :: dX, dY, dZ
        real(kind=r8) :: noX, noY, noZ

        dX = 2.0_r8 * Ranf() - 1.0_r8
        dY = 2.0_r8 * Ranf() - 1.0_r8
        dZ = 2.0_r8 * Ranf() - 1.0_r8

        noX = box_m(1,1) * dX + box_m(1,2) * dY + box_m(1,3) * dZ
        noY = box_m(2,1) * dX + box_m(2,2) * dY + box_m(2,3) * dZ
        noZ = box_m(3,1) * dX + box_m(3,2) * dY + box_m(3,3) * dZ

        xnew_ = xnew_ - xnew_(0) + noX  ! random position in simulation box
        ynew_ = ynew_ - ynew_(0) + noY
        znew_ = znew_ - znew_(0) + noZ

    end subroutine RandomTranslation


    subroutine ChangeBox (xNew, yNew, zNew, nBox, niBox, nVol)

        use mod_random_tool, only : RanMat, Ranf
        use mod_math_tool,   only : Inverse, Det

        implicit none

        integer :: i
        real(kind=r8) :: ranValue
        real(kind=r8), dimension(molNum_m,0:nAtoms_m), intent(out) :: xNew, yNew, zNew
        real(kind=r8), dimension(3,3), intent(out) :: nBox, niBox
        real(kind=r8), intent(out) :: nVol

        logical :: singularQ

        real(kind=r8), dimension(molNum_m) :: dx, dy, dz, dx0, dy0, dz0

        if (cubicBoxQ) then
            ranValue = 2.0_r8 * Ranf() - 1.0_r8
            do i = 1, 3
                nBox(i,i) = box_m(i,i) + dBoxMax * ranValue
            end do
        else
            nBox = box_m + dBoxMax * RanMat(3,3,'u-')
        end if
        call Inverse (nBox, niBox,singularQ,3)
        nVol = Det (nBox,3)

        dx0 = rX_m(:,0)
        dy0 = rY_m(:,0)
        dz0 = rZ_m(:,0)

        dx = iBox_m(1,1)*dx0 + iBox_m(1,2)*dy0 + iBox_m(1,3)*dz0
        dy =                   iBox_m(2,2)*dy0 + iBox_m(2,3)*dz0
        dz =                                     iBox_m(3,3)*dz0

        ! upper만 된다.

        dx = nBox(1,1)*dx + nBox(1,2)*dy + nBox(1,3)*dz
        dy =                nBox(2,2)*dy + nBox(2,3)*dz
        dz =                               nBox(3,3)*dz

        do i = 0, nAtoms_m

            xNew(:,i) = rX_m(:,i) - dx0 + dx
            yNew(:,i) = rY_m(:,i) - dy0 + dy
            zNew(:,i) = rZ_m(:,i) - dz0 + dz

        end do

        call GroupOldNew(0) ! calculate group position of overall molecules

    end subroutine ChangeBox


    subroutine IntraOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)

        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_

        if (.not.bondQ .and. .not.bendQ .and. .not.torsQ) stop 'why perform intra movement?'

        xOld_ = rX_m(i,:)
        yOld_ = rY_m(i,:)
        zOld_ = rZ_m(i,:)

        xNew_ = xOld_
        yNew_ = yOld_
        zNew_ = zOld_

        if (bondQ) call MoveBond(xNew_, yNew_, zNew_)    ! molecular vaibrational motion
        if (bendQ) call MoveBending(xNew_, yNew_, zNew_) ! molecular bending motion
        if (torsQ) call MoveTorsion(xNew_, yNew_, zNew_) ! molecular torsion motion

        call ReCenterize(xNew_, yNew_, zNew_) ! reset center of mass

        if (i == 1) then

            if (fixMol1Q) then    ! fix COM of molecule 1

                xNew_ = xNew_ - xNew_(0) + xOld_(0)
                yNew_ = yNew_ - yNew_(0) + yOld_(0)
                zNew_ = zNew_ - zNew_(0) + zOld_(0)

            else if (fixAtom1Q) then    ! fix position of atom 1 of molecule 1

                xNew_ = xNew_ - xNew_(1) + xOld_(1)
                yNew_ = yNew_ - yNew_(1) + yOld_(1)
                zNew_ = zNew_ - zNew_(1) + zOld_(1)

            end if ! fix kind

        end if ! i == 1
        call GroupOldNew(i) ! calculate group position of molecule i

    end subroutine IntraOldNew

    !********************
!   -------------------------------------------------------
    subroutine BondOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!   -------------------------------------------------------
        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
!   -------------------------------------------------------

        xOld_ = rX_m(i,:)
        yOld_ = rY_m(i,:)
        zOld_ = rZ_m(i,:)

        xNew_ = xOld_
        yNew_ = yOld_
        zNew_ = zOld_

        call MoveBond(xNew_, yNew_, zNew_) ! molecular vibrational motion
        call ReCenterize(xNew_, yNew_, zNew_) ! reset center of mass
        call GroupOldNew(i) ! calculate group position of molecule i

    end subroutine

!   -------------------------------------------------------
    subroutine BendOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!   -------------------------------------------------------
        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
!   -------------------------------------------------------

        xOld_ = rX_m(i,:)
        yOld_ = rY_m(i,:)
        zOld_ = rZ_m(i,:)

        xNew_ = xOld_
        yNew_ = yOld_
        zNew_ = zOld_

        call MoveBending(xNew_, yNew_, zNew_) ! molecular bending motion
        !call MoveTorsion(xNew_, yNew_, zNew_) ! molecular torsion motion
        call ReCenterize(xNew_, yNew_, zNew_) ! reset center of mass
        call GroupOldNew(i) ! calculate group position of molecule i

    end subroutine

!   -------------------------------------------------------
    subroutine TorsOldNew(i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
!   -------------------------------------------------------
        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
!   -------------------------------------------------------

        xOld_ = rX_m(i,:)
        yOld_ = rY_m(i,:)
        zOld_ = rZ_m(i,:)

        xNew_ = xOld_
        yNew_ = yOld_
        zNew_ = zOld_

        !call MoveBending(xNew_, yNew_, zNew_) ! molecular bending motion
        call MoveTorsion(xNew_, yNew_, zNew_) ! molecular torsion motion
        call ReCenterize(xNew_, yNew_, zNew_) ! reset center of mass
        call GroupOldNew(i) ! calculate group position of molecule i

    end subroutine

!   #######################################################

    subroutine SpecialOldNew(i, oldx_, oldy_, oldz_, newx_, newy_, newz_)

      use, intrinsic :: iso_c_binding

      implicit none

      integer(kind=c_int), parameter :: max_atom = 5
      integer(kind=c_int), intent(in) :: i
      real(kind=c_double), dimension(0:max_atom), intent(out) :: &
        oldx_, oldy_, oldz_, newx_, newy_, newz_

      oldx_ = rx_m(i, :)
      oldy_ = ry_m(i, :)
      oldz_ = rz_m(i, :)

      newx_ = oldx_
      newy_ = oldy_
      newz_ = oldz_

      call MoveSpecial(newx_, newy_, newz_)
      call ReCenterize(newx_, newy_, newz_)
      call GroupOldNew(i)

    end subroutine SpecialOldNew

!   -------------------------------------------------------
    subroutine MoveBond(xt, yt, zt)
!   -------------------------------------------------------
        use mod_random_tool, only : Ranf

        implicit none

        integer, parameter :: atomMax = 100 ! for speed up of subroutine.
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xt, yt, zt ! moved atomic position.

        logical :: mono1Q, mono2Q

        integer :: i, path
        integer :: bond, bType, a1, a2, ivib, ifx, nVib, avib, iv
        real(kind=r8) :: dbond, ratio, oldBondLength
        real(kind=r8), dimension(3) :: p, dBondVec
!   -------------------------------------------------------

        bond  = int(Ranf() * nBond + 1.0) ! random index-selection for positive index.
        dbond = dBondMax * (2.0_r8 * Ranf() - 1.0_r8) !

        a1 = index1InBond(bond)
        a2 = index2InBond(bond)

        mono1Q = monoBondQ(myAtomNum(a1))
        mono2Q = monoBondQ(myAtomNum(a2))

        if (mono1Q .and. mono2Q) then ! diatomic molecule

            if (Ranf() < 0.5_r8) then
                ivib = a1
                ifx  = a2
            else
                ivib = a2
                ifx  = a1
            end if

        else if (mono1Q .and. (.not.mono2Q)) then

            ivib = a1
            ifx  = a2

        else if ((.not.mono1Q) .and. mono2Q) then

            ivib = a2
            ifx  = a1

        else

            if (Ranf() < 0.5_r8) then
                ivib = a1
                ifx  = a2
            else
                ivib = a2
                ifx  = a1
            end if

        end if

        p(1) = xt(ivib) - xt(ifx)   ! bond vector
        p(2) = yt(ivib) - yt(ifx)
        p(3) = zt(ivib) - zt(ifx)

        oldBondLength = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
        ratio = dbond / oldBondLength

        dBondVec = ratio * p  ! make vector which is (bond new) - (bond old)

        do i = 1, 4

            if (molTree(ifx,i,1) /= ivib) cycle

            path = i
            exit

        end do

        nVib = molTree(ifx,path,0)

        do i = 1, nVib

            iv = molTree(ifx,path,i)

            xt(iv) = xt(iv) + dBondVec(1)
            yt(iv) = yt(iv) + dBondVec(2)
            zt(iv) = zt(iv) + dBondVec(3)

        end do

    end subroutine MoveBond

!   -------------------------------------------------------
    subroutine MoveBending(xt, yt, zt)
!   -------------------------------------------------------
        use mod_random_tool, only : Ranf

        implicit none

        integer, parameter  :: atomMax = 100 ! for speed of subroutine
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xt, yt, zt ! moved atomic position

        logical :: mono1Q, mono3Q

        integer :: i, path, nRotated
        integer :: bend, bType, a1, a2, a3, pivot, ir, ifx ! ir: index roatated, ifx: index fixed
        real(kind=r8) :: angle, cosa, sina, norm, ndr  ! ndr: n dot r
        real(kind=r8), dimension(3) :: p, nv, rcn, rt  ! nv: normal vector, rcn: r cross n, rt: rotated r
                                                       ! rt: relative position of rotated
        integer, dimension(atomMax) :: irot
        real(kind=r8), dimension(atomMax,3) :: r
!   -------------------------------------------------------

        bend  = int(Ranf() * nBend + 1.0)
        angle = dBendMax * (2.0_r8 * Ranf() - 1.0_r8) ! which is correct for detailed valance..?
        cosa  = cos(angle)
        sina  = sin(angle)

        a1 = Index1InBend(bend)    ! get atom index in bending(index of i = bend)
        a2 = Index2InBend(bend)
        a3 = Index3InBend(bend)

        pivot = a2

        mono1Q = monoBondQ(myAtomNum(a1)) ! true for mono bonding atom (H, F or Cl ...)
        mono3Q = monoBondQ(myAtomNum(a3))

        if (mono1Q .and. mono3Q) then ! select ramdomly

            if (Ranf() < 0.5_r8) then
                ir  = a1 ! roatated
                ifx = a3
            else
                ir  = a3
                ifx = a1
            end if

        else if (mono1Q .and. (.not.mono3Q)) then ! rotate a1 path

            ir  = a1
            ifx = a3

        else if ((.not.mono1Q) .and. mono3Q) then ! rotate a3 path

            ir  = a3
            ifx = a1

        else ! select ramdomly

            if (Ranf() < 0.5_r8) then
                ir  = a1 ! roatated
                ifx = a3
            else
                ir  = a3
                ifx = a1
            end if

        end if

        p(1) = xt(ifx) - xt(pivot)
        p(2) = yt(ifx) - yt(pivot)
        p(3) = zt(ifx) - zt(pivot)

        do i = 1, 4

            if (molTree(pivot,i,1) /= ir) cycle

            path = i
            exit

        end do

        nRotated = molTree(pivot,path,0)

        do i = 1, nRotated

            irot(i) = molTree(pivot,path,i)
            ir      = irot(i)
            r(i,1) = xt(ir) - xt(pivot)  ! relative position of atoms rotated
            r(i,2) = yt(ir) - yt(pivot)
            r(i,3) = zt(ir) - zt(pivot)

        end do

        nv(1) = r(1,2) * p(3) - r(1,3) * p(2)
        nv(2) = r(1,3) * p(1) - r(1,1) * p(3)
        nv(3) = r(1,1) * p(2) - r(1,2) * p(1)
        norm = dsqrt( nv(1)*nv(1) + nv(2)*nv(2) + nv(3)*nv(3) )
        nv(1) = nv(1) / norm
        nv(2) = nv(2) / norm
        nv(3) = nv(3) / norm

        do i = 1, nRotated

            ndr = nv(1) * r(i,1) + nv(2) * r(i,2) + nv(3) * r(i,3)
            rcn(1) = r(i,2) * nv(3) - r(i,3) * nv(2)
            rcn(2) = r(i,3) * nv(1) - r(i,1) * nv(3)
            rcn(3) = r(i,1) * nv(2) - r(i,2) * nv(1)
            rt(1) = nv(1)*ndr + ( r(i,1)-nv(1)*ndr ) * cosa + rcn(1) * sina
            rt(2) = nv(2)*ndr + ( r(i,2)-nv(2)*ndr ) * cosa + rcn(2) * sina
            rt(3) = nv(3)*ndr + ( r(i,3)-nv(3)*ndr ) * cosa + rcn(3) * sina
            ir = irot(i)
            xt(ir) = xt(pivot) + rt(1)
            yt(ir) = yt(pivot) + rt(2)
            zt(ir) = zt(pivot) + rt(3)

        end do ! nRotated

    end subroutine MoveBending

!   -------------------------------------------------------
    subroutine MoveTorsion(xt, yt, zt)
!   -------------------------------------------------------
        use mod_random_tool, only : Ranf

        implicit none

        integer, parameter  :: atomMax = 100
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xt, yt, zt ! moved atomic position

        integer :: i, path, nRotated
        integer :: tors, tType, a1, a2, a3, a4, pivot, ir, ifx ! ir: index roatated, ifx: index fixed
        real(kind=r8) :: angle, cosa, sina, norm, ndr  ! ndr: n dot r
        real(kind=r8), dimension(3) :: p, nv, rcn, rt  ! nv: normal vector, rcn: r cross n, rt: rotated r
                                                       ! rt: relative position of rotated
        integer, dimension(atomMax) :: irot
        real(kind=r8), dimension(atomMax,3) :: r
!   -------------------------------------------------------

        angle = dTorsMax * (2.0_r8 * Ranf() - 1.0_r8)
        if (randomTorsionQ) then
            if (ranf() > 0.95_r8) angle = angle * 9.0_r8   ! random torsional motion
        end if
        cosa  = cos(angle)
        sina  = sin(angle)
        tors  = int(Ranf() * nTorsion + 1.0)

        a1 = index1InTorsion(tors)
        a2 = index2InTorsion(tors)
        a3 = index3InTorsion(tors)
        a4 = index4InTorsion(tors)

        if (Ranf() < 0.5) then

            pivot = a3
            ir    = a4
            ifx   = a2

        else

            pivot = a2
            ir    = a1
            ifx   = a3

        end if

        nv(1) = xt(ifx) - xt(pivot)
        nv(2) = yt(ifx) - yt(pivot)
        nv(3) = zt(ifx) - zt(pivot)
        norm = dsqrt( nv(1)*nv(1) + nv(2)*nv(2) + nv(3)*nv(3) )
        nv(1) = nv(1) / norm
        nv(2) = nv(2) / norm
        nv(3) = nv(3) / norm

        do i = 1, 4
            if (molTree(ifx,i,1) /= pivot) cycle
            path = i
            exit
        end do

        nRotated = molTree(ifx,path,0)

        do i = 2, nRotated        ! i = 1 for pivot

            irot(i) = molTree(ifx,path,i)
            ir = irot(i)
            r(i,1) = xt(ir) - xt(pivot)
            r(i,2) = yt(ir) - yt(pivot)
            r(i,3) = zt(ir) - zt(pivot)

        end do

!     // Rotate atoms
        do i = 2, nRotated

            ndr = nv(1) * r(i,1) + nv(2) * r(i,2) + nv(3) * r(i,3)
            rcn(1) = r(i,2) * nv(3) - r(i,3) * nv(2)
            rcn(2) = r(i,3) * nv(1) - r(i,1) * nv(3)
            rcn(3) = r(i,1) * nv(2) - r(i,2) * nv(1)
            rt(1) = nv(1)*ndr + ( r(i,1)-nv(1)*ndr ) * cosa + rcn(1) * sina
            rt(2) = nv(2)*ndr + ( r(i,2)-nv(2)*ndr ) * cosa + rcn(2) * sina
            rt(3) = nv(3)*ndr + ( r(i,3)-nv(3)*ndr ) * cosa + rcn(3) * sina
            ir = irot(i)
            xt(ir) = xt(pivot) + rt(1)
            yt(ir) = yt(pivot) + rt(2)
            zt(ir) = zt(pivot) + rt(3)

        end do ! nRotated

    end subroutine MoveTorsion


    subroutine MoveSpecial(xt, yt, zt)

      use, intrinsic :: iso_c_binding
      use mod_random_tool, only : ranf

      implicit none

      integer, parameter :: max_atom = 5
      integer, parameter :: cd = c_double, ci = c_int
      real(kind=c_double), parameter :: TWO_PI = 6.283185307179586
      real(kind=c_double), dimension(0:max_atom), intent(out) :: xt, yt, zt

      integer(kind=c_int) :: i, i_axis
      real(kind=c_double) :: cosa = 0.0, sina = 0.0, ang = 0.0
      real(kind=c_double) :: norm = 0.0, ndr = 0.0 ! n dot r
      real(kind=c_double), dimension(3) :: nv, rcn, rt, r, pivot_r ! r cross n

      ! aint(3.0 * ranf()): 0 ~ 2.
      ! ang: 0(360), 120, 240.
      ang = TWO_PI / 3.0_cd * aint(3.0_cd * ranf())

      ! set sign of angle (rotational direction).
      if (ranf() > 0.5_cd) ang = -ang

      ! to reduce computaional cost.
      cosa = dcos(ang)
      sina = dsin(ang)

      ! select random pivot atom.
      i_axis = int(4.0_cd * ranf()) + 2 ! 2, 3, 4, 5 (2 is offset).

      ! make spece position vector of pivot axis atom.
      pivot_r = [xt(i_axis), yt(i_axis), zt(i_axis)]

      ! normalize axis vector.
      nv(1) = xt(1) - pivot_r(1)
      nv(2) = yt(1) - pivot_r(2)
      nv(3) = zt(1) - pivot_r(3)

      norm = dsqrt(dot_product(nv, nv))

      nv = nv / norm

      ! rotate each atom.
      do i = 2, 5
        ! neglect pivot atom
        if (i == i_axis) cycle
        ! make relative vector of atom i.
        ! origin of relative coordinates is pivot_r.
        r = [xt(i), yt(i), zt(i)] - pivot_r
        ! make dot and cross product.
        ndr = dot_product(nv, r)
        rcn = cross(r, nv)
        ! rotate relative vector.
        ! rt = nv * ndr + (r - nv * ndr) * cosa + rcn * sina
        ! Goldstein, 3rd edition, page 162.
        rt = r * cosa + nv * ndr * (1.0_cd - cosa) + rcn * sina
        ! move relative  vector to space axis.
        xt(i) = rt(1) + pivot_r(1)
        yt(i) = rt(2) + pivot_r(2)
        zt(i) = rt(3) + pivot_r(3)

      end do

      ! internal subroutines.
      contains

      ! define simple cross product routine.
      function cross(a, b) result(c)

        implicit none

        real(kind=c_double), dimension(3)  :: c
        real(kind=c_double), dimension(3), intent(in) ::  a, b

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)

      end function cross

    end subroutine MoveSpecial

!   --------------------------------------------------------
!   warning : for _rxi, _ryi and _rzi,
!             group coordinate is not updated.

    subroutine ShrinkSolute(ratio_, rxi_, ryi_, rzi_)
!   --------------------------------------------------------
        implicit none


        real(kind=r8), intent(in) :: ratio_
        real(kind=r8), optional, dimension(0:nAtoms_m), intent(inout) :: rxi_, ryi_, rzi_
        !integer, optional, intent(in) :: status

        ! parameters
        real(kind=r8), parameter :: maxRatio = 1.0_r8, minRatio = 0.001_r8

        ! static variables
        real(kind=r8), save :: savedRatio = 1.0_r8 ! initialized

        ! normal variables
        logical :: getOriginalMoleculeQ
        integer :: soli
        real(kind=r8) :: cx, cy, cz ! c = center
        real(kind=r8) :: factor, ratio

!   --------------------------------------------------------

        getOriginalMoleculeQ = .false.
        if (present(rxi_)) getOriginalMoleculeQ = .true.

        ratio = ratio_
        ! limit treatment
        if (ratio < minRatio) ratio = minRatio
        if (ratio > maxRatio) ratio = maxRatio

        ! save center of mass of molecule

        factor = ratio / savedRatio

        soli = soluteIndex

        if (getOriginalMoleculeQ) then

            cx = rxi_(0)
            cy = ryi_(0)
            cz = rzi_(0)

            rxi_ = rxi_ - cx  ! vector operation
            ryi_ = ryi_ - cy
            rzi_ = rzi_ - cz

            rxi_ = rxi_ * factor + cx ! vector operation
            ryi_ = ryi_ * factor + cy
            rzi_ = rzi_ * factor + cz

        else

            cx = rx_m(soli,0)
            cy = ry_m(soli,0)
            cz = rz_m(soli,0)

            rx_m(soli,:) = rx_m(soli,:) - cx ! vector operation
            ry_m(soli,:) = ry_m(soli,:) - cy
            rz_m(soli,:) = rz_m(soli,:) - cz

            rx_m(soli,:) = rx_m(soli,:) * factor + cx ! vector operation
            ry_m(soli,:) = ry_m(soli,:) * factor + cy
            rz_m(soli,:) = rz_m(soli,:) * factor + cz

            call CalculateGroupCenter(soluteIndex) ! update group coordinate

            savedRatio = ratio  ! save current ratio

        end if

        ! group coordinate is not updated.

    end subroutine ShrinkSolute

!   --------------------------------------------------------
    subroutine ShrinkOldNew(lambdaNew, xold_, yold_, zold_, &
                                       xnew_, ynew_, znew_)
!   --------------------------------------------------------
        implicit none

        real(kind=r8), intent(in) :: lambdaNew
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: &
        xold_, yold_, zold_, xnew_, ynew_, znew_

        integer :: soli
!   --------------------------------------------------------

        soli = soluteIndex

        xold_ = rx_m(soli,:)
        yold_ = ry_m(soli,:)
        zold_ = rz_m(soli,:)

        call ShrinkSolute(lambdaNew, xnew_, ynew_, znew_)
        call GroupOldNew(soluteIndex)

    end subroutine ShrinkOldNew

!   -------------------------------------------------------
    subroutine ReCenterize(xt, yt, zt)
!   -------------------------------------------------------
        implicit none

        integer :: a
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xt, yt, zt ! moved atomic position

        real(kind=r8) :: x0, y0, z0, wt
!   -------------------------------------------------------

        x0 = 0.0_r8
        y0 = 0.0_r8
        z0 = 0.0_r8

        do a = 1, nAtoms_m

            wt = mass(myType(a))
            x0 = x0 + xt(a) * wt
            y0 = y0 + yt(a) * wt
            z0 = z0 + zt(a) * wt

        end do

        xt(0) = x0 / molMass
        yt(0) = y0 / molMass
        zt(0) = z0 / molMass

    end subroutine ReCenterize
!   #######################################################

    subroutine SysCOMToOrigin()

        implicit none

        real(kind=r8) :: cX, cY, cZ

        cX = sum(rX_m(:,0)) / real(molNum_m) ! mass / mass = 1
        cY = sum(rY_m(:,0)) / real(molNum_m)
        cZ = sum(rZ_m(:,0)) / real(molNum_m)

        rX_m = rX_m - cX
        rY_m = rY_m - cY
        rZ_m = rZ_m - cZ

        call CalculateGroupCenter(0) ! recalculate group center

    end subroutine SysCOMToOrigin

    !********************










































    !********************

    subroutine TransPeriodicOldNew (i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)

        use mod_random_tool, only : Ranf

        implicit none

        integer, intent(in) :: i
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
        real(kind=r8) :: noX, noY, noZ
        real(kind=r8) :: dX, dY, dZ
        real(kind=r8) :: halfBox, Box

        halfBox = 0.5_r8 * box_m(1,1)
        Box     = box_m(1,1)

        xOld_ = rX_m(i,:)
        yOld_ = rY_m(i,:)
        zOld_ = rZ_m(i,:)

        dX = drMax*(2.0_r8*Ranf()-1.0_r8)
        dY = drMax*(2.0_r8*Ranf()-1.0_r8)
        dZ = drMax*(2.0_r8*Ranf()-1.0_r8)

        noX = box_m(1,1)*dX
        noY = box_m(2,2)*dY
        noZ = box_m(3,3)*dZ

        xNew_ = xOld_ + noX
        yNew_ = yOld_ + noY
        zNew_ = zOld_ + noZ

        if      (xNew_(0) >  halfBox) then; xNew_ = xNew_ - Box
        else if (xNew_(0) < -halfBox) then; xNew_ = xNew_ + Box; end if

        if      (yNew_(0) >  halfBox) then; yNew_ = yNew_ - Box
        else if (yNew_(0) < -halfBox) then; yNew_ = yNew_ + Box; end if

        if      (zNew_(0) >  halfBox) then; zNew_ = zNew_ - Box
        else if (zNew_(0) < -halfBox) then; zNew_ = zNew_ + Box; end if

    end subroutine TransPeriodicOldNew

    !********************

    subroutine ComboOldNew (i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_) ! 미완성

        use mod_random_tool, only : Ranf

        implicit none

        integer, intent(in) :: i
        integer :: xyz ! xyz : x축 : x, y축 : y, z축 : z
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
        real(kind=r8) :: noX, noY, noZ
        real(kind=r8) :: dX, dY, dZ
        real(kind=r8), dimension(0:nAtoms_m) :: fdX, fdY, fdZ
        real(kind=r8) :: cos_, sin_, ang, oX, oY, oZ ! dRotMax

        xOld_ = rX_m(i,:)
        yOld_ = rY_m(i,:)
        zOld_ = rZ_m(i,:)

        ! rotate

        oX = xOld_(0)
        oY = yOld_(0)
        oZ = zOld_(0)

        fdX = xOld_ - oX
        fdY = yOld_ - oY
        fdZ = zOld_ - oZ

        xyz  = int(3.0_r8*Ranf() + 1.0_r8)
        ang  = dRotMax*(2.0_r8*Ranf()-1.0_r8)
        cos_ = cos(ang)
        sin_ = sin(ang)

        if (xyz == 1) then

            xNew_ = fdX                      + oX
            yNew_ =     cos_*fdY - sin_*fdZ  + oY
            zNew_ =     sin_*fdY + cos_*fdZ  + oZ

        else if (xyz == 2) then

            xNew_ = cos_*fdX      -sin_*fdZ  + oX
            yNew_ =           fdY            + oY
            zNew_ = sin_*fdX      +cos_*fdZ  + oZ

        else if (xyz == 3) then

            xNew_ = cos_*fdX - sin_*fdY      + oX
            yNew_ = sin_*fdX + cos_*fdY      + oY
            zNew_ =                     fdZ  + oZ

        else; stop 'critical error occur in rotation'; end if

        ! translate

        dX = drMax*(2.0_r8*Ranf()-1.0_r8)
        dY = drMax*(2.0_r8*Ranf()-1.0_r8)
        dZ = drMax*(2.0_r8*Ranf()-1.0_r8)

        noX = box_m(1,1)*dX + box_m(1,2)*dY + box_m(1,3)*dZ
        noY = box_m(2,1)*dX + box_m(2,2)*dY + box_m(2,3)*dZ
        noZ = box_m(3,1)*dX + box_m(3,2)*dY + box_m(3,3)*dZ

        xNew_ = xNew_ + noX
        yNew_ = yNew_ + noY
        zNew_ = zNew_ + noZ

    end subroutine ComboOldNew ! 미완성

    !********************

    subroutine ComboPeriodicOldNew (i, xOld_, yOld_, zOld_, xNew_, yNew_, zNew_) ! 미완성

        use mod_random_tool, only : Ranf

        implicit none

        integer, intent(in) :: i
        integer :: xyz ! xyz : x축 : x, y축 : y, z축 : z
        real(kind=r8), dimension(0:nAtoms_m), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
        real(kind=r8) :: noX, noY, noZ
        real(kind=r8) :: dX, dY, dZ
        real(kind=r8), dimension(0:nAtoms_m) :: fdX, fdY, fdZ
        real(kind=r8) :: cos_, sin_, ang, oX, oY, oZ ! dRotMax
        real(kind=r8) :: halfBox, Box

        halfBox = 0.5_r8 * box_m(1,1)
        Box     = box_m(1,1)

        xOld_ = rX_m(i,:)
        yOld_ = rY_m(i,:)
        zOld_ = rZ_m(i,:)

        ! rotate

        oX = xOld_(0)
        oY = yOld_(0)
        oZ = zOld_(0)

        fdX = xOld_ - oX
        fdY = yOld_ - oY
        fdZ = zOld_ - oZ

        xyz  = int(3.0_r8*Ranf() + 1.0_r8)
        ang  = dRotMax*(2.0_r8*Ranf()-1.0_r8)
        cos_ = cos(ang)
        sin_ = sin(ang)

        if (xyz == 1) then

            xNew_ = fdX                      + oX
            yNew_ =     cos_*fdY - sin_*fdZ  + oY
            zNew_ =     sin_*fdY + cos_*fdZ  + oZ

        else if (xyz == 2) then

            xNew_ = cos_*fdX      -sin_*fdZ  + oX
            yNew_ =           fdY            + oY
            zNew_ = sin_*fdX      +cos_*fdZ  + oZ

        else if (xyz == 3) then

            xNew_ = cos_*fdX - sin_*fdY      + oX
            yNew_ = sin_*fdX + cos_*fdY      + oY
            zNew_ =                     fdZ  + oZ

        else; stop 'critical error occur in rotation'; end if

        ! translate

        dX = drMax*(2.0_r8*Ranf()-1.0_r8)
        dY = drMax*(2.0_r8*Ranf()-1.0_r8)
        dZ = drMax*(2.0_r8*Ranf()-1.0_r8)

        noX = box_m(1,1)*dX + box_m(1,2)*dY + box_m(1,3)*dZ
        noY = box_m(2,1)*dX + box_m(2,2)*dY + box_m(2,3)*dZ
        noZ = box_m(3,1)*dX + box_m(3,2)*dY + box_m(3,3)*dZ

        xNew_ = xNew_ + noX
        yNew_ = yNew_ + noY
        zNew_ = zNew_ + noZ

        if      (xNew_(0) >  halfBox) then; xNew_ = xNew_ - Box
        else if (xNew_(0) < -halfBox) then; xNew_ = xNew_ + Box; end if

        if      (yNew_(0) >  halfBox) then; yNew_ = yNew_ - Box
        else if (yNew_(0) < -halfBox) then; yNew_ = yNew_ + Box; end if

        if      (zNew_(0) >  halfBox) then; zNew_ = zNew_ - Box
        else if (zNew_(0) < -halfBox) then; zNew_ = zNew_ + Box; end if

    end subroutine ComboPeriodicOldNew ! 미완성

    !********************

    !********************

    subroutine ChangeCubicBox (xNew, yNew, zNew, nBox, niBox, nVol)

        use mod_random_tool, only : Ranf
        use mod_math_tool,   only : Inverse, Det

        implicit none

        integer :: i, j
        real(kind=r8), dimension(molNum_m,0:nAtoms_m), intent(out) :: xNew, yNew, zNew
        real(kind=r8), dimension(3,3), intent(out) :: nBox, niBox
        real(kind=r8), intent(out) :: nVol

        logical :: singularQ

        real(kind=r8), dimension(molNum_m) :: dx, dy, dz, dx0, dy0, dz0

        nBox(1,1) = box_m(1,1) + dBoxMax*(2.0_r8*Ranf()-1.0_r8)
        do i = 1, 3; do j = 1, 3;
            if (i==j) then
                nBox(i,j)  = nBox(1,1)
                niBox(i,j) = 1.0_r8 / nBox(1,1)
            else
                nBox(i,j)  = 0.0_r8
                niBox(i,j) = 0.0_r8
             end if
        end do; end do
        nVol = nBox(1,1)**3

        dx0 = rX_m(:,0)
        dy0 = rY_m(:,0)
        dz0 = rZ_m(:,0)

        dx = iBox_m(1,1)*dx0 + iBox_m(1,2)*dy0 + iBox_m(1,3)*dz0
        dy =                   iBox_m(2,2)*dy0 + iBox_m(2,3)*dz0
        dz =                                     iBox_m(3,3)*dz0

        ! diagnal만 된다.

        dx = nBox(1,1)*dx
        dy =                nBox(2,2)*dy
        dz =                               nBox(3,3)*dz

        do i = 0, nAtoms_m

            xNew(:,i) = rX_m(:,i) - dx0 + dx
            yNew(:,i) = rY_m(:,i) - dy0 + dy
            zNew(:,i) = rZ_m(:,i) - dz0 + dz

        end do

    end subroutine ChangeCubicBox

    !********************

end module mod_movement_information
