!----------------------------------------------------
!subroutine ReadInput()
!subroutine InitializeSimulation()
!subroutine MoveMolecule()
!subroutine SequenceMovement(i)
!subroutine ComboMovement(i)
!subroutine MoveSimulationBox()
!subroutine PrintSimulationResult()
!subroutine ResetForMainloop()
!subroutine FinalizeSimulation()
!----------------------------------------------------

!===================================================
 module mod_simulation_routine
!===================================================
    use mod_program_operation

    use mod_univ_const

    use mod_time_tool
    use mod_random_tool, only : RandomSettingAll, ranf

    use mod_system_information
    use mod_molecule_information
    use mod_config_information
    use mod_movement_information
    use mod_energy_core
    use mod_neighbor_list
    use mod_data_manipulation

    use mod_einstein_core

!===================================================
    implicit none

    integer, parameter, private :: r8 = selected_real_kind(15,300)

    public

    integer :: i, step
    character(len=31) :: filename
    logical :: acceptQ, adjustQ
    real(kind=r8) :: acE, nacE, uInter, uIntra, dbe, eOld, eNew, &
                     eOldIntra, eNewIntra, eTot, acIntra, nacIntra

    contains

!   ------------------------------------------------
    subroutine ReadInput()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        filename = 'system.txt'
        call GetSysInf (filename)
        open(1,file=trim(simulType)//'_data.txt') ! open datafile to write
        call PrintSysInf (1)
        call PrintSysInf

        filename = 'molecule.txt'
        call GetMolInf (filename)
        call PrintMolInf (1)
        call PrintMolInf

        filename = 'config.txt'
        call GetCoordinate (filename)
        call PrintConfigInf (1)
        call PrintConfigInf

        if (.not.configReadQ) then ! if no configuration, make initial configuration

            call MakeLatticePoint ! make crystal lattice for center of mass of molecule
            call GetMoleculeBody  ! get molecule frame to add lattice frame
            call GetLatticePoint  ! add molecule coordinate to lattice

        end if

        filename = trim('initial_config.txt')
        call WriteCoordinate (filename, hbox)
        call CalculateGroupCenter(0)
        !call PrintGroupCoordinate

        filename = 'movement.txt'
        call GetMoveInf (filename)
        call PrintMoveInf (1)
        call PrintMoveInf

        !filename = 'EinRotInf.txt'
        !call GetEinRotInf(filename)

        call InitCheckProb()

    end subroutine ReadInput


!   ------------------------------------------------
    subroutine InitializeSimulation()
!   ------------------------------------------------
        implicit none

        integer :: i = 1 ! why?
!   ------------------------------------------------

        ! neighbor list algorithm
        call InitializeNeighborList()
        !call PrintNeighborList()
        !call PrintNeighborList()

        call RandomSettingAll

        call ResetData()
        !write(*,*) 'energy error?'
        !call Sumup(eVdw, eCha, hBox, ihBox)
        call CalculateEnergy(0, 'OLD', eVdw, eCha)
        !write(*,*) 'no energy error'
        !write(*,*) 'intra error?'
        if (flexQ) call SumupIntra(eIntra)
        !write(*,*) 'no intra error'
        uInter = eVdw + eCha + EnergyLongRange(hBox)
        uIntra = eIntra
        eTot   = uInter + uIntra
        write(*,"(3(A,F10.3))") ' inter = ', uInter / dble(molNum) * KtokJ, &
                                ' intra = ', uIntra / dble(molNum) * KtokJ, &
                                ' LRC   = ', EnergyLongRange(hBox) / dble(molNum) * KtokJ
        write(1,"(3(A,F10.3))") ' inter = ', uInter / dble(molNum) * KtokJ, &
                                ' intra = ', uIntra / dble(molNum) * KtokJ, &
                                ' LRC   = ', EnergyLongRange(hBox) / dble(molNum) * KtokJ
        adjustQ = .true.
        acE     = uInter
        nacE    = 1.0_r8
        acBox   = hBox
        nacBox  = 1.0_r8
        call AccVolume()
        call SysCOMToOrigin ! change system center to origin

        call AccumulateData(0, uInter, uIntra, hBox)

    end subroutine InitializeSimulation

!   ------------------------------------------------
    subroutine MoveMolecule()
!   ------------------------------------------------
        implicit none

        integer :: i
!   ------------------------------------------------

        call ChangeToCombo(step)
        call CheckList()

        do i = 1, molNum

            if (sequenceQ) then ! sequence movement

                call SequenceMovement(i)

            else

                call ComboMovement(i)

            end if

        end do ! i: molecule number

        call MakeAngleHist()

    end subroutine MoveMolecule

!   ------------------------------------------------
    subroutine SequenceMovement(i)
!   ------------------------------------------------
        implicit none

        integer, intent(in) :: i
!   ------------------------------------------------

        call TransOldNew (i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)
        !call Energy (i, rxiOld, ryiOld, rziOld, eVdw, eCha, hBox, ihBox)
        call CalculateEnergy(i, 'OLD', eVdw, eCha)
        eOld = eVdw + eCha
        !call Energy (i, rxiNew, ryiNew, rziNew, eVdw, eCha, hBox, ihBox)
        call CalculateEnergy(i, 'NEW', eVdw, eCha)
        eNew = eVdw + eCha

        dbe  = (eNew - eOld) / temper
        call CheckAccept (dbe, acceptQ)

        if (acceptQ) then

            call UpdateCoordinate(i)
            !call CheckList()
            !rx(i,:) = rxinew
            !ry(i,:) = ryinew
            !rz(i,:) = rzinew

            eTot   = eTot - eOld + eNew
            uInter = uInter - eOld + eNew

            nTrans = nTrans + 1.0_r8

        end if

        tTrans = tTrans + 1.0_r8

        acE = acE + uInter; nacE = nacE + 1.0_r8
        call AccumulateData(step, uInter, uIntra, hBox)

        !orientational motion >>
!       -------------------------------------------------------------------------
        if (rotateQ) then  ! for nAtoms > 2 molecules
!       -------------------------------------------------------------------------
            call RotateOldNew (i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)
            !call Energy (i, rxiOld, ryiOld, rziOld, eVdw, eCha, hBox, ihBox)
            call CalculateEnergy(i, 'OLD', eVdw, eCha)
            eOld = eVdw + eCha
            !call Energy (i, rxiNew, ryiNew, rziNew, eVdw, eCha, hBox, ihBox)
            call CalculateEnergy(i, 'NEW', eVdw, eCha)
            eNew = eVdw + eCha

            dbe  = (eNew-eOld)/temper

            call CheckAccept (dbe, acceptQ)

            if (acceptQ) then

                call UpdateCoordinate(i)
                !call CheckList()
                !rx(i,:) = rxinew
                !ry(i,:) = ryinew
                !rz(i,:) = rzinew

                eTot   = eTot   - eOld + eNew
                uInter = uInter - eOld + eNew

                nRotate = nRotate + 1.0_r8

            end if

            tRotate = tRotate + 1.0_r8

            acE = acE + uInter; nacE = nacE + 1.0_r8
            call AccumulateData(step, uInter, uIntra, hBox)

        end if ! rotateQ
!       -------------------------------------------------------------------------


!       -------------------------------------------------------------------------
        if (flexQ) then ! only flexible molecules
!       -------------------------------------------------------------------------
            if (bondQ) then ! only for nBend > 0 molecules
!           ---------------------------------------------------------------------
                !bending motion >>
                call BondOldNew(i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)

                !call Energy (i, rxiOld, ryiOld, rziOld, eVdw, eCha, hBox, ihBox)
                call CalculateEnergy(i, 'OLD', eVdw, eCha)
                call EnergyIntra (rxiOld, ryiOld, rziOld, eOldIntra)
                eOld = eVdw + eCha
                !call Energy (i, rxiNew, ryiNew, rziNew, eVdw, eCha, hBox, ihBox)
                call CalculateEnergy(i, 'NEW', eVdw, eCha)
                call EnergyIntra (rxiNew, ryiNew, rziNew, eNewIntra)
                eNew = eVdw + eCha

                dbe  = (eNew + eNewIntra - eOld - eOldIntra) / temper

                call CheckAccept (dbe, acceptQ)

                if (acceptQ) then

                    call UpdateCoordinate(i)
                    !call CheckList()
                    !rx(i,:) = rxinew
                    !ry(i,:) = ryinew
                    !rz(i,:) = rzinew

                    uInter = uInter - eOld + eNew
                    uIntra = uIntra - eOldIntra + eNewIntra
                    eTot   = uInter + uIntra

                    nBond_ = nBond_ + 1.0_r8

                end if

                tBond = tBond + 1.0_r8

                acE = acE + uInter; nacE = nacE + 1.0_r8
                call AccumulateData(step, uInter, uIntra, hBox)

            end if ! bondQ

            if (bendQ) then ! only for nBend > 0 molecules
!           ---------------------------------------------------------------------
                !bending motion >>
                call BendOldNew(i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)

                !call Energy (i, rxiOld, ryiOld, rziOld, eVdw, eCha, hBox, ihBox)
                call CalculateEnergy(i, 'OLD', eVdw, eCha)
                call EnergyIntra (rxiOld, ryiOld, rziOld, eOldIntra)
                eOld = eVdw + eCha
                !call Energy (i, rxiNew, ryiNew, rziNew, eVdw, eCha, hBox, ihBox)
                call CalculateEnergy(i, 'NEW', eVdw, eCha)
                call EnergyIntra (rxiNew, ryiNew, rziNew, eNewIntra)
                eNew = eVdw + eCha

                dbe  = (eNew + eNewIntra - eOld - eOldIntra) / temper

                call CheckAccept (dbe, acceptQ)

                if (acceptQ) then

                    call UpdateCoordinate(i)
                    !call CheckList()
                    !rx(i,:) = rxinew
                    !ry(i,:) = ryinew
                    !rz(i,:) = rzinew

                    uInter = uInter - eOld + eNew
                    uIntra = uIntra - eOldIntra + eNewIntra
                    eTot   = uInter + uIntra

                    nBend_ = nBend_ + 1.0_r8

                end if

                tBend = tBend + 1.0_r8

                acE = acE + uInter; nacE = nacE + 1.0_r8
                call AccumulateData(step, uInter, uIntra, hBox)

            end if ! bendQ
!           ---------------------------------------------------------------------
            if (torsQ) then ! only for nTorsion > 0 molecules
!           ---------------------------------------------------------------------
                !torsion motion >>
                call TorsOldNew(i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)

                !call Energy (i, rxiOld, ryiOld, rziOld, eVdw, eCha, hBox, ihBox)
                call CalculateEnergy(i, 'OLD', eVdw, eCha)
                call EnergyIntra (rxiOld, ryiOld, rziOld, eOldIntra)
                eOld = eVdw + eCha
                !call Energy (i, rxiNew, ryiNew, rziNew, eVdw, eCha, hBox, ihBox)
                call CalculateEnergy(i, 'NEW', eVdw, eCha)
                call EnergyIntra (rxiNew, ryiNew, rziNew, eNewIntra)
                eNew = eVdw + eCha

                dbe  = (eNew + eNewIntra - eOld - eOldIntra) / temper

                call CheckAccept (dbe, acceptQ)

                if (acceptQ) then

                    call UpdateCoordinate(i)
                    !call CheckList()
                    !rx(i,:) = rxinew
                    !ry(i,:) = ryinew
                    !rz(i,:) = rzinew

                    uInter = uInter - eOld + eNew
                    uIntra = uIntra - eOldIntra + eNewIntra
                    eTot   = uInter + uIntra

                    nTors = nTors + 1.0_r8

                end if

                tTors = tTors + 1.0_r8

                acE = acE + uInter; nacE = nacE + 1.0_r8
                call AccumulateData(step, uInter, uIntra, hBox)

            end if ! torsQ
!           ---------------------------------------------------------------------
        end if ! end flex move
!       -------------------------------------------------------------------------
    end subroutine sequenceMovement


!   ------------------------------------------------
    subroutine ComboMovement(i)
!   ------------------------------------------------
        implicit none

        integer, intent(in) :: i
!   ------------------------------------------------

        call InterOldNew (i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)
        !call Energy (i, rxiOld, ryiOld, rziOld, eVdw, eCha, hBox, ihBox)
        call CalculateEnergy(i, 'OLD', eVdw, eCha)
        eOld = eVdw + eCha
        !call Energy (i, rxiNew, ryiNew, rziNew, eVdw, eCha, hBox, ihBox)
        call CalculateEnergy(i, 'NEW', eVdw, eCha)
        eNew = eVdw + eCha

        dbe  = (eNew - eOld) / temper
        call CheckAccept (dbe, acceptQ)

        if (acceptQ) then

            call UpdateCoordinate(i)

            eTot   = eTot - eOld + eNew
            uInter = uInter - eOld + eNew

            nInter = nInter + 1.0_r8

        end if

        tInter = tInter + 1.0_r8

        acE = acE + uInter; nacE = nacE + 1.0_r8
        call AccumulateData(step, uInter, uIntra, hBox)
        call AccEin(i)

        if (flexQ) then
           !intramolecular motion >>
!           ----------------------------------------------------------------------------------
            if (bondQ .or. bendQ .or. torsQ) then ! only for molecules having at least one intra movement
!           ----------------------------------------------------------------------------------
                !call ShrinkSolute(ranf())
                call IntraOldNew(i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)

                !call Energy (i, rxiOld, ryiOld, rziOld, eVdw, eCha, hBox, ihBox)
                call CalculateEnergy(i, 'OLD', eVdw, eCha)
                call EnergyIntra (rxiOld, ryiOld, rziOld, eOldIntra)
                eOld = eVdw + eCha
                !call Energy (i, rxiNew, ryiNew, rziNew, eVdw, eCha, hBox, ihBox)
                call CalculateEnergy(i, 'NEW', eVdw, eCha)
                call EnergyIntra (rxiNew, ryiNew, rziNew, eNewIntra)
                eNew = eVdw + eCha

                dbe  = (eNew + eNewIntra - eOld - eOldIntra) / temper

                call CheckAccept (dbe, acceptQ)

                if (acceptQ) then

                    call UpdateCoordinate(i)
                    !call CheckList()

                    !rx(i,:) = rxinew
                    !ry(i,:) = ryinew
                    !rz(i,:) = rzinew

                    uInter = uInter - eOld + eNew
                    uIntra = uIntra - eOldIntra + eNewIntra
                    eTot   = uInter + uIntra

                    nIntra = nIntra + 1.0_r8

                end if

                tIntra = tIntra + 1.0_r8

                acE = acE + uInter; nacE = nacE + 1.0_r8
                call AccumulateData(step, uInter, uIntra, hBox)
                call AccEin(i)

            end if ! bendQ or torsQ
!           -----------------------------------------------------------------------
        end if ! end flex move

    end subroutine ComboMovement




!   ------------------------------------------------
    subroutine MoveSimulationBox()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        if (.not.boxChangeQ) return

        if(mod(step,2) == 1) return ! change box each two step

        call ChangeBox (rxNew, ryNew, rzNew, newhBox, newihBox, newVolume)
        !call ChangeCubicBox (rxNew, ryNew, rzNew, newhBox, newihBox, newVolume)
        !call SumupCurrent (rxNew, ryNew, rzNew, eVdw, eCha, newhBox, newihBox)
        call CalculateEnergy(0, 'NEW', eVdw, eCha)
        eOld = uInter
        eNew = eVdw + eCha + EnergyLongRange (newhBox) ! new intra energy

        dbe  = beta*(eNew - eOld + BAR_A3toK*pressure*(newVolume-volume)) - &
               real(molNum)*log(newVolume/volume)
        call CheckAccept (dbe, acceptQ)

        if (acceptQ) then

            call UpdateCoordinate(0)
            !call CheckList()

            !rx = rxnew
            !ry = rynew
            !rz = rznew

            hBox  = newhBox
            ihBox = newihBox

            volume = newVolume

            eTot   = eTot   - eOld + eNew
            uInter = uInter - eOld + eNew

            nBox = nBox + 1.0_r8

        end if

        tBox = tBox + 1.0_r8

        acE   = acE + uInter; nacE   = nacE   + 1.0_r8
        acBox = acBox + hBox; nacBox = nacBox + 1.0_r8
        call AccumulateData(step, uInter, uIntra, hBox)
        !call AccEin(0) ! sampling for all molecules -> not used...

        call AccVolume

    end subroutine MoveSimulationBox

!   ------------------------------------------------
    subroutine PrintSimulationResult()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        filename = 'view.view'
        if (mod(step,visualStep)==0) call MakeConfigView (filename)

        if (mod(step,printStep) ==0) then

            call PrintTabulatedData()

            write(*,"('step = ',I10, '   energy = ', F20.5,' kJ/mol')")  step, &
            & acE/nacE/real(molNum)*KtokJ
            call PrintDensity
            call PrintAcceptanceRatio
            write(*,"('current lattice'/)")
            call PrintLatticeInf(hBox)
            write(*,"('average lattice'/)")
            if (boxChangeQ) call PrintLatticeInf(acBox/nacBox)

            write(1,"('step = ',I10, '   energy = ', F20.5,' kJ/mol')")  step, &
            & acE/nacE/real(molNum)*KtokJ
            call PrintDensity (1)
            call PrintAcceptanceRatio (1)
            write(1,"()")
            write(1,"('current lattice'/)")
            call PrintLatticeInf(hBox,1)
            write(1,"('average lattice'/)")
            if (boxChangeQ) call PrintLatticeInf(acBox/nacBox,1)

            filename = trim(currentConfig)
            call WriteCoordinate (filename, hbox)

            avBox = acBox / nacBox
            if (.not.boxChangeQ) avBox = hBox
            filename = trim(configResult)
            call WriteCoordinate (filename, avBox)

            call PrintRemainTime
            call PrintRemainTime (1)

            call WriteEinInf()
            call WriteAngleProb()
            !call WriteOrientationFactor(step)

        end if

    end subroutine PrintSimulationResult


!   ------------------------------------------------
    subroutine ResetForMainloop()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        if (step == preStep) then

            write(*,"('main step start'/)")
            write(1,"('main step start'/)")

            acE   = 0.0_r8; nacE   = 0.0_r8
            acBox = 0.0_r8; nacBox = 0.0_r8
            adjustQ = .false.

            ! re-calculate system energy.

            call SetupList()
            call CalculateEnergy(0, 'OLD', eVdw, eCha)
            if (flexQ) call SumupIntra(eIntra)
            write(*,"('check energy = ', 2F20.6)") eTot, eVdw + eCha + EnergyLongRange(hBox) + eIntra
            write(1,"('check energy = ', 2F20.6)") eTot, eVdw + eCha + EnergyLongRange(hBox) + eIntra

            eTot = eVdw + eCha + EnergyLongRange(hBox) + eIntra
            write(*,"('energy is changed to ', F20.6)") eTot
            write(1,"('energy is changed to ', F20.6)") eTot

            call SysCOMToOrigin() ! change system center to origin
            call AccVolume (.true.) ! reset = true
            call ResetData()
            call AccEin(1,.true.)
            call ResetAngleHist()

        end if

    end subroutine ResetForMainloop


!   ------------------------------------------------
    subroutine FinalizeSimulation()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        !call PrintNeighborList()
        !filename = trim('r_data.txt')
        !call WriteCoordinate (filename, hbox)
        !call PrintGroupCoordinate

        call WriteMoveInf

        !call Sumup(eVdw, eCha, hBox, ihBox)
        call SetupList()
        call CalculateEnergy(0, 'OLD', eVdw, eCha)
        if (flexQ) call SumupIntra(eIntra)
        write(*,"('check energy = ', 2F20.6)") eTot, eVdw + eCha + EnergyLongRange(hBox) + eIntra
        write(1,"('check energy = ', 2F20.6)") eTot, eVdw + eCha + EnergyLongRange(hBox) + eIntra

        if (boxChangeQ) hBox = acBox / nacBox ! for save

        filename = trim(configResult)
        call WriteCoordinate (filename, hBox)

        call WriteEinInf()

    end subroutine FinalizeSimulation

end module mod_simulation_routine
