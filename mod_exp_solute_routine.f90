!----------------------------------------------------
!subroutine ReadInput()
!subroutine InitializeSimulation()
!subroutine SequenceMovement(i)
!subroutine ComboMovement(i)
!subroutine MoveMolecule()
!subroutine MoveSimulationBox()
!subroutine PrintSimulationResult()
!subroutine ResetForMainloop()
!subroutine FinalizeSimulation()
!----------------------------------------------------

!=========================================
 module mod_exp_solute_routine
!=========================================
    use mod_program_operation

    use mod_univ_const

    use mod_time_tool
    use mod_random_tool,    only : RandomSettingAll
    use mod_character_tool, only : IntToCha

    use mod_system_information
    use mod_molecule_information
    use mod_config_information
    use mod_movement_information
    use mod_energy_core
    use mod_neighbor_list
    use mod_data_manipulation

    use mod_expans_solid
!===================================================
    implicit none

    integer, parameter, private :: r8 = selected_real_kind(15,300)

    public

    integer, parameter :: devide = 100, INFINITY = 1000000000
    integer :: maxCycle, freeCycle, mainCycle, binCycle
    logical :: freeQ, acceptQ ! free step for bin estimating step
    real(kind=r8) :: maxbeu, minbeu, beuOld, beuNew, acBeu, acBeuSq, nBeu, testBeu

    integer :: i, eeCycle
    character(len=31) :: filename
    real(kind=r8) :: uInter, uIntra, uEin, dbe, eOld, eNew, &
                     eOldIntra, eNewIntra, eTot, testEnergy, &
                     lambda, de, uRef, uTar       ! reference energy (ein crystal), Target energy (real crystal)

    real(kind=r8) :: dummyEnergy = 0.0_r8
!===================================================
    contains

!   ------------------------------------------------
    subroutine ReadInput()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------
        filename = 'expanded.txt'
        call GetExpSolidInf(filename)
        call AutoOperation
        filename = 'data'//IntToCha(partition)//'.txt'
        open(1,file=filename) ! open datafile to write

        filename = 'system.txt'
        call GetSysInf (filename)
        call PrintSysInf (1)
        call PrintSysInf

        filename = 'molecule.txt'
        call GetMolInf (filename)
        call PrintMolInf (1)
        call PrintMolInf

        filename = 'config.txt'
        !filename = 'nvt_configReault.txt'
        call GetCoordinate (filename)
        call PrintConfigInf (1)
        call PrintConfigInf

        filename = 'movement.txt'
        call GetMoveInf (filename)
        call PrintMoveInf (1)
        call PrintMoveInf

        call PrintExpSolidInf (1)
        call PrintExpSolidInf

        call OKsign

    end subroutine ReadInput






!   ------------------------------------------------
    subroutine InitializeSimulation()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        call InitializeNeighborList()
        call RandomSettingAll()

        call GetLambda(lambda)
        call ShrinkSolute(lambda)    ! initialize solute scale
        call CalculateGroupCenter(0) ! calculate center of groups

        call CalculateEnergy(0, 'OLD', eVdw, eCha, lambda)
        if (flexQ) call SumupIntra(eIntra)  ! independent for lambda

        sequenceQ = .false.   ! only combo move for EEM calculation.

        uInter = eVdw + eCha + EnergyLongRange(hBox) ! Van der Waals + charge + Long Range Correction.
        uIntra = eIntra                              ! Bond Vib + Bending Stratching + Torsion + Non-Bond interaction.
        eTot   = uInter + uIntra                     ! total energy of system.

        write(*,"(3(A,F10.3))") ' inter = ', uInter / dble(molNum) * KtokJ, &
                                ' intra = ', uIntra / dble(molNum) * KtokJ, &
                                ' LRC   = ', EnergyLongRange(hBox) / dble(molNum - 1) * kTokJ ! neglect solute
        write(1,"(3(A,F10.3))") ' inter = ', uInter / dble(molNum) * KtokJ, &
                                ' intra = ', uIntra / dble(molNum) * KtokJ, &
                                ' LRC   = ', EnergyLongRange(hBox) / dble(molNum - 1) * kTokJ

    end subroutine InitializeSimulation










!   ------------------------------------------------
    subroutine MoveMolecule()
!   ------------------------------------------------
        implicit none

        integer :: i
!   ------------------------------------------------
        do i = 1, molNum

            call CheckList()

            if (sequenceQ) then ! sequence movement, never used

                !call SequenceMovement(i)

            else

                call ComboMovement(i)

            end if

        end do ! i: molecule number

    end subroutine MoveMolecule






!   ------------------------------------------------
    subroutine ComboMovement(i)
!   ------------------------------------------------
        implicit none

        integer, intent(in) :: i
!   ------------------------------------------------

        call InterOldNew (i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)
        call GetLambda(lambda) ! get current lambda.
        call CalculateEnergy(i, 'OLD', eVdw, eCha, lambda) ! old energy
        eOld = eVdw + eCha
        call CalculateEnergy(i, 'NEW', eVdw, eCha, lambda) ! new energy
        eNew = eVdw + eCha

        dE   = eNew - eOld
        dbe  = beta * dE
        call CheckAccept (dbe, acceptQ)

        if (acceptQ) then

            call UpdateCoordinate(i)

            eTot   = eTot   - eOld + eNew
            uInter = uInter - eOld + eNew

            !call PRT(i)
            !call CalEne(i)

        end if

        if (flexQ) then
           !intramolecular motion >>
!           ----------------------------------------------------------------------------------
            if (bondQ .or. bendQ .or. torsQ) then ! only for molecules having at least one intra movement
!           ----------------------------------------------------------------------------------
                call IntraOldNew(i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)

                call GetLambda(lambda) ! get current lambda.

                call CalculateEnergy(i, 'OLD', eVdw, eCha, lambda)
                call EnergyIntra (rxiOld, ryiOld, rziOld, eOldIntra, i)
                eOld = eVdw + eCha

                call CalculateEnergy(i, 'NEW', eVdw, eCha, lambda)
                call EnergyIntra (rxiNew, ryiNew, rziNew, eNewIntra, i)
                eNew = eVdw + eCha

                dE   = eNew - eOld
                dbe  = beta * dE
                call CheckAccept (dbe, acceptQ)

                if (acceptQ) then

                    call UpdateCoordinate(i)

                    uInter = uInter - eOld      + eNew
                    uIntra = uIntra - eOldIntra + eNewIntra
                    eTot   = uInter + uIntra

                end if

            end if ! bendQ or torsQ
!           -----------------------------------------------------------------------
        end if ! end flex move

    end subroutine ComboMovement





!   ------------------------------------------------
    subroutine PrintSimulationResult(cycle_)
!   ------------------------------------------------
        implicit none

        logical, save :: iniQ = .true.
        integer, intent(in) :: cycle_
!   ------------------------------------------------

        if (iniQ) then

            write(*,"(4(A10))") 'cycle', 'u inter', 'u intra'
            write(1,"(4(A10))") 'cycle', 'u inter', 'u intra'

            iniQ = .false.

        end if

        filename = 'view'//IntToCha(partition)//'.view'
        if (mod(cycle_, visualStep) == 0) call MakeConfigView (filename)
        if (mod(cycle_, printStep)  == 0) then

            write(*,"(I10, 3(F10.3))") cycle_, &
                                       uInter / real(molNum) * KtokJ, &
                                       uIntra / real(molNum) * KtokJ


            write(1,"(I10, 3(F10.3))") cycle_, &
                                       uInter / real(molNum) * KtokJ, &
                                       uIntra / real(molNum) * KtokJ

        end if


    end subroutine PrintSimulationResult




!   ------------------------------------------------
    subroutine FinalizeSimulation()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        write(*,"(A20)") 'check energy'
        write(1,"(A20)") 'check energy'

        write(*,"(3(A20))") 'u inter', 'u intra', 'eTot'
        write(1,"(3(A20))") 'u inter', 'u intra', 'eTot'

        write(*,"(3(F20.3))") uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              eTot   / real(molNum) * KtokJ

        write(1,"(3(F20.3))") uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              eTot   / real(molNum) * KtokJ

        call GetLambda(lambda) ! get current lambda.
        call CalculateEnergy(0, 'OLD', eVdw, eCha, lambda)
        if (flexQ) call SumupIntra(eIntra)

        uInter = eVdw + eCha + EnergyLongRange(hBox) ! Van der Waals + charge + Long Range Correction.
        uIntra = eIntra                              ! Bond Vib + Bending Stratching + Torsion + Non-Bond interaction.
        eTot   = uInter + uIntra                     ! total energy of system... is it used...?

        write(*,"(3(F20.3))") uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              eTot   / real(molNum) * KtokJ

        write(1,"(3(F20.3))") uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              eTot   / real(molNum) * KtokJ

        call PrintFreeEnergy()

    end subroutine FinalizeSimulation


    subroutine CalEne(i)

        implicit none

        integer, intent(in) :: i

        real(kind=r8) :: uInter_, uIntra_, eTot_

        call GetLambda(lambda) ! get current lambda.
        call CalculateEnergy(0, 'OLD', eVdw, eCha, lambda)
        if (flexQ) call SumupIntra(eIntra)
        uInter_ = eVdw + eCha + EnergyLongRange(hBox) ! Van der Waals + charge + Long Range Correction.
        uIntra_ = eIntra                              ! Bond Vib + Bending Stratching + Torsion + Non-Bond interaction.
        eTot_   = uInter_ + uIntra_                   ! total energy of system... is it used...?

        write(*,"(I5, 3(F20.3))") i, uInter_ / real(molNum) * KtokJ, &
                                     uIntra_ / real(molNum) * KtokJ, &
                                     eTot_   / real(molNum) * KtokJ

        write(1,"(I5, 3(F20.3))") i, uInter_ / real(molNum) * KtokJ, &
                                     uIntra_ / real(molNum) * KtokJ, &
                                     eTot_   / real(molNum) * KtokJ

        call system('pause')

    end subroutine CalEne

    subroutine PRT(i)

        implicit none

        integer, intent(in) :: i

        write(*,"(I5, 3(F20.3))") i, uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              eTot   / real(molNum) * KtokJ

        write(1,"(I5, 3(F20.3))") i, uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              eTot   / real(molNum) * KtokJ

    end subroutine PRT

end module mod_exp_solute_routine



