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
 module mod_exp_simulation_routine
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

    use mod_einstein_core
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

!===================================================
    contains

!   ------------------------------------------------
    subroutine ReadInput()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------
        filename = 'expanded.txt'
        call GetExpSolidInf(filename)
        !call AutoOperation
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

        call GetEinInf()

        call PrintEinInf (1)
        call PrintEinInf

        call PrintExpSolidInf (1)
        call PrintExpSolidInf

        !call OKsign

        if (.not. isLastPartition()) call RunNextPartitionMain()

    end subroutine ReadInput






!   ------------------------------------------------
    subroutine InitializeSimulation()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        call InitializeNeighborList()
        call RandomSettingAll()
        call MolCenterToEinLattice()  ! avoid error in fixed one molecule method.
        call CalculateGroupCenter(0)   ! calculate center of groups
        call CalculateEnergy(0, 'OLD', eVdw, eCha)
        if (flexQ) call SumupIntra(eIntra)
        call EinSumup(uEin)

        sequenceQ = .false.   ! only combo move for EEM calculation.

        uInter = eVdw + eCha + EnergyLongRange(hBox) ! Van der Waals + charge + Long Range Correction.
        uIntra = eIntra                              ! Bond Vib + Bending Stratching + Torsion + Non-Bond interaction.
       !uEin                                         ! lattice harmonic vib + orientational vibration.
        eTot   = uInter + uIntra                     ! total energy of system... is it used...?

        write(*,"(3(A,F10.3))") ' inter = ', uInter / dble(molNum) * KtokJ, &
                                ' intra = ', uIntra / dble(molNum) * KtokJ, &
                                ' ein   = ', uEin   / dble(molNum) * kTokJ
        write(1,"(3(A,F10.3))") ' inter = ', uInter / dble(molNum) * KtokJ, &
                                ' intra = ', uIntra / dble(molNum) * KtokJ, &
                                ' ein   = ', uEin   / dble(molNum) * kTokJ

        ! bin searching initialize

        freeQ     = .true.
        mainCycle = 500
        freeCycle = 3000
        acBeu   = 0.
        acBeuSq = 0.
        nBeu    = 0.

        call MakeTarAndRef()

        beuOld  = beta * (xi(0) * uTar + (1.0-xi(0)) * uRef)
        beuNew  = beta * (xi(1) * uTar + (1.0-xi(1)) * uRef)
        testbeu = beuNew - beuOld
        maxbeu  = testbeu
        minbeu  = testbeu

    end subroutine InitializeSimulation





!   ------------------------------------------------
    subroutine SamplingBin()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        call MakeTarAndRef()

        beuOld  = beta * (xi(0) * uTar + (1.0-xi(0)) * uRef)
        beuNew  = beta * (xi(1) * uTar + (1.0-xi(1)) * uRef)

        testbeu = beuNew - beuOld
        if (testbeu > maxbeu) maxbeu = testbeu
        if (testbeu < minbeu) minbeu = testbeu
        acBeu   = acBeu   + testbeu
        acBeuSq = acBeuSq + testBeu**2
        nBeu    = nBeu    + 1.0

    end subroutine SamplingBin





!   ------------------------------------------------
    subroutine MakeBin()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        bin = 6.0 * sqrt(acBeuSq / nBeu - (acBeu / nBeu) ** 2) / dble(devide) ! 6 * sigma approximatley 99% of gaussian distribution
        if (bin > 1) bin = 0.1   ! for very large bin
        write(*,"(A15, F10.7)") 'bin size = ', bin
        write(1,"(A15, F10.7)") 'bin size = ', bin

    end subroutine MakeBin





!   ------------------------------------------------
    subroutine MakeTarAndRef()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        if (einlatticeType == 'mol' .or. einlatticeType == 'MOL') then

            uTar = uInter + uIntra
            uRef = uEin   + uIntra

        else if (einlatticeType == 'atom' .or. einlatticeType == 'ATOM') then  ! for all atom spring

            uTar = uInter + uIntra
            uRef = uEin

        end if

    end subroutine MakeTarAndRef






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
        integer :: i_special
!   ------------------------------------------------

        call InterOldNew (i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)
        call CalculateEnergy(i, 'OLD', eVdw, eCha)
        call EinEnergy (i, rxiOld, ryiOld, rziOld, eEinOld)
        eOld = eVdw + eCha
        call CalculateEnergy(i, 'NEW', eVdw, eCha)
        call EinEnergy (i, rxiNew, ryiNew, rziNew, eEinNew)
        eNew = eVdw + eCha

        dE    = eNew    - eOld
        deEin = eEinNew - eEinOld
        call GetLambda(lambda)   ! get current lambda

        dbe  = beta * (lambda * dE + (1.0 - lambda) * deEin)
        call CheckAccept (dbe, acceptQ)

        if (acceptQ) then
!        if (.false.) then ! freeze condition.
            call UpdateCoordinate(i)

            eTot   = eTot - eOld + eNew
            uInter = uInter - eOld + eNew
            uEin   = uEin - eEinOld + eEinNew

!            nInter = nInter + 1.0_r8

        end if

        ! repeat 1 special movement.
        do i_special = 1, 16

          call SpecialOldNew (i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)
          call EinEnergy (i, rxiOld, ryiOld, rziOld, eEinOld)
          call EinEnergy (i, rxiNew, ryiNew, rziNew, eEinNew)

          ! make energy difference.
          deEin = eEinNew - eEinOld

          ! get current coupling pramater.
          call GetLambda(lambda)

!          dbe  = beta * (1.0 - lambda) * deEin
          dbe  = 1.e3 * deEin

          call CheckAccept (dbe, acceptQ)

!          if (.false.) then
!
!            if (i == 50) then
!              write(*,"('[[¥ÄE_{Ein}]] = ', F10.4)") dbe
!              write(*,"('x_{old} = ', 5F10.3)") rxiold(1:5)
!              write(*,"('x_{new} = ', 5F10.3)") rxinew(1:5)
!              write(*,"('x_{upd} = ', 5F10.3)") rx(i,1:5)
!            end if
!!          end if

          if (acceptQ) then
            call UpdateCoordinate(i)
            uEin = uEin - eEinOld + eEinNew
          end if

        end do

        if (flexQ) then
           !intramolecular motion >>
!           ----------------------------------------------------------------------------------
            if (bondQ .or. bendQ .or. torsQ) then ! only for molecules having at least one intra movement
!           ----------------------------------------------------------------------------------
                call IntraOldNew(i, rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew)

                call CalculateEnergy(i, 'OLD', eVdw, eCha)
                call EnergyIntra (rxiOld, ryiOld, rziOld, eOldIntra)
                call EinEnergy (i, rxiOld, ryiOld, rziOld, eEinOld)
                eOld = eVdw + eCha

                call CalculateEnergy(i, 'NEW', eVdw, eCha)
                call EnergyIntra (rxiNew, ryiNew, rziNew, eNewIntra)
                call EinEnergy (i, rxiNew, ryiNew, rziNew, eEinNew)
                eNew = eVdw + eCha

                call MakeDeltaBetaEnergy()

                call CheckAccept (dbe, acceptQ)

                if (acceptQ) then

                    call UpdateCoordinate(i)

                    uInter = uInter - eOld      + eNew
                    uIntra = uIntra - eOldIntra + eNewIntra
                    uEin   = uEin   - eEinOld   + eEinNew
                    eTot   = uInter + uIntra

!                    nIntra = nIntra + 1.0_r8

                end if

!                tIntra = tIntra + 1.0_r8

!                acE = acE + uInter; nacE = nacE + 1.0_r8
                !call AccumulateData(step, uInter, uIntra, hBox)
                !call AcEinLattice (i)

            end if ! bendQ or torsQ
!           -----------------------------------------------------------------------
        end if ! end flex move

        contains

        subroutine MakeDeltaBetaEnergy()

            implicit none

            real(kind=r8) :: deInt

            if (einLatticeType == 'mol' .or. einLatticeType == 'MOL') then

                dE    = eNew    - eOld
                deEin = eEinNew - eEinOld
                deInt = eNewIntra - eOldIntra
                call GetLambda(lambda)   ! get current lambda

                dbe  = beta * (lambda * dE + (1.0 - lambda) * deEin + deInt)

            else if (einLatticeType == 'atom' .or. einLatticeType == 'ATOM') then

                dE    = eNew + eNewIntra - eOld - eOldIntra ! include intra energy
                deEin = eEinNew - eEinOld
                call GetLambda(lambda)   ! get current lambda

                dbe  = beta * (lambda * dE + (1.0 - lambda) * deEin)

            end if

        end subroutine MakeDeltaBetaEnergy

    end subroutine ComboMovement





!   ------------------------------------------------
    subroutine PrintSimulationResult(cycle_)
!   ------------------------------------------------
        implicit none

        logical, save :: iniQ = .true.
        integer, intent(in) :: cycle_
!   ------------------------------------------------

        if (iniQ) then

            write(*,"(4(A10))") 'cycle', 'u inter', 'u intra', 'u ein'
            write(1,"(4(A10))") 'cycle', 'u inter', 'u intra', 'u ein'

            iniQ = .false.

        end if

        filename = 'view'//IntToCha(partition)//'.view'
        if (mod(cycle_, visualStep) == 0) call MakeConfigView (filename)
        if (mod(cycle_, printStep)  == 0) then

            write(*,"(I10, 3(F10.3))") cycle_, &
                                       uInter / real(molNum) * KtokJ, &
                                       uIntra / real(molNum) * KtokJ, &
                                       uEin   / real(molNum) * KtokJ

            write(1,"(I10, 3(F10.3))") cycle_, &
                                       uInter / real(molNum) * KtokJ, &
                                       uIntra / real(molNum) * KtokJ, &
                                       uEin   / real(molNum) * KtokJ

        end if


    end subroutine PrintSimulationResult




!   ------------------------------------------------
    subroutine FinalizeSimulation()
!   ------------------------------------------------
        implicit none
!   ------------------------------------------------

        write(*,"(A20)") 'check energy'
        write(1,"(A20)") 'check energy'

        write(*,"(3(A20))") 'u inter', 'u intra', 'u ein'
        write(1,"(3(A20))") 'u inter', 'u intra', 'u ein'

        write(*,"(3(F20.3))") uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              uEin   / real(molNum) * KtokJ

        write(1,"(3(F20.3))") uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              uEin   / real(molNum) * KtokJ


        call CalculateEnergy(0, 'OLD', eVdw, eCha)
        if (flexQ) call SumupIntra(eIntra)
        call EinSumup(uEin)

        sequenceQ = .false.   ! only combo move for EEM calculation.

        uInter = eVdw + eCha + EnergyLongRange(hBox) ! Van der Waals + charge + Long Range Correction.
        uIntra = eIntra                              ! Bond Vib + Bending Stratching + Torsion + Non-Bond interaction.
       !uEin                                         ! lattice harmonic vib + orientational vibration.
        eTot   = uInter + uIntra                     ! total energy of system... is it used...?

        write(*,"(3(F20.3))") uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              uEin   / real(molNum) * KtokJ

        write(1,"(3(F20.3))") uInter / real(molNum) * KtokJ, &
                              uIntra / real(molNum) * KtokJ, &
                              uEin   / real(molNum) * KtokJ

        call PrintFreeEnergy()

    end subroutine FinalizeSimulation


end module mod_exp_simulation_routine



