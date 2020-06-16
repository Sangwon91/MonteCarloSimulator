subroutine RunConvention

    use mod_simulation_routine

    implicit none

    call ReadInput()
    call InitializeSimulation()

    do step = 1, preStep + mainStep

        call MoveMolecule()
        call MoveSimulationBox()
        call PrintSimulationResult()

        call ResetForMainloop() ! if end of prestep.

        call AdjustNVTMax()
        call AdjustNPTMax()

        call GetEndSignal()
        if (endSign) exit

        call CalRemainTime (step, preStep + mainStep)

    end do

    call FinalizeSimulation()

end subroutine RunConvention
