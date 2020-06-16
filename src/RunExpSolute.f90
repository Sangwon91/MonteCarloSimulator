!   -----------------------------------------------------------------
    subroutine RunExpSolute()
!   -----------------------------------------------------------------
        use mod_exp_solute_routine

        implicit none
!   -----------------------------------------------------------------

        eeCycle = 1

        call ReadInput()
        call InitializeSimulation()

        do while(.true.)

            call MoveMolecule()
            call ExpMCSolid(FLUID_TYPE, eTot, uInter, eeCycle)

            if (weightOK) then

                weightOK = .false.
                eeCycle  = 0

            end if

            if (.not.adjustModeQ) then

                if (eeCycle == mainStep) exit

            end if

            call PrintSimulationResult(eeCycle)

            call GetEndSignal
            if (endSign) exit

            eeCycle = eeCycle + 1

        end do

        call FinalizeSimulation()

    end subroutine RunExpSolute
