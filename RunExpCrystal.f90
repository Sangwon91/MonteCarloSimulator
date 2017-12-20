!   -----------------------------------------------------------------
    subroutine RunExpCrystal
!   -----------------------------------------------------------------
        use mod_exp_simulation_routine

        implicit none
!   -----------------------------------------------------------------

        call ReadInput()
        call InitializeSimulation()

        ! ==== bin search ====

        do binCycle = 1,  freeCycle + mainCycle

!            write(*, *) "[[assertion]]: the modifications are applied (bin cycle)."
            call MoveMolecule()
            if (binCycle >= freeCycle) call SamplingBin()

            call PrintSimulationResult(binCycle)

            call GetEndSignal
            if (endSign) go to 119

        end do ! bin-loop, to find bin of energy profile

        call MakeBin()

        ! ==== end of bin search ====

        eeCycle = 1

        do while(.true.)

            call MoveMolecule()

            call MakeTarAndRef()
            call ExpMCSolid(SOLID_TYPE, uTar, uRef, eeCycle)

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

 119    continue

        call FinalizeSimulation()

    end subroutine RunExpCrystal
