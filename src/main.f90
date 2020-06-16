!===== included module =====
include 'mod_time_tool.f90'
include 'mod_math_tool.f90'
include 'mod_random_tool.f90'
include 'mod_character_tool.f90'

include 'mod_module_set.f90'
include 'mod_simulation_routine.f90'
include 'mod_exp_simulation_routine.f90'
include 'mod_exp_solute_routine.f90'

!===== main program =====
program Main

    ! ---------------------- modules ---------------------------
    use mod_version_information
    use mod_time_tool
    use mod_system_information !, only : GetSysInf, simulType
    use mod_molecule_information, only : MolDeallocate
    use mod_config_information, only : ConfigDeallocate
    use mod_movement_information
    use mod_einstein_core, only : EinDeallocate, TurnOnCheckProb
    use mod_expans_solid, only : ExpSolDeallocate
    use mod_RDF_core, only : RDFDeallocate
    use mod_neighbor_list, only : DeallocateNeighborList
    use mod_program_operation
    ! ----------------------------------------------------------

    implicit none

    character(len=31) :: dc = 'system.txt'
    integer :: i
    character(len=32) arg

    call ProcessStateLog(RUNNING)

    !===== consol setting =====
    call system('color 0A')
    call system("mode con: cols=120 lines=2500")

    !===== get simulation type =====

    call GetSysInf(dc)

    do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (arg == '--version') then
        simulType = 'version'
        exit
      end if
      if (arg == '--check-prob') then
        call TurnOnCheckProb()
      end if
    end do

    !===== version reading  =====
    if (simulType == 'version') then

        call PrintVersionInformation()
        go to 1

    end if

    !===== start time calculation =====
    call cpu_time(time1) ! time calculation

    !===== choose simulation type =====
    select case (simulTypeENUM)

        case (NPT, NVT)

            call RunConvention

        case (EXP_CRYSTAL)

            call RunExpCrystal

        case (EXP_SOLUTE)

            call RunExpSolute

        case default

    end select

    !===== end time calculation =====
    call cpu_time(time2)

    call SecToHour (time2-time1)
    call SecToHour (time2-time1,1)

    !===== destruinclude 'mod_math_tool.f90'ctor =====
    call MolDeallocate
    call ConfigDeallocate
    call EinDeallocate
    call ExpSolDeallocate
    call RDFDeallocate
    call DeallocateNeighborList

    write(*,"('- END -')")
    write(1,"('- END -')")

    close(1)

1    continue

    call ProcessStateLog(DONE)

end program Main

!===== include some subroutine to use =====

include 'RunConvention.f90'
include 'RunExpCrystal.f90'
include 'RunExpSolute.f90'
