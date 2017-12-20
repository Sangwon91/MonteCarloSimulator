! subroutine GetSysInf(filename)
! subroutine PrintSysInf(unit_)

module mod_system_information

! system information module have many data
! 1. **
! 2. number of molecule
! 3. **
! 4. simulation type (NVT, NPT, expanded ensamble.. etc)
! 5. temperature, pressue
! 6. **
! 7. **
! 8. potential cutoff
! 9. print operation

    use mod_version_information

    implicit none

    public

    integer, parameter:: r8 = selected_real_kind(15,300)

    integer, parameter :: UNIT_SYS = 10
    integer :: molNum ! number of molecules
    logical :: ewaldQ, periodicQ, cubicBoxQ, neighborListQ, boxChangeQ, &
               fixCOMQ, fixAtom1Q, fixMol1Q, randomTorsionQ
    character(len=10) :: einLatticeType   ! mol (C.O.M of mol) or atom
    character(len=31) :: simulType
    character(len=61) :: configResult, currentConfig
    integer :: fixTrans, fixRot, fixBox, printStep, visualStep, sequenceStep = 1000
    integer :: mainStep, preStep



!    integer, parameter :: NVT         = 1, &
!                          NPT         = NVT + 1, &
!                          EXP_SOLUTE  = NPT + 1, &
!                          EXP_CRYSTAL = EXP_SOLUTE + 1

    enum, bind(c)

      enumerator :: NVT = 1
      enumerator :: NPT
      enumerator :: EXP_SOLUTE
      enumerator :: EXP_CRYSTAL

    end enum

    integer :: simulTypeENUM

    real(kind=r8) :: temper, pressure, rCut, rvCut, beta

    private :: r8

    contains

    subroutine GetSysInf (filename)

        implicit none

        character(len=31), intent(in) :: filename

        open(UNIT_SYS, file=filename)

            read(UNIT_SYS,*) simulType
            read(UNIT_SYS,*) molNum
            read(UNIT_SYS,*) temper
            read(UNIT_SYS,*) pressure  ! in unit of bar
            read(UNIT_SYS,*) preStep
            read(UNIT_SYS,*) mainStep
            read(UNIT_SYS,*) rCut
            !read(UNIT_SYS,*) fixTrans
            !read(UNIT_SYS,*) fixRot
            !read(UNIT_SYS,*) fixBox
            read(UNIT_SYS,*) printStep
            read(UNIT_SYS,*) visualStep
            !read(UNIT_SYS,*) sequenceStep
            read(UNIT_SYS,*) ewaldQ
            read(UNIT_SYS,*) periodicQ
            read(UNIT_SYS,*) cubicBoxQ
            read(UNIT_SYS,*) neighborListQ
            read(UNIT_SYS,*) rvCut
            read(UNIT_SYS,*) boxChangeQ
            read(UNIT_SYS,*) fixCOMQ
            read(UNIT_SYS,*) fixMol1Q
            read(UNIT_SYS,*) fixAtom1Q
            read(UNIT_SYS,*) einLatticeType
            read(UNIT_SYS,*) randomTorsionQ

        close(UNIT_SYS)

        select case (simulType)

            case ('nvt')

                simulTypeENUM = NVT

            case ('npt')

                simulTypeENUM = NPT

            case ('exp_solute')

                simulTypeENUM = EXP_SOLUTE

            case ('exp_crystal')

                simulTypeENUM = EXP_CRYSTAL

            case default

                stop 'invalid simulation type'

        end select

        if (.not. boxChangeQ) pressure = -1.0

        beta = 1.0_r8 / temper
        !configResult  = trim(simulType)//'_configReault.txt'
        !currentConfig = trim(simulType)//'_currentConfig.txt'
        configResult  = 'config.txt'
        currentConfig = 'current_config.txt'

    end subroutine GetSysInf

    !********************

    subroutine PrintSysInf (unit_)

        implicit none

        integer, optional, intent(in) :: unit_

        if (present(unit_)) then

            write(unit_,"('molecule simulation tool ver ', A5/)") trim(version)

            write(unit_,"(A/)")      'system information -> '

            write(unit_,"(A,A15)")   'simulation type    = ', simulType
            write(unit_,"(A,I15)")   'number of molecule = ', molNum
            write(unit_,"(A,F15.3)") 'temperature        = ', temper
            write(unit_,"(A,F15.3)") 'pressure           = ', pressure
            write(unit_,"(A,I15)")   'previous step      = ', preStep
            write(unit_,"(A,I15)")   'main step          = ', mainStep
            write(unit_,"(A,F15.3)") 'Rcut               = ', rCut
            !write(unit_,"(A,I15)")   'fix translation    = ', fixTrans
            !write(unit_,"(A,I15)")   'fix rotation       = ', fixRot
            !write(unit_,"(A,I15)")   'fix box            = ', fixBox
            write(unit_,"(A,I15)")   'print step         = ', printStep
            write(unit_,"(A,I15)")   'visual step        = ', visualStep
            !write(unit_,"(A,I15)")   'sequence step      = ', sequenceStep
            write(unit_,"(A,L15)")   'ewald summation    = ', ewaldQ
            write(unit_,"(A,L15)")   'periodic           = ', periodicQ
            write(unit_,"(A,L15)")   'cubic box          = ', cubicBoxQ
            write(unit_,"(A,L15)")   'neighbor list      = ', neighborListQ
            write(unit_,"(A,F15.3)") 'rvcut              = ', rvcut
            write(unit_,"(A,L15)")   'box change         = ', boxChangeQ
            write(unit_,"(A,L15)")   'fix C.O.M          = ', fixCOMQ
            write(unit_,"(A,L15)")   'fix mol 1          = ', fixMol1Q
            write(unit_,"(A,L15)")   'fix atom 1         = ', fixAtom1Q
            write(unit_,"(A,A15)")   'ein lattice type   = ', einLatticeType
            write(unit_,"(A,L15)")   'random torsion     = ', randomTorsionQ
            write(unit_,"()")

        else

            write(*,"('molecule simulation tool ver ', A5/)") trim(version)

            write(*,"(A/)")      'system information -> '

            write(*,"(A,A15)")   'simulation type    = ', simulType
            write(*,"(A,I15)")   'number of molecule = ', molNum
            write(*,"(A,F15.3)") 'temperature        = ', temper
            write(*,"(A,F15.3)") 'pressure           = ', pressure
            write(*,"(A,I15)")   'previous step      = ', preStep
            write(*,"(A,I15)")   'main step          = ', mainStep
            write(*,"(A,F15.3)") 'Rcut               = ', rCut
            !write(*,"(A,I15)")   'fix translation    = ', fixTrans
            !write(*,"(A,I15)")   'fix rotation       = ', fixRot
            !write(*,"(A,I15)")   'fix box            = ', fixBox
            write(*,"(A,I15)")   'print step         = ', printStep
            write(*,"(A,I15)")   'visual step        = ', visualStep
            !write(*,"(A,I15)")   'sequence step      = ', sequenceStep
            write(*,"(A,L15)")   'ewald summation    = ', ewaldQ
            write(*,"(A,L15)")   'periodic           = ', periodicQ
            write(*,"(A,L15)")   'cubic box          = ', cubicBoxQ
            write(*,"(A,L15)")   'neighbor list      = ', neighborListQ
            write(*,"(A,F15.3)") 'rvcut              = ', rvcut
            write(*,"(A,L15)")   'box change         = ', boxChangeQ
            write(*,"(A,L15)")   'fix C.O.M          = ', fixCOMQ
            write(*,"(A,L15)")   'fix mol 1          = ', fixMol1Q
            write(*,"(A,L15)")   'fix atom 1         = ', fixAtom1Q
            write(*,"(A,A15)")   'ein lattice type   = ', einLatticeType
            write(*,"(A,L15)")   'random torsion     = ', randomTorsionQ
            write(*,"()")

        end if

    end subroutine PrintSysInf

end module mod_system_information
