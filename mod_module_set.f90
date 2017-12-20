! ※ 프로그래밍 스타일 가이드
!
!============================================================

module mod_version_information

    implicit none

    integer, parameter :: r8 = selected_real_kind(15,300)

    public

    character(len=5), parameter :: version = '7.1.9'

    private :: r8

    contains

    subroutine PrintVersionInformation

    write(*,"(A, A5/)") 'current version : ', trim(version)

    write(*,*) 'ver 6.2.1 - configuration reading automatically adjusted'
    write(*,*) '          - nvt : read npt_configReault.txt (before : config.txt)'
    write(*,*) '          - eem : read nvt_configReault.txt (before : config.txt)'
    write(*,*) 'ver 6.2.2 - angle distribution bin is expanded to 1000'
    write(*,*) 'ver 6.2.3 - maximum rcut calculating in noncubic box'
    write(*,*) 'ver 6.3.0 - module of temper via expanded ensemble is added'
    write(*,*) '          - first line of expanded.txt is Tem(target) Tem(reference)'
    write(*,*) 'ver 6.3.1 - verseion information module is added'
    write(*,*) 'ver 6.3.2 - The orientational external field for linear molecule is'
    write(*,*) '            became [1 - cos(2w)]'
    write(*,*) 'ver 6.3.3 - angle distribution bin is reduced to 500'
    write(*,*) 'ver 6.3.4 - pass'
    write(*,*) 'ver 6.3.5 - in angle distribution module, int function is replace  to floor function'
    write(*,*) 'ver 7.0.0 - flexible module is added'
    write(*,*) 'ver 7.0.1 - nonbonded energy fixed.. and tabulated data is now supproted'
    write(*,*) 'ver 7.0.2 - vibrational motion is added, rescale bond length is deleted'
    write(*,*) 'ver 7.0.3 - random torsional motion is added with 5 % frequency with onoff shcame'
    write(*,*) 'ver 7.0.6 - fix1mol, fix1Atom are added, einstein module is more simplified'
    write(*,*) '            ein data accumulator is added, frame atom reader is added'
    write(*,*) 'ver 7.0.7 - adjust module fixed'
    write(*,*) 'ver 7.1.0 - Flexible Williams and Dreiding force field are now surpported.'
    write(*,*) 'ver 7.1.1 - Flexible Williams and Dreiding force field bug fixed - prevent too low distance.'
    write(*,*) 'ver 7.1.2 - process state log is added.'
    write(*,*) 'ver 7.1.3 - config style updated.'
    write(*,*) 'ver 7.1.4 - unitcell angle corrected'
    write(*,*) 'ver 7.1.5 - integer enumerator bug in intra energy is fixed'
    write(*,*) 'ver 7.1.6 - expanded ensemble open algorithm is changed.'
    write(*,*) 'ver 7.1.7 - some bugs are fixed.'
    write(*,*) 'ver 7.1.8 - enum, binc(c) is used for enumeration instead of '
    write(*,*) '            integer parameter for simulation type and molecule type.'
    write(*,*) 'ver 7.1.9 - type error in module of neighbor list is fixed'
    write(*,*) '            for nvt simulation... pressure be setted to -1'
    write(*,*) ' '
        call system("pause")

    end subroutine PrintVersionInformation

end module mod_version_information

!============================================================

module mod_program_operation

    implicit none

    private

    integer, parameter :: r8 = selected_real_kind(15,300)
    integer ,parameter :: READY = 1, &
                          RUNNING = READY + 1, &
                          DONE = RUNNING + 1

    logical            :: endSign

    public :: READY, RUNNING, DONE
    public :: endSign, GetEndSignal, ProcessStateLog

    contains

    subroutine GetEndSignal

        implicit none

        integer, parameter :: U_END = 11

        endSign = .false.

        open(U_END, file='end.end', status='old', err=1)

          read(U_END, *, end=2, err=3) endSign

        close(U_END, err=4)

1       continue
2       continue
3       continue
4       continue

    end subroutine GetEndSignal

    subroutine ProcessStateLog(state)

        implicit none

        integer, parameter :: unit_now = 23
        integer, intent(in) :: state

        open(unit=unit_now, file='process_state.log')

        select case (state)

            case (READY)

                write(unit_now,"('ready')")

            case (RUNNING)

                write(unit_now,"('running')")

            case (DONE)

                write(unit_now,"('done')")

            case default

        end select

        close(unit_now)

    end subroutine ProcessStateLog

end module mod_program_operation

!============================================================

module mod_periodic_table

    implicit none

    public

    character(len=2), parameter, dimension(20) :: atomName = &
    & (/'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
    &   'Na', 'Ma', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca'/)

end module mod_periodic_table

!============================================================

module mod_univ_const

    implicit none

    integer, parameter, private :: r8 = selected_real_kind(15,300) ! real type & 8 byte variable: r8

    public !무조건 외부로 사용할 상수들이라 defualt값이 public, for global variables.

    !real(kind=r8), parameter :: RGAS  = 1.0_r8        ! for reduced unit
    real(kind=r8), parameter :: RGAS  = 0.008314469_r8 ! gas const in unit of kJ/mol/K
    real(kind=r8), parameter :: KtokJ = RGAS           ! for energy unit K to kJ/NA (NA : avogadro number)
    real(kind=r8), parameter :: PI    = 3.141592653589793_r8, TWOPI = 2.0_r8 * PI
    real(kind=r8), parameter :: CHARGE_CONST = 167100.95662058413_r8 ! 1 / ( 4 PI E0 ) in Kelvin unit
    real(kind=r8), parameter :: RAD_TO_DEG = 180.0_r8/PI, DEG_TO_RAD = PI/180.0_r8
    real(kind=r8), parameter :: BAR_A3toK = 0.007242963171560301_r8  ! bar * A^3 / k = K
    !real(kind=r8), parameter :: BAR_A3toK = 1.0_r8    ! for reduced unit

end module mod_univ_const





!============================================================



module mod_filesystem

    implicit none

    public

    contains

!   ---------------------------------------------------------
    subroutine copy_file(in_file, out_file)
!   ---------------------------------------------------------
        implicit none

        character(len=*), intent(in) :: in_file, out_file

        integer :: in_unit = 111, out_unit = 222

        character(len=1000) :: line ! buffer
!   ---------------------------------------------------------
        open(unit=in_unit, file=in_file)
        open(unit=out_unit, file=out_file)
        !write(*,*) 'copy start'

        do while (.true.)

            read(in_unit, "(A1000)", err=100, end=101) line
            write(out_unit, "(A)") trim(line)

        end do

    100 continue
    101 continue

        close(out_unit)
        close(in_unit)

    end subroutine copy_file
!   ---------------------------------------------------------


end module mod_filesystem


!============================================================
include 'mod_system_information.f90'
include 'mod_molecule_information.f90'
include 'mod_config_information.f90'
include 'mod_neighbor_list.f90'
include 'mod_movement_information.f90'
include 'mod_energy_core.f90'
include 'mod_einstein_core.f90'
include 'mod_expans_solid.f90'
include 'mod_RDF_core.f90'
include 'mod_data_manipulation.f90'


!include 'mod_expans_liquidSep.module'
!include 'mod_expans_solid_temper.f90'  ! is not fully tested.
!include 'mod_expans_liquidTemper.f90'  ! is not fully tested.
!================================================================
