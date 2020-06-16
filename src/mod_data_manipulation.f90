!-----------------------------------------------------
!subroutine ResetData()
!subroutine AccumulateData(cycle_, e1_, e2_, box_)
!subroutine PrintTabulatedData()
!-----------------------------------------------------

module mod_data_manipulation

    use mod_univ_const
    use mod_system_information
    use mod_molecule_information, only : molMass

    implicit none

    private

    integer, parameter, private :: r8 = selected_real_kind(15,300)

    integer, parameter :: UNIT_DATA = 20
    real(kind=r8), parameter :: ZERO = 0.0_r8    

    real(kind=r8), dimension(20) :: dataSet
    integer :: nData

    integer :: nCycle
    real(kind=r8) :: acEInter, acEIntra, nacEInter, nacEIntra, nacBox, &
                     acVolume, nacVolume
    real(kind=r8), dimension(3,3) :: acBox

    public :: PrintTabulatedData, ResetData, AccumulateData

    contains    

!   -------------------------------------------------------
    subroutine ResetData()
!   -------------------------------------------------------
        implicit none
!   -------------------------------------------------------

        acEInter  = ZERO
        acEIntra  = ZERO
        nacEInter = ZERO 
        nacEIntra = ZERO
        acBox     = ZERO
        nacBox    = ZERO
        acVolume  = ZERO
        nacVolume = ZERO

    end subroutine ResetData





!   -------------------------------------------------------
    subroutine AccumulateData(cycle_, e1_, e2_, box_)
!   -------------------------------------------------------
        use mod_math_tool, only : Det        

        implicit none
  
        integer, intent(in) :: cycle_
        real(kind=r8), intent(in) :: e1_, e2_ 
        real(kind=r8), dimension(3,3), intent(in) :: box_
 
!   -------------------------------------------------------

        nCycle = cycle_

        nacEInter = nacEInter + 1
        nacEIntra = nacEIntra + 1
        nacBox    = nacBox    + 1
        nacVolume = nacVolume + 1

        acEInter = acEInter + e1_
        acEIntra = acEIntra + e2_
        acBox    = acBox  +  box_ 
        acVolume = acVolume + Det (box_,3)
       
    end subroutine AccumulateData






! -------------------------------------------------------
  subroutine PrintTabulatedData()
! -------------------------------------------------------
    implicit none

    logical, save :: initialQ = .true.

    real(kind=r8), parameter :: dTocm3 = 1.66054

    real(kind=r8) :: avEInter, avEIntra, avVolume, nMol, dens 
    real(kind=r8), dimension(3,3) :: mat, avBox
    real(kind=r8) norm1, norm2, norm3, ang12, ang23, ang31
! --------------------------------------------------------

    open(unit=UNIT_DATA, file='tabulated_data.txt', position = 'append')
        
    if (initialQ) then 

      ! not classified
      !open(unit=UNIT_DATA, file='tabulated_data.txt', position = 'append')
      write(UNIT_DATA,"()")
      write(UNIT_DATA,"()")
      write(UNIT_DATA, "(2A10, 2A20, 2A12, A15, 20A10)") &
      & 'cycle', 'N', 'EInter(kJ/mol)', 'EIntra(kJ/mol)', &
      & 'temper(K)', 'press(bar)', 'dens(g/cm^3)',   &
      & 'a(（)', 'b(（)', 'c(（)',                   &
      & 'メ(rad)' ,'モ(rad)', 'ャ(rad)' 
      write(UNIT_DATA,"()")
      initialQ = .false.

    end if

    avEInter = acEInter / nacEInter
    avEIntra = acEIntra / nacEIntra
    avBox    = acBox    / nacBox
    avVolume = acVolume / nacVolume
    nMol     = dble(molNum)

    avEInter = avEInter / nMol * KtokJ
    avEIntra = avEIntra / nMol * KtokJ

    mat = avBox

    dens  = nMol / avVolume * dTocm3 * molMass

    norm1 = sum(mat(:,1)*mat(:,1))**0.5_r8
    norm2 = sum(mat(:,2)*mat(:,2))**0.5_r8
    norm3 = sum(mat(:,3)*mat(:,3))**0.5_r8
        
    ang12 = acos(sum(mat(:,1)*mat(:,2))/norm1/norm2) * RAD_TO_DEG
    ang23 = acos(sum(mat(:,2)*mat(:,3))/norm2/norm3) * RAD_TO_DEG
    ang31 = acos(sum(mat(:,3)*mat(:,1))/norm3/norm1) * RAD_TO_DEG

    write(UNIT_DATA, "(2I10, 2F20.4, 2F12.4, F15.4, 20F10.4)") nCycle, molNum, avEInter, avEIntra, temper, pressure, dens, &
    norm1, norm2, norm3, ang23, ang31, ang12

    close(UNIT_DATA)
        
  end subroutine PrintTabulatedData

end module mod_data_manipulation