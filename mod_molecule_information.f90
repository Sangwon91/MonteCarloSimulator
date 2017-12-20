  ! imported by JAVA style: include 'module_name.module'

  !-------------------------------------------------------------------------------------------------------------
  !subroutine GetMolInf (fileName)        ! read input data
  !subroutine CheckMotions()              ! rotateQ, bondQ, etc...
  !subroutine ReadAtomType()              ! read atomic parameter (eps, sig, etc...)
  !subroutine ReadAtomIndex()             ! read type of atom in non-bond, vib, bend, tors motion
  !                                       ! make myTypeIn Series... (myTypeInBond, myTypeInBend, etc...)
  !subroutine ReadBond()                  ! read and make bond indices
  !subroutine ReadBondType()              ! read bond type
  !subroutine ReadBondParameter()         ! read bond parameters
  !subroutine ReadBendType()              ! read bend type
  !subroutine ReadBendParameter()         ! read bend parameters
  !subroutine ReadTorsionType()           ! read torsion type
  !subroutine ReadTorsionParameter()      ! read torsion parameters
  !subroutine ReadMoleculeCoordinate()    ! read molecule coordinate

  !subroutine ComToOrigin()               ! move molecule COM to origin
  !subroutine SearchBond()                ! make myBondType
  !subroutine SearchBending()             ! make myBendType, use after call SearchBond
  !subroutine SearchTorsion()             ! make myTorsionType, use after call SearchBending

  !subroutine MakeNeighborQ()             ! make neighborQ -> true for 1-2, 1-3 connection
  !subroutine MakeScalingQ()              ! make scalingQ  -> true for 1-4 connection
  !subroutine MakeAdjacentIndices()
  !subroutine MakeMolTree()
  !subroutine SetRelationship()
  !subroutine MakeMolView()               ! make view file for molecule
  !subroutine GroupingMolecule()
  !
  !subroutine DeallocateGroup()
  !subroutine MolDeallocate()
  !
  !subroutine PrintCheckedMotions(unit_)  ! print rotateQ, bondQ, etc...
  !subroutine PrintBondInformation(unit_)
  !subroutine PrintBendInformation(unit_)
  !subroutine PrintTorsionInformation(unit_)
  !subroutine PrintNeighborQ(unit_)
  !subroutine PrintScalingQ(unit_)
  !subroutine PrintAdjacentIndices(unit_)
  !subroutine PrintMolTree(unit_)
  !subroutine PrintRelationship(unit_)
  !subroutine PrintMolInf (unit_)
  !subroutine PrintInertia (unit_)
  !subroutine PrintGroup(unit_)
  !
  !function GiveBondType(b1, b2) result(bondType_)               ! search bond type at given atom indices
  !function GiveBendType(b1, b2, b3) result(bendType_)           ! search bend type at given atom indices
  !function GiveTorsionType(t1, t2, t3, t4) result(torsionType_) ! search torsion type at given atom indices
  !-------------------------------------------------------------------------------------------------------------

  module mod_molecule_information

  use mod_periodic_table

  implicit none

  private

  integer, parameter:: r8 = selected_real_kind(15,300)

  logical           :: flexQ = .false., rotateQ, bondQ, bendQ, torsQ

  integer, parameter:: UNIT_MOL = 10

  !integer, parameter :: OPLSAA = 1, &
  !  FW = OPLSAA + 1, & ! Flexible William potential
  !DREIDING = FW + 1  ! Dreiding

  enum, bind(c)

    enumerator :: OPLSAA = 100  ! OPLS-AA
    enumerator :: FW            ! Flexible William potential
    enumerator :: DREIDING      ! Dreiding

  end enum

  integer :: potenTypeENUM = 0

  integer, parameter :: MAX_BOND_TYPE = 20, &
    MAX_BEND_TYPE = 20, &
    MAX_TORS_TYPE = 20, &
    MAX_ATOM_TYPE = 10




  integer           :: nAtoms, maxCate ! nAtoms : number of atoms, maximum number of categoty
  integer           :: nType             ! number of atomic types
  integer           :: nBond, nBondType             ! number of atomic bonds
  integer           :: nBend, nBendType             ! number of atomic bends
  integer           :: nTorsion, nTorsionType         ! number of atomic torsions
  integer           :: frameAtom1, frameAtom2                 ! body frame의 기준이 되는 Atom
  logical           :: linearQ                                ! linear molecule = true, else false
  logical           :: allocateQ                              ! allocated array = true, esle false
  logical           :: getMolInfQ = .false.    ! molecule information을 얻었는지 체크.

  character(len=10)  :: potenType, mixingRule, moleculeType
  ! mixingRule a or A arithmetic mean, g or G geometric mean
  ! moleculeType: rigid, flexible or RIGID, FLEXIBLE

  integer, dimension(MAX_BOND_TYPE, 2) :: bondType
  integer, dimension(MAX_BEND_TYPE, 3) :: bendType
  integer, dimension(MAX_TORS_TYPE, 4) :: torsionType

  real(kind=r8), parameter :: scale14 = 0.5_r8

  real(kind=r8), dimension(MAX_BOND_TYPE) :: bondLength, kBond
  real(kind=r8), dimension(MAX_BEND_TYPE) :: kBend, theta0, theta0Deg, C_theta_DD, cos_theta0
  real(kind=r8), dimension(MAX_TORS_TYPE) :: v1, v2, v3, v4, v5, v6 ! for beding & torsion
  real(kind=r8), dimension(MAX_TORS_TYPE) :: d1, d2, d3, d4, d5, d6

  real(kind=r8)                  :: molMass
  real(kind=r8), dimension(3,3)  :: InertiaTensor
  real(kind=r8), dimension(20,4) :: molCate ! molecule category

  ! molCate(i,j) : 1. i = category index
  !                2. j = 1 -> epsilon
  !                       2 -> sigma
  !                       3 -> charge
  !                       4 -> atom number in category

  real(kind=r8), dimension(MAX_ATOM_TYPE, 6) :: atomType ! type of atom

  ! atomType(i,j): 1. i = type index
  !                 2. j = 1 -> atomic number
  !                        2 -> epsilon
  !                       3 -> sigma
  !                        4 -> charge
  !                        5 -> atomic mass

  real(kind=r8), dimension(20) :: atomTypeNum ! atomTypeNum(i): return number of atom in type i

  real(kind=r8), dimension(MAX_ATOM_TYPE) :: eps, sig, mass, charge, color

  real(kind=r8), dimension(MAX_ATOM_TYPE,MAX_ATOM_TYPE) :: epsij, sigij, chaij
  real(kind=r8), dimension(MAX_ATOM_TYPE,MAX_ATOM_TYPE) :: A_FW, B_FW, C_FW ! flexible william parameters.

  real(kind=r8), dimension(MAX_ATOM_TYPE,MAX_ATOM_TYPE) :: zeta_DD, D0_DD, R0_DD          ! Dreiding parameters.
  real(kind=r8), dimension(MAX_ATOM_TYPE,MAX_ATOM_TYPE) :: A_DD, B_DD, C_DD ! DD : Dreiding
  real(kind=r8), dimension(MAX_ATOM_TYPE,MAX_ATOM_TYPE) :: MIN_DISTANCE_DD

  ! ===== allocatable vaiables =====

  integer, dimension(:), allocatable :: myType, &      ! give my atomic Type index
  myAtomNum, &       ! give my atomic number
  myTypeInBond, &    ! to simplify ...
  myTypeInBend, &    ! to simplify input data..
  myTypeInTorsion    ! to simplify input data..
  integer, dimension(:,:), allocatable :: bondIndices, bendIndices, torsionIndices, adjacentIndices
  integer, dimension(:), allocatable :: myBondType, myBendType, myTorsionType
  integer, dimension(:,:,:), allocatable :: molTree

  integer, dimension(:), allocatable :: index1InBond, index2InBond, &
    &                                     index1InBend, index2InBend, index3InBend, &
    &          index1InTorsion,  index2InTorsion, index3InTorsion, index4InTorsion, &
    &                     bond12InBend, bond23InBend, direc12InBend, direc23InBend, &
    &                            bond12InTorsion, bond23InTorsion, bond34InTorsion, &
    &                         direc12InTorsion, direc23InTorsion, direc34InTorsion, &
    &       bend123InTorsion, bend234InTorsion, direc123InTorsion, direc234InTorsion

  logical, dimension(:,:), allocatable, public :: neighborQ, scalingQ

  real(kind=r8), dimension(:), allocatable :: molCrdX, molCrdY, molCrdZ    ! molecule structure coordinate



  !   -------------------------------------------------
  type Group  ! derived array for pointer operation
    !   -------------------------------------------------
    integer, dimension(:), pointer :: member !

  end type Group

  type(Group), dimension(:), allocatable :: localGroup

  integer :: nLocalGroup
  character(len=5) :: rcutType  ! rcutType can be 'mol', 'atom' or 'group'





  ! ===== public variables =====
  public :: nAtoms, color, molCrdX, molCrdY, molCrdZ, potenType, epsij, sigij, chaij, &
    frameAtom1, frameAtom2, maxCate, molCate, mixingRule, eps, sig, charge, molMass, linearQ
  public :: nType, nBond, nBend, nTorsion, moleculeType, atomType, myType, myAtomNum, &
    nBendType, nTorsionType, atomTypeNum
  public :: bondIndices, bendIndices, torsionIndices
  public :: index1InBend, index2InBend, index3InBend, molTree, mass
  public :: index1InTorsion,  index2InTorsion, index3InTorsion, index4InTorsion
  public :: scale14, index1InBond, index2InBond, bondLength, bond12InBend, bond23InBend, &
    myBendType, myTorsionType, direc12InBend, direc23InBend, bond12InTorsion, &
    bond23InTorsion, bond34InTorsion, bend123InTorsion, bend234InTorsion, &
    direc12InTorsion, direc23InTorsion, direc34InTorsion, myBondType, &
    direc123InTorsion, direc234InTorsion, &
    flexQ, rotateQ, bondQ, bendQ, torsQ, kBond
  public :: A_FW, B_FW, C_FW
  public :: A_DD, B_DD, C_DD, zeta_DD, R0_DD, D0_DD, MIN_DISTANCE_DD

  public :: potenTypeENUM, OPLSAA, FW, DREIDING

  public :: localGroup, nLocalGroup

  ! ===== public methods =====
  public :: GetMolInf, PrintMolInf, MakeMolView, MolDeallocate
  public :: kBend, theta0, v1, v2, v3, v4, v5, v6, &
    d1, d2, d3, d4, d5, d6, &
    C_theta_DD, cos_theta0

  ! ===== public derived type =====
  public :: Group

  contains


  !********************

  subroutine GetMolInf (fileName)    !  get molecule information

  character(len=31), intent(in) :: fileName

  integer :: i, j, dNum   ! dNum : dummy number
  logical :: echoQ = .true.

  linearQ = .false.

  open(UNIT_MOL, file=fileName)

  read(UNIT_MOL, *) linearQ       ! T: linear molecule, F: non-linear molecule.
  read(UNIT_MOL, *) potenType     ! OPLSAA, Flexible William and etc...

  ! potential type enumerator =============

  select case (potenType)

  case ('OPLSAA', 'oplsaa')

    potenTypeENUM = OPLSAA

  case ('FW', 'fw')

    potenTypeENUM = FW

  case ('DREIDING', 'dreiding')

    potenTypeENUM = DREIDING
    if (echoQ) write(*,*) DREIDING, potenTypeENUM

    case default

    stop 'invalid potential type'

  end select

  ! ======================================

  read(UNIT_MOL, *) moleculeType  ! RIGID, FLEXIBLE
  if (moleculeType == 'FLEXIBLE') flexQ = .true.
  if (echoQ) read(UNIT_MOL, *) rcutType      ! atom, mol, group

  read(UNIT_MOL, *) nType            ! number of atom type, C, F, N, etc...
  if (nType > MAX_ATOM_TYPE) stop 'atom type is larger than MAX_ATOM_TYPE'
  call ReadAtomType()                !!! have to be implementied
  if (echoQ) write(*,*) 'mol beg', 1

  read(UNIT_MOL, *) nAtoms        ! number of atoms in molecule.
  call ReadAtomIndex()            ! give atom a type.
  if (echoQ) write(*,*) 2

  read(UNIT_MOL, *) nBond            ! number of bond in molecule.
  call ReadBond()                    !!! have to be implementied
  if (echoQ) write(*,*) 3

  read(UNIT_MOL, *) nBondType
  if (nBondType > MAX_BOND_TYPE) stop 'bond type is larger than MAX_BOND_TYPE'
  call ReadBondType()
  if (echoQ) write(*,*) 4
  call ReadBondParameter()
  if (echoQ) write(*,*) 5

  read(UNIT_MOL, *) nBendType        ! number of type of bend in molecule.
  if (nBendType > MAX_BEND_TYPE) stop 'bend type is larger than MAX_BEND_TYPE'
  call ReadBendType()
  if (echoQ) write(*,*) 6
  call ReadBendParameter()
  if (echoQ) write(*,*) 7

  read(UNIT_MOL, *) nTorsionType    ! number of type of torsion in molecule.
  if (nTorsionType > MAX_TORS_TYPE) stop 'torsion type is larger than MAX_TORS_TYPE'
  call ReadTorsionType()
  if (echoQ) write(*,*) 8
  call ReadTorsionParameter()
  if (echoQ) write(*,*) 9

  call ReadMoleculeCoordinate()
  if (echoQ) write(*,*) 10

  read(UNIT_MOL, *) frameAtom1, frameAtom2
  if (echoQ) write(*,*) 11

  if (echoQ) write(*,"('searchBond()')")
  call SearchBond()
  if (echoQ) write(*,"('searchBending()')")
  call SearchBending()
  if (echoQ) write(*,"('searchTorsion()')")
  call SearchTorsion()

  if (echoQ) write(*,"('makeNeighborQ()')")
  call MakeNeighborQ()
  if (echoQ) write(*,"('searchScalingQ()')")
  call MakeScalingQ()

  if (echoQ) write(*,"('makeAdjacentIndices()')")
  call MakeAdjacentIndices()
  if (echoQ) write(*,"('makeMolTree()')")
  call MakeMolTree()
  !            write(*,"('my turn')")
  if (echoQ) write(*,"('SetRelationShip()')")
  call SetRelationship()
  if (echoQ) write(*,"('ComToOrigin()')")
  !            write(*,"('end my turn')")
  call ComToOrigin
  !            call RescaleBond()

  call CheckMotions()

  call GroupingMolecule()

  close(UNIT_MOL)

  getMolInfQ = .true.

  end subroutine GetMolInf

  !********************

  subroutine CheckMotions()

  implicit none

  rotateQ = .true.
  bondQ   = .true.
  bendQ   = .true.
  torsQ   = .true.

  if (nAtoms == 1  ) rotateQ = .false.
  if (nBond  == 0  ) bondQ   = .false.
  if (nBend  == 0  ) bendQ   = .false.
  if (nTorsion == 0) torsQ   = .false.

  end subroutine CheckMotions

  subroutine PrintCheckedMotions(unit_)

  implicit none

  integer, optional, intent(in) :: unit_
  character(len=10), parameter :: DASHED = '----------'

  if (present(unit_)) then
    write(unit_,"(4(A10))") DASHED, DASHED, DASHED, DASHED
    write(unit_,"(4(A10))") 'rotateQ', 'bondQ', 'bendQ', 'torsQ'
    write(unit_,"(4(A10))") DASHED, DASHED, DASHED, DASHED
    write(unit_,"(4(L10))")  rotateQ, bondQ, bendQ, torsQ
  else
    write(*,"(4(A10))") DASHED, DASHED, DASHED, DASHED
    write(*,"(4(A10))") 'rotateQ', 'bondQ', 'bendQ', 'torsQ'
    write(*,"(4(A10))") DASHED, DASHED, DASHED, DASHED
    write(*,"(4(L10))")  rotateQ, bondQ, bendQ, torsQ
  end if

  end subroutine PrintCheckedMotions

  !********************

  subroutine ReadAtomType()

  implicit none

  integer :: i, dum

  select case (potenTypeENUM)

  case (OPLSAA)

    do i = 1, nType
      ! atomic num, eps, sig, charge, mass
      read(UNIT_MOL, *) dum, atomType(i,1), atomType(i,2), &
        atomType(i,3), atomType(i,4), atomType(i,5)
    end do

  case (FW)

    do i = 1, nType
      ! atomic num, mass
      read(UNIT_MOL, *) dum, atomType(i,1), atomType(i,5)
    end do

    ! matrix a
    do i = 1, nType
      read(UNIT_MOL, *) A_FW(i, 1:nType)
    end do


    ! matrix b
    do i = 1, nType
      read(UNIT_MOL, *) B_FW(i, 1:nType)
    end do


    ! matrix c
    do i = 1, nType
      read(UNIT_MOL, *) C_FW(i, 1:nType)
    end do

  case (DREIDING)

    do i = 1, nType

      read(UNIT_MOL, *) dum, &
        atomType(i,1), & ! atomic number
      zeta_DD(i,i),    & ! force field parameter 1
      D0_DD(i,i),      & ! force field parameter 2
      R0_DD(i,i),      & ! force field parameter 3
      charge(i),     & ! ESP charge model
      mass(i)          ! atomic weight

    end do

  case default

    stop 'not surpported potential in readAtomType'

  end select

  end subroutine ReadAtomType

  !********************

  subroutine ReadAtomIndex()

  implicit none

  integer :: i, j, dum, type_

  atomTypeNum = 0

  select case (potenTypeENUM)    ! in this moment, only OPLS-AA..

  case (OPLSAA)

    if (allocated(myType)) deallocate(myType, myAtomNum)
    allocate(myType(nAtoms), myAtomNum(nAtoms),    &
      myTypeInBend(nAtoms), myTypeInTorsion(nAtoms), &
      myTypeInBond(nAtoms))

    do i = 1, nAtoms

      read(UNIT_MOL, *) dum, myType(i),  myTypeInBond(i), &
        myTypeInBend(i), myTypeInTorsion(i)
      type_ = myType(i)
      myAtomNum(i) = nint(atomType(type_, 1))
      atomTypeNum(type_) = atomTypeNum(type_) + 1
      eps(type_)         = atomType(type_, 2)
      sig(type_)         = atomType(type_, 3)
      charge(type_)      = atomType(type_, 4)
      mass(type_)        = atomType(type_, 5)
      color(type_)       = myAtomNum(i)    ! will be changed... in the future ..color index is not the atomic number.

    end do

    do i = 1, nType        ! make piar interaction parameters.
      do j = 1, nType

        sigij(i,j) = sqrt(sig(i)*sig(j)) ! geometric mean for OPLS-AA size parameter
        epsij(i,j) = sqrt(eps(i)*eps(j)) ! geometric mean
        chaij(i,j) = charge(i)*charge(j) ! charge integeraction.

      end do ! i
    end do ! j

    !call MolCategorize() ! will be removed...


  case (FW)

    if (allocated(myType)) deallocate(myType, myAtomNum)
    allocate(myType(nAtoms), myAtomNum(nAtoms),    &
      myTypeInBend(nAtoms), myTypeInTorsion(nAtoms), &
      myTypeInBond(nAtoms))

    do i = 1, nAtoms

      read(UNIT_MOL, *) dum, myType(i),  myTypeInBond(i), &
        myTypeInBend(i), myTypeInTorsion(i)
      type_ = myType(i)
      myAtomNum(i) = nint(atomType(type_, 1))
      atomTypeNum(type_) = atomTypeNum(type_) + 1
      !eps(type_)         = atomType(type_, 2)
      !sig(type_)         = atomType(type_, 3)
      !charge(type_)      = atomType(type_, 4)
      mass(type_)        = atomType(type_, 5)
      color(type_)       = myAtomNum(i)    ! will be changed... in the future ..color index is not the atomic number.

    end do


  case (DREIDING)

    if (allocated(myType)) deallocate(myType, myAtomNum)
    allocate (                   &
      myType(nAtoms),          &
      myAtomNum(nAtoms),       &
      myTypeInBend(nAtoms),    &
      myTypeInTorsion(nAtoms), &
      myTypeInBond(nAtoms)     &
      )

    do i = 1, nAtoms

      read(UNIT_MOL, *) dum,               &
        myType(i),         &
        myTypeInBond(i),   &
        myTypeInBend(i),   &
        myTypeInTorsion(i)


      type_ = myType(i)
      myAtomNum(i) = nint(atomType(type_, 1))
      atomTypeNum(type_) = atomTypeNum(type_) + 1
      color(type_) = myAtomNum(i)

    end do

    ! make A_DD, B_DD, C_DD

    do i = 1, nType

      A_DD(i,i) = D0_DD(i,i) * 6.0_r8 / (zeta_DD(i,i) - 6.0_r8) * exp(zeta_DD(i,i))
      B_DD(i,i) = D0_DD(i,i) * zeta_DD(i,i) / (zeta_DD(i,i) - 6.0_r8) * (R0_DD(i,i) ** 6)
      C_DD(i,i) = zeta_DD(i,i) / R0_DD(i,i)
      chaij(i,i) = charge(i) * charge(i) ! charge integeraction.

    end do

    do i = 1, nType - 1

      do j = i + 1, nType

        A_DD(i,j) = sqrt(A_DD(i,i) * A_DD(j,j))      ! geometric mean
        B_DD(i,j) = sqrt(B_DD(i,i) * B_DD(j,j))      ! geometric mean
        C_DD(i,j) = 0.5_r8 * (C_DD(i,i) + C_DD(j,j)) ! arithmetic mean

        ! symmetric condition.
        A_DD(j,i) = A_DD(i,j)
        B_DD(j,i) = B_DD(i,j)
        C_DD(j,i) = C_DD(i,j)

        ! charge interaction.
        chaij(i,j) = charge(i) * charge(j) ! charge integeraction.
        chaij(j,i) = chaij(i,j)

        zeta_DD(i,j) = sqrt(zeta_DD(i,i) * zeta_DD(j,j))
        zeta_DD(j,i) = zeta_DD(i,j)
        R0_DD(i,j) = sqrt(R0_DD(i,i) * R0_DD(j,j))
        R0_DD(j,i) = R0_DD(i,j)

      end do

    end do

    MIN_DISTANCE_DD = 1.2_r8 * (-0.0292458_r8 - 1.03031 / zeta_DD + 60.0109 / zeta_DD / zeta_DD ) * R0_DD ! parallel calculation.

    do i = 1, nType

      do j = 1, nType

        !write(*,"('min distance = ', 2I3, 2F10.4)") i, j, zeta_DD(i,j), MIN_DISTANCE_DD(i,j)
        !write(1,"('min distance = ', F10.4)") MIN_DISTANCE_DD(i,j)

      end do

    end do

    case default

    stop 'not surpported potential in readAtomIndex'

  end select

  end subroutine ReadAtomIndex

  !********************

  subroutine ReadBond()

  implicit none

  integer :: i, temp1, temp2

  if (allocated(bondIndices)) deallocate(bondIndices) ! for avoid runtime error.
  allocate(bondIndices(nBond,2))

  do i = 1, nBond
    read(UNIT_MOL, *) bondIndices(i,1), bondIndices(i,2)
    temp1 = bondIndices(i,1)
    temp2 = bondIndices(i,2)
    if (temp1 > temp2) then ! sort by index magnitude
      bondIndices(i,1) = temp2
      bondIndices(i,2) = temp1
    end if
  end do

  end subroutine ReadBond

  !********************
  !********************

  subroutine ReadBondType()

  implicit none

  integer :: i, dum, temp1, temp2

  do i = 1, nBondType

    read(UNIT_MOL, *) dum, bondType(i,1), bondType(i,2)
    temp1 = bendType(i,1)
    temp2 = bendType(i,2)

    if (temp1 > temp2) then ! sort by type num magnitude
      bendType(i,1) = temp2
      bendType(i,2) = temp1
    end if

  end do

  end subroutine ReadBondType

  !********************

  subroutine ReadBondParameter()

  implicit none

  integer :: i, dum
  select case (potenTypeENUM)

  case (OPLSAA, FW, DREIDING)

    do i = 1, nBondType

      read(UNIT_MOL, *) dum, kBond(i), bondLength(i) ! bond length indicates equilibrium bond length
      ! kBond [K] unit, bondLength [Å] unit
    end do

    case default

    stop 'not surpported potential in ReadBondParameter'

  end select

  end subroutine ReadBondParameter

  !********************

  subroutine SearchBond()

  implicit none

  integer :: bond, b1, b2

  if (allocated(myBondType)) deallocate(myBondType)
  allocate(myBondType(nBond))

  do bond = 1, nBond

    b1 = myTypeInBond(bondIndices(bond,1))
    b2 = myTypeInBond(bondIndices(bond,2))

    myBondType(bond) = GiveBondType(b1, b2)

  end do

  end subroutine SearchBond

  !********************

  function GiveBondType(b1, b2) result(bondType_)

  implicit none

  logical :: findBondTypeQ
  integer :: i
  integer :: temp1, temp2
  integer, intent(in) :: b1, b2
  integer :: b1_, b2_ ! copied variable for sorting
  integer :: bondType_

  findBondTypeQ = .false.

  b1_ = b1
  b2_ = b2

  if (b1_ > b2_) then ! sort for searching
    b1_ = b2
    b2_ = b1
  end if

  do i = 1, nBondType
    temp1 = bondType(i,1)  ! bondtype already sorted.
    temp2 = bondType(i,2)

    if (temp1 /= b1_) cycle
    if (temp2 /= b2_) cycle

    bondType_ = i
    findBondTypeQ = .true.
    exit ! from loop
  end do

  if (.not.findBondTypeQ) stop 'bond information is needed'

  end function GiveBondType

  !********************

  subroutine PrintBondInformation(unit_)

  implicit none

  integer, optional, intent(in) :: unit_

  integer :: i

  if (present(unit_)) then
    write(unit_, "(A, I5)") 'number of bond in molecule = ', nBond
    do i = 1, nBond
      write(unit_,"('[bond, type ',2(I3),']', 2(I5), ', length = ', F8.3)") i, myBondType(i), &
        bondIndices(i,1), bondIndices(i,2), bondLength(myBondType(i))
    end do
  else
    write(*, "(A, I5)") 'number of bond in molecule = ', nBond
    do i = 1, nBond
      write(*,"('[bond, type ',2(I3),']', 2(I5), ', length = ', F8.3)") i, myBondType(i), &
        bondIndices(i,1), bondIndices(i,2), bondLength(myBondType(i))
    end do
  end if

  end subroutine PrintBondInformation

  !********************
  subroutine ReadBendType()

  implicit none

  integer :: i, dum, temp1, temp2

  do i = 1, nBendType
    read(UNIT_MOL, *) dum, bendType(i,1), bendType(i,2), bendType(i,3)
    temp1 = bendType(i,1)
    temp2 = bendType(i,3)

    if (temp1 > temp2) then    ! atom types in bend are automatically sorted by magnitude.
      bendType(i,1) = temp2
      bendType(i,3) = temp1
    end if
  end do

  end subroutine ReadBendType

  !********************

  subroutine ReadBendParameter()

  use mod_univ_const, only : DEG_TO_RAD

  implicit none

  integer :: i, dum

  select case (potenTypeENUM)

  case (OPLSAA, FW)

    do i = 1, nBendType
      read(UNIT_MOL, *) dum, kBend(i), theta0Deg(i)
    end do

    theta0 = theta0Deg * DEG_TO_RAD  ! vector operation.

  case (DREIDING)

    do i = 1, nBendType

      read(UNIT_MOL, *) dum, kBend(i), theta0Deg(i)

      theta0(i)     = theta0Deg(i) * DEG_TO_RAD
      cos_theta0(i) = cos(theta0(i))
      C_theta_DD(i) = 0.5_r8 * kBend(i) / (sin(theta0(i)) ** 2)

    end do

    case default

    stop 'not surpported potential in readBondParameter'

  end select

  end subroutine ReadBendParameter

  !********************

  subroutine ReadTorsionType()

  implicit none

  integer :: i, dum, temp1, temp2, temp3, temp4

  do i = 1, nTorsionType
    read(UNIT_MOL, *) dum, torsionType(i,1), torsionType(i,2), &
      &                              torsionType(i,3), torsionType(i,4)
    temp1 = torsionType(i,1)
    temp2 = torsionType(i,2)
    temp3 = torsionType(i,3)
    temp4 = torsionType(i,4)

    if (temp1 > temp4) then    ! atom types in torsion are automatically sorted by 1-4 atoms magnitude.
      torsionType(i,1) = temp4
      torsionType(i,2) = temp3
      torsionType(i,3) = temp2
      torsionType(i,4) = temp1
    end if
  end do

  end subroutine ReadTorsionType

  !********************

  subroutine ReadTorsionParameter()

  implicit none

  integer :: i, dum

  select case (potenTypeENUM)

  case (OPLSAA, FW)

    do i = 1, nTorsionType

      read(UNIT_MOL, *) dum, v1(i), v2(i), v3(i), v4(i)

    end do

  case (DREIDING)

    do i = 1, nTorsionType

      read(UNIT_MOL, *) dum, v1(i), v2(i), v3(i), v4(i), v5(i), v6(i)
      read(UNIT_MOL, *) dum, d1(i), d2(i), d3(i), d4(i), d5(i), d6(i)

    end do

    case default

    stop 'not surpported potential in ReadTorsionParameter'

  end select

  end subroutine ReadTorsionParameter

  !********************

  subroutine ReadMoleculeCoordinate()

  implicit none

  integer :: i, dum

  if (allocated(molCrdX)) deallocate(molCrdX, molCrdY, molCrdZ)
  allocate(molCrdX(nAtoms), molCrdY(nAtoms), molCrdZ(nAtoms))

  do i = 1, nAtoms
    read(UNIT_MOL, *) dum, molCrdX(i), molCrdY(i), molCrdZ(i)
  end do

  end subroutine ReadMoleculeCoordinate
  !********************

  subroutine ComToOrigin()   ! change COM coordinate of molecule to origin

  real(kind=r8) :: comX, comY, comZ

  molMass = sum(mass(1:nType)*atomTypeNum(1:nType))

  comX    = sum(mass(myType(:))*molCrdX) / molMass
  comY    = sum(mass(myType(:))*molCrdY) / molMass
  comZ    = sum(mass(myType(:))*molCrdZ) / molMass

  molCrdX = molCrdX - comX
  molCrdY = molCrdY - comY
  molCrdZ = molCrdZ - comZ

  end subroutine ComToOrigin

  !*******************

  subroutine SearchBending()

  implicit none

  integer :: i, j, direction1, direction2, otherDirection2
  integer :: b1, b2, b3
  integer, dimension(-1:1) :: indexi   ! indexi(-1) or indexi(1)

  ! ===== count # of bending =====

  nBend = 0  ! initialize

  do i = 1, nBond - 1

    indexi(-1) = bondIndices(i,1) ! left direction
    indexi(1)  = bondIndices(i,2) ! right direction

    do j = i + 1, nBond

      do direction1 = -1, 1, 2 ! direction1 = -1 or 1

        do direction2 = 1, 2

          if (indexi(direction1) == bondIndices(j, direction2)) then ! become pivot

            nBend = nBend + 1 !  only count # of bending for array allocation.

          end if

        end do ! direction2 loop

      end do ! direction1 loop

    end do ! j loop

  end do ! i loop

  ! ===== end of counting =====

  if (allocated(bendIndices)) deallocate(bendIndices, myBendType)

  allocate(bendIndices(nBend,3), myBendType(nBend))

  nBend = 0  ! becomes zero to use as array index

  do i = 1, nBond - 1

    indexi(-1) = bondIndices(i,1) ! left direction
    indexi(1)  = bondIndices(i,2) ! right direction

    do j = i + 1, nBond          ! search other bond j which is well matched bond i

      do direction1 = -1, 1, 2 ! direction1 = -1 or 1

        do direction2 = 1, 2

          if (indexi(direction1) == bondIndices(j, direction2)) then ! become pivot

            nBend = nBend + 1 !  count # of bending for array index.
            bendIndices(nBend,1) = indexi(-direction1)
            bendIndices(nBend,2) = indexi(direction1)
            if (direction2 == 1) otherDirection2 = 2
            if (direction2 == 2) otherDirection2 = 1
            bendIndices(nBend,3) = bondIndices(j, otherDirection2)

            b1 = myTypeInBend(bendIndices(nBend,1)) ! my atom type
            b2 = myTypeInBend(bendIndices(nBend,2))
            b3 = myTypeInBend(bendIndices(nBend,3)) ! bend indices not sorted.

            myBendType(nBend) = GiveBendType(b1, b2, b3)

          end if

        end do ! direction2 loop

      end do ! direction1 loop

    end do ! j loop

  end do ! i loop

  end subroutine SearchBending

  !********************

  function GiveBendType(b1, b2, b3) result(bendType_)

  implicit none

  logical :: findBendTypeQ
  integer :: i
  integer :: temp1, temp2, temp3
  integer, intent(in) :: b1, b2 ,b3
  integer :: b1_, b2_, b3_ ! copied variable for sorting
  integer :: bendType_

  findBendTypeQ = .false.

  b1_ = b1
  b2_ = b2
  b3_ = b3

  if (b1_ > b3_) then ! sort for searching
    b1_ = b3
    b3_ = b1
  end if

  do i = 1, nBendType
    temp1 = bendType(i,1)
    temp2 = bendType(i,2)
    temp3 = bendType(i,3)

    if (temp1 /= b1_) cycle
    if (temp2 /= b2_) cycle
    if (temp3 /= b3_) cycle

    bendType_ = i
    findBendTypeQ = .true.
    exit ! from loop
  end do

  if (.not.findBendTypeQ) stop 'bend information is needed'

  end function GiveBendType

  !********************

  subroutine PrintBendInformation(unit_)

  implicit none

  integer, optional, intent(in) :: unit_

  integer :: i

  if (present(unit_)) then
    write(unit_, "(A, I5)") 'number of beding in molecule = ', nBend
    do i = 1, nBend
      write(unit_,"('[bend, type ',2(I3),']', 3(I5))") i, myBendType(i), &
        bendIndices(i,1), bendIndices(i,2), bendIndices(i,3)
    end do
  else
    write(*, "(A, I5)") 'number of beding in molecule = ', nBend
    do i = 1, nBend
      write(*,"('[bend, type ',2(I3),']', 3(I5))") i, myBendType(i), &
        bendIndices(i,1), bendIndices(i,2), bendIndices(i,3)
    end do
  end if

  end subroutine PrintBendInformation

  !*******************

  subroutine SearchTorsion()

  implicit none

  integer :: i, j, direction1, direction2
  integer :: t1, t2, t3, t4
  integer :: ipivot1, ipivot2, jpivot1, jpivot2

  ! ===== count # of torsion =====

  nTorsion = 0  ! initialize

  do i = 1, nBend - 1  ! for all bending group pair

    do direction1 = -1, 1, 2

      ipivot1 = bendIndices(i,2)
      ipivot2 = bendIndices(i,2+direction1)

      do j = i + 1, nBend

        do direction2 = -1, 1, 2

          jpivot1 = bendIndices(j,2+direction2)
          jpivot2 = bendIndices(j,2)

          if (ipivot1 /= jpivot1) cycle
          if (ipivot2 /= jpivot2) cycle

          ! if two beding group share torsion-axes

          nTorsion = nTorsion + 1

        end do ! direction2

      end do ! j-loop

    end do ! direction1

  end do ! i-loop

  ! ===== end of counting =====

  if (allocated(torsionIndices)) deallocate(torsionIndices, myTorsionType)
  allocate(torsionIndices(nTorsion,4), myTorsionType(nTorsion))

  nTorsion = 0  ! becomes zero to use as array index

  do i = 1, nBend - 1  ! for all bending group pair

    do direction1 = -1, 1, 2

      ipivot1 = bendIndices(i,2)
      ipivot2 = bendIndices(i,2+direction1)

      do j = i + 1, nBend

        do direction2 = -1, 1, 2

          jpivot1 = bendIndices(j,2+direction2)
          jpivot2 = bendIndices(j,2)

          if (ipivot1 /= jpivot1) cycle
          if (ipivot2 /= jpivot2) cycle

          ! if two beding group share torsion-axes

          nTorsion = nTorsion + 1
          torsionIndices(nTorsion,1) = bendIndices(i,2-direction1)
          torsionIndices(nTorsion,2) = bendIndices(i,2)
          torsionIndices(nTorsion,3) = bendIndices(j,2)
          torsionIndices(nTorsion,4) = bendIndices(j,2-direction2)

          t1 = myTypeInTorsion(torsionIndices(nTorsion,1))
          t2 = myTypeInTorsion(torsionIndices(nTorsion,2))
          t3 = myTypeInTorsion(torsionIndices(nTorsion,3))
          t4 = myTypeInTorsion(torsionIndices(nTorsion,4))

          myTorsionType(nTorsion) = GiveTorsionType(t1, t2, t3, t4)

        end do ! direction2

      end do ! j-loop

    end do ! direction1

  end do ! i-loop

  end subroutine SearchTorsion

  !******************************

  function GiveTorsionType(t1, t2, t3, t4) result(torsionType_)

  implicit none

  logical :: findTorsionTypeQ
  integer :: i
  integer :: temp1, temp2, temp3, temp4
  integer, intent(in) :: t1, t2, t3, t4
  integer :: t1_, t2_, t3_, t4_ ! copied variable for sorting
  integer :: torsionType_

  findTorsionTypeQ = .false.

  t1_ = t1
  t2_ = t2
  t3_ = t3
  t4_ = t4

  if (t1_ > t4_) then

    t1_ = t4
    t2_ = t3
    t3_ = t2
    t4_ = t1

  end if

  do i = 1, nTorsionType
    temp1 = torsionType(i,1)
    temp2 = torsionType(i,2)
    temp3 = torsionType(i,3)
    temp4 = torsionType(i,4)

    if (temp1 /= t1_) cycle
    if (temp2 /= t2_) cycle
    if (temp3 /= t3_) cycle
    if (temp4 /= t4_) cycle

    torsionType_ = i
    findTorsionTypeQ = .true.

    exit ! from loop
  end do

  if (.not.findTorsionTypeQ) stop 'torsion information is needed'

  end function GiveTorsionType

  !********************

  subroutine PrintTorsionInformation(unit_)

  implicit none

  integer, optional, intent(in) :: unit_

  integer :: i

  if (present(unit_)) then
    do i = 1, nTorsion
      write(unit_,"('[torsion, type ',2(I3),']', 4(I5))") i, myTorsionType(i), &
        torsionIndices(i,1), torsionIndices(i,2), torsionIndices(i,3), torsionIndices(i,4)
    end do
  else
    do i = 1, nTorsion
      write(*,"('[torsion, type ',2(I3),']', 4(I5))") i, myTorsionType(i), &
        torsionIndices(i,1), torsionIndices(i,2), torsionIndices(i,3), torsionIndices(i,4)
    end do
  end if

  end subroutine PrintTorsionInformation

  !********************

  subroutine MakeNeighborQ()

  implicit none

  integer :: i, j, k ,n

  if (allocated(neighborQ)) deallocate(neighborQ)
  allocate(neighborQ(nAtoms, nAtoms))

  ! initialize neighborQ to false
  do i = 1, nAtoms
    do j = 1, nAtoms
      neighborQ(i,j) = .false.
    end do
  end do

  do n = 1, nBend

    i = bendIndices(n,1)
    j = bendIndices(n,2)
    k = bendIndices(n,3)

    neighborQ(i,j) = .true.
    neighborQ(j,i) = .true.
    neighborQ(j,k) = .true.
    neighborQ(k,j) = .true.
    neighborQ(k,i) = .true.
    neighborQ(i,k) = .true.

  end do

  end subroutine MakeNeighborQ

  !********************

  subroutine PrintNeighborQ(unit_)

  implicit none

  integer, optional, intent(in) :: unit_

  integer :: i

  if (present(unit_)) then
    write(unit_, "(A)") 'neighbor matrix = '
    do i = 1, nAtoms
      write(unit_, *) neighborQ(i,:)
    end do

  else
    write(*, "(A)") 'neighbor matrix = '
    do i = 1, nAtoms
      write(*, *) neighborQ(i,:)
    end do

  end if

  end subroutine PrintNeighborQ

  !********************

  subroutine MakeScalingQ()

  implicit none

  integer :: i, j, n, end1, end2

  if (allocated(scalingQ)) deallocate(scalingQ)
  allocate(scalingQ(nAtoms, nAtoms))

  do i = 1, nAtoms
    do j = 1, nAtoms
      scalingQ(i,j) = .false.
    end do
  end do

  do n = 1, nTorsion

    end1 = torsionIndices(n,1)
    end2 = torsionIndices(n,4)

    scalingQ(end1,end2) = .true.
    scalingQ(end2,end1) = .true.

  end do

  end subroutine MakeScalingQ

  !********************

  subroutine PrintScalingQ(unit_)

  implicit none

  integer, optional, intent(in) :: unit_

  integer :: i

  if (present(unit_)) then
    write(unit_, "(A)") 'scaling matrix = '
    do i = 1, nAtoms
      write(unit_, *) scalingQ(i,:)
    end do

  else
    write(*, "(A)") 'scaling matrix = '
    do i = 1, nAtoms
      write(*, *) scalingQ(i,:)
    end do

  end if

  end subroutine PrintScalingQ

  !********************

  subroutine MakeAdjacentIndices()

  implicit none

  integer :: a, b, i, index, nAdjac, otherIndex

  if (allocated(adjacentIndices)) deallocate(adjacentIndices)
  allocate(adjacentIndices(nAtoms,0:4))

  do a = 1, nAtoms
    adjacentIndices(a,0) = 0
  end do

  do a = 1, nAtoms
    do b = 1, nBond
      do i = 1, 2

        index = bondIndices(b,i)
        if (a /= index) cycle
        ! if a = index
        adjacentIndices(a,0) = adjacentIndices(a,0) + 1
        nAdjac = adjacentIndices(a,0)
        if (i == 1) otherIndex = bondIndices(b,2)
        if (i == 2) otherIndex = bondIndices(b,1)
        adjacentIndices(a,nAdjac) = otherIndex

      end do
    end do ! b-loop (bond)
  end do ! a-loop (atom)

  end subroutine MakeAdjacentIndices

  !********************
  subroutine PrintAdjacentIndices(unit_)

  implicit none

  integer, optional, intent(in) :: unit_

  integer :: i, nAdjac

  if (present(unit_)) then
    write(unit_, "(A)") 'adjacent indices list'
    write(unit_, "(A6, A7, A8)") 'index', 'nAdjac', 'indices'
    do i = 1, nAtoms
      nAdjac = adjacentIndices(i,0)
      write(unit_,"(I6, I7, I8, 5(I3))") i, nAdjac, adjacentIndices(i,1:nAdjac)
    end do
  else
    write(*, "(A)") 'adjacent indices list'
    write(*, "(A6, A7, A8)") 'index', 'nAdjac', 'indices'
    do i = 1, nAtoms
      nAdjac = adjacentIndices(i,0)
      write(*,"(I6, I7, I8, 5(I3))") i, nAdjac, adjacentIndices(i,1:nAdjac)
    end do
  end if

  end subroutine PrintAdjacentIndices

  !********************

  subroutine MakeMolTree()

  implicit none

  logical :: forwardQ
  integer :: i, a, path, nIndices, startIndex, currentIndex, nextIndex, beforeIndex
  integer :: nForwardLevel, currentRestPath

  ! ===== allocatable variables =====
  integer, dimension(:), allocatable :: restPath, savedBeforeIndex

  if (allocated(molTree)) deallocate(molTree)
  allocate(molTree(nAtoms, 4, 0:nAtoms-1))
  allocate(restPath(nAtoms), savedBeforeIndex(0:nAtoms))

  do i = 1, nAtoms

    do a = 1, nAtoms
      restPath(a) = adjacentIndices(a,0)
    end do

    do path = 1, adjacentIndices(i,0)

      startIndex = adjacentIndices(i,path)
      currentIndex = startIndex

      nIndices = 1 ! because currentIndex.. at least has 1 index
      molTree(i,path,1) = startIndex

      forwardQ = .true.
      nForwardLevel = 0
      savedBeforeIndex(nForwardLevel) = i ! never return this index

      do while(.true.)  ! search all index contained this path

        if (forwardQ) then

          currentRestPath = restPath(currentIndex)
          if (currentRestPath == 0) then

            forwardQ = .false.
            cycle

          end if

          nextIndex   = adjacentIndices(currentIndex, currentRestPath)
          beforeIndex = savedBeforeIndex(nForwardLevel)

          if (nextIndex == beforeIndex) then ! self foling

            restPath(currentIndex) = restPath(currentIndex) - 1 ! give up path

          else

            nIndices = nIndices + 1
            molTree(i,path,nIndices) = nextIndex

            nForwardLevel = nForwardLevel + 1
            savedBeforeIndex(nForwardLevel) = currentIndex
            currentIndex = nextIndex

          end if

        else ! backword or stop loop

          if (startIndex == currentIndex .and. restPath(currentIndex) == 0) then
            molTree(i,path,0) = nIndices
            exit ! go to next path in i
          end if

          ! backword moving
          currentIndex = savedBeforeIndex(nForwardLevel)
          restPath(currentIndex) = restPath(currentIndex) - 1
          nForwardLevel = nForwardLevel - 1 ! reduce forward-level because of backword moving.
          if (restPath(currentIndex) > 0) forwardQ = .true. ! if there is anywhere to go.

        end if

      end do ! do while - loop

    end do ! path-loop

  end do ! i-loop: atoms

  deallocate(restPath, savedBeforeIndex)

  end subroutine MakeMolTree

  !********************

  subroutine PrintMolTree(unit_)

  implicit none

  integer, optional, intent(in) :: unit_
  integer :: a, i, nIndices, path

  if (present(unit_)) then
    write(unit_,"(A)") 'molecular tree data = '
    write(unit_,"(10(A5))") 'i', 'path', 'index'
    do a = 1, nAtoms
      do path = 1, 4
        nIndices = molTree(a,path,0)
        if (nIndices == 0) cycle
        write(unit_,"('[',2(I3),']',30(I3))") a, path, molTree(a,path,1:nIndices)
      end do ! path-loop
    end do ! a-loop: atoms
  else
    write(*,"(A)") 'molecular tree data = '
    write(*,"(10(A5))") 'i', 'path', 'index'
    do a = 1, nAtoms
      do path = 1, 4
        nIndices = molTree(a,path,0)
        if (nIndices == 0) cycle
        write(*,"('[',2(I3),']',30(I3))") a, path, molTree(a,path,1:nIndices)
      end do ! path-loop
    end do ! a-loop: atoms
  end if

  end subroutine PrintMolTree

  !********************

  subroutine SetRelationship()

  implicit none

  integer :: bond, bend, torsion, idx1, idx2, idx3, idx4

  !write(*,"(A)") '0'
  ! ====== allocate array used in energy calculation and move =============

  allocate(index1InBond(nBond), index2InBond(nBond))
  allocate(index1InBend(nBend), index2InBend(nBend), index3InBend(nBend))
  allocate(index1InTorsion(nTorsion), index2InTorsion(nTorsion), &
    &        index3InTorsion(nTorsion), index4InTorsion(nTorsion))
  allocate(bond12InBend(nBend), bond23InBend(nBend), &
    &       direc12InBend(nBend), direc23InBend(nBend))
  allocate(bond12InTorsion(nTorsion), bond23InTorsion(nTorsion), &
    &       bond34InTorsion(nTorsion), direc12InTorsion(nTorsion), &
    &      direc23InTorsion(nTorsion), direc34InTorsion(nTorsion), &
    &      bend123InTorsion(nTorsion), bend234InTorsion(nTorsion), &
    &      direc123InTorsion(nTorsion), direc234InTorsion(nTorsion))

  ! =======================================================================

  ! ==== set index ====

  do bond = 1, nBond
    index1InBond(bond) = bondIndices(bond,1)
    index2InBond(bond) = bondIndices(bond,2)
  end do
  !write(*,"(A)") '1'
  do bend = 1, nBend
    index1InBend(bend) = bendIndices(bend,1)
    index2InBend(bend) = bendIndices(bend,2)
    index3InBend(bend) = bendIndices(bend,3)
  end do
  !write(*,"(A)") '2'
  do torsion = 1, nTorsion
    index1InTorsion(torsion) = torsionIndices(torsion,1)
    index2InTorsion(torsion) = torsionIndices(torsion,2)
    index3InTorsion(torsion) = torsionIndices(torsion,3)
    index4InTorsion(torsion) = torsionIndices(torsion,4)
  end do
  !write(*,"(A)") '3'
  ! ===================

  ! ===== find bond in bend =====

  do bend = 1, nBend

    idx1 = index1InBend(bend)
    idx2 = index2InBend(bend)
    idx3 = index3InBend(bend)

    bond12InBend(bend) = 0
    do bond = 1, nBond

      if (idx1 == index1InBond(bond) .and. idx2 == index2InBond(bond)) then
        bond12InBend(bend)  = bond
        direc12InBend(bend) = 1
      else if (idx1 == index2InBond(bond) .and. idx2 == index1InBond(bond)) then
        bond12InBend(bend)  = bond
        direc12InBend(bend) = -1
      end if

    end do ! bond-loop
    if (bond12InBend(bend) == 0) stop 'relationship not found'


    bond23InBend(bend) = 0
    do bond = 1, nBond

      if (idx2 == index1InBond(bond) .and. idx3 == index2InBond(bond)) then
        bond23InBend(bend)  = bond
        direc23InBend(bend) = 1
      else if (idx2 == index2InBond(bond) .and. idx3 == index1InBond(bond)) then
        bond23InBend(bend)  = bond
        direc23InBend(bend) = -1
      end if

    end do ! bond-loop
    if (bond23InBend(bend) == 0) stop 'relationship not found'

  end do ! bend-loop
  !write(*,"(A)") '4'
  ! ===== find bond in torsion =====

  do torsion = 1, nTorsion

    idx1 = index1InTorsion(torsion)
    idx2 = index2InTorsion(torsion)
    idx3 = index3InTorsion(torsion)
    idx4 = index4InTorsion(torsion)

    bond12InTorsion(torsion) = 0
    do bond = 1, nBond

      if (idx1 == index1InBond(bond) .and. idx2 == index2InBond(bond)) then
        bond12InTorsion(torsion)  = bond
        direc12InTorsion(torsion) = 1
      else if (idx1 == index2InBond(bond) .and. idx2 == index1InBond(bond)) then
        bond12InTorsion(torsion)  = bond
        direc12InTorsion(torsion) = -1
      end if

    end do ! bond-loop
    if (bond12InTorsion(torsion) == 0) stop 'relationship not found'

    bond23InTorsion(torsion) = 0
    do bond = 1, nBond

      if (idx2 == index1InBond(bond) .and. idx3 == index2InBond(bond)) then
        bond23InTorsion(torsion)  = bond
        direc23InTorsion(torsion) = 1
      else if (idx2 == index2InBond(bond) .and. idx3 == index1InBond(bond)) then
        bond23InTorsion(torsion)  = bond
        direc23InTorsion(torsion) = -1
      end if

    end do ! bond-loop
    if (bond23InTorsion(torsion) == 0) stop 'relationship not found'

    bond34InTorsion(torsion) = 0
    do bond = 1, nBond

      if (idx3 == index1InBond(bond) .and. idx4 == index2InBond(bond)) then
        bond34InTorsion(torsion)  = bond
        direc34InTorsion(torsion) = 1
      else if (idx3 == index2InBond(bond) .and. idx4 == index1InBond(bond)) then
        bond34InTorsion(torsion)  = bond
        direc34InTorsion(torsion) = -1
      end if

    end do ! bond-loop
    if (bond34InTorsion(torsion) == 0) stop 'relationship not found'

  end do ! torsion-loop
  !write(*,"(A)") '5'
  ! ===== find bend in torsion =====

  do torsion = 1, nTorsion

    idx1 = index1InTorsion(torsion)
    idx2 = index2InTorsion(torsion)
    idx3 = index3InTorsion(torsion)
    idx4 = index4InTorsion(torsion)

    bend123InTorsion(torsion) = 0
    do bend = 1, nBend

      if (idx1 == index1InBend(bend) .and. idx2 == index2InBend(bend) .and. &
        &       idx3 == index3InBend(bend)) then

      bend123InTorsion(torsion) = bend
      direc123InTorsion(torsion) = 1

      else if (idx1 == index3InBend(bend) .and. idx2 == index2InBend(bend) .and. &
        &            idx3 == index1InBend(bend)) then

      bend123InTorsion(torsion) = bend
      direc123InTorsion(torsion) = -1

      end if

    end do ! bend-loop
    if (bend123InTorsion(torsion) == 0) stop 'relationship not found'

    bend234InTorsion(torsion) = 0
    do bend = 1, nBend

      if (idx2 == index1InBend(bend) .and. idx3 == index2InBend(bend) .and. &
        &       idx4 == index3InBend(bend)) then

      bend234InTorsion(torsion) = bend
      direc234InTorsion(torsion) = 1

      else if (idx2 == index3InBend(bend) .and. idx3 == index2InBend(bend) .and. &
        &            idx4 == index1InBend(bend)) then

      bend234InTorsion(torsion) = bend
      direc234InTorsion(torsion) = -1

      end if

    end do ! bend-loop
    if (bend234InTorsion(torsion) == 0) stop 'relationship not found'

  end do ! torsion-loop
  !write(*,"(A)") '6'
  end subroutine SetRelationship

  !********************

  subroutine PrintRelationship(unit_)

  implicit none

  integer, optional, intent(in) :: unit_

  integer :: i

  if (present(unit_)) then

    write(unit_,"(A)") 'bond'
    write(unit_,"(4(A4))") 'i', 'i1', 'i2'
    do i = 1, nBond
      write(unit_,"(4(I4))") i, index1InBond(i), index2InBond(i)
    end do

    write(unit_,"(A)") 'bend'
    write(unit_,"(/4(A4))") 'i', 'i1', 'i2', 'i3'
    do i = 1, nBend
      write(unit_,"(4(I4))") i, index1InBend(i), index2InBend(i), index3InBend(i)
      write(unit_,"(4(I4))") bond12InBend(i), bond23InBend(i)
      write(unit_,"(4(I4))") direc12InBend(i), direc23InBend(i)
    end do

    write(unit_,"(A)") 'torsion'
    write(unit_,"(/5(A4))") 'i', 'i1', 'i2', 'i3', 'i4'
    do i = 1, nTorsion
      write(unit_,"(/5(I4))") i, index1InTorsion(i), index2InTorsion(i), &
        &                          index3InTorsion(i), index4InTorsion(i)
      write(unit_,"(4(I4))") i, bond12InTorsion(i), bond23InTorsion(i), bond34InTorsion(i)
      write(unit_,"(4(I4))") i, direc12InTorsion(i), direc23InTorsion(i), direc34InTorsion(i)
    end do

    write(unit_,"(/5(A4))") 'i', 'i1', 'i2', 'i3', 'i4'
    do i = 1, nTorsion
      write(unit_,"(/5(I4))") i, index1InTorsion(i), index2InTorsion(i), &
        &                          index3InTorsion(i), index4InTorsion(i)
      write(unit_,"(4(I4))") i, bend123InTorsion(i), bend234InTorsion(i)
      write(unit_,"(4(I4))") i, direc123InTorsion(i), direc234InTorsion(i)
    end do

  else

    write(*,"(A)") 'bond'
    write(*,"(4(A4))") 'i', 'i1', 'i2'
    do i = 1, nBond
      write(*,"(4(I4))") i, index1InBond(i), index2InBond(i)
    end do

    write(*,"(A)") 'bend'
    write(*,"(/4(A4))") 'i', 'i1', 'i2', 'i3'
    do i = 1, nBend
      write(*,"(4(I4))") i, index1InBend(i), index2InBend(i), index3InBend(i)
      write(*,"(4(I4))") bond12InBend(i), bond23InBend(i)
      write(*,"(4(I4))") direc12InBend(i), direc23InBend(i)
    end do

    write(*,"(A)") 'torsion'
    write(*,"(/5(A4))") 'i', 'i1', 'i2', 'i3', 'i4'
    do i = 1, nTorsion
      write(*,"(/5(I4))") i, index1InTorsion(i), index2InTorsion(i), &
        &                          index3InTorsion(i), index4InTorsion(i)
      write(*,"(4(I4))") i, bond12InTorsion(i), bond23InTorsion(i), bond34InTorsion(i)
      write(*,"(4(I4))") i, direc12InTorsion(i), direc23InTorsion(i), direc34InTorsion(i)
    end do

    write(*,"(/5(A4))") 'i', 'i1', 'i2', 'i3', 'i4'
    do i = 1, nTorsion
      write(*,"(/5(I4))") i, index1InTorsion(i), index2InTorsion(i), &
        &                          index3InTorsion(i), index4InTorsion(i)
      write(*,"(4(I4))") i, bend123InTorsion(i), bend234InTorsion(i)
      write(*,"(4(I4))") i, direc123InTorsion(i), direc234InTorsion(i)
    end do

  end if

  end subroutine PrintRelationship

  !********************

  subroutine MolCategorize  ! not used

  implicit none

  logical :: findHomeQ = .false. ! 자신의 카테고리를 찾으면 .true.
  integer i, j

  maxCate = 1

  molCate(1,1) = eps(1)
  molCate(1,2) = sig(1)
  molCate(1,3) = charge(1)
  molCate(1,4) = 1.0_r8

  do i = 2, nAtoms; do j = 1, maxCate

    if (abs(molCate(j,1)-eps(i))    < 10e-5) then
      if (abs(molCate(j,2)-sig(i))    < 10e-5) then
        if (abs(molCate(j,3)-charge(i)) < 10e-5) then

          molCate(j,4) = molCate(j,4) + 1.0_r8
          findHomeQ = .true.

        end if; end if; end if;

      end do;

      if (.not.findHomeQ) then

        maxCate = maxCate + 1
        molCate(maxCate,1) = eps(i)
        molCate(maxCate,2) = sig(i)
        molCate(maxCate,3) = charge(i)
        molCate(maxCate,4) = 1.0_r8

      end if

      findHomeQ = .false.; end do

    end subroutine MolCategorize

    !********************

    subroutine PrintMolInf (unit_) ! print molecule properties

    integer :: i, ai
    integer, optional, intent(in) :: unit_

    if (.not.getMolInfQ) stop '입력받은 분자정보가 없습니다.'

    if (present(unit_)) then

      write(unit_,"(A/)")        'molecule information ->'
      write(unit_,"(A, L5)")     'linear molucule    : ', linearQ
      write(unit_,"(A, A10)")     'potential type     : ', potenType
      write(unit_,"(A, A10)")     'molecule type      : ', moleculeType
      write(unit_,"(A, I5)")     '# of atom type     : ', nType

      select case(potenTypeENUM)

      case (OPLSAA)

        write(unit_,"(10(A10))")   'type', 'atom', 'eps', 'sig', 'charge', 'mass'
        do i = 1, nType
          write(unit_,"(I10, A10, 4(F10.6))") i, atomName(nint(atomType(i,1))), atomType(i,2), atomType(i,3), &
            atomType(i,4), atomType(i,5)
        end do

      case (FW)

        write(unit_,"(10(A10))")   'type', 'atom', 'mass'

        do i = 1, nType
          write(unit_,"(I10, A10, 4(F10.6))") i, atomName(nint(atomType(i,1))), atomType(i,5)
        end do

        write(unit_,"(A10)") 'A matrix'

        do i = 1, nType
          write(unit_,"(20F20.6)") A_FW(i, 1:nType)
        end do

        write(unit_,"(A10)") 'B matrix'

        do i = 1, nType
          write(unit_,"(20F20.6)") B_FW(i, 1:nType)
        end do

        write(unit_,"(A10)") 'C matrix'

        do i = 1, nType
          write(unit_,"(20F20.6)") C_FW(i, 1:nType)
        end do

      case (DREIDING)

        write(unit_,"(10(A10))")     'type', 'atom', 'zeta', 'D0', 'R0', 'charge', 'mass'
        do i = 1, nType

          write(unit_,"(I10, A10, 4(F10.6))") &
            i, &
            atomName(nint(atomType(i,1))), &
            zeta_DD(i,i), &
            D0_DD(i,i), &
            R0_DD(i,i), &
            charge(i), &
            mass(i)

        end do

        case default

      end select

      write(unit_,"(A, I5)")     'number of atoms    : ', nAtoms
      write(unit_,"(A)") 'type of each atom'
      do i = 1, nAtoms
        write(unit_,"(2(I5))") i, myType(i)
      end do

      write(unit_,"(A, I5)")     '# of bond in molecule     : ', nBond
      write(unit_,"(A)") 'bond pair'
      do i = 1, nBond
        write(unit_,"(2(I5))") bondIndices(i,1), bondIndices(i,2)
      end do

      write(unit_,"(A)") 'index type'
      do i = 1, nBondType
        write(unit_,"(I6, 2(A3))") i, atomName(nint(atomType(bondType(i,1),1))), &
          atomName(nint(atomType(bondType(i,2),1)))
      end do


      write(unit_,"(A)") 'bond parameter'
      do i = 1, nBondType
        write(unit_,"(I6, 2(F20.6))") i, kBond(i), bondLength(i)
      end do

      call PrintBondInformation(unit_)

      write(unit_,"(A, I5)")     '# of bend type     : ', nBendType
      write(unit_,"(A)") 'type index'
      do i = 1, nBendType
        write(unit_,"(I6, 3(A3))") i, atomName(nint(atomType(bendType(i,1),1))), &
          atomName(nint(atomType(bendType(i,2),1))), atomName(nint(atomType(bendType(i,3),1)))
      end do

      write(unit_,"(A)") 'bend parameter'
      do i = 1, nBendType
        write(unit_,"(I6, 3(F20.6))") i, kBend(i), theta0Deg(i), theta0(i)
      end do

      call PrintBendInformation(unit_)

      write(unit_,"(A, I5)")     '# of torsion type     : ', nTorsionType
      write(unit_,"(A)") 'type index'
      do i = 1, nTorsionType
        write(unit_,"(I6, 4(A3))") i, atomName(nint(atomType(TorsionType(i,1),1))), &
          atomName(nint(atomType(TorsionType(i,2),1))), &
          atomName(nint(atomType(TorsionType(i,3),1))), &
          atomName(nint(atomType(TorsionType(i,4),1)))
      end do

      write(unit_,"(A)") 'torsion parameter'
      do i = 1, nTorsionType
        write(unit_,"(I6, 4(F20.6))") i, v1(i), v2(i), v3(i), v4(i)
      end do



      call PrintTorsionInformation(unit_)
      call PrintNeighborQ(unit_)
      call PrintScalingQ(unit_)
      call PrintAdjacentIndices(unit_)
      call PrintMolTree(unit_)
      call PrintRelationship(unit_)
      call PrintGroup(unit_)

      write(unit_,"(A, F10.3)")  'mass of molecule   : ', molMass

      write(unit_,"()")

      write(unit_,"('properties = '/, A5, 5A10/)") 'num', 'eps', 'sig', 'mass', 'charge', 'color'

      do i = 1, nAtoms
        ai = myType(i)
        write(unit_,"(I5, 5f10.3)") i, eps(ai), sig(ai), mass(ai), charge(ai), color(ai)
      end do

      write(unit_,"(/5A10/)") 'cate num', 'epsilon', 'sigma', 'charge', 'number'
      do i=1, nType
        write(unit_,"(I10, 4F10.4)") i, eps(i), sig(i), charge(i), mass(i)
      end do
      write(unit_,"()")
      call PrintInertia (unit_)

      call PrintCheckedMotions(unit_)
      write(unit_,"(A10, 2(I10))") 'frame atom ', frameAtom1, frameAtom2

    else

      write(*,"(A/)")        'molecule information ->'
      write(*,"(A, L5)")     'linear molucule    : ', linearQ
      write(*,"(A, A10)")     'potential type     : ', potenType
      write(*,"(A, A10)")     'molecule type      : ', moleculeType
      write(*,"(A, I5)")     '# of atom type     : ', nType
      select case(potenTypeENUM)

      case (OPLSAA)

        write(*,"(10(A10))")   'type', 'atom', 'eps', 'sig', 'charge', 'mass'
        do i = 1, nType
          write(*,"(I10, A10, 4(F10.6))") i, atomName(nint(atomType(i,1))), atomType(i,2), atomType(i,3), &
            atomType(i,4), atomType(i,5)
        end do

      case (FW)

        write(*,"(10(A10))")   'type', 'atom', 'mass'

        do i = 1, nType
          write(*,"(I10, A10, 4(F10.6))") i, atomName(nint(atomType(i,1))), atomType(i,5)
        end do

        write(*,"(A10)") 'A matrix'

        do i = 1, nType
          write(*,"(20F20.6)") A_FW(i, 1:nType)
        end do

        write(*,"(A10)") 'B matrix'

        do i = 1, nType
          write(*,"(20F20.6)") B_FW(i, 1:nType)
        end do

        write(*,"(A10)") 'C matrix'

        do i = 1, nType
          write(*,"(20F20.6)") C_FW(i, 1:nType)
        end do

      case (DREIDING)

        write(*,"(10(A10))")     'type', 'atom', 'zeta', 'D0', 'R0', 'charge', 'mass'
        do i = 1, nType

          write(*,"(I10, A10, 4(F10.6))") &
            i, &
            atomName(nint(atomType(i,1))), &
            zeta_DD(i,i), &
            D0_DD(i,i), &
            R0_DD(i,i), &
            charge(i), &
            mass(i)

        end do

        case default

      end select
      write(*,"(A, I5)")     'number of atoms    : ', nAtoms
      write(*,"(A)") 'type of each atom'
      do i = 1, nAtoms
        write(*,"(2(I5))") i, myType(i)
      end do

      write(*,"(A, I5)")     '# of bond in molecule     : ', nBond
      write(*,"(A)") 'bond pair'
      do i = 1, nBond
        write(*,"(2(I5))") bondIndices(i,1), bondIndices(i,2)
      end do

      write(*,"(A)") 'index type'
      do i = 1, nBondType
        write(*,"(I6, 2(A3))") i, atomName(nint(atomType(bondType(i,1),1))), &
          atomName(nint(atomType(bondType(i,2),1)))
      end do

      write(*,"(A)") 'bond parameter'
      do i = 1, nBondType
        write(*,"(I6, 2(F20.6))") i, kBond(i), bondLength(i)
      end do

      call PrintBondInformation()

      write(*,"(A, I5)")     '# of bend type     : ', nBendType
      write(*,"(A)") 'type index'
      do i = 1, nBendType
        write(*,"(I6, 3(A3))") i, atomName(nint(atomType(bendType(i,1),1))), &
          atomName(nint(atomType(bendType(i,2),1))), atomName(nint(atomType(bendType(i,3),1)))
      end do

      write(*,"(A)") 'bend parameter'
      do i = 1, nBendType
        write(*,"(I6, 3(F20.6))") i, kBend(i), theta0Deg(i), theta0(i)
      end do

      call PrintBendInformation()

      write(*,"(A, I5)")     '# of torsion type     : ', nTorsionType
      write(*,"(A)") 'type index'
      do i = 1, nTorsionType
        write(*,"(I6, 4(A3))") i, atomName(nint(atomType(TorsionType(i,1),1))), &
          atomName(nint(atomType(TorsionType(i,2),1))), &
          atomName(nint(atomType(TorsionType(i,3),1))), &
          atomName(nint(atomType(TorsionType(i,4),1)))
      end do

      write(*,"(A)") 'torsion parameter'
      do i = 1, nTorsionType
        write(*,"(I6, 4(F20.6))") i, v1(i), v2(i), v3(i), v4(i)
      end do



      call PrintTorsionInformation()
      call PrintNeighborQ()
      call PrintScalingQ()
      call PrintAdjacentIndices()
      call PrintMolTree()
      call PrintRelationship()
      call PrintGroup()

      write(*,"(A, F10.3)")  'mass of molecule   : ', molMass

      write(*,"()")

      write(*,"('properties = '/, A5, 5A10/)") 'num', 'eps', 'sig', 'mass', 'charge', 'color'

      do i = 1, nAtoms
        ai = myType(i)
        write(*,"(I5, 5f10.3)") i, eps(ai), sig(ai), mass(ai), charge(ai), color(ai)
      end do

      write(*,"(/5A10/)") 'cate num', 'epsilon', 'sigma', 'charge', 'number'
      do i=1, nType
        write(*,"(I10, 4F10.4)") i, eps(i), sig(i), charge(i), mass(i)
      end do
      write(*,"()")
      call PrintInertia

      call PrintCheckedMotions()

      write(*,"(A10, 2(I10))") 'frame atom ', frameAtom1, frameAtom2

    end if

    end subroutine PrintMolInf

    !********************

    subroutine MakeMolView()

    integer :: i, ai

    if (.not.getMolInfQ) stop '입력받은 분자정보가 없습니다.'

    open(UNIT_MOL, file='molview.view')

    write(UNIT_MOL,*) nAtoms

    do i = 1, nAtoms
      ai = myType(i)
      write(UNIT_MOL,*) color(ai), molCrdX(i), molCrdY(i), molCrdZ(i)
    end do

    write(UNIT_MOL,*) 0
    write(UNIT_MOL,*) 10

    close(UNIT_MOL)

    end subroutine MakeMolView

    !********************

    subroutine PrintInertia (unit_)

    implicit none

    integer :: i
    integer, optional, intent(in) :: unit_
    real(kind=r8) :: x, y, z, m, dsq

    InertiaTensor = 0.0_r8

    do i = 1, nAtoms

      x = molCrdX(i); y = molCrdY(i); z = molCrdZ(i); m = mass(myType(i))
      dsq = x**2 + y**2 + z**2

      InertiaTensor(1,1) = InertiaTensor(1,1) + m*(dsq - x**2)
      InertiaTensor(1,2) = InertiaTensor(1,2) - m*x*y
      InertiaTensor(1,3) = InertiaTensor(1,3) - m*x*z
      InertiaTensor(2,1) = InertiaTensor(2,1) - m*x*y
      InertiaTensor(2,2) = InertiaTensor(2,2) + m*(dsq - y**2)
      InertiaTensor(2,3) = InertiaTensor(2,3) - m*y*z
      InertiaTensor(3,1) = InertiaTensor(3,1) - m*x*z
      InertiaTensor(3,2) = InertiaTensor(3,2) - m*z*y
      InertiaTensor(3,3) = InertiaTensor(3,3) + m*(dsq - z**2)

    end do

    if (present(unit_)) then

      write(unit_,"('inertia tensor(g*A^2/mol) = '/)")
      write(unit_,"(3F20.6//,3F20.6//,3F20.6//)") Transpose(InertiaTensor)

    else

      write(*,"('inertia tensor(g*A^2/mol) = '/)")
      write(*,"(3F20.6//,3F20.6//,3F20.6//)") Transpose(InertiaTensor)

    end if

    end subroutine PrintInertia

    !********************
    subroutine RescaleBond  ! not used

    implicit none

    real(kind=r8), dimension(3) :: bondVec(3)
    real(kind=r8) :: currentNorm
    real(kind=r8) :: bRatio ! bondLength / CurrentNorm

    integer :: idx1, idx2, bond, i, path, n, nOther, nIndex

    do bond = 1, nBond

      idx1 = index1InBond(bond)
      idx2 = index2InBond(bond)
      bondVec(1) = molCrdX(idx2) - molCrdX(idx1)
      bondVec(2) = molCrdY(idx2) - molCrdY(idx1)
      bondVec(3) = molCrdZ(idx2) - molCrdZ(idx1)

      currentNorm = sqrt(bondVec(1)**2 + bondVec(2)**2 + bondVec(3)**2)
      bRatio = bondLength(myBondType(bond)) / currentNorm

      if (abs(bRatio - 1.0_r8) > 1.e-7) then
        write(*,"(I3, A, F10.3, A, F10.6)") bond, ' rescale bond length', &
          currentNorm, ' to ', bondLength(myBondType(bond) )
      end if

      do i = 1, 4
        if (molTree(idx1,i,1) == idx2) then
          path = i
          exit
        end if
      end do

      nOther = molTree(idx1,path,0)

      do n = 1, nOther

        nIndex = molTree(idx1,path,n)
        molCrdX(nIndex) = molCrdX(nIndex) - bondVec(1)
        molCrdY(nIndex) = molCrdY(nIndex) - bondVec(2)
        molCrdZ(nIndex) = molCrdZ(nIndex) - bondVec(3)

      end do ! n-loop

      bondVec = bondVec * bRatio ! rescale bond vector, using vector operation

      do n = 1, nOther

        nIndex = molTree(idx1,path,n)
        molCrdX(nIndex) = molCrdX(nIndex) + bondVec(1)
        molCrdY(nIndex) = molCrdY(nIndex) + bondVec(2)
        molCrdZ(nIndex) = molCrdZ(nIndex) + bondVec(3)

      end do ! n-loop

    end do ! bond-loop

    end subroutine RescaleBond

    !   -------------------------------------------------
    subroutine GroupingMolecule()
    !   -------------------------------------------------

    implicit none

    integer :: i, j

    !   -------------------------------------------------

    select case(rcutType)

    case('mol','MOL')
      call MolGrouping()
    case('atom','ATOM')
      call AtomGrouping()
    case('group','GROUP')

      call GroupGrouping() ! OPLSAA and FW but limited...
      !if (potenTypeENUM == OPLSAA) call GroupGrouping() ! charged molecule
      !if (potenTypeENUM == FW    ) call FWGrouping()    ! carbon center based

      case default

      stop 'what the groupping?'

    end select

    contains
    !       ------------------------
    subroutine MolGrouping()

    implicit none

    nLocalGroup = 1

    allocate(localGroup(nLocalGroup))
    allocate(localGroup(1)%member(0:nAtoms))

    localGroup(1)%member(0) = nAtoms

    do i = 1, nAtoms

      localGroup(1)%member(i) = i

    end do

    end subroutine MolGrouping

    !       ------------------------
    subroutine AtomGrouping()

    implicit none

    nLocalGroup = nAtoms

    allocate(localGroup(nLocalGroup))

    do i = 1, nLocalGroup
      allocate(localGroup(i)%member(0:1))
    end do
    do i = 1, nAtoms

      localGroup(i)%member(0) = 1
      localGroup(i)%member(1) = i

    end do

    end subroutine AtomGrouping

    !       ------------------------
    subroutine GroupGrouping() ! it is maybe changed..

    implicit none

    logical :: neutralQ = .false.

    integer :: a, nMem, nChecked, nAdjac ! number of adjacent atom

    integer :: nadj, adj

    real(kind=r8) :: myCharge, chargeSum, adjCharge

    nLocalGroup = 0
    nChecked    = 0

    do a = 1, nAtoms

      neutralQ = .false.
      if (adjacentIndices(a,0) == 1) cycle ! neglect single bond atom

      nLocalGroup = nLocalGroup + 1 ! count number of local group

    end do ! a-loop

    ! memory allocation
    allocate(localGroup(nLocalGroup))
    do i = 1, nLocalGroup
      allocate(localGroup(i)%member(0:5))
    end do

    nLocalGroup = 0

    ! indexing step
    do a = 1, nAtoms
      !call system('pause')
      neutralQ = .false.
      if (adjacentIndices(a,0) == 1) cycle ! neglect single bond atom

      nLocalGroup = nLocalGroup + 1
      !call system('pause')
      myCharge = charge(myType(a))

      if (abs(myCharge) < 1.e-5) neutralQ = .true.

      chargeSum = 0.0_r8
      !call system('pause')
      nAdjac = adjacentIndices(a,0)
      !call system('pause')
      localGroup(nLocalGroup)%member(0) = 1
      localGroup(nLocalGroup)%member(1) = a
      nChecked = nChecked + 1
      do nadj = 1, nAdjac
        !call system('pause')
        adj = adjacentIndices(a,nadj)
        if (neutralQ) then  ! for neutral case

          if (adjacentIndices(adj,0) /= 1) cycle ! neglect non-hydrogen

        else ! for charged case

          adjCharge = charge(myType(adj))
          if (myCharge * adjCharge >= -1.e-13) cycle ! neglect neutral or same chargeed atom

        end if

        localGroup(nLocalGroup)%member(0) = localGroup(nLocalGroup)%member(0) + 1
        nMem = localGroup(nLocalGroup)%member(0)
        localGroup(nLocalGroup)%member(nMem) = adj
        !write(*, "('group: ', I4, 'member: ', I4)") nLocalGroup, nMem, adj
        nChecked = nChecked + 1
      end do ! nadj-loop

    end do ! a-loop

    if (nChecked /= nAtoms) stop 'sum of group member /= # of atoms'

    end subroutine GroupGrouping
    !       --------------------------

    end subroutine GroupingMolecule



    !   ----------------------------------------
    subroutine PrintGroup(unit_)
    !   ----------------------------------------
    implicit none

    integer, intent(in), optional  :: unit_
    integer :: i, j, nMem
    integer, dimension(:), pointer :: mem
    !   ----------------------------------------

    if (present(unit_)) then
      write(unit_, "('number of group: ', I5)") nLocalGroup
      do i = 1, nLocalGroup
        mem => localGroup(i)%member
        nMem = mem(0)
        do j = 1, nMem
          write(unit_, "('group: ', I3, ', member: ', I3)") i, mem(j)
        end do ! j-loop
      end do ! i-loop
    else
      write(*, "('number of group: ', I5)") nLocalGroup
      do i = 1, nLocalGroup
        mem => localGroup(i)%member
        nMem = mem(0)
        do j = 1, nMem
          write(*, "('group: ', I3, ', member: ', I3)") i, mem(j)
        end do ! j-loop
      end do ! i-loop

    end if

    end subroutine PrintGroup








    !   ----------------------------------------
    subroutine DeallocateGroup()
    !   ----------------------------------------
    implicit none

    integer :: i
    !   ----------------------------------------

    do i = 1, nLocalGroup
      if (allocated(localGroup)) deallocate(localGroup(i)%member)
    end do

    if (allocated(localGroup)) deallocate(localGroup)

    end subroutine DeallocateGroup



    subroutine MolDeallocate() ! destructor in C++ or finalize in JAVA

    implicit none

    if (allocated(molCrdX))   deallocate(molCrdX, molCrdY, molCrdZ)
    if (allocated(myType))    deallocate(myType)
    if (allocated(myAtomNum)) deallocate(myAtomNum)
    if (allocated(bondIndices))    deallocate(bondIndices)
    if (allocated(bendIndices))    deallocate(bendIndices)
    if (allocated(torsionIndices)) deallocate(torsionIndices)
    if (allocated(myBendType))     deallocate(myBendType)
    if (allocated(myBendType))     deallocate(myBendType)
    if (allocated(myTorsionType))  deallocate(myTorsionType)
    if (allocated(neighborQ))      deallocate(neighborQ)
    if (allocated(scalingQ))       deallocate(scalingQ)
    if (allocated(adjacentIndices)) deallocate(adjacentIndices)
    if (allocated(molTree))         deallocate(molTree)
    if (allocated(index1InBond)) deallocate(index1InBond, index2InBond,         & ! array for efficient calculation of energy
    &                                     index1InBend, index2InBend, index3InBend, &
      &          index1InTorsion,  index2InTorsion, index3InTorsion, index4InTorsion, &
      &                     bond12InBend, bond23InBend, direc12InBend, direc23InBend, &
      &                            bond12InTorsion, bond23InTorsion, bond34InTorsion, &
      &                         direc12InTorsion, direc23InTorsion, direc34InTorsion, &
      &       bend123InTorsion, bend234InTorsion, direc123InTorsion, direc234InTorsion)
    if (allocated(myTypeInBend))    deallocate(myTypeInBend)
    if (allocated(myTypeInTorsion)) deallocate(myTypeInTorsion)
    if (allocated(myTypeInBond)) deallocate(myTypeInBond)
    call DeallocateGroup()

    !if (allocated()) deallocate()
    !if (allocated()) deallocate()


    write(*,"('Mol array is deallocated')")
    write(1,"('Mol array is deallocated')")

    end subroutine MolDeallocate

    !********************

  end module mod_molecule_information

  !============================================================
