  !--------------------------------------------------------------
  !subroutine GetExpSolidInf (filename)
  !
  !subroutine AutoOperation
  !subroutine OKsign          ! check input data reading
  !subroutine ExpMCSolid
  !subroutine ExpandArray(j,fb)
  !
  !subroutine GetLambda(lambda_)
  !
  !subroutine PrintExpSolidInf (unit_)
  !subroutine PrintFreeEnergy!
  !
  !subroutine ExpSolDeallocate
  !
  ! non-using subroutines
  !subroutine SwapIndex (i)
  !subroutine SpecialMove (i, delE)
  !--------------------------------------------------------------

  module mod_expans_solid

  use mod_univ_const,           only : CHARGE_CONST, KtokJ
  use mod_character_tool,       only : IntToCha
  !use mod_time_tool,            only : Delay
  use mod_random_tool,          only : Ranf, RandomSettingAll

  use mod_system_information,   only : molNum=>molNum, rCut=>rCut, beta=>beta, printStep=>printStep, temper=>temper, &
    simulTypeENUM => simulTypeENUM, EXP_SOLUTE => EXP_SOLUTE
  use mod_molecule_information, only : nAtoms=>nAtoms, eps=>eps, sig=>sig, cha=>charge, &
    mcX=>molCrdX, mcY=>molCrdY, mcZ=>molCrdZ, mixingRule=>mixingRule, &
    color=>color
  use mod_config_information,   only : rX => rX, rY => rY, rZ => rZ, &
    hBox => hBox, ihBox => ihBox, vol => volume, &
    rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew , &
    soluteIndex, UpdateCoordinate

  use mod_movement_information, only : ShrinkSolute, ShrinkOldNew
  use mod_energy_core,          only : CalculateEnergy

  implicit none

  private

  integer, parameter :: r8 = selected_real_kind(15,300)

  integer, parameter :: MXSTATE = 10000, MXBIN = 10000, MAXA = 120 ! Max Atom.. C70 never..
  integer, parameter :: FLUID_TYPE = 1, SOLID_TYPE = FLUID_TYPE + 1 ! enumerator for expanded ensemble type.

  logical :: adjustModeQ = .true., weightReadQ, weightOK, dimerSwapQ

  integer :: totalPartition, state = 1, repeatStep, MaxXi, partition, eachStep = 1

  integer       :: fmax, fmin, bmax, bmin
  integer       :: ixi, ieestep, icycle, ixiold
  real(kind=r8) :: xi0, xi1, bin, mxbin_1, ebinmax, ene, eneold, relax, maxdw, eReal, eEin, delEs
  real(kind=r8), allocatable, dimension(:)   :: xi, wxi, delw, delwOld, Pxi, nforxi, nbakxi, forxi, bakxi, oforbak
  real(kind=r8), allocatable, dimension(:,:) :: nbeufor, nbeubak

  public :: icycle, ene, eneold, partition, wxi, MaxXi, delw, delwOld, weightOK, ixi, delEs, eReal, eEin, xi, dimerSwapQ, &
    adjustModeQ, bin

  ! Enumerator.
  public :: FLUID_TYPE, SOLID_TYPE


  public :: GetExpSolidInf, &
    PrintExpSolidInf, &
    ExpMCSolid, &
    AutoOperation, &
    PrintFreeEnergy, &
    ExpSolDeallocate, &
    SpecialMove, &
    OKsign, &
    GetLambda, &
    isLastPartition, &
    RunNextPartitionMain, &
    get_ratio

  contains

  subroutine GetExpSolidInf (filename)

  implicit none

  character(len=31), intent(in) :: filename
  integer, parameter :: UNIT_EEM = 10, UNIT_PART = 11, UNIT_EEMOUT = 12
  integer            :: i, j, dummyPart, dmaxXi
  integer            :: openError
  real(kind=r8), dimension(0:MXSTATE) :: dweight

  real(kind=r8)      :: dummyUp, dummyDown, ls, le

  weightOK    = .false.
  adjustModeQ = .true.

  open(UNIT_EEM, file=filename)

  read(UNIT_EEM,*) totalPartition

  open(UNIT_PART, file='partition.part', status='old', iostat=openError)

  ! partition.part 파일이 없으면 첫번째 파티션입니다.
  if (openError >  0) then

    partition = 1

    ! eemout 파일의 타이틀을 작성합니다.
    open(UNIT_EEMOUT, file='eemout.txt', status='replace', position='append')

    write(UNIT_EEMOUT, "(2A30)") 'partition', 'delta free energy'

    close(UNIT_EEMOUT)

    call system("title 01")

  end if

  if (openError == 0) then

    read(UNIT_PART,*) partition
    close(UNIT_PART, status='delete')

  end if

  do i = 1, totalPartition

    read(UNIT_EEM,*) dummyPart, dmaxXi, ls, le

    if (dummyPart == partition) then

      maxXi = dmaxXi
      xi0   = ls
      xi1   = le

    end if

  end do
  ! ????
  if (MaxXi > MXSTATE) stop 'Max xi > Max State'

  read(UNIT_EEM,*) weightReadQ

  if (weightReadQ) then

    do i = 1, totalPartition

      read(UNIT_EEM,*) dummyPart!, dweight(0:MaxXi)

      if (dummyPart == partition) then

        backspace UNIT_EEM
        read(UNIT_EEM,*) dummyPart, dweight(0:MaxXi)

        wxi(0:MaxXi) = dweight(0:MaxXi)

      end if

    end do

    weightOK    = .true.
    adjustModeQ = .false.

  end if

  close(UNIT_EEM)

  ! initialization

  allocate(xi(0:MaxXi), wxi(0:MaxXi), delw(0:MaxXi), delwOld(0:MaxXi), &
    &Pxi(0:MaxXi), nforxi(0:MaxXi), nbakxi(0:MaxXi), forxi(0:MaxXi), &
    &bakxi(0:MaxXi), oforbak(0:MaxXi))

  allocate(nbeufor(0:MaxXi,0:0), nbeubak(0:MaxXi,0:0))

  do i = 0, MaxXi

    xi(i) =  xi0 + ( xi1 - xi0 ) * dble(i) / dble( MaxXi )
    if (.not.weightOK) wxi(i)= 0.0_r8

  end do

  ixi = 0

  do i = 0, maxXi - 1

    delw(i) = wxi(i+1) - wxi(i)

  end do

  ieestep = 0
  relax   = 1.0

  Pxi    = 0    ! probabilities
  nforxi = 0
  nbakxi = 0
  forxi  = 0
  bakxi  = 0

  oforbak = 0.0_r8     ! for what?
  nbeufor = 0.0_r8
  nbeubak = 0.0_r8

  fmax = 0
  fmin = 0
  bmax = 0
  bmin = 0

  dimerSwapQ = .false.

  if (simulTypeENUM == EXP_SOLUTE) bin = 0.1_r8

  end subroutine GetExpSolidInf

  logical Function isLastPartition()

    isLastPartition = (totalPartition == partition)

  end Function isLastPartition

  real(kind=r8) function get_ratio()

      real(kind=r8) min_, max_
      min_ = minval(Pxi(0:MaxXi))
      max_ = maxval(Pxi(0:MaxXi))
      if (min_ < 1.e-10) then
        get_ratio = -1.0
      else
        get_ratio =  max_ / min_
      end if

  end function get_ratio

  subroutine RunNextPartitionMain()

    implicit none

    integer, parameter :: UNIT_PART = 17

    character(len=6), parameter :: start = 'start '
    character(len=1), parameter :: ap    = char(34) ! char(34) = "
    character(len=4), parameter :: main  = 'main'
    character(len=51) :: command

    integer :: nextPartition

    nextPartition = partition + 1

    !다음 Main Program이 인식할 수 있도록 다음 파티션을 가리킨다.
    open(UNIT_PART, file='partition.part')

    write(UNIT_PART,"(I2)") nextPartition ! 현재 파티션의 다음 파티션 인덱스를 작성합니다.
    write(UNIT_PART,"(A)")  'end'

    close(UNIT_PART)

    ! 다음 Main을 실행하기 위해 System subroutine에 넣을 커멘드를 생성합니다.
    ! 예를들어 start "12" main 으로 생성됩니다.
    command = start//ap//IntToCha(nextPartition)//ap//' '//main

    write(*,*) 'check command: ', command
    !call system('pause')

    ! 다음 Partition 계산을 수행할 Main 프로그램을 실행합니다.
    call system(command)

  end subroutine RunNextPartitionMain

  subroutine GetLambda(lambda_)

  implicit none

  real(kind=r8), intent(out) :: lambda_

  lambda_ = xi(ixi)

  end subroutine GetLambda

  subroutine PrintExpSolidInf (unit_)

  implicit none

  integer, intent(in), optional :: unit_

  if (present(unit_)) then

    write(unit_,"(A/)") 'expanded ensemble information ->'

    write(unit_,"('max bin = ', I15)") MXBIN
    write(unit_,"('total partition number = ', I5)")  totalPartition
    write(unit_,"('current partition      = ', I5)")  partition
    write(unit_,"('max Xi                 = ', I5)")  maxXi
    write(unit_,"('lambda = ',11F10.6)") xi(0:MaxXi)
    write(unit_,"('weight = ',11F10.6)") wxi(0:MaxXi)
    write(unit_,"('Is exist initial weight? (T/F) = ', L2 )") weightReadQ
    write(unit_,"()")

  else

    write(*,"(A/)") 'expanded ensemble information ->'

    write(*,"('max bin = ', I15)") MXBIN
    write(*,"('total partition number = ', I5)")  totalPartition
    write(*,"('current partition      = ', I5)")  partition
    write(*,"('max Xi                 = ', I5)")  maxXi
    write(*,"('lambda = ',11F10.6)") xi(0:MaxXi)
    write(*,"('weight = ',11F10.6)") wxi(0:MaxXi)
    write(*,"('Is exist initial weight? (T/F) = ', L2 )") weightReadQ
    write(*,"()")

  end if

  end subroutine PrintExpSolidInf

  subroutine AutoOperation

  implicit none

  integer            :: checkPart
  integer, parameter :: unit_ = 11
  logical            :: readingQ
  character(len=6), parameter :: start = 'start '
  character(len=1), parameter :: ap    = char(34) ! char(34) = "
  character(len=4), parameter :: main  = 'main'
  character(len=51) :: command

  integer :: part_

  if (partition /= 1) return

  do part_ = 2, totalPartition

    command = start//ap//IntToCha(part_)//ap//' '//main
    !ex) start "12" main

    open(1, file='partition.part')

    write(1,"(I2)") part_
    write(1,"(A)")  'end'

    close(1)

    write(*,"('run partition = ', I2)") part_
    write(*, "(A50)") "command: "//command
    call system(command)

    !call Delay(1)
    !do while(.not.readingQ)
1   continue
    open(unit_, file='oksign.sign', status='old', err=1)

2   continue
3   continue
    read(unit_,*,err=2, end=3) readingQ, checkPart

4   continue
    close(unit_, status='delete', err=4)
    !end do
    if (.not.readingQ .or. part_ /= checkPart) call system('pause')

  end do

5 continue
  open(unit_, file='oksign.sign', err=5)
6 continue
  close(unit_, status='delete', err=6)

  open(1, file='eemout.txt', status='replace', position='append')

  write(1, "(2A30)") 'partition', 'delta free energy'

  partition = 1

  if (partition == 1) call system("title 01")

  end subroutine AutoOperation

  subroutine OKsign ! check input data reading

  implicit none

  integer, parameter :: unit_ = 10

  if (partition == 1) return

1 continue

  open(unit_, file='oksign.sign', status='replace', err=1)

3 continue
  write(unit_,"('T', I20)", err=3) partition

2 continue

  close(unit_, err=2)

  end subroutine OKsign

  subroutine ExpMCSolid(EXP_TYPE, eReal_, eEin_, icycle_)

  implicit none

  integer, intent(in) :: EXP_TYPE ! FLUID_TYPE or SOLID_TYPE
  logical, save       :: initialAdjustFQ = .true., initialAdjustBQ = .true.
  real(kind=r8), intent(inout) :: eReal_, eEin_
  integer, intent(in) :: icycle_

  logical                   :: overlapQ
  ! why so low? kind only 4byte?
  real(kind=r8), parameter  :: INFINITY = 1.e30, ZERO = 1.e-8

  integer                   :: eestep, i, j, iter, Imin, Imax, iarrmin, iarrmax
  real(kind=r8)             :: e, ec, arg, prob, paccfor, paccbak, probeu, beu, dw, dwn, farg, barg, nfarg, nbarg  ! energy of old configuration and state

  real(kind=r8)             :: muTot, Pmin, Pmax, eVDW, eQ

  icycle = icycle_
  if (mod(icycle,2) == 1) return

2 continue

  ixiold = ixi   ! save old configuration.
  !        if (EXP_TYPE == SOLID_TYPE) then
  !
  !            eneold = xi(ixiold) * eReal_ + (1.0-xi(ixiold)) * eEin_   ! energy of old state ixi
  !
  !        else if (EXP_TYPE == FLUID_TYPE) then
  !
  !        end if
  delEs  = 0.0

  ! choose random lambda transition.
  !----- begining of lambda ----------------
  if( ixiold == 0 ) then

    if( ranf() > 0.5 ) then                 ! stay

      ixi = 0

    else

      ixi = 1
      nforxi(ixiold) = nforxi(ixiold) + 1   ! number of forward trials

    endif
    !----- end of lambda --------------------
  else if( ixiold == MaxXi ) then

    if( ranf() > 0.5 ) then                 !  stay

      ixi = MaxXi

    else

      ixi = MaxXi - 1
      nbakxi(ixiold) = nbakxi(ixiold) + 1

    endif
    !----- other cases ----------------------
  else

    if( ranf() > 0.5 ) then ! forward test

      ixi = ixiold + 1
      nforxi(ixiold) = nforxi(ixiold) + 1

    else                    ! backward test

      ixi = ixiold - 1
      nbakxi(ixiold) = nbakxi(ixiold) + 1

    endif

  endif

  if ( ixi == ixiold ) go to 1             ! do nothing and bypass


  ! make delta energy
  if (EXP_TYPE == SOLID_TYPE) then

    eneold = xi(ixiold) * eReal_ + (1.0 - xi(ixiold)) * eEin_   ! energy of old state ixi
    ene    = xi(ixi)    * eReal_ + (1.0 - xi(ixi))    * eEin_

  else if (EXP_TYPE == FLUID_TYPE) then

    ! fluid delta energy is needed.
    call ShrinkOldNew(xi(ixi), rxiold, ryiold, rziold, &
      rxinew, ryinew, rzinew)

    call CalculateEnergy(soluteIndex, 'OLD', eVDW, eQ, xi(ixiold))
    eneold = eVDW + eQ
    call CalculateEnergy(soluteIndex, 'NEW', eVDW, eQ, xi(ixi))
    ene    = eVDW + eQ

  end if


  arg = - beta * (ene - eneold)

  if (adjustModeQ) then ! save test energy distribution.

    j = nint( (-arg) / bin )                      ! beta* del U = -arg, bin size = 0.1, range -10 ~ 10 in unit of kT

    if (EXP_TYPE == FLUID_TYPE) then

      if (j >=  MXBIN) j =  MXBIN ! because hardcore repulsion.
      if (j <= -MXBIN) j = -MXBIN ! not tested

    end if

    if( ixiold < ixi ) then ! forward direction.

      if (abs(j) /= abs(MXBIN)) then

        call ExpandArray (j,1)
        nbeufor(ixiold,j) = nbeufor(ixiold,j) + 1     ! count transition energy

      end if
      arg = arg + delw(ixiold)

    else ! ixiold > ixi    ! backward direction.

      if (abs(j) /= abs(MXBIN)) then

        call ExpandArray (j,2)
        nbeubak(ixiold,j) = nbeubak(ixiold,j) + 1     ! count transition energy

      end if
      arg = arg - delw(ixi)
      ! not delw(ixiold)
    endif

  else  ! no adjustment ! original method

    arg = arg + wxi(ixi) - wxi(ixiold)

  endif ! chk_adaptive


  if( arg > 0.d0 ) then

    prob = 1.d0

  else

    prob = dexp( arg )

  endif

  ! acceptance test
  if ( Ranf() < prob ) then ! accept

    ! lambda change
    if( ixiold < ixi ) then                 ! successful forward trial

      forxi(ixiold) = forxi(ixiold) + 1

    else if( ixiold > ixi ) then            ! successful backward trial

      bakxi(ixiold) = bakxi(ixiold) + 1

      !           else ! ixiold = ixi                     ! do noting

    endif

    !delEs = ene - eneold ! why ???

    if (EXP_TYPE == SOLID_TYPE) then

      ! do nothing

    else if (EXP_TYPE == FLUID_TYPE) then

      ! coordinate and energy update
      call ShrinkSolute(xi(ixi)) ! change saved ratio to xi(ixi)
      call UpdateCoordinate(soluteIndex) ! then update coordinate

      eReal_ = eReal_ + ene - eneold
      eEin_  = eEin_  + ene - eneold ! eEin = uInter in fluid case.

    end if

  else                      ! reject

    ixi = ixiold   ! back to origin state
    ene = eneold

    delEs = 0.0

  end if

  if (adjustModeQ) then  !new  !every 100 cycle
    if (mod(icycle,100) == 0) then ! adjust cycle

      do i = 0, MaxXi - 1

        if (nforxi(i) < 1 .or. nbakxi(i+1) < 1)      cycle  ! after at least 1 samplings
        if (nforxi(i) + nbakxi(i+1) < oforbak(i) + 10) cycle  ! do not renew weight before 10 times of transition

        dw = delw(i)
        relax = 1.

10000   do iter = 1, 500

          farg    = 0.
          nfarg   = 0.
          paccfor = 0.
          ! numerical integration.. square sum.
          do j = fmin, fmax

            if (nbeufor(i,j) == 0) cycle                      ! unvisited bin

            probeu = dble( nbeufor(i,j) ) / dble( nforxi(i) )

            !if( probeu < 0.001 ) cycle                         ! too small prob

            beu = dble(j) * bin
            arg = - beu + dw

            if (arg > ZERO) then

              paccfor = paccfor + probeu
              farg  = farg    + arg
              nfarg = nfarg   + 1.

            else

              paccfor = paccfor + exp(arg) * probeu
              farg  = farg    + arg
              nfarg = nfarg   + 1.

            endif

          end do ! j

          barg    = 0.
          nbarg   = 0.
          paccbak = 0.
          do j = bmin, bmax

            if (nbeubak(i+1,j) == 0) cycle  ! unvisited bin

            probeu = dble( nbeubak(i+1,j) ) / dble( nbakxi(i+1) )

            !if( probeu < 0.001 ) cycle

            beu = dble(j) * bin
            arg = - beu - dw

            if (arg > ZERO) then

              paccbak = paccbak + probeu
              barg    = barg    + arg
              nbarg   = nbarg   + 1.

            else

              paccbak = paccbak + exp(arg) * probeu
              barg    = barg    + arg
              nbarg   = nbarg   + 1.

            endif

          end do ! j

          if(paccfor < ZERO .or. paccbak < ZERO) then  ! very rare case if it occurs

            write(* ,*) 'warning: adaptive eemc, acc prob zero'
            write(1 ,*) 'warning: adaptive eemc, acc prob zero'
            write(* ,"(3i8,6f20.8)") i, nint(nforxi(i)), nint(nbakxi(i+1)), &
              paccfor, paccbak, &
              farg, barg, &
              nfarg, nbarg !check
            write(1 ,"(3i8,6f20.8)") i, nint(nforxi(i)), nint(nbakxi(i+1)), &
              paccfor, paccbak, &
              farg, barg, &
              nfarg, nbarg !check

            if (paccfor < ZERO) then

              if (nfarg < 1) exit
              write(*,"('mean farg = ', F20.6)") farg / nfarg
              write(1,"('mean farg = ', F20.6)") farg / nfarg
              dwn = dw - farg / nfarg

            end if

            if (paccbak < ZERO) then

              if (nbarg < 1) exit
              write(*,"('mean barg = ', F20.6)") barg / nbarg
              write(1,"('mean barg = ', F20.6)") barg / nbarg
              dwn = dw + barg / nbarg

            end if

          else

            dwn = dw - log(paccfor / paccbak) * relax

          endif

          if (abs(dwn-dw) < 0.001) then

            delw(i) = dwn
            oforbak(i) = nforxi(i) + nbakxi(i+1)
            exit  ! converged

          endif

          dw = dwn

        end do ! iter

        !       // error handling  --not happen usually
        if (iter > 500) then

          write(* ,"('--warning: weight not converged! ')")
          if (relax > 0.2) then

            relax = relax * 0.8
            write(* ,"('--relax factor = ',f6.4)") relax
            go to 10000   ! re calculate..

          else

            write(* ,"('--not converged! but will use current weight')")
            delw(i) = dwn
            oforbak(i) = nforxi(i) + nbakxi(i+1)

          endif

        endif
        !        // error handling

      end do ! i

      write(* ,"(/i9,9x,21(f10.4))") icycle,( delw(i), i= 0, MaxXi-1)
      write(1,"(i9,9x,21(f6.2)) ")   icycle,( delw(i), i= 0, MaxXi-1)
      write(*,"('fmin: ', I9,',fmax: ', I9,',bmin: ', I9,',bmax: ', I9)") fmin, fmax, bmin, bmax
      write(1,"('fmin: ', I9,',fmax: ', I9,',bmin: ', I9,',bmax: ', I9)") fmin, fmax, bmin, bmax

      maxdw = 0.0
      if ( minval(dabs(delw(0:maxXi-1))) > 0.0001 ) then ! non zero condition

        maxdw = maxval(dabs((delwOld(0:maxXi-1) - delw(0:maxXi-1))))

        if ( maxval(dabs(delw(0:maxXi-1))) < 0.1 ) then ! for too low value of weight
          ! use relative error
          maxdw = maxval(dabs((delwOld(0:maxXi-1) - delw(0:maxXi-1))/delw(0:maxXi-1)))

        end if

        pmin = minval(Pxi(0:MaxXi))
        pmax = maxval(Pxi(0:MaxXi))
        if (                                               &
          maxdw < 0.01 .and.                            &
          minval(Pxi(0:MaxXi)) > 500.0 .and.            &
          pmax / pmin < 2.0_r8                          &
          ) then ! 1% error .and. icycle > 10000

        adjustModeQ = .false.
        weightOK    = .true.
        Pxi = 0.0   ! initialize prob
        write(*,"('start hopping stage')")
        write(1,"('start hopping stage')")
        wxi(0) = 0.0
        do i = 1, maxXi

          wxi(i) = wxi(i-1) + delw(i-1)

        end do
        write(* ,"(/i9,9x,21(f15.4))") icycle,( wxi(i), i= 0, MaxXi)
        write(1 ,"(/i9,9x,21(f15.4))") icycle,( wxi(i), i= 0, MaxXi) ! unit 1 : monitoring file

        end if

        delwOld(0:maxXi-1) = delw(0:maxXi-1)

      end if



      !============== solid only ===================
      if (EXP_TYPE == SOLID_TYPE .or. EXP_TYPE == FLUID_TYPE) then

        pmin = minval(Pxi(0:MaxXi))
        pmax = maxval(Pxi(0:MaxXi))
        imin = minloc(Pxi(0:MaxXi),1) - 1
        imax = maxloc(Pxi(0:MaxXi),1) - 1

        if ( icycle > 1000 ) then
          if ( pmin == 0 .or. pmax / pmin > 2.0 ) then


            if (EXP_TYPE == SOLID_TYPE) then

              ixiold = ixi  ! save to old, current subensemble
              ixi = imin    ! change to minimum subensemble

              write(*,"('ixi is changed to minimum prob spot; ', I5)") ixi
              write(1,"('ixi is changed to minimum prob spot; ', I5)") ixi

            end if

            if (EXP_TYPE == FLUID_TYPE) then

              ! do nothing by pass

              !call ShrinkOldNew(xi(ixi), rxiold, ryiold, rziold, &
              !           rxinew, ryinew, rzinew)

              !call CalculateEnergy(soluteIndex, 'OLD', eVDW, eQ, xi(ixiold))
              !eneold = eVDW + eQ
              !call CalculateEnergy(soluteIndex, 'NEW', eVDW, eQ, xi(ixi))
              !ene    = eVDW + eQ
              ! coordinate and energy update
              !call ShrinkSolute(xi(ixi)) ! change saved ratio to xi(ixi)
              !call UpdateCoordinate(soluteIndex) ! then update coordinate

              !eReal_ = eReal_ + ene - eneold
              !eEin_  = eEin_  + ene - eneold ! eEin = uInter in fluid case.

            end if

          end if
        end if

      end if ! EXP_TYPE
      !============== solid only ===================

      write(*,"('max dw = ', F10.6, '   min Pxi, i = ', F6.0, I5, '   max Pxi, i = ', F6.0, I5)") &
        & maxdw, pmin, imin, pmax, imax
      write(1,"('max dw = ', F10.6, '   min Pxi, i = ', F6.0, I5, '   max Pxi, i = ', F6.0, I5)") &
        & maxdw, pmin, imin, pmax, imax

    end if ! adjust cycle

  end if ! chk_adaptive

1 continue

  Pxi(ixi) = Pxi(ixi) + 1 ! probabilities count

  if ( .not.adjustModeQ ) then

    if( mod(icycle,printStep)==0 ) then

      write(* ,"(/i9,9x,21(I12))") icycle,( nint(Pxi(i)), i= 0, MaxXi)
      write(1 ,"(/i9,9x,21(I12))") icycle,( nint(Pxi(i)), i= 0, MaxXi) ! unit 1 : monitoring file

      if ( minval(Pxi(0:maxXi)) > 0.001 ) then

        open(2, file='monitorMu'//IntToCha(partition)//'.txt', position='append')

        muTot = wxi(MaxXi) - wxi(0) - log(Pxi(MaxXi) / Pxi(0))
        write(* ,"(/'del mu(-) = ',9x,21(F12.5))") &
          & ( wxi(i)-wxi(i-1)-log(Pxi(i)/Pxi(i-1)), i= 1, MaxXi), temper*KtokJ*muTot
        write(1 ,"(/'del mu(-) = ',9x,21(F12.5))") &
          & ( wxi(i)-wxi(i-1)-log(Pxi(i)/Pxi(i-1)), i= 1, MaxXi), temper*KtokJ*muTot ! unit 1 : monitoring file
        write(2 ,"(I12, F12.5)") icycle, temper * KtokJ * muTot ! unit 1 : monitoring file

        close(2)

      end if

    end if

  end if

  end subroutine ExpMCSolid

  subroutine PrintFreeEnergy

  implicit none

  real(kind=r8) :: muTot

  muTot = 0.0_r8
  if (minval(Pxi) /= 0) muTot = wxi(MaxXi)-wxi(0)-log(Pxi(MaxXi)/Pxi(0))

  open(2, file='eemout.txt', position='append')

  write(2,"(I30, F30.6, F10.3)") partition, temper * KtokJ * muTot, get_ratio()

  close(2)

  end subroutine PrintFreeEnergy

  subroutine ExpandArray (j, fb) ! if j is larger than array, expand array to j

  logical             :: expQ
  logical, save       :: finiQ = .true., biniQ = .true.
  integer, intent(in) :: j
  integer, intent(in) :: fb
  integer, parameter  :: F = 1, B = 2 ! F = forword, B = backword
  integer             :: maxold, minold

  real(kind=r8), allocatable, dimension(:,:) :: temp

  expQ = .false.

  if (fb == F) then

    minold = fmin
    maxold = fmax

    if (j > fmax) then
      fmax = j
      expQ = .true.
    else if (j < fmin) then
      fmin = j
      expQ = .true.
    end if

    if (finiQ .and. expQ) then

      fmax = j
      fmin = j
      minold = j
      maxold = j
      finiQ = .false.

    end if

    if (expQ) then

      allocate(temp(0:MaxXi,minold:maxold))

      temp = nbeufor

      if (allocated(nbeufor)) deallocate(nbeufor)
      allocate(nbeufor(0:MaxXi,fmin:fmax))

      nbeufor = 0.0
      nbeufor(0:MaxXi,minold:maxold) = temp

    end if

  end if

  if (fb == B) then

    minold = bmin
    maxold = bmax

    if (j > bmax) then
      bmax = j
      expQ = .true.
    else if (j < bmin) then
      bmin = j
      expQ = .true.
    end if

    if (biniQ .and. expQ) then

      bmax = j
      bmin = j
      minold = j
      maxold = j
      biniQ = .false.

    end if

    if (expQ) then

      allocate(temp(0:MaxXi,minold:maxold))

      temp = nbeubak

      if (allocated(nbeubak)) deallocate(nbeubak)
      allocate(nbeubak(0:MaxXi,bmin:bmax))

      nbeubak = 0.0
      nbeubak(0:MaxXi,minold:maxold) = temp

    end if

  end if

  if (allocated(temp)) deallocate(temp)

  end subroutine ExpandArray

  subroutine SwapIndex (i)

  implicit none

  integer, intent(in) :: i
  real(kind=r8) :: dumX, dumY, dumZ

  dumX = rx(i,1)
  dumY = ry(i,1)
  dumZ = rz(i,1)

  rx(i,1) = rx(i,2)
  ry(i,1) = ry(i,2)
  rz(i,1) = rz(i,2)

  rx(i,2) = dumX
  ry(i,2) = dumY
  rz(i,2) = dumZ

  end subroutine SwapIndex

  subroutine SpecialMove (i, delE)

  use mod_einstein_core, only : EinEnergy

  implicit none

  integer, intent(in) :: i
  real(kind=r8), intent(out) :: delE

  real(kind=r8) :: oldE, newE, beu_

  call EinEnergy (i, rx(i,:), ry(i,:), rz(i,:), oldE)
  call SwapIndex (i)
  call EinEnergy (i, rx(i,:), ry(i,:), rz(i,:), newE)

  beu_ = beta * (1.0 - xi(ixi)) * (newE-oldE)
  if (newE < oldE) beu_ = 0
  if (exp(-beu_) > ranf()) then

    delE = newE-oldE

  else

    delE = 0.0
    call SwapIndex (i) ! back to origin

  end if

  end subroutine SpecialMove

  subroutine ExpSolDeallocate

  implicit none

  if (allocated(xi)) deallocate(xi, wxi, delw, delwOld, &
    &Pxi, nforxi, nbakxi, forxi, &
    &bakxi, oforbak)

  if (allocated(nbeufor) .and. allocated(nbeubak)) deallocate(nbeufor, nbeubak)

  write(*,"('ExpSol array is deallocated')")
  write(1,"('ExpSol array is deallocated')")

  end subroutine ExpSolDeallocate

  end module mod_expans_solid
