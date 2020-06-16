module mod_time_tool

    implicit none

    private

    integer, parameter :: r8 = selected_real_kind(15,300)
    real(kind=r8) :: time1, time2, rTime, rMin, rSec, remainTime ! r = remain

    public :: time1, time2
    public :: SecTohour, CalRemainTime, PrintRemainTime

    contains

    subroutine SecToHour (time_, unit_)

        implicit none

        real(kind=r8), intent(in) :: time_
        integer, optional, intent(in) :: unit_

        real(kind=r8) :: dtime
        integer :: hour_, min_, sec_

        dtime = time_

        hour_ = int(dtime / 3600.0_r8); dtime = dtime - real(3600*hour_)
        min_  = int(dtime / 60.0_r8  ); dtime = dtime - real(60*min_)
        sec_  = int(dtime)

        if (present(unit_)) then ;write(unit_,"(I5,'hour',I5,'min',I5,'sec')") hour_,min_,sec_
        else; write(*,"(I5,'hour',I5,'min',I5,'sec')") hour_,min_,sec_
        end if

    end subroutine SecToHour

    !********************

    subroutine CalRemainTime (curStep, maxStep)

        implicit none

        integer, intent(in) :: curStep, maxStep
        integer, parameter  :: sampleStep = 1000
        integer, save :: counter = 1
        integer       :: remainStep
        real(kind=r8), save :: countTime1, countTime2, runtime, runSpeed

        remainStep = maxStep - curStep

        counter = counter + 1

        call cpu_time(countTime2)
        if (counter == sampleStep) then

            runtime = countTime2 - countTime1
            runSpeed = dble(sampleStep) / runtime
            remainTime = dble(remainStep) / runSpeed
            call cpu_time(countTime1)
            counter = 1

        end if

    end subroutine CalRemainTime

    !********************

    subroutine PrintRemainTime (unit_)

        implicit none

        integer, intent(in), optional :: unit_

        real(kind=r8) :: dtime
        integer :: hour_, min_, sec_

        dtime = remainTime

        hour_ = int(dtime / 3600.0_r8); dtime = dtime - real(3600*hour_)
        min_  = int(dtime / 60.0_r8  ); dtime = dtime - real(60*min_)
        sec_  = int(dtime)

        if (present(unit_)) then
            write(unit_,"('remaining time : ',I5,'hour',I5,'min',I5,'sec')") hour_,min_,sec_
        else
            write(*,"('remaining time : ',I5,'hour',I5,'min',I5,'sec')") hour_,min_,sec_
        end if

    end subroutine PrintRemainTime

    !********************



    !********************

end module mod_time_tool

!============================================================
