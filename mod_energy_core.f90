!----------------------------------------------------------------------------------------------
!function EnergyLongRange (mat) result(dVlrc)         ! calculate LRC
!
!subroutine SumupIntra(eIntra_)                       ! calculate intra energy of all molecule
!subroutine EnergyIntra(xt, yt, zt, eIntra_)          ! calculate intra energy of molecule i
!subroutine CalculateEnergy(i, status, vdwE, chargeE) ! calculate inter energy of all molecule
!subroutine SpecificEnergy(i, status, vdwE, chargeE)  ! calculate inter energy of molecule i
!subroutine TotalEnergy(status, vdwE, chargeE)        ! calculate OLD or NEW energy
!subroutine CheckAccept (energy_, acceptQ)            ! accept
!
! non-using subroutines... will be deleted...
!subroutine Sumup (de, dec, dbox, dibox)
!subroutine SumupCurrent (rX_, rY_, rZ_, de, dec, dbox, dibox)
!subroutine Energy (i, xi_, yi_, zi_, de, dec, dbox, dibox)
!----------------------------------------------------------------------------------------------

module mod_energy_core

    use mod_univ_const, only : CHARGE_CONST, PI   ! to use univertial constant
    use mod_math_tool,  only : Det           ! to use erfc function
    use mod_system_information,   only : molN_e => molNum, rCut_e => rCut, &
                                         simulTypeENUM, EXP_SOLUTE


    use mod_molecule_information, only : sigij_e => sigij, epsij_e => epsij, &
                                         chaij_e => chaij, potenType_e => potenType, &
                                         maxC_e => nType, &
                                         nA_e => nAtoms, atomTypeNum, myType, &
                                         nBond, index1InBond, bondLength, index2InBond,&
                                         nBend, bond12InBend, bond23InBend, myBendType, &
                                         direc12InBend, direc23InBend, kBend, theta0, &
                                         nTorsion, bond12InTorsion, bond23InTorsion, bond34InTorsion, &
                                         bend123InTorsion, bend234InTorsion, myTorsionType, &
                                         direc12InTorsion, direc23InTorsion, direc34InTorsion, &
                                         v1, v2, v3, v4, neighborQ, scalingQ, scale14, myBondType, &
                                         direc123InTorsion, direc234InTorsion , &
                                         Group, localGroup, nLocalGroup, kBond, &
                                         A_FW, B_FW, C_FW, &
                                         A_DD, B_DD, C_DD, zeta_DD, R0_DD, D0_DD, &
                                         v5, v6, MIN_DISTANCE_DD, &
                                         C_theta_DD, cos_theta0, &
                                         d1, d2, d3, d4, d5, d6, &
                                         potenTypeENUM, OPLSAA, FW, DREIDING

    use mod_config_information,   only : rX_e => rX, rY_e => rY, rZ_e => rZ, &
                                         hBox_e => hBox, ihBox_e => ihBox, vol_e => volume, &   ! to use configuration information
                                         nGlobalGroup, grpX, grpY, grpZ, &
                                         grpXNew,  grpYNew,  grpZNew,  &
                                         grpiXOld, grpiYOld, grpiZOld, &
                                         grpiXNew, grpiYNew, grpiZNew, &
                                         myGroupType, globalGroupIndex, myMolIndex, &
                                         rxiOld, ryiOld, rziOld, rxiNew, ryiNew, rziNew, &
                                         newhBox, newihBox, rxNew, ryNew, rzNew, &
                                         soluteIndex

    use mod_movement_information, only : ShrinkSolute


    implicit none

    private

    integer, parameter :: r8 = selected_real_kind(15,300)
    !integer, parameter :: LJ = 1, EXP6 = 2
    real(kind=r8) :: eLrc, eVdw, eCha, eIntra ! long range, van der waals, charge energy

    public :: eLrc, eVdw, eCha, eIntra
    public :: EnergyLongRange, Sumup, SumupCurrent, Energy, CheckAccept, & !, &
              EnergyIntra, SumupIntra, CalculateEnergy
              !IPSSumup, IPSSumupCurrent, IPSEnergy

    contains

    !********************



!   -----------------------------------------------------------------------------------
    function EnergyLongRange (mat) result(dVlrc) ! have to be fixed
!   this function return LRC of system.
!   using simulation box and molecule information
!   -----------------------------------------------------------------------------------
        implicit none

        integer :: i, j, nSolvent
        real(kind=r8) :: dVlrc
        real(kind=r8) :: meanSig, meanEps
        real(kind=r8) :: Brc, Crc, repTerm, attTerm
        real(kind=r8) :: volumeLRC
        real(kind=r8), dimension(3,3), intent(in) :: mat ! simulation box matrix for volume calculation.

        dVlrc     = 0.0_r8
        volumeLRC = Det(mat,3)

!       ----------------------------------------------------------
        if (simulTypeENUM == EXP_SOLUTE) then ! for solute-insertion expanded ensemble simulation.

            nSolvent = molN_e - 1  ! neglect solute, treated as Long Range Correction.

        else

            nSolvent = molN_e      ! there is no solute in system.

        end if


!       ------------------------------------------------------------------------------------------
        select case (potenTypeENUM) ! select non-bonded potential type
!       ------------------------------------------------------------------------------------------

!           --------------------------------------------------------------------------------------
            case (OPLSAA) ! OPLS-AA
!           --------------------------------------------------------------------------------------
                do i=1, maxC_e

                   do j=1, maxC_e

                        meanSig = sigij_e(i,j)
                        meanEps = epsij_e(i,j)
                        dVlrc   = dVlrc + atomTypeNum(i) * atomTypeNum(j) * meanEps * &
                        meanSig ** 6 * (1.0_r8 / 3.0_r8 * (meanSig / rCut_e) ** 6 - 1.0_r8)

                    end do

                end do

                dVlrc = dVlrc * (8.0_r8 / 3.0_r8) * PI / &
                        volumeLRC / (rCut_e ** 3) * dble(nSolvent * nSolvent)

!           --------------------------------------------------------------------------------------
            case (FW) ! Flexible Williams
!           --------------------------------------------------------------------------------------
                do i = 1, maxC_e

                    do j = 1, maxC_e

                        Brc = B_FW(i,j) * rCut_e

                        repTerm = A_FW(i,j) * exp(-Brc) * (2.0_r8 + Brc * (2.0_r8 + Brc)) / (B_FW(i,j) ** 3)
                        attTerm = - C_FW(i,j) / 3.0_r8 / (Rcut_e ** 3)

                        dVLRC = dVLRC + atomTypeNum(i) * atomTypeNum(j) * (repTerm + attTerm)

                    end do

                end do

                dVLRC = dVLRC * 2.0_r8 * PI / volumeLRC * dble(nSolvent * nSolvent)

!           --------------------------------------------------------------------------------------
            case (DREIDING)
!           --------------------------------------------------------------------------------------

                do i = 1, maxC_e

                    do j = 1, maxC_e

                        Crc = C_DD(i,j) * rCut_e

                        repTerm = A_DD(i,j) * exp(-Crc) * (2.0_r8 + Crc * (2.0_r8 + Crc)) / (C_DD(i,j) ** 3)
                        attTerm = - B_DD(i,j) / 3.0_r8 / (Rcut_e ** 3)

                        dVLRC = dVLRC + atomTypeNum(i) * atomTypeNum(j) * (repTerm + attTerm)

                    end do

                end do

                dVLRC = dVLRC * 2.0_r8 * PI / volumeLRC * dble(nSolvent * nSolvent)

            case default

                stop 'invalid potential type'

        end select

    end function















!   ------------------------------------------------------------
    subroutine SumupIntra(eIntra_)
!   ------------------------------------------------------------
        implicit none

        real(kind=r8), intent(out) :: eIntra_
        real(kind=r8) :: de_

        integer :: mol
!   ------------------------------------------------------------

        eIntra_ = 0.0_r8
        do mol = 1, molN_e

            call EnergyIntra(rx_e(mol,:), ry_e(mol,:), rz_e(mol,:), de_, mol) ! intra energy of current
                                                      ! coordinate.
            eIntra_ = eIntra_ + de_                   ! optional variable 'mol' is used to check
                                                      ! the solute shrink
        end do

    end subroutine SumupIntra







!   ------------------------------------------------------------
    subroutine EnergyIntra(xt, yt, zt, eIntra_, index_)  ! soluteIndex is optional variable for
                                                         ! solute insertion expanded ensemble simulation.
!   ------------------------------------------------------------
        implicit none

        integer, parameter :: bondMax = 500, bendMax = 500, atomMax = 200

        real(kind=r8), dimension(0:nA_e), intent(inout) :: xt, yt, zt
        real(kind=r8), intent(out) :: eIntra_
        integer, optional, intent(in) :: index_

        integer :: i, j, ai, aj
        integer :: idx1, idx2, idx3, idx4, idx12, idx23, idx34, idx123, idx234 ! indices
        integer :: bond, bend, tors, bondType, bType, tType

        real(kind=r8), dimension(0:atomMax) :: dumxt, dumyt, dumzt
        real(kind=r8) :: arg, theta, th_th0
        real(kind=r8) :: eBond_, eBend_, eTors_, eNonBond_
        real(kind=r8) :: p12, p23, p34, p12p23, p23p34, p12p34
        real(kind=r8) :: cosphi, cosphi2, cos2phi, cos3phi, cos4phi, cos5phi, cos6phi, norm, coeff, t_cosphi
        real(kind=r8) :: rxij, ryij, rzij, rijsq, rij, s2, s6, ms, me, mc, eneV, eneQ, deV, deC
        real(kind=r8), dimension(bondMax,0:3) :: p     ! bond vector
        real(kind=r8), dimension(bendMax)     :: pp    ! bond vector inner product
        real(kind=r8) :: rij6, repTerm, attTerm, totalTerms
        real(kind=r8), parameter :: MIN_DISTANCE = 1.0_r8, INFINITY = 1.0e30
        logical :: soluteQ
!   ------------------------------------------------------------

        if (nBond > bondMax) stop 'bond max in Energyintra is too small'
        if (nBend > bendMax) stop 'bond max in Energyintra is too small'

        soluteQ = .false.
        if (simulTypeENUM == EXP_SOLUTE) then ! only for solute insertion ee simulation.

            if (present(index_)) then ! check index

                if (index_ == soluteIndex) soluteQ = .true.

            end if

        end if

        if (soluteQ) then

            dumxt(0:nA_e) = xt  ! save shrinked coordinate.
            dumyt(0:nA_e) = yt
            dumzt(0:nA_e) = zt
            call ShrinkSolute(1.0_r8, xt, yt, zt) ! reture to original molecule size.
                                                  ! savedRatio is not changed.
        end if

        !write(*,"('nBond = ', I5)") nBond
        eBond_ = 0.0_r8
        do bond = 1, nBond ! save bond vector & bond length

            idx1 = index1InBond(bond)
            idx2 = index2InBond(bond)
            !p(bond,0) = bondLength(myBondType(bond))
            p(bond,1) = xt(idx2) - xt(idx1)
            p(bond,2) = yt(idx2) - yt(idx1)
            p(bond,3) = zt(idx2) - zt(idx1)
            norm      = sqrt(p(bond,1)**2 + p(bond,2)**2 + p(bond,3)**2)
            !write(*,"('norm = ', F10.6)") norm
            p(bond,0) = norm

            bondType = myBondType(bond)

            if (potenTypeENUM == OPLSAA .or. potenTypeENUM == FW) then

                eBond_ = eBond_ + kBond(bondType) * (norm - bondLength(bondType)) ** 2

            else if (potenTypeENUM == DREIDING) then

                eBond_ = eBond_ + 0.5_r8 * kBond(bondType) * (norm - bondLength(bondType)) ** 2

            else

                stop 'invalid potential in intra energy > bond vibrational energy'

            end if

        end do ! bond-loop

        ! beding energy
        eBend_ = 0.0_r8
        do bend = 1, nBend

            idx12 = bond12InBend(bend)
            idx23 = bond23InBend(bend)
            ! save inner product of bonds in bend
            pp(bend) = (p(idx12,1) * p(idx23,1) + p(idx12,2) * p(idx23,2) + &
            p(idx12,3) * p(idx23,3)) * dble(direc12InBend(bend)*direc23InBend(bend))

            bType = myBendType(bend)
            arg   = -pp(bend) / p(idx12,0) / p(idx23,0) ! minus to get right value of theta
            if (arg >  1.0_r8) arg =  1.0_r8
            if (arg < -1.0_r8) arg = -1.0_r8

            if (potenTypeENUM == OPLSAA .or. potenTypeENUM == FW) then

                theta = acos(arg)
                th_th0 = theta - theta0(bType) ! in radian unit
                eBend_ = eBend_ + kBend(bType) * th_th0 * th_th0

            else if (potenTypeENUM == DREIDING) then

                th_th0 = arg - cos_theta0(bType) ! cos(th) - cos(th0)
                eBend_ = eBend_ + C_theta_DD(bType) * th_th0 * th_th0

            else

                stop 'invalid potential in intra energy > bend vibrational energy'

            end if

        end do ! bend-loop

        ! torsion energy
        eTors_ = 0.0_r8
        do tors = 1, nTorsion

            idx12  = bond12InTorsion(tors)
            idx23  = bond23InTorsion(tors)
            idx34  = bond34InTorsion(tors)
            idx123 = bend123InTorsion(tors)
            idx234 = bend234InTorsion(tors)

            p12 = p(idx12,0)
            p23 = p(idx23,0)
            p34 = p(idx34,0)
            p12p23 = pp(idx123)
            p23p34 = pp(idx234)
            p12p34 = (p(idx12,1) * p(idx34,1) + p(idx12,2) * p(idx34,2) + &
            p(idx12,3) * p(idx34,3)) * dble(direc12InTorsion(tors) * direc34InTorsion(tors))

            cosphi = (p12p23 * p23p34 - p12p34 * p23 * p23) / sqrt( &
                     (p12 * p12 * p23 * p23 - p12p23 * p12p23) * &
                     (p23 * p23 * p34 * p34 - p23p34 * p23p34))

            !cosphi2 = cosphi * cosphi
            t_cosphi = 2.0_r8 * cosphi
            cos2phi = t_cosphi * cosphi - 1.0_r8
            cos3phi = t_cosphi * cos2phi - cosphi
            !cos3phi = cosphi * (4.0_r8 * cosphi2 - 3.0_r8)

            !write(*,"(I3, A, 4F10.4)") tors, '  cos(phi) = ', cosphi, cos2phi, cos3phi, cos4phi

            tType = myTorsionType(tors)

            if (potenTypeENUM == OPLSAA) then

                cos4phi = t_cosphi * cos3phi - cos2phi
                !cos4phi = 8.0_r8 * cosphi2 * (cosphi2 - 1.0_r8) + 1.0_r8

                eTors_ = eTors_ + 0.5_r8 * (v1(tType) * (1.0_r8 + cosphi)  + &
                                            v2(tType) * (1.0_r8 - cos2phi) + &
                                            v3(tType) * (1.0_r8 + cos3phi) + &
                                            v4(tType) * (1.0_r8 - cos4phi))

            else if (potenTypeENUM == FW) then

                eTors_ = eTors_ + (v1(tType) * (1.0_r8 + cosphi)  + &
                                   v2(tType) * (1.0_r8 + cos2phi) + &
                                   v3(tType) * (1.0_r8 + cos3phi))

            else if (potenTypeENUM == DREIDING) then

                !cos5phi = cosphi  * (2.0_r8 * (cos4phi - cos2phi) + 1.0_r8)
                !cos6phi = cos2phi * (2.0_r8 *  cos4phi - 1.0_r8)
                cos5phi = t_cosphi * cos4phi - cos3phi
                cos6phi = t_cosphi * cos5phi - cos4phi

                eTors_ = eTors_ + 0.5_r8 * (v1(tType) * (1.0_r8 - d1(tType) * cosphi)  + &
                                            v2(tType) * (1.0_r8 - d2(tType) * cos2phi) + &
                                            v3(tType) * (1.0_r8 - d3(tType) * cos3phi) + &
                                            v4(tType) * (1.0_r8 - d4(tType) * cos4phi) + &
                                            v5(tType) * (1.0_r8 - d5(tType) * cos5phi) + &
                                            v6(tType) * (1.0_r8 - d6(tType) * cos6phi))
            end if

        end do ! torsion-loop
        !write(*,*) 'TORSION@@@@@@@@@@@@@@@@@@@@@@', eTors_
        ! non-bonded energy
        eneV = 0.0_r8
        eneQ = 0.0_r8
        eNonBond_ = 0.0_r8
        do i = 1, nA_e - 1

            ai = myType(i)

            do j = i + 1, nA_e

                if (neighborQ(i,j)) cycle

                if (potenTypeENUM == OPLSAA) then

                    aj = myType(j)

                    rxij = xt(j) - xt(i)
                    ryij = yt(j) - yt(i)
                    rzij = zt(j) - zt(i)

                    rijsq = rxij * rxij + ryij * ryij + rzij * rzij
                    rij   = sqrt(rijsq)

                    me    = epsij_e(ai,aj)
                    ms    = sigij_e(ai,aj)
                    mc    = chaij_e(ai,aj)

                    s2    = ms * ms / rijsq
                    s6    = s2 * s2 * s2

                    deV   = me * s6 * (s6 - 1.0_r8)
                    deC   = mc / rij

                    if (scalingQ(i,j)) then

                        deV = scale14 * deV
                        deC = scale14 * deC

                    end if

                    eneV  = eneV + deV
                    eneQ  = eneQ + deC

                else if (potenTypeENUM == FW) then

                    if (scalingQ(i,j)) cycle

                    aj = myType(j)

                    rxij = xt(j) - xt(i)
                    ryij = yt(j) - yt(i)
                    rzij = zt(j) - zt(i)

                    rijsq = rxij * rxij + ryij * ryij + rzij * rzij
                    rij   = sqrt(rijsq)
                    rij6  = rijsq * rijsq * rijsq

                    repTerm =  A_FW(ai,aj) * exp(-B_FW(ai,aj) * rij)
                    attTerm = -C_FW(ai,aj) / rij6
                    totalTerms = repTerm + attTerm

                    if (rij < MIN_DISTANCE) then

                        if (totalTerms < 0.0) then

                            eIntra_ = INFINITY
                            return

                        end if

                    end if

                    eneV = eneV + totalTerms

                else if (potenTypeENUM == DREIDING) then

                    aj = myType(j)

                    rxij = xt(j) - xt(i)
                    ryij = yt(j) - yt(i)
                    rzij = zt(j) - zt(i)

                    rijsq = rxij * rxij + ryij * ryij + rzij * rzij
                    rij   = sqrt(rijsq)
                    rij6  = rijsq * rijsq * rijsq

                    repTerm =  A_DD(ai,aj) * exp(-C_DD(ai,aj) * rij)
                    attTerm = -B_DD(ai,aj) / rij6
                    mc      =  chaij_e(ai,aj)
                    totalTerms = repTerm + attTerm

                    if (rij < MIN_DISTANCE_DD(ai,aj)) then

                        eIntra_ = INFINITY
                        return

                    end if

                    eneV = eneV + totalTerms
                    eneQ = eneQ + mc / rij

                else

                    stop 'invalid potential type in energy intra > nonbond'

                end if

            end do ! j-loop

        end do ! i-loop

        ! final treatment for efficient calculation.
        if (potenTypeENUM == OPLSAA) then

            eneV  = 4.0_r8 * eneV
            eneQ  = CHARGE_CONST * eneQ

        else if (potenTypeENUM == DREIDING) then

            eneQ = CHARGE_CONST * eneQ

        else if (potenTypeENUM == FW) then

        else

            stop 'invalid potential type in energy intra > final treatment'

        end if

        eNonBond_ = eneV + eneQ

        eIntra_ = eBond_ + eBend_ + eTors_ + eNonBond_

        if (soluteQ) then

            xt = dumxt(0:nA_e)   ! save shrinked coordinate.
            yt = dumyt(0:nA_e)   ! savedRatio is not changed.
            zt = dumzt(0:nA_e)

        end if

    end subroutine EnergyIntra









    !********************

    subroutine CheckAccept (energy_, acceptQ)
    ! check aceeptance using normal Metropolis scheme (non-biased).
        use mod_random_tool, only : Ranf

        implicit none

        real(kind=r8), intent(in) :: energy_
        logical, intent(out) :: acceptQ
        real(kind=r8) :: testE

        acceptQ = .false.

        testE = energy_

        if (testE < 75.0_r8) then

            if (testE < 0.0_r8) then

                acceptQ = .true.

            else if (exp(-testE) > Ranf()) then

                acceptQ = .true.

            end if

        end if

    end subroutine CheckAccept

    !********************







!   --------------------------------------------------------------------------------
    subroutine CalculateEnergy(i, status, vdwE, chargeE, lambda)

!   --------------------------------------------------------------------------------
        implicit none

        integer, intent(in) :: i
        character(len=3), intent(in) :: status
        real(kind=r8), intent(out) :: vdwE, chargeE
        real(kind=r8), intent(in), optional :: lambda ! solute calculation.

        logical :: allQ, exp_soluteQ
!   --------------------------------------------------------------------------------

        if (i == 0) allQ = .true.
        if (i /= 0) allQ = .false.
        exp_soluteQ = .false.
        if (present(lambda)) exp_soluteQ = .true.

        if (allQ) then  ! all energy in system calculation

            if (exp_soluteQ)      call TotalEnergy(status, vdwE, chargeE, lambda) ! for with solute calculation
            if (.not.exp_soluteQ) call TotalEnergy(status, vdwE, chargeE)

        else if (.not.allQ) then

            if (exp_soluteQ)      call SpecificEnergy(i, status, vdwE, chargeE, lambda) ! for with solute calculation
            if (.not.exp_soluteQ) call SpecificEnergy(i, status, vdwE, chargeE)

        else

            stop 'critical error in energy routine' ! never occurs

        end if

    end subroutine CalculateEnergy





!   ----------------------------------------------------------------------------------------------
    subroutine SpecificEnergy(i, status, vdwE, chargeE, lambda)
!   ----------------------------------------------------------------------------------------------
        use mod_neighbor_list

        implicit none

        integer, intent(in) :: i
        character(len=3), intent(in) :: status
        real(kind=r8), intent(out) :: vdwE, chargeE
        real(kind=r8), intent(in), optional :: lambda

        real(kind=r8), parameter :: ZERO = 0.0_r8, ONE = 1.0_r8, HALF = 0.5_r8, FOUR = 4.0_r8

        integer :: mg, ig, nl, jg, maxList, niMem, njMem, ma, na, ia, ja, iat, jat, j

        real(kind=r8) :: rcutsq
        integer, dimension(:), pointer :: list, imem, jmem
        real(kind=r8), dimension(:), pointer :: rxi, ryi, rzi, grpix, grpiy, grpiz
        real(kind=r8) :: gxij, gyij, gzij, pijsq, dxij, dyij, dzij, pxij, pyij, pzij, &
                         sxij, syij, szij, rxij, ryij, rzij, rijsq, rij
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33
        real(kind=r8) :: ms, me, mc, sr2, sr6
        real(kind=r8) :: rij6, repTerm, attTerm, totalTerms
        real(kind=r8), parameter :: MIN_DISTANCE = 1.0_r8, INFINITY = 1.0e30

        ! solute parameters and variables update Date : 2015-03-04
        real(kind=r8) :: TWO_SIX = 2.0_r8 ** (1.0_r8 / 6.0_r8)
        real(kind=r8) :: zeta, zeta4, lam, lam3, minDistance, ms_lam3
        logical :: soluteQ = .false., lambdaQ = .false.
!   ----------------------------------------------------------------------------------------------

        vdwE    = ZERO
        chargeE = ZERO
        rcutsq  = rcut_e * rcut_e

        soluteQ = .false.
        lambdaQ = .false.
        if (present(lambda)) then

            lambdaQ = .true.

        end if

        if (lambdaQ) then

            zeta = lambda - ONE
            if (zeta < ZERO) zeta = ZERO ! limit
            zeta4 = 0.25_r8 * zeta

            lam = lambda
            if (lam > ONE) lam = ONE ! limit
            lam3 = lam ** (1.0_r8 / 3.0_r8)

        end if


        if (status == 'OLD') then

            rxi => rxiOld
            ryi => ryiOld
            rzi => rziOld

            grpix => grpiXOld
            grpiy => grpiYOld
            grpiz => grpiZOld
            !write(*,*) "i'm old"

        end if

        if (status == 'NEW') then

            rxi => rxiNew
            ryi => ryiNew
            rzi => rziNew

            grpix => grpiXNew
            grpiy => grpiYNew
            grpiz => grpiZNew
            !write(*,*) "i'm new"

        end if

        b11  =  hbox_e(1,1)
        b12  =  hbox_e(1,2)
        b13  =  hbox_e(1,3)
        b22  =  hbox_e(2,2)
        b23  =  hbox_e(2,3)
        b33  =  hbox_e(3,3)

        ib11 = ihbox_e(1,1)
        ib12 = ihbox_e(1,2)
        ib13 = ihbox_e(1,3)
        ib22 = ihbox_e(2,2)
        ib23 = ihbox_e(2,3)
        ib33 = ihbox_e(3,3)

        do mg = 1, nLocalGroup

            ig = globalGroupIndex(i,mg)
            list => nlist(ig)%n
            maxList = list(0)

            imem => localGroup(mg)%member
            niMem = imem(0)

            do nl = 1, maxList

                jg = list(nl)

                gxij = grpiX(mg) - grpX(jg)
                gyij = grpiY(mg) - grpY(jg)
                gzij = grpiZ(mg) - grpZ(jg)

                sxij = ib11 * gxij + ib12 * gyij + ib13 * gzij
                syij =               ib22 * gyij + ib23 * gzij
                szij =                             ib33 * gzij

                if (sxij >  HALF) sxij = sxij - ONE
                if (sxij < -HALF) sxij = sxij + ONE

                if (syij >  HALF) syij = syij - ONE
                if (syij < -HALF) syij = syij + ONE

                if (szij >  HALF) szij = szij - ONE
                if (szij < -HALF) szij = szij + ONE

                pxij = b11 * sxij + b12 * syij + b13 * szij
                pyij =              b22 * syij + b23 * szij
                pzij =                           b33 * szij

                pijsq = pxij * pxij + pyij * pyij + pzij * pzij

                if (pijsq > rcutsq) cycle

                j    = myMolIndex(jg)
                jmem => localGroup(myGroupType(jg))%member
                njMem = jmem(0)
                dxij = gxij - pxij
                dyij = gyij - pyij
                dzij = gzij - pzij

                if (lambdaQ) then ! only for exp_solute simulation

                soluteQ = .false.  ! solute test
                    if (i == soluteIndex .or. j == soluteIndex) soluteQ = .true.

                end if

                do ma = 1, niMem

                    ia  = imem(ma)
                    iat = myType(ia)

                    do na = 1, njMem

                        ja  = jmem(na)
                        jat = myType(ja)

                        rxij = rxi(ia) - rx_e(j,ja) - dxij
                        ryij = ryi(ia) - ry_e(j,ja) - dyij
                        rzij = rzi(ia) - rz_e(j,ja) - dzij

                        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
                        rij   = sqrt(rijsq)

                        if (potenTypeENUM == OPLSAA) then

                            ms = sigij_e(iat,jat)
                            me = epsij_e(iat,jat)
                            mc = chaij_e(iat,jat)

                            if (soluteQ) then

                                ms_lam3 = ms * lam3    ! scaled sigma
                                minDistance = TWO_SIX * ms_lam3

                                sr2 = ms_lam3 * ms_lam3 / rijsq
                                sr6 = sr2 * sr2 * sr2

                                if (rij < minDistance) then

                                    vdwE = vdwE + me * (sr6 * (sr6 - ONE) + 0.25_r8 - zeta4)

                                else ! rij > minDistance

                                    vdwE = vdwE + zeta * me * sr6 * (sr6 - ONE)

                                end if

                                chargeE = chargeE + zeta * mc / rij

                            else

                                sr2 = ms * ms / rijsq
                                sr6 = sr2 * sr2 * sr2

                                vdwE    = vdwE + me * sr6 * (sr6 - ONE)
                                chargeE = chargeE + mc / rij

                            end if

                        else if (potenTypeENUM == FW) then

                            rij6 = rijsq * rijsq * rijsq

                            repTerm = A_FW(iat,jat) * exp(-B_FW(iat,jat) * rij)
                            attTerm = -C_FW(iat,jat) / rij6

                            totalTerms = repTerm + attTerm
                            if (rij < MIN_DISTANCE) then

                                if (totalTErms < ZERO) then

                                    VDWE = INFINITY
                                    return

                                end if

                            end if

                            VDWE = VDWE + totalTerms

                        else if (potenTypeENUM == DREIDING) then

                            if (rij < MIN_DISTANCE_DD(iat,jat)) then

                                VDWE    = INFINITY
                                chargeE = INFINITY
                                return

                            end if

                            rij6 = rijsq * rijsq * rijsq

                            repTerm =  A_DD(iat,jat) * exp(-C_DD(iat,jat) * rij)
                            attTerm = -B_DD(iat,jat) / rij6
                            totalTerms = repTerm + attTerm

                            mc = chaij_e(iat,jat)

                            VDWE    = VDWE + totalTerms
                            chargeE = chargeE + mc / rij

                        else

                            stop 'invalid potentype in energy core > inner loop'

                        end if

                    end do ! na-loop

                end do ! ma-loop

            end do ! nl-loop

        end do ! mg-loop

        if (potenTypeENUM == OPLSAA) then

            vdwE    = FOUR * vdwE
            chargeE = CHARGE_CONST * chargeE

        else if (potenTypeENUM == DREIDING) then

            chargeE = CHARGE_CONST * chargeE

        else if (potenTypeENUM == FW) then

        else

            stop 'invalid potentype in energy core > final treatment'

        end if

    end subroutine SpecificEnergy





!   --------------------------------------------------------------------------------
    subroutine TotalEnergy(status, vdwE, chargeE, lambda)
!   --------------------------------------------------------------------------------
        use mod_neighbor_list

        implicit none

        character(len=3) :: status ! 'OLD' or 'NEW'
        real(kind=r8), intent(out) :: vdwE, chargeE
          real(kind=r8), intent(in), optional :: lambda

        real(kind=r8), parameter :: ZERO = 0.0_r8, ONE = 1.0_r8, HALF = 0.5_r8, FOUR = 4.0_r8

        real(kind=r8), dimension(:),   pointer :: gx, gy, gz
        real(kind=r8), dimension(:,:), pointer :: rx, ry, rz, box, ibox

        real(kind=r8) :: gxij, gyij, gzij, pijsq, dxij, dyij, dzij, pxij, pyij, pzij, sxij, syij, &
                         szij, rxij, ryij, rzij, rijsq, rij
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33
        real(kind=r8) :: ms, me, mc, sr2, sr6, rcutsq
        real(kind=r8) :: rij6, repTerm, attTerm, totalTerms
        real(kind=r8), parameter :: MIN_DISTANCE = 1.0_r8, INFINITY = 1.0e30

        integer :: ig, jg, nl, i, j, ma, na, ia, ja, iat, jat, mg, ng    ! dummy index
        integer :: maxList, niMem, njMem
        integer, dimension(:), pointer :: list, imem, jmem

        ! solute parameters and variables update Date : 2015-03-04
        real(kind=r8) :: TWO_SIX = 2.0_r8 ** (1.0_r8 / 6.0_r8)
        real(kind=r8) :: zeta, zeta4, lam, lam3, minDistance, ms_lam3
        logical :: soluteQ = .false., lambdaQ = .false.
!   --------------------------------------------------------------------------------

        vdwE    = ZERO
        chargeE = ZERO
        rcutsq  = rcut_e * rcut_e

        soluteQ = .false.
        lambdaQ = .false.
        if (present(lambda)) then

            lambdaQ = .true.

        end if

        if (lambdaQ) then

            zeta = lambda - ONE
            if (zeta < ZERO) zeta = ZERO ! limit
            zeta4 = 0.25_r8 * zeta

            lam = lambda
            if (lam > ONE) lam = ONE ! limit
            lam3 = lam ** (1.0_r8 / 3.0_r8)

        end if

        if (status == 'OLD') then

            gx => grpX
            gy => grpY
            gz => grpZ

            rx => rx_e
            ry => ry_e
            rz => rz_e

            box  => hBox_e
            ibox => ihBox_e

        end if

        if (status == 'NEW') then

            gx => grpXNew
            gy => grpYNew
            gz => grpZNew

            rx => rxNew
            ry => ryNew
            rz => rzNew

            box  => newhBox
            ibox => newihBox

        end if

        b11  =  box(1,1)
        b12  =  box(1,2)
        b13  =  box(1,3)
        b22  =  box(2,2)
        b23  =  box(2,3)
        b33  =  box(3,3)

        ib11 = ibox(1,1)
        ib12 = ibox(1,2)
        ib13 = ibox(1,3)
        ib22 = ibox(2,2)
        ib23 = ibox(2,3)
        ib33 = ibox(3,3)

        do ig = 1, nGlobalGroup

            list => nlist(ig)%n
            maxList = list(0)
            i = myMolIndex(ig)

            imem => localGroup(myGroupType(ig))%member
            niMem = imem(0)

            do nl = 1, maxList

                jg = list(nl)
                if (ig >= jg) cycle ! i < j condition

                gxij = gx(ig) - gx(jg)
                gyij = gy(ig) - gy(jg)
                gzij = gz(ig) - gz(jg)

                sxij = ib11 * gxij + ib12 * gyij + ib13 * gzij
                syij =               ib22 * gyij + ib23 * gzij
                szij =                             ib33 * gzij

                if (sxij >  HALF) sxij = sxij - ONE
                if (sxij < -HALF) sxij = sxij + ONE

                if (syij >  HALF) syij = syij - ONE
                if (syij < -HALF) syij = syij + ONE

                if (szij >  HALF) szij = szij - ONE
                if (szij < -HALF) szij = szij + ONE

                pxij = b11 * sxij + b12 * syij + b13 * szij
                pyij =              b22 * syij + b23 * szij
                pzij =                           b33 * szij

                pijsq = pxij * pxij + pyij * pyij + pzij * pzij ! minimum image distance

                if (pijsq > rcutsq) cycle

                j = myMolIndex(jg)
                jmem => localGroup(myGroupType(jg))%member
                njMem = jmem(0)

                dxij = gxij - pxij
                dyij = gyij - pyij
                dzij = gzij - pzij

                if (lambdaQ) then ! only for exp_solute simulation

                    soluteQ = .false.  ! solute test
                    if (i == soluteIndex .or. j == soluteIndex) soluteQ = .true.

                end if

                do ma = 1, niMem

                    ia  = imem(ma)
                    iat = myType(ia)

                    do na = 1, njMem

                        ja  = jmem(na)
                        jat = myType(ja)

                        rxij = rx(i,ia) - rx(j,ja) - dxij
                        ryij = ry(i,ia) - ry(j,ja) - dyij
                        rzij = rz(i,ia) - rz(j,ja) - dzij

                        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
                        rij   = sqrt(rijsq)

                        if (potenTypeENUM == OPLSAA) then

                            ms = sigij_e(iat,jat)
                            me = epsij_e(iat,jat)
                            mc = chaij_e(iat,jat)

                            if (soluteQ) then ! for solute-solvect interaction.

                                ms_lam3 = ms * lam3    ! scaled sigma
                                minDistance = TWO_SIX * ms_lam3

                                sr2 = ms_lam3 * ms_lam3 / rijsq
                                sr6 = sr2 * sr2 * sr2

                                if (rij < minDistance) then

                                    vdwE = vdwE + me * (sr6 * (sr6 - ONE) + 0.25_r8 - zeta4)

                                else ! rij > minDistance

                                    vdwE = vdwE + zeta * me * sr6 * (sr6 - ONE)

                                end if

                                chargeE = chargeE + zeta * mc / rij

                            else ! for solvent

                                sr2 = ms * ms / rijsq
                                sr6 = sr2 * sr2 * sr2

                                vdwE    = vdwE + me * sr6 * (sr6 - ONE)
                                chargeE = chargeE + mc / rij

                            end if

                        else if (potenTypeENUM == FW) then

                            rij6 = rijsq * rijsq * rijsq
                            repTerm = A_FW(iat,jat) * exp(-B_FW(iat,jat) * rij)
                            attTerm = -C_FW(iat,jat) / rij6

                            totalTerms = repTerm + attTerm
                            if (rij < MIN_DISTANCE) then

                                if (totalTerms < ZERO) then

                                    VDWE = INFINITY
                                    return

                                end if

                            end if

                            VDWE = VDWE + totalTerms

                        else if (potenTypeENUM == DREIDING) then

                            if (rij < MIN_DISTANCE_DD(iat,jat)) then

                                VDWE    = INFINITY
                                chargeE = INFINITY
                                return

                            end if

                            rij6 = rijsq * rijsq * rijsq

                            repTerm =  A_DD(iat,jat) * exp(-C_DD(iat,jat) * rij)
                            attTerm = -B_DD(iat,jat) / rij6
                            totalTerms = repTerm + attTerm

                            mc = chaij_e(iat,jat)

                            VDWE    = VDWE + totalTerms
                            chargeE = chargeE + mc / rij

                        else

                            stop 'invalid potentype in energy core > inner loop'

                        end if

                    end do ! na-loop

                end do ! ma-loop

            end do ! nl-loop

        end do ! ig-loop

        if (potenTypeENUM == OPLSAA) then

            vdwE    = FOUR * vdwE
            chargeE = CHARGE_CONST * chargeE

        else if (potenTypeENUM == DREIDING) then

            chargeE = CHARGE_CONST * chargeE

        else if (potenTypeENUM == FW) then

        else

            stop 'invalid potentype in energy core > final treatment'

        end if

    end subroutine TotalEnergy






!   ---------------------------------------------------------------------------------
!   ---------------------------------------------------------------------------------
!   ---------------------------------------------------------------------------------
!   ---------------------------------------------------------------------------------
!   ---------------------------------------------------------------------------------
!   ---------------------------------------------------------------------------------















!   Dead routines...






   !********************

    subroutine Sumup (de, dec, dbox, dibox)

        implicit none

        integer :: i, j, ii, jj, ai, aj

        real(kind=r8), intent(out)    :: de, dec ! dummy energy, charge energy
        real(kind=r8), dimension(3,3), intent(in) :: dbox, dibox

        real(kind=r8) :: rcutsq
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33
        real(kind=r8) :: sr2, sr6, rxij, ryij, rzij, rij, rijsq, dxij, dyij, dzij, srxij, sryij, srzij, &
                         rrxij, rryij, rrzij
        real(kind=r8) :: ms, me, mc ! mean sigma, epsilon, charge

        de  = 0.0_r8
        dec = 0.0_r8
        rcutsq = rCut_e**2

        b11 = dbox(1,1)
        b12 = dbox(1,2)
        b13 = dbox(1,3)
        b22 = dbox(2,2)
        b23 = dbox(2,3)
        b33 = dbox(3,3)

        ib11 = dibox(1,1)
        ib12 = dibox(1,2)
        ib13 = dibox(1,3)
        ib22 = dibox(2,2)
        ib23 = dibox(2,3)
        ib33 = dibox(3,3)

        do i = 1, molN_e - 1

            do j = i+1, molN_e

                rrxij = rX_e(i,0) - rX_e(j,0)
                rryij = rY_e(i,0) - rY_e(j,0)
                rrzij = rZ_e(i,0) - rZ_e(j,0)

                dxij = ib11*rrxij + ib12*rryij + ib13*rrzij
                dyij =            + ib22*rryij + ib23*rrzij
                dzij =                         + ib33*rrzij

                if (dxij >  0.5_r8) dxij = dxij - 1.0_r8
                if (dxij < -0.5_r8) dxij = dxij + 1.0_r8

                if (dyij >  0.5_r8) dyij = dyij - 1.0_r8
                if (dyij < -0.5_r8) dyij = dyij + 1.0_r8

                if (dzij >  0.5_r8) dzij = dzij - 1.0_r8
                if (dzij < -0.5_r8) dzij = dzij + 1.0_r8

                srxij = b11*dxij + b12*dyij + b13*dzij
                sryij =            b22*dyij + b23*dzij
                srzij =                       b33*dzij

                rijsq = srxij*srxij + sryij*sryij + srzij*srzij

                if (rijsq < rcutsq) then

                    rrxij = rrxij - srxij
                    rryij = rryij - sryij
                    rrzij = rrzij - srzij

                    do ii = 1, nA_e; do jj = 1, nA_e

                    rxij = rX_e(i,ii) - rX_e(j,jj) - rrxij
                    ryij = rY_e(i,ii) - rY_e(j,jj) - rryij
                    rzij = rZ_e(i,ii) - rZ_e(j,jj) - rrzij

                    rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                    rij   = rijsq**0.5_r8

                    ai = myType(ii)
                    aj = myType(jj)

                    ms = sigij_e(ai,aj)
                    me = epsij_e(ai,aj)
                    mc = chaij_e(ai,aj)

                    sr2 = ms*ms / rijsq
                    sr6 = sr2*sr2*sr2

                    de  = de  + me * sr6 * (sr6-1.0_r8)
                    dec = dec + mc / rij

                    end do; end do;

                end if

            end do

        end do

    de  = 4.0_r8 * de
    dec = CHARGE_CONST * dec

    end subroutine Sumup

    !********************

    subroutine SumupCurrent (rX_, rY_, rZ_, de, dec, dbox, dibox)

        implicit none

        integer :: i, j, ii, jj, ai, aj

        real(kind=r8), intent(out)    :: de, dec ! dummy energy, charge energy
        real(kind=r8), dimension(3,3), intent(in) :: dbox, dibox
        real(kind=r8), dimension(molN_e,0:nA_e), intent(in) :: rX_, rY_, rZ_

        real(kind=r8) :: rcutsq
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33
        real(kind=r8) :: sr2, sr6, rxij, ryij, rzij, rij, rijsq, dxij, dyij, dzij, srxij, sryij, srzij, &
                         rrxij, rryij, rrzij
        real(kind=r8) :: ms, me, mc ! mean sigma, epsilon, charge

        de  = 0.0_r8
        dec = 0.0_r8
        rcutsq = rCut_e**2

        b11 = dbox(1,1)
        b12 = dbox(1,2)
        b13 = dbox(1,3)
        b22 = dbox(2,2)
        b23 = dbox(2,3)
        b33 = dbox(3,3)

        ib11 = dibox(1,1)
        ib12 = dibox(1,2)
        ib13 = dibox(1,3)
        ib22 = dibox(2,2)
        ib23 = dibox(2,3)
        ib33 = dibox(3,3)

        do i = 1, molN_e - 1

            do j = i+1, molN_e

                rrxij = rX_(i,0) - rX_(j,0)
                rryij = rY_(i,0) - rY_(j,0)
                rrzij = rZ_(i,0) - rZ_(j,0)

                dxij = ib11*rrxij + ib12*rryij + ib13*rrzij
                dyij =            + ib22*rryij + ib23*rrzij
                dzij =                         + ib33*rrzij

                if (dxij >  0.5_r8) dxij = dxij - 1.0_r8
                if (dxij < -0.5_r8) dxij = dxij + 1.0_r8

                if (dyij >  0.5_r8) dyij = dyij - 1.0_r8
                if (dyij < -0.5_r8) dyij = dyij + 1.0_r8

                if (dzij >  0.5_r8) dzij = dzij - 1.0_r8
                if (dzij < -0.5_r8) dzij = dzij + 1.0_r8

                srxij = b11*dxij + b12*dyij + b13*dzij
                sryij =            b22*dyij + b23*dzij
                srzij =                       b33*dzij

                rijsq = srxij*srxij + sryij*sryij + srzij*srzij

                if (rijsq < rcutsq) then

                    rrxij = rrxij - srxij
                    rryij = rryij - sryij
                    rrzij = rrzij - srzij

                    do ii = 1, nA_e; do jj = 1, nA_e

                    rxij = rX_(i,ii) - rX_(j,jj) - rrxij
                    ryij = rY_(i,ii) - rY_(j,jj) - rryij
                    rzij = rZ_(i,ii) - rZ_(j,jj) - rrzij

                    rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                    rij   = rijsq**0.5_r8

                    ai = myType(ii)
                    aj = myType(jj)

                    ms = sigij_e(ai,aj)
                    me = epsij_e(ai,aj)
                    mc = chaij_e(ai,aj)

                    sr2 = ms*ms / rijsq
                    sr6 = sr2*sr2*sr2

                    de  = de  + me * sr6 * (sr6-1.0_r8)
                    dec = dec + mc / rij

                    end do; end do;

                end if

            end do

        end do

    de  = 4.0_r8 * de
    dec = CHARGE_CONST * dec

    end subroutine SumupCurrent ! ¹Ì¿Ï¼º

    !********************

    subroutine Energy (i, xi_, yi_, zi_, de, dec, dbox, dibox)

        implicit none

        integer :: j, ii, jj, ai, aj

        integer, intent(in) :: i
        real(kind=r8), intent(out)    :: de, dec ! dummy energy, charge energy
        real(kind=r8), dimension(3,3), intent(in) :: dbox, dibox
        real(kind=r8), dimension(0:nA_e), intent(in) :: xi_, yi_, zi_

        real(kind=r8) :: rcutsq
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33
        real(kind=r8) :: sr2, sr6, rxij, ryij, rzij, rij, rijsq, dxij, dyij, dzij, srxij, sryij, srzij, &
                         rrxij, rryij, rrzij
        real(kind=r8) :: ms, me, mc ! mean sigma, epsilon, charge

        de  = 0.0_r8
        dec = 0.0_r8
        rcutsq = rCut_e**2

        b11 = dbox(1,1)
        b12 = dbox(1,2)
        b13 = dbox(1,3)
        b22 = dbox(2,2)
        b23 = dbox(2,3)
        b33 = dbox(3,3)

        ib11 = dibox(1,1)
        ib12 = dibox(1,2)
        ib13 = dibox(1,3)
        ib22 = dibox(2,2)
        ib23 = dibox(2,3)
        ib33 = dibox(3,3)

        do j = 1, molN_e

        if (i /= j) then ! outer if

            rrxij = xi_(0) - rX_e(j,0)
            rryij = yi_(0) - rY_e(j,0)
            rrzij = zi_(0) - rZ_e(j,0)

            dxij = ib11*rrxij + ib12*rryij + ib13*rrzij
            dyij =            + ib22*rryij + ib23*rrzij
            dzij =                         + ib33*rrzij

            if (dxij >  0.5_r8) dxij = dxij - 1.0_r8
            if (dxij < -0.5_r8) dxij = dxij + 1.0_r8

            if (dyij >  0.5_r8) dyij = dyij - 1.0_r8
            if (dyij < -0.5_r8) dyij = dyij + 1.0_r8

            if (dzij >  0.5_r8) dzij = dzij - 1.0_r8
            if (dzij < -0.5_r8) dzij = dzij + 1.0_r8

            srxij = b11*dxij + b12*dyij + b13*dzij
            sryij =            b22*dyij + b23*dzij
            srzij =                       b33*dzij

            rijsq = srxij*srxij + sryij*sryij + srzij*srzij

            if (rijsq < rcutsq) then

                rrxij = rrxij - srxij
                rryij = rryij - sryij
                rrzij = rrzij - srzij

                do ii = 1, nA_e; do jj = 1, nA_e

                rxij = xi_(ii) - rX_e(j,jj) - rrxij
                ryij = yi_(ii) - rY_e(j,jj) - rryij
                rzij = zi_(ii) - rZ_e(j,jj) - rrzij

                rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                rij   = rijsq**0.5_r8

                ai = myType(ii)
                aj = myType(jj)

                ms = sigij_e(ai,aj)
                me = epsij_e(ai,aj)
                mc = chaij_e(ai,aj)

                sr2 = ms*ms / rijsq
                sr6 = sr2*sr2*sr2

                de  = de  + me * sr6 * (sr6-1.0_r8)
                dec = dec + mc / rij

                end do; end do;

            end if

            end if ! outer if

        end do

    de  = 4.0_r8 * de
    dec = CHARGE_CONST * dec

    end subroutine Energy

    !********************









end module mod_energy_core



