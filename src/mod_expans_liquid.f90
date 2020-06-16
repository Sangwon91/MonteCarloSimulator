module mod_expans_liquid

    use mod_univ_const,           only : CHARGE_CONST, KtokJ
    use mod_character_tool,       only : IntToCha
    use mod_time_tool,            only : Delay   
    use mod_random_tool,          only : Ranf, RandomSettingAll   

    use mod_system_information,   only : molNum=>molNum, rCut=>rCut, beta=>beta, printStep=>printStep, temper=>temper
    use mod_molecule_information, only : nAtoms=>nAtoms, eps=>eps, sig=>sig, cha=>charge, &
                                         mcX=>molCrdX, mcY=>molCrdY, mcZ=>molCrdZ, mixingRule=>mixingRule, &
                                         color=>color
    use mod_config_information,   only : rX => rX, rY => rY, rZ => rZ, &
                                         hBox => hBox, ihBox => ihBox, vol => volume 

    use mod_movement_information, only : drMax=>drMax, dRotMax=>dRotMax, dBoxMax=>dBoxMax  

    implicit none
    
    private
    
    integer, parameter :: r8 = selected_real_kind(15,300)
    
    integer, parameter :: MXSTATE = 10, MXBIN = 151, MAXA = 30 ! Max Atom.. C70 never..
    
    logical :: adjustModeQ = .true., weightReadQ, weightOK

    integer :: totalPartition, state = 1, repeatStep, MaxXi, partition, eachStep = 1
    
    integer       :: ixi, ieestep, icycle, ixiold
    real(kind=r8) :: xi0, xi1, bin, mxbin_1, ebinmax, ene, eneold, relax, maxdw, delEs
    real(kind=r8), dimension(0:MAXA)     :: solX, solY, solZ, bl  ! solute (X,Y,Z), bond length
    real(kind=r8), dimension(0:MXSTATE)  :: xi, wxi, delw, delwOld, zeta, Pxi, nforxi, nbakxi, forxi, bakxi, oforbak
    real(kind=r8), dimension(0:MXSTATE,-MXBIN:MXBIN) :: nbeufor, nbeubak
    
    publiC :: icycle, solX, solY, solZ, ene, eneold, partition, wxi, MaxXi, delw, delwOld, weightOK, ixi, delEs
    public :: GetExpLiqInf, PrintExpLiqInf, SolEnergy, SolEnergyI, SolTransOldNew, SolRotateOldNew, ExpMCLiq, &
    &         GetMolecule, AutoOperation, MakeSolConfigView
              
    
    contains
    
    !********************
    
    subroutine GetExpLiqInf (filename)
    
        implicit none
        
        character(len=31), intent(in) :: filename
        !character(len=10)  :: dummyC
        integer, parameter :: UNIT_EEM = 10, UNIT_PART = 11
        integer            :: i, j, dummyPart, dmaxXi
        integer            :: openError
        real(kind=r8), dimension(0:MXSTATE) :: dweight

        real(kind=r8)      :: dummyUp, dummyDown, ls, le
 
        open(UNIT_EEM, file=filename)     

            read(UNIT_EEM,*) totalPartition           
            !read(UNIT_EEM,*) MaxXi
                  
            !if (MaxXi > MXSTATE) stop 'Max xi > Max State'   

            open(UNIT_PART, file='partition.part', status='old', iostat=openError)
            
            if (openError >  0) partition = 1
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
            
            if (MaxXi > MXSTATE) stop 'Max xi > Max State' 
            
            read(UNIT_EEM,*) weightReadQ

            if (weightReadQ) then
            
              do i = 1, totalPartition 
                
                read(UNIT_EEM,*) dummyPart, dweight(0:MaxXi)
                    
                if (dummyPart == partition) wxi = dweight

              end do  
            
              weightOK = .true.
            
            end if

        close(UNIT_EEM)
        
        ! initialization
        
        do i = 0, MaxXi
          xi(i) =  xi0 + ( xi1 - xi0 ) *dble(i)/dble( MaxXi )    
          wxi(i)= 0.0_r8
        end do
        
        ixi = 0

        zeta(0:MaxXi) = xi(0:MaxXi) - 1.0 
        where ( xi(0:MaxXi) > 1.0 ) 
          xi = 1.0               ! after LJ
        end where
        
        where ( zeta(0:MaxXi) < 0.0 )
          zeta = 0.0
        end where
        
        do i = 0, maxXi - 1
          delw(i) = wxi(i+1) - wxi(i)
	    end do   

        ieestep = 0
        relax   = 1.0
        
        adjustModeQ = .true.
        
        do i = 0, maxXi
          Pxi(i)    = 0    ! probabilities
          nforxi(i) = 0 
          nbakxi(i) = 0 
          forxi(i)  = 0 
          bakxi(i)  = 0
        end do     
        
        do i = 0, maxXi
	      do j = -MXBIN, MXBIN
	        nbeufor(i,j) = 0
	        nbeubak(i,j) = 0
	      end do
          oforbak(i) = 0     ! for what?
	    end do     

        bin = 0.1
        mxbin_1 = MXBIN - 1
        ebinmax = bin * dble(mxbin_1)
        
        weightOK = .false.
    
    !write(*,"('out partition = ', I2)") partition
        
    end subroutine GetExpLiqInf

    !********************
    
    subroutine PrintExpLiqInf (unit_)
    
        implicit none
        
        integer, intent(in), optional :: unit_
        
        if (present(unit_)) then
        
            write(unit_,"(A/)") 'expanded ensemble information ->'
            
            write(unit_,"('total partition number = ', I5)")  totalPartition
            write(unit_,"('current partition      = ', I5)")  partition
            write(unit_,"('max Xi                 = ', I5)")  maxXi
            write(unit_,"('lambda = ',11F10.6)") xi(0:MaxXi)
            write(unit_,"('zeta   = ',11F10.6)") zeta(0:MaxXi)
            write(unit_,"('weight = ',11F10.6)") wxi(0:MaxXi)
            write(unit_,"('Is exist initial weight? (T/F) = ', L2 )") weightReadQ
            write(unit_,"()")
        
        else
        
            write(*,"(A/)") 'expanded ensemble information ->'
            
            write(*,"('total partition number = ', I5)")  totalPartition
            write(*,"('current partition      = ', I5)")  partition
            write(*,"('max Xi                 = ', I5)")  maxXi
            write(*,"('lambda = ',11F10.6)") xi(0:MaxXi)
            write(*,"('zeta   = ',11F10.6)") zeta(0:MaxXi)
            write(*,"('weight = ',11F10.6)") wxi(0:MaxXi)
            write(*,"('Is exist initial weight? (T/F) = ', L2 )") weightReadQ
            write(*,"()")
        
        end if
    
    end subroutine PrintExpLiqInf
    
    !********************
    
    subroutine AutoOperation

        implicit none
        
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
            call system(command)
           
            call Delay(1)
        
        end do
        
        !call Delay(1)
        
        open(1, file='eemout.txt', status='replace', position='append')
        
        write(1, "(2A30)") 'partition', 'delta free energy'
        
        partition = 1
        
    end subroutine AutoOperation
    
    !********************
    
    subroutine GetMolecule ! only pure system

        implicit none
        
        integer :: i
        
        if ( nAtoms > MAXA ) stop 'error atom num > Max atom'
        
        solX(0:nAtoms) = rX(1,0:nAtoms)
        solY(0:nAtoms) = rY(1,0:nAtoms)
        solZ(0:nAtoms) = rZ(1,0:nAtoms)
        
        ! go molecule center to origin
        
        solX(0:nAtoms) = solX(0:nAtoms) - solX(0)
        solY(0:nAtoms) = solY(0:nAtoms) - solY(0)
        solZ(0:nAtoms) = solZ(0:nAtoms) - solZ(0)
        
        do i = 1, nAtoms ! save relative coordinate
        
          bl(i) = solX(i)**2 + solY(i)**2 + solZ(i)**2
          bl(i) = bl(i)**0.5
          
          if ( bl(i) < 10.e-7 ) bl(i) = -1.0  ! origin posotion atom
        
        end do
        
        call RandomSettingAll
        
        i = int(Ranf())
        
        solX(0:nAtoms) = solX(0:nAtoms) + (Ranf()-0.5_r8) * hBox(1,1)
        solY(0:nAtoms) = solY(0:nAtoms) + (Ranf()-0.5_r8) * hBox(2,2)
        solZ(0:nAtoms) = solZ(0:nAtoms) + (Ranf()-0.5_r8) * hBox(3,3)
        
        !write(*,*) hBox(1,1),hBox(2,2),hBox(3,3)

    end subroutine GetMolecule
    
    !********************
    
    subroutine Shrink (ixi_)
    
        implicit none
        
        integer, intent(in)            :: ixi_ ! current index of xi 
        
        integer :: i
        real(kind=r8), dimension(MAXA) :: cbl  ! current bond length
        real(kind=r8)                  :: scl, ox, oy, oz
        solX(1:nAtoms) = solX(1:nAtoms) - solX(0)
        solY(1:nAtoms) = solY(1:nAtoms) - solY(0)
        solZ(1:nAtoms) = solZ(1:nAtoms) - solZ(0)
        
        do i = 1, nAtoms
        
          cbl(i) = solX(i)**2 + solY(i)**2 + solZ(i)**2
          cbl(i) = cbl(i)**0.5
          
          if ( cbl(i) < 10.e-7 ) cbl(i) = -1.0  ! origin posotion atom
        
        end do    

        if ( xi(ixi_) > 0.001 ) then
		
          scl = xi(ixi_)
		  
        else
		
          scl = 0.001
		  
        end if

          do i = 1, nAtoms
        
            solX(i) = solX(i) / cbl(i) * bl(i) * scl
            solY(i) = solY(i) / cbl(i) * bl(i) * scl
            solZ(i) = solZ(i) / cbl(i) * bl(i) * scl
        
          end do

        solX(1:nAtoms) = solX(1:nAtoms) + solX(0)
        solY(1:nAtoms) = solY(1:nAtoms) + solY(0)
        solZ(1:nAtoms) = solZ(1:nAtoms) + solZ(0)
    
    end subroutine Shrink
    
    !********************
    
    subroutine SolEnergy (hx_, hY_, hZ_, de, dec, dBox, diBox) ! energy of solute, use current ixi
    
        implicit none
        
        integer :: j, ii, jj
        
        real(kind=r8), dimension(0:MAXA), intent(in) :: hX_, hY_, hZ_    ! solute coordinate
        real(kind=r8), intent(out) :: de, dec
        real(kind=r8), dimension(3,3), intent(in) :: dBox, diBox
        
        real(kind=r8) :: rcutsq
        real(kind=r8) :: b11,b22,b33, ib11,ib22,ib33
        real(kind=r8) :: sr2, sr6, rxij, ryij, rzij, rij, rijsq, dxij, dyij, dzij, srxij, sryij, srzij, &
                         rrxij, rryij, rrzij  
        real(kind=r8) :: ms, me, mc, xi1, xi3, zeta2 ! mean sigma, epsilon, charge    
        
        real(kind=r8), dimension(MAXA, MAXA) :: sigS, epsS, chaS ! sigma, epsilon, charge with solute

        de  = 0.0_r8
        dec = 0.0_r8
        rcutsq = rCut**2
        
        b11 = dbox(1,1)
        b22 = dbox(2,2)
        b33 = dbox(3,3)
        
        ib11 = dibox(1,1)
        ib22 = dibox(2,2)
        ib33 = dibox(3,3)     
        
        xi1   = xi(ixi)
        xi3   = xi1 ** 0.333333333333333333333333
        
        zeta2 = zeta(ixi) ** 0.5

        select case(mixingRule)
        
            case('a','A')
            
              do ii = 1, nAtoms; do jj = 1, nAtoms
                sigS(ii,jj) = xi3*(sig(ii) + sig(jj))*0.5_r8
              end do; end do    
                    
            case('g','G')
            
              do ii = 1, nAtoms; do jj = 1, nAtoms
                sigS(ii,jj) = xi3*(sig(ii) * sig(jj))**0.5_r8
              end do; end do   
            
            case default
                stop 'unclassified mixing rule'
        
        end select

        do ii = 1, nAtoms; do jj = 1, nAtoms
            epsS(ii,jj) = xi1*(eps(ii) * eps(jj))**0.5_r8
        end do; end do 

        do ii = 1, nAtoms; do jj = 1, nAtoms
            chaS(ii,jj) = zeta2*(cha(ii) * cha(jj))
        end do; end do   

        do j = 1, molNum

            rrxij = hX_(0) - rX(j,0)
            rryij = hY_(0) - rY(j,0)
            rrzij = hZ_(0) - rZ(j,0)
        
            dxij = ib11*rrxij
            dyij = ib22*rryij 
            dzij = ib33*rrzij
            
            if (dxij >  0.5_r8) dxij = dxij - 1.0_r8
            if (dxij < -0.5_r8) dxij = dxij + 1.0_r8
            
            if (dyij >  0.5_r8) dyij = dyij - 1.0_r8
            if (dyij < -0.5_r8) dyij = dyij + 1.0_r8
            
            if (dzij >  0.5_r8) dzij = dzij - 1.0_r8
            if (dzij < -0.5_r8) dzij = dzij + 1.0_r8
            
            srxij = b11*dxij 
            sryij = b22*dyij
            srzij = b33*dzij
        
            rijsq = srxij*srxij + sryij*sryij + srzij*srzij 

            if (rijsq < rcutsq) then

                rrxij = rrxij - srxij
                rryij = rryij - sryij
                rrzij = rrzij - srzij

                do ii = 1, nAtoms; do jj = 1, nAtoms
                
                rxij = hX_(ii) - rX(j,jj) - rrxij
                ryij = hY_(ii) - rY(j,jj) - rryij
                rzij = hZ_(ii) - rZ(j,jj) - rrzij           
           
                rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                
                ms = sigS(ii,jj)
                me = epsS(ii,jj)

                sr2 = ms*ms / rijsq
                sr6 = sr2*sr2*sr2   

                de  = de + me * sr6 * (sr6-1.0_r8)

                if ( zeta2 > 10.e-7 ) then ! if charge stage
                    
                   rij   = rijsq**0.5_r8
                    
                   mc  = chaS(ii,jj)
                   dec = dec + mc / rij  

                end if
                  
                end do; end do;
        
            end if
        
        end do
        
    de  = 4.0_r8 * de
    dec = CHARGE_CONST * dec          
   
    end subroutine SolEnergy
    
    !********************
    
    subroutine SolEnergyI (xOld_, yOld_, zOld_, de, dec, dBox, diBox)
    
        implicit none
        
        integer :: ii, jj
        
        real(kind=r8), dimension(0:nAtoms), intent(in) :: xOld_, yOld_, zOld_    
        real(kind=r8), intent(out) :: de, dec
        real(kind=r8), dimension(3,3), intent(in) :: dBox, diBox

        real(kind=r8) :: rcutsq
        real(kind=r8) :: b11,b22,b33, ib11,ib22,ib33
        real(kind=r8) :: sr2, sr6, rxij, ryij, rzij, rij, rijsq, dxij, dyij, dzij, srxij, sryij, srzij, &
                         rrxij, rryij, rrzij  
        real(kind=r8) :: ms, me, mc, xi1, xi3, zeta2 ! mean sigma, epsilon, charge    
        
        real(kind=r8), dimension(MAXA, MAXA) :: sigS, epsS, chaS ! sigma, epsilon, charge with solute

        de  = 0.0_r8
        dec = 0.0_r8
        rcutsq = rCut**2
        
        b11 = dbox(1,1)
        b22 = dbox(2,2)
        b33 = dbox(3,3)
        
        ib11 = dibox(1,1)
        ib22 = dibox(2,2)
        ib33 = dibox(3,3)     
        
        xi1   = xi(ixi)
        xi3   = xi1 ** 0.333333333333333333333333
        
        zeta2 = zeta(ixi) ** 0.5

        select case(mixingRule)
        
            case('a','A')
            
              do ii = 1, nAtoms; do jj = 1, nAtoms
                sigS(ii,jj) = xi3*(sig(ii) + sig(jj))*0.5_r8
              end do; end do    
                    
            case('g','G')
            
              do ii = 1, nAtoms; do jj = 1, nAtoms
                sigS(ii,jj) = xi3*(sig(ii) * sig(jj))**0.5_r8
              end do; end do   
            
            case default
                stop 'unclassified mixing rule'
        
        end select

        do ii = 1, nAtoms; do jj = 1, nAtoms
            epsS(ii,jj) = xi1*(eps(ii) * eps(jj))**0.5_r8
        end do; end do 

        do ii = 1, nAtoms; do jj = 1, nAtoms
            chaS(ii,jj) = zeta2*(cha(ii) * cha(jj))
        end do; end do     
 
        rrxij = SolX(0) - xOld_(0)
        rryij = SolY(0) - yOld_(0)
        rrzij = SolZ(0) - zOld_(0)
    
        dxij = ib11*rrxij
        dyij = ib22*rryij 
        dzij = ib33*rrzij
        
        if (dxij >  0.5_r8) dxij = dxij - 1.0_r8
        if (dxij < -0.5_r8) dxij = dxij + 1.0_r8
        
        if (dyij >  0.5_r8) dyij = dyij - 1.0_r8
        if (dyij < -0.5_r8) dyij = dyij + 1.0_r8
        
        if (dzij >  0.5_r8) dzij = dzij - 1.0_r8
        if (dzij < -0.5_r8) dzij = dzij + 1.0_r8
        
        srxij = b11*dxij 
        sryij = b22*dyij
        srzij = b33*dzij
    
        rijsq = srxij*srxij + sryij*sryij + srzij*srzij 

        if (rijsq < rcutsq) then

            rrxij = rrxij - srxij
            rryij = rryij - sryij
            rrzij = rrzij - srzij

            do ii = 1, nAtoms; do jj = 1, nAtoms
            
                rxij = SolX(ii) - xOld_(jj) - rrxij
                ryij = SolY(ii) - YOld_(jj) - rryij
                rzij = SolZ(ii) - ZOld_(jj) - rrzij           
           
                rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                
                ms = sigS(ii,jj)
                me = epsS(ii,jj)

                sr2 = ms*ms / rijsq
                sr6 = sr2*sr2*sr2   

                de  = de + me * sr6 * (sr6-1.0_r8)

                if ( zeta2 > 10.e-7 ) then ! if charge stage
                    
                   rij   = rijsq**0.5_r8
                    
                   mc  = chaS(ii,jj)
                   dec = dec + mc / rij  

                end if
              
            end do; end do;
    
        end if
        
    de  = 4.0_r8 * de
    dec = CHARGE_CONST * dec       

    end subroutine SolEnergyI
    
    !********************
    
    subroutine ExpMCliq

        implicit none
        
        real(kind=r8), parameter  :: Infinity = 10.e100, Zero = 0.0_r8
        
        integer                   :: eestep, i, j, iter
        real(kind=r8)             :: e, ec, arg, prob, paccfor, paccbak, probeu, beu, dw, dwn ! energy of old configuration and state
        
        real(kind=r8)             :: muTot
        !integer, parameter :: BYPASS = 2200
        
        call SolEnergy (solX, SolY, SolZ, e, ec, hBox, ihBox) ! get old energy

        ixiold = ixi
        eneold = e + ec   ! energy of state ixi
        delEs  = 0.0

        if( ixiold == 0 ) then        
          if( ranf() > 0.5 ) then                 ! stay
            ixi = 0      
          else 
            ixi = 1
            nforxi(ixiold) = nforxi(ixiold) + 1   ! number of forward trials
          endif
        elseif( ixiold == MaxXi ) then
          if( ranf() > 0.5 ) then                 !  stay
            ixi = MaxXi    
          else
            ixi = MaxXi - 1
            nbakxi(ixiold) = nbakxi(ixiold) + 1
          endif
        else
          if( ranf() > 0.5 ) then
            ixi = ixiold + 1
            nforxi(ixiold) = nforxi(ixiold) + 1
          else
            ixi = ixiold - 1
            nbakxi(ixiold) = nbakxi(ixiold) + 1
          endif
        endif     
        !write(*,"(I2, ' -> ', I2)") ixiold, ixi
        if ( ixi == ixiold ) go to 1             ! do nothing and bypass
        !write(*,*) 'shrink TEST'      
        call Shrink(ixi)                         ! shrink to new xi
        !write(*,*) 'shrink OK'
        !call system('pause')
        
        call SolEnergy (solX, solY, solZ, e, ec, hBox, ihBox) ! get new energy after shrinked
        
        ene = e + ec
        
        arg = - beta * ( ene - eneold )
        
        if( adjustModeQ )then

           if( -arg >  ebinmax ) then
             j = MXBIN                                     ! very large positive energy difference
           elseif( -arg < -ebinmax ) then
             j = -MXBIN                                    ! very large negative energy difference
           else  ! usual
             j = nint( (-arg) / bin )                      ! beta* del U = -arg, bin size = 0.1, range -10 ~ 10 in unit of kT
           endif
             if( iabs(j) > MXBIN ) stop 'bin index error'    ! was error spot 

           if( ixiold < ixi ) then
             nbeufor(ixiold,j) = nbeufor(ixiold,j) + 1     ! count transition energy
             arg = arg + delw(ixiold) 
           else ! ixiold > ixi
             nbeubak(ixiold,j) = nbeubak(ixiold,j) + 1     ! count transition energy
             arg = arg - delw(ixi)                                   ! not delw(ixiold) 
           endif

        else  ! no adjustment ! original method

           arg = arg + wxi(ixi) - wxi(ixiold) 

        endif ! chk_adaptive
        
        if( arg > 0.d0 ) then
          prob = 1.d0
        else  
          prob = dexp( arg ) 
        endif
        
        if ( Ranf() < prob ) then ! accept
        
          if( ixiold < ixi ) then                 ! successful forward trial 
            forxi(ixiold) = forxi(ixiold) + 1
          elseif( ixiold > ixi ) then             ! successful backward trial 
            bakxi(ixiold) = bakxi(ixiold) + 1
!         else ! ixiold = ixi                     ! do noting
          endif
          
          delEs = ene - eneold
        
        else                      ! reject
        
          ixi = ixiold   ! back to origin state
          ene = eneold
          
          call Shrink (ixi)
          
          delEs = 0.0 
        
        end if
        
        if( adjustModeQ .and. mod(icycle,100) == 0 ) then  !new  !every 100 cycle

        do i = 0, MaxXi-1
       
            if( nforxi(i) < 1 .or. nbakxi(i+1) < 1 )      cycle  ! after at least 1 samplings   
            if( nforxi(i) + nbakxi(i+1) < oforbak(i)+10 ) cycle  ! do not renew weight before 10 times of transition

            dw = delw(i)
            relax = 1.

10000   do iter = 1, 50

             paccfor = 0.
             do j = -MXBIN, MXBIN
             
                if( nbeufor(i,j) == 0 ) cycle                      ! unvisited bin
                probeu = dble( nbeufor(i,j) )/dble( nforxi(i) )
                if( probeu < 0.001 ) cycle                         ! too small prob
                  beu = dble(j) * bin
                if( j == -MXBIN ) beu = -Infinity                  ! not happen in most cases
                if( j ==  MXBIN ) beu =  Infinity                  ! hard-core repulsion
                  arg = - beu + dw
                if( arg > Zero ) then
                 paccfor = paccfor + probeu
                else
                 paccfor = paccfor + exp(arg) * probeu
                 
                  endif
             end do ! j

             paccbak = 0.
             do j = -MXBIN, MXBIN
             
                if( nbeubak(i+1,j) == 0 ) cycle  ! unvisited bin
                    probeu = dble( nbeubak(i+1,j) )/dble( nbakxi(i+1) )
                if( probeu < 0.001 ) cycle
                    beu = dble(j) * bin
                if( j == -MXBIN ) beu = -Infinity      
                if( j ==  MXBIN ) beu =  Infinity                ! hard-core repulsion
                    arg = - beu - dw
                if( arg > Zero ) then
                    paccbak = paccbak + probeu
                else
                    paccbak = paccbak + exp(arg) * probeu
                endif
                  
             end do ! j

             if( paccfor < Zero .or. paccbak < Zero ) then  ! very rare case if it occurs
                  write(* ,*)'warning: adaptive eemc, acc prob zero'          
                  write(13,*)'warning: adaptive eemc, acc prob zero'          
                  write(* ,"(4i8,2f7.3)")i,nforxi(i),nbakxi(i+1),paccfor,paccbak !check
                  write(13,"(4i8,2f7.3)")i,nforxi(i),nbakxi(i+1),paccfor,paccbak !check
                  dwn = dw  ! escaping the loop
             else      
               
        !	   // Update weight difference, usual case

                  dwn = dw - log( paccfor / paccbak ) * relax
             endif

!             if( chk_verbose == 'YES' ) 
!         :      write(*,"(4i8,f6.2,2f7.3)") i, nforxi(i), nbakxi(i+1),
!         :                                  iter, dwn, paccfor, paccbak     !check

             if( abs(dwn-dw) < 0.001 ) then
                delw(i) = dwn
                oforbak(i) = nforxi(i) + nbakxi(i+1)
                exit  ! converged
             endif

             dw = dwn
             
        end do ! iter
 
!       // error handling  --not happen usually
        if( iter > 50 ) then
          write(13,"('--warning: weight not converged! ')")
          write(* ,"('--warning: weight not converged! ')")
	      if( relax > 0.2 ) then
	        relax = relax * 0.8
            write(* ,"('--relax factor = ',f6.4)") relax
	      go to 10000
	      else
            write(* ,"('--not converged! but will use current weight')")
	        delw(i) = dwn
            oforbak(i) = nforxi(i) + nbakxi(i+1)
	      endif
	    endif
!        // error handling

         end do ! i

        !if( mod(icycle,printStep)==0 ) then
          write(* ,"(/i9,9x,21(f10.4))") icycle,( delw(i), i= 0, MaxXi-1)
          !write(34,"(i9,9x,21(f6.2)) ") icycle,( delw(i), i= 0, MaxXi-1)
        !endif
        if ( minval(dabs(delw(0:maxXi-1))) > 0.0001 ) then ! non zero condition
          maxdw = maxval(dabs((delwOld(0:maxXi-1) - delw(0:maxXi-1))/delw(0:maxXi-1)))
          
          write(*,"('max dw = ', F15.6, ' %')") maxdw*100.0
            
          if ( maxdw < 0.001 ) then ! 0.1% error
            adjustModeQ = .false.
            weightOK    = .true.
            write(*,"('start hopping stage')")
            wxi(0) = 0.0
            do i = 1, maxXi
              wxi(i) = wxi(i-1) + delw(i-1)
            end do
            write(* ,"(/i9,9x,21(f10.4))") icycle,( wxi(i), i= 0, MaxXi)
            write(1 ,"(/i9,9x,21(f10.4))") icycle,( wxi(i), i= 0, MaxXi) ! unit 1 : monitoring file            
          end if
        
          delwOld(0:maxXi-1) = delw(0:maxXi-1)
        
        end if

	    endif ! chk_adaptive
        
1   continue

    if ( xi(ixi) < 0.0001 ) call RandomTrans

    if ( .not.adjustModeQ ) then
      Pxi(ixi) = Pxi(ixi) + 1 ! probabilities count
      if( mod(icycle,printStep)==0 ) then
        write(* ,"(/i9,9x,21(I10))") icycle,( nint(Pxi(i)), i= 0, MaxXi)
        write(1 ,"(/i9,9x,21(I10))") icycle,( nint(Pxi(i)), i= 0, MaxXi) ! unit 1 : monitoring file
        
        if ( minval(Pxi(0:maxXi)) > 0.001 ) then
          open(2, file='monitorMu'//IntToCha(partition)//'.txt', position='append')
            muTot = wxi(MaxXi)-wxi(0)-log(Pxi(MaxXi)/Pxi(0))
            write(* ,"(/'del mu(-) = ',9x,21(F10.6))") ( wxi(i)-wxi(i-1)-log(Pxi(i)/Pxi(i-1)), i= 1, MaxXi), temper*KtokJ*muTot
            write(1 ,"(/'del mu(-) = ',9x,21(F10.6))") ( wxi(i)-wxi(i-1)-log(Pxi(i)/Pxi(i-1)), i= 1, MaxXi), temper*KtokJ*muTot ! unit 1 : monitoring file
            write(2 ,"(I15, F15.6)") icycle, temper*KtokJ*muTot ! unit 1 : monitoring file
          close(2)
        end if
      end if
    end if
 
    end subroutine ExpMCliq
    
    subroutine RandomTrans
    
        use mod_random_tool, only : Ranf
    
        implicit none
        
        solX(0:nAtoms) = solX(0:nAtoms) - solX(0)
        solY(0:nAtoms) = solY(0:nAtoms) - solY(0)
        solZ(0:nAtoms) = solZ(0:nAtoms) - solZ(0)

        solX(0:nAtoms) = solX(0:nAtoms) + (Ranf()-0.5_r8) * hBox(1,1)
        solY(0:nAtoms) = solY(0:nAtoms) + (Ranf()-0.5_r8) * hBox(2,2)
        solZ(0:nAtoms) = solZ(0:nAtoms) + (Ranf()-0.5_r8) * hBox(3,3)        
    
    end subroutine RandomTrans
    
    !********************
    
    subroutine SolTransOldNew (xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
    
        use mod_random_tool, only : Ranf
        
        implicit none
        
        real(kind=r8), dimension(0:nAtoms), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
        real(kind=r8) :: noX, noY, noZ
        real(kind=r8) :: dX, dY, dZ
        real(kind=r8) :: halfBox, Box
        
        halfBox = 0.5_r8 * hBox(1,1)
        Box     = hBox(1,1)

        xOld_ = SolX(0:nAtoms)
        yOld_ = SolY(0:nAtoms)
        zOld_ = SolZ(0:nAtoms)
        
        dX = drMax*(2.0_r8*Ranf()-1.0_r8)
        dY = drMax*(2.0_r8*Ranf()-1.0_r8)
        dZ = drMax*(2.0_r8*Ranf()-1.0_r8)
        
        noX = hBox(1,1)*dX 
        noY = hBox(2,2)*dY
        noZ = hBox(3,3)*dZ
        
        xNew_ = xOld_ + noX
        yNew_ = yOld_ + noY
        zNew_ = zOld_ + noZ
        
        if      (xNew_(0) >  halfBox) then; xNew_ = xNew_ - Box
        else if (xNew_(0) < -halfBox) then; xNew_ = xNew_ + Box; end if
        
        if      (yNew_(0) >  halfBox) then; yNew_ = yNew_ - Box
        else if (yNew_(0) < -halfBox) then; yNew_ = yNew_ + Box; end if
        
        if      (zNew_(0) >  halfBox) then; zNew_ = zNew_ - Box
        else if (zNew_(0) < -halfBox) then; zNew_ = zNew_ + Box; end if

    end subroutine SolTransOldNew 
   
    !********************
    
    subroutine SolRotateOldNew (xOld_, yOld_, zOld_, xNew_, yNew_, zNew_)
    
        use mod_random_tool, only : Ranf
    
        implicit none
        
        real(kind=r8), dimension(0:nAtoms), intent(out) :: xOld_, yOld_, zOld_, xNew_, yNew_, zNew_
        integer :: xyz ! xyz : xÃà : x, yÃà : y, zÃà : z
        real(kind=r8), dimension(0:nAtoms) :: dX, dY, dZ
        real(kind=r8) :: cos_, sin_, ang, oX, oY, oZ ! dRotMax
        
        xOld_ = solX(0:nAtoms)
        yOld_ = solY(0:nAtoms)
        zOld_ = solZ(0:nAtoms)
        
        oX = xOld_(0)
        oY = yOld_(0)
        oZ = zOld_(0)
        
        dX = xOld_ - oX
        dY = yOld_ - oY
        dZ = zOld_ - oZ
        
        xyz  = int(3.0_r8*Ranf() + 1.0_r8)
        ang  = dRotMax * (2.0_r8*Ranf()-1.0_r8)
        cos_ = cos(ang)
        sin_ = sin(ang)
        
        if (xyz == 1) then
        
            xNew_ = dX                     + oX
            yNew_ =     cos_*dY - sin_*dZ  + oY
            zNew_ =     sin_*dY + cos_*dZ  + oZ
         
        else if (xyz == 2) then
        
            xNew_ = cos_*dX      -sin_*dZ  + oX
            yNew_ =          dY            + oY
            zNew_ = sin_*dX      +cos_*dZ  + oZ
        
        else if (xyz == 3) then
        
            xNew_ = cos_*dX - sin_*dY      + oX
            yNew_ = sin_*dX + cos_*dY      + oY
            zNew_ =                    dZ  + oZ        
        
        else; stop 'critical error occur in rotation'; end if       
    
    end subroutine SolRotateOldNew
    
    !********************
    
    subroutine MakeSolConfigView (filename)
    
        use mod_math_tool, only : Det

        implicit none

        integer, parameter :: UNIT_VIEW = 10
        integer            :: i, j;
        character(len=31), intent(in) :: filename

        open(UNIT_VIEW,file=filename)
    
            write(UNIT_VIEW,*) (1+molNum)*nAtoms

            do i=1, molNum; do j=1, nAtoms
                write(UNIT_VIEW,"(4F20.5)") color(j), rx(i,j), ry(i,j), rz(i,j)
            end do; end do
            
            do i = 1, nAtoms
                write(UNIT_VIEW,"(4F20.5)") 8.0_r8, solX(i), solY(i), solZ(i)
            end do;
            
            write(UNIT_VIEW,*) 0
            write(UNIT_VIEW,*) Det (hBox,3)**(1.0_r8/3.0_r8)
            
        close(UNIT_VIEW)
    
    end subroutine MakeSolConfigView
    
    !********************

end module mod_expans_liquid