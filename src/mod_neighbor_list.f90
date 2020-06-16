module mod_neighbor_list

    use mod_system_information, only : rCut, rvCut, molNum, neighborListQ
    use mod_molecule_information
    use mod_config_information

    implicit none

    integer, parameter, private :: r8 = selected_real_kind(15,300)

    private ! default: private

!   -------------------------------------------------------
    type NeighborList  ! derived type of nlist
!   -------------------------------------------------------
        integer, dimension(:), pointer :: n ! neighbor
!   -------------------------------------------------------
    end type NeighborList

    type(NeighborList), dimension(:), allocatable :: nlist

    logical :: setupQ

    integer :: nUpdate
    integer, dimension(:), allocatable :: savedGrpX, savedGrpY, savedGrpZ

    public :: NeighborList, nlist

    public :: InitializeNeighborList, CheckList, PrintNeighborList, DeallocateNeighborList, SetupList



    contains

!   -------------------------------------------------------
    subroutine InitializeNeighborList()
!   -------------------------------------------------------
        implicit none

        integer :: ig, nl, m       

!   -------------------------------------------------------
 
        nUpdate = 0

        if (.not.allocated(nlist)) then 

            allocate(nlist(nGlobalGroup))
            do ig = 1, nGlobalGroup

                allocate(nlist(ig)%n(0:nGlobalGroup))

            end do ! ig-loop 

        end if

        do ig = 1, nGlobalGroup
  
            do nl = 1, nGlobalGroup

                nlist(ig)%n(nl) = 0 ! initialize to zero

            end do ! jg-loop

        end do ! ig-loop             

        allocate(savedGrpX(nGlobalGroup), savedGrpY(nGlobalGroup), savedGrpZ(nGlobalGroup))  

        savedGrpX = grpX ! save old configuration to check list
        savedGrpY = grpY
        savedGrpZ = grpZ
        
        
        setupQ = .true.
        call SetupList() ! setup linitial list

        if (.not.neighborListQ) then

            setupQ = .false.

            do ig = 1, nGlobalGroup
                
                m = myMolIndex(ig)
                nlist(ig)%n(0) = 0

                do nl = 1, nGlobalGroup

                    if (nl == ig) cycle
                    if (m == myMolIndex(nl)) cycle ! neglect inner group
                    nlist(ig)%n(0)  = nlist(ig)%n(0) + 1
                    nlist(ig)%n(nlist(ig)%n(0)) = nl ! save all pair and never update
                    
                end do ! jg-loop

            end do ! ig-loop      

        end if

        write(1,"('neighbor list is initialized')")
        !call PrintNeighborList()

    end subroutine InitializeNeighborList




!   -------------------------------------------------------
    subroutine SetupList()
!   -------------------------------------------------------
        implicit none

        integer :: im, in, ng, mg, ig, jg
        integer, dimension(:), pointer :: ilist, jlist

        real(kind=r8), parameter :: ONE = 1.0_r8, HALF = 0.5_r8
        real(kind=r8) :: rvcutsq
        ! g: group, d: reduced, p: peridic
        real(kind=r8) :: gxij, gyij, gzij, gijsq, dxij, dyij, dzij, pxij, pyij, pzij
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33
  
!   -------------------------------------------------------

        if (.not.setupQ) return

        rvCutsq = rvCut*rvCut

        b11  =  hbox(1,1)
        b12  =  hbox(1,2)
        b13  =  hbox(1,3)
        b22  =  hbox(2,2)
        b23  =  hbox(2,3)
        b33  =  hbox(3,3)
        
        ib11 = ihbox(1,1)
        ib12 = ihbox(1,2)
        ib13 = ihbox(1,3)
        ib22 = ihbox(2,2)
        ib23 = ihbox(2,3)
        ib33 = ihbox(3,3)

        do ig = 1, nGlobalGroup
            nlist(ig)%n(0) = 0
        end do ! ig-loop

        do im = 1, molNum - 1

            do in = im + 1, molNum
      
                do mg = 1, nLocalGroup
                    
                    ig = globalGroupIndex(im,mg)
                    ilist => nlist(ig)%n

                    do ng = 1, nLocalGroup

                        jg = globalGroupIndex(in,ng)
                        jlist => nlist(jg)%n

                        gxij = grpX(ig) - grpX(jg)                        
                        gyij = grpY(ig) - grpY(jg)
                        gzij = grpZ(ig) - grpZ(jg)  

                        dxij = ib11*gxij + ib12*gyij + ib13*gzij
                        dyij =             ib22*gyij + ib23*gzij
                        dzij =                         ib33*gzij

                        if (dxij >  HALF) dxij = dxij - ONE
                        if (dxij < -HALF) dxij = dxij + ONE  

                        if (dyij >  HALF) dyij = dyij - ONE
                        if (dyij < -HALF) dyij = dyij + ONE

                        if (dzij >  HALF) dzij = dzij - ONE
                        if (dzij < -HALF) dzij = dzij + ONE 
                          
                        pxij = b11*dxij + b12*dyij + b13*dzij
                        pyij =            b22*dyij + b23*dzij
                        pzij =                       b33*dzij  

                        gijsq = pxij*pxij + pyij*pyij + pzij*pzij    

                        if (gijsq > rvcutsq ) cycle ! neglect the region of r > rvcut               

                        ilist(0) = ilist(0) + 1 ! number of neighbors
                        ilist(ilist(0)) = jg

                        jlist(0) = jlist(0) + 1
                        jlist(jlist(0)) = ig   
                        
                    end do ! ng-loop 

                end do ! mg-loop

            end do ! in-loop

        end do ! im-loop

    end subroutine SetupList








!   -------------------------------------------------------
    subroutine CheckList()
!   -------------------------------------------------------
        implicit none

        logical updateQ
        
        integer :: m, mg, ig

        real(kind=r8), parameter :: saveFactor = 0.4_r8, ZERO = 0.0_r8, ONE = 1.0_r8, HALF = 0.5_r8
        real(kind=r8) :: dispLimit
        real(kind=r8) :: gxij, gyij, gzij, gijsq, dxij, dyij, dzij, pxij, pyij, pzij
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33

!   -------------------------------------------------------

        if (.not.neighborListQ) return

        updateQ = .false.

        dispLimit = saveFactor * (rvCut - rCut)

        b11  =  hbox(1,1)
        b12  =  hbox(1,2)
        b13  =  hbox(1,3)
        b22  =  hbox(2,2)
        b23  =  hbox(2,3)
        b33  =  hbox(3,3)
        
        ib11 = ihbox(1,1)
        ib12 = ihbox(1,2)
        ib13 = ihbox(1,3)
        ib22 = ihbox(2,2)
        ib23 = ihbox(2,3)
        ib33 = ihbox(3,3)

        do m = 1, molNum
 
            do mg = 1, nLocalGroup

                ig = globalGroupIndex(m,mg)
                
                ! check z direction displacement
                gzij = grpZ(ig) - savedGrpZ(ig) ! save to subindex 0
                if (gzij < ZERO) gzij = -gzij ! abs value
                dzij = ib33 * gzij            ! scaled coordinate
                if (dzij > HALF) dzij =  dzij - ONE
                if (dzij < ZERO) dzij = -dzij ! abs value
                gzij = b33 * dzij  
                if (gzij > dispLimit) updateQ = .true.
                if (updateQ) exit

                
                ! check y direction displacement
                gyij = grpY(ig) - savedGrpY(ig)
                if (gyij < ZERO) gyij = -gyij ! abs value
                dyij = ib22 * gyij + ib23 * gzij
                if (dyij > HALF) dyij =  dyij - ONE
                if (dyij < ZERO) dyij = -dyij ! abs value
                gyij = b22 * dyij + b23 * dzij
                if (gyij > dispLimit) updateQ = .true.
                if (updateQ) exit



                ! check x direction displacement
                gxij = grpX(ig) - savedGrpX(ig)
                if (gxij < ZERO) gxij = -gxij ! abs value
                dxij = ib11 * gxij + ib12 * gyij + ib13 * gzij
                if (dxij > HALF) dxij =  dxij - ONE
                if (dxij < ZERO) dxij = -dxij ! abs value
                gxij = b11 * dxij + b12 * dyij + b13 * dzij
                if (gxij > dispLimit) updateQ = .true.
                if (updateQ) exit


            end do ! mg-loop
                
            if (updateQ) exit

        end do ! m-loop (molecule)

        if (updateQ) then
            
            nUpdate = nUpdate + 1
            savedGrpX = grpX
            savedGrpY = grpY
            savedGrpZ = grpZ

            call SetupList()

        end if

    end subroutine CheckList





!   -------------------------------------------------------
    subroutine PrintNeighborList()
!   -------------------------------------------------------
        implicit none

        integer, parameter :: UNIT_ = 10, UNIT__ = 11

        integer, dimension(:), pointer :: list

        real(kind=r8), parameter :: ONE = 1.0_r8, HALF = 0.5_r8
        integer :: ig, maxList, nl, jg, mol
        real(kind=r8) :: gx, gy, gz
        real(kind=r8) :: gxij, gyij, gzij, gijsq, dxij, dyij, dzij, pxij, pyij, pzij
        real(kind=r8) :: b11,b12,b13,b22,b23,b33, ib11,ib12,ib13,ib22,ib23,ib33
!   -------------------------------------------------------

        b11  =  hbox(1,1)
        b12  =  hbox(1,2)
        b13  =  hbox(1,3)
        b22  =  hbox(2,2)
        b23  =  hbox(2,3)
        b33  =  hbox(3,3)
        
        ib11 = ihbox(1,1)
        ib12 = ihbox(1,2)
        ib13 = ihbox(1,3)
        ib22 = ihbox(2,2)
        ib23 = ihbox(2,3)
        ib33 = ihbox(3,3)


        open(UNIT_, file='neighbor_list.txt')
        open(UNIT__, file='all_list.txt')
        write(UNIT_, "(2(I20))") nUpdate
        write(UNIT__, "(2(I20))") nUpdate

        do ig = 1, nGlobalGroup
        
            list => nlist(ig)%n
            maxList = list(0)

            write(UNIT_, "(2(I5))") ig, maxList
            do nl = 1, maxList

                jg = list(nl)
                gxij = grpX(ig) - grpX(jg)                        
                gyij = grpY(ig) - grpY(jg)
                gzij = grpZ(ig) - grpZ(jg)  

                dxij = ib11*gxij + ib12*gyij + ib13*gzij
                dyij =             ib22*gyij + ib23*gzij
                dzij =                         ib33*gzij

                if (dxij >  HALF) dxij = dxij - ONE
                if (dxij < -HALF) dxij = dxij + ONE  

                if (dyij >  HALF) dyij = dyij - ONE
                if (dyij < -HALF) dyij = dyij + ONE

                if (dzij >  HALF) dzij = dzij - ONE
                if (dzij < -HALF) dzij = dzij + ONE 
                  
                pxij = b11*dxij + b12*dyij + b13*dzij
                pyij =            b22*dyij + b23*dzij
                pzij =                       b33*dzij  

                gijsq = pxij*pxij + pyij*pyij + pzij*pzij
                write(UNIT_,"(I3, F10.4)") jg, sqrt(gijsq)
                
            end do ! nl-loop
            write(UNIT_,"()")     

            mol = myMolIndex(ig)
            write(UNIT__, "(2(I5))") ig, nlist(ig)%n(0)
            do jg = 1, nGlobalGroup

                if (myMolIndex(jg) == mol) cycle ! neglect same molecule

                gxij = grpX(ig) - grpX(jg)                        
                gyij = grpY(ig) - grpY(jg)
                gzij = grpZ(ig) - grpZ(jg)  

                dxij = ib11*gxij + ib12*gyij + ib13*gzij
                dyij =             ib22*gyij + ib23*gzij
                dzij =                         ib33*gzij

                if (dxij >  HALF) dxij = dxij - ONE
                if (dxij < -HALF) dxij = dxij + ONE  

                if (dyij >  HALF) dyij = dyij - ONE
                if (dyij < -HALF) dyij = dyij + ONE

                if (dzij >  HALF) dzij = dzij - ONE
                if (dzij < -HALF) dzij = dzij + ONE 
                  
                pxij = b11*dxij + b12*dyij + b13*dzij
                pyij =            b22*dyij + b23*dzij
                pzij =                       b33*dzij  

                gijsq = pxij*pxij + pyij*pyij + pzij*pzij
                write(UNIT__,"(I3, F10.4)") jg, sqrt(gijsq)
                
            end do ! nl-loop
            write(UNIT__,"()") 

        end do ! ig-loop

        close(UNIT_)
        close(UNIT__)

    end subroutine PrintNeighborList




!   -------------------------------------------------------
    subroutine DeallocateNeighborList()
!   -------------------------------------------------------
        implicit none

        integer :: ig
!   -------------------------------------------------------

        do ig = 1, nGlobalGroup 
            if (allocated(nlist)) deallocate(nlist(ig)%n)
        end do ! ig-loop

        if (allocated(nlist)) deallocate(nlist)

        if(allocated(savedGrpX)) deallocate(savedGrpX, savedGrpY, savedGrpZ)

        write(*,"(A)") 'Neighbor list is deallocated'
        write(1,"(A)") 'Neighbor list is deallocated'

    end subroutine DeallocateNeighborList



    

end module mod_neighbor_list