module mod_character_tool

    implicit none
    
    private
    
    integer, parameter :: ASCII_a = 97,  ASCII_aa = 65, ASCII_0 = 48, &
                          ASCII_z = 122, ASCII_zz = 90
                          
    public :: IntToCha
    
    contains
    
    function IntToCha(num)
    
        implicit none
        
        integer, parameter    :: MAX = 2
        integer, intent(in)   :: num
        integer               :: i, count_, dnum
        integer, dimension(MAX) :: saveInt
        character(len=MAX)    :: IntToCha
        
        dnum    = num
        count_  = 0
        saveInt = 0

        do while (dnum > 0)
        
            count_          = count_ + 1
            saveInt(count_) = dnum - (dnum/10)*10
            dnum            = (dnum/10)

        end do

        IntToCha = char(ASCII_0 + saveInt(1))

        do i = 2, Max
        
            IntToCha = char(ASCII_0 + saveInt(i))//IntToCha
            
        end do
    
    end function IntToCha
    
    !************************************

end module mod_character_tool
