Subroutine isFaiNaN(i,flag)
use Phy
use mesh
implicit none
integer*8,intent(in) :: i
logical,intent(out)  :: flag
integer*8 :: j,Njb
    Njb = 4*N+4
    flag = .false.
    do j = 1,Njb
        if (isnan(Fai(j,i))) then
            flag = .true.
            exit
        endif
    enddo
    
end subroutine