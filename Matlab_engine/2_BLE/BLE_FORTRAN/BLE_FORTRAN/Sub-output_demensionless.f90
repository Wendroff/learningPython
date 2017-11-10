Subroutine output_demensionless(i)
use Phy
use mesh
implicit none
    integer*8,intent(in) :: i
    integer*8 :: j
    j = 1
    Fai(        j,i) = 0.0d0
    Fai(   N+1 +j,i) = 0.0d0
    Fai(2*(N+1)+j,i) = 0.0d0
    Fai(3*(N+1)+j,i) = 0.0d0
        
    do j = 1,N+1
        write(11,100) Fai(        j,i)
        write(12,100) Fai(   N+1 +j,i)
        write(13,100) Fai(2*(N+1)+j,i)
        write(14,100) Fai(3*(N+1)+j,i)
    enddo
100 FORMAT (E20.10)
end subroutine