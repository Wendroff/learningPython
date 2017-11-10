Subroutine output_tecout
use physicalvalue
use mesh
use OUTFLOW, only:S0
implicit none
integer*8 :: i,j
    open(11,file='output\tecout.dat')
    write(11,*) 'VARIABLES = s y Ro U V W T'
    write(11,*) 'ZONE I=',iend-i0-5,'J=',N+1
    do j = 1,N+1
    do i = 7,iend-i0+1
         write(11,100) Xout(i)+S0,Yout(i,j),DEN1(i,j),U1(i,j),V1(i,j),W1(i,j),T1(i,j)!спа©╦ы
    enddo
    enddo
    close(11)
    deallocate(U1,V1,W1,DEN1,T1,Xout,Yout)
100 FORMAT (F15.8,1X,8F20.13)
end subroutine