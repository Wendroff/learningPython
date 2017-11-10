Subroutine initial_mesh!(method_flag)
    use mesh
    implicit none
    !integer*8,intent(in)  :: method_flag
    integer*8 :: i,j
    real*16,dimension(:),allocatable :: z,dzy,dyz
    real*16  :: ymax,yi,a,b
    
    write(*,*) 'initializing mesh......'
    
    open(10,file='input/BLEgrid.dat')
    read(10,*)
    read(10,*) N,yi,ymax,i0,iend
    close(10)
    write(*,*) 'N   =',N
    write(*,*) 'yi  =',yi
    write(*,*) 'ymax=',ymax
    write(*,*) 'i0  =',i0
    
    allocate(y(N+1),z(N+1),dyz(N+1),dzy(N+1),LN(N+1,N+1),INN(N+1,N+1))
    do i = 1,N+1
        z(i)=cos(pi*real(N+1-i)/real(N))
        !z(2)=cos(pi*1.0_16)
    enddo
    a = ymax*yi/(ymax-2.0_16*yi);
    b = 1.0_16+2.0_16*a/ymax;
    
    do j=1,N+1;               
        y(j)   = a*(1.0_16+z(j))/(b-z(j));
        dzy(j) = (a+a*b)/(y(j)+a)**2;
        dyz(j) = 1.0_16/dzy(j)
    enddo
    
    ![LN,IN]=Dmat(N,z,y,dzy,dyz);
    call Dmat(z,dzy,dyz)
    !call Dmat_H
    
    
    
end Subroutine