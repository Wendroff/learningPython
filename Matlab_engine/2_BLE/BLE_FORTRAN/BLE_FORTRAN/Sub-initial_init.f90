Subroutine initial_init
    use Phy
    use mesh
    implicit none
    integer*8 :: i,k
    
    write(*,*) 'initializing inlet condition......'
    
    allocate(Fai(4*N+4,iend-i0+1))
    k=1;Fai = 0.0d0;
    k=1;Fai = 0.0d0;
    open(9,file='input/Vini.dat')
    do i = 1,N+1
        read(9,*) Fai(k,1)
        k = k+1
    enddo
    close(9)
    
    
    open(9,file='input/Uini.dat')
    do i = 1,N+1
        read(9,*) Fai(k,1)
        k = k+1
    enddo
    close(9)
    
    open(9,file='input/Wini.dat')
    do i = 1,N+1
        read(9,*) Fai(k,1)
        k = k+1
    enddo
    close(9)
    
    open(9,file='input/Tini.dat')
    do i = 1,N+1
        read(9,*) Fai(k,1)
        k = k+1
    enddo
    close(9)
end subroutine