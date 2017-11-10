Subroutine initial_outflow
    use OUTFLOW
    !use Referencevalue
    use mesh
    implicit none
    integer*8 :: i
    !real*8    :: S0
    
    write(*,*) 'initialzing the boundary condition at the edge of boundary layer'
    
    open(9,file='input/SOUT2.dat')
    !read(9,*) iend
    write(*,*) 'iend=',iend
    pause
    allocate(S(iend),DENe(iend),Ue(iend),We(iend),Te(iend))
    do i = 1,iend
        read(9,*) S(i)
        !S(i) = S(i)
    enddo
    close(9)
    S0 = S(i0)
    write(*,*) 'S0=',S0
    write(*,*) 'Moving the original point to the attachemt point'
    
    do i = 1,iend
        S(i) = S(i) - S0 !以驻点为原点
    enddo
    open(9,file='input/UOUT2.dat')
    do i = 1,iend
        read(9,*) Ue(i)
        !Ue(i) = Ue(i)
    enddo
    close(9)
    WRITE(*,"(1X,A8,F9.6)") "Max Ue=",maxval(Ue)
        
    open(9,file='input/WOUT2.dat')
    do i = 1,iend
        read(9,*) We(i)
        !We(i) = We(i)
    enddo
    close(9)
    WRITE(*,"(1X,A8,F9.6)") "Max We=",maxval(We)
    
    open(9,file='input/ROUT2.dat')
    do i = 1,iend
        read(9,*) DENe(i)
        !DENe(i) = DENe(i)
    enddo
    close(9)
    WRITE(*,"(1X,A8,F9.6)") "Max DENe=",maxval(DENe)
    
    open(9,file='input/TOUT2.dat')
    do i = 1,iend
        read(9,*) Te(i)
        !Te(i) = Te(i)
    enddo
    close(9)
    WRITE(*,"(1X,A8,F9.6)") "Max Te=",maxval(Te)
    
    
end subroutine