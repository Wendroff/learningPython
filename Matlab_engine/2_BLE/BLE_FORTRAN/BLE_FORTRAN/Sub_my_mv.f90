subroutine my_mv(A,b,C,n)
    implicit none
    
    integer*8,intent(in) :: n
    real*8,intent(in),dimension(n,n) :: A
    real*8,intent(in),dimension(n)   :: b
    real*8,intent(out),dimension(n)  :: C
    
    integer*8 i,j,k
    
    
    
    do i = 1,n
        C(i) = 0.0d0
        do j = 1,n
        
            C(i) = C(i) + A(i,j)*B(j)
        
        enddo
    enddo
    
    
end subroutine