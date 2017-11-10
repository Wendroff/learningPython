Subroutine my_rmm(A,B,C,l,m,n)
    implicit none
    integer*8,intent(in) :: l,m,n !A的行数，A的列数B的行数，B的列数
    real*8,intent(in)    :: A(l,m),B(m,n)
    real*8,intent(out)   :: C(l,n)
    integer*8 :: i,j,k
    
    c = 0.0d0
    do i = 1,l
        do j = 1,n
            do k = 1,m
                c(i,j) = c(i,j) + a(i,k)*b(k,j)
            enddo
        enddo
    enddo
            
    
    
    
    
end subroutine