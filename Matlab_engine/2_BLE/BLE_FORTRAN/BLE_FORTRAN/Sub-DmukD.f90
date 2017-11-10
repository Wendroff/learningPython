subroutine DmukD(mub,DmuD)
use mesh
implicit none
real*8,intent(in) :: mub(N+1) 
real*8,intent(out):: DmuD(N+1,N+1)
integer*8 :: j,k,e
    do j = 1,N+1
        do k = 1,N+1
            DmuD(j,k) = 0.0d0
        enddo
    enddo
    do j = 1,N+1
        do k = 1,N+1
            do e = 1,N+1
                DmuD(j,e) = DmuD(j,e) + LN(j,k)*mub(k)*LN(k,e)
            enddo
        enddo
    enddo

    
end subroutine