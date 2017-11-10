subroutine GenCSR3(N)
use M_CSR3
    integer*8,intent(in)    :: N
    !integer,allocatable   :: IA(:),JA(:),IndexA(:)
    !integer,intent(out)   :: num
    integer :: Jb(4*N+4,4*N+4)
    integer :: e,j
    
    Njb = 4*N+4
    num = (N+1)*(N-1)*8 + 6*(N-1) + 8 + (N+1)*2
    ALLOCATE(IA(Njb+1),JA(num),IndexA(num),a(num))!
    
    
    do j = 1,4*N+4
        do e = 1,4*N+4
            Jb(j,e) = 0
        enddo
    enddo
    
    do e = 2,N
        !DJb/DVt
        Jb(4*e-3,e)   = 1
        Jb(4*e-3,e-1) = 1
        Jb(4*e-2,e)   = 1
        Jb(4*e-1,e)   = 1
        Jb(4*e  ,e)   = 1
        !DJb/DUt
        Jb(4*e-2,e+N+1) = 1
        Jb(4*e-1,e+N+1) = 1
        Jb(4*e  ,e+N+1) = 1
        !DJb/DWt
        Jb(4*e-1,e+2*N+2) = 1
        !DJb/DTt
        Jb(4*e-2,e+3*N+3) = 1
        Jb(4*e  ,e+3*N+3) = 1
    enddo
        
        
    do e = 1,N+1
            
        !DJb/DUt
            
        do j = 2,N
            Jb(4*j-3,e+N+1) = 1
            Jb(4*j-2,e+N+1) = 1
            Jb(4*j  ,e+N+1) = 1
            
        !DJb/DWt
            
            
            Jb(4*j-1,e+2*N+2) = 1
            Jb(4*j  ,e+2*N+2) = 1
            
        !DJb/DTt
            
            
            Jb(4*j-2,e+3*N+3) = 1
            Jb(4*j-1,e+3*N+3) = 1
            Jb(4*j  ,e+3*N+3) = 1
        enddo
    enddo
    !±ß½çÌõ¼þ
    Jb(1,1) = 1
    Jb(2,N+2) = 1
    Jb(3,2*N+3) = 1
    do e = 1,N+1
        Jb(4,e+3*N+3) = 1
    enddo
        
    Jb(4*N+1,N+1) = 1; Jb(4*N+1,N) = 1
    Jb(4*N+2,2*N+2) = 1; Jb(4*N+3,3*N+3) = 1; Jb(4*N+4,4*N+4) = 1;
    do e = 1,N+1
        Jb(4*N+1,e+N+1) = 1
    enddo
    k = 0
    IA= 0
    JA= 0
    IndexA=0
    IA(1) = 1
    do i = 1,Njb
        IA(i+1) = IA(i)
        do j = 1,Njb
            if (Jb(i,j)==1) then
                k = k+1
                IA(i+1)  = IA(i+1)+1
                IndexA(k)= i
                JA(k)    = j
            endif
            !k = k+M(i,j)
        enddo
    enddo
    write(*,*) 'k  =',k
    
end subroutine