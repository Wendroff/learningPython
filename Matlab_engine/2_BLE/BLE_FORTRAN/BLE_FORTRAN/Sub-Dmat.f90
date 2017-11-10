Subroutine Dmat(z,dzy,dyz)
use mesh
implicit none
real*16,intent(in)    :: z(N+1),dzy(N+1),dyz(N+1)
integer*8 ::i,j,k
real*16    :: D(N+1,N+1),I1(N,N+1)
real*16    :: Cj(N+1),Ckj,Qik,ri1,ri!,temp
    
    do j = 2,N
        Cj(j) = 1.0_16
    enddo
    Cj(1)   = 2.0_16
    Cj(N+1) = 2.0_16
    
    do j=0,N
        do k=0,N

        
            if (j/=k) then
                !temp = 2.0_16*sin(pi+real(j+k)/2.0_16/real(N)*pi)*sin(real(j-k)/2.0_16/real(N)*pi)/dzy(j+1)
                if (mod(j+k,2)==0) then
                  D(j+1,k+1)=Cj(j+1)/Cj(k+1)/(z(j+1)-z(k+1))*dzy(j+1)
                  !D(j+1,k+1)=-Cj(j+1)/Cj(k+1)/temp
                else
                  D(j+1,k+1)=-Cj(j+1)/Cj(k+1)/(z(j+1)-z(k+1))*dzy(j+1)
                  !D(j+1,k+1)=Cj(j+1)/Cj(k+1)/temp
                endif
           
            elseif ((j==k).and.(j>=1).and.(j<=N-1)) then
                
                  D(j+1,k+1)=-z(k+1)/2.0_16/(1.0_16-z(k+1)**2)*dzy(j+1) 
              
            elseif ((j==0).and.(k==0)) then
                
                  D(j+1,k+1)=-(2.0_16*N**2+1.0_16)/6.0_16*dzy(j+1) 
              
            elseif ((k==N).and.(j==N))  then
                
                  D(j+1,k+1)=(2*N**2+1.0_16)/6.0_16*dzy(j+1);
              
            endif
         enddo
    enddo
    
    do j = 1,N
        do k = 1,N+1
            I1(j,k) = 0.0_16
        enddo
    enddo
    do j = 1,N+1
        do k = 1,N+1
            INN(j,k) = 0.0_16
        enddo
    enddo
    
    do i=0,N-1
        do j=0,N
        
            I1(i+1,j+1)=0.0_16;
        
            do k=0,N
            
                Ckj=2.0_16/Cj(k+1)/Cj(j+1)/N*cos(real(j*k)*pi/N);
            
                if ((k==0).and.(i>=0).and.(i<=N-1)) then
                
                    Qik=cos(real(i+1)*pi/real(N))-cos(real(i)*pi/real(N));
                
                elseif ((k==1).and.(i>=0).and.(i<=N-1)) then
                
                    Qik=-0.5_16*((sin(real(i+1)*pi/real(N)))**2-(sin(real(i)*pi/real(N)))**2 );
                
                elseif ((k>=2).and.(k<=N).and.(i>=0).and.(i<=N-1)) then
                
                    ri1 = real(i+1)*pi/N;
                    ri  = real(i)*pi/N;
                    Qik = 0.5_16*(cos(ri1*real(k+1))/real(k+1)+cos(ri1*real(1-k))/real(1-k)-cos(ri*real(k+1))/real(k+1)-cos(ri*real(1-k))/real(1-k));
                
                endif
                
                I1(i+1,j+1)  = I1(i+1,j+1)-Qik*Ckj*dyz(j+1);
                INN(i+2,j+1) = I1(i+1,j+1)
            enddo
        enddo
    enddo
    
    do j = 1,N+1
        do k = 1,N+1
            LN(j,k)  = D(j,k)
            INN(j,k) = 0.0_16
        enddo
    enddo
    do j = 1,N
        do k = 1,N+1
            INN(j+1,k) = I1(j,k)
            
        enddo
    enddo
    
    !call Test_output_Dmat
    
  end subroutine