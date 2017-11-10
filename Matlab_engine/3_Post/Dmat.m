function [LN,IN]=Dmat(N,z,y,dzy,dyz)
Cj=diag(eye(N+1));
Cj(1)=2;
Cj(N+1)=2;
D = zeros(N+1,N+1);
for j=0:N
    for k=0:N

        
        if j~=k
            
           D(j+1,k+1)=Cj(j+1)*(-1)^(j+k)/Cj(k+1)/(z(j+1)-z(k+1))*dzy(j+1);
           
        elseif (j==k)&&(j>=1)&&(j<=N-1)
                
              D(j+1,k+1)=-z(k+1)/2/(1-z(k+1)^2)*dzy(j+1);
              
        elseif j==0&&k==0
                
              D(j+1,k+1)=-(2*N^2+1)/6*dzy(j+1);
              
        elseif k==N&&j==N 
                
              D(j+1,k+1)=(2*N^2+1)/6*dzy(j+1);
              
        end
     end
end
  

                 
I1 = zeros(N,N+1);      
for i=0:N-1
    for j=0:N
        
        I1(i+1,j+1)=0;
        
        for k=0:N
            
                Ckj=2/Cj(k+1)/Cj(j+1)/N*cos(j*k*pi/N);
            
            if k==0&&i>=0&&i<=N-1
                
                Qik=cos((i+1)*pi/N)-cos(i*pi/N);
                
            elseif k==1&&i>=0&&i<=N-1
                
                Qik=-0.5*( (sin((i+1)*pi/N))^2-(sin(i*pi/N))^2 );
                
            elseif k>=2&&k<=N&&i>=0&&i<=N-1
                
                ri1=(i+1)*pi/N;
                ri=i*pi/N;
                Qik=0.5*( cos(ri1*(k+1))/(k+1)+cos(ri1*(1-k))/(1-k)-cos(ri*(k+1))/(k+1)-cos(ri*(1-k))/(1-k) );
                
            end
                
            I1(i+1,j+1)=I1(i+1,j+1)-Qik*Ckj*dyz(j+1);
        end
    end
end


LN=D;
IN=[zeros(1,N+1);
              I1;];
          
          