function   [Fai,H,Lbs]=meanorder4(Fai,i,N,LN,IN,y)

     global Pr r Rg T00
     global S DENe Ue We Te


      
     Fai(:,i)=Fai(:,i-1);
     V=Fai(1:N+1,:);
     U=Fai(N+1+1:2*(N+1),:);
     W=Fai(2*(N+1)+1:3*(N+1),:);
     T=Fai(3*(N+1)+1:4*(N+1),:);
     
     dx=S(i)-S(i-1);
     
     [mue,dmuex,dDenex,duex]=Finte(i,dx,DENe,Ue,Te);
     
     beta0=S(i)/mue*dmuex+S(i)/DENe(i)*dDenex+S(i)/Ue(i)*duex;
     
     beta1=S(i)/Ue(i)*duex;
     
     aae=sqrt(r*Rg*Te(i))
     Ue(i)
     Maue=Ue(i)/aae
     Mawe=We(i)/aae
%    pause
    
     a1=110.4/Te(i);
     a2=194/Te(i);
     
     for itr=1:1:100000000
         i 
         itr
         Vt=Fai(1:N+1,i);
         Ut=Fai(N+1+1:2*(N+1),i);
         Wt=Fai(2*(N+1)+1:3*(N+1),i);
         Tt=Fai(3*(N+1)+1:4*(N+1),i);
         DUt=LN*Ut;
         DWt=LN*Wt;
         
         for j=1:1:N+1

             mu=(1+a1)*Tt(j)^1.5/(Tt(j)+a1);
             
             mub(j)=(1+a1)*Tt(j)^1.5/(Tt(j)+a1)/Tt(j);
             
             k1(j)=(1+a2)*Tt(j)^1.5/(Tt(j)+a2);
             
             kb(j)=(1+a2)*Tt(j)^1.5/(Tt(j)+a2)/Tt(j);
         end
        
        %%%%%%%%%%%%%
         % x equation
         D22=S(i)*diag(Ut)*(50/24/dx)+diag(Vt)*LN+diag(beta1)*diag(Ut)*eye(N+1)-LN*diag(mub)*LN; 
         
         D24=-diag(beta1)*eye(N+1);
         
         % z equation
         D33=S(i)*diag(Ut)*(50/24/dx)+diag(Vt)*LN-LN*diag(mub)*LN; 
         
         % E equation
         D42=-(r-1)*Maue^2*diag(mub)*diag(DUt)*LN;
         
         D43=-(r-1)*Mawe^2*diag(mub)*diag(DWt)*LN;
         
         D44=S(i)*diag(Ut)*(50/24/dx)+diag(Vt)*LN-1/Pr*LN*diag(kb)*LN;
         %%%%%%%%%%%%%%
         DV0=eye(N+1);
         for j=2:1:N+1
             DV0(j,j-1)=-1;
         end
        %%%%%%%%%%%%%%%%%
        Lbs=[];
        F=[];
        for j=2:1:N
      
             Lbs=[Lbs;
                  DV0(j,:)       IN(j,:)*(S(i)*(50/24/dx)+0.5*(1+beta0))   zeros(1,N+1)   zeros(1,N+1);
                  zeros(1,N+1)   D22(j,:)                                  zeros(1,N+1)   D24(j,:);
                  zeros(1,N+1)   zeros(1,N+1)                              D33(j,:)       zeros(1,N+1);
                  zeros(1,N+1)   D42(j,:)                                  D43(j,:)       D44(j,:);  ];
             F=[F;
                IN(j,:)*S(i)*( 96/24/dx*U(:,i-1)-72/24/dx*U(:,i-2)+32/24/dx*U(:,i-3)-6/24/dx*U(:,i-4) );
                S(i)*Ut(j)*(   96/24/dx*U(j,i-1)-72/24/dx*U(j,i-2)+32/24/dx*U(j,i-3)-6/24/dx*U(j,i-4) );
                S(i)*Ut(j)*(   96/24/dx*W(j,i-1)-72/24/dx*W(j,i-2)+32/24/dx*W(j,i-3)-6/24/dx*W(j,i-4) );
                S(i)*Ut(j)*(   96/24/dx*T(j,i-1)-72/24/dx*T(j,i-2)+32/24/dx*T(j,i-3)-6/24/dx*T(j,i-4) );];
        end
  
        j=1;  %  wall
        
           Lbs=[1 zeros(1,N)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N+1);
                zeros(1,N+1)  1 zeros(1,N)   zeros(1,N+1)   zeros(1,N+1);
                zeros(1,N+1)  zeros(1,N+1)   1 zeros(1,N)   zeros(1,N+1);
                zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   LN(j,:);
                                                               Lbs;];
           F=[ 0;
               0;
               0;
               0;
               F;];
           
      
   
        j=N+1; % boundary layer
            Lbs=[ Lbs;             % by xu
                  DV0(j,:)      IN(j,:)*(S(i)*(50/24/dx)+0.5*(1+beta0)) zeros(1,N+1)   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N) 1                            zeros(1,N+1)   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N+1)                            zeros(1,N) 1   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N+1)                            zeros(1,N+1)   zeros(1,N) 1;
                  ];
   
            F=[F;
                IN(j,:)*S(i)*( 96/24/dx*U(:,i-1)-72/24/dx*U(:,i-2)+32/24/dx*U(:,i-3)-6/24/dx*U(:,i-4) );
                1;
                1;
                1;];

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         H=[];
         for j=2:1:N
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %  continiue
             RV1=zeros(1,N+1);  RU1=zeros(1,N+1);
             
             RV1(1,j-1)=-1;
             RV1(1,j)=1;
             
             RU1(1,j-1)=0.5*(y(j)-y(j-1))*(S(i)*50/24/dx+0.5*(1+beta0));     % order
             RU1(1,j)=0.5*(y(j)-y(j-1))*(S(i)*50/24/dx+0.5*(1+beta0));       % order
             
             %  x equation
             RU2=zeros(1,N+1);  RT2=zeros(1,N+1);
             
             RU2(1,j)=S(i)*50/24/dx*Ut(j)+beta1*Ut(j);                      % order
            
             RU2(1,j-1)=-Vt(j)/(y(j)-y(j-1));
             RU2(1,j)=RU2(1,j)+Vt(j)/(y(j)-y(j-1));
             
              
             RU2(1,j-1)=RU2(1,j-1)-(mub(j)+mub(j-1))/(y(j+1)-y(j-1))/(y(j)-y(j-1));
             RU2(1,j)=RU2(1,j)...
                      +1/(y(j+1)-y(j-1))*(mub(j+1)+mub(j))/(y(j+1)-y(j))...
                      +1/(y(j+1)-y(j-1))*(mub(j)+mub(j-1))/(y(j)-y(j-1));
                      
             RU2(1,j+1)=-(mub(j+1)+mub(j))/(y(j+1)-y(j-1))/(y(j+1)-y(j));

             
             RT2(1,j)=-beta1;
             
             %%  z equation
             RW3=zeros(1,N+1); 
                        
             RW3(1,j)=S(i)*50/24/dx*Ut(j);                                  % order

             RW3(1,j-1)=-Vt(j)/(y(j)-y(j-1));

             RW3(1,j)=RW3(1,j)+Vt(j)/(y(j)-y(j-1));

             RW3(1,j-1)=RW3(1,j-1)-(mub(j)+mub(j-1))/(y(j+1)-y(j-1))/(y(j)-y(j-1));
             RW3(1,j)=RW3(1,j)...
                      +1/(y(j+1)-y(j-1))*(mub(j+1)+mub(j))/(y(j+1)-y(j))...
                      +1/(y(j+1)-y(j-1))*(mub(j)+mub(j-1))/(y(j)-y(j-1));
             RW3(1,j+1)=-(mub(j+1)+mub(j))/(y(j+1)-y(j-1))/(y(j+1)-y(j));
      
             % E equation
             RU4=zeros(1,N+1);RW4=zeros(1,N+1);RT4=zeros(1,N+1);

             RU4(1,j-1)=-(r-1)*Maue^2*mub(j-1)*DUt(j)/(y(j)-y(j-1));
             RU4(1,j)=(r-1)*Maue^2*mub(j)*DUt(j)/(y(j)-y(j-1));
             
             RW4(1,j-1)=-(r-1)*Mawe^2*mub(j-1)*DWt(j)/(y(j)-y(j-1));
             RW4(1,j)=(r-1)*Mawe^2*mub(j)*DWt(j)/(y(j)-y(j-1));
           

             RT4(1,j)=S(i)*50/24/dx*Ut(j);                                  % order
             
             RT4(1,j-1)=-Vt(j)/(y(j)-y(j-1));
             RT4(1,j)=RT4(1,j)+Vt(j)/(y(j)-y(j-1));

             
             RT4(1,j-1)=RT4(1,j-1)-1/Pr*(kb(j)+kb(j-1))/(y(j+1)-y(j-1))/(y(j)-y(j-1));
             RT4(1,j)=RT4(1,j)...
                      +1/Pr*1/(y(j+1)-y(j-1))*(kb(j+1)+kb(j))/(y(j+1)-y(j))...
                      +1/Pr*1/(y(j+1)-y(j-1))*(kb(j)+kb(j-1))/(y(j)-y(j-1));
             RT4(1,j+1)=-1/Pr*(kb(j+1)+kb(j-1))/(y(j+1)-y(j-1))/(y(j+1)-y(j));
      
      
             H=[H;
                RV1          RU1            zeros(1,N+1) zeros(1,N+1);
                zeros(1,N+1) RU2            zeros(1,N+1) RT2;
                zeros(1,N+1) zeros(1,N+1)   RW3          zeros(1,N+1);
                zeros(1,N+1) RU4            RW4          RT4;];
         end
          j=1;
          H=[1 zeros(1,N)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N+1);
             zeros(1,N+1)  1 zeros(1,N)   zeros(1,N+1)   zeros(1,N+1);
             zeros(1,N+1)  zeros(1,N+1)   1 zeros(1,N)   zeros(1,N+1);
             zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   -1/y(2) 1/y(2) zeros(1,N-1);
                                                                   H;];

          j=N+1;
  
            RV1=zeros(1,N+1);
            RU1=zeros(1,N+1);

            RV1(1,j-1)=-1;
            RV1(1,j)=1;
            RU1(1,j-1)=0.5*(y(j)-y(j-1))*(S(i)*(50/24/dx)+0.5*(1+beta0));
            RU1(1,j)=0.5*(y(j)-y(j-1))*(S(i)*(50/24/dx)+0.5*(1+beta0));
      
            H=[H;
               RV1           RU1            zeros(1,N+1)   zeros(1,N+1);
               zeros(1,N+1)  zeros(1,N) 1   zeros(1,N+1)   zeros(1,N+1);
               zeros(1,N+1)  zeros(1,N+1)   zeros(1,N) 1   zeros(1,N+1);
               zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N) 1;
               ];

          omega=0.1;
%           (abs(eig(inv(H)*Lbs)))
%           pause
          Faixin=Fai(:,i)+omega*inv(H)*(F-Lbs*Fai(:,i));
  
          max(abs(F-Lbs*Fai(:,i)))
          if max(abs(F-Lbs*Fai(:,i)))<1e-8
             break
          end 
          Fai(:,i)=Faixin;
         
    end
    
    
    