function   [Fai]=mean2order(Fai,i,N,LN,IN,y)

global Pr r Rg T00;


global S DENe Ue We Te


      
     Fai(:,i)=Fai(:,i-1);
     V=Fai(1:N+1,:);
     U=Fai(N+1+1:2*(N+1),:);
     W=Fai(2*(N+1)+1:3*(N+1),:);
     T=Fai(3*(N+1)+1:4*(N+1),:);
     
     dx=S(i)-S(i-1);
     
     [mue,dmuex,dDenex,duex]=Finte(i,dx,DENe,Ue,Te)
     
     
     beta0=S(i)/mue*dmuex+S(i)/DENe(i)*dDenex+S(i)/Ue(i)*duex;
     
     beta1=S(i)/Ue(i)*duex;
     
     aae=sqrt(r*Rg*Te(i))
     Ue(i)
     Maue=Ue(i)/aae
     pause
     
     Mawe=We(i)/aae;
     
    
     a1=110.4/T00;
     a2=194/T00; 
     
     for itr=1:1:100000
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

         end
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        L0=eye(num+1);

         % x equation
         D22=S(i)/dx*diag(Ut)+diag(Vt)*LN+diag(beta1)*diag(Ut)*eye(N+1)-LN*diag(mub)*LN; 
         
         D24=-diag(beta1)*eye(N+1);
         
         % z equation
         D33=S(i)/dx*diag(Ut)+diag(Vt)*LN-LN*diag(mub)*LN; 
         
         % E equation
         D42=-(r-1)*Maue^2*diag(mub)*diag(DUt)*LN;
         
         D43=-(r-1)*Mawe^2*diag(mub)*diag(DWt)*LN;
         
         D44=S(i)/dx*diag(Ut)+diag(Vt)*LN-1/Pr*LN*diag(mub)*LN;
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
                  DV0(j,:)       IN(j,:)*(S(i)/dx+0.5*(1+beta0))   zeros(1,N+1)   zeros(1,N+1);
                  zeros(1,N+1)   D22(j,:)                          zeros(1,N+1)   D24(j,:);
                  zeros(1,N+1)   zeros(1,N+1)                      D33(j,:)       zeros(1,N+1);
                  zeros(1,N+1)   D42(j,:)                          D43(j,:)       D44(j,:);  ];
             F=[F;
                IN(j,:)*S(i)/dx*U(:,i-1);
                S(i)/dx*Ut(j)*U(j,i-1);
                S(i)/dx*Ut(j)*W(j,i-1);
                S(i)/dx*Ut(j)*T(j,i-1);];
        end
  
        j=1;  %  wall
        
           Lbs=[Lbs;
                1 zeros(1,N)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N+1);
                zeros(1,N+1)  1 zeros(1,N)   zeros(1,N+1)   zeros(1,N+1);
                zeros(1,N+1)  zeros(1,N+1)   1 zeros(1,N)   zeros(1,N+1);
                zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   LN(j,:);];
           F=[ 0;
               0;
               0;
               0;
               F;];
           
      
   
        j=N+1; % boundary layer
            Lbs=[ DV0(j,:)      IN(j,:)*(S(i)/dx+0.5*(1+beta0)) zeros(1,N+1)   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N) 1                    zeros(1,N+1)   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N+1)                    zeros(1,N) 1   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N+1)                    zeros(1,N+1)   zeros(1,N) 1;
                  Lbs;];
   
            F=[F;
                IN(j,:)*S(i)/dx*U(:,i-1);
                1;
                1;
                1;];
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         H=[];
         for j=2:1:N
             
             % initial
             % continiue
             RV1=zeros(1,N+1);  RU1=zeros(1,N+1);
            
             % x equation
    
             RU2=zeros(1,N+1);  RT2=zeros(1,N+1);
      
             % z equation
    
             RW3=zeros(1,N+1); 
             
             % E equation
             RU4=zeros(1,N+1);RW4=zeros(1,N+1);RT4=zeros(1,N+1);
             
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %  continiue
             RV1(1,j-1)=-1;
             RV1(1,j)=1;
             
             RU1(1,j-1)=0.5*(y(j)-y(j-1))*(S(i)/dx+0.5*(1+beta0));
             RU1(1,j)=0.5*(y(j)-y(j-1))*(S(i)/dx+0.5*(1+beta0));
             
             %  x equation
             RU2(1,j)=S(i)/dx*Ut(j)+beta1*Ut(j);
            
             RU2(1,j-1)=-Vt(j)/(y(j)-y(j-1));
             RU2(1,j)=RU2(1,j)+Vt(j)/(y(j)-y(j-1));
             
              
             RU2(1,j-1)=RU2(1,j-1)-(mub(j)+mub(j-1))/(y(j+1)-y(j-1))/(y(j)-y(j-1));
             RU2(1,j)=RU2(1,j)...
                      +1/(y(j+1)-y(j-1))*(mub(j+1)+mub(j))/(y(j+1)-y(j))...
                      +1/(y(j+1)-y(j-1))*(mub(j)+mub(j-1))/(y(j)-y(j-1));
                      
             RU2(1,j+1)=-(mub(j+1)+mub(j))/(y(j+1)-y(j-1))/(y(j+1)-y(j));

             
             RT2(1,j)=-beta1;
             
             %%  z equation
                        
             RW3(1,j)=S(i)/dx*Ut(j);

             RW3(1,j-1)=-Vt(j)/(y(j)-y(j-1));

             RW3(1,j)=RW3(1,j)+Vt(j)/(y(j)-y(j-1));

             RW3(1,j-1)=RW3(1,j-1)-(mub(j)+mub(j-1))/(y(j+1)-y(j-1))/(y(j)-y(j-1));
             RW3(1,j)=RW3(1,j)...
                      +1/(y(j+1)-y(j-1))*(mub(j+1)+mub(j))/(y(j+1)-y(j))...
                      +1/(y(j+1)-y(j-1))*(mub(j)+mub(j-1))/(y(j)-y(j-1));
             RW3(1,j+1)=-(mub(j+1)+mub(j))/(y(j+1)-y(j-1))/(y(j+1)-y(j));
      
             %%%%%%%%%%%%%%%%%%%%%%%%
             RU4(1,j-1)=-(r-1)*Maue^2*mub(j-1)*DUt(j)/(y(j)-y(j-1));
             RU4(1,j)=(r-1)*Maue^2*mub(j)*DUt(j)/(y(j)-y(j-1));
             
             RW4(1,j-1)=-(r-1)*Mawe^2*mub(j-1)*DWt(j)/(y(j)-y(j-1));
             RW4(1,j)=(r-1)*Mawe^2*mub(j)*DWt(j)/(y(j)-y(j-1));
           

             RT4(1,j)=S(i)/dx*Ut(j);
             
             RT4(1,j-1)=-Vt(j)/(y(j)-y(j-1));
             RT4(1,j)=RT4(1,j)+Vt(j)/(y(j)-y(j-1));

             
             RT4(1,j-1)=RT4(1,j-1)-1/Pr*(mub(j)+mub(j-1))/(y(j+1)-y(j-1))/(y(j)-y(j-1));
             RT4(1,j)=RT4(1,j)...
                      +1/Pr*1/(y(j+1)-y(j-1))*(mub(j+1)+mub(j))/(y(j+1)-y(j))...
                      +1/Pr*1/(y(j+1)-y(j-1))*(mub(j)+mub(j-1))/(y(j)-y(j-1));
             RT4(1,j+1)=-1/Pr*(mub(j+1)+mub(j-1))/(y(j+1)-y(j-1))/(y(j+1)-y(j));
      
      
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
             zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   -1 1 zeros(1,N-1);
                                                                   H;];

          j=N+1;
  
            RV1=zeros(1,N+1);
            RU1=zeros(1,N+1);

            RV1(1,j-1)=-1;
            RV1(1,j)=1;
            RU1(1,j-1)=0.5*(y(j)-y(j-1))*(S(i)/dx+0.5*(1+beta0));
            RU1(1,j)=0.5*(y(j)-y(j-1))*(S(i)/dx+0.5*(1+beta0));
      
            H=[RV1        RU1       zeros(1,N+1)   zeros(1,N+1);
               zeros(1,N+1)  zeros(1,N) 1   zeros(1,N+1)   zeros(1,N+1);
               zeros(1,N+1)  zeros(1,N+1)   zeros(1,N) 1   zeros(1,N+1);
               zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N) 1;
               H;];

          omega=0.2;
          abs(eig(inv(H)*Lbs))
          pause
          Faixin=Fai(:,i)+omega*inv(H)*(F-Lbs*Fai(:,i));
  
          max(abs(F-Lbs*Fai(:,i)))
          if max(abs(F-Lbs*Fai(:,i)))<1e-10
             break
          end 
          Fai(:,i)=Faixin;
          
          
     end