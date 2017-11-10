function   [Fai,H,Lbs]=meanorder5(Fai,i,N,LN,IN,y)

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
     
     for itr=1:1:10000
         disp(i); 
         disp(itr);
         Vt=Fai(1:N+1,i);
         Ut=Fai(N+1+1:2*(N+1),i);
         Wt=Fai(2*(N+1)+1:3*(N+1),i);
         Tt=Fai(3*(N+1)+1:4*(N+1),i);
         DUt=LN*Ut;
         DWt=LN*Wt;
         mub = ((1+a1)*ones(N+1,1)).*Tt.^1.5./(Tt+a1*ones(N+1,1))./Tt;
         kb  = ((1+a2)*ones(N+1,1)).*Tt.^1.5./(Tt+a2*ones(N+1,1))./Tt;
%          for j=1:1:N+1
% 
% %              mu=(1+a1)*Tt(j)^1.5/(Tt(j)+a1);
%              
%              mub(j)=(1+a1)*Tt(j)^1.5/(Tt(j)+a1)/Tt(j);
%              
% %              k1(j)=(1+a2)*Tt(j)^1.5/(Tt(j)+a2);
%              
%              kb(j)=(1+a2)*Tt(j)^1.5/(Tt(j)+a2)/Tt(j);
%          end
        
        %%%%%%%%%%%%%
         % x equation
         D22=S(i)*diag(Ut)*(274/120/dx)+diag(Vt)*LN+diag(beta1)*diag(Ut)*eye(N+1)-LN*diag(mub)*LN;   % order
         
         D24=-diag(beta1)*eye(N+1);
         
         % z equation
         D33=S(i)*diag(Ut)*(274/120/dx)+diag(Vt)*LN-LN*diag(mub)*LN; % order 
         
         % E equation
         D42=-(r-1)*Maue^2*diag(mub)*diag(DUt)*LN;
         
         D43=-(r-1)*Mawe^2*diag(mub)*diag(DWt)*LN;
         
         D44=S(i)*diag(Ut)*(274/120/dx)+diag(Vt)*LN-1/Pr*LN*diag(kb)*LN; % order
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
                  DV0(j,:)       IN(j,:)*(S(i)*(274/120/dx)+0.5*(1+beta0))   zeros(1,N+1)   zeros(1,N+1); % order
                  zeros(1,N+1)   D22(j,:)                                    zeros(1,N+1)   D24(j,:);
                  zeros(1,N+1)   zeros(1,N+1)                                D33(j,:)       zeros(1,N+1);
                  zeros(1,N+1)   D42(j,:)                                    D43(j,:)       D44(j,:);  ];
             F=[F;
                IN(j,:)*S(i)*( 600/120/dx*U(:,i-1)-600/120/dx*U(:,i-2)+400/120/dx*U(:,i-3)-150/120/dx*U(:,i-4)+24/120/dx*U(:,i-5) );
                S(i)*Ut(j)*(   600/120/dx*U(j,i-1)-600/120/dx*U(j,i-2)+400/120/dx*U(j,i-3)-150/120/dx*U(j,i-4)+24/120/dx*U(j,i-5) );
                S(i)*Ut(j)*(   600/120/dx*W(j,i-1)-600/120/dx*W(j,i-2)+400/120/dx*W(j,i-3)-150/120/dx*W(j,i-4)+24/120/dx*W(j,i-5) );
                S(i)*Ut(j)*(   600/120/dx*T(j,i-1)-600/120/dx*T(j,i-2)+400/120/dx*T(j,i-3)-150/120/dx*T(j,i-4)+24/120/dx*T(j,i-5) );];
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
                  DV0(j,:)      IN(j,:)*(S(i)*(274/120/dx)+0.5*(1+beta0)) zeros(1,N+1)   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N) 1                              zeros(1,N+1)   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N+1)                              zeros(1,N) 1   zeros(1,N+1);
                  zeros(1,N+1)  zeros(1,N+1)                              zeros(1,N+1)   zeros(1,N) 1;
                  ];
   
            F=[F;
                IN(j,:)*S(i)*( 600/120/dx*U(:,i-1)-600/120/dx*U(:,i-2)+400/120/dx*U(:,i-3)-150/120/dx*U(:,i-4)+24/120/dx*U(:,i-5) );
                1;
                1;
                1;];

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         H=zeros(4*N+4,4*N+4);
         for j=2:1:N
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %  continiue
             RV1=zeros(1,N+1);  RU1=zeros(1,N+1);
             
             RV1(1,j-1)=-1;
             RV1(1,j)=1;
             
             RU1(1,j-1)=0.5*(y(j)-y(j-1))*(S(i)*274/120/dx+0.5*(1+beta0));  % order
             RU1(1,j)=0.5*(y(j)-y(j-1))*(S(i)*274/120/dx+0.5*(1+beta0));    % order
             
             %  x equation
             RU2=zeros(1,N+1);  RT2=zeros(1,N+1);
             
             RU2(1,j)=S(i)*274/120/dx*Ut(j)+beta1*Ut(j);                    % order
            
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
                        
             RW3(1,j)=S(i)*274/120/dx*Ut(j);                                % order

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
           

             RT4(1,j)=S(i)*274/120/dx*Ut(j);                                % order
             
             RT4(1,j-1)=-Vt(j)/(y(j)-y(j-1));
             RT4(1,j)=RT4(1,j)+Vt(j)/(y(j)-y(j-1));

             
             RT4(1,j-1)=RT4(1,j-1)-1/Pr*(kb(j)+kb(j-1))/(y(j+1)-y(j-1))/(y(j)-y(j-1));
             RT4(1,j)=RT4(1,j)...
                      +1/Pr*1/(y(j+1)-y(j-1))*(kb(j+1)+kb(j))/(y(j+1)-y(j))...
                      +1/Pr*1/(y(j+1)-y(j-1))*(kb(j)+kb(j-1))/(y(j)-y(j-1));
             RT4(1,j+1)=-1/Pr*(kb(j+1)+kb(j-1))/(y(j+1)-y(j-1))/(y(j+1)-y(j));
      
      
%              H=[H;
%                 RV1          RU1            zeros(1,N+1) zeros(1,N+1);
%                 zeros(1,N+1) RU2            zeros(1,N+1) RT2;
%                 zeros(1,N+1) zeros(1,N+1)   RW3          zeros(1,N+1);
%                 zeros(1,N+1) RU4            RW4          RT4;];
            H(((j*4-3):(j*4)),:) = [RV1          RU1            zeros(1,N+1) zeros(1,N+1);
                                    zeros(1,N+1) RU2            zeros(1,N+1) RT2;
                                    zeros(1,N+1) zeros(1,N+1)   RW3          zeros(1,N+1);
                                    zeros(1,N+1) RU4            RW4          RT4;];
         end
          j=1;
%           H=[1 zeros(1,N)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N+1);
%              zeros(1,N+1)  1 zeros(1,N)   zeros(1,N+1)   zeros(1,N+1);
%              zeros(1,N+1)  zeros(1,N+1)   1 zeros(1,N)   zeros(1,N+1);
%              zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   -1/y(2) 1/y(2) zeros(1,N-1);
%                                                                    H;];
          H(((j*4-3):(j*4)),:) = [1 zeros(1,N)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N+1);
                                  zeros(1,N+1)  1 zeros(1,N)   zeros(1,N+1)   zeros(1,N+1);
                                  zeros(1,N+1)  zeros(1,N+1)   1 zeros(1,N)   zeros(1,N+1);
                                  zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   -1/y(2) 1/y(2) zeros(1,N-1)];

          j=N+1;
  
            RV1=zeros(1,N+1);
            RU1=zeros(1,N+1);

            RV1(1,j-1)=-1;
            RV1(1,j)=1;
            RU1(1,j-1)=0.5*(y(j)-y(j-1))*(S(i)*(274/120/dx)+0.5*(1+beta0));         % order
            RU1(1,j)=0.5*(y(j)-y(j-1))*(S(i)*(274/120/dx)+0.5*(1+beta0));           % order
      
%             H=[H;
%                RV1           RU1            zeros(1,N+1)   zeros(1,N+1);
%                zeros(1,N+1)  zeros(1,N) 1   zeros(1,N+1)   zeros(1,N+1);
%                zeros(1,N+1)  zeros(1,N+1)   zeros(1,N) 1   zeros(1,N+1);
%                zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N) 1;
%                ];

            H(((j*4-3):(j*4)),:) = [RV1           RU1            zeros(1,N+1)   zeros(1,N+1);
                                   zeros(1,N+1)  zeros(1,N) 1   zeros(1,N+1)   zeros(1,N+1);
                                   zeros(1,N+1)  zeros(1,N+1)   zeros(1,N) 1   zeros(1,N+1);
                                   zeros(1,N+1)  zeros(1,N+1)   zeros(1,N+1)   zeros(1,N) 1;
                                   ];
          omega=0.15;
%           (abs(eig(inv(H)*Lbs)))
%           pause
          Faixin=Fai(:,i)+omega*(H\(F-Lbs*Fai(:,i)));
  
          max(abs(F-Lbs*Fai(:,i)))
          if max(abs(F-Lbs*Fai(:,i)))<1e-8
             break
          end 
          Fai(:,i)=Faixin;
          if (any(isnan(Faixin))||any(isinf(Faixin)))
              disp('There are NaNs or Infs in Fai');break;
          end
    end
    
    
    