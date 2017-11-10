% 用于V4.0边界层程序的后处理
% 
% 
% 依旧只能处理流向均匀网格
% by Wendroff 2015-12-10



clear all
close all
format long
clc

global Pr r Rg ;

global S DENe Ue We Te


% Ma=3


Pr=0.72;
r=1.4;
Rg=287;
%-------NPSE无量纲化参考位置----------------------------
% index_reff = 10; %无量纲化参考位置的流向指标

[N,yi,ymax,~,i0] = A_read_BLEgrid;

% N=300;
vec=(N:-1:0)';
%%%%%%%%%%%%%%%%%%
z=cos(pi*vec/N);

% ymax=15;
% yi=4;
a=ymax*yi/(ymax-2*yi);
b=1+2*a/ymax;
y = zeros(N+1,1);
dzy = zeros(N+1,1);
for j=1:1:N+1;               %
    y(j)=a*(1+z(j))/(b-z(j));
    dzy(j)=(a+a*b)/(y(j)+a)^2;
end
dyz=1./dzy;

[LN,IN]=Dmat(N,z,y,dzy,dyz);
%%%%%%%%%%%%%%%%%%%%%%%%%%
% load iniatial state caculated by initial0.m
% load input/yini.dat
% load input/Vini.dat
% load input/Uini.dat
% load input/Wini.dat
% load input/Tini.dat
%  
% Vi0=Vini;
% Ui0=Uini;
% Wi0=Wini;
% Ti0=Tini;
% clear Vini Uini Wini Tini
load ../input/SOUT2.dat
load ../input/ROUT2.dat
load ../input/UOUT2.dat
load ../input/WOUT2.dat
load ../input/TOUT2.dat
load ../input/XcOUT.dat

load ../output/T_dimensionless.dat
load ../output/U_dimensionless.dat
load ../output/V_dimensionless.dat
load ../output/W_dimensionless.dat
n = length(V_dimensionless)/(N+1);
%只取其中部分点 
index = i0:1:length(SOUT2);
S0 = SOUT2(i0);
S=SOUT2(index) - S0*ones(length(index),1);
Ue=UOUT2(index);
We=WOUT2(index);
DENe=ROUT2(index);
Te=TOUT2(index);
xc=XcOUT(index);
% n = length(S); %number of discret point in x direction
% n = 40;    %用于调试程序
% S    = linspace(0,1,n);
% Ue   = 100*ones(1,n);
% We   = zeros(1,n);
% DENe = 1.21*ones(1,n);
% Te   = 300*ones(1,n);
% 
%      Fai(:,1)=[Vi0;Ui0;Wi0;Ti0;];
% %        
%      for i=2:2
%      
%      [Fai,H,Lbs]=meanorder1(Fai,i,N,LN,IN,y);
%      
%      end
% 
%      for i=3:3
%       
%      [Fai,H,Lbs]=meanorder2(Fai,i,N,LN,IN,y);
%      if (any(isnan(Fai(:,i)))||any(isinf(Fai(:,i))))
%         disp('There are NaNs or Infs in Fai');break;
%      end
%      end
% 
%      for i=4:4
%       
%      [Fai,H,Lbs]=meanorder3(Fai,i,N,LN,IN,y);
%      
%      end
% %      
% %      
%      for i=5:5
%       
%      [Fai,H,Lbs]=meanorder4(Fai,i,N,LN,IN,y);
%      
%      end
%      
%      
%      for i=6:n
%       
%      [Fai,H,Lbs]=meanorder5(Fai,i,N,LN,IN,y);
%      if (any(isnan(Fai(:,i)))||any(isinf(Fai(:,i))))
%         disp('There are NaNs or Infs in Fai');break;
%      end
%      end
%      n = i-1;
% output
%%
%data out put: in x-eta coordinate
% save all
yout  = zeros(N+1,n);
V     = zeros(N+1,n);
U1    = zeros(N+1,n);
W1    = zeros(N+1,n);
T1    = zeros(N+1,n);
%变量代换前的结果
%******************************************************************************************************
fid = fopen('../output/meanflow-dimensionless.dat', 'wt');
fprintf(fid,'VARIABLES = "s", "eta", "V", "U", "W", "T" \n  ZONE I=    %d  J = %d  F=POINT \n',n,(N+1));
for j=1:N+1
    for i=1:n
         fprintf(fid, '%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n', S(i), y(j) ,V_dimensionless((i-1)*(N+1)+j),U_dimensionless((i-1)*(N+1)+j),W_dimensionless((i-1)*(N+1)+j),T_dimensionless((i-1)*(N+1)+j));
%          eta(j,i) = y(j);
         V(j,i)  = V_dimensionless((i-1)*(N+1)+j);
         U1(j,i) = U_dimensionless((i-1)*(N+1)+j);
         W1(j,i) = W_dimensionless((i-1)*(N+1)+j);
         T1(j,i) = T_dimensionless((i-1)*(N+1)+j);
    end
end
fclose(fid);
%*****************************************************************************************************************
% transfer to x-y coordinate
% 

%
mue = zeros(n,1);
Re1 = zeros(n,1);
ReL = zeros(n,1);
L1 = zeros(n,1);
for i=2:n
    
    mue(i)=1.789*10^(-5)*(Te(i)/288)^1.5*(288+110.4)/(Te(i)+110.4);
    Re1(i)=DENe(i)*Ue(i)/mue(i);
    ReL(i)=sqrt(DENe(i)*Ue(i)*S(i)/mue(i));

    L1(i)=sqrt(mue(i)*S(i)/DENe(i)/Ue(i));
    INT=IN*T1(:,i);
      for j=2:1:N+1
          INT(j)=INT(j)+INT(j-1);
      end 
    yout(:,i)=L1(i)*INT; % results
    
%     U2(:,i)=U0(:,i)*Ue(i);
      
end 

dx=S(2)-S(1);
% pause
% U1 V1 W1 T1
%
% U1=U;
% W1=W;
% T1=T;
Ro1=1./T1;
Tx = zeros(N+1,n);
Tetax = zeros(N+1,n);
etax = zeros(N+1,n);
V1 = zeros(N+1,n);
V11= zeros(N+1,n);
for i=7:1:n

    Tx(:,i)=(274/120*T1(:,i)-600/120*T1(:,i-1)+600/120*T1(:,i-2)-400/120*T1(:,i-3)+150/120*T1(:,i-4)-24/120*T1(:,i-5) )/dx;
     
    INTx=IN*Tx(:,i);
    for j=2:1:N+1
        INTx(j)=INTx(j)+INTx(j-1);
    end
      
    dRex1=(274/120*Re1(i)-600/120*Re1(i-1)+600/120*Re1(i-2)-400/120*Re1(i-3)+150/120*Re1(i-4)-24/120*Re1(i-5) )/dx;
    
    Tetax(:,i)=-yout(:,i)/2/L1(i)*(1/S(i)-1/Re1(i)*dRex1)-INTx;
    youtx = (274/120*yout(:,i)-600/120*yout(:,i-1)+600/120*yout(:,i-2)-400/120*yout(:,i-3)+150/120*yout(:,i-4)-24/120*yout(:,i-5) )/dx;
    etax(:,i) = -youtx./(LN*yout(:,i));
%     etax(:,i)=Tetax(:,i)./T1(:,i); 
    V1(:,i)=T1(:,i).*V(:,i)/ReL(i)-1/ReL(i)*S(i)*U1(:,i).*T1(:,i).*etax(:,i);
    
end

Rox1 = zeros(N+1,n);
Ux1 = zeros(N+1,n);
Vx1 = zeros(N+1,n);
Wx1 = zeros(N+1,n);
Tx1 = zeros(N+1,n);
DRo1 = zeros(N+1,n);
DU1 = zeros(N+1,n);
DV1 = zeros(N+1,n);
DW1 = zeros(N+1,n);
DT1 = zeros(N+1,n);
DDRo1 = zeros(N+1,n);
DDU1 = zeros(N+1,n);
DDV1 = zeros(N+1,n);
DDW1 = zeros(N+1,n);
DDT1 = zeros(N+1,n);
for i=7:1:n
    
    Rox1(:,i)=(274/120*Ro1(:,i)-600/120*Ro1(:,i-1)+600/120*Ro1(:,i-2)-400/120*Ro1(:,i-3)+150/120*Ro1(:,i-4)-24/120*Ro1(:,i-5) )/dx;
    Ux1(:,i)=(274/120*U1(:,i)-600/120*U1(:,i-1)+600/120*U1(:,i-2)-400/120*U1(:,i-3)+150/120*U1(:,i-4)-24/120*U1(:,i-5) )/dx;
    Vx1(:,i)=(274/120*V1(:,i)-600/120*V1(:,i-1)+600/120*V1(:,i-2)-400/120*V1(:,i-3)+150/120*V1(:,i-4)-24/120*V1(:,i-5) )/dx;
    Wx1(:,i)=(274/120*W1(:,i)-600/120*W1(:,i-1)+600/120*W1(:,i-2)-400/120*W1(:,i-3)+150/120*W1(:,i-4)-24/120*W1(:,i-5) )/dx;
    Tx1(:,i)=(274/120*T1(:,i)-600/120*T1(:,i-1)+600/120*T1(:,i-2)-400/120*T1(:,i-3)+150/120*T1(:,i-4)-24/120*T1(:,i-5) )/dx;
    
    
    DRo1(:,i)=LN*Ro1(:,i);
    DU1(:,i)=LN*U1(:,i);
    DV1(:,i)=LN*V1(:,i);
    DW1(:,i)=LN*W1(:,i);
    DT1(:,i)=LN*T1(:,i);
    
    DDRo1(:,i)=LN*LN*Ro1(:,i);
    DDU1(:,i)=LN*LN*U1(:,i);
    DDV1(:,i)=LN*LN*V1(:,i);
    DDW1(:,i)=LN*LN*W1(:,i);
    DDT1(:,i)=LN*LN*T1(:,i);
end



%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%      x  y      %%%%
Rox2 = zeros(N+1,n);
Ux2 = zeros(N+1,n);
Vx2 = zeros(N+1,n);
Wx2 = zeros(N+1,n);
Tx2 = zeros(N+1,n);
DRo2 = zeros(N+1,n);
DU2 = zeros(N+1,n);
DV2 = zeros(N+1,n);
DW2 = zeros(N+1,n);
DT2 = zeros(N+1,n);
for i=7:1:n
   
     Rox2(:,i)=Rox1(:,i)+DRo1(:,i).*etax(:,i); %Roxout(:,i)=Rox1(:,i)+DRo1(:,i).*etax(:,i);
     Ux2(:,i)=Ux1(:,i)+DU1(:,i).*etax(:,i);    %Uxout(:,i)=Ux1(:,i)+DU1(:,i).*etax(:,i);
     Vx2(:,i)=Vx1(:,i)+DV1(:,i).*etax(:,i);    %Vxout(:,i)=Vx1(:,i)+DV1(:,i).*etax(:,i);
     Wx2(:,i)=Wx1(:,i)+DW1(:,i).*etax(:,i);    %Wxout(:,i)=Wx1(:,i)+DW1(:,i).*etax(:,i);
     Tx2(:,i)=Tx1(:,i)+DT1(:,i).*etax(:,i);    %Txout(:,i)=Tx1(:,i)+DT1(:,i).*etax(:,i);
    
    
    DRo2(:,i)=Ro1(:,i)/L1(i).*DRo1(:,i); % DRoout(:,i)=Ro1(:,i)/L1(i).*DRo1(:,i);
    DU2(:,i)=Ro1(:,i)/L1(i).*DU1(:,i);   % DUout(:,i)=Ro1(:,i)/L1(i).*DU1(:,i);
    DV2(:,i)=Ro1(:,i)/L1(i).*DV1(:,i);   % DVout(:,i)=Ro1(:,i)/L1(i).*DV1(:,i);
    DW2(:,i)=Ro1(:,i)/L1(i).*DW1(:,i);   % DWout(:,i)=Ro1(:,i)/L1(i).*DW1(:,i);
    DT2(:,i)=Ro1(:,i)/L1(i).*DT1(:,i);   % DTout(:,i)=Ro1(:,i)/L1(i).*DT1(:,i);
    
end


DDRo2 = zeros(N+1,n);
DDU2 = zeros(N+1,n);
DDV2 = zeros(N+1,n);
DDW2 = zeros(N+1,n);
DDT2 = zeros(N+1,n);
for i=7:1:n
    
    
    
    DDRo2(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDRo1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DRo1(:,i); % DDRoout(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDRo1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DRo1(:,i);
    DDU2(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDU1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DU1(:,i);    % DDUout(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDU1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DU1(:,i);
    DDV2(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDV1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DV1(:,i);    % DDVout(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDV1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DV1(:,i);
    DDW2(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDW1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DW1(:,i);    % DDWout(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDW1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DW1(:,i);
    DDT2(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDT1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DT1(:,i);    % DDTout(:,i)=Ro1(:,i)/L1(i).*Ro1(:,i)/L1(i).*DDT1(:,i)+1/L1(i)*Ro1(:,i)/L1(i).*DRo1(:,i).*DT1(:,i);
    
end


%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% x y dimensional %%%
Ro3 = zeros(N+1,n);
U3 = zeros(N+1,n);
V3 = zeros(N+1,n);
W3 = zeros(N+1,n);
T3 = zeros(N+1,n);
Rox3 = zeros(N+1,n);
Ux3 = zeros(N+1,n);
Vx3 = zeros(N+1,n);
Wx3 = zeros(N+1,n);
Tx3 = zeros(N+1,n);
DRo3 = zeros(N+1,n);
DU3 = zeros(N+1,n);
DV3 = zeros(N+1,n);
DW3 = zeros(N+1,n);
DT3 = zeros(N+1,n);
DDRo3 = zeros(N+1,n);
DDU3 = zeros(N+1,n);
DDV3 = zeros(N+1,n);
DDW3 = zeros(N+1,n);
DDT3 = zeros(N+1,n);
for i=7:1:n
    
    Ro3(:,i)=DENe(i)*Ro1(:,i);
    U3(:,i)=Ue(i)*U1(:,i);
    V3(:,i)=Ue(i)*V1(:,i);
    W3(:,i)=We(i)*W1(:,i);
    T3(:,i)=Te(i)*T1(:,i);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dRox=(274/120*DENe(i)-600/120*DENe(i-1)+600/120*DENe(i-2)-400/120*DENe(i-3)+150/120*DENe(i-4)-24/120*DENe(i-5) )/dx;
    Rox3(:,i)=DENe(i)*Rox2(:,i)+Ro1(:,i)*dRox;
    
    dUex=(274/120*Ue(i)-600/120*Ue(i-1)+600/120*Ue(i-2)-400/120*Ue(i-3)+150/120*Ue(i-4)-24/120*Ue(i-5) )/dx;
    Ux3(:,i)=Ue(i)*Ux2(:,i)+U1(:,i)*dUex;
    Vx3(:,i)=Ue(i)*Vx2(:,i)+V1(:,i)*dUex;

    dWex=(274/120*We(i)-600/120*We(i-1)+600/120*We(i-2)-400/120*We(i-3)+150/120*We(i-4)-24/120*We(i-5) )/dx;
    Wx3(:,i)=We(i)*Wx2(:,i)+W1(:,i)*dWex; 

    dTex=(274/120*Te(i)-600/120*Te(i-1)+600/120*Te(i-2)-400/120*Te(i-3)+150/120*Te(i-4)-24/120*Te(i-5) )/dx;
    Tx3(:,i)=Te(i)*Tx2(:,i)+T1(:,i)*dTex;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DRo3(:,i)=DENe(i)*DRo2(:,i);
    DU3(:,i)=Ue(i)*DU2(:,i);
    DV3(:,i)=Ue(i)*DV2(:,i);
    DW3(:,i)=We(i)*DW2(:,i);
    DT3(:,i)=Te(i)*DT2(:,i);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DDRo3(:,i)=DENe(i)*DDRo2(:,i);
    DDU3(:,i)=Ue(i)*DDU2(:,i);
    DDV3(:,i)=Ue(i)*DDV2(:,i);
    DDW3(:,i)=We(i)*DDW2(:,i);
    DDT3(:,i)=Te(i)*DDT2(:,i);
end

DUx3 = zeros(N+1,n);
DVx3 = zeros(N+1,n);
for i = 7:n
    DUx3(:,i)=(274/120*DU3(:,i)-600/120*DU3(:,i-1)+600/120*DU3(:,i-2)-400/120*DU3(:,i-3)+150/120*DU3(:,i-4)-24/120*DU3(:,i-5) )/dx+LN*DU3(:,i).*etax(:,i);
    DVx3(:,i)=(274/120*DV3(:,i)-600/120*DV3(:,i-1)+600/120*DV3(:,i-2)-400/120*DV3(:,i-3)+150/120*DV3(:,i-4)-24/120*DV3(:,i-5) )/dx+LN*DV3(:,i).*etax(:,i);
end

fid = fopen('../output/meanflow.dat', 'wt');
fprintf(fid,'VARIABLES = "s", "y", "Ro","U", "V", "W", "T" \n  ZONE I=    %d  J = %d  F=POINT \n',n-6,(N+1));
for j=1:N+1
    for i=7:n
         fprintf(fid, '%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n', S(i), yout(j,i) ,Ro3(j,i),U3(j,i),V3(j,i),W3(j,i),T3(j,i));
%          
    end
end
fclose(fid);
%
% FileNo = fopen('BASEFLOW/den3.dat','wt');
% 	for i=7:n
%         for j = 1:N+1
%             fprintf(FileNo,'%20.16f %20.16f %20.16f \n',Ro3(j,i),Rox3(j,i),DRo3(j,i));
%         end
% 	end 
% fclose(FileNo);
% 
% FileNo = fopen('BASEFLOW/U3.dat','wt');
% 	for i=7:n
%         for j = 1:N+1
%             fprintf(FileNo,'%20.16f %20.16f %20.16f %20.16f %20.16f \n',U3(j,i),Ux3(j,i),DU3(j,i),DUx3(j,i),DDU3(j,i));
%         end
% 	end 
% fclose(FileNo);
% 
% FileNo = fopen('BASEFLOW/V3.dat','wt');
% 	for i=7:n
%         for j = 1:N+1
%             fprintf(FileNo,'%20.16f %20.16f %20.16f %20.16f %20.16f \n',V3(j,i),Vx3(j,i),DV3(j,i),DVx3(j,i),DDV3(j,i));
%         end
% 	end 
% fclose(FileNo);
% 
% FileNo = fopen('BASEFLOW/W3.dat','wt');
% 	for i=7:n
%         for j = 1:N+1
%             fprintf(FileNo,'%20.16f %20.16f %20.16f %20.16f \n',W3(j,i),Wx3(j,i),DW3(j,i),DDW3(j,i));
%         end
% 	end 
% fclose(FileNo);
% 
% FileNo = fopen('BASEFLOW/T3.dat','wt');
% 	for i=7:n
%         for j = 1:N+1
%             fprintf(FileNo,'%20.16f %20.16f %20.16f %20.16f \n',T3(j,i),Tx3(j,i),DT3(j,i),DDT3(j,i));
%         end
% 	end 
% fclose(FileNo);
%
%差值到NPSE使用的网格上
NN = 301;                   %NPSE法向网格点数
% n_NPSE = n-5;               %NPSE流向网格点数
ymax=100;                    %定义远场位置
yi=5;                      %加密位置
dy=2/(NN-1);               	%Chebyshev区域间距
z=zeros(1,NN);   	%预定义z
dz=zeros(1,NN);   	%预定义z
for i=1:1:NN
    z(i)=-1+(i-1)*dy;
    
end

a=ymax*yi/(ymax-2*yi);    	%Chebyshev区域常数a
b=1+2*a/ymax;            	%Chebyshev区域常数b

ym=zeros(1,NN);    	%预定义ym
for i=1:1:NN;
    ym(i)=a*(1+z(i))/(b-z(i));
    dz(i)=(a+a*b)/(ym(i)+a)^2;
end
%注意此时ym为无量纲

%-----------------------------预定义速度场----------------------------------
U9=zeros(NN,n);Uy9=zeros(NN,n);Uyy9=zeros(NN,n);Ux9=zeros(NN,n);Uxy9=zeros(NN,n);
V9=zeros(NN,n);Vy9=zeros(NN,n);Vyy9=zeros(NN,n);Vx9=zeros(NN,n);Vxy9=zeros(NN,n);
T9=zeros(NN,n);Ty9=zeros(NN,n);Tyy9=zeros(NN,n);Tx9=zeros(NN,n);%Txy9=zeros(N,Nx); 
den9=zeros(NN,n);deny9=zeros(NN,n);denx9=zeros(NN,n);
W9=zeros(NN,n);Wy9=zeros(NN,n);Wyy9=zeros(NN,n);Wx9=zeros(NN,n);

%-------------------------------预定义常数----------------------------------
% infinity=30;                %定义远场（需要与子程序Guess一起改变）
%--------------------------------------------------------------------------

Y_Point=N-9;               %边界层解的流向网格点数(排除边界振荡部分)
extract=1:(N-9);
Y_Point_Ex=1000;           %附加点数
%************************参考量*******************************************
% x0    = S(index_reff);         %参考点的有量纲流向位置，即参考点到驻点的位置
% c      = 1.29;                   %垂直于前缘的弦长
% Ureff  = 19.5410;                %自由来流速度
% Treff  = 293.15;
% Roreff = 1.21;
% mureff = 1.8029E-5;
[Ureff,Treff,Roreff,mureff,c] = A_read_reference;
deta0  = sqrt(mureff*c/Ureff/Roreff);%参考长度
Re0    = Ureff*deta0/mureff;
Mareff = Ureff/sqrt(r*Rg*Treff);
%*************************************************************************

for i = 7:n
    ys=yout(:,i)./deta0;
    y_plus=linspace(ys(Y_Point),300,Y_Point_Ex+1)';            %增加10000个点,远场扩充到300
    Y_E=[ys(extract);y_plus(2:Y_Point_Ex+1)];
        %--将速度场等进行扩展,并无量纲化--
    U_E   =[U3(extract,i);U3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff;
    Uy_E  =[DU3(extract,i);DU3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0;
    Uyy_E =[DDU3(extract,i);DDU3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0^2;
    Ux_E  =[Ux3(extract,i);Ux3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0;
    Uxy_E =[DUx3(extract,i);DUx3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0^2;

%     V_E   =[V3(:,i);V3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff;
    Vy_E  =[DV3(extract,i);DV3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0;
    Vyy_E =[DDV3(extract,i);DDV3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0^2;
%     Vx_E  =[Vx3(:,i);Vx3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0;
    Vxy_E =[DVx3(extract,i);DVx3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0^2;
    
    W_E   =[W3(extract,i);W3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff;
    Wy_E  =[DW3(extract,i);DW3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0;
    Wyy_E =[DDW3(extract,i);DDW3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0^2;
    Wx_E  =[Wx3(extract,i);Wx3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0;
%     Wxy_E =[DWx3(:,i);DWx3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff*deta0^2;

    T_E=[T3(extract,i);T3(Y_Point,i)*ones(Y_Point_Ex,1)]/Treff;
    Ty_E=[DT3(extract,i);DT3(Y_Point,i)*ones(Y_Point_Ex,1)]/Treff*deta0;
    Tyy_E=[DDT3(extract,i);DDT3(Y_Point,i)*ones(Y_Point_Ex,1)]/Treff*deta0^2;
    Tx_E=[Tx3(extract,i);Tx3(Y_Point,i)*ones(Y_Point_Ex,1)]/Treff*deta0;
%     Txy_E=[DTx3;DTx3(Y_Point)*ones(Y_Point_Ex,1)]/Treff*deta0^2;

    Den_E=[Ro3(extract,i);Ro3(Y_Point,i)*ones(Y_Point_Ex,1)]/Roreff;
    Deny_E=[DRo3(extract,i);DRo3(Y_Point,i)*ones(Y_Point_Ex,1)]/Roreff*deta0;
    Denx_E=[Rox3(extract,i);Rox3(Y_Point,i)*ones(Y_Point_Ex,1)]/Roreff*deta0;
    
    U0=interp1(Y_E,U_E,ym,'spline');
    Uy0=interp1(Y_E,Uy_E,ym,'spline');
    Uyy0=interp1(Y_E,Uyy_E,ym,'spline');
    Ux0=interp1(Y_E,Ux_E,ym,'spline');
    Uxy0=interp1(Y_E,Uxy_E,ym,'spline');

    V0=interp1(ys,V3(:,i),ym,'linear','extrap')/Ureff;
    Vy0=interp1(Y_E,Vy_E,ym,'spline');
    Vyy0=interp1(Y_E,Vyy_E,ym,'spline');
%     Vx0=interp1(Y_E,Vx_E,ym,'spline');
    Vx0=interp1(ys,Vx3(:,i),ym,'linear','extrap')/Ureff*deta0;
    Vxy0=interp1(Y_E,Vxy_E,ym,'spline');

    T0=interp1(Y_E,T_E,ym,'spline');
    Ty0=interp1(Y_E,Ty_E,ym,'spline');
    Tyy0=interp1(Y_E,Tyy_E,ym,'spline');
    Tx0=interp1(Y_E,Tx_E,ym,'spline');
    
    W0=interp1(Y_E,W_E,ym,'spline');
    Wy0=interp1(Y_E,Wy_E,ym,'spline');
    Wyy0=interp1(Y_E,Wyy_E,ym,'spline');
    Wx0=interp1(Y_E,Wx_E,ym,'spline');
%     Txy0=interp1(Y_E,Txy_E,ym,'spline');

    Den0=interp1(Y_E,Den_E,ym,'spline');
    Deny0=interp1(Y_E,Deny_E,ym,'spline');
    Denx0=interp1(Y_E,Denx_E,ym,'spline');
    
     %最终结果赋值(待无量纲化)
    U9(:,i)=U0';
    Uy9(:,i)=Uy0';
    Uyy9(:,i)=Uyy0';
    Ux9(:,i)=Ux0';
    Uxy9(:,i)=Uxy0';
    
    V9(:,i)=V0';
    Vy9(:,i)=Vy0';
    Vyy9(:,i)=Vyy0';
    Vx9(:,i)=Vx0';
    Vxy9(:,i)=Vxy0';
    
    T9(:,i)=T0';
    Ty9(:,i)=Ty0';
    Tyy9(:,i)=Tyy0';
    Tx9(:,i)=Tx0';
    
    W9(:,i)=W0';
    Wy9(:,i)=Wy0';
    Wyy9(:,i)=Wyy0';
    Wx9(:,i)=Wx0';
%     Txy9(:,i)=Txy0';
    
    den9(:,i)=Den0';
    deny9(:,i)=Deny0';
    denx9(:,i)=Denx0';
    
end

den4=zeros(NN*n,3);
U4=zeros(NN*n,5);
V4=zeros(NN*n,5);
W4=zeros(NN*n,4);
T4=zeros(NN*n,4);

for i=7:1:n
%     disp(i)
    for j=1:1:NN
        den4(j+(i-1)*NN,:)=[den9(j,i),denx9(j,i),deny9(j,i)];
        U4(j+(i-1)*NN,:)=[U9(j,i),Ux9(j,i),Uy9(j,i),Uxy9(j,i),Uyy9(j,i)];
        V4(j+(i-1)*NN,:)=[V9(j,i),Vx9(j,i),Vy9(j,i),Vxy9(j,i),Vyy9(j,i)];
        T4(j+(i-1)*NN,:)=[T9(j,i),Tx9(j,i),Ty9(j,i),Tyy9(j,i)];
        W4(j+(i-1)*NN,:)=[W9(j,i),Wx9(j,i),Wy9(j,i),Wyy9(j,i)];
    end
end

if (exist('../output/BASEFLOW','dir')~=7)
    mkdir('../output/BASEFLOW');
end
fid1 = fopen('../output/BASEFLOW/den4.dat', 'wt');
fid2 = fopen('../output/BASEFLOW/U4.dat', 'wt');
fid3 = fopen('../output/BASEFLOW/V4.dat', 'wt');
fid4 = fopen('../output/BASEFLOW/W4.dat', 'wt');
fid5 = fopen('../output/BASEFLOW/T4.dat', 'wt');
for i=1:1:NN*n
    fprintf(fid1, '%17.15e %17.15e %17.15e\n', den4(i,:));
    fprintf(fid2, '%17.15e %17.15e %17.15e %17.15e %17.15e\n', U4(i,:));
    fprintf(fid3, '%17.15e %17.15e %17.15e %17.15e %17.15e\n', V4(i,:));
    fprintf(fid4, '%17.15e %17.15e %17.15e %17.15e\n', W4(i,:));
    fprintf(fid5, '%17.15e %17.15e %17.15e %17.15e \n', T4(i,:));
end
fclose(fid1);fclose(fid2);fclose(fid3);fclose(fid4);fclose(fid5);

%流向坐标文件输出（无量纲）
S=S + S0*ones(length(index),1);
% xc = S/c;
sout = S/deta0;
fid = fopen('../output/BASEFLOW/SOUT.dat', 'wt');
for i=1:1:n
    fprintf(fid, '%17.15e\n', sout(i));
end
fclose(fid);

fid = fopen('../output/BASEFLOW/xc.dat', 'wt');
for i=1:1:n
    fprintf(fid, '%17.15e\n', xc(i));
end
fclose(fid);

%----------tecplot output------------------------------

fid = fopen('../output/meanflow-tecout.dat','wt');
fprintf(fid,'VARIABLES = "s","x/c", "y", "DEN", "DENx","DENy","U","Ux","Uy","Uxy","Uyy","V","Vx","Vy","Vxy","Vyy","T","Tx","Ty","Tyy","W","Wx","Wy","Wyy" \n');
fprintf(fid,'ZONE T = DIFFERENCE I = %d J = %d F = POINT \n',NN,n-6);
for i = 7:n
    for j = 1:NN
        fprintf(fid,'%17.15e ',sout(i),xc(i),ym(j),den9(j,i),denx9(j,i),deny9(j,i),U9(j,i),Ux9(j,i),Uy9(j,i),Uxy9(j,i),Uyy9(j,i),V9(j,i),Vx9(j,i),Vy9(j,i),Vxy9(j,i),Vyy9(j,i),T9(j,i),Tx9(j,i),Ty9(j,i),Tyy9(j,i),W9(j,i),Wx9(j,i),Wy9(j,i),Wyy9(j,i));
        fprintf(fid,'\n');
    end
end
fclose(fid);
%测试用输出*****************************************************
fid1 = fopen('test/Utvsy.dat','wt');
fprintf(fid1,'variables = Ut y x/c\n ZONE T = MINE I = %d J = %d \n',NN,n-6);
fid2 = fopen('test/Wtvsy.dat','wt');
fprintf(fid2,'variables = Wt y x/c\n ZONE T = MINE I = %d J = %d \n',NN,n-6);
for i = 7:n
    for j = 1:NN
        theta = atan(W9(NN,i)/U9(NN,i));
        Utmax = sqrt(W9(NN,i)^2+U9(NN,i)^2);
        Ut = (U9(:,i)*cos(theta) + W9(:,i)*sin(theta))/Utmax;
        Wt = (U9(:,i)*sin(theta) - W9(:,i)*cos(theta))/Utmax;
        fprintf(fid1,'%17.15e ',Ut(j),ym(j)*deta0*1000,xc(i));
        fprintf(fid2,'%17.15e ',-Wt(j),ym(j)*deta0*1000,xc(i));
        fprintf(fid1,'\n');
        fprintf(fid2,'\n');
    end
end
fclose(fid1);fclose(fid2);
%生成NPSE控制文件************************************************
% fid_set = fopen('Setting.in','w');
% fprintf(fid_set,'%s\n','#Flow Varibles# ');
% fprintf(fid_set,'%s\n','Pr     r    Rg    Re0    Ma_reff    Te_reff    k0  L_reff(deta0) ');
% fprintf(fid_set,'%f %f %f %f %f %f XXXX %f \n\n',Pr,r,Re0,Rg,Mareff,Treff,deta0);
% fprintf(fid_set,'#Mesh# \n Ny		ymax	yi		iend0\n');
% fprintf(fid_set,'%d %d %d %d\n',NN,ymax,yi,n-5);
% fprintf(fid_set,'\n#Fourier Mode#\n m0		n0		Mode_m?\n\n\n#Calc#\ni0		iend	ndx		Amp\n');
% fprintf(fid_set,'xx %d 1 XXXX\n\n',n-5);
% fprintf(fid_set,'#Instability#\nomega0		beta0');
% fprintf(fid_set,'\n\n\n x0_reff Ureff Treff DENreff deta0 \n %f %f %f %f %f ',x0,Ureff,Treff,Roreff,deta0);
% fclose(fid_set);
%生成LST控制文件************************************************
if (exist('../output/LSToutput','dir')~=7)
    mkdir('../output/LSToutput');
end
fid_set = fopen('../output/LSToutput/input.dat','w');
fprintf(fid_set,' Ureff  Treff  Roreff  mureff  Re0  deta0  beta   omega   i \n');
fprintf(fid_set,'%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f XXXX XXXX XXXX',Ureff,Treff,Roreff,mureff,Re0,deta0);
fclose(fid_set);
%生成LST需要的网格文件******************************************
fid_set = fopen('../output/LSToutput/Setgrid.dat','w');
fprintf(fid_set,'Ny   ymax   yi  iend\n');
fprintf(fid_set,'%d %f %f %d',NN,ymax,yi,n);
fclose(fid_set);
save all;
% 
% 
% 
