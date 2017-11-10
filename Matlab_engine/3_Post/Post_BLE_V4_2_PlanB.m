% 用于V4.0边界层程序的后处理
% 此版本采用先插值在求导数的方式，精度可能有所降低
% 目前仅仅用作验证导数计算的第二套方案
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
%读入计算结果*************************************
load ../input/SOUT2.dat
load ../input/ROUT2.dat
load ../input/UOUT2.dat
load ../input/WOUT2.dat
load ../input/TOUT2.dat

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
% V11= zeros(N+1,n);
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


Ro3 = zeros(N+1,n);
U3 = zeros(N+1,n);
V3 = zeros(N+1,n);
W3 = zeros(N+1,n);
T3 = zeros(N+1,n);

for i=7:1:n
    
    Ro3(:,i)=DENe(i)*Ro1(:,i);
    U3(:,i)=Ue(i)*U1(:,i);
    V3(:,i)=Ue(i)*V1(:,i);
    W3(:,i)=We(i)*W1(:,i);
    T3(:,i)=Te(i)*T1(:,i);
    

end



fid = fopen('../output/meanflow-dimensional.dat', 'wt');
fprintf(fid,'VARIABLES = "s", "y", "Ro","U", "V", "W", "T" \n  ZONE I=    %d  J = %d  F=POINT \n',n-6,(N+1));
for j=1:N+1
    for i=7:n
         fprintf(fid, '%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f\n', S(i), yout(j,i) ,Ro3(j,i),U3(j,i),V3(j,i),W3(j,i),T3(j,i));
%          
    end
end
fclose(fid);

%差值到NPSE使用的网格上
NN = 301;                   %NPSE法向网格点数
% n_NPSE = n-5;             %NPSE流向网格点数
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
DINm = Dmat_dis(ym);
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
[Ureff,Treff,Roreff,mureff,c] = A_read_reference;
deta0  = sqrt(mureff*c/Ureff/Roreff);%参考长度
S=S + S0*ones(length(index),1);
xc = S/c;
sout = S/deta0; %流向无量纲位置
dx   = sout(2) - sout(1);

Re0    = Ureff*deta0/mureff;
Mareff = Ureff/sqrt(r*Rg*Treff);
%*************************************************************************

for i = 7:n
    ys=yout(:,i)./deta0;
    y_plus=linspace(ys(Y_Point),300,Y_Point_Ex+1)';            %增加10000个点,远场扩充到300
    Y_E=[ys(extract);y_plus(2:Y_Point_Ex+1)];
        %--将速度场等进行扩展,并无量纲化--
    U_E   =[U3(extract,i);U3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff;
    W_E   =[W3(extract,i);W3(Y_Point,i)*ones(Y_Point_Ex,1)]/Ureff;
    T_E=[T3(extract,i);T3(Y_Point,i)*ones(Y_Point_Ex,1)]/Treff;
    Den_E=[Ro3(extract,i);Ro3(Y_Point,i)*ones(Y_Point_Ex,1)]/Roreff;
    
    U0=interp1(Y_E,U_E,ym,'spline');
    V0=interp1(ys,V3(:,i),ym,'linear','extrap')/Ureff;
    T0=interp1(Y_E,T_E,ym,'spline');
    W0=interp1(Y_E,W_E,ym,'spline');
    Den0=interp1(Y_E,Den_E,ym,'spline');
    
    
     %最终结果赋值
    U9(:,i)   =U0';
    Uy9(:,i)  =DINm*U0';
    Uyy9(:,i) =DINm*DINm*U0';
    
    
    V9(:,i)   =V0';
    Vy9(:,i)  =DINm*V0';
    Vyy9(:,i) =DINm*Vy9(:,i);
    
    
    T9(:,i)   =T0';
    Ty9(:,i)  =DINm*T0';
    Tyy9(:,i) =DINm*Ty9(:,i);
    
    
    W9(:,i)   =W0';
    Wy9(:,i)  =DINm*W0';
    Wyy9(:,i) =DINm*Wy9(:,i);
    
%     Txy9(:,i)=Txy0';
    
    den9(:,i) =Den0';
    deny9(:,i)=DINm*Den0';
    
    
end

for i = 12:n
    Ux9(:,i)  =(274/120*U9(:,i)-600/120*U9(:,i-1)+600/120*U9(:,i-2)-400/120*U9(:,i-3)+150/120*U9(:,i-4)-24/120*U9(:,i-5) )/dx;
    Uxy9(:,i) =(274/120*Uy9(:,i)-600/120*Uy9(:,i-1)+600/120*Uy9(:,i-2)-400/120*Uy9(:,i-3)+150/120*Uy9(:,i-4)-24/120*Uy9(:,i-5) )/dx;
    Vx9(:,i)  =(274/120*V9(:,i)-600/120*V9(:,i-1)+600/120*V9(:,i-2)-400/120*V9(:,i-3)+150/120*V9(:,i-4)-24/120*V9(:,i-5) )/dx;
    Vxy9(:,i) =(274/120*Vy9(:,i)-600/120*Vy9(:,i-1)+600/120*Vy9(:,i-2)-400/120*Vy9(:,i-3)+150/120*Vy9(:,i-4)-24/120*Vy9(:,i-5) )/dx;
    Tx9(:,i)  =(274/120*T9(:,i)-600/120*T9(:,i-1)+600/120*T9(:,i-2)-400/120*T9(:,i-3)+150/120*T9(:,i-4)-24/120*T9(:,i-5) )/dx;
    Wx9(:,i)  =(274/120*W9(:,i)-600/120*W9(:,i-1)+600/120*W9(:,i-2)-400/120*W9(:,i-3)+150/120*W9(:,i-4)-24/120*W9(:,i-5) )/dx;
    denx9(:,i)=(274/120*den9(:,i)-600/120*den9(:,i-1)+600/120*den9(:,i-2)-400/120*den9(:,i-3)+150/120*den9(:,i-4)-24/120*den9(:,i-5) )/dx;
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
% fid1 = fopen('test/Utvsy.dat','wt');
% fprintf(fid1,'variables = Ut y x/c\n ZONE T = MINE I = %d J = %d \n',NN,n-6);
% fid2 = fopen('test/Wtvsy.dat','wt');
% fprintf(fid2,'variables = Wt y x/c\n ZONE T = MINE I = %d J = %d \n',NN,n-6);
% for i = 7:n
%     for j = 1:NN
%         theta = atan(W9(NN,i)/U9(NN,i));
%         Utmax = sqrt(W9(NN,i)^2+U9(NN,i)^2);
%         Ut = (U9(:,i)*cos(theta) + W9(:,i)*sin(theta))/Utmax;
%         Wt = (U9(:,i)*sin(theta) - W9(:,i)*cos(theta))/Utmax;
%         fprintf(fid1,'%17.15e ',Ut(j),ym(j)*deta0*1000,xc(i));
%         fprintf(fid2,'%17.15e ',Wt(j),ym(j)*deta0*1000,xc(i));
%         fprintf(fid1,'\n');
%         fprintf(fid2,'\n');
%     end
% end
% fclose(fid1);fclose(fid2);
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
fprintf(fid_set,'Ny   ymax   yi  \n');
fprintf(fid_set,'%d %f %f',NN,ymax,yi);
fclose(fid_set);
save all;
% 
% 
% 
