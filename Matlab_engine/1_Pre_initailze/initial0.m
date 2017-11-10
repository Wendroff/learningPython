

%  2010 1 18  by  xugl
%  2015/12/5 Wendroff 进行适当修改，用于V4.0版边界层计算程序
clc
clear all
close all
format long

global Pr r Mawe Maue Te;
global infinity




Pr=0.72;
r=1.4;
Rg=287;

[N,yi,ymax,yini,i0] = A_read_BLEgrid;
load ../input/UOUT2.dat
load ../input/WOUT2.dat
load ../input/TOUT2.dat
Ue = UOUT2(i0);
We = WOUT2(i0);
Te = TOUT2(i0);
clear UOUT2 WOUT2 TOUT2
% Te=293.14999999999998000 ;
%    2.731098743318289e+02
% Ue= 0.01002579609534197 ;
Maue=Ue/sqrt(r*Rg*Te);
% We= 0;
Mawe=We/sqrt(r*Rg*Te);
infinity=20;

num=500;

sinit=bvpinit(linspace(0,infinity,100),@guess);
sol=bvp4c(@bla1,@blbc,sinit);
t0=linspace(0,infinity,num);
% y=bvpval(sol,t0)';
y=deval(sol,t0)';



y0=t0';
U00=y(:,1);
V00=y(:,2);
W00=y(:,3);
T00=y(:,4);



Vini=interp1(y0,V00,yini,'spline');
Uini=interp1(y0,U00,yini,'spline');
Wini=interp1(y0,W00,yini,'spline');
Tini=interp1(y0,T00,yini,'spline');
% xxx(:,1)=t0;
% xxx(:,2)=U00;
% xxx(:,3)=T00;
%%
% plot(t0,T00);
plot(T00,t0);
%%
% fid = fopen('yini.dat', 'wt');
%     for j=1:1:num
%          fprintf(fid, ' %20.16f\n', yini(j));
%     end
% fclose(fid);

fid = fopen('../input/Uini.dat', 'wt');
    for j=1:1:N+1
         fprintf(fid, ' %20.16f\n', Uini(j));
    end
fclose(fid);

fid = fopen('../input/Vini.dat', 'wt');
    for j=1:1:N+1
         fprintf(fid, ' %20.16f\n', Vini(j));
    end
fclose(fid);

fid = fopen('../input/Wini.dat', 'wt');
    for j=1:1:N+1
         fprintf(fid, ' %20.16f\n', Wini(j));
    end
fclose(fid);


fid = fopen('../input/Tini.dat', 'wt');
    for j=1:1:N+1
         fprintf(fid, ' %20.16f\n', Tini(j));
    end
fclose(fid);

%%
subplot(2,2,1);
plot(Uini,yini,'linewidth',4);
title('Uini','FontName','Times New Roman','FontSize',20);
ylabel('Y','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);
grid on
subplot(2,2,2);
plot(Vini,yini,'linewidth',4);
title('Vini','FontName','Times New Roman','FontSize',20);
ylabel('Y','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);
grid on
subplot(2,2,3);
plot(Wini,yini,'linewidth',4);
title('Wini','FontName','Times New Roman','FontSize',20);
ylabel('Y','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);
grid on
subplot(2,2,4);
plot(Tini,yini,'linewidth',4);
title('Tini','FontName','Times New Roman','FontSize',20);
ylabel('Y','FontName','Times New Roman','FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);
grid on