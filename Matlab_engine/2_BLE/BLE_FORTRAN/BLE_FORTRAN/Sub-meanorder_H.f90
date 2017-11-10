Subroutine meanorder(i)
    use Phy
    use mesh
    use OUTFLOW
    use Referencevalue
    
    implicit none
    integer*8  :: i !计算Fai（i）
    !integer*8  :: try !=1计算失败
    !real*8  :: duex1
    integer*8  :: itr !迭代步数
    real*8     :: dx,mue,dmuex,dDenex,duex,beta,omega,res,res1!,res_jb
    real*8     :: Lbsu(4*N+4),F(4*N+4),dL(4*N+4),dphi(4*N+4),Jb(4*N+4,4*N+4),Jb_real(4*N+4,4*N+4)
    real*8     :: phi(4*N+4)
    real*8     :: k0,k1
    real*8     :: beta0,beta1,aae,Maue,Mawe,a1,a2
    integer    :: j,k!,L_GS
    integer    :: Njb!,L,IL,M
    !real*8     :: V(N+1),U(N+1),W(N+1),T(N+1)
    !
    res1 = 1.0
    Njb = 4*N+4
    dx=S(i0+i-1)-S(i0+i-2)
    if(i==2) then
        call Finte1order(i,dx,mue,dmuex,dDenex,duex)
    elseif (i==3) then
        call Finte2order(i,dx,mue,dmuex,dDenex,duex)
    elseif (i==4) then
        call Finte3order(i,dx,mue,dmuex,dDenex,duex)
    elseif (i==5) then
        call Finte4order(i,dx,mue,dmuex,dDenex,duex)
    elseif (i>5) then
        call Finte5order(i,dx,mue,dmuex,dDenex,duex)
    else
        write(*,*) 'meanorder_H do not fit i=',i
    endif
    
    beta0=S(i0+i-1)/mue*dmuex+S(i0+i-1)/DENe(i0+i-1)*dDenex+S(i0+i-1)/Ue(i0+i-1)*duex
     
    beta1=S(i0+i-1)/Ue(i0+i-1)*duex;
     
    aae=sqrt(r*Rg*Te(i0+i-1))
    
    Maue=Ue(i0+i-1)/aae
    Mawe=We(i0+i-1)/aae
    a1  =110.4/Te(i0+i-1)
    a2  =194.0/Te(i0+i-1)
    !k1  = -Ue(i0+i-2)-2.0*S(i0+i-2)*duex1
    !k0  = -Ue(i0+i-1)-2.0*S(i0+i-1)*duex
    !
    !a1  = Ue(i0+i-1)/Ue(i0+i-2)
    !a2  = We(i0+i-1)/We(i0+i-2)
    !a3  = Te(i0+i-1)/Te(i0+i-2)
    do j = 1,N+1
        Fai(j*4-3,i)=Fai(j*4-3,i-1)! + y(j)*(k0-k1)
        Fai(j*4-2,i)=Fai(j*4-2,i-1)!*a1
        Fai(j*4-1,i)=Fai(j*4-1,i-1)!*a2
        Fai(j*4  ,i)=Fai(j*4  ,i-1)!*a3
        !V=Fai(1:N+1,:);
        !U=Fai(N+1+1:2*(N+1),:);
        !W=Fai(2*(N+1)+1:3*(N+1),:);
        !T=Fai(3*(N+1)+1:4*(N+1),:);
        !phi(j) = Fai(j,i-1)
        phi(j*4-3)=Fai(j*4-3,i-1)
        phi(j*4-2)=Fai(j*4-2,i-1)!*a1
        phi(j*4-1)=Fai(j*4-1,i-1)!*a2
        phi(j*4  )=Fai(j*4  ,i-1)!*a3
    enddo
    !if (try==0) then
        omega = 0.5
    !else
    !    omega = 0.01
    !endif
     ![mue,dmuex,dDenex,duex]=Finte(i,dx,DENe,Ue,Te);
    

    !if (Compress_flag) then
    !    beta=S(i0+i-1)/r/Ma/Ma/Pe*dPex
    !else
    !    beta=-Ue(i0+i-1)/Te(i0+i-1)*S(i0+i-1)*duex
    !endif
    
    do itr = 1,50
        write(*,*) 'i=',i,'itr=',itr
        call Lbs_H(i,phi,Lbsu,F,Jb,beta0,beta1,a1,a2,Maue,Mawe,dx)
        !call real_Jb(i,phi,Lbsu,Jb_real,beta0,beta1,a1,a2,Maue,Mawe,dx)
        !open(12,file='output/Jb.dat')
        do j = 1,Njb
            dL(j) = (F(j) - Lbsu(j))*omega
            !do k = 1,Njb
            !    Jb(j,k) = Jb(j,k) - Jb_real(j,k)
            !    write(12,*) j,k,Jb(j,k),Jb_real(j,k)
            !enddo
        enddo
        !close(12)
        res = maxval(abs(dL))/omega
        !res_jb = maxval(abs(Jb))
        write(*,*) 'res=',res
        !write(*,*) 'res_jb=',res_jb
        if (res1 < res) then
            if (omega>0.01) then
                omega = omega/2.0
                write(*,*) 'omega=',omega
            endif
        elseif ((res1 > res) .and. (res<1.0)) then
            if (omega<0.8) then
                omega = omega*1.1
                write(*,*) 'omega=',omega
            endif
        endif
        res1 = res
        if (res<1e-9) then
            do j = 1,Njb
                Fai(j,i) = phi(j)
                
            enddo
            !try = 0
            !duex1 = duex
            return
        endif
        
        !call Jb1order(i,Jb,dx,beta)
        call My_pardiso(Jb,dL,Njb)
        !call AGAUS(Jb,dL,Njb,dphi,L_GS,JS_GS)
        !call AGGJE(Jb,Njb,dL,L_GS,JS_GS)
        !L = 6
        !IL= 2*L+1
        !M = 1
        !call ABAND(H,DL,NJB,L,IL,M,L_GS)
        !if (L_GS == 0) then
        !    write(*,*) 'Jb is singular'
        !    try = 1
        !    !pause
        !    return
        !endif
        do j = 1,Njb
            phi(j) = phi(j) + dL(j)
        enddo
        !do j = 1,Njb
        !    phi(j) = Faixin(j) 
        !enddo
    enddo
    !call F1order(i,phi,F,dx)
    !call Test_Jb(i,phi,beta,dx,Jb)
    !try = -1
    write(*,*) 'res is not below 1e-8'
    do j = 1,Njb
        Fai(j,i) = phi(j)
    enddo
    
    
end subroutine
    
    
    
subroutine Lbs_H(i,phii,Lbsu,F,Jb,beta0,beta1,a1,a2,Maue,Mawe,dx)
use Phy
use mesh
use OUTFLOW
use Referencevalue
!use HH
implicit none
integer*8,intent(in) :: i !Fai的指标
real*8,intent(in)    :: phii(4*N+4)
real*8,intent(in)    :: beta0,beta1,a1,a2,Maue,Mawe,dx !该值仅与x有关，提前算好
real*8,intent(out)   :: Lbsu(4*N+4),F(4*N+4)!,HHH(4*N+4,13)!方程左/右端计算出来的向量
!real*8 :: H(4*N+4,4*N+4)
real*8,intent(out) :: Jb(4*N+4,4*N+4)
real*8 :: Vt(N+1),Ut(N+1),Wt(N+1),Tt(N+1),DUt(N+1),DWt(N+1),DTt(N+1),INNUt(N+1)
real*8 :: mu,k,mub(N+1),kb(N+1),muDUt(N+1),muDWt(N+1),DmuDUt(N+1),DmuDWt(N+1),kDTt(N+1),DkDTt(N+1),DmubDt(N+1),DkbDT(N+1)
real*8,allocatable :: U1(:),W1(:),T1(:)
real*8 ::U_b(N+1),W_b(N+1),T_b(N+1),INNU_b(N+1)!上游位置的物理量和差分
real*8 :: DmuD(N+1,N+1),DkD(N+1,N+1)
real*8 :: t0  
integer*8  :: j,e
integer*8  :: Njb,l,n_b
!计算Fai（：，i）的LBS（U） 注意，输出的时候为[V0,U0,W0,T0,V1,U1...WN,TN]
    if (i==2) then
        t0  = 1.0_8/dx
        n_b = 1
    elseif (i==3) then
        t0  = 3.0_8/2.0_8/dx
        n_b = 2
    elseif (i==4) then
        t0  = 11.0_8/6.0_8/dx
        n_b = 3
    elseif (i==5) then
        t0  = 50.0_8/24.0_8/dx
        n_b = 4
    elseif (i>5) then
        t0  = 274.0_8/120.0_8/dx
        n_b = 5
    else
        write(*,*) 'meanorder_H do not fit i=',i
    endif
        allocate(U1(n_b),W1(n_b),T1(n_b))
        do j = 1,N+1
            Vt(j) = phii(j)
            Ut(j) = phii(N+1+j)
            Wt(j) = phii(2*(N+1)+j)
            Tt(j) = phii(3*(N+1)+j)
            do e = 1,n_b
                U1(e) = Fai(N+1+j    ,i-e)
                W1(e) = Fai(2*(N+1)+j,i-e)
                T1(e) = Fai(3*(N+1)+j,i-e)
            enddo
            if (i==2) then
                U_b(j) = -U1(1)/dx
                W_b(j) = -W1(1)/dx
                T_b(j) = -T1(1)/dx
            elseif (i==3) then
                call dx2orderC(U1,dx,U_b(j))
                call dx2orderC(W1,dx,W_b(j))
                call dx2orderC(T1,dx,T_b(j))
            elseif (i==4) then
                call dx3orderC(U1,dx,U_b(j))
                call dx3orderC(W1,dx,W_b(j))
                call dx3orderC(T1,dx,T_b(j))
            elseif (i==5) then
                call dx4orderC(U1,dx,U_b(j))
                call dx4orderC(W1,dx,W_b(j))
                call dx4orderC(T1,dx,T_b(j))
            elseif (i>5) then
                call dx5orderC(U1,dx,U_b(j))
                call dx5orderC(W1,dx,W_b(j))
                call dx5orderC(T1,dx,T_b(j))
            else
                write(*,*) 'meanorder_H do not fit i=',i
            endif
            
        enddo
        do j = 1,4*N+4
            Lbsu(j)  = 0.0_8
            F(j)     = 0.0_8
        enddo
        call my_mv(INN,U_b,INNU_b,N+1)
        call my_mv(LN,Ut,DUt,N+1)!DUt=LN*Ut;
        call my_mv(LN,Wt,DWt,N+1)!DWt=LN*Wt;
        call my_mv(LN,Tt,DTt,N+1)!DWt=LN*Wt;
        call my_mv(INN,Ut,INNUt,N+1)
        
        do j=1,N+1
            !call sutherland(mu,Tt(j))!mu=(1+a1)*Tt(j)^1.5/(Tt(j)+a1);
            mub(j)   = (1.0_8+a1)*Tt(j)**1.5_8/(Tt(j)+a1)/Tt(j)
            DmubDT(j)= (1.0_8+a1)*(a1-Tt(j))*0.5_8/(a1+Tt(j))**2/Tt(j)**0.5_8
            muDUt(j) = mub(j)*DUt(j)
            muDWt(j) = mub(j)*DWt(j)
            !call sutherlandK(k,Tt(j))!k1(j)=(1+a2)*Tt(j)^1.5/(Tt(j)+a2);
            kb(j)    = (1.0_8+a2)*Tt(j)**1.5_8/(Tt(j)+a2)/Tt(j)
            DkbDT(j)= (1.0_8+a2)*(a2-Tt(j))*0.5_8/(a2+Tt(j))**2/Tt(j)**0.5_8
            kDTt(j)  = Kb(j)*DTt(j)
        enddo
        call my_mv(LN,muDUt,DmuDUt,N+1)
        call my_mv(LN,muDWt,DmuDWt,N+1)
        call my_mv(LN,kDTt,DkDTt,N+1)
        
        do j = 2,N
            Lbsu(j*4-3) = Vt(j) - Vt(j-1) + (S(i0+i-1)*t0+0.5_8*(1.0_8+beta0))*INNUt(j)
            Lbsu(j*4-2) = S(i0+i-1)*Ut(j)*(t0*Ut(j)+U_b(j)) + Vt(j)*DUt(j) - DmuDUt(j) - beta1*(Tt(j)-Ut(j)*Ut(j))
            Lbsu(j*4-1) = S(i0+i-1)*Ut(j)*(t0*Wt(j)+W_b(j)) + Vt(j)*DWt(j) - DmuDWt(j)
            Lbsu(j*4  ) = S(i0+i-1)*Ut(j)*(t0*Tt(j)+T_b(j)) + Vt(j)*DTt(j) - DkDTt(j)/Pr - (r-1.0_8)*Maue**2*mub(j)*DUt(j)**2 - (r-1.0_8)*Mawe**2*mub(j)*DWt(j)**2
            F(j*4-3)    = -S(i0+i-1)*INNU_b(j)
        enddo
        
        Lbsu(1)    = Vt(1)
        Lbsu(2)    = Ut(1)
        Lbsu(3)    = Wt(1)
        Lbsu(4)    = DTt(1)
        Lbsu(4*N+1)= Vt(N+1) - Vt(N) + (S(i0+i-1)*t0+0.5d0*(1.0d0+beta0))*INNUt(N+1)
        Lbsu(4*N+2)= Ut(N+1)
        Lbsu(4*N+3)= Wt(N+1)
        Lbsu(4*N+4)= Tt(N+1)
        
        j = N+1
        F(j*4-3) = -S(i0+i-1)*INNU_b(j)
        F(j*4-2) = 1.0_8
        F(j*4-1) = 1.0_8
        F(j*4  ) = 1.0_8
!计算H矩阵
        Njb = 4*N+4
        !do j = 1,Njb
        !    do e = 1,13
        !        HHH(j,e) = 0.0
        !    enddo
        !enddo
        !do j = 1,4*N+4
        !    do e = 1,4*N+4
        !        H(j,e) = 0.0
        !    enddo
        !enddo
        !call my_mv(LN_H,Ut,DUt,N+1)!DUt=LN_H*Ut;
        !call my_mv(LN_H,Wt,DWt,N+1)!DWt=LN_H*Wt;
        !call my_mv(LN_H,Tt,DTt,N+1)!DWt=LN_H*Wt;
        !
        !do j=1,N+1
        !    call sutherland(mu,Tt(j))!mu=(1+a1)*Tt(j)^1.5/(Tt(j)+a1);
        !    !call Fun_DmubDT(DmubDT(j),Tt(j))
        !    mub(j)   = mu/Tt(j);
        !    !call Test_Fun_DmubDt(Tt(j))
        !    !muDUt(j) = mub(j)*DUt(j)
        !    !muDWt(j) = mub(j)*DWt(j)
        !    call sutherlandK(k,Tt(j))!k1(j)=(1+a2)*Tt(j)^1.5/(Tt(j)+a2);
        !    !call Fun_DkbDT(DkbDT(j),Tt(j))
        !    kb(j)    = k*Tt(j);
        !    !kDTt(j)  = Kb(j)*DTt(j)
        !enddo
        !call DmukD_H(mub,DmuD)
        !call DmukD_H(Kb,DkD)
        !
        !do e = 2,N
        !    !DH/DVt
        !    H(4*e-3,4*e-3) = 1.0
        !    H(4*e-3,4*e-7) = -1.0
        !    H(4*e-2,4*e-3) = DUt(e)
        !    H(4*e-1,4*e-3) = DWt(e)
        !    H(4*e  ,4*e-3) = DTt(e)
        !    !DH/DUt
        !    H(4*e-2,4*e-2) = 2.0*S(i0+i-1)*(t0*Ut(e)+U_b(e)) + 2.0*S(i0+i-1)*t0*Ut(e)
        !    H(4*e-1,4*e-2) = 2.0*S(i0+i-1)*(t0*Wt(e)+W_b(e))
        !    H(4*e  ,4*e-2) = 2.0*S(i0+i-1)*(t0*Tt(e)+T_b(e)) - 2.0*(r-1.0)*Ma*Ma*beta*Tt(e)
        !    !DH/DWt
        !    H(4*e-1,4*e-1) = 2.0*S(i0+i-1)*t0*Ut(e)
        !    !DH/DTt
        !    H(4*e-2,4*e  ) = 2.0*beta
        !    H(4*e  ,4*e  ) = 2.0*S(i0+i-1)*t0*Ut(e) - 2.0*(r-1.0)*Ma*Ma*beta*Ut(e) - (r-1.0)*Ma*Ma*DmubDT(e)*(DUt(e)**2+DWt(e)**2)
        !enddo
        !
        !
        !do e = 1,N+1
        !    
        !    !DH/DUt
        !    
        !    do j = 2,N
        !        H(4*j-3,4*e-2) = H(4*j-3,4*e-2) + INN_H(j,e)*(2.0*S(i0+i-1)*t0 + 1.0)
        !        H(4*j-2,4*e-2) = H(4*j-2,4*e-2) + Vt(j)*LN_H(j,e) - DmuD(j,e)
        !        H(4*j  ,4*e-2) = H(4*j  ,4*e-2) - 2.0*(r-1.0)*ma*ma*mub(j)*DUt(j)*LN_H(j,e)
        !    
        !    !DH/DWt
        !    
        !    
        !        H(4*j-1,4*e-1) = H(4*j-1,4*e-1) + Vt(j)*LN_H(j,e) - DmuD(j,e)
        !        H(4*j  ,4*e-1) = - 2.0*(r-1.0)*ma*ma*mub(j)*DWt(j)*LN_H(j,e)
        !    
        !    !DH/DTt
        !    
        !    
        !        H(4*j-2,4*e  ) = H(4*j-2,4*e  ) - LN_H(j,e)*DmubDT(e)*DUt(e)
        !        H(4*j-1,4*e  ) = H(4*j-1,4*e  ) - LN_H(j,e)*DmubDT(e)*DWt(e)
        !        H(4*j  ,4*e  ) = H(4*j  ,4*e  ) + Vt(j)*LN_H(j,e) - DkD(j,e)/Pr - (LN_H(j,e)*DkbDT(e)*DTt(e))/Pr 
        !    enddo
        !enddo
        !!边界条件
        !H(1,1) = 1.0;
        !H(2,2) = 1.0;
        !H(3,3) = 1.0; 
        !do e = 1,N+1
        !    H(4,4*e  ) = LN_H(1,e)
        !enddo
        !
        !H(4*N+1,4*N+1) = 1.0;H(4*N+1,4*N-3) = -1.0;
        !H(4*N+2,4*N+2) = 1.0; H(4*N+3,4*N+3) = 1.0; H(4*N+4,4*N+4) = 1.0;
        !do e = 1,N+1
        !    H(4*N+1,4*e-2) = INN_H(N+1,e)*(2.0*S(i0+i-1)*t0 + 1.0)
        !enddo
        !
        !!测试
        !!if (i==5) then
        !    !call Test_output_M(H,Njb)
        !!endif
        !!****
        !
        !l = 6
        !do e = 1,Njb
        !    if(e<=l) then
        !        do j = 1,(e+l)
        !            HHH(e,j) = H(e,j)
        !        enddo
        !    elseif((Njb-e)<=l) then
        !        do j = e-l,Njb
        !            HHH(e,j-e+l+1) = H(e,j)
        !        enddo
        !    else
        !        do j = e-l,e+l
        !            HHH(e,j-e+l+1) = H(e,j)
        !        enddo
        !    endif
        !enddo
!计算Jb矩阵
        do j = 1,4*N+4
            do e = 1,4*N+4
                Jb(j,e) = 0.0_8
            enddo
        enddo
        !do j=1,N+1
        !    call sutherland(mu,Tt(j))!mu=(1+a1)*Tt(j)^1.5/(Tt(j)+a1);
        !    call Fun_DmubDT(DmubDT(j),Tt(j))
        !    mub(j)   = mu/Tt(j);
        !    !call Test_Fun_DmubDt(Tt(j))
        !    !muDUt(j) = mub(j)*DUt(j)
        !    !muDWt(j) = mub(j)*DWt(j)
        !    call sutherlandK(k,Tt(j))!k1(j)=(1+a2)*Tt(j)^1.5/(Tt(j)+a2);
        !    call Fun_DkbDT(DkbDT(j),Tt(j))
        !    kb(j)    = k*Tt(j);
        !    !kDTt(j)  = Kb(j)*DTt(j)
        !enddo
        call DmukD(mub,DmuD)
        call DmukD(Kb,DkD)
        
        do e = 2,N
            !DJb/DVt
            Jb(4*e-3,e) = 1.0_8
            Jb(4*e-3,e-1) = -1.0_8
            Jb(4*e-2,e) = DUt(e)
            Jb(4*e-1,e) = DWt(e)
            Jb(4*e  ,e) = DTt(e)
            !DJb/DUt
            Jb(4*e-2,e+N+1) = S(i0+i-1)*(t0*Ut(e)+U_b(e)) + S(i0+i-1)*t0*Ut(e) + 2.0_8*beta1*Ut(e)
            Jb(4*e-1,e+N+1) = S(i0+i-1)*(t0*Wt(e)+W_b(e))
            Jb(4*e  ,e+N+1) = S(i0+i-1)*(t0*Tt(e)+T_b(e))! - 2.0_8*(r-1.0_8)*Ma*Ma*beta*Tt(e)
            !DJb/DWt
            Jb(4*e-1,e+2*N+2) = S(i0+i-1)*t0*Ut(e)
            !DJb/DTt
            Jb(4*e-2,e+3*N+3) = -beta1
            Jb(4*e  ,e+3*N+3) = S(i0+i-1)*t0*Ut(e) - (r-1.0_8)*DmubDT(e)*((Maue*DUt(e))**2+(Maue*DWt(e))**2) !- 2.0_8*(r-1.0_8)*Ma*Ma*beta*Ut(e)
        enddo
        
        
        do e = 1,N+1
            
            !DJb/DUt
            
            do j = 2,N
                Jb(4*j-3,e+N+1) = Jb(4*j-3,e+N+1) + INN(j,e)*(S(i0+i-1)*t0+0.5_8*(1.0_8+beta0))
                Jb(4*j-2,e+N+1) = Jb(4*j-2,e+N+1) + Vt(j)*LN(j,e) - DmuD(j,e)
                Jb(4*j  ,e+N+1) = Jb(4*j  ,e+N+1) - 2.0_8*(r-1.0_8)*Maue*Maue*mub(j)*DUt(j)*LN(j,e)
            
            !DJb/DWt
            
            
                Jb(4*j-1,e+2*N+2) = Jb(4*j-1,e+2*N+2) + Vt(j)*LN(j,e) - DmuD(j,e)
                Jb(4*j  ,e+2*N+2) = - 2.0_8*(r-1.0_8)*Mawe*Mawe*mub(j)*DWt(j)*LN(j,e)
            
            !DJb/DTt
            
            
                Jb(4*j-2,e+3*N+3) = Jb(4*j-2,e+3*N+3) - LN(j,e)*DmubDT(e)*DUt(e)
                Jb(4*j-1,e+3*N+3) = Jb(4*j-1,e+3*N+3) - LN(j,e)*DmubDT(e)*DWt(e)
                Jb(4*j  ,e+3*N+3) = Jb(4*j  ,e+3*N+3) + Vt(j)*LN(j,e) - (LN(j,e)*DkbDT(e)*DTt(e))/Pr - DkD(j,e)/Pr
            enddo
        enddo
        !边界条件
        Jb(1,1) = 1.0_8;
        Jb(2,N+2) = 1.0_8;
        Jb(3,2*N+3) = 1.0_8; 
        do e = 1,N+1
            Jb(4,e+3*N+3) = LN(1,e)
        enddo
        
        Jb(4*N+1,N+1) = 1.0_8;Jb(4*N+1,N) = -1.0_8;
        Jb(4*N+2,2*N+2) = 1.0_8; Jb(4*N+3,3*N+3) = 1.0_8; Jb(4*N+4,4*N+4) = 1.0_8;
        do e = 1,N+1
            Jb(4*N+1,e+N+1) = INN(N+1,e)*(S(i0+i-1)*t0+0.5_8*(1.0_8+beta0))
        enddo
    
        deallocate(U1,W1,T1)
end subroutine
    
Subroutine dx2orderC(y1,dx,U_b)
implicit none 
    real*8,intent(in) :: y1(2) !y_i和y_i-1~y_i-5
    real*8,intent(in) :: dx
    real*8,intent(out):: U_b !输出的效果为：t0*Ui + U_b = dU/dx
    
    U_b = ( - 4.0_8*y1(1) + 1.0_8*y1(2) )/2.0_8/dx
    
    end subroutine

Subroutine dx3orderC(y1,dx,U_b)
implicit none 
    real*8,intent(in) :: y1(3) !y_i和y_i-1~y_i-5
    real*8,intent(in) :: dx
    real*8,intent(out):: U_b !输出的效果为：t0*Ui + U_b = dU/dx
    
    U_b = ( - 18.0_8*y1(1) + 9.0_8*y1(2) -2.0_8*y1(3))/6.0_8/dx
    
    end subroutine   
    
Subroutine dx4orderC(y1,dx,U_b)
implicit none 
    real*8,intent(in) :: y1(4) !y_i和y_i-1~y_i-5
    real*8,intent(in) :: dx
    real*8,intent(out):: U_b !输出的效果为：t0*Ui + U_b = dU/dx
    
    U_b = ( - 96.0_8*y1(1) + 72.0_8*y1(2) - 32.0_8*y1(3) + 6.0_8*y1(4) )/24.0_8/dx
    
    end subroutine    
    
Subroutine dx5orderC(y1,dx,U_b)
implicit none 
    real*8,intent(in) :: y1(5) !y_i和y_i-1~y_i-5
    real*8,intent(in) :: dx
    real*8,intent(out):: U_b !输出的效果为：t0*Ui + U_b = dU/dx
    
    U_b = ( - 600.0_8*y1(1) + 600.0_8*y1(2) - 400.0_8*y1(3) + 150.0_8*y1(4) - 24.0_8*y1(5))/120.0_8/dx
    
end subroutine    