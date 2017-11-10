Subroutine Finte1order(i,dx,mue,dmuex,dDenex,duex)
    use OUTFLOW
    use Referencevalue
    use mesh
    implicit none
    integer*8,intent(in) :: i !对应于Fai的指标
    real*8,intent(in)    :: dx
    real*8,intent(out)   :: mue,dmuex,dDenex,duex
    real*8               :: mu1,mu

    mu=1.789e-5_8*(Te(i0+i-1)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-1)+110.4_8)
    mue=mu
    mu1=1.789e-5_8*(Te(i0+i-2)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-2)+110.4_8)
    dmuex=(mu-mu1)/dx
     
    dDenex=(DENe(i0+i-1)-DENe(i0+i-2))/dx
    
    duex=(Ue(i0+i-1)-Ue(i0+i-2))/dx
    !call sutherland(mue,Te(i0+i-1))
    !call sutherland(mu1,Te(i0+i-2))
    !dmuex=(mue-mu1)/dx;
    ! 
    !dDenex=(DENe(i0+i-1)-DENe(i0+i-2))/dx;
    !
    !duex=(Ue(i0+i-1)-Ue(i0+i-2))/dx;
    
    !Pe  = DENe(i0+i-1)*Te(i0+i-1)
    !Pe1 = DENe(i0+i-2)*Te(i0+i-2)
    !dPex= (Pe-Pe1)/dx
end subroutine
    
!subroutine sutherlandK(k,T)
!use Referencevalue    
!implicit none
!    real*8,intent(in)  :: T    !注意这里是用Treff=288.15K，MUreff=1.78E-5无量纲化后的
!    real*8,intent(out) :: k
!    real*8             :: C
!    !C = 194.4/Treff
!    !k = T**1.5*(C+1.0_8)/(C+T)
!    
!end subroutine
!    
!subroutine sutherland(mu,T)
!use Referencevalue    
!implicit none
!    real*8,intent(in)  :: T    !注意这里是用Treff=288.15K，MUreff=1.78E-5无量纲化后的
!    real*8,intent(out) :: mu
!    real*8             :: C
!    !C = 110.4/Treff
!    !mu = T**1.5*(C+1.0_8)/(C+T)
!    mu = 1.789e-5_8*(T/288.0_8)**1.5_8*(288.0_8+110.4_8)/(T+110.4_8)
!end subroutine    
!    
!subroutine Fun_DmubDT(DmubDT,T)
!use Referencevalue 
!implicit none
!    real*8,intent(in) :: T
!    real*8,intent(out):: DmubDT
!    real*8 :: C
!    C = 110.4_8/Treff
!    DmubDT = (C+1.0_8)*(C-T)/2.0_8/(C+T)/(C+T)/(T**0.5_8)
!    
!end subroutine
!
!subroutine Fun_DkbDT(DkbDT,T)
!use Referencevalue 
!implicit none
!    real*8,intent(in) :: T
!    real*8,intent(out):: DkbDT
!    real*8 :: C
!    C = 194.4_8/Treff
!    DkbDT = T**1.5_8*(C+1.0_8)*(5.0_8*C+3.0_8*T)/2.0_8/(C+T)**2
!    
!end subroutine
    
Subroutine Finte2order(i,dx,mue,dmuex,dDenex,duex)
use OUTFLOW
use Referencevalue
use mesh
implicit none
integer*8,intent(in) :: i !对应于Fai的指标
real*8,intent(in)    :: dx
real*8,intent(out)   :: mue,dmuex,dDenex,duex
real*8               :: mu1,mu2,mu
integer*4            :: j
    
    mu=1.789e-5_8*(Te(i0+i-1)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-1)+110.4_8);
    mue=mu;
    mu1=1.789e-5_8*(Te(i0+i-2)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-2)+110.4_8);
    mu2=1.789e-5_8*(Te(i0+i-3)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-3)+110.4_8);
    dmuex=(1.5_8*mu-2.0_8*mu1+0.5_8*mu2)/dx;
     
    dDenex=(1.5_8*DENe(i0+i-1)-2.0_8*DENe(i0+i-2)+0.5_8*DENe(i0+i-3))/dx;
     
    duex=(1.5_8*Ue(i0+i-1)-2.0_8*Ue(i0+i-2)+0.5_8*Ue(i0+i-3))/dx;
    
    !call sutherland(mue,Te(i0+i-1))
    !do j = 1,2
    !    call sutherland(mu1(j),Te(i0+i-1-j))
    !    Dene1(j)= DENe(i0+i-1-j)
    !    ue1(j)  = Ue(i0+i-1-j)
    !    Pe1(j)  = DENe(i0+i-1-j)*Te(i0+i-1-j)
    !enddo
    !call dx2order(mue,mu1,dx,dmuex)
    !call dx2order(Ue(i0+i-1),ue1,dx,duex)
    !call dx2order(DENe(i0+i-1),DENe1,dx,dDenex)
    !Pe  = DENe(i0+i-1)*Te(i0+i-1)
    !call dx2order(Pe,Pe1,dx,dPex)
    !dmuex=(mue-mu1)/dx;
    
    !dDenex=(DENe(i0+i-1)-DENe(i0+i-2))/dx;
    
    !duex=(Ue(i0+i-1)-Ue(i0+i-2))/dx;
    
    
    !Pe1 = DENe(i0+i-2)*Te(i0+i-2)
    !dPex= (Pe-Pe1)/dx
end subroutine    
    
    
!Subroutine dx2order(y0,y1,dx,dydx)
!implicit none 
!    real*8,intent(in) :: y0,y1(2) !y_i和y_i-1~y_i-5
!    real*8,intent(in) :: dx
!    real*8,intent(out):: dydx
!    
!    dydx = ( 3.0_8*y0 - 4.0_8*y1(1) + 1.0_8*y1(2) )/2.0_8/dx
!end subroutine    
    
Subroutine Finte3order(i,dx,mue,dmuex,dDenex,duex)
    use OUTFLOW
    use Referencevalue
    use mesh
    implicit none
    integer*8,intent(in) :: i !对应于Fai的指标
    real*8,intent(in)    :: dx
    real*8,intent(out)   :: mue,dmuex,dDenex,duex
    real*8               :: mu1,mu2,mu3,mu
    integer*4            :: j

    mu=1.789e-5_8*(Te(i0+i-1)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-1)+110.4_8)
    mue=mu
    mu1=1.789e-5_8*(Te(i0+i-2)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-2)+110.4_8)
    mu2=1.789e-5_8*(Te(i0+i-3)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-3)+110.4_8)
    mu3=1.789e-5_8*(Te(i0+i-4)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-4)+110.4_8)
    dmuex=(11.0_8*mu-18.0_8*mu1+9.0_8*mu2-2.0_8*mu3)/dx/6.0_8
     
    dDenex=(11.0_8*DENe(i0+i-1)-18.0_8*DENe(i0+i-2)+9.0_8*DENe(i0+i-3)-2.0_8*DENe(i0+i-4))/dx/6.0_8
     
    duex=(11.0_8*Ue(i0+i-1)-18.0_8*Ue(i0+i-2)+9.0_8*Ue(i0+i-3)-2.0_8*Ue(i0+i-4))/dx/6.0_8
    
    !call sutherland(mue,Te(i0+i-1))
    !do j = 1,3
    !    call sutherland(mu1(j),Te(i0+i-1-j))
    !    Dene1(j)= DENe(i0+i-1-j)
    !    ue1(j)  = Ue(i0+i-1-j)
    !    Pe1(j)  = DENe(i0+i-1-j)*Te(i0+i-1-j)
    !enddo
    !call dx3order(mue,mu1,dx,dmuex)
    !call dx3order(Ue(i0+i-1),ue1,dx,duex)
    !call dx3order(DENe(i0+i-1),DENe1,dx,dDenex)
    !Pe  = DENe(i0+i-1)*Te(i0+i-1)
    !call dx3order(Pe,Pe1,dx,dPex)
    !dmuex=(mue-mu1)/dx;
    
    !dDenex=(DENe(i0+i-1)-DENe(i0+i-2))/dx;
    
    !duex=(Ue(i0+i-1)-Ue(i0+i-2))/dx;
    
    
    !Pe1 = DENe(i0+i-2)*Te(i0+i-2)
    !dPex= (Pe-Pe1)/dx
end subroutine    
    
    
    
!Subroutine dx3order(y0,y1,dx,dydx)
!implicit none 
!    real*8,intent(in) :: y0,y1(3) !y_i和y_i-1~y_i-5
!    real*8,intent(in) :: dx
!    real*8,intent(out):: dydx
!    
!    dydx = ( 11.0_8*y0 - 18.0_8*y1(1) + 9.0_8*y1(2) -2.0_8*y1(3))/6.0_8/dx
!end subroutine
    
Subroutine Finte4order(i,dx,mue,dmuex,dDenex,duex)
    use OUTFLOW
    use Referencevalue
    use mesh
    implicit none
    integer*8,intent(in) :: i !对应于Fai的指标
    real*8,intent(in)    :: dx
    real*8,intent(out)   :: mue,dmuex,dDenex,duex
    real*8               :: mu1,mu2,mu3,mu4,mu
    integer*4            :: j

    mu=1.789e-5_8*(Te(i0+i-1)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-1)+110.4_8)
    mue=mu;
    mu1=1.789e-5_8*(Te(i0+i-2)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-2)+110.4_8)
    mu2=1.789e-5_8*(Te(i0+i-3)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-3)+110.4_8)
    mu3=1.789e-5_8*(Te(i0+i-4)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-4)+110.4_8)
    mu4=1.789e-5_8*(Te(i0+i-5)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-5)+110.4_8)
    dmuex=(50.0_8*mu-96.0_8*mu1+72.0_8*mu2-32.0_8*mu3+6.0_8*mu4)/dx/24.0_8
     
    dDenex=(50.0_8*DENe(i0+i-1)-96.0_8*DENe(i0+i-2)+72.0_8*DENe(i0+i-3)-32.0_8*DENe(i0+i-4)+6.0_8*DENe(i0+i-5) )/dx/24.0_8
     
    duex=(50.0_8*Ue(i0+i-1)-96.0_8*Ue(i0+i-2)+72.0_8*Ue(i0+i-3)-32.0_8*Ue(i0+i-4)+6.0_8*Ue(i0+i-5))/dx/24.0_8
    
    !call sutherland(mue,Te(i0+i-1))
    !do j = 1,4
    !    call sutherland(mu1(j),Te(i0+i-1-j))
    !    Dene1(j)= DENe(i0+i-1-j)
    !    ue1(j)  = Ue(i0+i-1-j)
    !    Pe1(j)  = DENe(i0+i-1-j)*Te(i0+i-1-j)
    !enddo
    !call dx4order(mue,mu1,dx,dmuex)
    !call dx4order(Ue(i0+i-1),ue1,dx,duex)
    !call dx4order(DENe(i0+i-1),DENe1,dx,dDenex)
    !Pe  = DENe(i0+i-1)*Te(i0+i-1)
    !call dx4order(Pe,Pe1,dx,dPex)
    !dmuex=(mue-mu1)/dx;
    
    !dDenex=(DENe(i0+i-1)-DENe(i0+i-2))/dx;
    
    !duex=(Ue(i0+i-1)-Ue(i0+i-2))/dx;
    
    
    !Pe1 = DENe(i0+i-2)*Te(i0+i-2)
    !dPex= (Pe-Pe1)/dx
end subroutine    
    
    
    
!Subroutine dx4order(y0,y1,dx,dydx)
!implicit none 
!    real*8,intent(in) :: y0,y1(4) !y_i和y_i-1~y_i-5
!    real*8,intent(in) :: dx
!    real*8,intent(out):: dydx
!    
!    dydx = ( 50.0_8*y0 - 96.0_8*y1(1) + 72.0_8*y1(2) - 32.0_8*y1(3) + 6.0_8*y1(4) )/24.0_8/dx
!end subroutine

    
Subroutine Finte5order(i,dx,mue,dmuex,dDenex,duex)
    use OUTFLOW
    use Referencevalue
    use mesh
    implicit none
    integer*8,intent(in) :: i !对应于Fai的指标
    real*8,intent(in)    :: dx
    real*8,intent(out)   :: mue,dmuex,dDenex,duex!,Pe,dPex
    real*8               :: mu1,mu2,mu3,mu4,mu5,mu
    integer*4            :: j

    mu=1.789e-5_8*(Te(i0+i-1)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-1)+110.4_8)
    mue=mu;
    mu1=1.789e-5_8*(Te(i0+i-2)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-2)+110.4_8)
    mu2=1.789e-5_8*(Te(i0+i-3)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-3)+110.4_8)
    mu3=1.789e-5_8*(Te(i0+i-4)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-4)+110.4_8)
    mu4=1.789e-5_8*(Te(i0+i-5)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-5)+110.4_8)
    mu5=1.789e-5_8*(Te(i0+i-6)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-6)+110.4_8)
    dmuex=(274.0_8*mu-600.0_8*mu1+600.0_8*mu2-400.0_8*mu3+150.0_8*mu4-24.0_8*mu5)/dx/120.0_8
     
    dDenex=(274.0_8*DENe(i0+i-1)-600.0_8*DENe(i0+i-2)+600.0_8*DENe(i0+i-3)-400.0_8*DENe(i0+i-4)+150.0_8*DENe(i0+i-5)-24.0_8*DENe(i0+i-6))/dx/120.0_8
     
    duex=(274.0_8*Ue(i0+i-1)-600.0_8*Ue(i0+i-2)+600.0_8*Ue(i0+i-3)-400.0_8*Ue(i0+i-4)+150.0_8*Ue(i0+i-5)-24.0_8*Ue(i0+i-6) )/dx/120
    
    !call sutherland(mue,Te(i0+i-1))
    !do j = 1,5
    !    call sutherland(mu1(j),Te(i0+i-1-j))
    !    Dene1(j)= DENe(i0+i-1-j)
    !    ue1(j)  = Ue(i0+i-1-j)
    !    !Pe1(j)  = DENe(i0+i-1-j)*Te(i0+i-1-j)
    !enddo
    !call dx5order(mue,mu1,dx,dmuex)
    !call dx5order(Ue(i0+i-1),ue1,dx,duex)
    !call dx5order(DENe(i0+i-1),DENe1,dx,dDenex)
    !Pe  = DENe(i0+i-1)*Te(i0+i-1)
    !call dx5order(Pe,Pe1,dx,dPex)
    !dmuex=(mue-mu1)/dx;
    
    !dDenex=(DENe(i0+i-1)-DENe(i0+i-2))/dx;
    
    !duex=(Ue(i0+i-1)-Ue(i0+i-2))/dx;
    
    
    !Pe1 = DENe(i0+i-2)*Te(i0+i-2)
    !dPex= (Pe-Pe1)/dx
end subroutine    