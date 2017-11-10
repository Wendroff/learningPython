subroutine Transferxy
use physicalvalue
use mesh
use Phy
use OUTFLOW
implicit none
real*8 :: mue(i0:iend),Re1(i0:iend),ReL(i0:iend),L1(i0:iend),INTT(N+1),LNyout(N+1)
real*8 :: Vt(iend-i0+1,N+1),Ut(N+1),Wt(N+1),Tt(N+1),Tx(7:iend-i0+1,N+1),Tetax(7:iend-i0+1,N+1),etax(7:iend-i0+1,N+1)
real*8 :: INTx(N+1)
real*8 :: dx,dRex1,youtx
integer*8 :: i,j
    allocate(U1(iend-i0+1,N+1),V1(iend-i0+1,N+1),W1(iend-i0+1,N+1),DEN1(iend-i0+1,N+1),T1(iend-i0+1,N+1),Xout(iend-i0+1),Yout(iend-i0+1,N+1))
    mue(i0) = 0d0
    Re1(i0) = 0d0
    ReL(i0) = 0d0
    L1(i0)  = 0d0
    do i=1,iend-i0+1
        do j=1,N+1
            Vt(i,j) = Fai(        j,i)
            Ut(j)   = Fai(   N+1 +j,i)
            Wt(j)   = Fai(2*(N+1)+j,i)
            Tt(j)   = Fai(3*(N+1)+j,i)
        enddo
        mue(i0+i-1)=1.789e-5_8*(Te(i0+i-1)/288.0_8)**1.5_8*(288.0_8+110.4_8)/(Te(i0+i-1)+110.4_8)
        Re1(i0+i-1)=DENe(i0+i-1)*Ue(i0+i-1)/mue(i0+i-1)
        ReL(i0+i-1)=sqrt(DENe(i0+i-1)*Ue(i0+i-1)*S(i0+i-1)/mue(i0+i-1))
        L1(i0+i-1) =sqrt(mue(i0+i-1)*S(i0+i-1)/DENe(i0+i-1)/Ue(i0+i-1))
        !INT=IN*T1(:,i);
        call my_mv(INN,Tt,INTT,N+1)
        do j=2,N+1
            INTT(j) =INTT(j)+INTT(j-1)
        enddo
        
        do j = 1,N+1
            yout(i,j)=L1(i0+i-1)*INTT(j)
            U1(i,j)  =Ut(j)
            T1(i,j)  =Tt(j)
            W1(i,j)  =Wt(j)
            !V1(i,j)  =0d0
        enddo
        xout(i) = S(i0+i-1)
    enddo
    deallocate(Fai)
    
    do i=7,iend-i0+1
        dx = xout(i) - xout(i-1)
        do j = 1,N+1
            Tx(i,j)=(274.0_8*T1(i,j)-600.0_8*T1(i-1,j)+600.0_8*T1(i-2,j)-400.0_8*T1(i-3,j)+150.0_8*T1(i-4,j)-24.0_8*T1(i-5,j))/dx/120.0_8
        enddo
        !INTx=IN*Tx(i,j);
        call my_mv(INN,Tx(i,:),INTx,N+1)
        do j=2,N+1
            INTx(j)=INTx(j)+INTx(j-1)
        enddo
      
        dRex1=(274.0_8*Re1(i0+i-1)-600.0_8*Re1(i0+i-2)+600.0_8*Re1(i0+i-3)-400.0_8*Re1(i0+i-4)+150.0_8*Re1(i0+i-5)-24.0_8*Re1(i0+i-6))/dx/120.0_8
        do j = 1,N+1
            Tetax(i,j)=-yout(i,j)/2.0_8/L1(i+i0-1)*(1.0_8/S(i+i0-1)-1.0_8/Re1(i+i0-1)*dRex1)-INTx(j)
            youtx = (274.0_8*yout(i,j)-600.0_8*yout(i-1,j)+600.0_8*yout(i-2,j)-400.0_8*yout(i-3,j)+150.0_8*yout(i-4,j)-24.0_8*yout(i-5,j))/dx/120.0_8
            call my_mv(LN,yout(i,:),LNyout,N+1)
            etax(i,j) = -youtx/LNyout(j)

            V1(i,j)=T1(i,j)*Vt(i,j)/ReL(i+i0-1)-1.0_8/ReL(i+i0-1)*S(i+i0-1)*U1(i,j)*T1(i,j)*etax(i,j)
        enddo
    enddo
    
    do i=7,iend-i0+1
        do j = 1,N+1
            Den1(i,j)=1.0_8/T1(i,j)*Dene(i+i0-1) !不能和T1有量纲化交换顺序
            U1(i,j)  =U1(i,j)*Ue(i+i0-1)
            T1(i,j)  =T1(i,j)*Te(i+i0-1)
            W1(i,j)  =W1(i,j)*We(i+i0-1)
            V1(i,j)  =V1(i,j)*Ue(i+i0-1)
        enddo
    enddo
end subroutine