!初步完成于2015年9月8日 —— 王哲夫
!BLE_FORTRAN Version 4.0
!    
!    
!    
program BLE_FORTRAN
use OUTFLOW
use mesh
implicit none
integer*8 :: i
logical   :: flag
    
    call initial_mesh!(method_flag)
    call initial_Ref
    call initial_outflow
    call initial_init
    call GenCSR3(N)
    call mkl_set_NUM_THREADS(4)
    
    open(11,file='output\V_dimensionless.dat')
    open(12,file='output\U_dimensionless.dat')
    open(13,file='output\W_dimensionless.dat')
    open(14,file='output\T_dimensionless.dat')
    
    i = 1
    call output_demensionless(i)
    do i = 2,iend-i0+1
        call meanorder(i)
        call isFaiNaN(i,flag)
        if (flag) then
            iend = i
            exit
        endif
        call output_demensionless(i)
    enddo
    close(11)
    close(12)
    close(13)
    close(14)
    call TransferXY
    call output_tecout
    
    deallocate(y,LN,INN)
end program