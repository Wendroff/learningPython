Subroutine initial_Ref
    use Referencevalue
    implicit none
    
    open(9,file='input/Ref.dat')
    read(9,*)
    read(9,*) Ureff,Treff,DENreff,MUreff,Cordlenthg
    write(*,*) 'Ureff  =',Ureff
    write(*,*) 'Treff  =',Treff
    write(*,*) 'DENreff=',DENreff
    write(*,*) 'MUreff =',MUreff
    write(*,*) 'C      =',Cordlenthg
    Pr = 0.72d0;r = 1.4d0;Rg=287d0;
    close(9)
    Ma = Ureff/sqrt(r*Rg*Treff)
    write(*,*) 'Ma     =',Ma
    deta0 = sqrt(MUreff*Cordlenthg/DENreff/Ureff)
    write(*,*) 'deta0  =',deta0
    
    
end subroutine