        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 02 09:30:06 2015
        MODULE MY_RMM__genmod
          INTERFACE 
            SUBROUTINE MY_RMM(A,B,C,L,M,N)
              INTEGER(KIND=8), INTENT(IN) :: N
              INTEGER(KIND=8), INTENT(IN) :: M
              INTEGER(KIND=8), INTENT(IN) :: L
              REAL(KIND=8), INTENT(IN) :: A(L,M)
              REAL(KIND=8), INTENT(IN) :: B(M,N)
              REAL(KIND=8), INTENT(OUT) :: C(L,N)
            END SUBROUTINE MY_RMM
          END INTERFACE 
        END MODULE MY_RMM__genmod
