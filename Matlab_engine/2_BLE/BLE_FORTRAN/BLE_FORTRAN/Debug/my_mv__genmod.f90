        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 02 09:30:06 2015
        MODULE MY_MV__genmod
          INTERFACE 
            SUBROUTINE MY_MV(A,B,C,N)
              INTEGER(KIND=8), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: A(N,N)
              REAL(KIND=8), INTENT(IN) :: B(N)
              REAL(KIND=8), INTENT(OUT) :: C(N)
            END SUBROUTINE MY_MV
          END INTERFACE 
        END MODULE MY_MV__genmod
