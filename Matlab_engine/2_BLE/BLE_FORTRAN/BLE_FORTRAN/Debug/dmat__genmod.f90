        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 02 09:30:13 2015
        MODULE DMAT__genmod
          INTERFACE 
            SUBROUTINE DMAT(Z,DZY,DYZ)
              USE MESH
              REAL(KIND=16), INTENT(IN) :: Z(N+1)
              REAL(KIND=16), INTENT(IN) :: DZY(N+1)
              REAL(KIND=16), INTENT(IN) :: DYZ(N+1)
            END SUBROUTINE DMAT
          END INTERFACE 
        END MODULE DMAT__genmod
