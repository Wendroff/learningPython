        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 02 09:30:15 2015
        MODULE DMUKD__genmod
          INTERFACE 
            SUBROUTINE DMUKD(MUB,DMUD)
              USE MESH
              REAL(KIND=8), INTENT(IN) :: MUB(N+1)
              REAL(KIND=8), INTENT(OUT) :: DMUD(N+1,N+1)
            END SUBROUTINE DMUKD
          END INTERFACE 
        END MODULE DMUKD__genmod
