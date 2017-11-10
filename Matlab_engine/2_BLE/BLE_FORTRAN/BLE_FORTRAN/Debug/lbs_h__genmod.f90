        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep 07 10:46:44 2015
        MODULE LBS_H__genmod
          INTERFACE 
            SUBROUTINE LBS_H(I,PHII,LBSU,F,JB,BETA0,BETA1,A1,A2,MAUE,   &
     &MAWE,DX)
              USE MESH
              INTEGER(KIND=8), INTENT(IN) :: I
              REAL(KIND=8), INTENT(IN) :: PHII(4*N+4)
              REAL(KIND=8), INTENT(OUT) :: LBSU(4*N+4)
              REAL(KIND=8), INTENT(OUT) :: F(4*N+4)
              REAL(KIND=8), INTENT(OUT) :: JB(4*N+4,4*N+4)
              REAL(KIND=8), INTENT(IN) :: BETA0
              REAL(KIND=8), INTENT(IN) :: BETA1
              REAL(KIND=8), INTENT(IN) :: A1
              REAL(KIND=8), INTENT(IN) :: A2
              REAL(KIND=8), INTENT(IN) :: MAUE
              REAL(KIND=8), INTENT(IN) :: MAWE
              REAL(KIND=8), INTENT(IN) :: DX
            END SUBROUTINE LBS_H
          END INTERFACE 
        END MODULE LBS_H__genmod
