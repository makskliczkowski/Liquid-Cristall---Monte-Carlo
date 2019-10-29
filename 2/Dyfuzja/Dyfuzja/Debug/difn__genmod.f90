        !COMPILER-GENERATED INTERFACE MODULE: Thu Mar 22 14:25:30 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIFN__genmod
          INTERFACE 
            FUNCTION DIFN(L,K,N,A,NX,NY,PX,PY)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: L
              INTEGER(KIND=4) :: K
              INTEGER(KIND=4) :: A(L,L)
              INTEGER(KIND=4) :: NX(0:L)
              INTEGER(KIND=4) :: NY(0:L)
              INTEGER(KIND=4) :: PX(0:L)
              INTEGER(KIND=4) :: PY(0:L)
              REAL(KIND=4) :: DIFN
            END FUNCTION DIFN
          END INTERFACE 
        END MODULE DIFN__genmod
