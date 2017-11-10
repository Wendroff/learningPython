Subroutine My_pardiso(Jb,b,N)
use M_CSR3
    implicit none
    integer  :: N
    REAL*8   :: Jb(N,N)
    INTEGER*8 pt(64)

!C..     All other variables

    INTEGER     maxfct, mnum, mtype, phase,  nrhs, error, msglvl
    INTEGER     iparm(64)
    real*8  ::   b(n)
    real*8  ::   x(n)
    integer           ::       ii,j,k

    INTEGER i, idum(1)
    real*8  ddum(1)
    
    
    
    do i = 1,num
        a(i) = Jb(IndexA(i),JA(i))
    enddo
    maxfct = 1
    mnum   = 1
    nrhs   = 1
    
    !C
!C  .. Setup Pardiso control parameters und initialize the solvers
!C     internal adress pointers. This is only necessary for the FIRST
!C     call of the PARDISO solver.
!C
      do i = 1, 64
         iparm(i) = 0
      end do
    
      iparm(1) = 1 ! no solver default
      iparm(2) = 2 ! fill-in reordering from METIS
      iparm(3) = 4 ! numbers of processors, value of OMP_NUM_THREADS
      iparm(4) = 0 ! no iterative-direct algorithm
      iparm(5) = 0 ! no user fill-in reducing permutation
      iparm(6) = 1 ! =0 solution on the first n compoments of x
      iparm(7) = 0 ! not in use
      iparm(8) = 2 ! numbers of iterative refinement steps
      iparm(9) = 0 ! not in use
      iparm(10) = 13 ! perturbe the pivot elements with 1E-13
      iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
      iparm(12) = 0 ! not in use
      iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric).
      iparm(14) = 0 ! Output: number of perturbed pivots
      iparm(15) = 0 ! not in use
      iparm(16) = 0 ! not in use
      iparm(17) = 0 ! not in use
      iparm(18) = -1 ! Output: number of nonzeros in the factor LU
      iparm(19) = -1 ! Output: Mflops for LU factorization
      iparm(20) = 0 ! Output: Numbers of CG Iterations
      error     = 0 ! initialize error flag
      msglvl    = 0 ! print statistical information
      mtype     = 11! real unsymmetric matrix

!C.. Initiliaze the internal solver memory pointer. This is only
!C necessary for the FIRST call of the PARDISO solver.
      do i = 1, 64
         pt(i) = 0
      end do

!C..   Reordering and Symbolic Factorization, This step also allocates
!C     all memory that is necessary for the factorization
!C
      phase     = 11      ! only reordering and symbolic factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,idum, nrhs, iparm, msglvl, ddum, ddum, error)

      !WRITE(*,*) 'Reordering completed ... '

      IF (error .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', error
        STOP 1
      END IF

      !WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
      !WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

!C.. Factorization.
      phase     = 22  ! only factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,idum, nrhs, iparm, msglvl, ddum, ddum, error)

      !WRITE(*,*) 'Factorization completed ...  '
      IF (error .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', error
        STOP 1
      ENDIF

!C.. Back substitution and iterative refinement
      phase     = 33  ! only factorization
      iparm(8)  = 3   ! max numbers of iterative refinement steps
      
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,idum, nrhs, iparm, msglvl, b, x, error)
      !WRITE(*,*) 'Solve completed ... '

!C.. Termination and release of memory
      phase     = -1           ! release internal memory
      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,idum, nrhs, iparm, msglvl, ddum, ddum, error)
      
      return
end subroutine