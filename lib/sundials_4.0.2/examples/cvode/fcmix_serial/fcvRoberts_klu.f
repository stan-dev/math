C     ----------------------------------------------------------------
C     Programmer(s): Ting Yan @ SMU
C          Based on cvRoberts_klu.c and modified to Fortran 77
C     ----------------------------------------------------------------
C     SUNDIALS Copyright Start
C     Copyright (c) 2002-2019, Lawrence Livermore National Security
C     and Southern Methodist University.
C     All rights reserved.
C
C     See the top-level LICENSE and NOTICE files for details.
C
C     SPDX-License-Identifier: BSD-3-Clause
C     SUNDIALS Copyright End
C     ----------------------------------------------------------------
C     FCVODE Example Problem: Robertson kinetics
C
C     The following is a simple example problem, with the coding
C     needed for its solution by CVODE. The problem is from chemical
C     kinetics, and consists of the following three rate equations:
C
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C     dy3/dt = 3.e7*y2**2
C
C     on the interval from t = 0.0 to t = 4.e10, with initial
C     conditions:
C
C     y1 = 1.0, y2 = y3 = 0.
C
C     The problem is stiff. While integrating the system, we also
C     enable the root finding feature to find the points at which
C     y1 = 1.e-4 or at which y3 = 0.01. The following coding solves
C     this problem with CVODE, using the Fortran/C interface routine
C     package. This solution uses the BDF method, Newton iteration with
C     the the KLU sparse direct linear solver, and a user-supplied
C     Jacobian routine.
C     It uses a scalar relative tolerance and a vector absolute
C     tolerance. Output is printed in decades from t = .4 to t = 4.e10.
C     Run statistics (optional outputs) are printed at the end.
C     ----------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER*4 IER, LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNJE, LNGE
      INTEGER*4 METH, ITOL, ITASK, JOUT, NOUT, IERROOT
      INTEGER*4 INFO(2)
      INTEGER*4 I
C The following declaration specification should match C type long int
      INTEGER*8 NEQ, IPAR, IOUT(22), IVAL, MXNLI, MXETF, NNZ
      INTEGER*4 SPARSETYPE
      DOUBLE PRECISION RTOL, T, T0, TOUT, H0, NLCONV
      DOUBLE PRECISION Y(3), ATOL(3), ROUT(6), RPAR
C
      DATA LNST/3/, LNFE/4/, LNETF/5/, LNCF/6/, LNNI/7/, LNSETUP/8/,
     1     LNGE/12/, LNJE/17/
C
C     Problem constants
C     number of equations
      NEQ = 3
C     initial time
      T0 = 0.0D0
C     initial y components
      Y(1) = 1.0D0
      Y(2) = 0.0D0
      Y(3) = 0.0D0
C     basic integration method, 2 for BDF
      METH = 2
C     type for absolute tolerance, 2 for array
      ITOL = 2
C     scalar relative tolerance
      RTOL = 1.0D-4
C     vector absolute tolerance components
      ATOL(1) = 1.0D-6
      ATOL(2) = 1.0D-12
      ATOL(3) = 1.0D-4
C     first output time
      TOUT = 0.4D0
C     call CVODE in normal mode
      ITASK = 1
      JOUT = 0
C     number of output times
      NOUT = 12
C
      WRITE(6, 10) NEQ
 10   FORMAT('Klu example problem:'//
     1       ' Robertson kinetics, NEQ = ', I2//)

C     create serial vector
      CALL FNVINITS(1, NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
 20     FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
        STOP
      ENDIF

C     initialize sparse matrix module
C     maximum number of nonzeros in the sparse Jac
      NNZ = NEQ * NEQ
C     CSC
      SPARSETYPE = 0
      CALL FSUNSPARSEMATINIT(1, NEQ, NEQ, NNZ, SPARSETYPE, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,25) IER
 25     FORMAT(///' SUNDIALS_ERROR: FSUNSPARSEMATINIT returned IER = ',
     1         I5)
        STOP
      ENDIF

C     initialize KLU sparse direct linear solver module
      CALL FSUNKLUINIT(1, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,28) IER
 28     FORMAT(///' SUNDIALS_ERROR: FSUNKLUINIT returned IER = ', I5)
        STOP
      ENDIF

C     Call FCVMALLOC to create the solver memory and specify the 
C     Backward Differentiation Formula
      CALL FCVMALLOC(T0, Y, METH, ITOL, RTOL, ATOL,
     1               IOUT, ROUT, IPAR, RPAR, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,30) IER
 30     FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
        STOP
      ENDIF

C     Set the FCVODE input
C     max no. of internal steps before t_out
      IVAL = 1000
      CALL FCVSETIIN('MAX_NSTEPS', IVAL, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,31) IER
 31     FORMAT(///' SUNDIALS_ERROR: FCVSETIIN returned IER = ', I5)
        STOP
      ENDIF

C     max no. of error test failures
      MXETF = 20
      CALL FCVSETIIN('MAX_ERRFAIL', MXETF, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,31) IER
        STOP
      ENDIF
      
C     initial step size
      H0 = 1.0D-4 * RTOL
      CALL FCVSETRIN('INIT_STEP', H0, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,32) IER
 32     FORMAT(///' SUNDIALS_ERROR: FCVSETRIN returned IER = ', I5)
        STOP
      ENDIF
      
C     coefficient in the nonlinear convergence test
      NLCONV = 1.0D-4
      CALL FCVSETRIN('NLCONV_COEF', NLCONV, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,32) IER
        STOP
      ENDIF

C     Call FCVROOTINIT to specify the root function g with 2 components
      CALL FCVROOTINIT(2, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,35) IER
 35     FORMAT(///' SUNDIALS_ERROR: FCVROOTINIT returned IER = ', I5)
        CALL FCVFREE
        STOP
      ENDIF

C     attach the matrix and linear solver modules to CVLs interface
      CALL FCVLSINIT(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,40) IER
 40     FORMAT(///' SUNDIALS_ERROR: FCVLSINIT returned IER = ', I5)
        CALL FCVFREE
        STOP
      ENDIF

C     indicate a sparse Jacobian function is provided
      CALL FCVSPARSESETJAC(IER)
      IF (IER .NE. 0) THEN
        WRITE(6,45) IER
 45     FORMAT(///' SUNDIALS_ERROR: FCVSPARSESETJAC returned IER = ',
     1         I5)
        CALL FCVFREE
        STOP
      ENDIF

C     In loop, call FCVODE, print results, and test for error.
      DO WHILE(JOUT .LT. NOUT)
C
        CALL FCVODE(TOUT, T, Y, ITASK, IER)
C
        WRITE(6,50) T, Y(1), Y(2), Y(3)
 50     FORMAT('At t = ', E12.4, '   y = ', 3E14.6)
C
        IF (IER .LT. 0) THEN
          WRITE(6,60) IER, IOUT(15)
 60       FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /,
     1           '                 Linear Solver returned IER = ', I5)
          CALL FCVROOTFREE
          CALL FCVFREE
          STOP
        ENDIF
C
        IF (IER .EQ. 2) THEN
          CALL FCVROOTINFO(2, INFO, IERROOT)
          IF (IERROOT .LT. 0) THEN
            WRITE(6,65) IERROOT
 65         FORMAT(///' SUNDIALS_ERROR: FCVROOTINFO returned IER = ',
     1             I5)
            CALL FCVROOTFREE
            CALL FCVFREE
            STOP
          ENDIF
          WRITE(6,70) (INFO(I), I = 1, 2)
 70       FORMAT(5X, 'Above is a root, INFO() = ', 2I3)
        ENDIF
C
        IF (IER .EQ. 0) THEN
          TOUT = TOUT * 10.0D0
          JOUT = JOUT + 1
        ENDIF
C
      ENDDO

C     obtain a derivative of the solution
      CALL FCVDKY(T, 1, Y, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,80) IER
 80     FORMAT(///' SUNDIALS_ERROR: FCVDKY returned IER = ', I4)
        CALL FCVROOTFREE
        CALL FCVFREE
        STOP
      ENDIF
      WRITE(6,85) Y(1), Y(2), Y(3)
 85   FORMAT(/'Final value of ydot = ', 3E14.6)

C     output run statistics
      WRITE(6,90) IOUT(LNST), IOUT(LNFE), IOUT(LNJE), IOUT(LNSETUP),
     1            IOUT(LNNI), IOUT(LNCF), IOUT(LNETF), IOUT(LNGE)
 90   FORMAT(//'Final statistics:'//
     1       ' No. steps = ', I4, '   No. f-s = ', I4,
     2       '   No. J-s = ', I4, '   No. LU-s = ', I4/
     3       ' No. nonlinear iterations = ', I4/
     4       ' No. nonlinear convergence failures = ', I4/
     5       ' No. error test failures = ', I4/
     6       ' No. root function evals = ', I4)
C
      CALL FCVROOTFREE
      CALL FCVFREE
C
      STOP
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
C     Fortran routine for right-hand side function
      IMPLICIT NONE
C
C The following declaration specification should match C type long int
      INTEGER*8 IPAR(*)
      INTEGER*4 IER
      DOUBLE PRECISION T, Y(*), YDOT(*), RPAR(*)
C
      YDOT(1) = -0.04D0 * Y(1) + 1.0D4 * Y(2) * Y(3)
      YDOT(3) = 3.0D7 * Y(2) * Y(2)
      YDOT(2) = -YDOT(1) - YDOT(3)
C
      IER = 0
C
      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVROOTFN(T, Y, G, IPAR, RPAR, IER)
C     Fortran routine for root finding
      IMPLICIT NONE
C
      DOUBLE PRECISION T, Y(*), G(*), RPAR(*)
C The following declaration specification should match C type long int
      INTEGER*8 IPAR(*)
      INTEGER*4 IER
C
      G(1) = Y(1) - 1.0D-4
      G(2) = Y(3) - 1.0D-2
C
      IER = 0

      RETURN
      END

C     ----------------------------------------------------------------

      SUBROUTINE FCVSPJAC(T, Y, FY, N, NNZ, JDATA, JRVALS,
     1                    JCPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)
C     Fortran routine for user-supplied CSC format Jacobian
      IMPLICIT NONE
C
C The following declaration specification should match C type long int
      INTEGER*8 N, NNZ, IPAR(*)
      INTEGER*8 JRVALS(NNZ), JCPTRS(N+1)
      INTEGER*4 IER
      DOUBLE PRECISION T, Y(*), FY(*), H, RPAR(*)
      DOUBLE PRECISION JDATA(NNZ)
      DOUBLE PRECISION WK1(*), WK2(*), WK3(*)
C
      DOUBLE PRECISION Y1, Y2, Y3
C
      Y1 = Y(1)
      Y2 = Y(2)
      Y3 = Y(3)
      JCPTRS(1) = 0
      JCPTRS(2) = 3
      JCPTRS(3) = 6
      JCPTRS(4) = 9

      JDATA(1) = -0.04D0
      JRVALS(1) = 0
      JDATA(2) = 0.04D0
      JRVALS(2) = 1
      JDATA(3) = 0.0D0
      JRVALS(3) = 2

      JDATA(4) = 1.0D4 * Y3
      JRVALS(4) = 0
      JDATA(5) = -1.0D4 * Y3 - 6.0D7 * Y2
      JRVALS(5) = 1
      JDATA(6) = 6.0D7 * Y2
      JRVALS(6) = 2
      
      JDATA(7) = 1.0D4 * Y2
      JRVALS(7) = 0
      JDATA(8) = -1.0D4 * Y2
      JRVALS(8) = 1
      JDATA(9) = 0.0D0
      JRVALS(9) = 2
C
      IER = 0
C
      RETURN
      END
