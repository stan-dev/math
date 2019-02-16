C     --------------------------------------------------------------------
C     SUNDIALS Copyright Start
C     Copyright (c) 2002-2019, Lawrence Livermore National Security
C     and Southern Methodist University.
C     All rights reserved.
C
C     See the top-level LICENSE and NOTICE files for details.
C
C     SPDX-License-Identifier: BSD-3-Clause
C     SUNDIALS Copyright End
C     --------------------------------------------------------------------
C     FCVODE Example Problem: Robertson kinetics, Lapack linear solver
C                             with dense user Jacobian.
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
C     package. This solution uses the BDF method and a user-supplied
C     Jacobian routine, and prints results at t = .4, 4., ..., 4.e10.
C     It uses ITOL = 2 and ATOL much smaller for y2 than y1 or y3
C     because y2 has much smaller values. At the end of the run,
C     various counters of interest are printed.
C
C     Note that this problem should only work with SUNDIALS configured
C     to use 'realtype' as 'double' and 'sunindextype' as '32bit'
C     --------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER*4 IER, LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNJE, LNGE
      INTEGER*4 METH, ITOL, ITASK, JOUT, NOUT, IERROOT
      INTEGER*4 INFO(2)
      INTEGER*4 I, NEQ
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR, IOUT(25)
      DOUBLE PRECISION RTOL, T, T0, TOUT
      DOUBLE PRECISION Y(3), ATOL(3), ROUT(10), RPAR
C
      DATA LNST/3/, LNFE/4/, LNETF/5/, LNCF/6/, LNNI/7/, LNSETUP/8/,
     1     LNGE/12/, LNJE/17/
C
      NEQ = 3
      T0 = 0.0D0
      Y(1) = 1.0D0
      Y(2) = 0.0D0
      Y(3) = 0.0D0
      METH = 2
      ITOL = 2
      RTOL = 1.0D-4
      ATOL(1) = 1.0D-8
      ATOL(2) = 1.0D-14
      ATOL(3) = 1.0D-6
      TOUT = 0.4D0
      ITASK = 1
      JOUT = 0
      NOUT = 12
C
      WRITE(6,10) NEQ
 10   FORMAT('Dense example problem:'//
     1       ' Robertson kinetics, NEQ = ', I2//)

C     create serial vector
      CALL FNVINITS(1, NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,20) IER
 20     FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
        STOP
      ENDIF

C     initialize dense matrix module
      CALL FSUNDENSEMATINIT(1, NEQ, NEQ, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,25) IER
 25     FORMAT(///' SUNDIALS_ERROR: FSUNDENSEMATINIT returned IER = ',
     1         I5)
        STOP
      ENDIF

C     initialize LAPACK dense linear solver module
      CALL FSUNLAPACKDENSEINIT(1, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,28) IER
 28     FORMAT(///' SUNDIALS_ERROR: FSUNLAPACKDENSEINIT returned IER = '
     1         , I5)
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

C     indicate a dense Jacobian function is provided
      CALL FCVDENSESETJAC(1, IER)
      IF (IER .NE. 0) THEN
        WRITE(6,45) IER
 45     FORMAT(///' SUNDIALS_ERROR: FCVDENSESETJAC returned IER = ', I5)
        CALL FCVFREE
        STOP
      ENDIF

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
C
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
C
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

      SUBROUTINE FCVDJAC(N, T, Y, FY, JAC, H, IPAR, RPAR, 
     1                   V1, V2, V3, IER)
C     Fortran routine for dense user-supplied Jacobian
      IMPLICIT NONE
C
C The following declaration specification should match C type long int
      INTEGER*8 N, IPAR(*)
      INTEGER*4 IER
      DOUBLE PRECISION T, Y(*), FY(*), JAC(N,*), H, RPAR(*)
      DOUBLE PRECISION V1(*), V2(*), V3(*)
C
      DOUBLE PRECISION Y1, Y2, Y3
C
      Y1 = Y(1)
      Y2 = Y(2)
      Y3 = Y(3)
      JAC(1,1) = -0.04D0
      JAC(1,2) = 1.0D4 * Y3
      JAC(1,3) = 1.0D4 * Y2
      JAC(2,1) =  0.04D0
      JAC(2,2) = -1.0D4 * Y3 - 6.0D7 * Y2
      JAC(2,3) = -1.0D4 * Y2
      JAC(3,3) = 0.0D0
      JAC(3,2) = 6.0D7 * Y2
      JAC(3,3) = 0.0D0
C
      IER = 0
C
      RETURN
      END
