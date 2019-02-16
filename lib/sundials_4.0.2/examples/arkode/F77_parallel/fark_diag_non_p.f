C     ----------------------------------------------------------------
C     Programmer(s): Daniel R. Reynolds @ SMU
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
C     Diagonal ODE example. Nonstiff case: alpha = 10/NEQ.
C     ----------------------------------------------------------------
C
C     Include MPI-Fortran header file for MPI_COMM_WORLD, MPI types.
C
      IMPLICIT NONE
C
      INCLUDE "mpif.h"
      INTEGER*8 NLOCAL
      PARAMETER (NLOCAL=2)
C
      INTEGER*4 IER, MYPE, NPES, NOUT, LNST, LNST_ATT, LNFE, LNFI, LNNI
      INTEGER*4 LNCF, LNETF, METH, IATOL, ITASK, JOUT
      INTEGER*8 NEQ, I, NST, NST_ATT, NFE, NFI, NNI, NCFN, NETF
      INTEGER*8 IOUT(35), IPAR(2)
      DOUBLE PRECISION Y(128), ROUT(6), RPAR(1)
      DOUBLE PRECISION ATOL, RTOL, DTOUT, T, ALPHA, TOUT
      DOUBLE PRECISION ERMAX, ERRI, GERMAX
C
      DATA ATOL/1.0D-10/, RTOL/1.0D-5/, DTOUT/0.1D0/, NOUT/10/
      DATA LNST/3/, LNST_ATT/6/, LNFE/7/, LNFI/8/, LNNI/11/, LNCF/12/, 
     1     LNETF/10/
C
C     Get NPES and MYPE.  Requires initialization of MPI.
      CALL MPI_INIT(IER)
      IF (IER .NE. 0) THEN
         WRITE(6,5) IER
 5       FORMAT(///' MPI_ERROR: MPI_INIT returned IER = ', I5)
         STOP
      ENDIF
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPES, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,6) IER
 6       FORMAT(///' MPI_ERROR: MPI_COMM_SIZE returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYPE, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,7) IER
 7       FORMAT(///' MPI_ERROR: MPI_COMM_RANK returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
C     Set input arguments.
      NEQ = NPES * NLOCAL
      T = 0.0D0
      METH = 1
      IATOL = 1
      ITASK = 1
c     Set parameter ALPHA
      ALPHA  = 10.0D0 / NEQ
C
C     Load IPAR and RPAR
      IPAR(1) = NLOCAL
      IPAR(2) = MYPE
      RPAR(1) = ALPHA
C
      DO 10 I = 1, NLOCAL
  10    Y(I) = 1.0D0
C
      IF (MYPE .EQ. 0) THEN
         WRITE(6,11) NEQ, NLOCAL, ALPHA
  11     FORMAT('Diagonal test problem:'//' NEQ = ', I3, /
     1        ' NLOCAL = ', I3, /
     2        ' parameter alpha = ', F8.3)
         WRITE(6,12)
  12     FORMAT(' ydot_i = -alpha*i * y_i (i = 1,...,NEQ)')
         WRITE(6,13) RTOL, ATOL
  13     FORMAT(' RTOL, ATOL = ', 2E10.1)
         WRITE(6,14)
  14     FORMAT(' Method is ERK')
         WRITE(6,15) NPES
  15     FORMAT(' Number of processors = ', I3//)
      ENDIF
C
      CALL FNVINITP(MPI_COMM_WORLD, 4, NLOCAL, NEQ, IER)
C
      IF (IER .NE. 0) THEN
         WRITE(6,20) IER
  20     FORMAT(///' SUNDIALS_ERROR: FNVINITP returned IER = ', I5)
         CALL MPI_FINALIZE(IER)
         STOP
      ENDIF
C
      CALL FARKMALLOC(T, Y, METH, IATOL, RTOL, ATOL,
     1                IOUT, ROUT, IPAR, RPAR, IER)
C
      IF (IER .NE. 0) THEN
         WRITE(6,30) IER
  30     FORMAT(///' SUNDIALS_ERROR: FARKMALLOC returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
C     Loop through tout values, call solver, print output, test for failure.
      TOUT = DTOUT
      DO 70 JOUT = 1, NOUT
C
        CALL FARKODE(TOUT, T, Y, ITASK, IER)
C
        IF (MYPE .EQ. 0) WRITE(6,40) T, IOUT(LNST), IOUT(LNST_ATT), 
     &       IOUT(LNFE), IOUT(LNFI)
  40    FORMAT(' t = ', D10.2, 5X, 'steps = ', I5,
     &         '  (attempted = ', I5, '),  fe = ', I5, 
     &         '  fi = ', I5)
C
        IF (IER .NE. 0) THEN
           WRITE(6,60) IER, IOUT(16)
  60       FORMAT(///' SUNDIALS_ERROR: FARKODE returned IER = ', I5, /,
     &            '                 Linear Solver returned IER = ', I5)
           CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
           STOP
        ENDIF
C
        TOUT = TOUT + DTOUT
  70    CONTINUE
C
C     Get max. absolute error in the local vector.
      ERMAX = 0.0D0
      DO 75 I = 1, NLOCAL
        ERRI  = Y(I) - EXP(-ALPHA * (MYPE * NLOCAL + I) * T)
        ERMAX = MAX(ERMAX, ABS(ERRI))
  75    CONTINUE
C     Get global max. error from MPI_REDUCE call.
      CALL MPI_REDUCE(ERMAX, GERMAX, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
     1                0, MPI_COMM_WORLD, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,80) IER
  80     FORMAT(///' MPI_ERROR: MPI_REDUCE returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
      IF (MYPE .EQ. 0) WRITE(6,85) GERMAX
  85  FORMAT(/'Max. absolute error is ', E10.2/)
C
C     Print final statistics.
      NST = IOUT(LNST)
      NST_ATT = IOUT(LNST_ATT)
      NFE = IOUT(LNFE)
      NFI = IOUT(LNFI)
      NNI = IOUT(LNNI)
      NCFN = IOUT(LNCF)
      NETF = IOUT(LNETF)
      IF (MYPE .EQ. 0) WRITE (6,90) NST, NST_ATT, NFE, NFI, NNI, NCFN, 
     &     NETF
  90  FORMAT(/'Final statistics:'//
     &       ' number of steps = ', I5/
     &       ' number of step attempts = ', I5/
     &       ' number of fe evals. = ', I5/
     &       ' number of fi evals. = ', I5/
     &       ' number of nonlinear iters. = ', I5/
     &       ' number of nonlinear conv. failures = ', I3/
     &       ' number of error test failures = ', I3)
C
C     Free the memory and finalize MPI.
      CALL FARKFREE
      CALL MPI_FINALIZE(IER)
      IF (IER .NE. 0) THEN
         WRITE(6,95) IER
 95      FORMAT(///' MPI_ERROR: MPI_FINALIZE returned IER = ', I5)
         STOP
      ENDIF
C
      STOP
      END
C
C     ------------------------------------------------------------------------
C
      SUBROUTINE FARKIFUN(T, Y, YDOT, IPAR, RPAR, IER)
C     Routine for right-hand side function fi
C
      IMPLICIT NONE
C
      INTEGER*4 IER
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION T, Y(*), YDOT(*), RPAR(*)
C
      INTEGER*8 NLOCAL, I
C
      NLOCAL = IPAR(1)
C
      DO I = 1, NLOCAL
         YDOT(I) = 0.D0
      ENDDO
C
      IER = 0
C
      RETURN
      END
C     ------------------------------------------------------------------------
C
      SUBROUTINE FARKEFUN(T, Y, YDOT, IPAR, RPAR, IER)
C     Routine for right-hand side function fe
C
      IMPLICIT NONE
C
      INTEGER*4 IER
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION T, Y(*), YDOT(*), RPAR(*)
C
      INTEGER*8 MYPE
      INTEGER*8 NLOCAL, I
      DOUBLE PRECISION ALPHA
C
      NLOCAL = IPAR(1)
      MYPE = IPAR(2)
      ALPHA = RPAR(1)
C
      DO I = 1, NLOCAL
         YDOT(I) = -ALPHA * (MYPE * NLOCAL + I) * Y(I)
      ENDDO
C
      IER = 0
C
      RETURN
      END
