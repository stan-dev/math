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
C     Diagonal ODE example. Stiff case, with BDF/SPGMR, diagonal
C     preconditioner. Solved with preconditioning on left, then with
C     preconditioning on right.
C     --------------------------------------------------------------------
C
C     Include MPI-Fortran header file for MPI_COMM_WORLD, MPI types.

      IMPLICIT NONE
C
      INCLUDE "mpif.h"
C
C The following declaration specification should match C type long int.
      INTEGER*8 NLOCAL, NEQ, I, IOUT(25), IPAR(2)
      PARAMETER (NLOCAL=10)
C
      INTEGER*4 LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNPE, LNLI, LNPS
      INTEGER*4 LNCFL, NOUT, MYPE, NPES, IER, METH, IATOL
      INTEGER*4 ITASK, IPRE, IGS, JOUT
      INTEGER*8 NST, NFE, NPSET, NPE, NPS, NNI, NLI
      INTEGER*8 NCFL, NETF, NCFN
      DOUBLE PRECISION Y(1024), ROUT(10), RPAR(1)
      DOUBLE PRECISION ATOL, DTOUT, T, ALPHA, RTOL, TOUT, ERMAX, ERRI
      DOUBLE PRECISION GERMAX, AVDIM
C
      DATA ATOL/1.0D-10/, RTOL/1.0D-5/, DTOUT/0.1D0/, NOUT/10/
      DATA LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/, LNNI/7/, LNSETUP/8/, 
     1     LNPE/20/, LNLI/22/, LNPS/21/, LNCFL/23/
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
      METH = 2
      IATOL = 1
      ITASK = 1
      IPRE = 1
      IGS = 1
C     Set parameter alpha
      ALPHA  = 10.0D0
C
C     Load IPAR and RPAR
      IPAR(1) = NLOCAL
      IPAR(2) = MYPE
      RPAR(1) = ALPHA
C
C     Do remaining initializations for first case: IPRE = 1 (prec. on left).
C
      DO 10 I = 1, NLOCAL
  10    Y(I) = 1.0D0
C
      IF (MYPE .EQ. 0) THEN
        WRITE(6,11) NEQ, ALPHA
  11    FORMAT('Diagonal test problem:'//' NEQ = ', I3,
     1         ' parameter alpha = ', F8.3)
        WRITE(6,12)
  12    FORMAT(' ydot_i = -alpha*i * y_i (i = 1,...,NEQ)')
        WRITE(6,13) RTOL, ATOL
  13    FORMAT(' RTOL, ATOL = ', 2E10.1)
        WRITE(6,14)
  14    FORMAT(' Method is BDF/NEWTON/SPGMR'/
     1         ' Diagonal preconditioner uses approximate Jacobian')
        WRITE(6,15) NPES
  15    FORMAT(' Number of processors = ', I3)
        WRITE(6,16)
  16    FORMAT(//'Preconditioning on left'/)
        ENDIF
C
      CALL FNVINITP(MPI_COMM_WORLD, 1, NLOCAL, NEQ, IER)
C
      IF (IER .NE. 0) THEN
         WRITE(6,20) IER
  20     FORMAT(///' SUNDIALS_ERROR: FNVINITP returned IER = ', I5)
         CALL MPI_FINALIZE(IER)
         STOP
      ENDIF
C
C     initialize SPGMR linear solver module
      call FSUNSPGMRINIT(1, IPRE, 0, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,25) IER
 25      FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRINIT IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      call FSUNSPGMRSETGSTYPE(1, IGS, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,27) IER
 27      FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRSETGSTYPE IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      CALL FCVMALLOC(T, Y, METH, IATOL, RTOL, ATOL,
     1               IOUT, ROUT, IPAR, RPAR, IER)
C
      IF (IER .NE. 0) THEN
         WRITE(6,30) IER
  30     FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
C     attach linear solver module to CVLs interface
      CALL FCVLSINIT(IER)
      IF (IER .NE. 0) THEN
         WRITE(6,32) IER
 32      FORMAT(///' SUNDIALS_ERROR: FCVLSINIT returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C     
C     attach preconditioner to CVLs interface
      CALL FCVLSSETPREC(1, IER)
C
C     Loop through tout values, call solver, print output, test for failure.
      TOUT = DTOUT
      DO 70 JOUT = 1, NOUT
C
        CALL FCVODE(TOUT, T, Y, ITASK, IER)
C
        IF (MYPE .EQ. 0) WRITE(6,40) T, IOUT(LNST), IOUT(LNFE)
  40    FORMAT(' t = ', E10.2, 5X, 'no. steps = ', I5,
     &         '   no. f-s = ', I5)
C
        IF (IER .NE. 0) THEN
           WRITE(6,60) IER, IOUT(15)
  60       FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /,
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
  75    ERMAX = MAX(ERMAX, ABS(ERRI))
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
      NFE = IOUT(LNFE)
      NPSET = IOUT(LNSETUP)
      NPE = IOUT(LNPE)
      NPS = IOUT(LNPS)
      NNI = IOUT(LNNI)
      NLI = IOUT(LNLI)
      AVDIM = DBLE(NLI) / DBLE(NNI)
      NCFN = IOUT(LNCF)
      NCFL = IOUT(LNCFL)
      NETF = IOUT(LNETF)
      IF (MYPE .EQ. 0)
     1  WRITE (6,90) NST, NFE, NPSET, NPE, NPS, NNI, NLI, AVDIM, NCFN,
     &               NCFL, NETF
  90  FORMAT(/'Final statistics:'//
     &       ' number of steps        = ', I5, 5X,
     &       'number of f evals.     =', I5/
     &       ' number of prec. setups = ', I5/
     &       ' number of prec. evals. = ', I5, 5X,
     &       'number of prec. solves = ', I5/
     &       ' number of nonl. iters. = ', I5, 5X,
     &       'number of lin. iters.  = ', I5/
     &       ' average Krylov subspace dimension (NLI/NNI)  = ', F8.4/
     &       ' number of conv. failures.. nonlinear = ', I3,
     &       '  linear = ', I3/
     &       ' number of error test failures = ', I3)
C
C     Re-initialize to run second case: IPRE = 2 (prec. on right).
      IPRE = 2
      T = 0.0D0
      DO 110 I = 1, NLOCAL
 110    Y(I) = 1.0D0
C
      IF (MYPE .EQ. 0)  WRITE(6,111) 
 111  FORMAT(//60('-')///'Preconditioning on right'/)
C
      CALL FCVREINIT(T, Y, IATOL, RTOL, ATOL, IER)
C
      IF (IER .NE. 0) THEN
         WRITE(6,130) IER
 130     FORMAT(///' SUNDIALS_ERROR: FCVREINIT returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      CALL FSUNSPGMRSETPRECTYPE(1, IPRE, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,140) IER
 140     FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRSETPRECTYPE IER = ',I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
C     Loop through tout values, call solver, print output, test for failure.
      TOUT = DTOUT
      DO 170 JOUT = 1, NOUT
C
        CALL FCVODE(TOUT, T, Y, ITASK, IER)
C
        IF (MYPE .EQ. 0) WRITE(6,40) T, IOUT(LNST), IOUT(LNFE)
C
        IF (IER .NE. 0) THEN
           WRITE(6,60) IER
           CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
           STOP
        ENDIF
C
        TOUT = TOUT + DTOUT
 170    CONTINUE
C
C     Get max. absolute error in the local vector.
      ERMAX = 0.0D0
      DO 175 I = 1, NLOCAL
        ERRI  = Y(I) - EXP(-ALPHA * (MYPE * NLOCAL + I) * T)
 175    ERMAX = MAX(ERMAX, ABS(ERRI))
C     Get global max. error from MPI_REDUCE call.
      CALL MPI_REDUCE(ERMAX, GERMAX, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
     1                0, MPI_COMM_WORLD, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,80) IER
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
      IF (MYPE .EQ. 0) WRITE(6,85) GERMAX
C     
C     Print final statistics.
      NST = IOUT(LNST)
      NFE = IOUT(LNFE)
      NPSET = IOUT(LNSETUP)
      NPE = IOUT(LNPE)
      NPS = IOUT(LNPS)
      NNI = IOUT(LNNI)
      NLI = IOUT(LNLI)
      AVDIM = DBLE(NLI) / DBLE(NNI)
      NCFN = IOUT(LNCF)
      NCFL = IOUT(LNCFL)
      NETF = IOUT(LNETF)
      IF (MYPE .EQ. 0)
     1  WRITE (6,90) NST, NFE, NPSET, NPE, NPS, NNI, NLI, AVDIM, NCFN,
     &               NCFL, NETF
C
C     Free the memory and finalize MPI.
      CALL FCVFREE
      CALL MPI_FINALIZE(IER)
      IF (IER .NE. 0) THEN
         WRITE(6,195) IER
 195     FORMAT(///' MPI_ERROR: MPI_FINALIZE returned IER = ', I5)
         STOP
      ENDIF
C
      STOP
      END
C
C     ------------------------------------------------------------------------
C
      SUBROUTINE FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
C     Routine for right-hand side function f
      IMPLICIT NONE
C
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*), MYPE, I, NLOCAL
      INTEGER*4 IER
      DOUBLE PRECISION T, Y(*), YDOT(*), RPAR(*)
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
C
C     ------------------------------------------------------------------------
C
      SUBROUTINE FCVPSOL(T, Y, FY, R, Z, GAMMA, DELTA, LR,
     &                   IPAR, RPAR, IER)
C     Routine to solve preconditioner linear system
C     This routine uses a diagonal preconditioner P = I - gamma*J,
C     where J is a diagonal approximation to the true Jacobian, given by:
C     J = diag(0, 0, 0, -4*alpha, ..., -N*alpha).
C     The vector r is copied to z, and the inverse of P (restricted to the
C     local vector segment) is applied to the vector z.
      IMPLICIT NONE
C
C The following declaration specification should match C type long int.
      INTEGER*8 IPAR(*), NLOCAL, I, MYPE, ISTART, IBASE
      INTEGER*4 IER, LR
      DOUBLE PRECISION T, Y(*), FY(*), R(*), Z(*)
      DOUBLE PRECISION GAMMA, DELTA, RPAR(*)
      DOUBLE PRECISION PSUBI, ALPHA
C
      NLOCAL = IPAR(1)
      MYPE = IPAR(2)
      ALPHA = RPAR(1)
C
      DO I = 1, NLOCAL
         Z(I) = R(I)
      ENDDO
C
      IBASE = MYPE * NLOCAL
      ISTART = MAX(1, 4 - IBASE)
      DO I = ISTART, NLOCAL
         PSUBI = 1.0D0 + GAMMA * ALPHA * (IBASE + I)
         Z(I) = Z(I) / PSUBI
      ENDDO
C
      RETURN
      END
C
C     ------------------------------------------------------------------------
C
      SUBROUTINE FCVPSET(T, Y, FY, JOK, JCUR, GAMMA, H,
     &                   IPAR, RPAR, IER)
C     Empty subroutine. Not needed for the preconditioner, but required
C     by the FCVODE module.
      INTEGER*4 IER, JOK, JCUR
      DOUBLE PRECISION T, Y(*), FY(*), GAMMA, H
C The following declaration specification should match C type long int
      INTEGER*8 IPAR(*)
      DOUBLE PRECISION RPAR(*)

      IER = 0
      RETURN
      END
