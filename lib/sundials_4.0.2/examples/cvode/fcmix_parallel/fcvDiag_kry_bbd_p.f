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
C     Diagonal ODE example. Stiff case, with diagonal preconditioner.
C     Uses FCVODE interfaces and FCVBBD interfaces.
C     Solves problem twice -- with left and right preconditioning.
C     --------------------------------------------------------------------
C
C     Include MPI-Fortran header file for MPI_COMM_WORLD, MPI types.

      IMPLICIT NONE
C
      INCLUDE "mpif.h"
C
C The following declaration specification should match C type long int.
      INTEGER*8 NLOCAL, NEQ, I, IOUT(25), IPAR(2), MUDQ, MLDQ, MU, ML
      PARAMETER (NLOCAL=10)   
C
      INTEGER*4 NOUT, LNST, LNFE, LNSETUP, LNNI, LNCF, LNETF, LNPE
      INTEGER*4 LNLI, LNPS, LNCFL, MYPE, IER, NPES, METH
      INTEGER*4 LLENRW, LLENIW, LLENRWLS, LLENIWLS
      INTEGER*4 IATOL, ITASK, IPRE, IGS, JOUT
C The following declaration specification should match C type long int.
      INTEGER*8 LENRWBBD, LENIWBBD, NGEBBD
      INTEGER*8 NST, NFE, NPSET, NPE, NPS, NNI, NLI, NCFN, NCFL, NETF
      INTEGER*8 LENRW, LENIW, LENRWLS, LENIWLS
      DOUBLE PRECISION Y(1024), ROUT(10), RPAR(1)
      DOUBLE PRECISION ALPHA, TOUT, ERMAX, AVDIM
      DOUBLE PRECISION ATOL, ERRI, RTOL, GERMAX, DTOUT, T
C     
      DATA ATOL/1.0D-10/, RTOL/1.0D-5/, DTOUT/0.1D0/, NOUT/10/
      DATA LLENRW/1/, LLENIW/2/, LNST/3/, LNFE/4/, LNETF/5/,  LNCF/6/,
     1     LNNI/7/, LNSETUP/8/, LLENRWLS/13/, LLENIWLS/14/,
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
      DO I = 1, NLOCAL
         Y(I) = 1.0D0
      ENDDO
C     
      IF (MYPE .EQ. 0) THEN
         WRITE(6,15) NEQ, ALPHA, RTOL, ATOL, NPES
 15      FORMAT('Diagonal test problem:'//' NEQ = ', I3, /
     &          ' parameter alpha = ', F8.3/
     &          ' ydot_i = -alpha*i * y_i (i = 1,...,NEQ)'/
     &          ' RTOL, ATOL = ', 2E10.1/
     &          ' Method is BDF/NEWTON/SPGMR'/
     &          ' Preconditioner is band-block-diagonal, using CVBBDPRE'
     &          /' Number of processors = ', I3/)
      ENDIF
C     
      CALL FNVINITP(MPI_COMM_WORLD, 1, NLOCAL, NEQ, IER)
C     
      IF (IER .NE. 0) THEN
         WRITE(6,20) IER
 20      FORMAT(///' SUNDIALS_ERROR: FNVINITP returned IER = ', I5)
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
     &               IOUT, ROUT, IPAR, RPAR, IER)
C     
      IF (IER .NE. 0) THEN
         WRITE(6,30) IER
 30      FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
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
      MUDQ = 0
      MLDQ = 0
      MU = 0
      ML = 0
      CALL FCVBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, 0.0D0, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,35) IER
 35      FORMAT(///' SUNDIALS_ERROR: FCVBBDINIT returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      IF (MYPE .EQ. 0) WRITE(6,38)
 38   FORMAT(/'Preconditioning on left'/)
C     
C     Looping point for cases IPRE = 1 and 2.
C     
 40   CONTINUE
C     
C     Loop through tout values, call solver, print output, test for failure.
      TOUT = DTOUT
      DO 60 JOUT = 1, NOUT
C     
         CALL FCVODE(TOUT, T, Y, ITASK, IER)
C     
         IF (MYPE .EQ. 0) WRITE(6,45) T, IOUT(LNST), IOUT(LNFE)
 45      FORMAT(' t = ', E10.2, 5X, 'no. steps = ', I5,
     &          '   no. f-s = ', I5)
C     
         IF (IER .NE. 0) THEN
            WRITE(6,50) IER, IOUT(15)
 50         FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /,
     &             '                 Linear Solver returned IER = ', I5)
            CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
            STOP
         ENDIF
C     
         TOUT = TOUT + DTOUT
 60   CONTINUE
C     
C     Get max. absolute error in the local vector.
      ERMAX = 0.0D0
      DO 65 I = 1, NLOCAL
         ERRI  = Y(I) - EXP(-ALPHA * (MYPE * NLOCAL + I) * T)
         ERMAX = MAX(ERMAX, ABS(ERRI))
 65   CONTINUE
C     Get global max. error from MPI_REDUCE call.
      CALL MPI_REDUCE(ERMAX, GERMAX, 1, MPI_DOUBLE_PRECISION, MPI_MAX,
     &                0, MPI_COMM_WORLD, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,70) IER
 70      FORMAT(///' MPI_ERROR: MPI_REDUCE returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
      IF (MYPE .EQ. 0) WRITE(6,75) GERMAX
 75   FORMAT(/'Max. absolute error is', E10.2/)
C     
C     Print final statistics.
      IF (MYPE .EQ. 0) THEN
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
         LENRW = IOUT(LLENRW)
         LENIW = IOUT(LLENIW)
         LENRWLS = IOUT(LLENRWLS)
         LENIWLS = IOUT(LLENIWLS)
         WRITE(6,80) NST, NFE, NPSET, NPE, NPS, NNI, NLI, AVDIM, NCFN,
     &               NCFL, NETF, LENRW, LENIW, LENRWLS, LENIWLS
 80      FORMAT(/'Final statistics:'//
     &          ' number of steps        = ', I5, 4X,
     &          ' number of f evals.     = ', I5/
     &          ' number of prec. setups = ', I5/
     &          ' number of prec. evals. = ', I5, 4X,
     &          ' number of prec. solves = ', I5/
     &          ' number of nonl. iters. = ', I5, 4X,
     &          ' number of lin. iters.  = ', I5/
     &          ' average Krylov subspace dimension (NLI/NNI) = ',F8.4/
     &          ' number of conv. failures.. nonlinear = ', I3,
     &          '  linear = ', I3/
     &          ' number of error test failures = ', I3/
     &          ' main solver real/int workspace sizes   = ',2I5/
     &          ' linear solver real/int workspace sizes = ',2I5)
         CALL FCVBBDOPT(LENRWBBD, LENIWBBD, NGEBBD)
         WRITE(6,82) LENRWBBD, LENIWBBD, NGEBBD
 82      FORMAT('In CVBBDPRE:'/
     &          ' real/int local workspace = ', 2I5/
     &          ' number of g evals. = ', I5)
      ENDIF
C     
C     If IPRE = 1, re-initialize T, Y, and the solver, and loop for
C     case IPRE = 2.  Otherwise jump to final block.
      IF (IPRE .EQ. 2) GO TO 99
C
      T = 0.0D0
      DO I = 1, NLOCAL
         Y(I) = 1.0D0
      ENDDO         
C
      CALL FCVREINIT(T, Y, IATOL, RTOL, ATOL, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,91) IER
 91      FORMAT(///' SUNDIALS_ERROR: FCVREINIT returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      IPRE = 2
C
      CALL FCVBBDREINIT(MUDQ, MLDQ, 0.0D0, IER) 
      IF (IER .NE. 0) THEN
         WRITE(6,92) IER
 92      FORMAT(///' SUNDIALS_ERROR: FCVBBDREINIT returned IER = ', I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      CALL FSUNSPGMRSETPRECTYPE(1, IPRE, IER)
      IF (IER .NE. 0) THEN
         WRITE(6,93) IER
 93      FORMAT(///' SUNDIALS_ERROR: FSUNSPGMRSETPRECTYPE IER = ',I5)
         CALL MPI_ABORT(MPI_COMM_WORLD, 1, IER)
         STOP
      ENDIF
C
      IF (MYPE .EQ. 0) WRITE(6,95)
 95   FORMAT(//60('-')///'Preconditioning on right'/)
      GO TO 40
C     
C     Free the memory and finalize MPI.
 99   CALL FCVFREE
      CALL MPI_FINALIZE(IER)
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
      SUBROUTINE FCVGLOCFN(NLOC, T, YLOC, GLOC, IPAR, RPAR, IER)
C     Routine to define local approximate function g, here the same as f. 
      IMPLICIT NONE
C
C The following declaration specification should match C type long int.
      INTEGER*8 NLOC, IPAR(*)
      INTEGER*4 IER
      DOUBLE PRECISION T, YLOC(*), GLOC(*), RPAR(*)
C
      CALL FCVFUN(T, YLOC, GLOC, IPAR, RPAR, IER)
C
      RETURN
      END
C
C     ------------------------------------------------------------------------
C
      SUBROUTINE FCVCOMMFN(NLOC, T, YLOC, IPAR, RPAR, IER)
C     Routine to perform communication required for evaluation of g.
      INTEGER*8 NLOC, IPAR(*)
      INTEGER*4 IER
      DOUBLE PRECISION T, YLOC(*), RPAR(*)
      IER = 0
      RETURN
      END
