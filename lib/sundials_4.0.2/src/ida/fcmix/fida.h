/*---------------------------------------------------------------
 * Programmer(s): Aaron Collier and Radu Serban @ LLNL
 *                Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This is the header file for FIDA, the Fortran interface to
 * the IDA package.
 *--------------------------------------------------------------*/

/*=============================================================================
                   FIDA Interface Package

  The FIDA Interface Package is a package of C functions which support
  the use of the IDA solver, for the solution of DAE systems, in a
  mixed Fortran/C setting.  While IDA is written in C, it is assumed
  here that the user's calling program and user-supplied problem-defining
  routines are written in Fortran.  This package provides the necessary
  interface to IDA for any acceptable NVECTOR implementation.

  A summary of the user-callable functions, with the corresponding
  IDA functions, are as follows:

    Fortran                      IDA
    ---------------------        --------------------------------
    FNVINITS                     N_VNew_Serial
    FNVINITP                     N_VNew_Parallel
    FNVINITOMP                   N_VNew_OpenMP
    FNVINITPTS                   N_VNew_Pthreads

    FSUNBANDMATINIT              SUNBandMatrix
    FSUNDENSEMATINIT             SUNDenseMatrix
    FSUNSPARSEMATINIT            SUNSparseMatrix

    FSUNBANDLINSOLINIT           SUNBandLinearSolver
    FSUNDENSELINSOLINIT          SUNDenseLinearSolver
    FSUNKLUINIT                  SUNKLU
    FSUNKLUREINIT                SUNKLUReinit
    FSUNLAPACKBANDINIT           SUNLapackBand
    FSUNLAPACKDENSEINIT          SUNLapackDense
    FSUNPCGINIT                  SUNPCG
    FSUNSPBCGSINIT               SUNSPBCGS
    FSUNSPFGMRINIT               SUNSPFGMR
    FSUNSPGMRINIT                SUNSPGMR
    FSUNSPTFQMRINIT              SUNSPTFQMR
    FSUNSUPERLUMTINIT            SUNSuperLUMT

    FIDAMALLOC                   IDACreate, IDASetUserData and IDAInit
    FIDAREINIT                   IDAReInit

    FIDASETIIN                   IDASet* (integer arguments)
    FIDASETRIN                   IDASet* (real arguments)
    FIDASETVIN                   IDASet* (vector arguments)

    FIDATOLREINIT                IDASetTolerances

    FIDACALCIC                   IDACalcIC

    FIDAEWTSET                   IDAWFtolerances

    FIDALSINIT                   IDASetLinearSolver
    FIDALSSETEPSLIN              IDASetEpsLin
    FIDALSSETINCREMENTFACTOR     IDASetIncrementFactor
    FIDALSSETJAC                 IDASetJacTimes
    FIDALSSETPREC                IDASetPreconditioner
    FIDADENSESETJAC              IDASetJacFn
    FIDABANDSETJAC               IDASetJacFn
    FIDASPARSESETJAC             IDASetJacFn

    FIDANLSINIT                  IDASetNonlinearSolver

    FIDASOLVE                    IDASolve, IDAGet*, and IDA*Get*

    FIDAGETDKY                   IDAGetDky

    FIDAGETERRWEIGHTS            IDAGetErrWeights

    FIDAGETESTLOCALERR           IDAGetEstLocalErrors

    FIDAFREE                     IDAFree
    ---------------------        --------------------------------

  The user-supplied functions, each listed with the corresponding interface
  function which calls it (and its type within IDA), are as follows:

    Fortran:           Interface Fcn:           IDA Type:
    -------------      ------------------       -----------------------
    FIDARESFUN         FIDAresfn                IDAResFn
    FIDADJAC           FIDADenseJac             IDALsJacFn
    FIDABJAC           FIDABandJac              IDALsJacFn
    FIDASPJAC          FIDASparseJac            IDALsJacFn
    FIDAPSET           FIDAPSet                 IDALsPrecSetupFn
    FIDAPSOL           FIDAPSol                 IDALsPrecSolveFn
    FIDAJTSETUP        FIDAJTSetup              IDALsJacTimesSetupFn
    FIDAJTIMES         FIDAJtimes               IDALsJacTimesVecFn
    FIDAEWT            FIDAEwtSet               IDAEwtFn
    -------------      ------------------       -----------------------

  In contrast to the case of direct use of IDA, the names of all user-supplied
  routines here are fixed, in order to maximize portability for the resulting
  mixed-language program.

  Important note on portability:
  In this package, the names of the interface functions, and the names of
  the Fortran user routines called by them, appear as dummy names
  which are mapped to actual values by a series of definitions, in this
  and other header files.

  =============================================================================

                   Usage of the FIDA Interface Package

  The usage of FIDA requires calls to a few different interface
  functions, depending on the method options selected, and one or more
  user-supplied routines which define the problem to be solved.  These
  function calls and user routines are summarized separately below.

  Some details are omitted, and the user is referred to the user documents
  on IDA for more complete documentation.  Information on the
  arguments of any given user-callable interface routine, or of a given
  user-supplied function called by an interface function, can be found in
  the documentation on the corresponding function in the IDA package.

  The number labels on the instructions below end with s for instructions
  that are specific to use with the serial/OpenMP/PThreads NVector package,
  and end with p are specific to use with the N_VParallel package.

 -----------------------------------------------------------------------------

                               Data Types

 Throughout this documentation, we will refer to data types according to
 their usage in SUNDIALS.  The equivalent types to these may vary,
 depending on your computer architecture and on how SUNDIALS was compiled.
 A Fortran user should take care that all arguments passed through this
 Fortran/C interface are declared of the appropriate type.

 Integers: SUNDIALS uses 'int', 'long int' and 'sunindextype' types.  At
 compilation, SUNDIALS allows the configuration of the 'index' type, that
 accepts values of 32-bit signed and 64-bit signed.  This choice dictates
 the size of a SUNDIALS 'sunindextype' variable.
   int      -- equivalent to an INTEGER or INTEGER*4 in Fortran
   long int -- equivalent to an INTEGER*8 in Fortran (Linux/UNIX/OSX), or
               equivalent to an INTEGER in Windows
   sunindextype -- this will depend on the SUNDIALS configuration:
               32-bit -- equivalent to an INTEGER or INTEGER*4 in Fortran
               64-bit -- equivalent to an INTEGER*8 in Fortran

 Real numbers:  At compilation, SUNDIALS allows the configuration option
 '--with-precision', that accepts values of 'single', 'double' or
 'extended' (the default is 'double').  This choice dictates the size of a
 SUNDIALS 'realtype' variable.  The corresponding Fortran types for these
 'realtype' sizes are:
   single   -- equivalent to a REAL or REAL*4 in Fortran
   double   -- equivalent to a DOUBLE PRECISION or REAL*8 in Fortran
   extended -- equivalent to a REAL*16 in Fortran

  -----------------------------------------------------------------------------

  (1) User-supplied residual routine: FIDARESFUN

      The user must in all cases supply the following Fortran routine

        SUBROUTINE FIDARESFUN(T, Y, YP, R, IPAR, RPAR, IER)

      It must set the R array to F(t,y,y'), the residual of the DAE system.

      The arguments are:
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives [realtype, input]
        R    -- array containing DAE residuals [realtype, output]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful,
                  >0 if a recoverable error occurred,
                  <0 if an unrecoverable error ocurred.

  (2s) Optional user-supplied dense Jacobian approximation routine: FIDADJAC

      As an option when using the Dense or LapackDense linear solvers, the
      user may supply a routine that computes a dense approximation of the
      system Jacobian J = dF/dy' + c_j*dF/dy. If supplied, it must have the
      following form:

        SUBROUTINE FIDADJAC(NEQ, T, Y, YP, R, DJAC, CJ, EWT, H,
       1                    IPAR, RPAR, WK1, WK2, WK3, IER)

      This routine must compute the Jacobian and store it columnwise in DJAC.

      The arguments are:
        NEQ  -- number of rows in the matrix [long int, input]
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives [realtype, input]
        R    -- array containing DAE residuals [realtype, input]
        DJAC -- 2D array containing the jacobian entries [realtype of size
                (NEQ,NEQ), output]
        CJ   -- scalar in the system Jacobian proportional to inverse step
                size [realtype, input]
        EWT  -- array containing error weight vector [realtype, input]
        H    -- current step size [realtype, input]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        WK*  -- array containing temporary workspace of same size as Y
                [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful,
                  >0 if a recoverable error occurred,
                  <0 if an unrecoverable error ocurred.

  (2s) Optional user-supplied band Jacobian approximation routine: FIDABJAC

      As an option when using the Band or LapackBand linear solvers, the
      user may supply a routine that computes a band approximation of the
      system Jacobian J = dF/dy' + c_j*dF/dy. If supplied, it must have the
      following form:

        SUBROUTINE FIDABJAC(NEQ, MU, ML, MDIM, T, Y, YP, R, CJ, BJAC,
       1                    EWT, H, IPAR, RPAR, WK1, WK2, WK3, IER)

      This routine must load the MDIM by N array BJAC with the Jacobian
      matrix at the current (t,y,y') in band form.  Store in BJAC(k,j)
      the Jacobian element J(i,j) with k = i - j + MU + 1
      (k = 1 ... ML+MU+1) and j = 1 ... N.

      The arguments are:
        NEQ  -- number of rows in the matrix [long int, input]
        MU   -- upper half-bandwidth of the matrix [long int, input]
        ML   -- lower half-bandwidth of the matrix [long int, input]
        MDIM -- leading dimension of BJAC array [long int, input]
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives [realtype, input]
        R    -- array containing DAE residuals [realtype, input]
        DJAC -- 2D array containing the jacobian entries [realtype of size
                (NEQ,NEQ), output]
        CJ   -- scalar in the system Jacobian proportional to inverse step
                size [realtype, input]
        EWT  -- array containing error weight vector [realtype, input]
        H    -- current step size [realtype, input]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        WK*  -- array containing temporary workspace of same size as Y
                [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful,
                  >0 if a recoverable error occurred,
                  <0 if an unrecoverable error ocurred.

  (2s) User-supplied sparse Jacobian approximation routine: FIDASPJAC

      When using the KLU or SuperLUMT linear solvers, the user *must* supply
      a routine that computes a compressed-sparse-column [or
      compressed-sparse-row] approximation of the system Jacobian
      J = dF/dy' + c_j*dF/dy.  If supplied, it must have the following form:

        SUBROUTINE FIDASPJAC(T, CJ, Y, YP, R, N, NNZ, JDATA, JRVALS,
       1                     JCPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)

      It must load the N by N compressed sparse column [row] matrix with
      storage for NNZ nonzeros, stored in the arrays JDATA (nonzero values),
      JRVALS (row [column] indices for each nonzero), JCOLPTRS (indices for
      start of each column [row]), with the Jacobian matrix in CSC [CSR]
      form (see sunmatrix_sparse.h for more information).

      The arguments are:
          T    -- current time [realtype, input]
          CJ   -- scalar in the system Jacobian proportional
                  to inverse step size [realtype, input]
          Y    -- array containing state variables [realtype, input]
          YP   -- array containing state derivatives [realtype, input]
          R    -- array containing system residual F(T, Y, YP) [realtype, input]
          N    -- number of matrix rows/columns in Jacobian [int, input]
          NNZ  -- allocated length of nonzero storage [int, input]
          JDATA -- nonzero values in Jacobian
                  [realtype of length NNZ, output]
          JRVALS -- row [column] indices for each nonzero in Jacobian
                   [int of length NNZ, output]
          JCPTRS -- pointers to each Jacobian column [row] in preceding arrays
                  [int of length N+1, output]
          H    -- current step size [realtype, input]
          IPAR -- array containing integer user data that was passed to
                  FIDAMALLOC [long int, input]
          RPAR -- array containing real user data that was passed to
                  FIDAMALLOC [realtype, input]
          WK*  -- array containing temporary workspace of same size as Y
                  [realtype, input]
          IER  -- return flag [int, output]:
                     0 if successful,
                    >0 if a recoverable error occurred,
                    <0 if an unrecoverable error ocurred.

      NOTE: this may ONLY be used if SUNDIALS has been configured with
      sunindextype set to 64-bit integers.

  (2) Optional user-supplied Jacobian-vector product setup routine:
      FIDAJTSETUP

      As an option when using the IDALS linear solver interface with a 
      matrix-free linear solver module, the user may supply a routine that 
      computes the product of the system Jacobian J = dF/dy' + c_j*dF/dy 
      and a given vector v, as well as a routine to set up any user data 
      structures in preparation for the matrix-vector product.  If a 
      'setup' routine is supplied, it must have the following form:

        SUBROUTINE FIDAJTSETUP(T, Y, YP, R, CJ, EWT, H, IPAR, RPAR, IER)

      It must perform any relevant preparations for subsequent calls to the
      user-provided FIDAJTIMES routine (see below).

      The arguments are:
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives [realtype, input]
        R    -- array containing DAE residuals [realtype, input]
        CJ   -- scalar in the system Jacobian proportional to inverse step
                size [realtype, input]
        EWT  -- array containing error weight vector [realtype, input]
        H    -- current step size [realtype, input]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful,
                   nonzero if an error.

  (2) Optional user-supplied Jacobian-vector product routine: FIDAJTIMES

      As an option when using the IDALS linear solver interface with a 
      matrix-free linear solver module, the user may supply a routine 
      that computes the product of the system Jacobian 
      J = dF/dy' + c_j*dF/dy and a given vector v.  If supplied, it must 
      have the following form:

         SUBROUTINE FIDAJTIMES(T, Y, YP, R, V, FJV, CJ, EWT, H,
        1                      IPAR, RPAR, WK1, WK2, IER)

      This routine must compute the product vector Jv, where the vector v
      is stored in V, and store the product in FJV.

      The arguments are:
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives [realtype, input]
        R    -- array containing DAE residuals [realtype, input]
        V    -- array containing vector to multiply [realtype, input]
        FJV  -- array containing product vector [realtype, output]
        CJ   -- scalar in the system Jacobian proportional to inverse step
                size [realtype, input]
        EWT  -- array containing error weight vector [realtype, input]
        H    -- current step size [realtype, input]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful,
                   nonzero if an error.

  (3) Optional user-supplied preconditioner setup/solve routines: FIDAPSET
      and FIDAPSOL

      As an option when using the IDALS linear solver interface and an 
      iterative linear solver module, the user may supply routines to 
      setup and apply the preconditioner.  If supplied, these must have 
      the following form:

        SUBROUTINE FIDAPSET(T, Y, YP, R, CJ, EWT, H, IPAR, RPAR, IER)

      This routine must perform any evaluation of Jacobian-related data and
      preprocessing needed for the solution of the preconditioner linear
      systems by FIDAPSOL.

      The arguments are:
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives [realtype, input]
        R    -- array containing DAE residuals [realtype, input]
        CJ   -- scalar in the system Jacobian proportional to inverse step
                size [realtype, input]
        EWT  -- array containing error weight vector [realtype, input]
        H    -- current step size [realtype, input]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful,
                   nonzero if an error.

      The user-supplied routine FIDAPSOL must have the form:

         SUBROUTINE FIDAPSOL(T, Y, YP, R, RV, ZV, CJ, DELTA, EWT,
        1                    IPAR, RPAR, IER)

      This routine must solve the preconditioner linear system Pz = r,
      where r = RV is input, and store the solution z in ZV.

      The arguments are:
        T    -- current time [realtype, input]
        Y    -- array containing state variables [realtype, input]
        YP   -- array containing state variable derivatives [realtype, input]
        R    -- array containing DAE residuals [realtype, input]
        RV   -- right-hand side array [realtype, input]
        ZV   -- solution array [realtype, output]
        CJ   -- scalar in the system Jacobian proportional to inverse step
                size [realtype, input]
        DELTA -- desired residual tolerance [realtype, input]
        EWT  -- array containing error weight vector [realtype, input]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful,
                   nonzero if an error.

  (4) Optional user-supplied error weight vector routine: FIDAEWT

      As an option to providing the relative and absolute tolerances, the
      user may supply a routine that computes the weights used in the WRMS
      norms.  If supplied, it must have the following form:

        SUBROUTINE FIDAEWT(Y, EWT, IPAR, RPAR, IER)

      It must store the error weights in EWT, given the current solution
      vector Y.

      The arguments are:
        Y    -- array containing state variables [realtype, input]
        EWT  -- array containing the error weight vector [realtype, output]
        IPAR -- array containing integer user data that was passed to
                FIDAMALLOC [long int, input]
        RPAR -- array containing real user data that was passed to
                FIDAMALLOC [realtype, input]
        IER  -- return flag [int, output]:
                   0 if successful,
                   nonzero if an error.

  -----------------------------------------------------------------------------

  (5) Initialization:  FNVINITS / FNVINITP / FNVINITOMP / FNVINITPTS,
                       FSUNBANDMATINIT / FSUNDENSEMATINIT /
                          FSUNSPARSEMATINIT,
                       FSUNBANDLINSOLINIT / FSUNDENSELINSOLINIT /
                          FSUNKLUINIT / FSUNKLUREINIT / FSUNKLUSETORDERING /
                          FSUNLAPACKBANDINIT / FSUNLAPACKDENSEINIT /
                          FSUNPCGINIT / FSUNSPBCGSINIT / FSUNSPFGMRINIT /
                          FSUNSPGMRINIT / FSUNSPTFQMRINIT / FSUNSUPERLUMTINIT /
                          FSUNSUPERLUMTSETORDERING,
                       FIDAMALLOC,
                       FIDALSINIT
                       FIDAREINIT,
                       FIDATOLREINIT,
                       FIDACALCIC,

      NOTE: the initialization order is important!  It *must* proceed as
      shown: vector, matrix (if used), linear solver (if used), IDA,
      IDALS, reinit.

  (5.1) To initialize the a vector specification for storing the solution
      data, the user must make one of the following calls:

        (serial)
           CALL FNVINITS(2, NEQ, IER)
        (MPI parallel)
           CALL FNVINITP(COMM, 2, NLOCAL, NGLOBAL, IER)
        (OpenMP threaded)
           CALL FNVINITOMP(2, NEQ, NUM_THREADS, IER)
        (PThreads threaded)
           CALL FNVINITPTS(2, NEQ, NUM_THREADS, IER)

      In each of these, one argument is an int containing the IDA solver
      ID (2).

      The other arguments are:
         NEQ = size of vectors [long int, input]
         COMM = the MPI communicator [int, input]
         NLOCAL = local size of vectors on this processor
            [long int, input]
         NGLOBAL = the system size, and the global size of vectors (the sum
            of all values of NLOCAL) [long int, input]
         NUM_THREADS = number of threads
         IER = return completion flag [int, output]:
                  0 = success,
                 -1 = failure.

  (5.2) To initialize a band/dense/sparse matrix structure for
      storing the system Jacobian and for use within a direct linear solver,
      the user must make one of the following calls:

           CALL FSUNBANDMATINIT(2, N, MU, ML, SMU, IER)
           CALL FSUNDENSEMATINIT(2, M, N, IER)
           CALL FSUNSPARSEMATINIT(2, M, N, NNZ, SPARSETYPE, IER)

      In each of these, one argument is an int containing the IDA solver
      ID (2).

      The other arguments are:

         M = the number of rows of the matrix [long int, input]
         N = the number of columns of the matrix [long int, input]
         MU = the number of upper bands (diagonal not included) in a banded
            matrix [long int, input]
         ML = the number of lower bands (diagonal not included) in a banded
            matrix [long int, input]
         SMU = the number of upper bands to store (diagonal not included)
            for factorization of a banded matrix [long int, input]
         NNZ = the storage size (upper bound on the number of nonzeros) for
            a sparse matrix [long int, input]
         SPARSETYPE = integer denoting use of CSC (0) vs CSR (1) storage
            for a sparse matrix [int, input]
         IER = return completion flag [int, output]:
                  0 = success,
                 -1 = failure.

  (5.3) To initialize a linear solver structure for solving linear systems
      arising from solution to the DAE, the user must make one of the
      following calls:

           CALL FSUNBANDLINSOLINIT(2, IER)
           CALL FSUNDENSELINSOLINIT(2, IER)
           CALL FSUNKLUINIT(2, IER)
           CALL FSUNLAPACKBANDINIT(2, IER)
           CALL FSUNLAPACKDENSEINIT(2, IER)
           CALL FSUNPCGINIT(2, PRETYPE, MAXL, IER)
           CALL FSUNSPBCGSINIT(2, PRETYPE, MAXL, IER)
           CALL FSUNSPFGMRINIT(2, PRETYPE, MAXL, IER)
           CALL FSUNSPGMRINIT(2, PRETYPE, MAXL, IER)
           CALL FSUNSPTFQMRINIT(2, PRETYPE, MAXL, IER)
           CALL FSUNSUPERLUMTINIT(2, NUM_THREADS, IER)

      Or once these have been initialized, their solver parameters may be
      modified via calls to the functions

           CALL FSUNKLUSETORDERING(2, ORD_CHOICE, IER)
           CALL FSUNSUPERLUMTSETORDERING(2, ORD_CHOICE, IER)

           CALL FSUNPCGSETPRECTYPE(2, PRETYPE, IER)
           CALL FSUNPCGSETMAXL(2, MAXL, IER)
           CALL FSUNSPBCGSSETPRECTYPE(2, PRETYPE, IER)
           CALL FSUNSPBCGSSETMAXL(2, MAXL, IER)
           CALL FSUNSPFGMRSETGSTYPE(2, GSTYPE, IER)
           CALL FSUNSPFGMRSETPRECTYPE(2, PRETYPE, IER)
           CALL FSUNSPGMRSETGSTYPE(2, GSTYPE, IER)
           CALL FSUNSPGMRSETPRECTYPE(2, PRETYPE, IER)
           CALL FSUNSPTFQMRSETPRECTYPE(2, PRETYPE, IER)
           CALL FSUNSPTFQMRSETMAXL(2, MAXL, IER)

      In all of the above, one argument is an int containing the IDA solver
      ID (2).

      The other arguments are:

         NNZ = the storage size (upper bound on the number of nonzeros) for
            a sparse matrix [long int, input]
         ORD_CHOICE = integer denoting ordering choice (see
            SUNKLUSetOrdering and SUNSuperLUMTSetOrdering documentation
            for details) [int, input]
         PRETYPE = type of preconditioning to perform (0=none, 1=left,
            2=right, 3=both) [int, input]
         MAXL = maximum Krylov subspace dimension [int, input]
         GSTYPE = choice of Gram-Schmidt orthogonalization algorithm
            (0=modified, 1=classical) [int, input]
         IER = return completion flag [int, output]:
                   0 = success,
                  -1 = failure.

  (5.4) To set various problem and solution parameters and allocate
      internal memory, make the following call:

         CALL FIDAMALLOC(T0, Y0, YP0, IATOL, RTOL, ATOL,
        1                IOUT, ROUT, IPAR, RPAR, IER)

      The arguments are:
         T0 = initial value of t [realtype, input]
         Y0 = array of initial conditions for y(t0) [realtype, input]
         YP0 = array of initial conditions for y'(t0) [realtype, input]
         IATOL = type for absolute tolerance ATOL [int, input]:
                   1 = scalar,
                   2 = array,
                   3 = user-supplied function; the user must supply a routine
                       FIDAEWT to compute the error weight vector.
         RTOL = scalar relative tolerance [realtype, input]
         ATOL = scalar or array absolute tolerance [realtype, input]
         IOUT = array of length at least 21 for integer optional outputs
                [long int, output]
         ROUT = array of length at least 6 for real optional outputs
                [realtype, output]
         IPAR = array of user integer data [long int, input/output]
         RPAR = array with user real data [realtype, input/output]
         IER  = return completion flag [int, output]:
                   0 = SUCCESS,
                  -1 = failure (see printed message for failure details).

      The user data arrays IPAR and RPAR are passed unmodified to all
      subsequent calls to user-provided routines. Modifications to either
      array inside a user-provided routine will be propagated. Using these
      two arrays, the user can dispense with Common blocks to pass data
      betwen user-provided routines.

      The optional outputs are:
            LENRW   = IOUT( 1) -> IDAGetWorkSpace
            LENIW   = IOUT( 2) -> IDAGetWorkSpace
            NST     = IOUT( 3) -> IDAGetNumSteps
            NRE     = IOUT( 4) -> IDAGetNumResEvals
            NETF    = IOUT( 5) -> IDAGetNumErrTestFails
            NCFN    = IOUT( 6) -> IDAGetNumNonlinSolvConvFails
            NNI     = IOUT( 7) -> IDAGetNumNonlinSolvIters
            NSETUPS = IOUT( 8) -> IDAGetNumLinSolvSetups
            KLAST   = IOUT( 9) -> IDAGetLastOrder
            KCUR    = IOUT(10) -> IDAGetCurrentOrder
            NBCKTRK = IOUT(11) -> IDAGetNumBacktrackOps
            NGE     = IOUT(12) -> IDAGetNumGEvals

            HINUSED = ROUT( 1) -> IDAGetActualInitStep
            HLAST   = ROUT( 2) -> IDAGetLastStep
            HCUR    = ROUT( 3) -> IDAGetCurrentStep
            TCUR    = ROUT( 4) -> IDAGetCurrentTime
            TOLSFAC = ROUT( 5) -> IDAGetTolScaleFactor
            UNITRND = ROUT( 6) -> UNIT_ROUNDOFF
      See the IDA manual for details.

  (5.5) To attach the linear solver created in step (5.3) to the 
      IDALS interface, using the command:

        CALL FIDALSINIT(IER)

      The arguments are:
            IER  = return completion flag [int, output]:
                   0 = SUCCESS,
                  -1 = failure (see printed message for failure details).

  (5.6) If the user program includes the FIDAEWT routine for the evaluation
      of the error weights, the following call must be made

        CALL FIDAEWTSET(FLAG, IER)

      with FLAG = 1 to specify that FIDAEWT is provided and should be used;
      FLAG = 0 resets to the default EWT formulation.
      The return flag IER is 0 if successful, and nonzero otherwise.

  (5.7) If the user program includes the FIDABJAC routine for the
      evaluation of the band approximation to the Jacobian, then following
      the call to FIDALSINIT, the following call must be made

        CALL FIDABANDSETJAC(FLAG, IER)

      with the int FLAG=1 to specify that FIDABJAC is provided and should be
      used; FLAG=0 specifies a reset to the internal finite difference
      Jacobian approximation.  The int return flag IER=0 if successful,
      nonzero otherwise.

      If the user program includes the FIDADJAC routine for the evaluation
      of the dense approximation to the Jacobian, then after the call to
      FIDALSINIT, the following call must be made

        CALL FIDADENSESETJAC(FLAG, IER)

      with the int FLAG=1 to specify that FIDADJAC is provided and should be
      used; FLAG=0 specifies a reset to the internal finite difference
      Jacobian approximation.  The int return flag IER=0 if successful, and
      nonzero otherwise.

      When using a sparse matrix and linear solver the user must provide the
      FIDASPJAC routine for the evaluation of the sparse approximation to
      the Jacobian.  To indicate that this routine has been provided, after
      the call to FIDALSINIT, the following call must be made

        CALL FIDASPARSESETJAC(IER)

      The int return flag IER=0 if successful, and nonzero otherwise.

  (5.8) If the user program includes the FIDAJTSETUP and FIDAJTIMES
      routines for setup of a Jacobian-times-vector product, then after 
      creating the IDALS interface, the following call must be made:

        CALL FIDALSSETJAC(FLAG, IER)

      with the int FLAG=1 to specify that FIDAJTSETUP and FIDAJTIMES are
      provided and should be used; FLAG=0 specifies a reset to the internal
      finite difference approximation to this product).  The int return
      flag IER=0 if successful, and nonzero otherwise.

  (5.9) If the user program includes the FIDAPSET and FIDAPSOL routines
      for supplying a preconditioner to an iterative linear solver, then
      after creating the IDALS interface, the following call must be made

        CALL FIDALSSETPREC(FLAG, IER)

      with the int FLAG=1.  If FLAG=0 then preconditioning with these
      routines will be disabled. The return flag IER=0 if successful,
      nonzero otherwise.

  (5.10) If the user wishes to use one of IDA's built-in preconditioning
      module, FIDABBD, then that should be initialized after creating the
      IDALS interface using the call

        CALL FIDABBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, DQRELY, IER)

      Detailed explanation of the inputs to these functions, as well as any
      requirements of user-supplied functions on which these preconditioning
      modules rely, may be found in the header file fidabbd.h.

  (5.11) To set various integer optional inputs, make the folowing call:

        CALL FIDASETIIN(KEY, VALUE, IER)

      to set the integer input VALUE to the optional input specified by the
      quoted character string KEY.  VALUE must be a Fortran integer of size
      commensurate with a C "long int".  KEY must be one of the following:
      MAX_ORD, MAX_NSTEPS, MAX_ERRFAIL, MAX_NITERS, MAX_CONVFAIL,
      SUPPRESS_ALG, MAX_NSTEPS_IC, MAX_NITERS_IC, MAX_NJE_IC, LS_OFF_IC.
      The int return flag IER is 0 if successful, and nonzero otherwise.

  (5.12) To set various real optional inputs, make the folowing call:

        CALL FIDASETRIN(KEY, VALUE, IER)

      to set the realtype value VALUE to the optional input specified by the
      quoted character string KEY.  VALUE must be a Fortran real-valued
      number of size commensurate with the SUNDIALS "realtype".  KEY must
      one of the following: INIT_STEP, MAX_STEP, MIIN_STEP, STOP_TIME,
      NLCONV_COEF. The int return flag IER is 0 if successful, and nonzero
      otherwise.

  (5.13) To set the vector of variable IDs or the vector of constraints,
      make the following call:

        CALL FIDASETVIN(KEY, ARRAY, IER)

      where ARRAY is an array of realtype and the quoted character string
      KEY is one of: ID_VEC or CONSTR_VEC.  The int return flag IER is 0
      if successful, and nonzero otherwise.

  (5.14) To re-initialize the FIDA solver for the solution of a new problem
      of the same size as one already solved, make the following call:

        CALL FIDAREINIT(T0, Y0, YP0, IATOL, RTOL, ATOL, ID, CONSTR, IER)

      The arguments have the same names and meanings as those of FIDAMALLOC.
      FIDAREINIT performs the same initializations as FIDAMALLOC, but does
      no memory allocation for IDA data structures, using instead the
      existing internal memory created by the previous FIDAMALLOC call.
      The subsequent calls to attach the linear system solver is only needed
      if the matrix or linear solver objects have been re-created.

  (5.15) To modify the tolerance parameters, make the following call:

        CALL FIDATOLREINIT(IATOL, RTOL, ATOL, IER)

      The arguments have the same names and meanings as those of FIDAMALLOC.
      FIDATOLREINIT simply calls IDASetTolerances with the given arguments.

  (5.16) To compute consistent initial conditions for an index-one DAE system,
      make the following call:

        CALL FIDACALCIC(ICOPT, TOUT, IER)

      The arguments are:
         ICOPT = specifies the option [int, input]:
                 1 = IDA_YP_YDP_INIT
                 2 = IDA_Y_INIT
                 (See user guide for additional details)
         TOUT  = the first value of t at which a solution will
                 be requested from FIDASOLVE [realtype, input].
         IER   = return completion flag [int, output].

  (5.17) The FSUNKLU solver will reuse much of the factorization information
      from one solve to the next.  If at any time the user wants to force a
      full refactorization or if the number of nonzeros in the Jacobian
      matrix changes, the user should make the call

         CALL FSUNKLUREINIT(2, NNZ, REINIT_TYPE, IER)

      The arguments are:
         NNZ = the maximum number of nonzeros [int; input]
         REINIT_TYPE = 1 or 2.  For a value of 1, the matrix will be
           destroyed and a new one will be allocated with NNZ nonzeros.
           For a value of 2, only symbolic and numeric factorizations will
           be completed.

  -----------------------------------------------------------------------------

  (6) Optional outputs from the IDALS linear solver interface (stored in the
      IOUT array that was passed to FIDAMALLOC)

         LENRWLS = IOUT(13) -> IDAGetLinWorkSpace
         LENIWLS = IOUT(14) -> IDAGetLinWorkSpace
         LSTF    = IOUT(15) -> IDAGetLastLinFlag
         NRELS   = IOUT(16) -> IDAGetNumLinResEvals
         NJE     = IOUT(17) -> IDAGetNumJacEvals
         NJTS    = IOUT(18) -> IDAGetJTSetupEvals
         NJT     = IOUT(19) -> IDAGetJtimesEvals
         NPE     = IOUT(20) -> IDAGetPrecEvals
         NPS     = IOUT(21) -> IDAGetPrecSolves
         NLI     = IOUT(22) -> IDAGetLinIters
         NLCF    = IOUT(23) -> IDAGetLinConvFails

      See the IDA manual for more detailed descriptions of any of the
      above.

  -----------------------------------------------------------------------------

  (7) The solver: FIDASOLVE

      To solve the DAE system, make the following call:

         CALL FIDASOLVE(TOUT, TRET, Y, YP, ITASK, IER)

      The arguments are:
        TOUT = next value of t at which a solution is desired [realtype, input]
        TRET = value of t reached by the solver [realtype, output]
        Y = array containing state variables on output [realtype, output]
        YP = array containing state derivatives on output [realtype, output]
        ITASK = task indicator [int, input]:
                   1 = normal mode (overshoot TOUT and interpolate)
                   2 = one-step mode (return after each internal step taken)
                   3 = normal tstop mode (like 1, but integration never
                       proceeds past TSTOP, which must be specified through a
                       call to FIDASETRIN using the key 'STOP_TIME')
                   4 = one step tstop (like 2, but integration never goes
                       past TSTOP)
        IER = completion flag [int, output]:
                   0 = success,
                   1 = tstop return,
                   2 = root return,
                   negative values are failure modes (see IDA manual).
      The current values of the optional outputs are immediately available in
      the IOUT and ROUT arrays.

  -----------------------------------------------------------------------------

  (8) Getting current solution derivative: FIDAGETDKY

      To obtain interpolated values of y and y' for any value of t in the
      last internal step taken by IDA, make the following call:

         CALL FIDAGETDKY(T, K, DKY, IER)

      The arguments are:
        T = time at which solution derivative is desired, within the interval
            [TCUR-HU,TCUR], [realtype, input].
        K = derivative order (0 .le. K .le. QU) [int, input]
        DKY = array containing computed K-th derivative of y [realtype, output]
        IER = return flag [int, output]: 0=success, <0 = illegal argument.

  -----------------------------------------------------------------------------

  (9) Get the current error weight vector: FIDAGETERRWEIGHTS

      To obtain the current error weight vector, make the following call:

        CALL FIDAGETERRWEIGHTS(EWT, IER)

      The arguments are:
        EWT = array containing the error weight vector [realtype, output]
        IER = return flag [int, output]: 0=success, nonzero if an error.

  -----------------------------------------------------------------------------

  (10) Get an estimate of the local error: FIDAGETESTLOCALERR

      To obtain the current error estimate vector, make the following call:

        CALL FIDAGETESTLOCALERR(ELE, IER)

      The arguments are:
        ELE = array with the estimated local error vector [realtype, output]
        IER = return flag [int, output]: 0=success, nonzero if an error.

  -----------------------------------------------------------------------------

  (11) Memory freeing: FIDAFREE

      To the free the internal memory created by the calls to FIDAMALLOC,
      FIDALSINIT, the generic linear solver and matrix modules,
      and FNVINIT*, make the following call:

        CALL FIDAFREE

 =============================================================================*/

#ifndef _FIDA_H
#define _FIDA_H

#include <ida/ida.h>                         /* definition of type IDAResFn */
#include <sundials/sundials_linearsolver.h>  /* definition of type SUNLinearSolver */
#include <sundials/sundials_matrix.h>        /* definition of type SUNMatrix */
#include <sundials/sundials_nvector.h>       /* definition of type N_Vector */
#include <sundials/sundials_types.h>         /* definition of type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)

#define FIDA_MALLOC         SUNDIALS_F77_FUNC(fidamalloc, FIDAMALLOC)
#define FIDA_REINIT         SUNDIALS_F77_FUNC(fidareinit, FIDAREINIT)
#define FIDA_SETIIN         SUNDIALS_F77_FUNC(fidasetiin, FIDASETIIN)
#define FIDA_SETRIN         SUNDIALS_F77_FUNC(fidasetrin, FIDASETRIN)
#define FIDA_SETVIN         SUNDIALS_F77_FUNC(fidasetvin, FIDASETVIN)
#define FIDA_TOLREINIT      SUNDIALS_F77_FUNC(fidatolreinit, FIDATOLREINIT)
#define FIDA_SOLVE          SUNDIALS_F77_FUNC(fidasolve, FIDASOLVE)
#define FIDA_FREE           SUNDIALS_F77_FUNC(fidafree, FIDAFREE)
#define FIDA_CALCIC         SUNDIALS_F77_FUNC(fidacalcic, FIDACALCIC)
#define FIDA_LSINIT         SUNDIALS_F77_FUNC(fidalsinit, FIDALSINIT)
#define FIDA_LSSETEPSLIN    SUNDIALS_F77_FUNC(fidalssetepslin, FIDALSSETEPSLIN)
#define FIDA_LSSETINCREMENTFACTOR SUNDIALS_F77_FUNC(fidalssetincrementfactor, FIDALSSETINCREMENTFACTOR)
#define FIDA_BANDSETJAC     SUNDIALS_F77_FUNC(fidabandsetjac, FIDABANDSETJAC)
#define FIDA_BJAC           SUNDIALS_F77_FUNC(fidabjac, FIDABJAC)
#define FIDA_DENSESETJAC    SUNDIALS_F77_FUNC(fidadensesetjac, FIDADENSESETJAC)
#define FIDA_DJAC           SUNDIALS_F77_FUNC(fidadjac, FIDADJAC)
#define FIDA_SPARSESETJAC   SUNDIALS_F77_FUNC(fidasparsesetjac, FIDASPARSESETJAC)
#define FIDA_SPJAC          SUNDIALS_F77_FUNC(fidaspjac, FIDASPJAC)
#define FIDA_LSSETJAC       SUNDIALS_F77_FUNC(fidalssetjac, FIDALSSETJAC)
#define FIDA_JTSETUP        SUNDIALS_F77_FUNC(fidajtsetup, FIDAJTSETUP)
#define FIDA_JTIMES         SUNDIALS_F77_FUNC(fidajtimes, FIDAJTIMES)
#define FIDA_LSSETPREC      SUNDIALS_F77_FUNC(fidalssetprec, FIDALSSETPREC)
#define FIDA_PSET           SUNDIALS_F77_FUNC(fidapset, FIDAPSET)
#define FIDA_PSOL           SUNDIALS_F77_FUNC(fidapsol, FIDAPSOL)
#define FIDA_RESFUN         SUNDIALS_F77_FUNC(fidaresfun, FIDARESFUN)
#define FIDA_EWTSET         SUNDIALS_F77_FUNC(fidaewtset, FIDAEWTSET)
#define FIDA_EWT            SUNDIALS_F77_FUNC(fidaewt, FIDAEWT)
#define FIDA_GETDKY         SUNDIALS_F77_FUNC(fidagetdky, FIDAGETDKY)
#define FIDA_GETERRWEIGHTS  SUNDIALS_F77_FUNC(fidageterrweights, FIDAGETERRWEIGHTS)
#define FIDA_GETESTLOCALERR SUNDIALS_F77_FUNC(fidagetestlocalerr, FIDAGETESTLOCALERR)
#define FIDA_NLSINIT        SUNDIALS_F77_FUNC(fidanlsinit, FIDANLSINIT)

/*** DEPRECATED ***/  
#define FIDA_DLSINIT        SUNDIALS_F77_FUNC(fidadlsinit, FIDADLSINIT)
#define FIDA_SPILSINIT      SUNDIALS_F77_FUNC(fidaspilsinit,FIDASPILSINIT)
#define FIDA_SPILSSETEPSLIN SUNDIALS_F77_FUNC(fidaspilssetepslin, FIDASPILSSETEPSLIN)
#define FIDA_SPILSSETINCREMENTFACTOR SUNDIALS_F77_FUNC(fidaspilssetincrementfactor, FIDASPILSSETINCREMENTFACTOR)
#define FIDA_SPILSSETJAC    SUNDIALS_F77_FUNC(fidaspilssetjac, FIDASPILSSETJAC)
#define FIDA_SPILSSETPREC   SUNDIALS_F77_FUNC(fidaspilssetprec, FIDASPILSSETPREC)
/******************/  
  
#else

#define FIDA_MALLOC         fidamalloc_
#define FIDA_REINIT         fidareinit_
#define FIDA_SETIIN         fidasetiin_
#define FIDA_SETRIN         fidasetrin_
#define FIDA_SETVIN         fidasetvin_
#define FIDA_TOLREINIT      fidatolreinit_
#define FIDA_SOLVE          fidasolve_
#define FIDA_FREE           fidafree_
#define FIDA_CALCIC         fidacalcic_
#define FIDA_LSINIT         fidalsinit_
#define FIDA_LSSETEPSLIN    fidalssetepslin_
#define FIDA_LSSETINCREMENTFACTOR fidalssetincrementfactor_
#define FIDA_BANDSETJAC     fidabandsetjac_
#define FIDA_BJAC           fidabjac_
#define FIDA_DENSESETJAC    fidadensesetjac_
#define FIDA_DJAC           fidadjac_
#define FIDA_SPARSESETJAC   fidasparsesetjac_
#define FIDA_SPJAC          fidaspjac_
#define FIDA_LSSETJAC       fidalssetjac_
#define FIDA_JTSETUP        fidajtsetup_
#define FIDA_JTIMES         fidajtimes_
#define FIDA_LSSETPREC      fidalssetprec_
#define FIDA_PSET           fidapset_
#define FIDA_PSOL           fidapsol_
#define FIDA_RESFUN         fidaresfun_
#define FIDA_EWTSET         fidaewtset_
#define FIDA_EWT            fidaewt_
#define FIDA_GETDKY         fidagetdky_
#define FIDA_GETERRWEIGHTS  fidageterrweights_
#define FIDA_GETESTLOCALERR fidagetestlocalerr_
#define FIDA_NLSINIT        fidanlsinit_

/*** DEPRECATED ***/  
#define FIDA_DLSINIT        fidadlsinit_
#define FIDA_SPILSINIT      fidaspilsinit_
#define FIDA_SPILSSETEPSLIN fidaspilssetepslin_
#define FIDA_SPILSSETINCREMENTFACTOR fidaspilssetincrementfactor_
#define FIDA_SPILSSETJAC    fidaspilssetjac_
#define FIDA_SPILSSETPREC   fidaspilssetprec_
/******************/  
  
#endif

/* Type for user data */

typedef struct {
  realtype *rpar;
  long int *ipar;
} *FIDAUserData;

/* Prototypes of exported functions */

void FIDA_MALLOC(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 long int *iout, realtype *rout,
                 long int *ipar, realtype *rpar, int *ier);
void FIDA_REINIT(realtype *t0, realtype *yy0, realtype *yp0,
                 int *iatol, realtype *rtol, realtype *atol,
                 int *ier);

void FIDA_SETIIN(char key_name[], long int *ival, int *ier);
void FIDA_SETRIN(char key_name[], realtype *rval, int *ier);
void FIDA_SETVIN(char key_name[], realtype *vval, int *ier);

void FIDA_TOLREINIT(int *iatol, realtype *rtol, realtype *atol, int *ier);
void FIDA_CALCIC(int *icopt, realtype *tout1, int *ier);

void FIDA_LSINIT(int *ier);
void FIDA_LSSETEPSLIN(realtype *eplifac, int *ier);
void FIDA_LSSETINCREMENTFACTOR(realtype *dqincfac, int *ier);
void FIDA_LSSETJAC(int *flag, int *ier);
void FIDA_LSSETPREC(int *flag, int *ier);
void FIDA_DENSESETJAC(int *flag, int *ier);
void FIDA_BANDSETJAC(int *flag, int *ier);
void FIDA_SPARSESETJAC(int *ier);

/*** DEPRECATED ***/  
void FIDA_DLSINIT(int *ier);
void FIDA_SPILSINIT(int *ier);
void FIDA_SPILSSETEPSLIN(realtype *eplifac, int *ier);
void FIDA_SPILSSETINCREMENTFACTOR(realtype *dqincfac, int *ier);
void FIDA_SPILSSETJAC(int *flag, int *ier);
void FIDA_SPILSSETPREC(int *flag, int *ier);
/******************/  
  
void FIDA_NLSINIT(int *ier);

void FIDA_SOLVE(realtype *tout, realtype *tret, realtype *yret,
                realtype *ypret, int *itask, int *ier);

void FIDA_FREE(void);
void FIDA_EWTSET(int *flag, int *ier);
void FIDA_GETDKY(realtype *t, int *k, realtype *dky, int *ier);
void FIDA_GETERRWEIGHTS(realtype *eweight, int *ier);
void FIDA_GETESTLOCALERR(realtype *ele, int *ier);

/* Prototypes: Functions Called by the IDA Solver */

int FIDAresfn(realtype t, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data);

int FIDADenseJac(realtype t, realtype c_j, N_Vector yy, N_Vector yp,
                 N_Vector rr, SUNMatrix Jac, void *user_data,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDABandJac(realtype t, realtype c_j, N_Vector yy, N_Vector yp,
                N_Vector rr, SUNMatrix Jac, void *user_data,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDASparseJac(realtype t, realtype c_j, N_Vector y, N_Vector yp,
		  N_Vector rr, SUNMatrix Jac, void *user_data,
		  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

int FIDAJTSetup(realtype t, N_Vector y, N_Vector yp, N_Vector r,
                realtype c_j, void *user_data);

int FIDAJtimes(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
               N_Vector v, N_Vector Jv,
               realtype c_j, void *user_data,
               N_Vector vtemp1, N_Vector vtemp2);

int FIDAPSet(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
             realtype c_j, void *user_data);

int FIDAPSol(realtype t, N_Vector yy, N_Vector yp, N_Vector rr,
             N_Vector rvec, N_Vector zvec,
             realtype c_j, realtype delta, void *user_data);

int FIDAEwtSet(N_Vector yy, N_Vector ewt, void *user_data);

void FIDANullMatrix();
void FIDANullNonlinSol();

/* Declarations for global variables shared amongst various routines */
extern N_Vector F2C_IDA_vec;                 /* defined in FNVECTOR module */
extern N_Vector F2C_IDA_ypvec;               /* defined in fida.c */
extern N_Vector F2C_IDA_ewtvec;              /* defined in fida.c */
extern SUNMatrix F2C_IDA_matrix;             /* defined in FSUNMATRIX module */
extern SUNLinearSolver F2C_IDA_linsol;       /* defined in FSUNLINSOL module */
extern SUNNonlinearSolver F2C_IDA_nonlinsol; /* defined in FSUNNONLINSOL module */
extern void *IDA_idamem;                     /* defined in fida.c */
extern long int *IDA_iout;                   /* defined in fida.c */
extern realtype *IDA_rout;                   /* defined in fida.c */
extern int IDA_nrtfn;                        /* defined in fida.c */

#ifdef __cplusplus
}
#endif

#endif
