/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds and Ting Yan @ SMU
 *    Alan C. Hindmarsh, Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for FCVODE, the Fortran interface to
 * the CVODE package.                                            
 * -----------------------------------------------------------------
 */

/*
 * =============================================================================
 *
 *                FCVODE Interface Package
 *
 * The FCVODE Interface Package is a package of C functions which support
 * the use of the CVODE solver, for the solution of ODE systems 
 * dy/dt = f(t,y), in a mixed Fortran/C setting.  While CVODE is written
 * in C, it is assumed here that the user's calling program and
 * user-supplied problem-defining routines are written in Fortran. 
 * This package provides the necessary interface to CVODE for both the
 * serial and the parallel NVECTOR implementations.
 * 
 * A summary of the user-callable functions, with the corresponding 
 * CVODE functions, are as follows:
 * 
 *   Fortran                    CVODE
 *   ---------------------      --------------------------------
 *   FNVINITS                   N_VNew_Serial
 *   FNVINITP                   N_VNew_Parallel
 *   FNVINITOMP                 N_VNew_OpenMP
 *   FNVINITPTS                 N_VNew_Pthreads
 *
 *   FSUNBANDMATINIT            SUNBandMatrix
 *   FSUNDENSEMATINIT           SUNDenseMatrix
 *   FSUNSPARSEMATINIT          SUNSparseMatrix
 *
 *   FSUNBANDLINSOLINIT         SUNBandLinearSolver
 *   FSUNDENSELINSOLINIT        SUNDenseLinearSolver
 *   FSUNKLUINIT                SUNKLU
 *   FSUNKLUREINIT              SUNKLUReinit
 *   FSUNLAPACKBANDINIT         SUNLapackBand
 *   FSUNLAPACKDENSEINIT        SUNLapackDense
 *   FSUNPCGINIT                SUNPCG
 *   FSUNSPBCGSINIT             SUNSPBCGS
 *   FSUNSPFGMRINIT             SUNSPFGMR
 *   FSUNSPGMRINIT              SUNSPGMR
 *   FSUNSPTFQMRINIT            SUNSPTFQMR
 *   FSUNSUPERLUMTINIT          SUNSuperLUMT
 *
 *   FCVMALLOC                  CVodeCreate, CVodeSetUserData,
 *                                 and CVodeInit
 *   FCVREINIT                  CVReInit
 *   FCVSETIIN                  CVodeSet* (integer arguments)
 *   FCVSETRIN                  CVodeSet* (real arguments)
 *   FCVSETVIN                  CVodeSet* (vector arguments)
 *   FCVEWTSET                  CVodeWFtolerances
 *
 *   FCVLSINIT                  CVodeSetLinearSolver
 *   FCVLSSETEPSLIN             CVodeSetEpsLin
 *   FCVLSSETJAC                CVodeSetJacTimes
 *   FCVLSSETPREC               CVodeSetPreconditioner
 *   FCVDENSESETJAC             CVodeSetJacFn
 *   FCVBANDSETJAC              CVodeSetJacFn
 *   FCVSPARSESETJAC            CVodeSetJacFn
 *
 *   FCVDIAG                    CVDiag
 *
 *   FCVNLSINIT                 CVSetNonlinearSolver
 *
 *   FCVODE                     CVode, CVodeGet*, and CV*Get*
 *   FCVDKY                     CVodeGetDky
 * 
 *   FCVGETERRWEIGHTS           CVodeGetErrWeights
 *   FCVGETESTLOCALERR          CVodeGetEstLocalErrors
 *
 *   FCVFREE                    CVodeFree
 *   ---------------------      --------------------------------
 * 
 * The user-supplied functions, each listed with the corresponding interface
 * function which calls it (and its type within CVODE), are as follows:
 *
 *   Fortran:           Interface Fcn:           CVODE Type:
 *   -------------      ------------------       -----------------------
 *   FCVFUN             FCVf                     CVRhsFn
 *   FCVDJAC            FCVDenseJac              CVLsJacFn
 *   FCVBJAC            FCVBandJac               CVLsJacFn
 *   FCVSPJAC           FCVSparseJac             CVLsJacFn
 *   FCVPSET            FCVPSet                  CVLsPrecSetupFn
 *   FCVPSOL            FCVPSol                  CVLsPrecSolveFn
 *   FCVJTSETUP         FCVJTSetup               CVLsJacTimesSetupFn
 *   FCVJTIMES          FCVJtimes                CVLsJacTimesVecFn
 *   FCVEWT             FCVEwtSet                CVEwtFn
 *   -------------      ------------------       -----------------------
 *
 * In contrast to the case of direct use of CVODE, and of most Fortran ODE
 * solvers, the names of all user-supplied routines here are fixed, in
 * order to maximize portability for the resulting mixed-language program.
 * 
 * Important note on portability.
 * In this package, the names of the interface functions, and the names of
 * the Fortran user routines called by them, appear as dummy names
 * which are mapped to actual values by a series of definitions, in this
 * and other header files.
 * 
 * =============================================================================
 * 
 *                  Usage of the FCVODE Interface Package
 * 
 * The usage of FCVODE requires calls to a variety of interface
 * functions, depending on the method options selected, and one or more
 * user-supplied routines which define the problem to be solved.  These
 * function calls and user routines are summarized separately below.
 * 
 * Some details are omitted, and the user is referred to the user documents
 * on CVODE for more complete documentation.  Information on the
 * arguments of any given user-callable interface routine, or of a given
 * user-supplied function called by an interface function, can be found in
 * the documentation on the corresponding function in the CVODE package.
 * 
 * The number labels on the instructions below end with s for instructions
 * that are specific to use with the serial/OpenMP/PThreads package; similarly 
 * those that end with p are specific to use with the N_VParallel package.
 *
 * -----------------------------------------------------------------------------
 *
 *                               Data Types
 *
 * Throughout this documentation, we will refer to data types according to 
 * their usage in SUNDIALS.  The equivalent types to these may vary, 
 * depending on your computer architecture and on how SUNDIALS was compiled.  
 * A Fortran user should take care that all arguments passed through this 
 * Fortran/C interface are declared of the appropriate type.
 * 
 * Integers: SUNDIALS uses 'int', 'long int' and 'sunindextype' types.  At 
 * compilation, SUNDIALS allows the configuration of the 'index' type, that 
 * accepts values of 32-bit signed and 64-bit signed.  This choice dictates 
 * the size of a SUNDIALS 'sunindextype' variable.
 *   int      -- equivalent to an INTEGER or INTEGER*4 in Fortran
 *   long int -- equivalent to an INTEGER*8 in Fortran (Linux/UNIX/OSX), or 
 *               equivalent to an INTEGER in Windows
 *   sunindextype -- this will depend on the SUNDIALS configuration:
 *               32-bit -- equivalent to an INTEGER or INTEGER*4 in Fortran
 *               64-bit -- equivalent to an INTEGER*8 in Fortran
 *	      
 * Real numbers:  At compilation, SUNDIALS allows the configuration option 
 * '--with-precision', that accepts values of 'single', 'double' or 
 * 'extended' (the default is 'double').  This choice dictates the size of a 
 * SUNDIALS 'realtype' variable.  The corresponding Fortran types for these 
 * 'realtype' sizes are:
 *   single   -- equivalent to a REAL or REAL*4 in Fortran
 *   double   -- equivalent to a DOUBLE PRECISION or REAL*8 in Fortran
 *   extended -- equivalent to a REAL*16 in Fortran
 *
 * -----------------------------------------------------------------------------
 *
 * (1) User-supplied right-hand side routine: FCVFUN
 *
 *   The user must in all cases supply the following Fortran routine
 *
 *       SUBROUTINE FCVFUN (T, Y, YDOT, IPAR, RPAR, IER)
 *
 *   It must set the YDOT array to f(t,y), the right-hand side of the ODE 
 *   system, as function of T = t and the array Y = y.
 *
 *   The arguments are:
 *       Y    -- array containing state variables [realtype, input]
 *       YDOT -- array containing state derivatives [realtype, output]
 *       IPAR -- array containing integer user data that was passed to
 *               FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful, 
 *                 >0 if a recoverable error occurred,
 *                 <0 if an unrecoverable error ocurred.
 * 
 * (2s) Optional user-supplied dense Jacobian approximation routine: FCVDJAC
 *
 *   As an option when using the DENSE or LAPACKDENSE linear solvers, the user may
 *   supply a routine that computes a dense approximation of the system Jacobian 
 *   J = df/dy. If supplied, it must have the following form:
 *
 *       SUBROUTINE FCVDJAC(NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)
 *
 *   Typically this routine will use only NEQ, T, Y, and DJAC. It must compute
 *   the Jacobian and store it columnwise in DJAC.
 *
 *   The arguments are:
 *       NEQ  -- number of rows in the matrix [long int, input]
 *       T    -- current time [realtype, input]
 *       Y    -- array containing state variables [realtype, input]
 *       FY   -- array containing state derivatives [realtype, input]
 *       DJAC -- 2D array containing the jacobian entries [realtype of size
 *               (NEQ,NEQ), output]
 *       H    -- current step size [realtype, input]
 *       IPAR -- array containing integer user data that was passed to
 *               FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       WK*  -- array containing temporary workspace of same size as Y 
 *               [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful, 
 *                 >0 if a recoverable error occurred,
 *                 <0 if an unrecoverable error ocurred.
 * 
 * (2s) Optional user-supplied band Jacobian approximation routine: FCVBJAC
 *
 *   As an option when using the BAND or LAPACKBAND linear solvers, the user 
 *   may supply a routine that computes a band approximation of the system 
 *   Jacobian J = df/dy. If supplied, it must have the following form:
 *
 *       SUBROUTINE FCVBJAC(NEQ, MU, ML, MDIM, T, Y, FY, BJAC, H,
 *      1                   IPAR, RPAR, WK1, WK2, WK3, IER)
 *
 *   Typically this routine will use only NEQ, MU, ML, T, Y, and BJAC. 
 *   It must load the MDIM by N array BJAC with the Jacobian matrix at the
 *   current (t,y) in band form.  Store in BJAC(k,j) the Jacobian element J(i,j)
 *   with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and j = 1 ... N.
 *
 *   The arguments are:
 *       NEQ  -- number of rows in the matrix [long int, input]
 *       MU   -- upper half-bandwidth of the matrix [long int, input]
 *       ML   -- lower half-bandwidth of the matrix [long int, input]
 *       MDIM -- leading dimension of BJAC array [long int, input]
 *       T    -- current time [realtype, input]
 *       Y    -- array containing state variables [realtype, input]
 *       FY   -- array containing state derivatives [realtype, input]
 *       BJAC -- 2D array containing the jacobian entries [realtype of size
 *               (MDIM,NEQ), output]
 *       H    -- current step size [realtype, input]
 *       IPAR -- array containing integer user data that was passed to
 *               FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       WK*  -- array containing temporary workspace of same size as Y 
 *               [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful, 
 *                 >0 if a recoverable error occurred,
 *                 <0 if an unrecoverable error ocurred.
 *
 * 
 * (2s) User-supplied sparse Jacobian approximation routine: FCVSPJAC
 *
 *   Required when using the KLU or SuperLUMT linear solvers, the 
 *   user must supply a routine that computes a compressed-sparse-column [or 
 *   compressed-sparse-row] approximation of the system Jacobian 
 *   J = dfi(t,y)/dy.  If supplied, it must have the following form:
 *
 *       SUBROUTINE FCVSPJAC(T, Y, FY, N, NNZ, JDATA, JRVALS, JCPTRS, 
 *      1                    H, IPAR, RPAR, WK1, WK2, WK3, IER)
 *
 *   This routine must load the N by N compressed sparse column [or row] matrix 
 *   with storage for NNZ nonzeros, stored in the arrays JDATA (nonzero
 *   values), JRVALS (row [or column] indices for each nonzero), JCOLPTRS (indices 
 *   for start of each column [or row]), with the Jacobian matrix at the current
 *   (t,y) in CSC [or CSR] form (see sunmatrix_sparse.h for more information).
 *
 *   The arguments are:
 *         T    -- current time [realtype, input]
 *         Y    -- array containing state variables [realtype, input]
 *         FY   -- array containing state derivatives [realtype, input]
 *         N    -- number of matrix rows/columns in Jacobian [int, input]
 *         NNZ  -- allocated length of nonzero storage [int, input]
 *        JDATA -- nonzero values in Jacobian
 *                 [realtype of length NNZ, output]
 *       JRVALS -- row [or column] indices for each nonzero in Jacobian
 *                 [int of length NNZ, output]
 *       JCPTRS -- pointers to each Jacobian column [or row] in preceding arrays
 *                 [int of length N+1, output]
 *         H    -- current step size [realtype, input]
 *         IPAR -- array containing integer user data that was passed to
 *                 FCVMALLOC [long int, input]
 *         RPAR -- array containing real user data that was passed to
 *                 FCVMALLOC [realtype, input]
 *         WK*  -- array containing temporary workspace of same size as Y 
 *                 [realtype, input]
 *         IER  -- return flag [int, output]:
 *                    0 if successful, 
 *                   >0 if a recoverable error occurred,
 *                   <0 if an unrecoverable error ocurred.
 *
 *   NOTE: this may ONLY be used if SUNDIALS has been configured with 
 *   long int set to 64-bit integers.
 * 
 * (2) Optional user-supplied Jacobian-vector product setup routine: 
 *   FCVJTSETUP
 *
 *   As an option when using the CVLS linear solver interface with a 
 *   matrix-free linear solver, the user may supply a routine that computes 
 *   the product of the system Jacobian J = dfi(t,y)/dy and a given vector v, 
 *   as well as a routine to set up any user data structures in preparation 
 *   for the matrix-vector product.  If a 'setup' routine is supplied, it 
 *   must have the following form:
 *
 *       SUBROUTINE FCVJTSETUP(T, Y, FY, H, IPAR, RPAR, IER)
 *
 *   Typically this routine will use only T and Y.  It must perform any 
 *   relevant preparations for subsequent calls to the user-provided 
 *   FCVJTIMES routine (see below).  
 *
 *   The arguments are:
 *       T    -- current time [realtype, input]
 *       Y    -- array containing state variables [realtype, input]
 *       FY   -- array containing state derivatives [realtype, input]
 *       H    -- current step size [realtype, input]
 *       IPAR -- array containing integer user data that was passed to
 *               FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful, 
 *                  nonzero if an error.
 * 
 * (2) Optional user-supplied Jacobian-vector product routine: FCVJTIMES
 *
 *   As an option when using the SP* linear solver, the user may supply
 *   a routine that computes the product of the system Jacobian J = df/dy and 
 *   a given vector v.  If supplied, it must have the following form:
 *
 *       SUBROUTINE FCVJTIMES (V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)
 *
 *   Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
 *   compute the product vector Jv where the vector v is stored in V, and store
 *   the product in FJV.  
 * 
 *   The arguments are:
 *       V    -- array containing vector to multiply [realtype, input]
 *       JV   -- array containing product vector [realtype, output]
 *       T    -- current time [realtype, input]
 *       Y    -- array containing state variables [realtype, input]
 *       FY   -- array containing state derivatives [realtype, input]
 *       H    -- current step size [realtype, input]
 *       IPAR -- array containing integer user data that was passed to
 *               FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       WORK -- array containing temporary workspace of same size as Y 
 *               [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful, 
 *                  nonzero if an error.
 *
 * (3) Optional user-supplied preconditioner setup/solve routines: FCVPSET 
 *   and FCVPSOL
 *
 *     As an option when using the CVLS linear solver interface and an 
 *     iterative linear solver, the user may supply routines to setup and
 *     apply the preconditioner.  If supplied, these must have the 
 *     following form:
 *
 *       SUBROUTINE FCVPSET(T,Y,FY,JOK,JCUR,GAMMA,H,IPAR,RPAR,IER)
 *
 *     This routine must set up the preconditioner P to be used in the 
 *     subsequent call to FCVPSOL.  The preconditioner (or the product of 
 *     the left and right preconditioners if using both) should be an 
 *     approximation to the matrix  A = I - GAMMA*J  (J = Jacobian),
 *
 *     The arguments are:
 *       T = current time [realtype, input]
 *       Y = current state variable array [realtype, input]
 *       FY = current state variable derivative array [realtype, input]
 *       JOK = flag indicating whether Jacobian-related data needs to be 
 *           recomputed [int, input]:
 *                  0 = recompute, 
 *		  1 = reuse with the current value of GAMMA
 *       JCUR = return flag to denote if Jacobian data was recomputed
 *           [realtype, output], 1=yes, 0=no
 *       GAMMA = Jacobian scaling factor [realtype, input]
 *       H = current time step [realtype, input]
 *       IPAR = array of user integer data [long int, input/output]
 *       RPAR = array with user real data [realtype, input/output]
 *       IER  = return completion flag [int, output]:
 *                  0 = SUCCESS,
 *                 >0 = recoverable failure
 *                 <0 = non-recoverable failure
 *
 *     The user-supplied routine FCVPSOL must have the form:
 *
 *       SUBROUTINE FCVPSOL(T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,IER)
 *
 *     Typically this routine will use only T, Y, GAMMA, R, LR, and Z.  It
 *     must solve the preconditioner linear system Pz = r.  The preconditioner
 *     (or the product of the left and right preconditioners if both are 
 *     nontrivial) should be an approximation to the matrix  I - GAMMA*J  
 *     (J = Jacobian).
 *
 *     The arguments are:
 *       T = current time [realtype, input]
 *       Y = current state variable array [realtype, input]
 *       FY = current state variable derivative array [realtype, input]
 *       R = right-hand side array [realtype, input]
 *       Z = solution array [realtype, output]
 *       GAMMA = Jacobian scaling factor [realtype, input]
 *       DELTA = desired residual tolerance [realtype, input]
 *       LR = flag denoting to solve the right or left preconditioner system
 *                  1 = left preconditioner
 *                  2 = right preconditioner
 *       IPAR = array of user integer data [long int, input/output]
 *       RPAR = array with user real data [realtype, input/output]
 *       IER  = return completion flag [int, output]:
 *                  0 = SUCCESS,
 *                 >0 = recoverable failure
 *                 <0 = non-recoverable failure
 *
 * (4) Optional user-supplied error weight vector routine: FCVEWT
 *
 *   As an option to providing the relative and absolute tolerances, the user
 *   may supply a routine that computes the weights used in the WRMS norms.
 *   If supplied, it must have the following form:
 *
 *       SUBROUTINE FCVEWT (Y, EWT, IPAR, RPAR, IER)
 *
 *   It must store the error weights in EWT, given the current solution vector Y.
 *
 *   The arguments are:
 *       Y    -- array containing state variables [realtype, input]
 *       EWT  -- array containing the error weight vector [realtype, output]
 *       IPAR -- array containing integer user data that was passed to
 *               FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful, 
 *                  nonzero if an error.
 *
 *
 * -----------------------------------------------------------------------------
 *
 * (5) Initialization:  FNVINITS / FNVINITP / FNVINITOMP / FNVINITPTS, 
 *                      FSUNBANDMATINIT / FSUNDENSEMATINIT / 
 *                         FSUNSPARSEMATINIT,
 *                      FSUNBANDLINSOLINIT / FSUNDENSELINSOLINIT / 
 *                         FSUNKLUINIT / FSUNKLUREINIT /
 *                         FSUNKLUSETORDERING / FSUNLAPACKBANDINIT / 
 *                         FSUNLAPACKDENSEINIT / FSUNPCGINIT / 
 *                         FSUNSPBCGSINIT / FSUNSPFGMRINIT / FSUNSPGMRINIT / 
 *                         FSUNSPTFQMRINIT / FSUNSUPERLUMTINIT /
 *                         FSUNSUPERLUMTSETORDERING,
 *                      FCVMALLOC,
 *                      FCVLSINIT,
 *                      FCVREINIT
 *
 *   NOTE: the initialization order is important!  It *must* proceed as 
 *   shown: vector, matrix (if used), linear solver (if used), CVode, 
 *   CVLs, reinit.
 * 
 * (5.1s) To initialize the a vector specification for storing the solution 
 *   data, the user must make one of the following calls:
 *
 *       (serial)   
 *          CALL FNVINITS(1, NEQ, IER)
 *       (MPI parallel)
 *          CALL FNVINITP(COMM, 1, NLOCAL, NGLOBAL, IER)
 *       (OpenMP threaded)
 *          CALL FNVINITOMP(1, NEQ, NUM_THREADS, IER)
 *       (PThreads threaded)
 *          CALL FNVINITPTS(1, NEQ, NUM_THREADS, IER)
 *
 *   In each of these, one argument is an int containing the CVODE solver 
 *   ID (1). 
 *
 *   The other arguments are:
 *        NEQ = size of vectors [long int, input]
 *        COMM = the MPI communicator [int, input]
 *        NLOCAL = local size of vectors on this processor 
 *           [long int, input]
 *        NGLOBAL = the system size, and the global size of vectors (the sum 
 *           of all values of NLOCAL) [long int, input]
 *        NUM_THREADS = number of threads
 *        IER = return completion flag [int, output]:
 *	          0 = success, 
 *		 -1 = failure.
 *
 * (5.2) To initialize a band/dense/sparse matrix structure for 
 *   storing the system Jacobian and for use within a direct linear solver,
 *   the user must make one of the following calls:
 * 
 *          CALL FSUNBANDMATINIT(1, N, MU, ML, SMU, IER)
 *          CALL FSUNDENSEMATINIT(1, M, N, IER)
 *          CALL FSUNSPARSEMATINIT(1, M, N, NNZ, SPARSETYPE, IER)
 *
 *   In each of these, one argument is an int containing the CVODE solver 
 *   ID (1). 
 *
 *   The other arguments are:
 *
 *        M = the number of rows of the matrix [long int, input]
 *        N = the number of columns of the matrix [long int, input]
 *        MU = the number of upper bands (diagonal not included) in a banded 
 *           matrix [long int, input]
 *        ML = the number of lower bands (diagonal not included) in a banded 
 *           matrix [long int, input]
 *        SMU = the number of upper bands to store (diagonal not included) 
 *           for factorization of a banded matrix [long int, input]
 *        NNZ = the storage size (upper bound on the number of nonzeros) for 
 *           a sparse matrix [long int, input]
 *        SPARSETYPE = integer denoting use of CSC (0) vs CSR (1) storage 
 *           for a sparse matrix [int, input]
 *        IER = return completion flag [int, output]:
 *	          0 = success, 
 *		 -1 = failure.
 *
 * (5.3) To initialize a linear solver structure for solving linear systems 
 *   arising from implicit treatment of the IVP, the user must make 
 *   one of the following calls:
 *
 *          CALL FSUNBANDLINSOLINIT(1, IER)
 *          CALL FSUNDENSELINSOLINIT(1, IER)
 *          CALL FSUNKLUINIT(1, IER)
 *          CALL FSUNLAPACKBANDINIT(1, IER)
 *          CALL FSUNLAPACKDENSEINIT(1, IER)
 *          CALL FSUNPCGINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSPBCGSINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSPFGMRINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSPGMRINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSPTFQMRINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSUPERLUMTINIT(1, NUM_THREADS, IER)
 *
 *   Or once these have been initialized, their solver parameters may be
 *   modified via calls to the functions
 *
 *          CALL FSUNKLUSETORDERING(1, ORD_CHOICE, IER)
 *          CALL FSUNSUPERLUMTSETORDERING(1, ORD_CHOICE, IER)
 *
 *          CALL FSUNPCGSETPRECTYPE(1, PRETYPE, IER)
 *          CALL FSUNPCGSETMAXL(1, MAXL, IER)
 *          CALL FSUNSPBCGSSETPRECTYPE(1, PRETYPE, IER)
 *          CALL FSUNSPBCGSSETMAXL(1, MAXL, IER)
 *          CALL FSUNSPFGMRSETGSTYPE(1, GSTYPE, IER)
 *          CALL FSUNSPFGMRSETPRECTYPE(1, PRETYPE, IER)
 *          CALL FSUNSPGMRSETGSTYPE(1, GSTYPE, IER)
 *          CALL FSUNSPGMRSETPRECTYPE(1, PRETYPE, IER)
 *          CALL FSUNSPTFQMRSETPRECTYPE(1, PRETYPE, IER)
 *          CALL FSUNSPTFQMRSETMAXL(1, MAXL, IER)
 *
 *   In all of the above, one argument is an int containing the CVODE solver 
 *   ID (1). 
 *
 *   The other arguments are:
 *
 *        NNZ = the storage size (upper bound on the number of nonzeros) for 
 *           a sparse matrix [long int, input]
 *        ORD_CHOICE = integer denoting ordering choice (see 
 *           SUNKLUSetOrdering and SUNSuperLUMTSetOrdering documentation 
 *           for details) [int, input]
 *        PRETYPE = type of preconditioning to perform (0=none, 1=left, 
 *           2=right, 3=both) [int, input]
 *        MAXL = maximum Krylov subspace dimension [int, input]
 *        GSTYPE = choice of Gram-Schmidt orthogonalization algorithm 
 *           (0=modified, 1=classical) [int, input]
 *        IER = return completion flag [int, output]:
 *	          0 = success, 
 *		 -1 = failure.
 *
 *
 * (5.4) To set various problem and solution parameters and allocate
 *   internal memory, make the following call:
 *
 *       CALL FCVMALLOC(T0, Y0, METH, IATOL, RTOL, ATOL,
 *      1               IOUT, ROUT, IPAR, RPAR, IER)
 *
 *   The arguments are:
 *      T0     = initial value of t [realtype, input]
 *      Y0     = array of initial conditions [realtype, input]
 *      METH   = flag denoting basic integration method [int, input]:
 *                   1 = Adams (nonstiff), 
 *                   2 = BDF (stiff)
 *      IATOL  = flag denoting type for absolute tolerance ATOL [int, input]: 
 *                   1 = scalar, 
 *                   2 = array.
 *                   3 = user-supplied function; the user must supply a routine 
 *                       FCVEWT to compute the error weight vector.
 *      RTOL   = scalar relative tolerance [realtype, input]
 *      ATOL   = scalar or array absolute tolerance [realtype, input]
 *      IOUT   = array of length 21 for integer optional outputs 
 *               [long int, output]
 *      ROUT   = array of length 6 for real optional outputs [realtype, output]
 *      IPAR   = array with user integer data [long int, input/output]
 *      RPAR   = array with user real data [realtype, input/output]
 *	IER    = return completion flag [int, output]:
 *                  0 = SUCCESS,
 *                 -1 = failure (see printed message for failure details).
 *
 *   The user data arrays IPAR and RPAR are passed unmodified to all subsequent
 *   calls to user-provided routines. Modifications to either array inside a
 *   user-provided routine will be propagated. Using these two arrays, the user
 *   can dispense with Common blocks to pass data betwen user-provided routines.
 * 
 *   The optional outputs are:
 *           LENRW   = IOUT( 1) from CVodeGetWorkSpace
 *           LENIW   = IOUT( 2) from CVodeGetWorkSpace
 *           NST     = IOUT( 3) from CVodeGetNumSteps
 *           NFE     = IOUT( 4) from CVodeGetNumRhsEvals
 *           NETF    = IOUT( 5) from CVodeGetNumErrTestFails
 *           NCFN    = IOUT( 6) from CVodeGetNumNonlinSolvConvFails
 *           NNI     = IOUT( 7) from CVodeGetNumNonlinSolvIters
 *           NSETUPS = IOUT( 8) from CVodeGetNumLinSolvSetups
 *           QU      = IOUT( 9) from CVodeGetLastOrder
 *           QCUR    = IOUT(10) from CVodeGetCurrentOrder
 *           NOR     = IOUT(11) from CVodeGetNumStabLimOrderReds
 *           NGE     = IOUT(12) from CVodeGetNumGEvals
 *
 *           H0U     = ROUT( 1) from CVodeGetActualInitStep
 *           HU      = ROUT( 2) from CVodeGetLastStep
 *           HCUR    = ROUT( 3) from CVodeGetCurrentStep
 *           TCUR    = ROUT( 4) from CVodeGetCurrentTime
 *           TOLSF   = ROUT( 5) from CVodeGetTolScaleFactor
 *           UROUND  = ROUT( 6) from UNIT_ROUNDOFF
 *   See the CVODE manual for details. 
 *
 * (5.5) If a linear solver was created in step (5.3) then it must be 
 *   attached to CVode.  If the user called any one of FSUNBANDLINSOLINIT, 
 *   FSUNDENSELINSOLINIT, FSUNKLUINIT, FSUNLAPACKBANDINIT, 
 *   FSUNLAPACKDENSEINIT, FSUNSUPERLUMTINIT, FSUNPCGINIT, 
 *   FSUNSPBCGSINIT, FSUNSPFGMRINIT, FSUNSPGMRINIT, or FSUNSPTFQMRINIT, 
 *   then this must be attached to the CVLS interface using the command:
 *
 *       CALL FCVLSINIT(IER)
 *
 *   The arguments are:
 *	IER  = return completion flag [int, output]:
 *                  0 = SUCCESS,
 *                 -1 = failure (see printed message for failure details).
 *
 * (5.5) If the user instead wishes to use a diagonal approximate Jacobian for 
 *   solving the Newton systems, then it must be created and attached to CVode.  
 *   This choice is appropriate when the Jacobian can be well approximated by
 *   a diagonal matrix.  The user must make the call:
 *       CALL FCVDIAG(IER)
 *
 *   The arguments are:
 *	IER  = return completion flag [int, output]:
 *                  0 = SUCCESS,
 *                 -1 = failure (see printed message for failure details).
 *
 * (5.6) If the user program includes the FCVEWT routine for the evaluation 
 *   of the error weights, the following call must be made
 *
 *       CALL FCVEWTSET(FLAG, IER)
 *
 *   with FLAG = 1 to specify that FCVEWT is provided and 
 *   should be used; FLAG = 0 resets to the default EWT formulation.
 *   The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * (5.7) If the user program includes the FCVBJAC routine for the 
 *   evaluation of the band approximation to the Jacobian, then following 
 *   the call to FCVLSINIT, the following call must be made 
 *
 *       CALL FCVBANDSETJAC(FLAG, IER)
 *
 *   with the int FLAG=1 to specify that FCVBJAC is provided and should be 
 *   used; FLAG=0 specifies a reset to the internal finite difference 
 *   Jacobian approximation.  The int return flag IER=0 if successful, 
 *   nonzero otherwise.
 * 
 *   If the user program includes the FCVDJAC routine for the evaluation 
 *   of the dense approximation to the Jacobian, then after the call to 
 *   FCVLSINIT, the following call must be made 
 *
 *       CALL FCVDENSESETJAC(FLAG, IER)
 *
 *   with the int FLAG=1 to specify that FCVDJAC is provided and should be 
 *   used; FLAG=0 specifies a reset to the internal finite difference 
 *   Jacobian approximation.  The int return flag IER=0 if successful, and 
 *   nonzero otherwise.
 * 
 *   When using a sparse matrix and linear solver the user must provide the
 *   FCVSPJAC routine for the evaluation of the sparse approximation to 
 *   the Jacobian.  To indicate that this routine has been provided, after 
 *   the call to FCVLSINIT, the following call must be made 
 *
 *       CALL FCVSPARSESETJAC(IER)
 *
 *   The int return flag IER=0 if successful, and nonzero otherwise.
 *
 * (5.8) If the user program includes the FCVJTSETUP and FCVJTIMES 
 *   routines for setup of a Jacobian-times-vector product (for use with 
 *   the CVLS interface), then after creating the CVLS interface, 
 *   the following call must be made:
 *
 *       CALL FCVLSSETJAC(FLAG, IER)
 *
 *   with the int FLAG=1 to specify that FCVJTSETUP and FCVJTIMES are 
 *   provided and should be used; FLAG=0 specifies a reset to the internal 
 *   finite difference approximation to this product).  The int return 
 *   flag IER=0 if successful, and nonzero otherwise.
 * 
 * (5.9) If the user program includes the FCVPSET and FCVPSOL routines 
 *   for supplying a preconditioner to an iterative linear solver, then 
 *   after creating the CVLS interface, the following call must be made
 *
 *       CALL FCVLSSETPREC(FLAG, IER)
 *
 *   with the int FLAG=1.  If FLAG=0 then preconditioning with these 
 *   routines will be disabled. The return flag IER=0 if successful, 
 *   nonzero otherwise.
 *
 * (5.10) If the user wishes to use one of CVode's built-in preconditioning 
 *   modules, FCVBP or FCVBBD, then that should be initialized after 
 *   creating the CVLS interface using one of the calls
 *
 *       CALL FCVBPINIT(NEQ, MU, ML, IER)
 *       CALL FCVBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, DQRELY, IER)
 *
 *   Detailed explanation of the inputs to these functions, as well as any 
 *   requirements of user-supplied functions on which these preconditioning 
 *   modules rely, may be found in the header files for each module, 
 *   fcvbp.h or fcvbbd.h, respectively.
 *
 *
 *
 * (5.11) To re-initialize the CVODE solver for the solution of a new problem
 *   of the same size as one already solved, make the following call:
 *
 *       CALL FCVREINIT(T0, Y0, IATOL, RTOL, ATOL, IER)
 *
 *   The arguments have the same names and meanings as those of FCVMALLOC,
 *   except that METH has been omitted from the argument list
 *   (being unchanged for the new problem).  
 *   FCVREINIT performs the same initializations as FCVMALLOC, but does no memory 
 *   allocation, using instead the existing internal memory created by the
 *   previous FCVMALLOC call.  The subsequent calls to 
 *   attach the linear system solver is only needed if that object has been re-created.
 * 
 * (5.12) The SUNKLU solver will reuse much of the factorization information 
 *   from one solve to the next.  If at any time the user wants to force a 
 *   full refactorization or if the number of nonzeros in the Jacobian 
 *   matrix changes, the user should make the call
 *
 *        CALL FSUNKLUREINIT(4, NNZ, REINIT_TYPE, IER)
 *
 *   The arguments are:
 *        NNZ = the maximum number of nonzeros [int; input]
 *        REINIT_TYPE = 1 or 2.  For a value of 1, the matrix will be 
 *          destroyed and a new one will be allocated with NNZ nonzeros.  
 *          For a value of 2, only symbolic and numeric factorizations will 
 *          be completed. 
 * 
 * (5.13) To set various integer optional inputs, make the folowing call:
 *
 *       CALL FCVSETIIN(KEY, VALUE, IER)
 *
 *   to set the integer value VAL to the optional input specified by the
 *   quoted character string KEY.  VALUE must be a Fortran integer of size 
 *   commensurate with a C "long int".
 *   KEY must be one of the following: MAX_ORD, MAX_NSTEPS, MAX_ERRFAIL, 
 *   MAX_NITERS, MAX_CONVFAIL, HNIL_WARNS, STAB_LIM.  The int return flag 
 *   IER is 0 if successful, and <0 otherwise.
 *
 * (5.14) To set various real optional inputs, make the folowing call:
 *
 *       CALL FCVSETRIN(KEY, VALUE, IER)
 *
 *   to set the real value VAL to the optional input specified by the
 *   quoted character string KEY.  VALUE must be a Fortran real-valued 
 *   number of size commensurate with the SUNDIALS "realtype".  KEY must 
 *   be one of the following: INIT_STEP, MAX_STEP, MIN_STEP, STOP_TIME,
 *   NLCONV_COEF.  The int return flag IER is 0 if successful, and <0 otherwise.
 *
 * (5.15) To set the vector of constraints, make the following call:
 *
 *      CALL CVSETVIN(KEY, ARRAY, IER)
 *
 *    where ARRAY is an array of realtype and the quoted character string
 *    KEY is CONSTR_VEC.  The int return flag IER is 0 if successful, and
 *    nonzero otherwise.
 * 
 * -----------------------------------------------------------------------------
 *
 * (6) Optional outputs from CVLS linear solvers (stored in the 
 *   IOUT array that was passed to FCVMALLOC)
 *
 *   Optional outputs specific to the CVLS interface:
 *        LENRWLS  = IOUT(13) from CVodeGetLinWorkSpace (realtype space)
 *        LENIWLS  = IOUT(14) from CVodeGetLinWorkSpace (integer space)
 *        LSTF     = IOUT(15) from CVodeGetLastLinFlag
 *        NFELS    = IOUT(16) from CVodeGetNumLinRhsEvals
 *        NJE      = IOUT(17) from CVodeGetNumJacEvals
 *        NJTS     = IOUT(18) from CVodeGetNumJTSetupEvals
 *        NJTV     = IOUT(19) from CVodeGetNumJtimesEvals
 *        NPE      = IOUT(20) from CVodeGetNumPrecEvals
 *        NPS      = IOUT(21) from CVodeGetNumPrecSolves
 *        NLI      = IOUT(22) from CVodeGetNumLinIters
 *        NCFL     = IOUT(23) from CVodeGetNumLinConvFails
 *
 *   Optional outputs specific to the DIAG case are:
 *        LENRWLS  = IOUT(13) from CVDiagGetWorkSpace
 *        LENIWLS  = IOUT(14) from CVDiagGetWorkSpace
 *        LSTF     = IOUT(15) from CVDiagGetLastFlag
 *        NFELS    = IOUT(16) from CVDiagGetNumRhsEvals
 * 
 *   See the CVODE manual for more detailed descriptions of any of the 
 *   above.
 *
 * -----------------------------------------------------------------------------
 *
 * (7) The integrator: FCVODE
 *
 *   Carrying out the integration is accomplished by making calls as follows:
 *
 *       CALL FCVODE (TOUT, T, Y, ITASK, IER)
 *
 *   The arguments are:
 *      TOUT  = next value of t at which a solution is desired [realtype, input]
 *      T     = value of t reached by the solver on output [realtype, output]
 *      Y     = array containing the computed solution on output [realtype, output]
 *      ITASK = task indicator [int, input]: 
 *              1 = normal mode (overshoot TOUT and interpolate)
 *              2 = one-step mode (return after each internal step taken)
 *              3 = normal tstop mode (like 1, but integration never proceeds past 
 *                  TSTOP, which must be specified through a call to FCVSETRIN
 *                  using the key 'STOP_TIME')
 *              4 = one step tstop (like 2, but integration never goes past TSTOP)
 *      IER   = completion flag [int, output]: 
 *              0 = success, 
 *              1 = tstop return, 
 *              2 = root return, 
 *              values -1 ... -10 are various failure modes (see CVODE manual).
 *   The current values of the optional outputs are available in IOUT and ROUT.
 * 
 * -----------------------------------------------------------------------------
 *
 * (8) Computing solution derivatives: FCVDKY
 *
 *   To obtain a derivative of the solution, of order up to the current method
 *   order, make the following call:
 *
 *       CALL FCVDKY (T, K, DKY, IER)
 *
 *   The arguments are:
 *      T   = value of t at which solution derivative is desired, in
 *             [TCUR-HU,TCUR], [realtype, input].
 *      K   = derivative order (0 .le. K .le. QU) [int, input]
 *      DKY = array containing computed K-th derivative of y [realtype, output]
 *      IER = return flag [int, output]: = 0 for success, <0 for illegal argument
 * 
 * -----------------------------------------------------------------------------
 *
 * (9) Get the current error weight vector: FCVGETERRWEIGHTS
 *
 *   To obtain the current error weight vector, make the following call:
 *
 *       CALL FCVGETERRWEIGHTS(EWT, IER)
 *
 *   The arguments are:
 *       EWT = array containing the error weight vector [realtype, output]
 *       IER = return flag [int, output]: 0=success, nonzero if an error.
 *
 * -----------------------------------------------------------------------------
 *
 * (10) Memory freeing: FCVFREE 
 *
 *   To free the internal memory created by the calls to FCVMALLOC, 
 *   FCVLSINIT and FNVINIT*, make the call
 *
 *       CALL FCVFREE
 * 
 * =============================================================================
 */

#ifndef _FCVODE_H
#define _FCVODE_H

/* header files  */
#include <cvode/cvode.h>
#include <sundials/sundials_linearsolver.h>  /* definition of type SUNLinearSolver */
#include <sundials/sundials_matrix.h>        /* definition of type SUNMatrix */
#include <sundials/sundials_nvector.h>       /* definition of type N_Vector */
#include <sundials/sundials_types.h>         /* definition of type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */

#if defined(SUNDIALS_F77_FUNC)

#define FCV_MALLOC         SUNDIALS_F77_FUNC(fcvmalloc, FCVMALLOC)
#define FCV_REINIT         SUNDIALS_F77_FUNC(fcvreinit, FCVREINIT)
#define FCV_SETIIN         SUNDIALS_F77_FUNC(fcvsetiin, FCVSETIIN)
#define FCV_SETRIN         SUNDIALS_F77_FUNC(fcvsetrin, FCVSETRIN)
#define FCV_SETVIN         SUNDIALS_F77_FUNC(fcvsetvin, FCVSETVIN)
#define FCV_EWTSET         SUNDIALS_F77_FUNC(fcvewtset, FCVEWTSET)
#define FCV_LSINIT         SUNDIALS_F77_FUNC(fcvlsinit, FCVLSINIT)
#define FCV_LSSETJAC       SUNDIALS_F77_FUNC(fcvlssetjac, FCVLSSETJAC)
#define FCV_LSSETPREC      SUNDIALS_F77_FUNC(fcvlssetprec, FCVLSSETPREC)
#define FCV_LSSETEPSLIN    SUNDIALS_F77_FUNC(fcvlssetepslin, FCVLSSETEPSLIN)
#define FCV_DENSESETJAC    SUNDIALS_F77_FUNC(fcvdensesetjac, FCVDENSESETJAC)
#define FCV_BANDSETJAC     SUNDIALS_F77_FUNC(fcvbandsetjac, FCVBANDSETJAC)
#define FCV_SPARSESETJAC   SUNDIALS_F77_FUNC(fcvsparsesetjac, FCVSPARSESETJAC)  
#define FCV_DIAG           SUNDIALS_F77_FUNC(fcvdiag, FCVDIAG)
#define FCV_CVODE          SUNDIALS_F77_FUNC(fcvode, FCVODE)
#define FCV_DKY            SUNDIALS_F77_FUNC(fcvdky, FCVDKY)
#define FCV_FREE           SUNDIALS_F77_FUNC(fcvfree, FCVFREE)
#define FCV_FUN            SUNDIALS_F77_FUNC(fcvfun, FCVFUN)
#define FCV_DJAC           SUNDIALS_F77_FUNC(fcvdjac, FCVDJAC)
#define FCV_BJAC           SUNDIALS_F77_FUNC(fcvbjac, FCVBJAC)
#define FCV_SPJAC          SUNDIALS_F77_FUNC(fcvspjac, FCVSPJAC)
#define FCV_PSOL           SUNDIALS_F77_FUNC(fcvpsol, FCVPSOL)
#define FCV_PSET           SUNDIALS_F77_FUNC(fcvpset, FCVPSET)
#define FCV_JTSETUP        SUNDIALS_F77_FUNC(fcvjtsetup, FCVJTSETUP)
#define FCV_JTIMES         SUNDIALS_F77_FUNC(fcvjtimes, FCVJTIMES)
#define FCV_EWT            SUNDIALS_F77_FUNC(fcvewt, FCVEWT)
#define FCV_GETERRWEIGHTS  SUNDIALS_F77_FUNC(fcvgeterrweights, FCVGETERRWEIGHTS)
#define FCV_GETESTLOCALERR SUNDIALS_F77_FUNC(fcvgetestlocalerr, FCVGETESTLOCALERR)
#define FCV_NLSINIT        SUNDIALS_F77_FUNC(fcvnlsinit, FCVNLSINIT)

/*---DEPRECATED---*/
#define FCV_DLSINIT        SUNDIALS_F77_FUNC(fcvdlsinit, FCVDLSINIT)
#define FCV_DLSSETJAC      SUNDIALS_F77_FUNC(fcvdlssetjac, FCVDLSSETJAC)
#define FCV_SPILSINIT      SUNDIALS_F77_FUNC(fcvspilsinit, FCVSPILSINIT)
#define FCV_SPILSSETPREC   SUNDIALS_F77_FUNC(fcvspilssetprec, FCVSPILSSETPREC)
/*----------------*/

#else

#define FCV_MALLOC         fcvmalloc_
#define FCV_REINIT         fcvreinit_
#define FCV_SETIIN         fcvsetiin_
#define FCV_SETRIN         fcvsetrin_
#define FCV_SETVIN         fcvsetvin_
#define FCV_EWTSET         fcvewtset_
#define FCV_LSINIT         fcvlsinit_
#define FCV_LSSETJAC       fcvlssetjac_
#define FCV_LSSETPREC      fcvlssetprec_
#define FCV_LSSETEPSLIN    fcvlssetepslin_
#define FCV_DENSESETJAC    fcvdensesetjac_
#define FCV_BANDSETJAC     fcvbandsetjac_
#define FCV_SPARSESETJAC   fcvsparsesetjac_
#define FCV_DIAG           fcvdiag_
#define FCV_CVODE          fcvode_
#define FCV_DKY            fcvdky_
#define FCV_FREE           fcvfree_
#define FCV_FUN            fcvfun_
#define FCV_DJAC           fcvdjac_
#define FCV_BJAC           fcvbjac_
#define FCV_SPJAC          fcvspjac_
#define FCV_PSOL           fcvpsol_
#define FCV_PSET           fcvpset_
#define FCV_JTSETUP        fcvjtsetup_
#define FCV_JTIMES         fcvjtimes_
#define FCV_EWT            fcvewt_
#define FCV_GETERRWEIGHTS  fcvgeterrweights_
#define FCV_GETESTLOCALERR fcvgetestlocalerr_
#define FCV_NLSINIT        fcvnlsinit_

/*---DEPRECATED---*/
#define FCV_DLSINIT        fcvdlsinit_
#define FCV_SPILSINIT      fcvspilsinit_
#define FCV_SPILSINIT      fcvspilsinit_
#define FCV_SPILSSETPREC   fcvspilssetprec_
/*----------------*/
  
#endif

  /* Type for user data */

  typedef struct {
    realtype *rpar;
    long int *ipar;
  } *FCVUserData;

  /* Prototypes of exported functions */

  void FCV_MALLOC(realtype *t0, realtype *y0,
                  int *meth, int *iatol,
                  realtype *rtol, realtype *atol,
                  long int *iout, realtype *rout,
                  long int *ipar, realtype *rpar,
                  int *ier);

  void FCV_REINIT(realtype *t0, realtype *y0,
                  int *iatol, realtype *rtol, realtype *atol,
                  int *ier);

  void FCV_SETIIN(char key_name[], long int *ival, int *ier);

  void FCV_SETRIN(char key_name[], realtype *rval, int *ier);

  void FCV_SETVIN(char key_name[], realtype *vval, int *ier);
  void FCV_EWTSET(int *flag, int *ier);

  void FCV_LSINIT(int *ier);
  void FCV_LSSETJAC(int *flag, int *ier);
  void FCV_LSSETPREC(int *flag, int *ier);
  void FCV_LSSETEPSLIN(realtype *eplifac, int *ier);
  void FCV_DENSESETJAC(int *flag, int *ier);
  void FCV_BANDSETJAC(int *flag, int *ier);
  void FCV_SPARSESETJAC(int *ier);

/*---DEPRECATED---*/
  void FCV_DLSINIT(int *ier);
  void FCV_DLSSETJAC(int *flag, int *ier);
  void FCV_SPILSINIT(int *ier);
  void FCV_SPILSSETPREC(int *flag, int *ier);
  void FCV_SPILSSETEPSLIN(realtype *eplifac, int *ier);
/*----------------*/
  
  void FCV_DIAG(int *ier);

  void FCV_NLSINIT(int *ier);

  void FCV_CVODE(realtype *tout, realtype *t, realtype *y, int *itask, int *ier);

  void FCV_DKY(realtype *t, int *k, realtype *dky, int *ier);

  void FCV_GETERRWEIGHTS(realtype *eweight, int *ier);
  void FCV_GETESTLOCALERR(realtype *ele, int *ier);

  void FCV_FREE(void);


  /* Prototypes: Functions Called by the CVODE Solver */
  
  int FCVf(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  
  int FCVDenseJac(realtype t, N_Vector y, N_Vector fy, 
                  SUNMatrix J, void *user_data,
                  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FCVBandJac(realtype t, N_Vector y, N_Vector fy,
                 SUNMatrix J, void *user_data,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FCVSparseJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
		   void *user_data, N_Vector vtemp1,
		   N_Vector vtemp2, N_Vector vtemp3);

  int FCVPSet(realtype tn, N_Vector y, N_Vector fy, booleantype jok,
              booleantype *jcurPtr, realtype gamma, void *user_data);
  
  int FCVPSol(realtype tn, N_Vector y, N_Vector fy, 
              N_Vector r, N_Vector z,
              realtype gamma, realtype delta,
              int lr, void *user_data);
  
  int FCVJTSetup(realtype t, N_Vector y, N_Vector fy, void *user_data);
  
  int FCVJtimes(N_Vector v, N_Vector Jv, realtype t, 
                N_Vector y, N_Vector fy,
                void *user_data, N_Vector work);
  
  int FCVEwtSet(N_Vector y, N_Vector ewt, void *user_data);

  void FCVNullMatrix();
  void FCVNullLinsol();
  void FCVNullNonlinSol();

  /* Declarations for global variables shared amongst various routines */

  extern N_Vector F2C_CVODE_vec;                 /* defined in FNVECTOR module      */
  extern SUNMatrix F2C_CVODE_matrix;             /* defined in FSUNMATRIX module    */
  extern SUNLinearSolver F2C_CVODE_linsol;       /* defined in FSUNLINSOL module    */
  extern SUNNonlinearSolver F2C_CVODE_nonlinsol; /* defined in FSUNNONLINSOL module */

  extern void *CV_cvodemem;        /* defined in fcvode.c */
  extern long int *CV_iout;        /* defined in fcvode.c */
  extern realtype *CV_rout;        /* defined in fcvode.c */
  extern int CV_nrtfn;             /* defined in fcvode.c */
  extern int CV_ls;                /* defined in fcvode.c */

  /* Linear solver IDs */

  enum { CV_LS_STD = 0, CV_LS_DIAG = 1 };

#ifdef __cplusplus
}
#endif

#endif
