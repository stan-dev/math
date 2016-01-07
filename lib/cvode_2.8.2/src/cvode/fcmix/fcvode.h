/*
 * -----------------------------------------------------------------
 * $Revision: 4423 $
 * $Date: 2015-03-08 17:23:10 -0700 (Sun, 08 Mar 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, Radu Serban and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
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
 * The user-callable functions, with the corresponding CVODE functions,
 * are as follows:
 * 
 *   FNVINITS and FNVINITP interface to N_VNew_Serial and
 *               N_VNew_Parallel, respectively
 *   FNVINITOMP               N_VNew_OpenMP
 *   FNVINITPTS               N_VNew_Pthreads
 * 
 *   FCVMALLOC  interfaces to CVodeCreate, CVodeSetUserData, and CVodeInit
 * 
 *   FCVREINIT  interfaces to CVReInit
 * 
 *   FCVSETIIN and FCVSETRIN interface to CVodeSet*
 *
 *   FCVEWTSET  interfaces to CVodeWFtolerances
 * 
 *   FCVDIAG    interfaces to CVDiag
 * 
 *   FCVDENSE   interfaces to CVDense
 *   FCVDENSESETJAC   interfaces to CVDenseSetJacFn
 * 
 *   FCVBAND    interfaces to CVBand
 *   FCVBANDSETJAC    interfaces to CVBandSetJacFn
 *
 *   FCVLAPACKDENSE   interfaces to CVLapackDense
 *   FCVLAPACKBAND    interfaces to CVLapackBand
 *   FCVLAPACKDENSESETJAC  interfaces to CVLapackSetJacFn
 *   FCVLAPACKBANDSETJAC   interfaces to CVLapackSetJacFn
 *
 *   FCVKLU         interfaces to CVKLU
 *   FCVKLUReinit   interfaces to CVKLUReinit
 *   FCVSUPERLUMT   interfaces to CVSuperLUMT
 *
 *   FCVSPGMR and FCVSPGMRREINIT interface to CVSpgmr and CVSpilsSet*
 *   FCVSPBCG, FCVSPBCGREINIT interface to CVSpbcg and CVSpilsSet*
 *   FCVSPTFQMR, FCVSPTFQMRREINIT interface to CVSptfqmr and CVSpilsSet*
 *
 *   FCVSPILSSETJAC   interfaces to CVSpilsSetJacTimesVecFn
 *   FCVSPILSSETPREC  interfaces to CVSpilsSetPreconditioner
 * 
 *   FCVODE     interfaces to CVode, CVodeGet*, and CV*Get*
 * 
 *   FCVDKY     interfaces to CVodeGetDky
 * 
 *   FCVGETERRWEIGHTS  interfaces to CVodeGetErrWeights
 *
 *   FCVGETESTLOCALERR  interfaces to CVodeGetEstLocalErrors
 *
 *   FCVFREE    interfaces to CVodeFree
 * 
 * The user-supplied functions, each listed with the corresponding interface
 * function which calls it (and its type within CVODE), are as follows:
 *   FCVFUN    is called by the interface function FCVf of type CVRhsFn
 *   FCVDJAC   is called by the interface fn. FCVDenseJac of type CVDenseJacFn
 *   FCVBJAC   is called by the interface fn. FCVBandJac of type CVBandJacFn
 *   FCVLDJAC  is called by the interface fn. FCVLapackDenseJac of type CVLapackJacFn
 *   FCVLBJAC  is called by the interface fn. FCVLapackBandJac of type CVLapackJacFn
 *   FCVPSOL   is called by the interface fn. FCVPSol of type CVSpilsPrecSolveFn
 *   FCVPSET   is called by the interface fn. FCVPSet of type CVSpilsPrecSetupFn
 *   FCVJTIMES is called by interface fn. FCVJtimes of type CVSpilsJacTimesVecFn
 *   FCVSPJAC  is called by interface fn. FCVSparseJac of type CVSlsSparseJacFn
 *   FCVEWT    is called by interface fn. FCVEwtSet of type CVEwtFn
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
 * The usage of FCVODE requires calls to five or more interface
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
 * that apply to the serial version of CVODE only, and end with p for
 * those that apply to the parallel version only.
 *
 * -----------------------------------------------------------------------------
 *
 * (1) User-supplied right-hand side routine: FCVFUN
 * The user must in all cases supply the following Fortran routine
 *       SUBROUTINE FCVFUN (T, Y, YDOT, IPAR, RPAR, IER)
 *       DIMENSION Y(*), YDOT(*), IPAR(*), RPAR(*)
 * It must set the YDOT array to f(t,y), the right-hand side of the ODE 
 * system, as function of T = t and the array Y = y.  Here Y and YDOT
 * are distributed vectors. IPAR and RPAR are arrays of integer and real user 
 * data, respectively as passed to FCVMALLOC.
 * On return, set IER = 0 if successful, IER > 0 if a recoverable error occurred,
 * and IER < 0 if an unrecoverable error ocurred.
 * 
 * (2s) Optional user-supplied dense Jacobian approximation routine: FCVDJAC
 * As an option when using the DENSE linear solver, the user may supply a
 * routine that computes a dense approximation of the system Jacobian 
 * J = df/dy. If supplied, it must have the following form:
 *       SUBROUTINE FCVDJAC (NEQ, T, Y, FY, DJAC, H, IPAR, RPAR, WK1, WK2, WK3, IER)
 *       DIMENSION Y(*), FY(*), DJAC(NEQ,*), IPAR(*), RPAR(*), WK1(*), WK2(*), WK3(*)
 * Typically this routine will use only NEQ, T, Y, and DJAC. It must compute
 * the Jacobian and store it columnwise in DJAC.
 * IPAR and RPAR are user (integer and real) arrays passed to FCVMALLOC.
 * On return, set IER = 0 if successful, IER > 0 if a recoverable error occurred,
 * and IER < 0 if an unrecoverable error ocurred.
 * 
 * (3s) Optional user-supplied band Jacobian approximation routine: FCVBJAC
 * As an option when using the BAND linear solver, the user may supply a
 * routine that computes a band approximation of the system Jacobian 
 * J = df/dy. If supplied, it must have the following form:
 *       SUBROUTINE FCVBJAC (NEQ, MU, ML, MDIM, T, Y, FY, BJAC, H,
 *      1                    IPAR, RPAR, WK1, WK2, WK3, IER)
 *       DIMENSION Y(*), FY(*), BJAC(MDIM,*), IPAR(*), RPAR(*), WK1(*), WK2(*), WK3(*)
 * Typically this routine will use only NEQ, MU, ML, T, Y, and BJAC. 
 * It must load the MDIM by N array BJAC with the Jacobian matrix at the
 * current (t,y) in band form.  Store in BJAC(k,j) the Jacobian element J(i,j)
 * with k = i - j + MU + 1 (k = 1 ... ML+MU+1) and j = 1 ... N.
 * IPAR and RPAR are user (integer and real) arrays passed to FCVMALLOC.
 * On return, set IER = 0 if successful, IER > 0 if a recoverable error occurred,
 * and IER < 0 if an unrecoverable error ocurred.
 * 
 (4s) User-supplied sparse Jacobian approximation routine: FCVSPJAC

     Required when using the CVKLU or CVSuperLUMT linear solvers, the 
     user must supply a routine that computes a compressed-sparse-column 
     approximation of the system Jacobian J = dfi(t,y)/dy.  If supplied, 
     it must have the following form:

       SUBROUTINE FCVSPJAC(T, Y, FY, N, NNZ, JDATA, JRVALS, 
      &                     JCPTRS, H, IPAR, RPAR, WK1, WK2, WK3, IER)

     This routine must load the N by N compressed sparse column matrix 
     with storage for NNZ nonzeros, stored in the arrays JDATA (nonzero
     values), JRVALS (row indices for each nonzero), JCOLPTRS (indices 
     for start of each column), with the Jacobian matrix at the current
     (t,y) in CSC form (see sundials_sparse.h for more information).

     The arguments are:
         T    -- current time [realtype, input]
         Y    -- array containing state variables [realtype, input]
         FY   -- array containing state derivatives [realtype, input]
         N    -- number of matrix rows/columns in Jacobian [int, input]
         NNZ  -- allocated length of nonzero storage [int, input]
        JDATA -- nonzero values in Jacobian
                 [realtype of length NNZ, output]
       JRVALS -- row indices for each nonzero in Jacobian
                  [int of length NNZ, output]
       JCPTRS -- pointers to each Jacobian column in preceding arrays
                 [int of length N+1, output]
         H    -- current step size [realtype, input]
         IPAR -- array containing integer user data that was passed to
                 FCVMALLOC [long int, input]
         RPAR -- array containing real user data that was passed to
                 FCVMALLOC [realtype, input]
         WK*  -- array containing temporary workspace of same size as Y 
                 [realtype, input]
         IER  -- return flag [int, output]:
                    0 if successful, 
                   >0 if a recoverable error occurred,
                   <0 if an unrecoverable error ocurred.

 * (5s) Optional user-supplied Lapack dense Jacobian routine: FCVLDJAC
 * See the description for FCVDJAC. NOTE: the dense Jacobian matrix
 * is NOT set to zero before calling the user's FCVLDJAC.
 *
 * (6s) Optional user-supplied Lapack band Jacobian routine: FCVLBJAC
 * See the description for FCVBJAC. NOTE: the band Jacobian matrix
 * is NOT set to zero before calling the user's FCVLBJAC.
 *
 * (7) Optional user-supplied Jacobian-vector product routine: FCVJTIMES
 * As an option when using the SP* linear solver, the user may supply
 * a routine that computes the product of the system Jacobian J = df/dy and 
 * a given vector v.  If supplied, it must have the following form:
 *       SUBROUTINE FCVJTIMES (V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)
 *       DIMENSION V(*), FJV(*), Y(*), FY(*), IPAR(*), RPAR(*), WORK(*)
 * Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
 * compute the product vector Jv where the vector v is stored in V, and store
 * the product in FJV.  On return, set IER = 0 if FCVJTIMES was successful,
 * and nonzero otherwise.
 * IPAR and RPAR are user (integer and real) arrays passed to FCVMALLOC.
 * 
 * (8) Optional user-supplied error weight vector routine: FCVEWT
 * As an option to providing the relative and absolute tolerances, the user
 * may supply a routine that computes the weights used in the WRMS norms.
 * If supplied, it must have the following form:
 *       SUBROUTINE FCVEWT (Y, EWT, IPAR, RPAR, IER)
 *       DIMENSION Y(*), EWT(*), IPAR(*), RPAR(*)
 * It must store the error weights in EWT, given the current solution vector Y.
 * On return, set IER = 0 if successful, and nonzero otherwise.
 * IPAR and RPAR are user (integer and real) arrays passed to FCVMALLOC.
 *
 * -----------------------------------------------------------------------------
 *
 * (9) Initialization:  FNVINITS/FNVINITP/FNVINITOMP/FNVINITPTS, 
 *                      FCVMALLOC, FCVREINIT
 * 
 * (9.1s) To initialize the serial machine environment, the user must make
 * the following call:
 *        CALL FNVINITS (1, NEQ, IER)
 * where the first argument is the CVODE solver ID. The other arguments are:
 * NEQ     = size of vectors
 * IER     = return completion flag. Values are 0 = success, -1 = failure.
 * 
 * (9.1p) To initialize the distributed memory parallel machine environment, 
 * the user must make the following call:
 *        CALL FNVINITP (1, NLOCAL, NGLOBAL, IER)
 * The arguments are:
 * NLOCAL  = local size of vectors on this processor
 * NGLOBAL = the system size, and the global size of vectors (the sum 
 *           of all values of NLOCAL)
 * IER     = return completion flag. Values are 0 = success, -1 = failure.
 * Note: If MPI was initialized by the user, the communicator must be
 * set to MPI_COMM_WORLD.  If not, this routine initializes MPI and sets
 * the communicator equal to MPI_COMM_WORLD.
 * 
 (9.1omp) To initialize the openMP threaded vector kernel, 
          the user must make the following call:

          CALL FNVINITOMP (1, NEQ, NUM_THREADS, IER)

        The arguments are:
          NEQ = size of vectors
          NUM_THREADS = number of threads
          IER = return completion flag. Values are 0 = success, -1 = failure.

 (9.1pts) To initialize the Pthreads threaded vector kernel, 
          the user must make the following call:

          CALL FNVINITOMP (1, NEQ, NUM_THREADS, IER)

        The arguments are:
          NEQ = size of vectors
          NUM_THREADS = number of threads
          IER = return completion flag. Values are 0 = success, -1 = failure.

 * (9.2) To set various problem and solution parameters and allocate
 * internal memory, make the following call:
 *       CALL FCVMALLOC(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL,
 *      1               IOUT, ROUT, IPAR, RPAR, IER)
 * The arguments are:
 * T0     = initial value of t
 * Y0     = array of initial conditions
 * METH   = basic integration method: 1 = Adams (nonstiff), 2 = BDF (stiff)
 * ITMETH = nonlinear iteration method: 1=functional iteration, 2=Newton iter.
 * IATOL  = type for absolute tolerance ATOL: 1 = scalar, 2 = array.
 *          If IATOL = 3, then the user must supply a routine FCVEWT to compute
 *          the error weight vector.
 * RTOL   = relative tolerance (scalar)
 * ATOL   = absolute tolerance (scalar or array)
 * IOUT   = array of length 21 for integer optional outputs
 *          (declare as INTEGER*4 or INTEGER*8 according to C type long int)
 * ROUT   = array of length 6 for real optional outputs
 * IPAR   = array with user integer data
 *          (declare as INTEGER*4 or INTEGER*8 according to C type long int)
 * RPAR   = array with user real data
 * IER    = return completion flag.  Values are 0 = SUCCESS, and -1 = failure.
 *          See printed message for details in case of failure.
 *
 * The user data arrays IPAR and RPAR are passed unmodified to all subsequent
 * calls to user-provided routines. Modifications to either array inside a
 * user-provided routine will be propagated. Using these two arrays, the user
 * can dispense with Common blocks to pass data betwen user-provided routines.
 * 
 * The optional outputs are:
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
 * See the CVODE manual for details. 
 *
 * If the user program includes the FCVEWT routine for the evaluation of the 
 * error weights, the following call must be made
 *       CALL FCVEWTSET(FLAG, IER)
 * with FLAG = 1 to specify that FCVEWT is provided.
 * The return flag IER is 0 if successful, and nonzero otherwise.
 *
 * (9.3) To re-initialize the CVODE solver for the solution of a new problem
 * of the same size as one already solved, make the following call:
 *       CALL FCVREINIT(T0, Y0, IATOL, RTOL, ATOL, IER)
 * The arguments have the same names and meanings as those of FCVMALLOC,
 * except that METH and ITMETH  have been omitted from the argument list 
 * (being unchanged for the new problem).  
 * FCVREINIT performs the same initializations as FCVMALLOC, but does no memory 
 * allocation, using instead the existing internal memory created by the
 * previous  FCVMALLOC call.  The call to specify the linear system solution
 * method may or may not be needed; see paragraph (7) below.
 * 
 * (9.4) To set various integer optional inputs, make the folowing call:
 *       CALL FCVSETIIN(KEY, VALUE, IER)
 * to set the integer value VAL to the optional input specified by the
 * quoted character string KEY.
 * KEY is one of the following: MAX_ORD, MAX_NSTEPS, MAX_ERRFAIL, MAX_NITERS, 
 * MAX_CONVFAIL, HNIL_WARNS, STAB_LIM.
 *
 * To set various real optional inputs, make the folowing call:
 *       CALL FCVSETRIN(KEY, VALUE, IER)
 * to set the real value VAL to the optional input specified by the
 * quoted character string KEY.
 * KEY is one of the following: INIT_STEP, MAX_STEP, MIN_STEP, STOP_TIME,
 * NLCONV_COEF.
 *
 * FCVSETIIN and FCVSETRIN return IER = 0 if successful and IER < 0 if an 
 * error occured.
 *
 * -----------------------------------------------------------------------------
 *
 * (10) Specification of linear system solution method.
 * In the case of a stiff system, the implicit BDF method involves the solution
 * of linear systems related to the Jacobian J = df/dy of the ODE system.
 * CVODE presently includes four choices for the treatment of these systems,
 * and the user of FCVODE must call a routine with a specific name to make the
 * desired choice.
 * 
 * (10.1) Diagonal approximate Jacobian.
 * This choice is appropriate when the Jacobian can be well approximated by
 * a diagonal matrix.  The user must make the call:
 *       CALL FCVDIAG(IER)
 * IER is an error return flag: 0 = success, negative value = error.
 * There is no additional user-supplied routine.  
 *
 * Optional outputs specific to the DIAG case are:
 *        LENRWLS = IOUT(13) from CVDiagGetWorkSpace
 *        LENIWLS = IOUT(14) from CVDiagGetWorkSpace
 *        LSTF    = IOUT(15) from CVDiagGetLastFlag
 *        NFELS   = IOUT(16) from CVDiagGetNumRhsEvals
 * See the CVODE manual for descriptions.
 * 
 * (10.2s) DENSE treatment of the linear system.
 * The user must make the call
 *       CALL FCVDENSE(NEQ, IER)
 * The argument is:
 * IER = error return flag: 0 = success , negative value = an error occured
 * 
 * If the user program includes the FCVDJAC routine for the evaluation of the 
 * dense approximation to the Jacobian, the following call must be made
 *       CALL FCVDENSESETJAC(FLAG, IER)
 * with FLAG = 1 to specify that FCVDJAC is provided.  (FLAG = 0 specifies
 * using the internal finite differences approximation to the Jacobian.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 * 
 * Optional outputs specific to the DENSE case are:
 *        LENRWLS = IOUT(13) from CVDenseGetWorkSpace
 *        LENIWLS = IOUT(14) from CVDenseGetWorkSpace
 *        LSTF    = IOUT(15) from CVDenseGetLastFlag
 *        NFELS   = IOUT(16) from CVDenseGetNumRhsEvals
 *        NJED    = IOUT(17) from CVDenseGetNumJacEvals
 * See the CVODE manual for descriptions.
 * 
 * (10.3s) BAND treatment of the linear system
 * The user must make the call
 *       CALL FCVBAND(NEQ, MU, ML, IER)
 * The arguments are:
 * MU  = upper bandwidth
 * ML  = lower bandwidth
 * IER = error return flag: 0 = success , negative value = an error occured
 * 
 * If the user program includes the FCVBJAC routine for the evaluation of the 
 * band approximation to the Jacobian, the following call must be made
 *       CALL FCVBANDSETJAC(FLAG, IER)
 * with FLAG = 1 to specify that FCVBJAC is provided.  (FLAG = 0 specifies
 * using the internal finite differences approximation to the Jacobian.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 * 
 * Optional outputs specific to the BAND case are:
 *        LENRWLS = IOUT(13) from CVBandGetWorkSpace
 *        LENIWLS = IOUT(14) from CVBandGetWorkSpace
 *        LSTF    = IOUT(15) from CVBandGetLastFlag
 *        NFELS   = IOUT(16) from CVBandGetNumRhsEvals
 *        NJEB    = IOUT(17) from CVBandGetNumJacEvals
 * See the CVODE manual for descriptions.
 *
 * (10.4s) LAPACK dense treatment of the linear system
 * The user must make the call
 *       CALL FCVLAPACKDENSE(NEQ, IER)
 * and, optionally
 *       CALL FCVLAPACKDENSESETJAC(FLAG, IER)
 * with FLAG=1 if the user provides the function FCVLDJAC. 
 * See (9.2s) for more details.
 *
 * (10.5s) LAPACK band treatment of the linear system
 * The user must make the call
 *       CALL FCVLAPACKBAND(NEQ, IER)
 * and, optionally
 *       CALL FCVLAPACKBANDSETJAC(FLAG, IER)
 * with FLAG=1 if the user provides the function FCVLBJAC. 
 * See (9.3s)
 *
  (10.6s) SPARSE treatment of the linear system using the KLU solver.

     The user must make the call

       CALL FCVKLU(NEQ, NNZ, ORDERING, IER)

     The arguments are:
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	ORDERING = the matrix ordering desired, possible values
	   come from the KLU package (0 = AMD, 1 = COLAMD) [int; input]
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
 
     The CVODE KLU solver will reuse much of the factorization information from one
     nonlinear iteration to the next.  If at any time the user wants to force a full
     refactorization or if the number of nonzeros in the Jacobian matrix changes, the
     user should make the call

       CALL FCVKLUREINIT(NEQ, NNZ, REINIT_TYPE)

     The arguments are:
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	REINIT_TYPE = 1 or 2.  For a value of 1, the matrix will be destroyed and 
          a new one will be allocated with NNZ nonzeros.  For a value of 2, 
	  only symbolic and numeric factorizations will be completed. 
 
     When using FCVKLU, the user is required to supply the FCVSPJAC 
     routine for the evaluation of the sparse approximation to the 
     Jacobian, as discussed above with the other user-supplied routines.
 
     Optional outputs specific to the KLU case are:
        LSTF    = IOUT(14) from CVSlsGetLastFlag
        NJES    = IOUT(16) from CVSlsGetNumJacEvals
     See the CVODE manual for descriptions.
 
 (10.7s) SPARSE treatment of the linear system using the SuperLUMT solver.

     The user must make the call

       CALL FCVSUPERLUMT(NTHREADS, NEQ, NNZ, ORDERING, IER)

     The arguments are:
        NTHREADS = desired number of threads to use [int; input]
        NEQ = the problem size [int; input]
        NNZ = the maximum number of nonzeros [int; input]
	ORDERING = the matrix ordering desired, possible values
	   come from the SuperLU_MT package [int; input]
           0 = Natural
           1 = Minimum degree on A^T A
           2 = Minimum degree on A^T + A
           3 = COLAMD
	IER = error return flag [int, output]: 
	         0 = success, 
		 negative = error.
	 
     At this time, there is no reinitialization capability for the SUNDIALS 
     interfaces to the SuperLUMT solver.

     When using FCVSUPERLUMT, the user is required to supply the FCVSPJAC 
     routine for the evaluation of the sparse approximation to the 
     Jacobian, as discussed above with the other user-supplied routines.
 
     Optional outputs specific to the SUPERLUMT case are:
        LSTF    = IOUT(14) from CVSlsGetLastFlag
        NJES    = IOUT(16) from CVSlsGetNumJacEvals
     See the CVODE manual for descriptions.
 
 * (10.8) SPGMR treatment of the linear systems.
 * For the Scaled Preconditioned GMRES solution of the linear systems,
 * the user must make the following call:
 *       CALL FCVSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)              
 * The arguments are:
 * IPRETYPE = preconditioner type: 
 *              0 = none 
 *              1 = left only
 *              2 = right only
 *              3 = both sides
 * IGSTYPE  = Gram-schmidt process type: 
 *              1 = modified G-S
 *              2 = classical G-S.
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * DELT     = linear convergence tolerance factor; 0.0 indicates default.
 * IER      = error return flag: 0 = success; negative value = an error occured
 * 
 * 
 * Optional outputs specific to the SPGMR case are:
 *        LENRWLS = IOUT(13) from CVSpgmrGetWorkSpace
 *        LENIWLS = IOUT(14) from CVSpgmrGetWorkSpace
 *        LSTF    = IOUT(15) from CVSpgmrGetLastFlag
 *        NFELS   = IOUT(16) from CVSpgmrGetRhsEvals
 *        NJTV    = IOUT(17) from CVSpgmrGetJtimesEvals
 *        NPE     = IOUT(18) from CVSpgmrGetPrecEvals
 *        NPS     = IOUT(19) from CVSpgmrGetPrecSolves
 *        NLI     = IOUT(20) from CVSpgmrGetLinIters
 *        NCFL    = IOUT(21) from CVSpgmrGetConvFails
 * See the CVODE manual for descriptions.
 * 
 * If a sequence of problems of the same size is being solved using the
 * SPGMR linear solver, then following the call to FCVREINIT, a call to the
 * FCVSPGMRREINIT routine is needed if any of IPRETYPE, IGSTYPE, DELT is
 * being changed.  In that case, call FCVSPGMRREINIT as follows:
 *       CALL FCVSPGMRREINIT(IPRETYPE, IGSTYPE, DELT, IER)              
 * The arguments have the same meanings as for FCVSPGMR.  If MAXL is being
 * changed, then call FCVSPGMR instead.
 * 
 * (10.9) SPBCG treatment of the linear systems.
 * For the Scaled Preconditioned Bi-CGSTAB solution of the linear systems,
 * the user must make the following call:
 *       CALL FCVSPBCG(IPRETYPE, MAXL, DELT, IER)              
 * The arguments are:
 * IPRETYPE = preconditioner type: 
 *              0 = none 
 *              1 = left only
 *              2 = right only
 *              3 = both sides
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * DELT     = linear convergence tolerance factor; 0.0 indicates default.
 * IER      = error return flag: 0 = success; negative value = an error occured
 * 
 * Optional outputs specific to the SPBCG case are:
 *        LENRWLS = IOUT(13) from CVSpbcgGetWorkSpace
 *        LENIWLS = IOUT(14) from CVSpbcgGetWorkSpace
 *        LSTF    = IOUT(15) from CVSpbcgGetLastFlag
 *        NFELS   = IOUT(16) from CVSpbcgGetRhsEvals
 *        NJTV    = IOUT(17) from CVSpbcgGetJtimesEvals
 *        NPE     = IOUT(18) from CVSpbcgGetPrecEvals
 *        NPS     = IOUT(19) from CVSpbcgGetPrecSolves
 *        NLI     = IOUT(20) from CVSpbcgGetLinIters
 *        NCFL    = IOUT(21) from CVSpbcgGetConvFails
  * See the CVODE manual for descriptions.
 * 
 * If a sequence of problems of the same size is being solved using the
 * SPBCG linear solver, then following the call to FCVREINIT, a call to the
 * FCVSPBCGREINIT routine is needed if any of its arguments is
 * being changed.  The call is:
 *       CALL FCVSPBCGREINIT(IPRETYPE, MAXL, DELT, IER)              
 * The arguments have the same meanings as for FCVSPBCG.
 *
 * (10.10) SPTFQMR treatment of the linear systems.
 * For the Scaled Preconditioned TFQMR solution of the linear systems,
 * the user must make the following call:
 *       CALL FCVSPTFQMR(IPRETYPE, MAXL, DELT, IER)              
 * The arguments are:
 * IPRETYPE = preconditioner type: 
 *              0 = none 
 *              1 = left only
 *              2 = right only
 *              3 = both sides
 * MAXL     = maximum Krylov subspace dimension; 0 indicates default.
 * DELT     = linear convergence tolerance factor; 0.0 indicates default.
 * IER      = error return flag: 0 = success; negative value = an error occured
 * 
 * Optional outputs specific to the SPTFQMR case are:
 *        LENRWLS = IOUT(13) from CVSptfqmrGetWorkSpace
 *        LENIWLS = IOUT(14) from CVSptfqmrGetWorkSpace
 *        LSTF    = IOUT(15) from CVSptfqmrGetLastFlag
 *        NFELS   = IOUT(16) from CVSptfqmrGetRhsEvals
 *        NJTV    = IOUT(17) from CVSptfqmrGetJtimesEvals
 *        NPE     = IOUT(18) from CVSptfqmrGetPrecEvals
 *        NPS     = IOUT(19) from CVSptfqmrGetPrecSolves
 *        NLI     = IOUT(20) from CVSptfqmrGetLinIters
 *        NCFL    = IOUT(21) from CVSptfqmrGetConvFails
 * See the CVODE manual for descriptions.
 *
 * If a sequence of problems of the same size is being solved using the
 * SPTFQMR linear solver, then following the call to FCVREINIT, a call to the
 * FCVSPTFQMRREINIT routine is needed if any of its arguments is
 * being changed.  The call is:
 *       CALL FCVSPTFQMRREINIT(IPRETYPE, MAXL, DELT, IER)              
 * The arguments have the same meanings as for FCVSPTFQMR.
 *
 * (10.11) Usage of user-supplied routines for the Krylov solvers
 *
 * If the user program includes the FCVJTIMES routine for the evaluation of the 
 * Jacobian vector product, the following call must be made
 *       CALL FCVSPILSSETJAC(FLAG, IER)
 * with FLAG = 1 to specify that FCVJTIMES is provided.  (FLAG = 0 specifies
 * using and internal finite difference approximation to this product.)
 * The return flag IER is 0 if successful, and nonzero otherwise.
 * 
 * Usage of the user-supplied routines FCVPSOL and FCVPSET for solution of the 
 * preconditioner linear system requires the following call:
 *       CALL FCVSPILSSETPREC(FLAG, IER)
 * with FLAG = 1. The return flag IER is 0 if successful, nonzero otherwise.
 * The user-supplied routine FCVPSOL must have the form:
 *       SUBROUTINE FCVPSOL (T,Y,FY,R,Z,GAMMA,DELTA,LR,IPAR,RPAR,VT,IER)
 *       DIMENSION Y(*), FY(*), VT(*), R(*), Z(*), IPAR(*), RPAR(*)
 * Typically this routine will use only NEQ, T, Y, GAMMA, R, LR, and Z.  It
 * must solve the preconditioner linear system Pz = r, where r = R is input, 
 * and store the solution z in Z.  Here P is the left preconditioner if LR = 1
 * and the right preconditioner if LR = 2.  The preconditioner (or the product
 * of the left and right preconditioners if both are nontrivial) should be an 
 * approximation to the matrix I - GAMMA*J (I = identity, J = Jacobian).
 * IPAR and RPAR are user (integer and real) arrays passed to FCVMALLOC.
 * On return, set IER = 0 if successful, IER > 0 if a recoverable error occurred,
 * and IER < 0 if an unrecoverable error ocurred.
 *
 * -----------------------------------------------------------------------------
 *
 * (11) The integrator: FCVODE
 * Carrying out the integration is accomplished by making calls as follows:
 *       CALL FCVODE (TOUT, T, Y, ITASK, IER)
 * The arguments are:
 * TOUT  = next value of t at which a solution is desired (input)
 * T     = value of t reached by the solver on output
 * Y     = array containing the computed solution on output
 * ITASK = task indicator: 1 = normal mode (overshoot TOUT and interpolate)
 *         2 = one-step mode (return after each internal step taken)
 *         3 = normal tstop mode (like 1, but integration never proceeds past 
 *             TSTOP, which must be specified through a call to FCVSETRIN
 *             using the key 'STOP_TIME')
 *         4 = one step tstop (like 2, but integration never goes past TSTOP)
 * IER   = completion flag: 0 = success, 1 = tstop return, 2 = root return, 
 *         values -1 ... -10 are various failure modes (see CVODE manual).
 * The current values of the optional outputs are available in IOUT and ROUT.
 * 
 * -----------------------------------------------------------------------------
 *
 * (12) Computing solution derivatives: FCVDKY
 * To obtain a derivative of the solution, of order up to the current method
 * order, make the following call:
 *       CALL FCVDKY (T, K, DKY, IER)
 * The arguments are:
 * T   = value of t at which solution derivative is desired, in [TCUR-HU,TCUR].
 * K   = derivative order (0 .le. K .le. QU)
 * DKY = array containing computed K-th derivative of y on return
 * IER = return flag: = 0 for success, < 0 for illegal argument.
 * 
 * -----------------------------------------------------------------------------
 *
 * (13) Memory freeing: FCVFREE 
 * To free the internal memory created by the calls to FCVMALLOC and
 * FNVINITS or FNVINITP, make the call
 *       CALL FCVFREE
 * 
 * =============================================================================
 */

#ifndef _FCVODE_H
#define _FCVODE_H

/* header files  */
#include <cvode/cvode.h>
#include <sundials/sundials_direct.h>  /* definition of type DlsMat   */
#include <sundials/sundials_sparse.h>  /* definition of type SlsMat   */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */

#if defined(SUNDIALS_F77_FUNC)

#define FCV_MALLOC         SUNDIALS_F77_FUNC(fcvmalloc, FCVMALLOC)
#define FCV_REINIT         SUNDIALS_F77_FUNC(fcvreinit, FCVREINIT)
#define FCV_SETIIN         SUNDIALS_F77_FUNC(fcvsetiin, FCVSETIIN)
#define FCV_SETRIN         SUNDIALS_F77_FUNC(fcvsetrin, FCVSETRIN)
#define FCV_EWTSET         SUNDIALS_F77_FUNC(fcvewtset, FCVEWTSET)
#define FCV_DIAG           SUNDIALS_F77_FUNC(fcvdiag, FCVDIAG)
#define FCV_DENSE          SUNDIALS_F77_FUNC(fcvdense, FCVDENSE)
#define FCV_DENSESETJAC    SUNDIALS_F77_FUNC(fcvdensesetjac, FCVDENSESETJAC)
#define FCV_BAND           SUNDIALS_F77_FUNC(fcvband, FCVBAND)
#define FCV_BANDSETJAC     SUNDIALS_F77_FUNC(fcvbandsetjac, FCVBANDSETJAC)
#define FCV_LAPACKDENSE    SUNDIALS_F77_FUNC(fcvlapackdense, FCVLAPACKDENSE)
#define FCV_LAPACKDENSESETJAC   SUNDIALS_F77_FUNC(fcvlapackdensesetjac, FCVLAPACKDENSESETJAC)
#define FCV_LAPACKBAND     SUNDIALS_F77_FUNC(fcvlapackband, FCVLAPACKBAND)
#define FCV_LAPACKBANDSETJAC    SUNDIALS_F77_FUNC(fcvlapackbandsetjac, FCVLAPACKBANDSETJAC)
#define FCV_KLU            SUNDIALS_F77_FUNC(fcvklu, FCVKLU)
#define FCV_KLUREINIT      SUNDIALS_F77_FUNC(fcvklureinit, FCVKLUREINIT)
#define FCV_SUPERLUMT      SUNDIALS_F77_FUNC(fcvsuperlumt, FCVSUPERLUMT)
#define FCV_SPTFQMR        SUNDIALS_F77_FUNC(fcvsptfqmr, FCVSPTFQMR)
#define FCV_SPTFQMRREINIT  SUNDIALS_F77_FUNC(fcvsptfqmrreinit, FCVSPTFQMRREINIT)
#define FCV_SPBCG          SUNDIALS_F77_FUNC(fcvspbcg, FCVSPBCG)
#define FCV_SPBCGREINIT    SUNDIALS_F77_FUNC(fcvspbcgreinit, FCVSPBCGREINIT)
#define FCV_SPGMR          SUNDIALS_F77_FUNC(fcvspgmr, FCVSPGMR)
#define FCV_SPGMRREINIT    SUNDIALS_F77_FUNC(fcvspgmrreinit, FCVSPGMRREINIT)
#define FCV_SPILSSETJAC    SUNDIALS_F77_FUNC(fcvspilssetjac, FCVSPILSSETJAC)
#define FCV_SPILSSETPREC   SUNDIALS_F77_FUNC(fcvspilssetprec, FCVSPILSSETPREC)
#define FCV_CVODE          SUNDIALS_F77_FUNC(fcvode, FCVODE)
#define FCV_DKY            SUNDIALS_F77_FUNC(fcvdky, FCVDKY)
#define FCV_FREE           SUNDIALS_F77_FUNC(fcvfree, FCVFREE)
#define FCV_FUN            SUNDIALS_F77_FUNC(fcvfun, FCVFUN)
#define FCV_DJAC           SUNDIALS_F77_FUNC(fcvdjac, FCVDJAC)
#define FCV_BJAC           SUNDIALS_F77_FUNC(fcvbjac, FCVBJAC)
#define FCV_PSOL           SUNDIALS_F77_FUNC(fcvpsol, FCVPSOL)
#define FCV_PSET           SUNDIALS_F77_FUNC(fcvpset, FCVPSET)
#define FCV_JTIMES         SUNDIALS_F77_FUNC(fcvjtimes, FCVJTIMES)
#define FCV_EWT            SUNDIALS_F77_FUNC(fcvewt, FCVEWT)
#define FCV_GETERRWEIGHTS  SUNDIALS_F77_FUNC(fcvgeterrweights, FCVGETERRWEIGHTS)
#define FCV_GETESTLOCALERR SUNDIALS_F77_FUNC(fcvgetestlocalerr, FCVGETESTLOCALERR)

#else

#define FCV_MALLOC         fcvmalloc_
#define FCV_REINIT         fcvreinit_
#define FCV_SETIIN         fcvsetiin_
#define FCV_SETRIN         fcvsetrin_
#define FCV_EWTSET         fcvewtset_
#define FCV_DIAG           fcvdiag_
#define FCV_DENSE          fcvdense_
#define FCV_DENSESETJAC    fcvdensesetjac_
#define FCV_BAND           fcvband_
#define FCV_BANDSETJAC     fcvbandsetjac_
#define FCV_LAPACKDENSE    fcvlapackdense_
#define FCV_LAPACKDENSESETJAC   fcvlapackdensesetjac_
#define FCV_LAPACKBAND     fcvlapackband_
#define FCV_LAPACKBANDSETJAC    fcvlapackbandsetjac_
#define FCV_KLU            fcvklu_
#define FCV_KLUREINIT      fcvklureinit_
#define FCV_SUPERLUMT      fcvsuperlumt_
#define FCV_SPTFQMR        fcvsptfqmr_
#define FCV_SPTFQMRREINIT  fcvsptfqmrreinit_
#define FCV_SPBCG          fcvspbcg_
#define FCV_SPBCGREINIT    fcvspbcgreinit_
#define FCV_SPGMR          fcvspgmr_
#define FCV_SPGMRREINIT    fcvspgmrreinit_
#define FCV_SPILSSETJAC    fcvspilssetjac_
#define FCV_SPILSSETPREC   fcvspilssetprec_
#define FCV_CVODE          fcvode_
#define FCV_DKY            fcvdky_
#define FCV_FREE           fcvfree_
#define FCV_FUN            fcvfun_
#define FCV_DJAC           fcvdjac_
#define FCV_BJAC           fcvbjac_
#define FCV_PSOL           fcvpsol_
#define FCV_PSET           fcvpset_
#define FCV_JTIMES         fcvjtimes_
#define FCV_EWT            fcvewt_
#define FCV_GETERRWEIGHTS  fcvgeterrweights_
#define FCV_GETESTLOCALERR fcvgetestlocalerr_

#endif

  /* Type for user data */

  typedef struct {
    realtype *rpar;
    long int *ipar;
  } *FCVUserData;

  /* Prototypes of exported functions */

  void FCV_MALLOC(realtype *t0, realtype *y0,
                  int *meth, int *itmeth, int *iatol,
                  realtype *rtol, realtype *atol,
                  long int *iout, realtype *rout,
                  long int *ipar, realtype *rpar,
                  int *ier);

  void FCV_REINIT(realtype *t0, realtype *y0,
                  int *iatol, realtype *rtol, realtype *atol,
                  int *ier);

  void FCV_SETIIN(char key_name[], long int *ival, int *ier);

  void FCV_SETRIN(char key_name[], realtype *rval, int *ier);

  void FCV_EWTSET(int *flag, int *ier);

  void FCV_DIAG(int *ier);

  void FCV_DENSE(long int *neq, int *ier);
  void FCV_DENSESETJAC(int *flag, int *ier);

  void FCV_BAND(long int *neq, long int *mupper, long int *mlower, int *ier);
  void FCV_BANDSETJAC(int *flag, int *ier);

  void FCV_LAPACKDENSE(int *neq, int *ier);
  void FCV_LAPACKDENSESETJAC(int *flag, int *ier);
  void FCV_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier);
  void FCV_LAPACKBANDSETJAC(int *flag, int *ier);

  void FCV_KLU(int *neq, int *nnz, int *ordering, int *ier);
  void FCV_KLUREINIT(int *neq, int *nnz, int *reinit_type, int *ier);
  void FCV_SUPERLUMT(int *nthreads, int *neq, int *nnz, int *ordering, int *ier);

  void FCV_SPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier);
  void FCV_SPGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier);

  void FCV_SPBCG(int *pretype, int *maxl, realtype *delt, int *ier);
  void FCV_SPBCGREINIT(int *pretype, int *maxl, realtype *delt, int *ier);

  void FCV_SPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier);
  void FCV_SPTFQMRREINIT(int *pretype, int *maxl, realtype *delt, int *ier);

  void FCV_SPILSSETJAC(int *flag, int *ier);
  void FCV_SPILSSETPREC(int *flag, int *ier);
  
  void FCV_CVODE(realtype *tout, realtype *t, realtype *y, int *itask, int *ier);

  void FCV_DKY(realtype *t, int *k, realtype *dky, int *ier);

  void FCV_GETERRWEIGHTS(realtype *eweight, int *ier);
  void FCV_GETESTLOCALERR(realtype *ele, int *ier);

  void FCV_FREE(void);


  /* Prototypes: Functions Called by the CVODE Solver */
  
  int FCVf(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  
  int FCVDenseJac(long int N, realtype t, 
                  N_Vector y, N_Vector fy, 
                  DlsMat J, void *user_data,
                  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FCVBandJac(long int N, long int mupper, long int mlower,
                 realtype t, N_Vector y, N_Vector fy,
                 DlsMat J, void *user_data,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FCVLapackDenseJac(long int N, realtype t,
                        N_Vector y, N_Vector fy, 
                        DlsMat Jac, void *user_data,
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  int FCVLapackBandJac(long int N, long int mupper, long int mlower,
                       realtype t, N_Vector y, N_Vector fy, 
                       DlsMat Jac, void *user_data,
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  int FCVPSet(realtype tn, N_Vector y,N_Vector fy, booleantype jok,
              booleantype *jcurPtr, realtype gamma, void *user_data,
              N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  
  int FCVPSol(realtype tn, N_Vector y, N_Vector fy, 
              N_Vector r, N_Vector z,
              realtype gamma, realtype delta,
              int lr, void *user_data, N_Vector vtemp);
  
  int FCVJtimes(N_Vector v, N_Vector Jv, realtype t, 
                N_Vector y, N_Vector fy,
                void *user_data, N_Vector work);
  
  int FCVEwtSet(N_Vector y, N_Vector ewt, void *user_data);

  /* Declarations for global variables shared amongst various routines */

  extern N_Vector F2C_CVODE_vec;   /* defined in FNVECTOR module */

  extern void *CV_cvodemem;        /* defined in fcvode.c */
  extern long int *CV_iout;        /* defined in fcvode.c */
  extern realtype *CV_rout;        /* defined in fcvode.c */
  extern int CV_nrtfn;             /* defined in fcvode.c */
  extern int CV_ls;                /* defined in fcvode.c */

  /* Linear solver IDs */

  enum { CV_LS_DENSE = 1, CV_LS_BAND = 2, CV_LS_DIAG = 3,
         CV_LS_LAPACKDENSE = 4, CV_LS_LAPACKBAND = 5,
	 CV_LS_KLU = 6, CV_LS_SUPERLUMT = 7, 
	 CV_LS_SPGMR = 8, CV_LS_SPBCG = 9, CV_LS_SPTFQMR = 10 };

#ifdef __cplusplus
}
#endif

#endif
