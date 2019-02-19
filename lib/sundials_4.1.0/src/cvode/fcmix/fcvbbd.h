/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *    Alan Hindmarsh, Radu Serban and Aaron Collier @ LLNL
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
 *---------------------------------------------------------------
 * This is the Fortran interface include file for the BBD
 * preconditioner (CVBBDPRE)
 * -----------------------------------------------------------------
 */

/*
 * ==============================================================================
 *
 *                   FCVBBD Interface Package
 *
 * The FCVBBD Interface Package is a package of C functions which,
 * together with the FCVODE Interface Package, support the use of the
 * CVODE solver and MPI-parallel N_Vector module, along with the CVBBDPRE
 * preconditioner module, for the solution of ODE systems in a mixed
 * Fortran/C setting.  The combination of CVODE and CVBBDPRE solves systems
 * dy/dt = f(t,y) using a Krylov iterative linear solver via the CVSPILS
 * interface, and with a preconditioner that is block-diagonal with banded blocks.
 * While CVODE and CVBBDPRE are written in C, it is assumed here that the user's
 * calling program and user-supplied problem-defining routines are written in
 * Fortran.
 *
 * The user-callable functions in this package, with the corresponding
 * CVODE and CVBBDPRE functions, are as follows:
 *
 *   Fortran               CVODE
 *   --------------        ---------------------------
 *   FCVBBDININT           CVBBDPrecInit
 *   FCVBBDREINIT          CVBBDPrecReInit
 *   FCVBBDOPT             (accesses optional outputs)
 *   --------------        ---------------------------
 *
 * In addition to the Fortran right-hand side function FCVFUN, the
 * user-supplied functions used by this package, are listed below,
 * each with the corresponding interface function which calls it (and its
 * type within CVBBDPRE or CVODE):
 *
 *   Fortran           CVODE            Type
 *   --------------    -----------      -----------------
 *   FCVLOCFN          FCVgloc          CVLocalFn
 *   FCVCOMMF          FCVcfn           CVCommFn
 *   FCVJTSETUP(*)     FCVJTSetup       CVSpilsJTSetupFn
 *   FCVJTIMES(*)      FCVJtimes        CVSpilsJtimesFn
 *   --------------    -----------      -----------------
 *   (*) = optional
 *
 * Important notes on portability:
 *
 * The names of all user-supplied routines here are fixed, in order to
 * maximize portability for the resulting mixed-language program.
 *
 * In this package, the names of the interface functions, and the names of
 * the Fortran user routines called by them, appear as dummy names
 * which are mapped to actual values by a series of definitions in the
 * header file fcvbbd.h.
 *
 * ==============================================================================
 *
 *               Usage of the FCVODE/FCVBBD Interface Packages
 *
 * The usage of the combined interface packages FCVODE and FCVBBD requires
 * calls to a variety of interface functions, and three or more user-supplied
 * routines which define the problem to be solved and indirectly define
 * the preconditioner.  These function calls and user routines are
 * summarized separately below.
 *
 * Some details are omitted, and the user is referred to the CVODE user document
 * for more complete information.
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
 *       YDOT -- array containing state derivatives [realtype,
 *               output]
 *       IPAR -- array containing integer user data that was passed
 *               to FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful,
 *                 >0 if a recoverable error occurred,
 *                 <0 if an unrecoverable error ocurred.
 *
 * (2) User-supplied routines to define preconditoner: FCVLOCFN and FCVCOMMF
 *
 *   The routines in the CVBBDPRE module provide a preconditioner matrix
 *   for CVODE that is block-diagonal with banded blocks.  The blocking
 *   corresponds to the distribution of the dependent variable vector y
 *   among the processors.  Each preconditioner block is generated from the
 *   Jacobian of the local part (on the current processor) of a given
 *   function g(t,y) approximating f(t,y).  The blocks are generated by a
 *   difference quotient scheme on each processor independently, utilizing
 *   an assumed banded structure with given half-bandwidths.  A separate
 *   pair of half-bandwidths defines the band matrix retained.
 *
 * (2.1) Local approximate function FCVLOCFN.
 *
 *   The user must supply a subroutine of the form
 *
 *       SUBROUTINE FCVLOCFN (NLOC, T, YLOC, GLOC, IPAR, RPAR, IER)
 *
 *   To compute the function g(t,y) which approximates the right-hand side
 *   function f(t,y).  This function is to be computed locally, i.e. without
 *   interprocess communication.  (The case where g is mathematically
 *   identical to f is allowed.)
 *
 *   The arguments are:
 *       NLOC -- local problem size [long int, input]
 *       T    -- current time [realtype, input]
 *       YLOC -- array containing local state variables
 *               [realtype, input]
 *       GLOC -- array containing local state derivatives
 *               [realtype, output]
 *       IPAR -- array containing integer user data that was passed
 *               to FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful,
 *                 >0 if a recoverable error occurred,
 *                 <0 if an unrecoverable error ocurred.
 *
 * (2.2) Communication function FCVCOMMF.
 *
 *   The user must also supply a subroutine of the form
 *
 *       SUBROUTINE FCVCOMMF (NLOC, T, YLOC, IPAR, RPAR, IER)
 *
 *   which is to perform all interprocess communication necessary to
 *   evaluate the approximate right-hand side function g described above.
 *   This function takes as input the local vector length NLOC, the
 *   independent variable value T = t, and the local real dependent
 *   variable array YLOC.  It is expected to save communicated data in
 *   work space defined by the user, and made available to CVLOCFN.
 *   Each call to the FCVCOMMF is preceded by a call to FCVFUN with the same
 *   (t,y) arguments.  Thus FCVCOMMF can omit any communications done by
 *   FCVFUN if relevant to the evaluation of g.
 *
 *   The arguments are:
 *       NLOC -- local problem size [long int, input]
 *       T    -- current time [realtype, input]
 *       YLOC -- array containing local state variables
 *               [realtype, input]
 *       IPAR -- array containing integer user data that was passed
 *               to FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful,
 *                 >0 if a recoverable error occurred,
 *                 <0 if an unrecoverable error ocurred.
 *
 * (3) Optional user-supplied Jacobian-vector setup and product
 *   functions: FCVJTSETUP and FCVJTIMES
 *
 *   As an option, the user may supply a routine that computes the product
 *   of the system Jacobian J = df/dy and a given vector v.  If supplied, a
 *   'setup' routine to prepare any user data structures must exist, and
 *   have the form:
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
 *   The accompanying Jacobian matrix-vector product routine must
 *   have the following form:
 *
 *       SUBROUTINE FCVJTIMES (V, FJV, T, Y, FY, EWT, IPAR, RPAR, WORK, IER)
 *
 *   Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
 *   compute the product vector Jv, where the vector v is stored in V, and store
 *   the product in FJV.
 *
 *     The arguments are:
 *       V    -- vector to multiply [realtype, input]
 *       FJV  -- product vector [realtype, output]
 *       T    -- current time [realtype, input]
 *       Y    -- state variables [realtype, input]
 *       FY   -- state derivatives [realtype, input]
 *       H    -- current step size [realtype, input]
 *       IPAR -- array containing integer user data that was passed
 *               to FCVMALLOC [long int, input]
 *       RPAR -- array containing real user data that was passed to
 *               FCVMALLOC [realtype, input]
 *       WORK -- array containing temporary workspace of same size
 *               as Y [realtype, input]
 *       IER  -- return flag [int, output]:
 *                  0 if successful,
 *                  nonzero if an error.
 *
 * (4) Initialization:  FNVINITP, generic iterative linear solver
 *   initialization, FCVMALLOC, FCVSPILSINIT, and FCVBBDINIT.
 *
 * (4.1) To initialize the parallel vector specification, the user must make
 *   the following call:
 *
 *        CALL FNVINITP(COMM, 1, NLOCAL, NGLOBAL, IER)
 *
 *   where the second argument is an int containing the CVODE
 *   solver ID (1). The other arguments are:
 *        COMM = the MPI communicator [int, input]
 *        NLOCAL = local vector size on this processor
 *           [long int, input]
 *        NGLOBAL = system size, and the global size of vectors
 *           (the sum of all values of NLOCAL) [long int, input]
 *        IER = return completion flag [int, ouptut].
 *                  0 = success,
 *                 -1 = failure.
 *
 * (4.2) To initialize a generic iterative linear solver structure for
 *   solving linear systems arising from implicit or IMEX treatment
 *   of the IVP, the user must make one of the following calls:
 *
 *          CALL FSUNPCGINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSPBCGSINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSPFGMRINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSPGMRINIT(1, PRETYPE, MAXL, IER)
 *          CALL FSUNSPTFQMRINIT(1, PRETYPE, MAXL, IER)
 *
 *   In each of these, one argument is an int containing the CVODE solver
 *   ID (1).
 *
 *     The other arguments are:
 *
 *        PRETYPE = type of preconditioning to perform (0=none, 1=left,
 *           2=right, 3=both) [int, input]
 *        MAXL = maximum Krylov subspace dimension [int, input]
 *        IER = return completion flag [int, output]:
 *	          0 = success,
 *		 -1 = failure.
 *
 * (4.3) To set various problem and solution parameters and allocate
 * internal memory for CVODE, make the following call:
 *
 *       CALL FCVMALLOC(T0, Y0, METH, IATOL, RTOL, ATOL,
 *      1               IOUT, ROUT, IPAR, RPAR, IER)
 *
 * The arguments are:
 *      T0     = initial value of t [realtype, input]
 *      Y0     = array of initial conditions [realtype, input]
 *      METH   = flag denoting integration method [int, input]:
 *                  1 = Adams (nonstiff),
 *                  2 = BDF (stiff)
 *      IATOL  = flag denoting type for absolute tolerance ATOL [int, input]:
 *                  1 = scalar,
 *                  2 = array
 *      RTOL   = scalar relative tolerance [realtype, input]
 *      ATOL   = scalar or array absolute tolerance [realtype, input]
 *      IOUT   = array of length at least 21 for integer optional outputs
 *               [long int, output]
 *      ROUT   = array of length at least 6 for real optional outputs
 *               [realtype, output]
 *      IPAR   = array with user integer data [long int, in/out]
 *      RPAR   = array with user real data [realtype, in/out]
 *      IER    = return completion flag [int, output]:
 *                  0 = success,
 *                 -1 = failure (see printed message for details).
 *
 *   The user data arrays IPAR and RPAR are passed unmodified to
 *   all subsequent calls to user-provided routines. Changes to
 *   either array inside a user-provided routine will be
 *   propagated. Using these two arrays, the user can dispense
 *   with COMMON blocks to pass data betwen user-provided
 *   routines.
 *
 * (4.3) Create the CVSPILS interface to attach the generic
 *   iterative linear solver to CVode, by making the following call:
 *
 *       CALL FCVSPILSINIT(IER)
 *
 *   The arguments are:
 *	IER = error return flag [int, output]:
 *	       0 = success;
 *	      <0 = an error occured
 *
 * (4.4) To allocate memory and initialize data associated with the CVBBDPRE
 *   preconditioner, make the following call:
 *
 *       CALL FCVBBDINIT(NLOCAL, MUDQ, MLDQ, MU, ML, DQRELY, IER)
 *
 * The arguments are:
 *        NLOCAL = local vector size on this process
 *             [long int, input]
 *        MUDQ = upper half-bandwidth to be used in the computation
 *             of the local Jacobian blocks by difference
 *             quotients.  These may be smaller than the true
 *             half-bandwidths of the Jacobian of the local block
 *             of g, when smaller values may provide greater
 *             efficiency [long int, input]
 *        MLDQ = lower half-bandwidth to be used in the computation
 *             of the local Jacobian blocks by difference
 *             quotients [long int, input]
 *        MU = upper half-bandwidth of the band matrix that is
 *             retained as an approximation of the local Jacobian
 *             block (may be smaller than MUDQ) [long int, input]
 *        ML = lower half-bandwidth of the band matrix that is
 *             retained as an approximation of the local Jacobian
 *             block (may be smaller than MLDQ) [long int, input]
 *        DQRELY = relative increment factor in y for difference
 *             quotients [realtype, input]
 *                    0.0 = default (sqrt(unit roundoff))
 *        IER = return completion flag [int, output]:
 *                    0 = success
 *                   <0 = an error occurred
 *
 * (4.5) To specify whether the Krylov linear solver should use the
 *   supplied FCVJTIMES or the internal finite difference approximation,
 *   make the call
 *
 *        CALL FCVSPILSSETJAC(FLAG, IER)
 *
 *   with the int FLAG=1 to specify that FCVJTSETUP and FCVJTIMES
 *   are provided (FLAG=0 specifies to use and internal finite
 *   difference approximation to this product).  The int return
 *   flag IER=0 if successful, and nonzero otherwise.
 *
 * (5) Re-initialization: FCVREINIT, FCVBBDREINIT
 *
 *   If a sequence of problems of the same size is being solved using the
 *   Krylov linear solver in combination with the CVBBDPRE preconditioner,
 *   then the CVODE package can be reinitialized for the second and
 *   subsequent problems so as to avoid further memory allocation.  First,
 *   in place of the call to FCVMALLOC, make the following call:
 *
 *       CALL FCVREINIT(T0, Y0, IATOL, RTOL, ATOL, IER)
 *
 *   The arguments have the same names and meanings as those of FCVMALLOC, except
 *   that METH has been omitted from the argument list (being unchanged
 *   for the new problem).  FCVREINIT performs the same initializations as
 *   FCVMALLOC, but does no memory allocation, using instead the existing
 *   internal memory created by the previous FCVMALLOC call.
 *
 *   If there is no change in any of the linear solver or
 *   preconditioner arguments, then no additional calls are
 *   necessary.
 *
 *   Following the call to FCVREINIT, if there is no change in any of the
 *   linear solver arguments, but the user wishes to modify the values of
 *   MUDQ, MLDQ or DQRELY from the previous call to FCVBBDINIT, then a user
 *   may call:
 *
 *      CALL FCVBBDREINIT(MUDQ, MLDQ, DQRELY, IER)
 *
 *   This reinitializes the BBD preconditioner, but without reallocating
 *   its memory. The arguments of the have the same names and meanings as
 *   FCVBBDINIT.
 *
 *   However, if there is a change in any of the linear solver
 *   arguments or other preconditioner arguments, then a call to
 *   FSUNPCGINIT, FSUNSPBCGSINIT, FSUNSPFGMRINIT, FSUNSPGMRINIT,
 *   or FSUNSPTFQMRINIT is required; in this case the linear
 *   solver memory is reallocated.  Following this call, the
 *   CVSPILS interface must also be reconstructed using another
 *   call to FCVSPILSINIT (interface memory is freed and
 *   reallocated), as well as a subsequent call to FCVBBDINIT.
 *
 *
 * (6) The integrator: FCVODE
 *
 *   Carrying out the integration is accomplished by making calls as follows:
 *
 *       CALL FCVODE (TOUT, T, Y, ITASK, IER)
 *
 *   The arguments are:
 *      TOUT  = next value of t at which a solution is desired [realtype, input]
 *      T     = value of t reached by the solver [realtype, output]
 *      Y     = array containing the computed solution [realtype, output]
 *      ITASK = task indicator [int, input]:
 *              1 = normal mode (overshoot TOUT and interpolate)
 *              2 = one-step mode (return after each internal step taken)
 *              3 = normal mode with TSTOP check
 *              4 = one-step mode with TSTOP check
 *      IER   = completion flag [int, output]:
 *              0 = success,
 *              1 = TSTOP return,
 *              2 = root return,
 *              negative values are various failure modes (see CVODE User Guide).
 *    The current values of the optional outputs are available in IOUT and ROUT.
 *
 * (7) Optional outputs: FCVBBDOPT
 *
 *   Optional outputs specific to the CVSPILS solver interface are
 *        LENRWLS = IOUT(13) from CVSpilsGetWorkSpace
 *        LENIWLS = IOUT(14) from CVSpilsGetWorkSpace
 *        LSTF    = IOUT(15) from CVSpilsGetLastFlag
 *        NFELS   = IOUT(16) from CVSpilsGetNumRhsEvals
 *        NJTV    = IOUT(17) from CVSpilsGetNumJtimesEvals
 *        NPE     = IOUT(18) from CVSpilsGetNumPrecEvals
 *        NPS     = IOUT(19) from CVSpilsGetNumPrecSolves
 *        NLI     = IOUT(20) from CVSpilsGetNumLinIters
 *        NCFL    = IOUT(21) from CVSpilsGetNumConvFails
 *   See the CVODE manual for descriptions.
 *
 *   To obtain the optional outputs associated with the CVBBDPRE module, make
 *   the following call:
 *
 *       CALL FCVBBDOPT (LENRWBBD, LENIWBBD, NGEBBD)
 *
 *   The arguments returned are:
 *      LENRWBBD = length of real preconditioner work space, in realtype words.
 *                 This size is local to the current processor [long int, output]
 *      LENIWBBD = length of integer preconditioner work space, in integer words.
 *                 This size is local to the current processor [long int, output]
 *      NGEBBD   = number of g(t,y) evaluations (calls to CVLOCFN) so far
 *                 [long int, output]
 *
 * (8) Computing solution derivatives: FCVDKY
 *
 *   To obtain a derivative of the solution (optionally), of order up to
 *   the current method order, make the following call:
 *
 *       CALL FCVDKY (T, K, DKY)
 *
 *   The arguments are:
 *       T = time at which solution derivative is desired, within
 *           the interval [TCUR-HU,TCUR], [realtype, input].
 *       K = derivative order (0 .le. K .le. QU) [int, input]
 *       DKY = array containing computed K-th derivative of y
 *           [realtype, output]
 *       IER = return flag [int, output]:
 *                    0 = success
 *                   <0 = illegal argument.
 *
 * (9) Memory freeing: FCVFREE
 *
 *   To the free the internal memory created by the calls to FNVINIT*,
 *   FCVMALLOC, FCVSPILSINIT and FCVBBDINIT, make the following call:
 *
 *       CALL FCVFREE
 *
 * ==============================================================================
 */

#ifndef _FCVBBD_H
#define _FCVBBD_H

/* header files  */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */

#if defined(SUNDIALS_F77_FUNC)

#define FCV_BBDINIT    SUNDIALS_F77_FUNC(fcvbbdinit, FCVBBDINIT)
#define FCV_BBDREINIT  SUNDIALS_F77_FUNC(fcvbbdreinit, FCVBBDREINIT)
#define FCV_BBDOPT     SUNDIALS_F77_FUNC(fcvbbdopt, FCVBBDOPT)
#define FCV_GLOCFN     SUNDIALS_F77_FUNC(fcvglocfn, FCVGLOCFN)
#define FCV_COMMFN     SUNDIALS_F77_FUNC(fcvcommfn, FCVCOMMFN)

#else

#define FCV_BBDINIT    fcvbbdinit_
#define FCV_BBDREINIT  fcvbbdreinit_
#define FCV_BBDOPT     fcvbbdopt_
#define FCV_GLOCFN     fcvglocfn_
#define FCV_COMMFN     fcvcommfn_

#endif

/* Prototypes of exported functions */

void FCV_BBDINIT(long int *Nloc, long int *mudq,
                 long int *mldq, long int *mu,
                 long int *ml, realtype* dqrely, int *ier);
void FCV_BBDREINIT(long int *mudq, long int *mldq,
                   realtype* dqrely, int *ier);
void FCV_BBDOPT(long int *lenrwbbd, long int *leniwbbd,
                long int *ngebbd);

/* Prototypes: Functions Called by the CVBBDPRE Module */

int FCVgloc(long int Nloc, realtype t, N_Vector yloc,
            N_Vector gloc, void *user_data);

int FCVcfn(long int Nloc, realtype t, N_Vector y,
           void *user_data);

#ifdef __cplusplus
}
#endif

#endif
