/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban and Aaron Collier @ LLNL
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
 * This is the Fortran interface include file for the BAND
 * preconditioner (CVBANDPRE).
 * -----------------------------------------------------------------
 */

/*
 * ==============================================================================
 *
 *                  FCVBP Interface Package
 *
 * The FCVBP Interface Package is a package of C functions which,
 * together with the FCVODE Interface Package, support the use of the
 * CVODE solver and serial, OpenMP or PThreads vector module with the
 * CVBANDPRE preconditioner module, for the solution of ODE systems in
 * a mixed Fortran/C setting.  The combination of CVODE and CVBANDPRE solves
 * systems dy/dt = f(t,y) using a Krylov iterative linear solver with a banded
 * difference quotient Jacobian-based preconditioner.
 *
 * The user-callable functions in this package, with the corresponding
 * CVODE and CVBBDPRE functions, are as follows:
 *
 *   Fortran              CVODE
 *   -------------        ---------------------------
 *   FCVBPINIT            CVBandPrecInit
 *   FCVBPOPT             (accesses optional outputs)
 *   -------------        ---------------------------
 *
 * In addition to the Fortran right-hand side function FCVFUN, the
 * user may (optionally) supply routines FCVJTSETUP and FCVJTIMES which
 * are called by the interface function FCVJTSetup of type CVSpilsJTSetupFn
 * and the interface function FCVJtimes of type CVSpilsJtimesFn.
 *
 * Important notes on portability.
 *
 * The names of all user-supplied routines here are fixed, in order to
 * maximize portability for the resulting mixed-language program.
 *
 * In this package, the names of the interface functions, and the names of
 * the Fortran user routines called by them, appear as dummy names
 * which are mapped to actual values by a series of definitions in the
 * header file fcvbp.h.
 *
 * ==============================================================================
 *
 *               Usage of the FCVODE/FCVBP Interface Packages
 *
 * The usage of the combined interface packages FCVODE and FCVBP requires
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
 *
 * (2) Optional user-supplied Jacobian-vector setup and product
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
 * (3) Initialization:  FNVINITP, generic iterative linear solver
 *   initialization, FCVMALLOC, FCVSPILSINIT, and FCVBPINIT.
 *
 * (3.1) To initialize the vector specification, the user must make
 *   one of the following calls:
 *
 *       (serial)
 *          CALL FNVINITS(4, NEQ, IER)
 *       (OpenMP threaded)
 *          CALL FNVINITOMP(4, NEQ, NUM_THREADS, IER)
 *       (PThreads threaded)
 *          CALL FNVINITPTS(4, NEQ, NUM_THREADS, IER)
 *
 *   where the first argument is an int containing the CVODE
 *   solver ID (4). The other arguments are:
 *        NEQ = size of vectors [long int, input]
 *        NUM_THREADS = number of threads
 *        IER = return completion flag [int, output]:
 *	          0 = success,
 *		 -1 = failure.
 *
 * (3.2) To initialize a generic iterative linear solver structure for
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
 *
 * (3.2) To set various problem and solution parameters and allocate
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
 * (3.3) Create the CVSPILS interface to attach the generic
 *   iterative linear solver to CVode, by making the following call:
 *
 *       CALL FCVSPILSINIT(IER)
 *
 *   The arguments are:
 *	IER = error return flag [int, output]:
 *	       0 = success;
 *	      <0 = an error occured
 *
 * (3.4) To allocate memory and initialize data associated with the CVBANDPRE
 *   preconditioner, make the following call:
 *
 *       CALL FCVBPINIT(NEQ, MU, ML, IER)
 *
 *   The arguments are:
 *        NEQ = problem size [long int, input]
 *        MU = upper half-bandwidth of the band matrix that is retained as
 *             an approximation of the Jacobian [long int, input]
 *        ML = lower half-bandwidth of the band matrix approximant
 *             to the Jacobian [long int, input]
 *        IER = return completion flag [int, output]:
 *                    0 = success
 *                   <0 = an error occurred
 *
 *
 * (3.5) To specify whether the Krylov linear solver should use the
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
 * (4) The integrator: FCVODE
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
 * (5) Optional outputs: FCVBPOPT
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
 *   To obtain the optional outputs associated with the CVBANDPRE module, make
 *   the following call:
 *
 *       CALL FCVBPOPT(LENRWBP, LENIWBP, NFEBP)
 *
 *   The arguments returned are:
 *      LENRWBP = length of real preconditioner work space, in realtype words.
 *                This size is local to the current processor.
 *      LENIWBP = length of integer preconditioner work space, in integer words.
 *                This size is local to the current processor.
 *      NFEBP   = number of f(t,y) evaluations for CVBANDPRE
 *
 * (6) Computing solution derivatives: FCVDKY
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
 * (7) Memory freeing: FCVFREE
 *
 *   To the free the internal memory created by the calls to FNVINIT*,
 *   FCVMALLOC, FCVSPILSINIT and FCVBPINIT, make the following call:
 *
 *       CALL FCVFREE
 *
 * ==============================================================================
 */

#ifndef _FCVBP_H
#define _FCVBP_H

/* header files  */
#include <sundials/sundials_nvector.h> /* definition of type N_Vector */
#include <sundials/sundials_types.h>   /* definition of type realtype */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* Definitions of interface function names */

#if defined(SUNDIALS_F77_FUNC)

#define FCV_BPINIT    SUNDIALS_F77_FUNC(fcvbpinit, FCVBPINIT)
#define FCV_BPOPT     SUNDIALS_F77_FUNC(fcvbpopt, FCVBPOPT)

#else

#define FCV_BPINIT    fcvbpinit_
#define FCV_BPOPT     fcvbpopt_

#endif

/* Prototypes of exported function */
void FCV_BPINIT(long int *N, long int *mu,
                long int *ml, int *ier);
void FCV_BPOPT(long int *lenrwbp, long int *leniwbp,
               long int *nfebp);

#ifdef __cplusplus
}
#endif

#endif
