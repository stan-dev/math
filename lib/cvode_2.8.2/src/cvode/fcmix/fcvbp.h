/*
 * -----------------------------------------------------------------
 * $Revision: 4378 $
 * $Date: 2015-02-19 10:55:14 -0800 (Thu, 19 Feb 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
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
 * CVODE solver (serial version) with the CVBANDPRE preconditioner module,
 * for the solution of ODE systems in a mixed Fortran/C setting.  The
 * combination of CVODE and CVBANDPRE solves systems dy/dt = f(t,y) with the
 * SPGMR (scaled preconditioned GMRES), SPTFQMR (scaled preconditioned TFQMR),
 *  or SPBCG (scaled preconditioned Bi-CGSTAB) method for the linear systems
 * that arise, and with a banded difference quotient Jacobian-based preconditioner.
 * 
 * The user-callable functions in this package, with the corresponding
 * CVODE and CVBBDPRE functions, are as follows: 
 *   FCVBPINIT    interfaces to CVBandPrecInit
 *   FCVBPOPT     accesses optional outputs
 * 
 * In addition to the Fortran right-hand side function FCVFUN, the
 * user may (optionally) supply a routine FCVJTIMES which is called by 
 * the interface function FCVJtimes of type CVSpilsJtimesFn.
 * (The names of all user-supplied routines here are fixed, in order to
 * maximize portability for the resulting mixed-language program.)
 * 
 * Important note on portability.
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
 * calls to seven to ten interface functions, and one or two user-supplied
 * routines which define the problem to be solved and indirectly define
 * the preconditioner.  These function calls and user routines are
 * summarized separately below.
 * 
 * Some details are omitted, and the user is referred to the CVODE user document 
 * for more complete information.
 * 
 * (1) User-supplied right-hand side routine: FCVFUN
 * The user must in all cases supply the following Fortran routine
 *       SUBROUTINE FCVFUN (T, Y, YDOT, IPAR, RPAR, IER)
 *       DIMENSION Y(*), YDOT(*), IPAR(*), RPAR(*)
 * It must set the YDOT array to f(t,y), the right-hand side of the ODE
 * system, as function of T = t and the array Y = y.  Here Y and YDOT
 * are distributed vectors.
 * 
 * (2) Optional user-supplied Jacobian-vector product routine: FCVJTIMES
 * As an option, the user may supply a routine that computes the product
 * of the system Jacobian J = df/dy and a given vector v.  If supplied, it
 * must have the following form:
 *       SUBROUTINE FCVJTIMES (V, FJV, T, Y, FY, EWT, IPAR, RPAR, WORK, IER)
 *       DIMENSION V(*), FJV(*), Y(*), FY(*), EWT(*), IPAR(*), RPAR(*), WORK(*)
 * Typically this routine will use only NEQ, T, Y, V, and FJV.  It must
 * compute the product vector Jv, where the vector v is stored in V, and store
 * the product in FJV.  On return, set IER = 0 if FCVJTIMES was successful,
 * and nonzero otherwise.
 * 
 * (3) Initialization:  FNVINITS, FCVMALLOC, FCVBPINIT.
 * 
 * (3.1) To initialize the serial vector specification, the user must make 
 * the following call:
 *        CALL FNVINITS(NEQ, IER)
 * where NEQ is the problem size and IER is a return completion flag.
 * Possible values for IER are 0 = success, -1 = failure.
 * 
 * (3.2) To set various problem and solution parameters and allocate
 * internal memory for CVODE, make the following call:
 *       CALL FCVMALLOC(T0, Y0, METH, ITMETH, IATOL, RTOL, ATOL,
 *      1               IOUT, ROUT, IPAR, RPAR, IER)
 * The arguments are:
 * T0     = initial value of t
 * Y0     = array of initial conditions
 * METH   = basic integration method: 1 = Adams (nonstiff), 2 = BDF (stiff)
 * ITMETH = nonlinear iteration method: 1 = functional iteration, 2 = Newton iter.
 * IATOL  = type for absolute tolerance ATOL: 1 = scalar, 2 = array
 * RTOL   = relative tolerance (scalar)
 * ATOL   = absolute tolerance (scalar or array)
 * IOUT   = array of length 21 for integer optional outputs
 *          (declare as INTEGER*4 or INTEGER*8 according to C type long int)
 * ROUT   = array of length 6 for real optional outputs
 * IPAR   = array with user integer data
 *          (declare as INTEGER*4 or INTEGER*8 according to C type long int)
 * RPAR   = array with user real data
 * IER    = return completion flag.  Values are 0 = success, and -1 = failure.
 *          See printed message for details in case of failure.
 * 
 * (3.3) To allocate memory and initialize data associated with the CVBANDPRE
 * preconditioner, make the following call:
 *       CALL FCVBPINIT(NEQ, MU, ML, IER)
 * The arguments are:
 * NEQ       = problem size
 * MU, ML    = upper and lower half-bandwidths of the band matrix that 
 *             is retained as an approximation of the Jacobian.
 * IER       = return completion flag: IER=0: success, IER<0: and error occurred
 *
 * (3.4A) To specify the SPGMR linear solver with the CVBANDPRE preconditioner,
 * make the following call
 *       CALL FCVSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)
 * The arguments are:
 * IPRETYPE  = preconditioner type: 
 *            0 = none
 *            1 = left only
 *            2 = right only
 *            3 = both sides.
 * IGSTYPE   = Gram-schmidt process type: 0 = modified G-S, 1 = classical G-S.
 * MAXL      = maximum Krylov subspace dimension; 0 indicates default.
 * DELT      = linear convergence tolerance factor; 0.0 indicates default.
 * IER       = return completion flag: IER=0: success, IER<0: ans error occurred
 *
 * (3.4B) To specify the SPBCG linear solver with the CVBANDPRE preconditioner,
 * make the following call
 *       CALL FCVSPBCG(IPRETYPE, MAXL, DELT, IER)
 * The arguments are:
 * IPRETYPE  = preconditioner type: 
 *            0 = none
 *            1 = left only
 *            2 = right only
 *            3 = both sides.
 * MAXL      = maximum Krylov subspace dimension; 0 indicates default.
 * DELT      = linear convergence tolerance factor; 0.0 indicates default.
 * IER       = return completion flag: IER=0: success, IER<0: ans error occurred
 *
 * (3.4C) To specify the SPTFQMR linear solver with the CVBANDPRE preconditioner,
 * make the following call
 *       CALL FCVSPTFQMR(IPRETYPE, MAXL, DELT, IER)
 * The arguments are:
 * IPRETYPE  = preconditioner type: 
 *            0 = none
 *            1 = left only
 *            2 = right only
 *            3 = both sides.
 * MAXL      = maximum Krylov subspace dimension; 0 indicates default.
 * DELT      = linear convergence tolerance factor; 0.0 indicates default.
 * IER       = return completion flag: IER=0: success, IER<0: ans error occurred
 *
 * (3.5) To specify whether the Krylov linear solver (GMRES, Bi-CGSTAB, or TFQMR) 
 * should use the supplied FCVJTIMES or the internal finite difference approximation, 
 * make the call
 *        CALL FCVSPILSSETJAC(FLAG, IER)
 * where FLAG=0 for finite differences approxaimtion or
 *       FLAG=1 to use the supplied routine FCVJTIMES
 *
 * (4) The integrator: FCVODE
 * Carrying out the integration is accomplished by making calls as follows:
 *       CALL FCVODE (TOUT, T, Y, ITASK, IER)
 * The arguments are:
 * TOUT  = next value of t at which a solution is desired (input)
 * T     = value of t reached by the solver on output
 * Y     = array containing the computed solution on output
 * ITASK = task indicator: 1 = normal mode (overshoot TOUT and interpolate);
 *         2 = one-step mode (return after each internal step taken);
 *         3 = normal mode with TSTOP; 4 = one-step mode with TSTOP.
 * IER   = completion flag: 0 = success, 1 = TSTOP return, 2 = root return,
 *         negative values are various failure modes (see CVODE User Guide).
 * The current values of the optional outputs are available in IOUT and ROUT.
 * 
 * (5) Optional outputs: FCVBPOPT
 * Optional outputs specific to the SP* solver are LRW, LIW, LFLG, NFELS, NJTV,
 * NPE, NPS, NLI, NCFL, stored in IOUT(13)...IOUT(21).
 * To obtain the optional outputs associated with the CVBANDPRE module, make
 * the following call:
 *       CALL FCVBPOPT(LENRWBP, LENIWBP, NFEBP)
 * The arguments returned are:
 * LENRWBP = length of real preconditioner work space, in realtype words.
 *           This size is local to the current processor.
 * LENIWBP = length of integer preconditioner work space, in integer words.
 *           This size is local to the current processor.
 * NFEBP   = number of f(t,y) evaluations for CVBANDPRE
 * 
 * (6) Computing solution derivatives: FCVDKY
 * To obtain a derivative of the solution (optionally), of order up to
 * the current method order, make the following call:
 *       CALL FCVDKY (T, K, DKY)
 * The arguments are:
 * T   = value of t at which solution derivative is desired
 * K   = derivative order (0 .le. K .le. QU)
 * DKY = array containing computed K-th derivative of y on return
 * 
 * (7) Memory freeing: FCVFREE
 *   To the free the internal memory created by the calls to FNVINITS,
 * FCVMALLOC, and FCVBPINIT, make the following call:
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
void FCV_BPINIT(long int *N, long int *mu, long int *ml, int *ier);
void FCV_BPOPT(long int *lenrwbp, long int *leniwbp, long int *nfebp);

#ifdef __cplusplus
}
#endif

#endif
