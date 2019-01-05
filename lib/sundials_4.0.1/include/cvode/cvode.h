/*
 * -----------------------------------------------------------------
 * $Revision$
 * $Date$
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban
 *                and Dan Shumaker @ LLNL
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
 * This is the interface file for the main CVODE integrator.
 * -----------------------------------------------------------------
 *
 * CVODE is used to solve numerically the ordinary initial value
 * problem:
 *
 *                 y' = f(t,y),
 *                 y(t0) = y0,
 *
 * where t0, y0 in R^N, and f: R x R^N -> R^N are given.
 *
 * -----------------------------------------------------------------
 */

#ifndef _CVODE_H
#define _CVODE_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <cvode/cvode_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * =================================================================
 *              C V O D E     C O N S T A N T S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Enumerations for inputs to CVodeCreate and CVode.
 * -----------------------------------------------------------------
 * Symbolic constants for the lmm parameter to CVodeCreate
 * and the input parameter itask to CVode, are given below.
 *
 * lmm:   The user of the CVODE package specifies whether to use the
 *        CV_ADAMS (Adams-Moulton) or CV_BDF (Backward Differentiation
 *        Formula) linear multistep method. The BDF method is
 *        recommended for stiff problems, and the CV_ADAMS method is
 *        recommended for nonstiff problems.
 *
 * itask: The itask input parameter to CVode indicates the job
 *        of the solver for the next user step. The CV_NORMAL
 *        itask is to have the solver take internal steps until
 *        it has reached or just passed the user specified tout
 *        parameter. The solver then interpolates in order to
 *        return an approximate value of y(tout). The CV_ONE_STEP
 *        option tells the solver to just take one internal step
 *        and return the solution at the point reached by that step.
 * -----------------------------------------------------------------
 */

/* lmm */
#define CV_ADAMS 1
#define CV_BDF   2

/* itask */
#define CV_NORMAL         1
#define CV_ONE_STEP       2

/*
 * ----------------------------------------
 * CVODE return flags
 * ----------------------------------------
 */

#define CV_SUCCESS               0
#define CV_TSTOP_RETURN          1
#define CV_ROOT_RETURN           2

#define CV_WARNING              99

#define CV_TOO_MUCH_WORK        -1
#define CV_TOO_MUCH_ACC         -2
#define CV_ERR_FAILURE          -3
#define CV_CONV_FAILURE         -4

#define CV_LINIT_FAIL           -5
#define CV_LSETUP_FAIL          -6
#define CV_LSOLVE_FAIL          -7
#define CV_RHSFUNC_FAIL         -8
#define CV_FIRST_RHSFUNC_ERR    -9
#define CV_REPTD_RHSFUNC_ERR    -10
#define CV_UNREC_RHSFUNC_ERR    -11
#define CV_RTFUNC_FAIL          -12
#define CV_NLS_INIT_FAIL        -13
#define CV_NLS_SETUP_FAIL       -14
#define CV_CONSTR_FAIL          -15

#define CV_MEM_FAIL             -20
#define CV_MEM_NULL             -21
#define CV_ILL_INPUT            -22
#define CV_NO_MALLOC            -23
#define CV_BAD_K                -24
#define CV_BAD_T                -25
#define CV_BAD_DKY              -26
#define CV_TOO_CLOSE            -27
#define CV_VECTOROP_ERR         -28

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type : CVRhsFn
 * -----------------------------------------------------------------
 * The f function which defines the right hand side of the ODE
 * system y' = f(t,y) must have type CVRhsFn.
 * f takes as input the independent variable value t, and the
 * dependent variable vector y.  It stores the result of f(t,y)
 * in the vector ydot.  The y and ydot arguments are of type
 * N_Vector.
 * (Allocation of memory for ydot is handled within CVODE)
 * The user_data parameter is the same as the user_data
 * parameter set by the user through the CVodeSetUserData routine.
 * This user-supplied pointer is passed to the user's f function
 * every time it is called.
 *
 * A CVRhsFn should return 0 if successful, a negative value if
 * an unrecoverable error occured, and a positive value if a 
 * recoverable error (e.g. invalid y values) occured. 
 * If an unrecoverable occured, the integration is halted. 
 * If a recoverable error occured, then (in most cases) CVODE
 * will try to correct and retry.
 * -----------------------------------------------------------------
 */

typedef int (*CVRhsFn)(realtype t, N_Vector y,
		       N_Vector ydot, void *user_data);

/*
 * -----------------------------------------------------------------
 * Type : CVRootFn
 * -----------------------------------------------------------------
 * A function g, which defines a set of functions g_i(t,y) whose
 * roots are sought during the integration, must have type CVRootFn.
 * The function g takes as input the independent variable value
 * t, and the dependent variable vector y.  It stores the nrtfn
 * values g_i(t,y) in the realtype array gout.
 * (Allocation of memory for gout is handled within CVODE.)
 * The user_data parameter is the same as that passed by the user
 * to the CVodeSetUserData routine.  This user-supplied pointer is
 * passed to the user's g function every time it is called.
 *
 * A CVRootFn should return 0 if successful or a non-zero value
 * if an error occured (in which case the integration will be halted).
 * -----------------------------------------------------------------
 */

typedef int (*CVRootFn)(realtype t, N_Vector y, realtype *gout, void *user_data);

/*
 * -----------------------------------------------------------------
 * Type : CVEwtFn
 * -----------------------------------------------------------------
 * A function e, which sets the error weight vector ewt, must have
 * type CVEwtFn.
 * The function e takes as input the current dependent variable y.
 * It must set the vector of error weights used in the WRMS norm:
 * 
 *   ||y||_WRMS = sqrt [ 1/N * sum ( ewt_i * y_i)^2 ]
 *
 * Typically, the vector ewt has components:
 * 
 *   ewt_i = 1 / (reltol * |y_i| + abstol_i)
 *
 * The user_data parameter is the same as that passed by the user
 * to the CVodeSetUserData routine.  This user-supplied pointer is
 * passed to the user's e function every time it is called.
 * A CVEwtFn e must return 0 if the error weight vector has been
 * successfuly set and a non-zero value otherwise.
 * -----------------------------------------------------------------
 */

typedef int (*CVEwtFn)(N_Vector y, N_Vector ewt, void *user_data);

/*
 * -----------------------------------------------------------------
 * Type : CVErrHandlerFn
 * -----------------------------------------------------------------
 * A function eh, which handles error messages, must have type
 * CVErrHandlerFn.
 * The function eh takes as input the error code, the name of the
 * module reporting the error, the error message, and a pointer to
 * user data, the same as that passed to CVodeSetUserData.
 * 
 * All error codes are negative, except CV_WARNING which indicates 
 * a warning (the solver continues).
 *
 * A CVErrHandlerFn has no return value.
 * -----------------------------------------------------------------
 */

typedef void (*CVErrHandlerFn)(int error_code, 
			       const char *module, const char *function, 
			       char *msg, void *user_data); 

/*
 * =================================================================
 *          U S E R - C A L L A B L E   R O U T I N E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Function : CVodeCreate
 * -----------------------------------------------------------------
 * CVodeCreate creates an internal memory block for a problem to
 * be solved by CVODE.
 *
 * lmm   is the type of linear multistep method to be used.
 *       The legal values are CV_ADAMS and CV_BDF (see previous
 *       description).
 *
 * If successful, CVodeCreate returns a pointer to initialized
 * problem memory. This pointer should be passed to CVodeInit.
 * If an initialization error occurs, CVodeCreate prints an error
 * message to standard err and returns NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void *CVodeCreate(int lmm);

/*
 * -----------------------------------------------------------------
 * Integrator optional input specification functions
 * -----------------------------------------------------------------
 * The following functions can be called to set optional inputs
 * to values other than the defaults given below:
 *
 * Function                |  Optional input / [ default value ]
 * -----------------------------------------------------------------
 *                         |
 * CVodeSetErrHandlerFn    | user-provided ErrHandler function.
 *                         | [internal]
 *                         |
 * CVodeSetErrFile         | the file pointer for an error file
 *                         | where all CVODE warning and error
 *                         | messages will be written if the default
 *                         | internal error handling function is used. 
 *                         | This parameter can be stdout (standard 
 *                         | output), stderr (standard error), or a 
 *                         | file pointer (corresponding to a user 
 *                         | error file opened for writing) returned 
 *                         | by fopen.
 *                         | If not called, then all messages will
 *                         | be written to the standard error stream.
 *                         | [stderr]
 *                         |
 * CVodeSetUserData        | a pointer to user data that will be
 *                         | passed to the user's f function every
 *                         | time f is called.
 *                         | [NULL]
 *                         |
 * CVodeSetMaxOrd          | maximum lmm order to be used by the
 *                         | solver.
 *                         | [12 for Adams , 5 for BDF]
 *                         |
 * CVodeSetMaxNumSteps     | maximum number of internal steps to be
 *                         | taken by the solver in its attempt to
 *                         | reach tout.
 *                         | [500]
 *                         |
 * CVodeSetMaxHnilWarns    | maximum number of warning messages
 *                         | issued by the solver that t+h==t on the
 *                         | next internal step. A value of -1 means
 *                         | no such messages are issued.
 *                         | [10]
 *                         |
 * CVodeSetStabLimDet      | flag to turn on/off stability limit
 *                         | detection (SUNTRUE = on, SUNFALSE = off).
 *                         | When BDF is used and order is 3 or
 *                         | greater, CVsldet is called to detect
 *                         | stability limit.  If limit is detected,
 *                         | the order is reduced.
 *                         | [SUNFALSE]
 *                         |
 * CVodeSetInitStep        | initial step size.
 *                         | [estimated by CVODE]
 *                         |
 * CVodeSetMinStep         | minimum absolute value of step size
 *                         | allowed.
 *                         | [0.0]
 *                         |
 * CVodeSetMaxStep         | maximum absolute value of step size
 *                         | allowed.
 *                         | [infinity]
 *                         |
 * CVodeSetStopTime        | the independent variable value past
 *                         | which the solution is not to proceed.
 *                         | [infinity]
 *                         |
 * CVodeSetMaxErrTestFails | Maximum number of error test failures
 *                         | in attempting one step.
 *                         | [7]
 *                         |
 * CVodeSetMaxNonlinIters  | Maximum number of nonlinear solver
 *                         | iterations at one solution.
 *                         | [3]
 *                         |
 * CVodeSetMaxConvFails    | Maximum number of convergence failures
 *                         | allowed in attempting one step.
 *                         | [10]
 *                         |
 * CVodeSetNonlinConvCoef  | Coefficient in the nonlinear
 *                         | convergence test.
 *                         | [0.1]
 *                         |
 * CVodeSetConstraints	   | an N_Vector defining inequality
 *                         | constraints for each component of the
 *                         | solution vector y. If a given element
 *                         | of this vector has values +2 or -2,
 *                         | then the corresponding component of y
 *                         | will be constrained to be > 0.0 or
 *                         | < 0.0, respectively, while if it is +1
 *                         | or -1, the y component is constrained
 *                         | to be >= 0.0 or <= 0.0, respectively.
 *                         | If a component of constraints is 0.0,
 *                         | then no constraint is imposed on the
 *                         | corresponding component of y.
 *                         | The presence of a non-NULL constraints
 *                         | vector that is not 0.0 (ZERO) in all
 *                         | components will cause constraint
 *                         | checking to be performed.
 *                         |
 * -----------------------------------------------------------------
 *                            |
 * CVodeSetRootDirection      | Specifies the direction of zero
 *                            | crossings to be monitored
 *                            | [both directions]
 *                            |
 * CVodeSetNoInactiveRootWarn | disable warning about possible
 *                            | g==0 at beginning of integration
 *                            | 
 * -----------------------------------------------------------------

 * -----------------------------------------------------------------
 * Return flag:
 *   CV_SUCCESS   if successful
 *   CV_MEM_NULL  if the cvode memory is NULL
 *   CV_ILL_INPUT if an argument has an illegal value
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeSetErrHandlerFn(void *cvode_mem, CVErrHandlerFn ehfun, void *eh_data);
SUNDIALS_EXPORT int CVodeSetErrFile(void *cvode_mem, FILE *errfp);
SUNDIALS_EXPORT int CVodeSetUserData(void *cvode_mem, void *user_data);
SUNDIALS_EXPORT int CVodeSetMaxOrd(void *cvode_mem, int maxord);
SUNDIALS_EXPORT int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps);
SUNDIALS_EXPORT int CVodeSetMaxHnilWarns(void *cvode_mem, int mxhnil);
SUNDIALS_EXPORT int CVodeSetStabLimDet(void *cvode_mem, booleantype stldet);
SUNDIALS_EXPORT int CVodeSetInitStep(void *cvode_mem, realtype hin);
SUNDIALS_EXPORT int CVodeSetMinStep(void *cvode_mem, realtype hmin);
SUNDIALS_EXPORT int CVodeSetMaxStep(void *cvode_mem, realtype hmax);
SUNDIALS_EXPORT int CVodeSetStopTime(void *cvode_mem, realtype tstop);
SUNDIALS_EXPORT int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef);
SUNDIALS_EXPORT int CVodeSetMaxNonlinIters(void *cvode_mem, int maxcor);
SUNDIALS_EXPORT int CVodeSetMaxConvFails(void *cvode_mem, int maxncf);
SUNDIALS_EXPORT int CVodeSetNonlinConvCoef(void *cvode_mem, realtype nlscoef);
SUNDIALS_EXPORT int CVodeSetConstraints(void *cvode_mem, N_Vector constraints);

SUNDIALS_EXPORT int CVodeSetRootDirection(void *cvode_mem, int *rootdir);
SUNDIALS_EXPORT int CVodeSetNoInactiveRootWarn(void *cvode_mem);

/*
 * -----------------------------------------------------------------
 * Function : CVodeInit
 * -----------------------------------------------------------------
 * CVodeInit allocates and initializes memory for a problem to
 * to be solved by CVODE.
 *
 * cvode_mem is pointer to CVODE memory returned by CVodeCreate.
 *
 * f       is the name of the C function defining the right-hand
 *         side function in y' = f(t,y).
 *
 * t0      is the initial value of t.
 *
 * y0      is the initial condition vector y(t0).
 *
 * Return flag:
 *  CV_SUCCESS if successful
 *  CV_MEM_NULL if the cvode memory was NULL
 *  CV_MEM_FAIL if a memory allocation failed
 *  CV_ILL_INPUT f an argument has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0);

/*
 * -----------------------------------------------------------------
 * Function : CVodeReInit
 * -----------------------------------------------------------------
 * CVodeReInit re-initializes CVode for the solution of a problem,
 * where a prior call to CVodeInit has been made with the same
 * problem size N. CVodeReInit performs the same input checking
 * and initializations that CVodeInit does.
 * But it does no memory allocation, assuming that the existing
 * internal memory is sufficient for the new problem.
 *
 * The use of CVodeReInit requires that the maximum method order,
 * maxord, is no larger for the new problem than for the problem
 * specified in the last call to CVodeInit.  This condition is
 * automatically fulfilled if the multistep method parameter lmm
 * is unchanged (or changed from CV_ADAMS to CV_BDF) and the default
 * value for maxord is specified.
 *
 * All of the arguments to CVodeReInit have names and meanings
 * identical to those of CVodeInit.
 *
 * The return value of CVodeReInit is equal to CV_SUCCESS = 0 if
 * there were no errors; otherwise it is a negative int equal to:
 *   CV_MEM_NULL      indicating cvode_mem was NULL (i.e.,
 *                    CVodeCreate has not been called).
 *   CV_NO_MALLOC     indicating that cvode_mem has not been
 *                    allocated (i.e., CVodeInit has not been
 *                    called).
 *   CV_ILL_INPUT     indicating an input argument was illegal
 *                    (including an attempt to increase maxord).
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0);

/*
 * -----------------------------------------------------------------
 * Functions : CVodeSStolerances
 *             CVodeSVtolerances
 *             CVodeWFtolerances
 * -----------------------------------------------------------------
 *
 * These functions specify the integration tolerances. One of them
 * MUST be called before the first call to CVode.
 *
 * CVodeSStolerances specifies scalar relative and absolute tolerances.
 * CVodeSVtolerances specifies scalar relative tolerance and a vector
 *   absolute tolerance (a potentially different absolute tolerance 
 *   for each vector component).
 * CVodeWFtolerances specifies a user-provides function (of type CVEwtFn)
 *   which will be called to set the error weight vector.
 *
 * The tolerances reltol and abstol define a vector of error weights,
 * ewt, with components
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol)      (in the SS case), or
 *   ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (in the SV case).
 * This vector is used in all error and convergence tests, which
 * use a weighted RMS norm on all error-like vectors v:
 *    WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
 * where N is the problem dimension.
 *
 * The return value of these functions is equal to CV_SUCCESS = 0 if
 * there were no errors; otherwise it is a negative int equal to:
 *   CV_MEM_NULL      indicating cvode_mem was NULL (i.e.,
 *                    CVodeCreate has not been called).
 *   CV_NO_MALLOC     indicating that cvode_mem has not been
 *                    allocated (i.e., CVodeInit has not been
 *                    called).
 *   CV_ILL_INPUT     indicating an input argument was illegal
 *                    (e.g. a negative tolerance)
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol);
SUNDIALS_EXPORT int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol);
SUNDIALS_EXPORT int CVodeWFtolerances(void *cvode_mem, CVEwtFn efun);

/*
 * -----------------------------------------------------------------
 * Function : CVodeRootInit
 * -----------------------------------------------------------------
 * CVodeRootInit initializes a rootfinding problem to be solved
 * during the integration of the ODE system.  It must be called
 * after CVodeCreate, and before CVode.  The arguments are:
 *
 * cvode_mem = pointer to CVODE memory returned by CVodeCreate.
 *
 * nrtfn     = number of functions g_i, an int >= 0.
 *
 * g         = name of user-supplied function, of type CVRootFn,
 *             defining the functions g_i whose roots are sought.
 *
 * If a new problem is to be solved with a call to CVodeReInit,
 * where the new problem has no root functions but the prior one
 * did, then call CVodeRootInit with nrtfn = 0.
 *
 * The return value of CVodeRootInit is CV_SUCCESS = 0 if there were
 * no errors; otherwise it is a negative int equal to:
 *   CV_MEM_NULL     indicating cvode_mem was NULL, or
 *   CV_MEM_FAIL     indicating a memory allocation failed.
 *                    (including an attempt to increase maxord).
 *   CV_ILL_INPUT   indicating nrtfn > 0 but g = NULL.
 * In case of an error return, an error message is also printed.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g);

/*
 * -----------------------------------------------------------------
 * Function : CVode
 * -----------------------------------------------------------------
 * CVode integrates the ODE over an interval in t.
 * If itask is CV_NORMAL, then the solver integrates from its
 * current internal t value to a point at or beyond tout, then
 * interpolates to t = tout and returns y(tout) in the user-
 * allocated vector yout. If itask is CV_ONE_STEP, then the solver
 * takes one internal time step and returns in yout the value of
 * y at the new internal time. In this case, tout is used only
 * during the first call to CVode to determine the direction of
 * integration and the rough scale of the t variable. If tstop is
 * enabled (through a call to CVodeSetStopTime), then CVode returns
 * the solution at tstop. Once the integrator returns at a tstop
 * time, any future testing for tstop is disabled (and can be 
 * reenabled only though a new call to CVodeSetStopTime).
 * The time reached by the solver is placed in (*tret). The
 * user is responsible for allocating the memory for this value.
 *
 * cvode_mem is the pointer to CVODE memory returned by
 *           CVodeCreate.
 *
 * tout  is the next time at which a computed solution is desired.
 *
 * yout  is the computed solution vector. In CV_NORMAL mode with no
 *       errors and no roots found, yout=y(tout).
 *
 * tret  is a pointer to a real location. CVode sets (*tret) to
 *       the time reached by the solver and returns
 *       yout=y(*tret).
 *
 * itask is CV_NORMAL or CV_ONE_STEP. These two modes are described above.
 *
 * Here is a brief description of each return value:
 *
 * CV_SUCCESS:      CVode succeeded and no roots were found.
 *
 * CV_ROOT_RETURN:  CVode succeeded, and found one or more roots.
 *                  If nrtfn > 1, call CVodeGetRootInfo to see
 *                  which g_i were found to have a root at (*tret).
 *
 * CV_TSTOP_RETURN: CVode succeeded and returned at tstop.
 *
 * CV_MEM_NULL:     The cvode_mem argument was NULL.
 *
 * CV_NO_MALLOC:    cvode_mem was not allocated.
 *
 * CV_ILL_INPUT:    One of the inputs to CVode is illegal. This
 *                  includes the situation when a component of the
 *                  error weight vectors becomes < 0 during
 *                  internal time-stepping.  It also includes the
 *                  situation where a root of one of the root
 *                  functions was found both at t0 and very near t0.
 *                  The ILL_INPUT flag will also be returned if the
 *                  linear solver routine CV--- (called by the user
 *                  after calling CVodeCreate) failed to set one of
 *                  the linear solver-related fields in cvode_mem or
 *                  if the linear solver's init routine failed. In
 *                  any case, the user should see the printed
 *                  error message for more details.
 *
 * CV_TOO_MUCH_WORK: The solver took mxstep internal steps but
 *                  could not reach tout. The default value for
 *                  mxstep is MXSTEP_DEFAULT = 500.
 *
 * CV_TOO_MUCH_ACC: The solver could not satisfy the accuracy
 *                  demanded by the user for some internal step.
 *
 * CV_CONSTR_FAIL:  The inequality constraints were violated,
 *                  and the solver was unable to recover.
 *
 * CV_ERR_FAILURE:  Error test failures occurred too many times
 *                  (= MXNEF = 7) during one internal time step or
 *                  occurred with |h| = hmin.
 *
 * CV_CONV_FAILURE: Convergence test failures occurred too many
 *                  times (= MXNCF = 10) during one internal time
 *                  step or occurred with |h| = hmin.
 *
 * CV_LINIT_FAIL:   The linear solver's initialization function 
 *                  failed.
 *
 * CV_LSETUP_FAIL:  The linear solver's setup routine failed in an
 *                  unrecoverable manner.
 *
 * CV_LSOLVE_FAIL:  The linear solver's solve routine failed in an
 *                  unrecoverable manner.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVode(void *cvode_mem, realtype tout, N_Vector yout,
			  realtype *tret, int itask);

/*
 * -----------------------------------------------------------------
 * Function : CVodeGetDky
 * -----------------------------------------------------------------
 * CVodeGetDky computes the kth derivative of the y function at
 * time t, where tn-hu <= t <= tn, tn denotes the current
 * internal time reached, and hu is the last internal step size
 * successfully used by the solver. The user may request
 * k=0, 1, ..., qu, where qu is the order last used. The
 * derivative vector is returned in dky. This vector must be
 * allocated by the caller. It is only legal to call this
 * function after a successful return from CVode.
 *
 * cvode_mem is the pointer to CVODE memory returned by
 *           CVodeCreate.
 *
 * t   is the time at which the kth derivative of y is evaluated.
 *     The legal range for t is [tn-hu,tn] as described above.
 *
 * k   is the order of the derivative of y to be computed. The
 *     legal range for k is [0,qu] as described above.
 *
 * dky is the output derivative vector [((d/dy)^k)y](t).
 *
 * The return value for CVodeGetDky is one of:
 *
 *   CV_SUCCESS:  CVodeGetDky succeeded.
 *
 *   CV_BAD_K:    k is not in the range 0, 1, ..., qu.
 *
 *   CV_BAD_T:    t is not in the interval [tn-hu,tn].
 *
 *   CV_BAD_DKY:  The dky argument was NULL.
 *
 *   CV_MEM_NULL: The cvode_mem argument was NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky);

/*
 * -----------------------------------------------------------------
 * Integrator optional output extraction functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the main integrator.
 * -----------------------------------------------------------------
 * CVodeGetWorkSpace returns the CVODE real and integer workspaces
 * CVodeGetNumSteps returns the cumulative number of internal
 *                  steps taken by the solver
 * CVodeGetNumRhsEvals returns the number of calls to the user's
 *                     f function
 * CVodeGetNumLinSolvSetups returns the number of calls made to
 *                          the linear solver's setup routine
 * CVodeGetNumErrTestFails returns the number of local error test
 *                         failures that have occured
 * CVodeGetLastOrder returns the order used during the last
 *                   internal step
 * CVodeGetCurrentOrder returns the order to be used on the next
 *                      internal step
 * CVodeGetNumStabLimOrderReds returns the number of order
 *                             reductions due to stability limit
 *                             detection
 * CVodeGetActualInitStep returns the actual initial step size
 *                        used by CVODE
 * CVodeGetLastStep returns the step size for the last internal
 *                  step
 * CVodeGetCurrentStep returns the step size to be attempted on
 *                     the next internal step
 * CVodeGetCurrentTime returns the current internal time reached
 *                     by the solver
 * CVodeGetTolScaleFactor returns a suggested factor by which the
 *                        user's tolerances should be scaled when
 *                        too much accuracy has been requested for
 *                        some internal step
 * CVodeGetErrWeights returns the current error weight vector.
 *                    The user must allocate space for eweight.
 * CVodeGetEstLocalErrors returns the vector of estimated local
 *                        errors. The user must allocate space
 *                        for ele.
 * CVodeGetNumGEvals returns the number of calls to the user's
 *                   g function (for rootfinding)
 * CVodeGetRootInfo returns the indices for which g_i was found to 
 *                  have a root. The user must allocate space for 
 *                  rootsfound. For i = 0 ... nrtfn-1, 
 *                  rootsfound[i] = 1 if g_i has a root, and = 0 if not.
 *
 * CVodeGet* return values:
 *   CV_SUCCESS   if succesful
 *   CV_MEM_NULL  if the cvode memory was NULL
 *   CV_NO_SLDET  if stability limit was not turned on
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeGetWorkSpace(void *cvode_mem, long int *lenrw, long int *leniw);
SUNDIALS_EXPORT int CVodeGetNumSteps(void *cvode_mem, long int *nsteps);
SUNDIALS_EXPORT int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevals);
SUNDIALS_EXPORT int CVodeGetNumLinSolvSetups(void *cvode_mem, long int *nlinsetups);
SUNDIALS_EXPORT int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails);
SUNDIALS_EXPORT int CVodeGetLastOrder(void *cvode_mem, int *qlast);
SUNDIALS_EXPORT int CVodeGetCurrentOrder(void *cvode_mem, int *qcur);
SUNDIALS_EXPORT int CVodeGetNumStabLimOrderReds(void *cvode_mem, long int *nslred);
SUNDIALS_EXPORT int CVodeGetActualInitStep(void *cvode_mem, realtype *hinused);
SUNDIALS_EXPORT int CVodeGetLastStep(void *cvode_mem, realtype *hlast);
SUNDIALS_EXPORT int CVodeGetCurrentStep(void *cvode_mem, realtype *hcur);
SUNDIALS_EXPORT int CVodeGetCurrentTime(void *cvode_mem, realtype *tcur);
SUNDIALS_EXPORT int CVodeGetTolScaleFactor(void *cvode_mem, realtype *tolsfac);
SUNDIALS_EXPORT int CVodeGetErrWeights(void *cvode_mem, N_Vector eweight);
SUNDIALS_EXPORT int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele);
SUNDIALS_EXPORT int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals);
SUNDIALS_EXPORT int CVodeGetRootInfo(void *cvode_mem, int *rootsfound);

/*
 * -----------------------------------------------------------------
 * As a convenience, the following functions provides the
 * optional outputs in one group.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeGetIntegratorStats(void *cvode_mem, long int *nsteps,
					    long int *nfevals, long int *nlinsetups,
					    long int *netfails, int *qlast,
					    int *qcur, realtype *hinused, realtype *hlast,
					    realtype *hcur, realtype *tcur);

/*
 * -----------------------------------------------------------------
 * Nonlinear solver optional output extraction functions
 * -----------------------------------------------------------------
 * The following functions can be called to get optional outputs
 * and statistics related to the nonlinear solver.
 * -----------------------------------------------------------------
 * CVodeGetNumNonlinSolvIters returns the number of nonlinear
 *                            solver iterations performed.
 * CVodeGetNumNonlinSolvConvFails returns the number of nonlinear
 *                                convergence failures.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters);
SUNDIALS_EXPORT int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails);

/*
 * -----------------------------------------------------------------
 * As a convenience, the following function provides the
 * nonlinear solver optional outputs in a group.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters,
					    long int *nncfails);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a CVODE return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *CVodeGetReturnFlagName(long int flag);

/*
 * -----------------------------------------------------------------
 * Function : CVodeFree
 * -----------------------------------------------------------------
 * CVodeFree frees the problem memory cvode_mem allocated by
 * CVodeCreate and CVodeInit. Its only argument is the pointer
 * cvode_mem returned by CVodeCreate.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void CVodeFree(void **cvode_mem);

/*
 * -----------------------------------------------------------------
 * Function : CVodeSetNonlinearSolver
 * -----------------------------------------------------------------
 * CVodeNonlinearSolver attaches a nonlinear solver object to the
 * CVODE memory.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVodeSetNonlinearSolver(void *cvode_mem, SUNNonlinearSolver NLS);

#ifdef __cplusplus
}
#endif

#endif
