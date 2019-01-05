/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Southern Methodist University and Lawrence Livermore
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * Header file for the Explicit Runge Kutta time step module for ARKode.
 *--------------------------------------------------------------*/

#ifndef _ERKSTEP_H
#define _ERKSTEP_H

#include <sundials/sundials_nvector.h>
#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ERKSTEP Constants
  ===============================================================*/

/* Default Butcher tables for each order */
#define DEFAULT_ERK_2           HEUN_EULER_2_1_2
#define DEFAULT_ERK_3           BOGACKI_SHAMPINE_4_2_3
#define DEFAULT_ERK_4           ZONNEVELD_5_3_4
#define DEFAULT_ERK_5           CASH_KARP_6_4_5
#define DEFAULT_ERK_6           VERNER_8_5_6
#define DEFAULT_ERK_8           FEHLBERG_13_7_8


/*===============================================================
  ERKSTEP Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  ERKStepCreate

  This creates an internal memory block for a problem to be
  solved by ARKode, using the ERK time step module.  If successful,
  it returns a pointer to initialized problem memory. This
  pointer should be passed to all ERKStep-related routines.
  If an initialization error occurs, this routine will print an
  error message to standard err and return NULL.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void* ERKStepCreate(ARKRhsFn f, realtype t0,
                                    N_Vector y0);


/*---------------------------------------------------------------
  ERKStepResize

  This re-initializes the ERKStep memory for a problem with a
  changing vector size.  It is assumed that the problem dynamics
  before and after the vector resize will be comparable, so that
  all time-stepping heuristics prior to calling ERKStepResize
  remain valid after the call.  If instead the dynamics should be
  re-calibrated, the ERKStep memory structure should be deleted
  with a call to ERKStepFree, and re-created with a call to
  ERKStepCreate.

  To aid in the vector-resize operation, the user can supply a
  vector resize function, that will take as input an N_Vector with
  the previous size, and return as output a corresponding vector
  of the new size.  If this function (of type ARKVecResizeFn) is
  not supplied (i.e. is set to NULL), then all existing N_Vectors
  will be destroyed and re-cloned from the input vector.

  In the case that the dynamical time scale should be modified
  slightly from the previous time scale, an input "hscale" is
  allowed, that will re-scale the upcoming time step by the
  specified factor.  If a value <= 0 is specified, the default of
  1.0 will be used.

  Other arguments:
    arkode_mem       Existing ERKStep memory data structure.
    ynew             The newly-sized solution vector, holding
                     the current dependent variable values.
    t0               The current value of the independent
                     variable.
    resize_data      User-supplied data structure that will be
                     passed to the supplied resize function.

  The return value of ERKStepResize is equal to ARK_SUCCESS = 0 if
  there were no errors; otherwise it is a negative int equal to:
    ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
                     ERKStepCreate has not been called).
    ARK_NO_MALLOC    indicating that arkode_mem has not been
                     allocated.
    ARK_ILL_INPUT    indicating an input argument was illegal
                     (including an error from the supplied
                     resize function).
  In case of an error return, an error message is also printed.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepResize(void *arkode_mem, N_Vector ynew,
                                  realtype hscale, realtype t0,
                                  ARKVecResizeFn resize,
                                  void *resize_data);


/*---------------------------------------------------------------
  ERKStepReInit

  This re-initializes the ERK time step module.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepReInit(void* arkode_mem, ARKRhsFn f,
                                  realtype t0, N_Vector y0);


/*---------------------------------------------------------------
  ERKStepSStolerances, ERKStepSVtolerances, ERKStepWFtolerances

  These specify the integration tolerances. One of them SHOULD be
  called before the first call to ERKStepEvolve; otherwise
  default values of reltol=1e-4 and abstol=1e-9 will be used,
  which may be entirely incorrect for a specific problem.

  ERKStepSStolerances specifies scalar relative and absolute
    tolerances.
  ERKStepSVtolerances specifies scalar relative tolerance and a
    vector absolute tolerance (a potentially different absolute
    tolerance for each vector component).
  ERKStepWFtolerances specifies a user-provided function (of type
    ARKEwtFn) which will be called to set the error weight vector.

  The tolerances reltol and abstol define a vector of error
  weights, ewt, with components
    ewt[i] = 1/(reltol*abs(y[i]) + abstol)      (in SS case), or
    ewt[i] = 1/(reltol*abs(y[i]) + abstol[i])   (in SV case).
  This vector is used in all error and convergence tests, which
  use a weighted RMS norm on all error-like vectors v:
     WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*ewt[i])^2 ),
  where N is the problem dimension.

  The return value of these functions is equal to ARK_SUCCESS=0
  if there were no errors; otherwise it is a negative int equal
  to:
    ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
                     ERKStepCreate has not been called).
    ARK_NO_MALLOC    indicating that arkode_mem has not been
                     allocated.
    ARK_ILL_INPUT    indicating an input argument was illegal
                     (e.g. a negative tolerance)
  In case of an error return, an error message is also printed.
  --------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepSStolerances(void *arkode_mem,
                                        realtype reltol,
                                        realtype abstol);
SUNDIALS_EXPORT int ERKStepSVtolerances(void *arkode_mem,
                                        realtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int ERKStepWFtolerances(void *arkode_mem,
                                        ARKEwtFn efun);


/*---------------------------------------------------------------
  ERKStepRootInit

  This initializes a rootfinding problem to be solved during the
  integration of the ODE system.  It must be called after
  ERKStepCreate, and before ERKStepEvolve.  The arguments are:

  arkode_mem = pointer to ARKode memory returned by ERKStepCreate.

  nrtfn      = number of functions g_i, an integer >= 0.

  g          = name of user-supplied function, of type ARKRootFn,
               defining the functions g_i whose roots are sought.

  If a new problem is to be solved with a call to ERKStepReInit,
  where the new problem has no root functions but the prior one
  did, then call ERKStepRootInit again with nrtfn = 0.

  The return value is ARK_SUCCESS = 0 if there
  were no errors; otherwise it is a negative int equal to:
    ARK_MEM_NULL    indicating arkode_mem was NULL, or
    ARK_MEM_FAIL    indicating a memory allocation failed.
                    (including an attempt to increase maxord).
    ARK_ILL_INPUT   indicating nrtfn > 0 but g = NULL.
  In case of an error return, an error message is also printed.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepRootInit(void *arkode_mem, int nrtfn,
                                    ARKRootFn g);


/*---------------------------------------------------------------
  ERKStep optional input specification functions -- ALL of these
  must be called AFTER ERKStepCreate.
  -----------------------------------------------------------------
  The following functions can be called to set optional inputs
  to values other than the defaults given below.

  Function                   |  Optional input / [ default value ]
  -----------------------------------------------------------------
  ERKStepSetDefaults         | resets all optional inputs to ERKStep
                             | default values.  Does not change
                             | problem-defining function pointers or
                             | user_data pointer.  Also leaves alone
                             | any data structures/options related
                             | to root-finding (those can be reset
                             | using ARKodeRootInit).
                             | [internal]
                             |
  ERKStepSetOrder            | method order to be used by the solver.
                             | [4]
                             |
  ERKStepSetDenseOrder       | polynomial order to be used for dense
                             | output.  Allowed values are between 0
                             | and min(q,5) (where q is the order of
                             | the integrator)
                             | [3]
                             |
  ERKStepSetErrHandlerFn     | user-provided ErrHandler function.
                             | [internal]
                             |
  ERKStepSetErrFile          | the file pointer for an error file
                             | where all ERKStep warning and error
                             | messages will be written if the
                             | default internal error handling
                             | function is used. This parameter can
                             | be stdout (standard output), stderr
                             | (standard error), or a file pointer
                             | (corresponding to a user error file
                             | opened for writing) returned by fopen.
                             | If not called, then all messages will
                             | be written to stderr.
                             | [stderr]
                             |
  ERKStepSetUserData         | a pointer to user data that will be
                             | passed to the user's f function every
                             | time f is called.
                             | [NULL]
                             |
  ERKStepSetDiagnostics      | the file pointer for a diagnostics file
                             | where all ERKStep adaptivity and solver
                             | information is written.  This parameter can
                             | be stdout or stderr, though the preferred
                             | approach is to specify a file pointer
                             | (corresponding to a user diagnostics file
                             | opened for writing) returned by fopen.  If
                             | not called, or if called with a NULL file
                             | pointer, all diagnostics output is disabled.
                             | NOTE: when run in parallel, only one process
                             | should set a non-NULL value for this pointer,
                             | since statistics from all processes would be
                             | identical.
                             | [NULL]
                             |
  ERKStepSetTable            | specifies to use a customized Butcher
                             | table for the explicit portion of the
                             | system.  This automatically calls
                             | ERKStepSetExplicit
                             | [determined by ARKode based on order]
                             |
  ERKStepSetTableNum         | specifies to use a built-in Butcher
                             | table for the explicit portion of the
                             | system.  The integer argument should
                             | match an existing method in
                             | ARKodeButcherTable_LoadERK() within the file
                             | arkode_butcher.c.  Error-checking is
                             | performed to ensure that the table
                             | exists, and is not implicit.  This
                             | automatically calls ERKStepSetExplicit
                             | [determined by ARKode based on order]
                             |
  ERKStepSetMaxNumSteps      | maximum number of internal steps to be
                             | taken by the solver in its attempt to
                             | reach tout.
                             | [500]
                             |
  ERKStepSetMaxHnilWarns     | maximum number of warning messages
                             | issued by the solver that t+h==t on
                             | the next internal step. A value of -1
                             | means no such messages are issued.
                             | [10]
                             |
  ERKStepSetInitStep         | initial step size.
                             | [estimated internally]
                             |
  ERKStepSetMinStep          | minimum absolute value of step size
                             | allowed.
                             | [0.0]
                             |
  ERKStepSetMaxStep          | maximum absolute value of step size
                             | allowed.
                             | [infinity]
                             |
  ERKStepSetStopTime         | the independent variable value past
                             | which the solution is not to proceed.
                             | [infinity]
                             |
  ERKStepSetFixedStep        | specifies to use a fixed step size
                             | throughout integration
                             | [off]
                             |
  ERKStepSetCFLFraction      | safety factor to use for explicitly
                             | stable steps
                             | [0.5]
                             |
  ERKStepSetSafetyFactor     | safety factor to use for error-based
                             | step adaptivity
                             | [0.96]
                             |
  ERKStepSetErrorBias        | error bias factor to use in error-based
                             | step adaptivity
                             | [1.5]
                             |
  ERKStepSetMaxGrowth        | maximum growth factor for successive
                             | time steps (not including the first step).
                             | [20.0]
                             |
  ERKStepSetMaxFirstGrowth   | maximum growth factor for first step.
                             | [10000.0]
                             |
  ERKStepSetMaxEFailGrowth   | maximum growth factor after an error failure.
                             | [0.3]
                             |
  ERKStepSetSmallNumEFails   | maximum number of error failures before
                             | MaxFailGrowth factor is used.
                             | [2]
                             |
  ERKStepSetFixedStepBounds  | step growth interval to force retention of
                             | the same step size
                             | [1.0 1.5]
                             |
  ERKStepSetAdaptivityMethod | Method to use for time step adaptivity
                             | [0]
                             |
  ERKStepSetAdaptivityFn     | user-provided time step adaptivity
                             | function.
                             | [internal]
                             |
  ERKStepSetStabilityFn      | user-provided explicit time step
                             | stability function.
                             | [internal]
                             |
  ERKStepSetMaxErrTestFails  | Maximum number of error test failures
                             | in attempting one step.
                             | [7]
  -----------------------------------------------------------------
  ERKStepSetRootDirection      | Specifies the direction of zero
                               | crossings to be monitored
                               | [both directions]
                               |
  ERKStepSetNoInactiveRootWarn | disable warning about possible
                               | g==0 at beginning of integration
  -----------------------------------------------------------------
  Return flag:
    ARK_SUCCESS   if successful
    ARK_MEM_NULL  if the arkode memory is NULL
    ARK_ILL_INPUT if an argument has an illegal value
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int ERKStepSetOrder(void *arkode_mem, int maxord);
SUNDIALS_EXPORT int ERKStepSetDenseOrder(void *arkode_mem, int dord);
SUNDIALS_EXPORT int ERKStepSetTable(void *arkode_mem,
                                    ARKodeButcherTable B);
SUNDIALS_EXPORT int ERKStepSetTableNum(void *arkode_mem, int itable);
SUNDIALS_EXPORT int ERKStepSetCFLFraction(void *arkode_mem,
                                          realtype cfl_frac);
SUNDIALS_EXPORT int ERKStepSetSafetyFactor(void *arkode_mem,
                                           realtype safety);
SUNDIALS_EXPORT int ERKStepSetErrorBias(void *arkode_mem,
                                        realtype bias);
SUNDIALS_EXPORT int ERKStepSetMaxGrowth(void *arkode_mem,
                                        realtype mx_growth);
SUNDIALS_EXPORT int ERKStepSetFixedStepBounds(void *arkode_mem,
                                              realtype lb, realtype ub);
SUNDIALS_EXPORT int ERKStepSetAdaptivityMethod(void *arkode_mem,
                                               int imethod,
                                               int idefault, int pq,
                                               realtype *adapt_params);
SUNDIALS_EXPORT int ERKStepSetAdaptivityFn(void *arkode_mem,
                                           ARKAdaptFn hfun,
                                           void *h_data);
SUNDIALS_EXPORT int ERKStepSetMaxFirstGrowth(void *arkode_mem,
                                             realtype etamx1);
SUNDIALS_EXPORT int ERKStepSetMaxEFailGrowth(void *arkode_mem,
                                             realtype etamxf);
SUNDIALS_EXPORT int ERKStepSetSmallNumEFails(void *arkode_mem,
                                             int small_nef);
SUNDIALS_EXPORT int ERKStepSetStabilityFn(void *arkode_mem,
                                          ARKExpStabFn EStab,
                                          void *estab_data);
SUNDIALS_EXPORT int ERKStepSetMaxErrTestFails(void *arkode_mem,
                                              int maxnef);
SUNDIALS_EXPORT int ERKStepSetMaxNumSteps(void *arkode_mem,
                                          long int mxsteps);
SUNDIALS_EXPORT int ERKStepSetMaxHnilWarns(void *arkode_mem,
                                           int mxhnil);
SUNDIALS_EXPORT int ERKStepSetInitStep(void *arkode_mem,
                                       realtype hin);
SUNDIALS_EXPORT int ERKStepSetMinStep(void *arkode_mem,
                                      realtype hmin);
SUNDIALS_EXPORT int ERKStepSetMaxStep(void *arkode_mem,
                                      realtype hmax);
SUNDIALS_EXPORT int ERKStepSetStopTime(void *arkode_mem,
                                       realtype tstop);
SUNDIALS_EXPORT int ERKStepSetFixedStep(void *arkode_mem,
                                        realtype hfixed);

SUNDIALS_EXPORT int ERKStepSetRootDirection(void *arkode_mem,
                                            int *rootdir);
SUNDIALS_EXPORT int ERKStepSetNoInactiveRootWarn(void *arkode_mem);

SUNDIALS_EXPORT int ERKStepSetErrHandlerFn(void *arkode_mem,
                                           ARKErrHandlerFn ehfun,
                                           void *eh_data);
SUNDIALS_EXPORT int ERKStepSetErrFile(void *arkode_mem,
                                      FILE *errfp);
SUNDIALS_EXPORT int ERKStepSetUserData(void *arkode_mem,
                                       void *user_data);
SUNDIALS_EXPORT int ERKStepSetDiagnostics(void *arkode_mem,
                                          FILE *diagfp);

SUNDIALS_EXPORT int ERKStepSetPostprocessStepFn(void *arkode_mem,
                                                ARKPostProcessStepFn ProcessStep);


/*---------------------------------------------------------------
  ERKStepEvolve

  This integrates the ODE over an interval in t.

  ERKStepEvolve may be run in one of two modes (ARK_NORMAL or
  ARK_ONE_STEP), as determined by the itask argument:

  If itask is ARK_NORMAL, then the solver integrates from its
  current internal t value to a point at or beyond tout, then
  interpolates to t = tout and returns y(tout) in the user-
  allocated vector yout.  This interpolation is typically less
  accurate than the full time step solutions produced by the
  solver, since the interpolating polynomial relies on the
  internal stage solutions, that may have reduced accuracy in
  comparison with the full time step solutions.  If the user
  wishes that this returned value have full method accuracy, they
  may issue a call to ERKStepSetStopTime before the call to
  ERKStepEvolve to specify a fixed stop time to end the time step
  and return to the user.  Once the integrator returns at a tstop
  time, any future testing for tstop is disabled (and can be
  reenabled only though a new call to ERKStepSetStopTime).

  If itask is ARK_ONE_STEP, then the solver takes one internal
  time step and returns in yout the value of y at the new internal
  time. In this case, tout is used only during the first call to
  ERKStepEvolve to determine the direction of integration and the
  rough scale of the t variable.  As with the ARK_NORMAL mode, a
  user may specify a specific stop time for output of this step,
  assuming that the requested step is smaller than the step taken
  by the method.

  The time reached by the solver is placed in (*tret). The
  user is responsible for allocating the memory for this value.

  arkode_mem is the pointer to ERKStep memory returned by
             ERKStepCreate.

  tout  is the next time at which a computed solution is desired.

  yout  is the computed solution vector. In ARK_NORMAL mode with no
        errors and no roots found, yout=y(tout).

  tret  is a pointer to a real location. ERKStepEvolve sets (*tret)
        to the time reached by the solver and returns yout=y(*tret).

  itask is ARK_NORMAL or ARK_ONE_STEP, as described above.

  Here is a brief description of each return value:

  ARK_SUCCESS:      ERKStepEvolve succeeded and no roots were found.

  ARK_ROOT_RETURN:  ERKStepEvolve succeeded, and found one or more roots.
                    If nrtfn > 1, call ERKStepGetRootInfo to see
                    which g_i were found to have a root at (*tret).

  ARK_TSTOP_RETURN: ERKStepEvolve succeeded and returned at tstop.

  ARK_MEM_NULL:     The arkode_mem argument was NULL.

  ARK_NO_MALLOC:    arkode_mem was not allocated.

  ARK_ILL_INPUT:    One of the inputs is illegal. This
                    includes the situation when a component of the
                    error weight vectors becomes < 0 during
                    internal time-stepping.  It also includes the
                    situation where a root of one of the root
                    functions was found both at t0 and very near t0.
                    The ILL_INPUT flag will also be returned if the
                    linear solver routine ARK--- (called by the user
                    after calling ERKStepCreate) failed to set one of
                    the linear solver-related fields in arkode_mem or
                    if the linear solver's init routine failed. In
                    any case, the user should see the printed
                    error message for more details.

  ARK_TOO_MUCH_WORK: The solver took mxstep internal steps but
                    could not reach tout prior to the maximum number
                    of steps (ark_mxstep).

  ARK_TOO_MUCH_ACC: The solver could not satisfy the accuracy
                    demanded by the user for some internal step.

  ARK_ERR_FAILURE:  Error test failures occurred too many times
                    (= ark_maxnef) during one internal time step
                    or occurred with |h| = hmin.

  ARK_CONV_FAILURE: Convergence test failures occurred too many
                    times (= ark_maxncf) during one internal time
                    step or occurred with |h| = hmin.

  ARK_LINIT_FAIL:   The linear solver's initialization function
                    failed.

  ARK_LSETUP_FAIL:  The linear solver's setup routine failed in an
                    unrecoverable manner.

  ARK_LSOLVE_FAIL:  The linear solver's solve routine failed in an
                    unrecoverable manner.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepEvolve(void *arkode_mem, realtype tout,
                                  N_Vector yout, realtype *tret,
                                  int itask);


/*---------------------------------------------------------------
  ERKStepGetDky

  This computes the kth derivative of the y function at time t,
  where tn-hu <= t <= tn, tn denotes the current internal time
  reached, and hu is the last internal step size successfully
  used by the solver. The user may request k=0, 1, ..., d, where
  d = min(5,q), with q the order of accuracy for the time
  integration method. The derivative vector is returned in dky.
  This vector must be allocated by the caller. It is only legal
  to call this function after a successful return from
  ERKStepEvolve.

  arkode_mem is the pointer to ERKStep memory returned by
             ERKStepCreate.

  t   is the time at which the kth derivative of y is evaluated.
      The legal range for t is [tn-hu,tn] as described above.

  k   is the order of the derivative of y to be computed. The
      legal range for k is [0,min(q,3)] as described above.

  dky is the output derivative vector [((d/dy)^k)y](t).

  The return value for ERKStepGetDky is one of:

    ARK_SUCCESS:  ERKStepGetDky succeeded.

    ARK_BAD_K:    k is not in the range 0, 1, ..., s-1.

    ARK_BAD_T:    t is not in the interval [tn-hu,tn].

    ARK_BAD_DKY:  The dky argument was NULL.

    ARK_MEM_NULL: The arkode_mem argument was NULL.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepGetDky(void *arkode_mem, realtype t,
                                  int k, N_Vector dky);


/*---------------------------------------------------------------
  Optional outputs from the ERKStep module:  the following
  functions can be called to get optional outputs and statistics
  related to the main integrator.

  ERKStepGetNumExpSteps returns the cumulative number of stability
                        limited steps taken by the solver

  ERKStepGetNumAccSteps returns the cumulative number of accuracy
                        limited steps taken by the solver

  ERKStepGetNumStepAttempts returns the total number of steps
                            attempted by the solver

  ERKStepGetNumRhsEvals returns the number of calls to the user's
                        f functions

  ERKStepGetNumErrTestFails returns the number of local error test
                            failures that have occured

  ERKStepGetCurrentButcherTable returns the Butcher table
                                currently in use

  ERKStepGetEstLocalErrors returns the vector of estimated local
                           errors. The user must allocate space
                           for ele.

  ERKStepGetNumSteps returns the cumulative number of internal
                     steps taken by the solver

  ERKStepGetActualInitStep returns the actual initial step size
                           used by ERKStep

  ERKStepGetLastStep returns the step size for the last internal
                     step

  ERKStepGetCurrentStep returns the step size to be attempted on
                        the next internal step

  ERKStepGetCurrentTime returns the current internal time reached
                        by the solver

  ERKStepGetTolScaleFactor returns a suggested factor by which the
                           user's tolerances should be scaled when
                           too much accuracy has been requested for
                           some internal step

  ERKStepGetErrWeights returns the current error weight vector.
                       The user must allocate space for eweight.

  ERKStepGetWorkSpace returns the ERKStep real and integer workspaces

  ERKStepGetNumGEvals returns the number of calls to the user's
                      g function (for rootfinding)

  ERKStepGetRootInfo returns the indices for which g_i was found to
                     have a root. The user must allocate space for
                     rootsfound. For i = 0 ... nrtfn-1,
                     rootsfound[i] = 1 if g_i has a root, and = 0
                     if not.

  The return value of ERKStepGet* is one of:
     ARK_SUCCESS   if successful
     ARK_MEM_NULL  if the ERKStep memory structure was NULL
     ARK_LMEM_NULL if a linear solver memory structure was NULL
  ---------------------------------------------------------------*/

SUNDIALS_EXPORT int ERKStepGetNumExpSteps(void *arkode_mem,
                                          long int *expsteps);
SUNDIALS_EXPORT int ERKStepGetNumAccSteps(void *arkode_mem,
                                          long int *accsteps);
SUNDIALS_EXPORT int ERKStepGetNumStepAttempts(void *arkode_mem,
                                              long int *step_attempts);
SUNDIALS_EXPORT int ERKStepGetNumRhsEvals(void *arkode_mem,
                                          long int *nfevals);
SUNDIALS_EXPORT int ERKStepGetNumErrTestFails(void *arkode_mem,
                                              long int *netfails);
SUNDIALS_EXPORT int ERKStepGetCurrentButcherTable(void *arkode_mem,
                                                  ARKodeButcherTable *B);
SUNDIALS_EXPORT int ERKStepGetEstLocalErrors(void *arkode_mem,
                                             N_Vector ele);
SUNDIALS_EXPORT int ERKStepGetWorkSpace(void *arkode_mem,
                                        long int *lenrw,
                                        long int *leniw);
SUNDIALS_EXPORT int ERKStepGetNumSteps(void *arkode_mem,
                                       long int *nsteps);
SUNDIALS_EXPORT int ERKStepGetActualInitStep(void *arkode_mem,
                                             realtype *hinused);
SUNDIALS_EXPORT int ERKStepGetLastStep(void *arkode_mem,
                                       realtype *hlast);
SUNDIALS_EXPORT int ERKStepGetCurrentStep(void *arkode_mem,
                                          realtype *hcur);
SUNDIALS_EXPORT int ERKStepGetCurrentTime(void *arkode_mem,
                                          realtype *tcur);
SUNDIALS_EXPORT int ERKStepGetTolScaleFactor(void *arkode_mem,
                                             realtype *tolsfac);
SUNDIALS_EXPORT int ERKStepGetErrWeights(void *arkode_mem,
                                         N_Vector eweight);
SUNDIALS_EXPORT int ERKStepGetNumGEvals(void *arkode_mem,
                                        long int *ngevals);
SUNDIALS_EXPORT int ERKStepGetRootInfo(void *arkode_mem,
                                       int *rootsfound);

/*---------------------------------------------------------------
  As a convenience, the following functions provide the
  optional outputs in grouped form.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepGetTimestepperStats(void *arkode_mem,
                                               long int *expsteps,
                                               long int *accsteps,
                                               long int *step_attempts,
                                               long int *nfevals,
                                               long int *netfails);
SUNDIALS_EXPORT int ERKStepGetStepStats(void *arkode_mem,
                                        long int *nsteps,
                                        realtype *hinused,
                                        realtype *hlast,
                                        realtype *hcur,
                                        realtype *tcur);

/*---------------------------------------------------------------
  The following function returns the name of the constant
  associated with a ERKStep return flag
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT char *ERKStepGetReturnFlagName(long int flag);

/*---------------------------------------------------------------
  Function : ERKStepWriteParameters
  -----------------------------------------------------------------
  ERKStepWriteParameters outputs all timestepper module parameters
  to the provided file pointer.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepWriteParameters(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
  Function : ERKStepWriteButcher
  -----------------------------------------------------------------
  ERKStepWriteButcher outputs the Butcher table to the
  provided file pointer.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ERKStepWriteButcher(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
  ERKStepFree

  This frees the problem memory arkode_mem allocated by
  ERKStepCreate. Its only argument is the pointer arkode_mem
  returned by ERKStepCreate.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void ERKStepFree(void **arkode_mem);

/*---------------------------------------------------------------
  ERKStepPrintMem

  This routine outputs the memory from the ERKStep structure and
  the main ARKode infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void ERKStepPrintMem(void* arkode_mem, FILE* outfile);


#ifdef __cplusplus
}
#endif

#endif
