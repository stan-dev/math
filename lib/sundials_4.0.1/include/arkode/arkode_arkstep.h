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
 * Header file for the Additive Runge Kutta time step module for ARKode.
 *--------------------------------------------------------------*/

#ifndef _ARKSTEP_H
#define _ARKSTEP_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <arkode/arkode.h>
#include <arkode/arkode_ls.h>
#include <arkode/arkode_butcher_erk.h>
#include <arkode/arkode_butcher_dirk.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARKSTEP Constants
  ===============================================================*/

/* Default Butcher tables for each method/order */

/*    explicit */
#define DEFAULT_ERK_2           HEUN_EULER_2_1_2
#define DEFAULT_ERK_3           BOGACKI_SHAMPINE_4_2_3
#define DEFAULT_ERK_4           ZONNEVELD_5_3_4
#define DEFAULT_ERK_5           CASH_KARP_6_4_5
#define DEFAULT_ERK_6           VERNER_8_5_6
#define DEFAULT_ERK_8           FEHLBERG_13_7_8

/*    implicit */
#define DEFAULT_DIRK_2          SDIRK_2_1_2
#define DEFAULT_DIRK_3          ARK324L2SA_DIRK_4_2_3
#define DEFAULT_DIRK_4          SDIRK_5_3_4
#define DEFAULT_DIRK_5          ARK548L2SA_DIRK_8_4_5

/*    ImEx */
#define DEFAULT_ARK_ETABLE_3    ARK324L2SA_ERK_4_2_3
#define DEFAULT_ARK_ETABLE_4    ARK436L2SA_ERK_6_3_4
#define DEFAULT_ARK_ETABLE_5    ARK548L2SA_ERK_8_4_5
#define DEFAULT_ARK_ITABLE_3    ARK324L2SA_DIRK_4_2_3
#define DEFAULT_ARK_ITABLE_4    ARK436L2SA_DIRK_6_3_4
#define DEFAULT_ARK_ITABLE_5    ARK548L2SA_DIRK_8_4_5


/*===============================================================
  ARKSTEP Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  ARKStepCreate

  This creates an internal memory block for a problem to be
  solved by ARKode, using the ARK time step module.  If successful,
  it returns a pointer to initialized problem memory. This
  pointer should be passed to all ARKStep-related routines.
  If an initialization error occurs, this routine will print an
  error message to standard err and return NULL.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void* ARKStepCreate(ARKRhsFn fe, ARKRhsFn fi,
                                    realtype t0, N_Vector y0);



/*---------------------------------------------------------------
  ARKStepResize

  This re-initializes the ARKStep memory for a problem with a
  changing vector size.  It is assumed that the problem dynamics
  before and after the vector resize will be comparable, so that
  all time-stepping heuristics prior to calling ARKStepResize
  remain valid after the call.  If instead the dynamics should be
  re-calibrated, the ARKStep memory structure should be deleted
  with a call to ARKStepFree, and re-created with a call to
  ARKStepCreate.

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
    arkode_mem       Existing ARKStep memory data structure.
    ynew             The newly-sized solution vector, holding
                     the current dependent variable values.
    t0               The current value of the independent
                     variable.
    resize_data      User-supplied data structure that will be
                     passed to the supplied resize function.

  The return value of ARKStepResize is equal to ARK_SUCCESS = 0 if
  there were no errors; otherwise it is a negative int equal to:
    ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
                     ARKStepCreate has not been called).
    ARK_NO_MALLOC    indicating that arkode_mem has not been
                     allocated.
    ARK_ILL_INPUT    indicating an input argument was illegal
                     (including an error from the supplied
                     resize function).
  In case of an error return, an error message is also printed.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepResize(void *arkode_mem, N_Vector ynew,
                                  realtype hscale, realtype t0,
                                  ARKVecResizeFn resize,
                                  void *resize_data);


/*---------------------------------------------------------------
  ARKStepReInit

  This re-initializes the ARK time step module.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepReInit(void* arkode_mem, ARKRhsFn fe,
                                  ARKRhsFn fi, realtype t0, N_Vector y0);


/*---------------------------------------------------------------
  ARKStepSStolerances, ARKStepSVtolerances, ARKStepWFtolerances

  These specify the integration tolerances. One of them SHOULD be
  called before the first call to ARKStepEvolve; otherwise
  default values of reltol=1e-4 and abstol=1e-9 will be used,
  which may be entirely incorrect for a specific problem.

  ARKStepSStolerances specifies scalar relative and absolute
    tolerances.
  ARKStepSVtolerances specifies scalar relative tolerance and a
    vector absolute tolerance (a potentially different absolute
    tolerance for each vector component).
  ARKStepWFtolerances specifies a user-provided function (of type
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
                     ARKStepCreate has not been called).
    ARK_NO_MALLOC    indicating that arkode_mem has not been
                     allocated.
    ARK_ILL_INPUT    indicating an input argument was illegal
                     (e.g. a negative tolerance)
  In case of an error return, an error message is also printed.
  --------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepSStolerances(void *arkode_mem,
                                        realtype reltol,
                                        realtype abstol);
SUNDIALS_EXPORT int ARKStepSVtolerances(void *arkode_mem,
                                        realtype reltol,
                                        N_Vector abstol);
SUNDIALS_EXPORT int ARKStepWFtolerances(void *arkode_mem,
                                        ARKEwtFn efun);


/*---------------------------------------------------------------
  ARKStepResStolerance, ARKStepResVtolerance, ARKStepResFtolerance

  These functions specify the absolute residual tolerance.
  Specification of the absolute residual tolerance is only
  necessary for problems with non-identity mass matrices in which
  the units of the solution vector y dramatically differ from the
  units of My, where M is the user-supplied mass matrix.  If this
  occurs, one of these routines SHOULD be called before the first
  call to ARKStepEvolve; otherwise the default value of
  rabstol=1e-9 will be used, which may be entirely incorrect for
  a specific problem.

  ARKStepResStolerances specifies a scalar residual tolerance.

  ARKStepResVtolerances specifies a vector residual tolerance
    (a potentially different absolute residual tolerance for
    each vector component).

  ARKStepResFtolerances specifies a user-provides function (of
    type ARKRwtFn) which will be called to set the residual
    weight vector.

  The tolerances reltol (defined for both the solution and
  residual) and rabstol define a vector of residual weights,
  rwt, with components
    rwt[i] = 1/(reltol*abs(My[i]) + abstol)     (in S case), or
    rwt[i] = 1/(reltol*abs(My[i]) + abstol[i])  (in V case).
  This vector is used in all solver convergence tests, which
  use a weighted RMS norm on all residual-like vectors v:
     WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v[i]*rwt[i])^2 ),
  where N is the problem dimension.

  The return value of these functions is equal to ARK_SUCCESS=0
  if there were no errors; otherwise it is a negative int equal
  to:
    ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
                     ARKStepCreate has not been called).
    ARK_NO_MALLOC    indicating that arkode_mem has not been
                     allocated.
    ARK_ILL_INPUT    indicating an input argument was illegal
                     (e.g. a negative tolerance)
  In case of an error return, an error message is also printed.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepResStolerance(void *arkode_mem,
                                         realtype rabstol);
SUNDIALS_EXPORT int ARKStepResVtolerance(void *arkode_mem,
                                         N_Vector rabstol);
SUNDIALS_EXPORT int ARKStepResFtolerance(void *arkode_mem,
                                         ARKRwtFn rfun);


/*---------------------------------------------------------------
  ARKStepSetLinearSolver and ARKStepSetMassLinearSolver

  These routines are required inputs for problems that require
  a linear solver for either the implicit stage solves or 
  non-identity mass matrix solves.

  ARKStepSetLinearSolver specifies the SUNLinearSolver object 
    that ARKode should use.  The 'LS' argument must be non-NULL, 
    but A can be NULL if the solver requires no SUNMatrix object.

  ARKStepSetMassLinearSolver specifies the SUNLinearSolver object 
    that ARKode should use for mass-matrix systems.  The 'LS' 
    argument must be non-NULL, but M can be NULL if the solver 
    requires no SUNMatrix object.

  NOTE: when solving an implicit or IMEX IVP with non-identity 
  mass matrix and iterative linear solver, both the system and
  mass solvers must be have the same type (i.e. you cannot 
  combine a direct system solver with an iterative mass matrix 
  solver, etc.).
  
  The return value is one of:
     ARKLS_SUCCESS   if successful
     ARKLS_MEM_NULL  if the ARKode memory was NULL
     ARKLS_ILL_INPUT if the linear solver memory was NULL
---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepSetLinearSolver(void *arkode_mem, 
                                           SUNLinearSolver LS,
                                           SUNMatrix A);

SUNDIALS_EXPORT int ARKStepSetMassLinearSolver(void *arkode_mem, 
                                               SUNLinearSolver LS,
                                               SUNMatrix M,
                                               booleantype time_dep);
  

/*---------------------------------------------------------------
  ARKStepRootInit

  This initializes a rootfinding problem to be solved during the
  integration of the ODE system.  It must be called after
  ARKStepCreate, and before ARKStepEvolve.  The arguments are:

  arkode_mem = pointer to ARKode memory returned by ARKStepCreate.

  nrtfn      = number of functions g_i, an integer >= 0.

  g          = name of user-supplied function, of type ARKRootFn,
               defining the functions g_i whose roots are sought.

  If a new problem is to be solved with a call to ARKStepReInit,
  where the new problem has no root functions but the prior one
  did, then call ARKStepRootInit again with nrtfn = 0.

  The return value is ARK_SUCCESS = 0 if there
  were no errors; otherwise it is a negative int equal to:
    ARK_MEM_NULL    indicating arkode_mem was NULL, or
    ARK_MEM_FAIL    indicating a memory allocation failed.
                    (including an attempt to increase maxord).
    ARK_ILL_INPUT   indicating nrtfn > 0 but g = NULL.
  In case of an error return, an error message is also printed.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepRootInit(void *arkode_mem, int nrtfn,
                                    ARKRootFn g);


/*---------------------------------------------------------------
  ARKStep optional input specification functions -- ALL of these
  must be called AFTER ARKStepCreate.  They can be called to set 
  optional inputs to non-default values.

  Function   [default value]
    Description
  -----------------------------------------------------------------
  ARKStepSetDefaults   [internal]
    resets all optional inputs to ARKStep default values.  Does not
    change problem-defining function pointers or user_data pointer.
    Also leaves alone any data structures/options related to
    root-finding (those can be reset using ARKodeRootInit).
  
  ARKStepSetOptimalParams   [internal]
    sets all adaptivity and solver parameters to our 'best guess'
    values, for a given integration method (ERK, DIRK, ARK) and a
    given method order.  Should only be called after the method
    order and integration method have been set.

  ARKStepSetOrder  [4]
    method order to be used by the solver.

  ARKStepSetDenseOrder  [3]
    polynomial order to be used for dense output.  Allowed values
    are between 0 and min(q,5) (where q is the order of the
    integrator)

  ARKStepSetNonlinearSolver  [SUNNonlinsol_Newton]
    attaches a non-default SUNNonlinearSolver object to be used in
    solving for implicit stage solutions.

  ARKStepSetErrHandlerFn  [internal]
    user-provided ErrHandler function.

  ARKStepSetErrFile   [stderr]
    the file pointer for an error file where all ARKStep warning
    and error messages will be written if the default internal
    error handling function is used. This parameter can be stdout
    (standard output), stderr (standard error), or a file pointer
    (corresponding to a user error file opened for writing)
    returned by fopen.  If not called, then all messages will
    be written to stderr.

  ARKStepSetUserData  [NULL]
    a pointer to user data that will be passed to the user's fi
    and fe functions every time they are called.

  ARKStepSetDiagnostics  [NULL]
    the file pointer for a diagnostics file where all ARKStep
    adaptivity and solver information is written.  This parameter
    can be stdout or stderr, though the preferred approach is to
    specify a file pointer corresponding to a user diagnostics file
    opened for writing) returned by fopen.  If not called, or if
    called with a NULL file pointer, all diagnostics output is
    disabled.
    NOTE: when run in parallel, only one process should set a
    non-NULL value for this pointer, since statistics from all
    processes would be identical.

  ARKStepSetLinear  [SUNFALSE]
    specifies that the implicit portion of the problem is linear,
    and to tighten the linear solver tolerances while taking only
    one Newton iteration.

  ARKStepSetNonlinear  [SUNTRUE]
    specifies that the implicit portion of the problem is nonlinear.
    Used to undo a previous call to ARKStepSetLinear

  ARKStepSetExplicit   [SUNFALSE]
    specifies that implicit portion of problem is disabled, and
    to use an explicit RK method.

  ARKStepSetImplicit   [SUNFALSE]
    specifies that explicit portion of problem is disabled, and to
    use an implicit RK method.

  ARKStepSetImEx   [SUNTRUE]
    specifies that problem has both implicit and explicit parts,
    and to use an ARK method.

  ARKStepSetTables  [determined by ARKode based on order]
    specifies to use customized Butcher tables for the IMEX system.

  ARKStepSetTableNum  [determined by ARKode based on order]
    specifies to use a built-in Butcher tables for the ImEx system.
    The integer arguments should match existing methods in
    ARKodeButcherTable_LoadERK() and ARKodeButcherTable_LoadDIRK().
    Error-checking is performed to ensure that the tables exist.

  ARKStepSetMaxNumSteps  [500]
    maximum number of internal steps to be taken by the solver in
    its attempt to reach tout.

  ARKStepSetMaxHnilWarns  [10]
    maximum number of warning messages issued by the solver that
    t+h==t on the next internal step. A value of -1 means no such
    messages are issued.

  ARKStepSetInitStep  [estimated internally]
    initial step size.

  ARKStepSetMinStep  [0.0]
    minimum absolute value of step size allowed.

  ARKStepSetMaxStep  [infinity]
    maximum absolute value of step size allowed.

  ARKStepSetStopTime  [infinity]
    the independent variable value past which the solution is
    not to proceed.

  ARKStepSetFixedStep  [off]
    specifies to use a fixed step size throughout integration

  ARKStepSetCFLFraction  [0.5]
    safety factor to use for explicitly stable steps

  ARKStepSetSafetyFactor  [0.96]
    safety factor to use for error-based step adaptivity

  ARKStepSetErrorBias  [1.5]
    error bias factor to use in error-based step adaptivity

  ARKStepSetMaxGrowth  [20.0]
    maximum growth factor for successive time steps (not
    including the first step).

  ARKStepSetMaxFirstGrowth  [10000.0]
    maximum growth factor for first step.

  ARKStepSetMaxEFailGrowth  [0.3]
    maximum growth factor after an error failure.

  ARKStepSetSmallNumEFails  [2]
    maximum number of error failures before MaxFailGrowth factor
    is used.

  ARKStepSetMaxCFailGrowth  [0.25]
    maximum growth factor after a convergence failure.

  ARKStepSetFixedStepBounds  [1.0 1.5]
    step growth interval to force retention of the same step size

  ARKStepSetAdaptivityMethod  [0]
    Method to use for time step adaptivity

  ARKStepSetAdaptivityFn  [internal]
    user-provided time step adaptivity function.

  ARKStepSetNonlinCRDown  [0.3]
    user-provided nonlinear convergence rate constant.

  ARKStepSetNonlinRDiv  [2.3]
    user-provided nonlinear divergence ratio.

  ARKStepSetDeltaGammaMax  [0.2]
    user-provided linear setup decision constant.

  ARKStepSetMaxStepsBetweenLSet  [20]
    user-provided linear setup decision constant.

  ARKStepSetPredictorMethod  [0]
    Method to use for predicting implicit solutions.

  ARKStepSetStabilityFn  [internal]
    user-provided explicit time step stability function.

  ARKStepSetMaxErrTestFails  [7]
    Maximum number of error test failures in attempting one step.

  ARKStepSetMaxNonlinIters  [3]
    Maximum number of nonlinear solver iterations at one stage solution.

  ARKStepSetMaxConvFails  [10]
    Maximum number of convergence failures allowed in attempting one step.
 
  ARKStepSetNonlinConvCoef  [0.1]
    Coefficient in the nonlinear convergence test.

  ARKStepSetRootDirection  [both directions]
    Specifies the direction of zero crossings to be monitored

  ARKStepSetNoInactiveRootWarn  [enabled]
    disable warning about possible g==0 at beginning of integration
  -----------------------------------------------------------------
  Return flag:
  ARK_SUCCESS   if successful
  ARK_MEM_NULL  if the arkode memory is NULL
  ARK_ILL_INPUT if an argument has an illegal value
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int ARKStepSetOptimalParams(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetOrder(void *arkode_mem, int maxord);
SUNDIALS_EXPORT int ARKStepSetDenseOrder(void *arkode_mem, int dord);
SUNDIALS_EXPORT int ARKStepSetNonlinearSolver(void *arkode_mem,
                                              SUNNonlinearSolver NLS);
SUNDIALS_EXPORT int ARKStepSetLinear(void *arkode_mem, int timedepend);
SUNDIALS_EXPORT int ARKStepSetNonlinear(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetExplicit(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetImplicit(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetImEx(void *arkode_mem);
SUNDIALS_EXPORT int ARKStepSetTables(void *arkode_mem, int q, int p,
                                     ARKodeButcherTable Bi,
                                     ARKodeButcherTable Be);
SUNDIALS_EXPORT int ARKStepSetTableNum(void *arkode_mem,
                                       int itable, int etable);
SUNDIALS_EXPORT int ARKStepSetCFLFraction(void *arkode_mem,
                                          realtype cfl_frac);
SUNDIALS_EXPORT int ARKStepSetSafetyFactor(void *arkode_mem,
                                           realtype safety);
SUNDIALS_EXPORT int ARKStepSetErrorBias(void *arkode_mem,
                                        realtype bias);
SUNDIALS_EXPORT int ARKStepSetMaxGrowth(void *arkode_mem,
                                        realtype mx_growth);
SUNDIALS_EXPORT int ARKStepSetFixedStepBounds(void *arkode_mem,
                                              realtype lb, realtype ub);
SUNDIALS_EXPORT int ARKStepSetAdaptivityMethod(void *arkode_mem,
                                               int imethod,
                                               int idefault, int pq,
                                               realtype *adapt_params);
SUNDIALS_EXPORT int ARKStepSetAdaptivityFn(void *arkode_mem,
                                           ARKAdaptFn hfun,
                                           void *h_data);
SUNDIALS_EXPORT int ARKStepSetMaxFirstGrowth(void *arkode_mem,
                                             realtype etamx1);
SUNDIALS_EXPORT int ARKStepSetMaxEFailGrowth(void *arkode_mem,
                                             realtype etamxf);
SUNDIALS_EXPORT int ARKStepSetSmallNumEFails(void *arkode_mem,
                                             int small_nef);
SUNDIALS_EXPORT int ARKStepSetMaxCFailGrowth(void *arkode_mem,
                                             realtype etacf);
SUNDIALS_EXPORT int ARKStepSetNonlinCRDown(void *arkode_mem,
                                           realtype crdown);
SUNDIALS_EXPORT int ARKStepSetNonlinRDiv(void *arkode_mem,
                                         realtype rdiv);
SUNDIALS_EXPORT int ARKStepSetDeltaGammaMax(void *arkode_mem,
                                            realtype dgmax);
SUNDIALS_EXPORT int ARKStepSetMaxStepsBetweenLSet(void *arkode_mem,
                                                  int msbp);
SUNDIALS_EXPORT int ARKStepSetPredictorMethod(void *arkode_mem,
                                              int method);
SUNDIALS_EXPORT int ARKStepSetStabilityFn(void *arkode_mem,
                                          ARKExpStabFn EStab,
                                          void *estab_data);
SUNDIALS_EXPORT int ARKStepSetMaxErrTestFails(void *arkode_mem,
                                              int maxnef);
SUNDIALS_EXPORT int ARKStepSetMaxNonlinIters(void *arkode_mem,
                                             int maxcor);
SUNDIALS_EXPORT int ARKStepSetMaxConvFails(void *arkode_mem,
                                           int maxncf);
SUNDIALS_EXPORT int ARKStepSetNonlinConvCoef(void *arkode_mem,
                                             realtype nlscoef);
SUNDIALS_EXPORT int ARKStepSetMaxNumSteps(void *arkode_mem,
                                          long int mxsteps);
SUNDIALS_EXPORT int ARKStepSetMaxHnilWarns(void *arkode_mem,
                                           int mxhnil);
SUNDIALS_EXPORT int ARKStepSetInitStep(void *arkode_mem,
                                       realtype hin);
SUNDIALS_EXPORT int ARKStepSetMinStep(void *arkode_mem,
                                      realtype hmin);
SUNDIALS_EXPORT int ARKStepSetMaxStep(void *arkode_mem,
                                      realtype hmax);
SUNDIALS_EXPORT int ARKStepSetStopTime(void *arkode_mem,
                                       realtype tstop);
SUNDIALS_EXPORT int ARKStepSetFixedStep(void *arkode_mem,
                                        realtype hfixed);

SUNDIALS_EXPORT int ARKStepSetRootDirection(void *arkode_mem,
                                            int *rootdir);
SUNDIALS_EXPORT int ARKStepSetNoInactiveRootWarn(void *arkode_mem);

SUNDIALS_EXPORT int ARKStepSetErrHandlerFn(void *arkode_mem,
                                           ARKErrHandlerFn ehfun,
                                           void *eh_data);
SUNDIALS_EXPORT int ARKStepSetErrFile(void *arkode_mem,
                                      FILE *errfp);
SUNDIALS_EXPORT int ARKStepSetUserData(void *arkode_mem,
                                       void *user_data);
SUNDIALS_EXPORT int ARKStepSetDiagnostics(void *arkode_mem,
                                          FILE *diagfp);

SUNDIALS_EXPORT int ARKStepSetPostprocessStepFn(void *arkode_mem,
                                                ARKPostProcessStepFn ProcessStep);

/*---------------------------------------------------------------
  ARKStep linear solver interface optional input specification 
  functions -- ALL of these must be called AFTER 
  ARKStepSetLinearSolver and/or ARKStepSetMassLinearSolver.  

  ARKStepSetJacFn specifies the Jacobian approximation routine to
    be used when constructing J.  By default, a difference 
    quotient approximation is used for dense/band SUNMatrix 
    objects; for all other matrix types passed to 
    ARKStepSetLinearSolver, this routine must be user-supplied).  
    See "arkode/arkode_ls.h" for a complete description of the 
    function prototype.

  ARKStepSetMassFn specifies the mass matrix approximation routine 
    to be used for constructing M.  This MUST be supplied if the 
    SUNMatrix M supplied to ARKStepSetMassLinearSolver was 
    non-NULL; otherwise it defaults to NULL.  See 
    "arkode/arkode_ls.h" for a complete description of the 
    function prototype.

  ARKStepSetMaxStepsBetweenJac specifies the maximum number of time
    steps to wait before recomputation of the Jacobian or 
    recommendation to update the preconditioner.  This differs from 
    the ARKStepSetMaxStepsBetweenLSet, which merely indicates the 
    frequency with which the linear solver setup routine is called.  
    Default value is 50.

  ARKStepSetEpsLin specifies the factor by which the tolerance on
    the nonlinear iteration is multiplied to get a tolerance on 
    the linear iteration.  Default value is 0.05.

  ARKStepSetMassEpsLin specifies the factor by which the tolerance
    on the nonlinear iteration is multiplied to get a tolerance 
    on the mass matrix linear iteration.  Default value is 0.05.

  ARKStepSetPreconditioner specifies the PrecSetup and PrecSolve 
    functions.  Default is NULL for both arguments (no 
    preconditioning).  See "arkode/arkode_ls.h" for a complete 
    description of the function prototype.

  ARKStepSetMassPreconditioner specifies the mass matrix MPrecSetup
    and MPrecSolve functions.  Default is NULL for both arguments
    (no preconditioning).  See "arkode/arkode_ls.h" for a complete 
    description of the function prototype.

  ARKStepSetJacTimes specifies the jtsetup and jtimes functions. 
    Default is to use an internal finite difference approximation
    routine (no setup).  See "arkode/arkode_ls.h" for a complete 
    description of the function prototype.

  ARKStepSetMassTimes specifies the mtsetup and mtimes functions. 
    Note that there do not exist built-in finite-difference 
    approximation routines for this.  Hence, either 
    ARKStepSetMassLinearSolver must have been called with non-NULL 
    SUNMatrix M, or this function MUST be called with non-NULL 
    'mtimes'.  See "arkode/arkode_ls.h" for a complete 
    description of the function prototype.

  The return value of ARKStepSet* is one of:
    ARKLS_SUCCESS      if successful
    ARKLS_MEM_NULL     if the arkode memory was NULL
    ARKLS_LMEM_NULL    if the linear solver memory was NULL
    ARKLS_MASSMEM_NULL if the mass matrix solver memory was NULL
    ARKLS_ILL_INPUT    if an input has an illegal value
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepSetJacFn(void *arkode_mem, ARKLsJacFn jac);
SUNDIALS_EXPORT int ARKStepSetMassFn(void *arkode_mem, ARKLsMassFn mass);
SUNDIALS_EXPORT int ARKStepSetMaxStepsBetweenJac(void *arkode_mem,
                                                 long int msbj);
SUNDIALS_EXPORT int ARKStepSetEpsLin(void *arkode_mem, realtype eplifac);
SUNDIALS_EXPORT int ARKStepSetMassEpsLin(void *arkode_mem, realtype eplifac);
SUNDIALS_EXPORT int ARKStepSetPreconditioner(void *arkode_mem, 
                                             ARKLsPrecSetupFn psetup,
                                             ARKLsPrecSolveFn psolve);
SUNDIALS_EXPORT int ARKStepSetMassPreconditioner(void *arkode_mem, 
                                                 ARKLsMassPrecSetupFn psetup,
                                                 ARKLsMassPrecSolveFn psolve);
SUNDIALS_EXPORT int ARKStepSetJacTimes(void *arkode_mem, 
                                       ARKLsJacTimesSetupFn jtsetup,
                                       ARKLsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int ARKStepSetMassTimes(void *arkode_mem, 
                                        ARKLsMassTimesSetupFn msetup,
                                        ARKLsMassTimesVecFn mtimes,
                                        void *mtimes_data);
  

/*---------------------------------------------------------------
  ARKStepEvolve

  This integrates the ODE over an interval in t.

  ARKStepEvolve may be run in one of two modes (ARK_NORMAL or
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
  may issue a call to ARKStepSetStopTime before the call to
  ARKStepEvolve to specify a fixed stop time to end the time step
  and return to the user.  Once the integrator returns at a tstop
  time, any future testing for tstop is disabled (and can be
  reenabled only though a new call to ARKStepSetStopTime).

  If itask is ARK_ONE_STEP, then the solver takes one internal
  time step and returns in yout the value of y at the new internal
  time. In this case, tout is used only during the first call to
  ARKStepEvolve to determine the direction of integration and the
  rough scale of the t variable.  As with the ARK_NORMAL mode, a
  user may specify a specific stop time for output of this step,
  assuming that the requested step is smaller than the step taken
  by the method.

  The time reached by the solver is placed in (*tret). The
  user is responsible for allocating the memory for this value.

  arkode_mem is the pointer to ARKStep memory returned by
             ARKStepCreate.

  tout  is the next time at which a computed solution is desired.

  yout  is the computed solution vector. In ARK_NORMAL mode with no
        errors and no roots found, yout=y(tout).

  tret  is a pointer to a real location. ARKStepEvolve sets (*tret)
        to the time reached by the solver and returns yout=y(*tret).

  itask is ARK_NORMAL or ARK_ONE_STEP, as described above.

  Here is a brief description of each return value:

  ARK_SUCCESS:      ARKStepEvolve succeeded and no roots were found.

  ARK_ROOT_RETURN:  ARKStepEvolve succeeded, and found one or more roots.
                    If nrtfn > 1, call ARKStepGetRootInfo to see
                    which g_i were found to have a root at (*tret).

  ARK_TSTOP_RETURN: ARKStepEvolve succeeded and returned at tstop.

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
                    after calling ARKStepCreate) failed to set one of
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
SUNDIALS_EXPORT int ARKStepEvolve(void *arkode_mem, realtype tout,
                                  N_Vector yout, realtype *tret,
                                  int itask);


/*---------------------------------------------------------------
  ARKStepGetDky

  This computes the kth derivative of the y function at time t,
  where tn-hu <= t <= tn, tn denotes the current internal time
  reached, and hu is the last internal step size successfully
  used by the solver. The user may request k=0, 1, ..., d, where
  d = min(5,q), with q the order of accuracy for the time
  integration method. The derivative vector is returned in dky.
  This vector must be allocated by the caller. It is only legal
  to call this function after a successful return from
  ARKStepEvolve.

  arkode_mem is the pointer to ARKStep memory returned by
             ARKStepCreate.

  t   is the time at which the kth derivative of y is evaluated.
      The legal range for t is [tn-hu,tn] as described above.

  k   is the order of the derivative of y to be computed. The
      legal range for k is [0,min(q,3)] as described above.

  dky is the output derivative vector [((d/dy)^k)y](t).

  The return value for ARKStepGetDky is one of:

    ARK_SUCCESS:  ARKStepGetDky succeeded.

    ARK_BAD_K:    k is not in the range 0, 1, ..., s-1.

    ARK_BAD_T:    t is not in the interval [tn-hu,tn].

    ARK_BAD_DKY:  The dky argument was NULL.

    ARK_MEM_NULL: The arkode_mem argument was NULL.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepGetDky(void *arkode_mem, realtype t,
                                  int k, N_Vector dky);


/*---------------------------------------------------------------
  Optional outputs from the ARKStep module:  the following
  functions can be called to get optional outputs and statistics
  related to the main integrator.

  ARKStepGetNumExpSteps returns the cumulative number of stability
    limited steps taken by the solver

  ARKStepGetNumAccSteps returns the cumulative number of accuracy
    limited steps taken by the solver

  ARKStepGetNumStepAttempts returns the total number of steps
    attempted by the solver

  ARKStepGetNumRhsEvals returns the number of calls to the user's
    fe and fi functions

  ARKStepGetNumLinSolvSetups returns the number of calls made to
    the linear solver's setup routine

  ARKStepGetNumErrTestFails returns the number of local error test
    failures that have occured

  ARKStepGetCurrentButcherTables returns the explicit and implicit
    Butcher tables currently in use

  ARKStepGetEstLocalErrors returns the vector of estimated local
    errors. The user must allocate space for ele.

  ARKStepGetNumSteps returns the cumulative number of internal
    steps taken by the solver

  ARKStepGetActualInitStep returns the actual initial step size
    used by ARKStep

  ARKStepGetLastStep returns the step size for the last internal
    step

  ARKStepGetCurrentStep returns the step size to be attempted on
    the next internal step

  ARKStepGetCurrentTime returns the current internal time reached
    by the solver

  ARKStepGetTolScaleFactor returns a suggested factor by which the
    user's tolerances should be scaled when too much accuracy has 
    been requested for some internal step

  ARKStepGetErrWeights returns the current error weight vector.
    The user must allocate space for eweight.

  ARKStepGetResWeights returns the current residual weight vector.
    The user must allocate space for rweight.

  ARKStepGetWorkSpace returns the ARKStep real and integer workspaces

  ARKStepGetNumGEvals returns the number of calls to the user's
    g function (for rootfinding)

  ARKStepGetRootInfo returns the indices for which g_i was found to
    have a root. The user must allocate space for rootsfound. For 
    i = 0 ... nrtfn-1, rootsfound[i] = 1 if g_i has a root, 
    and = 0  if not.

  The return value of ARKStepGet* is one of:
     ARK_SUCCESS   if successful
     ARK_MEM_NULL  if the ARKode or ARKStep memory structures
                   were NULL
     ARK_LMEM_NULL if a linear solver memory structure was NULL
  ---------------------------------------------------------------*/

SUNDIALS_EXPORT int ARKStepGetNumExpSteps(void *arkode_mem,
                                          long int *expsteps);
SUNDIALS_EXPORT int ARKStepGetNumAccSteps(void *arkode_mem,
                                          long int *accsteps);
SUNDIALS_EXPORT int ARKStepGetNumStepAttempts(void *arkode_mem,
                                              long int *step_attempts);
SUNDIALS_EXPORT int ARKStepGetNumRhsEvals(void *arkode_mem,
                                          long int *nfe_evals,
                                          long int *nfi_evals);
SUNDIALS_EXPORT int ARKStepGetNumLinSolvSetups(void *arkode_mem,
                                               long int *nlinsetups);
SUNDIALS_EXPORT int ARKStepGetNumErrTestFails(void *arkode_mem,
                                              long int *netfails);
SUNDIALS_EXPORT int ARKStepGetCurrentButcherTables(void *arkode_mem,
                                                   ARKodeButcherTable *Bi,
                                                   ARKodeButcherTable *Be);
SUNDIALS_EXPORT int ARKStepGetEstLocalErrors(void *arkode_mem,
                                             N_Vector ele);
SUNDIALS_EXPORT int ARKStepGetWorkSpace(void *arkode_mem,
                                        long int *lenrw,
                                        long int *leniw);
SUNDIALS_EXPORT int ARKStepGetNumSteps(void *arkode_mem,
                                       long int *nsteps);
SUNDIALS_EXPORT int ARKStepGetActualInitStep(void *arkode_mem,
                                             realtype *hinused);
SUNDIALS_EXPORT int ARKStepGetLastStep(void *arkode_mem,
                                       realtype *hlast);
SUNDIALS_EXPORT int ARKStepGetCurrentStep(void *arkode_mem,
                                          realtype *hcur);
SUNDIALS_EXPORT int ARKStepGetCurrentTime(void *arkode_mem,
                                          realtype *tcur);
SUNDIALS_EXPORT int ARKStepGetTolScaleFactor(void *arkode_mem,
                                             realtype *tolsfac);
SUNDIALS_EXPORT int ARKStepGetErrWeights(void *arkode_mem,
                                         N_Vector eweight);
SUNDIALS_EXPORT int ARKStepGetResWeights(void *arkode_mem,
                                         N_Vector rweight);
SUNDIALS_EXPORT int ARKStepGetNumGEvals(void *arkode_mem,
                                        long int *ngevals);
SUNDIALS_EXPORT int ARKStepGetRootInfo(void *arkode_mem,
                                       int *rootsfound);

/*---------------------------------------------------------------
  As a convenience, the following functions provides many of the
  above optional outputs in grouped form.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepGetTimestepperStats(void *arkode_mem,
                                               long int *expsteps,
                                               long int *accsteps,
                                               long int *step_attempts,
                                               long int *nfe_evals,
                                               long int *nfi_evals,
                                               long int *nlinsetups,
                                               long int *netfails);
SUNDIALS_EXPORT int ARKStepGetStepStats(void *arkode_mem,
                                        long int *nsteps,
                                        realtype *hinused,
                                        realtype *hlast,
                                        realtype *hcur,
                                        realtype *tcur);

/*---------------------------------------------------------------
  Nonlinear solver optional output extraction functions
  -----------------------------------------------------------------
  The following functions can be called to get optional outputs
  and statistics related to the nonlinear solver.
  -----------------------------------------------------------------
  ARKStepGetNumNonlinSolvIters returns the number of nonlinear
    solver iterations performed.

  ARKStepGetNumNonlinSolvConvFails returns the number of nonlinear
    convergence failures.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepGetNumNonlinSolvIters(void *arkode_mem,
                                                long int *nniters);
SUNDIALS_EXPORT int ARKStepGetNumNonlinSolvConvFails(void *arkode_mem,
                                                    long int *nncfails);

/*---------------------------------------------------------------
  As a convenience, the following function provides the
  nonlinear solver optional outputs in a group.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepGetNonlinSolvStats(void *arkode_mem,
                                             long int *nniters,
                                             long int *nncfails);

/*---------------------------------------------------------------
  Linear solver optional output extraction functions
  -----------------------------------------------------------------
  The following functions can be called to get optional outputs
  and statistics related to the linear solvers.
  -----------------------------------------------------------------
  ARKStepGetLinWorkSpace returns the real and integer workspace 
    used by the ARKLs interface.

  ARKStepGetNumJacEvals returns the number of calls made to the
    Jacobian evaluation routine jac.
  
  ARKStepGetNumPrecEvals returns the number of preconditioner
    evaluations, i.e. the number of calls made
    to PrecSetup with jok==SUNFALSE.
  
  ARKStepGetNumPrecSolves returns the number of calls made to
    PrecSolve.
  
  ARKStepGetNumLinIters returns the number of linear iterations.
  
  ARKStepGetNumLinConvFails returns the number of linear
    convergence failures.
  
  ARKStepGetNumJTSetupEvals returns the number of calls to jtsetup.
  
  ARKStepGetNumJtimesEvals returns the number of calls to jtimes.
  
  ARKStepGetNumLinRhsEvals returns the number of calls to the user
    f routine due to perform finite difference Jacobian times 
    vector evaluations or Jacobian matrix approximations.
  
  ARKStepGetLastLinFlag returns the last error flag set by any of
    the ARKLS interface functions.
  
  ARKStepGetMassWorkSpace returns the real and integer workspace used
    by the mass matrix ARKLs interface.
  
  ARKStepGetNumMassSetups returns the number of calls made to the
    mass matrix solver setup routine
  
  ARKStepGetNumMassMult returns the number of calls to either the 
    internal or user-supplied mtimes routine.
  
  ARKStepGetNumMassSolves returns the number of calls made to the
    mass matrix solver 'solve' routine
  
  ARKStepGetNumMassPrecEvals returns the number of mass matrix 
    preconditioner evaluations, i.e. the number of calls made to 
    MPrecSetup.
  
  ARKStepGetNumMassPrecSolves returns the number of calls made to
    MPrecSolve.
  
  ARKStepGetNumMassIters returns the number of mass matrix solver
    iterations.
  
  ARKStepGetNumMassConvFails returns the number of mass matrix solver
    convergence failures.
  
  ARKStepGetNumMTSetups returns the number of calls to mtsetup.
  
  ARKStepGetLastMassFlag returns the last error flag set by any of
    the ARKLs interface functions on the mass matrix solve.
  
  ARKStepGetLinReturnFlagName returns the name of the constant 
    associated with a ARKLS return flag.
  
  The return value of the above routines is one of:
     ARKLS_SUCCESS   if successful
     ARKLS_MEM_NULL  if the arkode memory was NULL
     ARKLS_LMEM_NULL if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepGetLinWorkSpace(void *arkode_mem, 
                                           long int *lenrwLS, 
                                           long int *leniwLS);
SUNDIALS_EXPORT int ARKStepGetNumJacEvals(void *arkode_mem, 
                                          long int *njevals);
SUNDIALS_EXPORT int ARKStepGetNumPrecEvals(void *arkode_mem, 
                                           long int *npevals);
SUNDIALS_EXPORT int ARKStepGetNumPrecSolves(void *arkode_mem, 
                                            long int *npsolves);
SUNDIALS_EXPORT int ARKStepGetNumLinIters(void *arkode_mem, 
                                          long int *nliters);
SUNDIALS_EXPORT int ARKStepGetNumLinConvFails(void *arkode_mem, 
                                              long int *nlcfails);
SUNDIALS_EXPORT int ARKStepGetNumJTSetupEvals(void *arkode_mem,
                                              long int *njtsetups);
SUNDIALS_EXPORT int ARKStepGetNumJtimesEvals(void *arkode_mem, 
                                             long int *njvevals);
SUNDIALS_EXPORT int ARKStepGetNumLinRhsEvals(void *arkode_mem, 
                                             long int *nfevalsLS); 
SUNDIALS_EXPORT int ARKStepGetLastLinFlag(void *arkode_mem, 
                                          long int *flag);

SUNDIALS_EXPORT int ARKStepGetMassWorkSpace(void *arkode_mem, 
                                            long int *lenrwMLS, 
                                            long int *leniwMLS);
SUNDIALS_EXPORT int ARKStepGetNumMassSetups(void *arkode_mem, 
                                            long int *nmsetups);
SUNDIALS_EXPORT int ARKStepGetNumMassMult(void *arkode_mem, 
                                          long int *nmvevals);
SUNDIALS_EXPORT int ARKStepGetNumMassSolves(void *arkode_mem, 
                                            long int *nmsolves);
SUNDIALS_EXPORT int ARKStepGetNumMassPrecEvals(void *arkode_mem, 
                                               long int *nmpevals);
SUNDIALS_EXPORT int ARKStepGetNumMassPrecSolves(void *arkode_mem, 
                                                long int *nmpsolves);
SUNDIALS_EXPORT int ARKStepGetNumMassIters(void *arkode_mem, 
                                           long int *nmiters);
SUNDIALS_EXPORT int ARKStepGetNumMassConvFails(void *arkode_mem, 
                                               long int *nmcfails);
SUNDIALS_EXPORT int ARKStepGetNumMTSetups(void *arkode_mem,
                                          long int *nmtsetups);
SUNDIALS_EXPORT int ARKStepGetLastMassFlag(void *arkode_mem, 
                                           long int *flag);

SUNDIALS_EXPORT char *ARKStepGetLinReturnFlagName(long int flag);
  
  
/*---------------------------------------------------------------
  The following function returns the name of the constant
  associated with a ARKStep return flag
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT char *ARKStepGetReturnFlagName(long int flag);

/*---------------------------------------------------------------
  Function : ARKStepWriteParameters
  -----------------------------------------------------------------
  ARKStepWriteParameters outputs all timestepper module parameters
  to the provided file pointer.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepWriteParameters(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
  Function : ARKStepWriteButcher
  -----------------------------------------------------------------
  ARKStepWriteButcher outputs the Butcher tables to the
  provided file pointer.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int ARKStepWriteButcher(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
  ARKStepFree

  This frees the problem memory arkode_mem allocated by
  ARKStepCreate. Its only argument is the pointer arkode_mem
  returned by ARKStepCreate.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void ARKStepFree(void **arkode_mem);

/*---------------------------------------------------------------
  ARKStepPrintMem

  This routine outputs the memory from the ARKStep structure and
  the main ARKode infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void ARKStepPrintMem(void* arkode_mem, FILE* outfile);

#ifdef __cplusplus
}
#endif

#endif
