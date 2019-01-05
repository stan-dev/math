/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------------------
 * Header file for the MRI time step module for ARKode.
 * ---------------------------------------------------------------------------*/

#ifndef _MRISTEP_H
#define _MRISTEP_H

#include <sundials/sundials_nvector.h>
#include <arkode/arkode.h>
#include <arkode/arkode_butcher_erk.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  MRISTEP Constants
  ===============================================================*/

/* Default Butcher tables for each order */
#define DEFAULT_MRI_STABLE_3  KNOTH_WOLKE_3_3
#define DEFAULT_MRI_FTABLE_3  KNOTH_WOLKE_3_3

/*===============================================================
  MRISTEP Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  MRIStepCreate

  This creates an internal memory block for a problem to be
  solved by ARKode, using the MRI time step module.  If successful,
  it returns a pointer to initialized problem memory. This
  pointer should be passed to all MRIStep-related routines.
  If an initialization error occurs, this routine will print an
  error message to standard err and return NULL.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void* MRIStepCreate(ARKRhsFn fs, ARKRhsFn ff,
                                    realtype t0, N_Vector y0);


/*---------------------------------------------------------------
  MRIStepResize

  This re-initializes the MRIStep memory for a problem with a
  changing vector size.  It is assumed that the problem dynamics
  before and after the vector resize will be comparable, so that
  the current step sizes remain valid after the call.  If instead
  the step sizes need to be updated a call to MRIStepSetFixedStep
  should be made after resizing.

  To aid in the vector-resize operation, the user can supply a
  vector resize function, that will take as input an N_Vector with
  the previous size, and return as output a corresponding vector
  of the new size.  If this function (of type ARKVecResizeFn) is
  not supplied (i.e. is set to NULL), then all existing N_Vectors
  will be destroyed and re-cloned from the input vector.

  Other arguments:
    arkode_mem       Existing MRIStep memory data structure.
    ynew             The newly-sized solution vector, holding
                     the current dependent variable values.
    t0               The current value of the independent
                     variable.
    resize_data      User-supplied data structure that will be
                     passed to the supplied resize function.

  The return value of MRIStepResize is equal to ARK_SUCCESS = 0 if
  there were no errors; otherwise it is a negative int equal to:
    ARK_MEM_NULL     indicating arkode_mem was NULL (i.e.,
                     MRIStepCreate has not been called).
    ARK_NO_MALLOC    indicating that arkode_mem has not been
                     allocated.
    ARK_ILL_INPUT    indicating an input argument was illegal
                     (including an error from the supplied
                     resize function).
  In case of an error return, an error message is also printed.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int MRIStepResize(void *arkode_mem, N_Vector ynew,
                                  realtype t0, ARKVecResizeFn resize,
                                  void *resize_data);


/*---------------------------------------------------------------
  MRIStepReInit

  This re-initializes the MRI time step module.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int MRIStepReInit(void* arkode_mem, ARKRhsFn fs, ARKRhsFn ff,
                                  realtype t0, N_Vector y0);


/*---------------------------------------------------------------
  MRIStepRootInit

  This initializes a rootfinding problem to be solved during the
  integration of the ODE system.  It must be called after
  MRIStepCreate, and before MRIStepEvolve.  The arguments are:

  arkode_mem = pointer to ARKode memory returned by MRIStepCreate.

  nrtfn      = number of functions g_i, an integer >= 0.

  g          = name of user-supplied function, of type ARKRootFn,
               defining the functions g_i whose roots are sought.

  If a new problem is to be solved with a call to MRIStepReInit,
  where the new problem has no root functions but the prior one
  did, then call MRIStepRootInit again with nrtfn = 0.

  The return value is ARK_SUCCESS = 0 if there
  were no errors; otherwise it is a negative int equal to:
    ARK_MEM_NULL    indicating arkode_mem was NULL, or
    ARK_MEM_FAIL    indicating a memory allocation failed.
                    (including an attempt to increase maxord).
    ARK_ILL_INPUT   indicating nrtfn > 0 but g = NULL.
  In case of an error return, an error message is also printed.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int MRIStepRootInit(void *arkode_mem, int nrtfn,
                                    ARKRootFn g);


/*---------------------------------------------------------------
  MRIStep optional input specification functions -- ALL of these
  must be called AFTER MRIStepCreate.
  -----------------------------------------------------------------
  The following functions can be called to set optional inputs
  to values other than the defaults given below.

  Function                   |  Optional input / [ default value ]
  -----------------------------------------------------------------
  MRIStepSetDefaults         | resets all optional inputs to MRIStep
                             | default values.  Does not change
                             | problem-defining function pointers or
                             | user_data pointer.  Also leaves alone
                             | any data structures/options related
                             | to root-finding (those can be reset
                             | using ARKodeRootInit).
                             | [internal]
                             |
  MRIStepSetDenseOrder       | polynomial order to be used for dense
                             | output.  Allowed values are between 0
                             | and min(q,5) (where q is the order of
                             | the integrator)
                             | [3]
                             |
  MRIStepSetErrHandlerFn     | user-provided ErrHandler function.
                             | [internal]
                             |
  MRIStepSetErrFile          | the file pointer for an error file
                             | where all MRIStep warning and error
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
  MRIStepSetUserData         | a pointer to user data that will be
                             | passed to the user's f function every
                             | time f is called.
                             | [NULL]
                             |
  MRIStepSetDiagnostics      | the file pointer for a diagnostics file
                             | where all MRIStep adaptivity and solver
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
  MRIStepSetTables           | specifies to use customized Butcher
                             | tables for the slow and fast portions of the
                             | system.
                             | [KNOTH_WOLKE_3_3 / KNOTH_WOLKE_3_3]
                             |
  MRIStepSetTableNum         | specifies to use built-in Butcher
                             | tables for the slow and fast portions of
                             | the system.  The integer argument should
                             | match an existing method in
                             | ARKodeButcherTable_LoadERK() within the
                             | file arkode_butcher_erk.c.  Error-checking
                             | is performed to ensure that these tables
                             | exist, and are not implicit.
                             | [KNOTH_WOLKE_3_3 / KNOTH_WOLKE_3_3]
                             |
  MRIStepSetMaxNumSteps      | maximum number of internal steps to be
                             | taken by the solver in its attempt to
                             | reach tout.
                             | [500]
                             |
  MRIStepSetMaxHnilWarns     | maximum number of warning messages
                             | issued by the solver that t+h==t on
                             | the next internal step. A value of -1
                             | means no such messages are issued.
                             | [10]
                             |
  MRIStepSetStopTime         | the independent variable value past
                             | which the solution is not to proceed.
                             | [infinity]
                             |
  MRIStepSetFixedStep        | specifies to use a fixed step size
                             | throughout integration
                             | [off]
                             |
  -----------------------------------------------------------------
  MRIStepSetRootDirection      | Specifies the direction of zero
                               | crossings to be monitored
                               | [both directions]
                               |
  MRIStepSetNoInactiveRootWarn | disable warning about possible
                               | g==0 at beginning of integration
  -----------------------------------------------------------------
  Return flag:
    ARK_SUCCESS   if successful
    ARK_MEM_NULL  if the arkode memory is NULL
    ARK_ILL_INPUT if an argument has an illegal value
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int MRIStepSetDefaults(void* arkode_mem);
SUNDIALS_EXPORT int MRIStepSetDenseOrder(void *arkode_mem, int dord);
SUNDIALS_EXPORT int MRIStepSetTables(void *arkode_mem,
                                     int q,
                                     ARKodeButcherTable Bs,
                                     ARKodeButcherTable Bf);
SUNDIALS_EXPORT int MRIStepSetTableNum(void *arkode_mem, int istable,
                                       int iftable);
SUNDIALS_EXPORT int MRIStepSetMaxNumSteps(void *arkode_mem,
                                          long int mxsteps);
SUNDIALS_EXPORT int MRIStepSetMaxHnilWarns(void *arkode_mem,
                                           int mxhnil);
SUNDIALS_EXPORT int MRIStepSetStopTime(void *arkode_mem,
                                       realtype tstop);
SUNDIALS_EXPORT int MRIStepSetFixedStep(void *arkode_mem,
                                        realtype hsfixed,
                                        realtype hffixed);

SUNDIALS_EXPORT int MRIStepSetRootDirection(void *arkode_mem,
                                            int *rootdir);
SUNDIALS_EXPORT int MRIStepSetNoInactiveRootWarn(void *arkode_mem);

SUNDIALS_EXPORT int MRIStepSetErrHandlerFn(void *arkode_mem,
                                           ARKErrHandlerFn ehfun,
                                           void *eh_data);
SUNDIALS_EXPORT int MRIStepSetErrFile(void *arkode_mem,
                                      FILE *errfp);
SUNDIALS_EXPORT int MRIStepSetUserData(void *arkode_mem,
                                       void *user_data);
SUNDIALS_EXPORT int MRIStepSetDiagnostics(void *arkode_mem,
                                          FILE *diagfp);

SUNDIALS_EXPORT int MRIStepSetPostprocessStepFn(void *arkode_mem,
                                                ARKPostProcessStepFn ProcessStep);


/*---------------------------------------------------------------
  MRIStepEvolve

  This integrates the ODE over an interval in t.

  MRIStepEvolve may be run in one of two modes (ARK_NORMAL or
  ARK_ONE_STEP), as determined by the itask argument:

  If itask is ARK_NORMAL, then the solver integrates from its
  current internal t value to a point at or beyond tout, then
  interpolates to t = tout and returns y(tout) in the user-
  allocated vector yout.  This interpolation is typically less
  accurate than the full time step solutions produced by the
  solver.  If the user wishes that this returned value have full
  method accuracy, they may issue a call to MRIStepSetStopTime
  before the call to MRIStepEvolve to specify a fixed stop time
  to end the time step and return to the user.  Once the
  integrator returns at a tstop time, any future testing for
  tstop is disabled (and can be reenabled only though a new call
  to MRIStepSetStopTime).

  If itask is ARK_ONE_STEP, then the solver takes one internal
  time step and returns in yout the value of y at the new internal
  time. In this case, tout is used only during the first call to
  MRIStepEvolve to determine the direction of integration and the
  rough scale of the t variable.  As with the ARK_NORMAL mode, a
  user may specify a specific stop time for output of this step,
  assuming that the requested step is smaller than the step taken
  by the method.

  The time reached by the solver is placed in (*tret). The
  user is responsible for allocating the memory for this value.

  arkode_mem is the pointer to MRIStep memory returned by
             MRIStepCreate.

  tout  is the next time at which a computed solution is desired.

  yout  is the computed solution vector. In ARK_NORMAL mode with no
        errors and no roots found, yout=y(tout).

  tret  is a pointer to a real location. MRIStepEvolve sets (*tret)
        to the time reached by the solver and returns yout=y(*tret).

  itask is ARK_NORMAL or ARK_ONE_STEP, as described above.

  Here is a brief description of each return value:

  ARK_SUCCESS:      MRIStepEvolve succeeded and no roots were found.

  ARK_ROOT_RETURN:  MRIStepEvolve succeeded, and found one or more roots.
                    If nrtfn > 1, call MRIStepGetRootInfo to see
                    which g_i were found to have a root at (*tret).

  ARK_TSTOP_RETURN: MRIStepEvolve succeeded and returned at tstop.

  ARK_MEM_NULL:     The arkode_mem argument was NULL.

  ARK_NO_MALLOC:    arkode_mem was not allocated.

  ARK_ILL_INPUT:    One of the inputs is illegal. This
                    includes the situation where a root of one of
                    the root functions was found both at t0 and very
                    near t0. In any case, the user should see the
                    printed error message for more details.

  ARK_INNERSTEP_FAIL:  The inner stepper or a call to an inner
                       stepper function failed. Use the function
                       MRIStepGetLastInnerStepFlag to get the
                       exact return flag.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int MRIStepEvolve(void *arkode_mem, realtype tout,
                                  N_Vector yout, realtype *tret,
                                  int itask);


/*---------------------------------------------------------------
  MRIStepGetDky

  This computes the kth derivative of the y function at time t,
  where tn-hu <= t <= tn, tn denotes the current internal time
  reached, and hu is the last internal step size successfully
  used by the solver. The user may request k=0, 1, ..., d, where
  d = min(5,q), with q the order of accuracy for the time
  integration method. The derivative vector is returned in dky.
  This vector must be allocated by the caller. It is only legal
  to call this function after a successful return from
  MRIStepEvolve.

  arkode_mem is the pointer to MRIStep memory returned by
             MRIStepCreate.

  t   is the time at which the kth derivative of y is evaluated.
      The legal range for t is [tn-hu,tn] as described above.

  k   is the order of the derivative of y to be computed. The
      legal range for k is [0,min(q,5)] as described above.

  dky is the output derivative vector [((d/dy)^k)y](t).

  The return value for MRIStepGetDky is one of:

    ARK_SUCCESS:  MRIStepGetDky succeeded.

    ARK_BAD_K:    k is not in the range 0, 1, ..., s-1.

    ARK_BAD_T:    t is not in the interval [tn-hu,tn].

    ARK_BAD_DKY:  The dky argument was NULL.

    ARK_MEM_NULL: The arkode_mem argument was NULL.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int MRIStepGetDky(void *arkode_mem, realtype t,
                                  int k, N_Vector dky);


/*---------------------------------------------------------------
  Optional outputs from the MRIStep module:  the following
  functions can be called to get optional outputs and statistics
  related to the main integrator.

  MRIStepGetNumRhsEvals returns the number of calls to the user's
                        f functions

  MRIStepGetCurrentButcherTables returns the Butcher tables
                                 currently in use

  MRIStepGetWorkSpace returns the MRIStep real and integer workspaces

  MRIStepGetNumSteps returns the cumulative number of slow and fast
                     internal steps taken by the solver

  MRIStepGetLastStep returns the step size for the last internal
                     step

  MRIStepGetCurrentTime returns the current internal time reached
                        by the solver

  MRIStepGetNumGEvals returns the number of calls to the user's
                      g function (for rootfinding)

  MRIStepGetRootInfo returns the indices for which g_i was found to
                     have a root. The user must allocate space for
                     rootsfound. For i = 0 ... nrtfn-1,
                     rootsfound[i] = 1 if g_i has a root, and = 0
                     if not.

  MRIStepGetLastInnerStepFlag returns the last flag returned from
                              the inner stepper.

  The return value of MRIStepGet* is one of:
     ARK_SUCCESS   if successful
     ARK_MEM_NULL  if the MRIStep memory structure was NULL
     ARK_LMEM_NULL if a linear solver memory structure was NULL
  ---------------------------------------------------------------*/

SUNDIALS_EXPORT int MRIStepGetNumRhsEvals(void *arkode_mem,
                                          long int *nfs_evals,
                                          long int *nff_evals);
SUNDIALS_EXPORT int MRIStepGetCurrentButcherTables(void *arkode_mem,
                                                   ARKodeButcherTable *Bs,
                                                   ARKodeButcherTable *Bf);
SUNDIALS_EXPORT int MRIStepGetWorkSpace(void *arkode_mem,
                                        long int *lenrw,
                                        long int *leniw);
SUNDIALS_EXPORT int MRIStepGetNumSteps(void *arkode_mem,
                                       long int *nssteps, long int *nfsteps);
SUNDIALS_EXPORT int MRIStepGetLastStep(void *arkode_mem,
                                       realtype *hlast);
SUNDIALS_EXPORT int MRIStepGetCurrentTime(void *arkode_mem,
                                          realtype *tcur);
SUNDIALS_EXPORT int MRIStepGetNumGEvals(void *arkode_mem,
                                        long int *ngevals);
SUNDIALS_EXPORT int MRIStepGetRootInfo(void *arkode_mem,
                                       int *rootsfound);
SUNDIALS_EXPORT int MRIStepGetLastInnerStepFlag(void *arkode_mem, int *flag);

/*---------------------------------------------------------------
  The following function returns the name of the constant
  associated with a MRIStep return flag
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT char *MRIStepGetReturnFlagName(long int flag);

/*---------------------------------------------------------------
  Function : MRIStepWriteParameters
  -----------------------------------------------------------------
  MRIStepWriteParameters outputs all timestepper module parameters
  to the provided file pointer.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int MRIStepWriteParameters(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
  Function : MRIStepWriteButcher
  -----------------------------------------------------------------
  MRIStepWriteButcher outputs the Butcher tables to the
  provided file pointer.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int MRIStepWriteButcher(void *arkode_mem, FILE *fp);

/*---------------------------------------------------------------
  MRIStepFree

  This frees the problem memory arkode_mem allocated by
  MRIStepCreate. Its only argument is the pointer arkode_mem
  returned by MRIStepCreate.
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void MRIStepFree(void **arkode_mem);

/*---------------------------------------------------------------
  MRIStepPrintMem

  This routine outputs the memory from the MRIStep structure and
  the main ARKode infrastructure to a specified file pointer
  (useful when debugging).
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT void MRIStepPrintMem(void* arkode_mem, FILE* outfile);


#ifdef __cplusplus
}
#endif

#endif
