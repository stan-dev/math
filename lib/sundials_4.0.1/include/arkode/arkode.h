/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and
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
 * This is the interface file for the main ARKode infrastructure.
 *---------------------------------------------------------------
 * ARKode is used to numerically solve the ordinary initial value
 * problems using one-step methods.  Users do not call ARKode
 * infrastructure routines directly; they instead interact with
 * one of the time stepping modules built on top of ARKode.
 * These time step modules define their supported problem types,
 * solver options, etc.
 *
 * This file serves to define constants and provide function
 * prototypes for use across ARKode-based time integration
 * modules.
 *--------------------------------------------------------------*/

#ifndef _ARKODE_H
#define _ARKODE_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
                         ARKode Constants
===============================================================*/

/*---------------------------------------------------------------
  Enumerations for inputs to *StepCreate and *StepEvolve.
  -----------------------------------------------------------------
  itask: The itask input parameter to ARKode indicates the job
         of the solver for the next user step. The ARK_NORMAL
         itask is to have the solver take internal steps until
         it has reached or just passed the user specified tout
         parameter. The solver then interpolates in order to
         return an approximate value of y(tout). The ARK_ONE_STEP
         option tells the solver to just take one internal step
         and return the solution at the point reached by that
         step.  If fixed time steps are desired, then the user
         can run in ARK_ONE_STEP mode, and set hmin and hmax to
         be the desired time step -- this will effectively
         disable all temporal error control, and the user may
         need to increase the iteration counts for the nonlinear
         and linear solvers to guarantee convergence for larger
         step sizes.
  ---------------------------------------------------------------*/

/* itask */
#define ARK_NORMAL         1
#define ARK_ONE_STEP       2


/* ARKode return flags */
#define ARK_SUCCESS               0
#define ARK_TSTOP_RETURN          1
#define ARK_ROOT_RETURN           2

#define ARK_WARNING              99

#define ARK_TOO_MUCH_WORK        -1
#define ARK_TOO_MUCH_ACC         -2
#define ARK_ERR_FAILURE          -3
#define ARK_CONV_FAILURE         -4

#define ARK_LINIT_FAIL           -5
#define ARK_LSETUP_FAIL          -6
#define ARK_LSOLVE_FAIL          -7
#define ARK_RHSFUNC_FAIL         -8
#define ARK_FIRST_RHSFUNC_ERR    -9
#define ARK_REPTD_RHSFUNC_ERR    -10
#define ARK_UNREC_RHSFUNC_ERR    -11
#define ARK_RTFUNC_FAIL          -12
#define ARK_LFREE_FAIL           -13
#define ARK_MASSINIT_FAIL        -14
#define ARK_MASSSETUP_FAIL       -15
#define ARK_MASSSOLVE_FAIL       -16
#define ARK_MASSFREE_FAIL        -17
#define ARK_MASSMULT_FAIL        -18

#define ARK_MEM_FAIL             -20
#define ARK_MEM_NULL             -21
#define ARK_ILL_INPUT            -22
#define ARK_NO_MALLOC            -23
#define ARK_BAD_K                -24
#define ARK_BAD_T                -25
#define ARK_BAD_DKY              -26
#define ARK_TOO_CLOSE            -27

#define ARK_POSTPROCESS_FAIL     -28
#define ARK_VECTOROP_ERR         -29

#define ARK_NLS_INIT_FAIL        -30
#define ARK_NLS_SETUP_FAIL       -31
#define ARK_NLS_SETUP_RECVR      -32
#define ARK_NLS_OP_ERR           -33

#define ARK_INNERSTEP_FAIL       -34

#define ARK_UNRECOGNIZED_ERROR   -99


/*===============================================================
                          FUNCTION TYPES
  ===============================================================*/

/*-----------------------------------------------------------------
  Type : ARKRhsFn
  -----------------------------------------------------------------
  The f functions which define the right hand side of ODE
  systems, f(t,y), must have type ARKRhsFn.  f takes as input the
  independent variable value t, and the dependent variable vector
  y.  It stores the result of f(t,y) in the vector ydot.  The y
  and ydot arguments are of type N_Vector.
  (Allocation of memory for ydot is handled within ARKode)
  The user_data parameter is the same as the user_data
  parameter set by the user through the *StepSetUserData routine.
  This user-supplied pointer is passed to the user's f()
  function every time it is called.

  An ARKRhsFn should return 0 if successful, a negative value if
  an unrecoverable error occured, and a positive value if a
  recoverable error (e.g. invalid y values) occured.
  If an unrecoverable occured, the integration is halted.
  If a recoverable error occured, then (in most cases) ARKode
  will try to correct and retry.
  ---------------------------------------------------------------*/
typedef int (*ARKRhsFn)(realtype t, N_Vector y,
                        N_Vector ydot, void *user_data);

/*-----------------------------------------------------------------
  Type : ARKRootFn
  -----------------------------------------------------------------
  A function g, which defines a set of functions g_i(t,y) whose
  roots are sought during the integration, must have type
  ARKRootFn. The function g takes as input the independent
  variable value t, and the dependent variable vector y.  It
  stores the nrtfn values g_i(t,y) in the realtype array gout.
  (Allocation of memory for gout is handled within ARKode.)
  The user_data parameter is the same as that passed by the user
  to the *StepSetUserData routine.  This user-supplied pointer
  is passed to the user's g function every time it is called.

  An ARKRootFn should return 0 if successful or a non-zero value
  if an error occured (in which case the integration will be
  halted).
  ---------------------------------------------------------------*/
typedef int (*ARKRootFn)(realtype t, N_Vector y,
                         realtype *gout, void *user_data);

/*-----------------------------------------------------------------
  Type : ARKEwtFn
  -----------------------------------------------------------------
  A function e, which sets the error weight vector ewt, must have
  type ARKEwtFn.  The function e takes as input the current
  dependent variable y. It must set the vector of error weights
  used in the WRMS norm:

    ||y||_WRMS = sqrt [ 1/N * sum ( ewt_i * y_i)^2 ]

  Typically, the vector ewt has components:

    ewt_i = 1 / (reltol * |y_i| + abstol_i)

  The user_data parameter is the same as that passed by the user
  to the *StepSetUserData routine.  This user-supplied pointer
  is passed to the user's e function every time it is called.
  An ARKEwtFn e must return 0 if the error weight vector has been
  successfuly set and a non-zero value otherwise.
  ---------------------------------------------------------------*/
typedef int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void *user_data);

/*-----------------------------------------------------------------
  Type : ARKRwtFn
  -----------------------------------------------------------------
  A function r, which sets the residual weight vector rwt, must
  have type ARKRwtFn.  The function r takes as input the current
  dependent variable y.  It must set the vector of residual
  weights used in the WRMS norm:

    ||v||_WRMS = sqrt [ 1/N * sum ( rwt_i * v_i)^2 ]

  Typically, the vector rwt has components:

    rwt_i = 1 / (reltol * |(M*y)_i| + rabstol_i)

  The user_data parameter is the same as that passed by the user
  to the *StepSetUserData routine.  This user-supplied pointer
  is passed to the user's r function every time it is called.
  An ARKRwtFn e must return 0 if the residual weight vector has
  been successfuly set and a non-zero value otherwise.
  ---------------------------------------------------------------*/
typedef int (*ARKRwtFn)(N_Vector y, N_Vector rwt, void *user_data);

/*-----------------------------------------------------------------
  Type : ARKErrHandlerFn
  -----------------------------------------------------------------
  A function eh, which handles error messages, must have type
  ARKErrHandlerFn.  The function eh takes as input the error code,
  the name of the module reporting the error, the error message,
  and a pointer to user data, the same as that passed to
  *StepSetUserData.

  All error codes are negative, except ARK_WARNING which indicates
  a warning (the solver continues).

  An ARKErrHandlerFn has no return value.
  ---------------------------------------------------------------*/
typedef void (*ARKErrHandlerFn)(int error_code, const char *module,
                                const char *function, char *msg,
                                void *user_data);

/*-----------------------------------------------------------------
  Type : ARKAdaptFn
  -----------------------------------------------------------------
  A function which sets the new time step h, must have type
  ARKAdaptFn.  The function takes as input the current dependent
  variable y, the current time t, the last 3 step sizes h,
  the last 3 error estimates, the method order q, the embedding 
  order p, and a pointer to user data. The function must set the 
  scalar step size for the upcoming time step.  This value will 
  subsequently be bounded by the user-supplied values for the 
  minimum and maximum allowed time step, and the time step 
  satisfying the explicit stability restriction.  The user_data 
  parameter is the same as that passed by the user to the 
  *StepSetUserData routine.  This user-supplied pointer is 
  passed to the function every time it is called.

  An ARKAdaptFn must return 0 if the new time step has been
  successfuly set and a non-zero value otherwise.
  ---------------------------------------------------------------*/
typedef int (*ARKAdaptFn)(N_Vector y, realtype t, realtype h1,
                          realtype h2, realtype h3,
                          realtype e1, realtype e2,
                          realtype e3, int q, int p,
                          realtype *hnew, void *user_data);

/*-----------------------------------------------------------------
  Type : ARKExpStabFn
  -----------------------------------------------------------------
  A function which returns the time step satisfying the stability
  restriction for an explicit portion of the ODE.  The function
  takes as input the current dependent variable y, the current
  time t, and a pointer to user data. The function must set the
  scalar step size satisfying the stability restriction for the
  upcoming time step.  This value will subsequently be bounded by
  the user-supplied values for the minimum and maximum allowed
  time step, and the accuracy-based time step.  The user_data
  parameter is the same as that passed by the user to the
  *StepSetUserData routine.  This user-supplied pointer is passed
  to the function every time it is called.

  If this function is not supplied (NULL), or if it returns a
  negative time step size, then ARKode will assume that there is
  no explicit stability restriction on the time step size.

  An ARKExpStabFn must return 0 if the step size limit has been
  successfuly set and a non-zero value otherwise.
  ---------------------------------------------------------------*/
typedef int (*ARKExpStabFn)(N_Vector y, realtype t,
                            realtype *hstab, void *user_data);

/*-----------------------------------------------------------------
  Type : ARKVecResizeFn
  -----------------------------------------------------------------
  When calling *StepResize, the user may specify a vector resize
  function to be used to convert any existing N_Vectors in the
  ARKode memory structure to the new problem size.  This would
  typically be used if there is a user-supplied N_Vector module
  that allows dynamic resizing of the vector data structures
  without the need to delete/allocate memory on each call.

  The default behavior will be to delete the vector memory and
  re-clone from the new vector; if this is the desired behavior
  then specification of the ARKVecResizeFn is not recommended.

  The first argument, 'y', is the vector to be resized.

  The second argument, 'ytemplate', is the user-provided vector
  with the "new" size, that may be used as a template.

  The third argument, 'user_data', is a user-provided data
  structure to *StepResize, in case additional data is
  necessary for the resize operation.

  An ARKVecResizeFn should return 0 if successful, and a nonzero
  value if an error occurred.
  ---------------------------------------------------------------*/
typedef int (*ARKVecResizeFn)(N_Vector y, N_Vector ytemplate,
                              void *user_data);

/*-----------------------------------------------------------------
  Type : ARKPostProcessStepFn
  -----------------------------------------------------------------
  A function that is used to process the results of each timestep
  solution, in preparation for subsequent steps.  A routine of
  this type is designed for tasks such as inter-processor
  communication, computation of derived quantities, etc..

  IF THIS IS USED TO MODIFY ANY OF THE ACTIVE STATE DATA, THEN ALL
  THEORETICAL GUARANTEES OF SOLUTION ACCURACY AND STABILITY ARE
  LOST.

  Inputs:
    t          current time of ARKode solution
    y          current ARKode solution N_Vector for processing
    user_data  the structure passed by the user to the
               ARKodeSetUserData routine.
  ---------------------------------------------------------------*/
typedef int (*ARKPostProcessStepFn)(realtype t, N_Vector y,
                                    void *user_data);


#ifdef __cplusplus
}
#endif

#endif
