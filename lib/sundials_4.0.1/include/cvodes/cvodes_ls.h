/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Radu Serban @ LLNL
 *-----------------------------------------------------------------
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
 *-----------------------------------------------------------------
 * This is the header file for CVODES' linear solver interface, 
 * CVSLS.
 *
 * Part I contains type definitions and functions for using CVSLS
 * on forward problems (IVP integration and/or FSA)
 *
 * Part II contains type definitions and functions for using CVSLS
 * on adjoint (backward) problems
 *-----------------------------------------------------------------*/

#ifndef _CVSLS_H
#define _CVSLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*=================================================================
  CVLS Constants
  =================================================================*/

#define CVLS_SUCCESS          0
#define CVLS_MEM_NULL        -1
#define CVLS_LMEM_NULL       -2
#define CVLS_ILL_INPUT       -3
#define CVLS_MEM_FAIL        -4
#define CVLS_PMEM_NULL       -5
#define CVLS_JACFUNC_UNRECVR -6
#define CVLS_JACFUNC_RECVR   -7
#define CVLS_SUNMAT_FAIL     -8
#define CVLS_SUNLS_FAIL      -9

/* Return values for the adjoint module */

#define CVLS_NO_ADJ          -101
#define CVLS_LMEMB_NULL      -102


/*=================================================================
  PART I - forward problems
  =================================================================*/


/*-----------------------------------------------------------------
  CVSLS user-supplied function prototypes
  -----------------------------------------------------------------*/
  
/*-----------------------------------------------------------------
  Type: CVLsJacFn
 
  A Jacobian approximation function Jac must be of type CVLsJacFn.
  Its parameters are:
 
  Jac is the SUNMatrix matrix that will be loaded by a CVLsJacFn 
  with an approximation to the Jacobian 
      matrix J = (df_i/dy_j) at the point (t,y). 
 
  t   is the current value of the independent variable.
 
  y   is the current value of the dependent variable vector,
      namely the predicted value of y(t).
 
  fy  is the vector f(t,y).
 
  user_data is a pointer to user data - the same as the user_data
      parameter passed to CVodeSetUserdata.
 
  tmp1, tmp2, and tmp3 are pointers to memory allocated for
  vectors of length N which can be used by a CVLsJacFn
  as temporary storage or work space.
 
  A CVLsJacFn should return 0 if successful, a positive 
  value if a recoverable error occurred, and a negative value if 
  an unrecoverable error occurred.
 
  NOTE: See the relevant SUNMatrix implementation header files and
      documentation for mechanisms to inquire about matrix 
      dimensions, and for efficient ways to set matrix entries.
                                                                 
  NOTE: If the user's Jacobian routine needs other quantities,   
      they are accessible as follows: hcur (the current stepsize)
      and ewt (the error weight vector) are accessible through   
      CVodeGetCurrentStep and CVodeGetErrWeights, respectively 
      (see cvode.h). The unit roundoff is available as 
      UNIT_ROUNDOFF defined in sundials_types.h.
  -----------------------------------------------------------------*/
typedef int (*CVLsJacFn)(realtype t, N_Vector y, N_Vector fy, 
                         SUNMatrix Jac, void *user_data,
                         N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  
/*-----------------------------------------------------------------
  Type : CVLsPrecSetupFn

  The user-supplied preconditioner setup function PrecSetup and
  the user-supplied preconditioner solve function PrecSolve
  together must define left and right preconditoner matrices
  P1 and P2 (either of which may be trivial), such that the
  product P1*P2 is an approximation to the Newton matrix
  M = I - gamma*J.  Here J is the system Jacobian J = df/dy,
  and gamma is a scalar proportional to the integration step
  size h.  The solution of systems P z = r, with P = P1 or P2,
  is to be carried out by the PrecSolve function, and PrecSetup
  is to do any necessary setup operations.

  The user-supplied preconditioner setup function PrecSetup
  is to evaluate and preprocess any Jacobian-related data
  needed by the preconditioner solve function PrecSolve.
  This might include forming a crude approximate Jacobian,
  and performing an LU factorization on the resulting
  approximation to M.  This function will not be called in
  advance of every call to PrecSolve, but instead will be called
  only as often as necessary to achieve convergence within the
  Inexact Newton iteration.  If the PrecSolve function needs no
  preparation, the PrecSetup function can be NULL.

  For greater efficiency, the PrecSetup function may save
  Jacobian-related data and reuse it, rather than generating it
  from scratch.  In this case, it should use the input flag jok
  to decide whether to recompute the data, and set the output
  flag *jcurPtr accordingly.

  Each call to the PrecSetup function is preceded by a call to
  the RhsFn f with the same (t,y) arguments.  Thus the PrecSetup
  function can use any auxiliary data that is computed and
  saved by the f function and made accessible to PrecSetup.

  A function PrecSetup must have the prototype given below.
  Its parameters are as follows:

  t       is the current value of the independent variable.

  y       is the current value of the dependent variable vector,
           namely the predicted value of y(t).

  fy      is the vector f(t,y).

  jok     is an input flag indicating whether Jacobian-related
          data needs to be recomputed, as follows:
            jok == SUNFALSE means recompute Jacobian-related data
                   from scratch.
            jok == SUNTRUE  means that Jacobian data, if saved from
                   the previous PrecSetup call, can be reused
                   (with the current value of gamma).
          A Precset call with jok == SUNTRUE can only occur after
          a call with jok == SUNFALSE.

  jcurPtr is a pointer to an output integer flag which is
          to be set by PrecSetup as follows:
          Set *jcurPtr = SUNTRUE if Jacobian data was recomputed.
          Set *jcurPtr = SUNFALSE if Jacobian data was not recomputed,
                         but saved data was reused.

  gamma   is the scalar appearing in the Newton matrix.

  user_data  is a pointer to user data - the same as the user_data
          parameter passed to the CVodeSetUserData function.

  NOTE: If the user's preconditioner needs other quantities,
        they are accessible as follows: hcur (the current stepsize)
        and ewt (the error weight vector) are accessible through
        CVodeGetCurrentStep and CVodeGetErrWeights, respectively).
        The unit roundoff is available as UNIT_ROUNDOFF defined in
        sundials_types.h.

  Returned value:
  The value to be returned by the PrecSetup function is a flag
  indicating whether it was successful.  This value should be
    0   if successful,
    > 0 for a recoverable error (step will be retried),
    < 0 for an unrecoverable error (integration is halted).
  -----------------------------------------------------------------*/
typedef int (*CVLsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
                               booleantype jok, booleantype *jcurPtr,
                               realtype gamma, void *user_data);

  
/*-----------------------------------------------------------------
  Type : CVLsPrecSolveFn

  The user-supplied preconditioner solve function PrecSolve
  is to solve a linear system P z = r in which the matrix P is
  one of the preconditioner matrices P1 or P2, depending on the
  type of preconditioning chosen.

  A function PrecSolve must have the prototype given below.
  Its parameters are as follows:

  t      is the current value of the independent variable.

  y      is the current value of the dependent variable vector.

  fy     is the vector f(t,y).

  r      is the right-hand side vector of the linear system.

  z      is the output vector computed by PrecSolve.

  gamma  is the scalar appearing in the Newton matrix.

  delta  is an input tolerance for use by PSolve if it uses
         an iterative method in its solution.  In that case,
         the residual vector Res = r - P z of the system
         should be made less than delta in weighted L2 norm,
         i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta.
         Note: the error weight vector ewt can be obtained
         through a call to the routine CVodeGetErrWeights.

  lr     is an input flag indicating whether PrecSolve is to use
         the left preconditioner P1 or right preconditioner
         P2: lr = 1 means use P1, and lr = 2 means use P2.

  user_data  is a pointer to user data - the same as the user_data
          parameter passed to the CVodeSetUserData function.

  Returned value:
  The value to be returned by the PrecSolve function is a flag
  indicating whether it was successful.  This value should be
    0 if successful,
    positive for a recoverable error (step will be retried),
    negative for an unrecoverable error (integration is halted).
  -----------------------------------------------------------------*/
typedef int (*CVLsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
                               N_Vector r, N_Vector z, realtype gamma,
                               realtype delta, int lr, void *user_data);

  
/*---------------------------------------------------------------
  Type: CVLsJacTimesSetupFn

  The user-supplied Jacobian-times-vector product setup function
  JacTimesSetup and the user-supplied Jacobian-times-vector
  product function JTimes together must generate the product
  J*v for v, where J is the Jacobian df/dy, or an approximation
  to it, and v is a given vector.

  Each call to the JacTimesSetup function is preceded by a call
  to the RhsFn fi with the same (t,y) arguments.  Thus the
  JacTimesSetup function can use any auxiliary data that is
  computed and saved by the f function and made accessible to
  JacTimesSetup.

  A function JacTimesSetup must have the prototype given below.
  Its parameters are as follows:

  t       is the current value of the independent variable.

  y       is the current value of the dependent variable vector,
           namely the predicted value of y(t).

  fy      is the vector f(t,y).

  user_data  is a pointer to user data - the same as the user_data
          parameter passed to the CVodeSetUserData function.

  Returned value:
  The value to be returned by the JacTimesSetup function is a flag
  indicating whether it was successful.  This value should be
    0   if successful,
    > 0 for a recoverable error (step will be retried),
    < 0 for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*CVLsJacTimesSetupFn)(realtype t, N_Vector y,
                                   N_Vector fy, void *user_data);


/*-----------------------------------------------------------------
  Type : CVLsJacTimesVecFn

  The user-supplied function jtimes is to generate the product
  J*v for given v, where J is the Jacobian df/dy, or an
  approximation to it, and v is a given vector. It should return
  0 if successful a positive value for a recoverable error or
  a negative value for an unrecoverable failure.

  A function jtimes must have the prototype given below. Its
  parameters are as follows:

    v        is the N_Vector to be multiplied by J.

    Jv       is the output N_Vector containing J*v.

    t        is the current value of the independent variable.

    y        is the current value of the dependent variable
             vector.

    fy       is the vector f(t,y).

    user_data   is a pointer to user data, the same as the user_data
             parameter passed to the CVodeSetUserData function.

    tmp      is a pointer to memory allocated for an N_Vector
             which can be used by Jtimes for work space.
  -----------------------------------------------------------------*/
typedef int (*CVLsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t,
                                 N_Vector y, N_Vector fy,
                                 void *user_data, N_Vector tmp);

  
/*-----------------------------------------------------------------
  CVSLS Exported functions
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Required inputs to the CVLS linear solver interface

  CVodeSetLinearSolver specifies the SUNLinearSolver object that 
  CVode should use.  This is required if CVode is solving a problem 
  with a nonlinear solver that requires an inner linear solver 
  (e.g. Newton and not fixed-point).  The 'LS' argument must be 
  non-NULL, but A can be NULL if the solver requires no SUNMatrix 
  object.

  The return value is one of:
     CVLS_SUCCESS   if successful
     CVLS_MEM_NULL  if the CVODE memory was NULL
     CVLS_MEM_FAIL  if a memory allocation request failed
     CVLS_ILL_INPUT if the linear solver memory was NULL
 ---------------------------------------------------------------*/
SUNDIALS_EXPORT int CVodeSetLinearSolver(void *cvode_mem,
                                         SUNLinearSolver LS,
                                         SUNMatrix A);


/*-----------------------------------------------------------------
  Optional inputs to the CVLS linear solver interface
  -----------------------------------------------------------------
  CVodeSetJacFn specifies the Jacobian approximation routine to
    be used when constructing J.  By default, a difference 
    quotient approximation is used for dense/band SUNMatrix 
    objects; for all other matrix types passed to 
    CVodeSetLinearSolver, this routine must be user-supplied).
 
  CVodeSetMaxStepsBetweenJac specifies the maximum number of time
    steps to wait before recomputation of the Jacobian or 
    recommendation to update the preconditioner.  This differs from 
    the CVodeSetMaxStepsBetweenLSet, which merely indicates the 
    frequency with which the linear solver setup routine is called.  
    Default value is 50.

  CVodeSetEpsLin specifies the factor by which the tolerance on
    the nonlinear iteration is multiplied to get a tolerance on 
    the linear iteration.  Default value is 0.05.

  CVodeSetPreconditioner specifies the pset and psolve
    functions for user-supplied preconditioning.  Default is NULL 
    for both arguments (no preconditioning)

  CVodeSetJacTimes specifies the jtsetup and jtimes functions.
    Default is to use an internal finite difference approximation 
    routine with no extra jtsetup.

  The return value of these routines is one of:
     CVLS_SUCCESS   if successful
     CVLS_MEM_NULL  if the cvode memory was NULL
     CVLS_LMEM_NULL if the linear solver memory was NULL
     CVLS_ILL_INPUT if an input has an illegal value
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT int CVodeSetJacFn(void *cvode_mem, CVLsJacFn jac);
SUNDIALS_EXPORT int CVodeSetMaxStepsBetweenJac(void *cvode_mem,
                                               long int msbj);
SUNDIALS_EXPORT int CVodeSetEpsLin(void *cvode_mem, realtype eplifac);
SUNDIALS_EXPORT int CVodeSetPreconditioner(void *cvode_mem,
                                           CVLsPrecSetupFn pset,
                                           CVLsPrecSolveFn psolve);
SUNDIALS_EXPORT int CVodeSetJacTimes(void *cvode_mem,
                                     CVLsJacTimesSetupFn jtsetup,
                                     CVLsJacTimesVecFn jtimes);

/*-----------------------------------------------------------------
  Optional outputs from the CVLS linear solver interface
  -----------------------------------------------------------------
  CVodeGetLinWorkSpace returns the real and integer workspace used
     by the CVLS module.

  CVodeGetNumJacEvals returns the number of calls made to the
     Jacobian evaluation routine jac.

  CVodeGetNumPrecEvals returns the number of preconditioner
     evaluations, i.e. the number of calls made to PrecSetup with 
     jok==SUNFALSE.

  CVodeGetNumPrecSolves returns the number of calls made to
     PrecSolve.

  CVodeGetNumLinIters returns the number of linear iterations.

  CVodeGetNumLinConvFails returns the number of linear
     convergence failures.

  CVodeGetNumJTSetupEvals returns the number of calls to jtsetup.

  CVodeGetNumJtimesEvals returns the number of calls to jtimes.

  CVodeGetNumLinRhsEvals returns the number of calls to the user
     f routine due to finite difference Jacobian times vector 
     evaluations or Jacobian matrix approximations.

  CVodeGetLastLinFlag returns the last error flag set by any of
     the CVLS interface functions.

  The return value of these routines is one of:
     CVLS_SUCCESS   if successful
     CVLS_MEM_NULL  if the cvode memory was NULL
     CVLS_LMEM_NULL if the linear solver memory was NULL
 -----------------------------------------------------------------*/

SUNDIALS_EXPORT int CVodeGetLinWorkSpace(void *cvode_mem,
                                         long int *lenrwLS,
                                         long int *leniwLS);
SUNDIALS_EXPORT int CVodeGetNumJacEvals(void *cvode_mem,
                                        long int *njevals);
SUNDIALS_EXPORT int CVodeGetNumPrecEvals(void *cvode_mem,
                                         long int *npevals);
SUNDIALS_EXPORT int CVodeGetNumPrecSolves(void *cvode_mem,
                                          long int *npsolves);
SUNDIALS_EXPORT int CVodeGetNumLinIters(void *cvode_mem,
                                        long int *nliters);
SUNDIALS_EXPORT int CVodeGetNumLinConvFails(void *cvode_mem,
                                            long int *nlcfails);
SUNDIALS_EXPORT int CVodeGetNumJTSetupEvals(void *cvode_mem,
                                              long int *njtsetups);
SUNDIALS_EXPORT int CVodeGetNumJtimesEvals(void *cvode_mem,
                                           long int *njvevals);
SUNDIALS_EXPORT int CVodeGetNumLinRhsEvals(void *cvode_mem,
                                           long int *nfevalsLS);
SUNDIALS_EXPORT int CVodeGetLastLinFlag(void *cvode_mem,
                                        long int *flag);

/*-----------------------------------------------------------------
  The following function returns the name of the constant
  associated with a CVLS return flag
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT char *CVodeGetLinReturnFlagName(long int flag);

  

/*=================================================================
  PART II - backward problems
  =================================================================*/

/*-----------------------------------------------------------------
  CVSLS user-supplied function prototypes
  -----------------------------------------------------------------*/
  
/*-----------------------------------------------------------------
  Type: CVLsJacFnB
  -----------------------------------------------------------------
  A Jacobian approximation function jacB for the adjoint
  (backward) problem must have the prototype given below. 
  -----------------------------------------------------------------*/
typedef int (*CVLsJacFnB)(realtype t, N_Vector y, N_Vector yB,
                          N_Vector fyB, SUNMatrix JB,
                          void *user_dataB, N_Vector tmp1B,
                          N_Vector tmp2B, N_Vector tmp3B);

/*-----------------------------------------------------------------
  Type: CVLsJacFnBS
  -----------------------------------------------------------------
  A Jacobian approximation function jacBS for the adjoint
  (backward) problem, sensitivity-dependent case,  must have the
  prototype given below. 
  -----------------------------------------------------------------*/
typedef int (*CVLsJacFnBS)(realtype t, N_Vector y, N_Vector *yS,
                           N_Vector yB, N_Vector fyB, SUNMatrix JB,
                           void *user_dataB, N_Vector tmp1B,
                           N_Vector tmp2B, N_Vector tmp3B);
  
/*-----------------------------------------------------------------
  Type : CVLsPrecSetupFnB
  -----------------------------------------------------------------
  A function PrecSetupB for the adjoint (backward) problem must 
  have the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*CVLsPrecSetupFnB)(realtype t, N_Vector y, N_Vector yB,
                                N_Vector fyB, booleantype jokB,
                                booleantype *jcurPtrB,
                                realtype gammaB, void *user_dataB);


/*----------------------------------------------------------------
  Type : CVLsPrecSetupFnBS
  -----------------------------------------------------------------
  A function PrecSetupBS for the adjoint (backward) problem must 
  have the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*CVLsPrecSetupFnBS)(realtype t, N_Vector y,
                                 N_Vector *yS, N_Vector yB,
                                 N_Vector fyB, booleantype jokB,
                                 booleantype *jcurPtrB,
                                 realtype gammaB, void *user_dataB);


/*-----------------------------------------------------------------
  Type : CVLsPrecSolveFnB
  -----------------------------------------------------------------
  A function PrecSolveB for the adjoint (backward) problem  must 
  have the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*CVLsPrecSolveFnB)(realtype t, N_Vector y, N_Vector yB, 
                                N_Vector fyB, N_Vector rB, 
                                N_Vector zB, realtype gammaB,
                                realtype deltaB, int lrB,
                                void *user_dataB);

/*-----------------------------------------------------------------
  Type : CVLsPrecSolveFnBS
  -----------------------------------------------------------------
  A function PrecSolveBS for the adjoint (backward) problem  must 
  have the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*CVLsPrecSolveFnBS)(realtype t, N_Vector y, N_Vector *yS,
                                 N_Vector yB, N_Vector fyB,
                                 N_Vector rB, N_Vector zB,
                                 realtype gammaB, realtype deltaB,
                                 int lrB, void *user_dataB);

/*-----------------------------------------------------------------
  Type : CVLsJacTimesSetupFnB
  -----------------------------------------------------------------
  A function jtsetupB for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*CVLsJacTimesSetupFnB)(realtype t, N_Vector y, N_Vector yB,
                                    N_Vector fyB, void *jac_dataB);

/*-----------------------------------------------------------------
  Type : CVLsJacTimesSetupFnBS
  -----------------------------------------------------------------
  A function jtsetupBS for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*CVLsJacTimesSetupFnBS)(realtype t, N_Vector y,
                                     N_Vector *yS, N_Vector yB,
                                     N_Vector fyB, void *jac_dataB);

/*-----------------------------------------------------------------
  Type : CVLsJacTimesVecFnB
  -----------------------------------------------------------------
  A function jtimesB for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*CVLsJacTimesVecFnB)(N_Vector vB, N_Vector JvB, realtype t,
                                  N_Vector y, N_Vector yB, N_Vector fyB,
                                  void *jac_dataB, N_Vector tmpB);

/*-----------------------------------------------------------------
  Type : CVLsJacTimesVecFnBS
  -----------------------------------------------------------------
  A function jtimesBS for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*CVLsJacTimesVecFnBS)(N_Vector vB, N_Vector JvB,
                                   realtype t, N_Vector y, N_Vector *yS,
                                   N_Vector yB, N_Vector fyB,
                                   void *jac_dataB, N_Vector tmpB);

/*-----------------------------------------------------------------
  CVSLS Exported functions
  -----------------------------------------------------------------*/

/*---------------------------------------------------------------
  Required input for the CVSLS linear solver interface:

  CVodeSetLinearSolverB specifies the SUNLinearSolver object that 
  should be used for the backwards integration.  The 'which' 
  argument is the int returned by CVodeCreateB.

  The return value is one of:
     CVLS_SUCCESS   if successful
     CVLS_MEM_NULL  if the cvode memory was NULL
     CVLS_ILL_INPUT if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int CVodeSetLinearSolverB(void *cvode_mem,
                                          int which,
                                          SUNLinearSolver LS,
                                          SUNMatrix A);

/*-----------------------------------------------------------------
  Each CVodeSet***B or CVodeSet***BS function below links the
  main CVODES integrator with the corresponding CVSLS
  optional input function for the backward integration.
  The 'which' argument is the int returned by CVodeCreateB.
  -----------------------------------------------------------------*/

SUNDIALS_EXPORT int CVodeSetJacFnB(void *cvode_mem, int which,
                                   CVLsJacFnB jacB);
SUNDIALS_EXPORT int CVodeSetJacFnBS(void *cvode_mem, int which,
                                    CVLsJacFnBS jacBS);
  
SUNDIALS_EXPORT int CVodeSetEpsLinB(void *cvode_mem, int which,
                                    realtype eplifacB);

SUNDIALS_EXPORT int CVodeSetPreconditionerB(void *cvode_mem, int which, 
                                            CVLsPrecSetupFnB psetB,
                                            CVLsPrecSolveFnB psolveB);
SUNDIALS_EXPORT int CVodeSetPreconditionerBS(void *cvode_mem, int which, 
                                             CVLsPrecSetupFnBS psetBS,
                                             CVLsPrecSolveFnBS psolveBS);

SUNDIALS_EXPORT int CVodeSetJacTimesB(void *cvode_mem, int which, 
                                      CVLsJacTimesSetupFnB jtsetupB,
                                      CVLsJacTimesVecFnB jtimesB);
SUNDIALS_EXPORT int CVodeSetJacTimesBS(void *cvode_mem, int which, 
                                       CVLsJacTimesSetupFnBS jtsetupBS,
                                       CVLsJacTimesVecFnBS jtimesBS);

#ifdef __cplusplus
}
#endif

#endif
