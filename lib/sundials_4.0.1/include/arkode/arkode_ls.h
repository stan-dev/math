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
 * Header file for ARKode's linear solver interface.
 *--------------------------------------------------------------*/

#ifndef _ARKLS_H
#define _ARKLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  ARKLS Constants
  ===============================================================*/

/* ARKLS return values */
#define ARKLS_SUCCESS           0
#define ARKLS_MEM_NULL         -1
#define ARKLS_LMEM_NULL        -2
#define ARKLS_ILL_INPUT        -3
#define ARKLS_MEM_FAIL         -4
#define ARKLS_PMEM_NULL        -5
#define ARKLS_MASSMEM_NULL     -6
#define ARKLS_JACFUNC_UNRECVR  -7
#define ARKLS_JACFUNC_RECVR    -8
#define ARKLS_MASSFUNC_UNRECVR -9
#define ARKLS_MASSFUNC_RECVR   -10
#define ARKLS_SUNMAT_FAIL      -11
#define ARKLS_SUNLS_FAIL       -12


/*===============================================================
  ARKLS user-supplied function prototypes
  ===============================================================*/

/*---------------------------------------------------------------
  Type: ARKLsJacFn

  A Jacobian approximation function Jac must be of type 
  ARKLsJacFn. Its parameters are:

  Jac is the SUNMatrix that will be loaded by a ARKLsJacFn 
      with an approximation to the Jacobian matrix 
      J = (df_i/dy_j) at the point (t,y). 

  t   is the current value of the independent variable.

  y   is the current value of the dependent variable vector,
      namely the predicted value of y(t).

  fy  is the vector f(t,y).

  user_data is a pointer to user data - the same as the user_data
      parameter passed to ARKodeSetUserData.

  tmp1, tmp2, and tmp3 are pointers to memory allocated for
      vectors of length N which can be used by a ARKLsJacFn
      as temporary storage or work space.

  A ARKLsJacFn should return 0 if successful, a positive 
  value if a recoverable error occurred, and a negative value if 
  an unrecoverable error occurred.

  NOTE: See the relevant SUNMatrix implementation header files
      and documentation for mechanisms to inquire about matrix 
      dimensions, and for efficient ways to set matrix entries.

  NOTE: If the user's Jacobian routine needs other quantities,   
      they are accessible as follows: hcur (the current stepsize)
      and ewt (the error weight vector) are accessible through   
      ARKodeGetCurrentStep and ARKodeGetErrWeights, respectively 
      (see arkode.h). The unit roundoff is available as 
      UNIT_ROUNDOFF defined in sundials_types.h.

---------------------------------------------------------------*/
typedef int (*ARKLsJacFn)(realtype t, N_Vector y, N_Vector fy, 
                          SUNMatrix Jac, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


/*---------------------------------------------------------------
  Type: ARKLsMassFn

  A mass matrix approximation function Mass must be of type 
  ARKLsMassFn. Its parameters are:

  t   is the current value of the independent variable.

  M   is the SUNMatrix that will be loaded by a ARKLsMassFn 
      with an approximation to the mass matrix.

  user_data is a pointer to user data - the same as the user_data
      parameter passed to ARKodeSetUserData.

  tmp1, tmp2, and tmp3 are pointers to memory allocated for
  vectors of length N which can be used by a ARKLsMassFn
  as temporary storage or work space.

  A ARKLsMassFn should return 0 if successful, and a 
  negative value if an unrecoverable error occurred.

  ---------------------------------------------------------------*/
typedef int (*ARKLsMassFn)(realtype t, SUNMatrix M, void *user_data, 
                           N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  

/*---------------------------------------------------------------
  Type: ARKLsPrecSetupFn

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
          parameter passed to the ARKodeSetUserData function.

  NOTE: If the user's preconditioner needs other quantities,
        they are accessible as follows: hcur (the current stepsize)
        and ewt (the error weight vector) are accessible through
        ARKodeGetCurrentStep and ARKodeGetErrWeights, respectively).
        The unit roundoff is available as UNIT_ROUNDOFF defined in
        sundials_types.h.

  Returned value:
  The value to be returned by the PrecSetup function is a flag
  indicating whether it was successful.  This value should be
    0   if successful,
    > 0 for a recoverable error (step will be retried),
    < 0 for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*ARKLsPrecSetupFn)(realtype t, N_Vector y, 
                                N_Vector fy, booleantype jok, 
                                booleantype *jcurPtr,
                                realtype gamma, void *user_data);
  

/*---------------------------------------------------------------
  Type: ARKLsPrecSolveFn

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
         through a call to the routine ARKodeGetErrWeights.

  lr     is an input flag indicating whether PrecSolve is to use
         the left preconditioner P1 or right preconditioner
         P2: lr = 1 means use P1, and lr = 2 means use P2.

  user_data  is a pointer to user data - the same as the user_data
         parameter passed to the ARKodeSetUserData function.

  Returned value:
  The value to be returned by the PrecSolve function is a flag
  indicating whether it was successful.  This value should be
    0 if successful,
    positive for a recoverable error (step will be retried),
    negative for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*ARKLsPrecSolveFn)(realtype t, N_Vector y, 
                                N_Vector fy, N_Vector r, 
                                N_Vector z, realtype gamma, 
                                realtype delta, int lr, 
                                void *user_data);


/*---------------------------------------------------------------
  Type: ARKLsJacTimesSetupFn

  The user-supplied Jacobian-times-vector product setup function 
  JacTimesSetup and the user-supplied Jacobian-times-vector 
  product function JTimes together must generate the product
  J*v for v, where J is the Jacobian df/dy, or an approximation 
  to it, and v is a given vector. It should return 0 if 
  successful a positive value for a recoverable error or a
  negative value for an unrecoverable failure.

  Each call to the JacTimesSetup function is preceded by a call 
  to the RhsFn fi with the same (t,y) arguments.  Thus the 
  JacTimesSetup function can use any auxiliary data that is 
  computed and saved by the fi function and made accessible to 
  JacTimesSetup.

  A function JacTimesSetup must have the prototype given below.
  Its parameters are as follows:

  t       is the current value of the independent variable.

  y       is the current value of the dependent variable vector,
          namely the predicted value of y(t).

  fy      is the vector f(t,y).

  user_data  is a pointer to user data - the same as the user_data
          parameter passed to the ARKodeSetUserData function.

  Returned value:
  The value to be returned by the JacTimesSetup function is a flag
  indicating whether it was successful.  This value should be
    0   if successful,
    > 0 for a recoverable error (step will be retried),
    < 0 for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*ARKLsJacTimesSetupFn)(realtype t, N_Vector y, 
                                    N_Vector fy, void *user_data);


/*---------------------------------------------------------------
  Type: ARKLsJacTimesVecFn

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
           parameter passed to the ARKodeSetUserData function.

  tmp      is a pointer to memory allocated for an N_Vector
           which can be used by Jtimes for work space.
  ---------------------------------------------------------------*/
typedef int (*ARKLsJacTimesVecFn)(N_Vector v, N_Vector Jv, 
                                  realtype t, N_Vector y, 
                                  N_Vector fy, void *user_data, 
                                  N_Vector tmp);


/*---------------------------------------------------------------
  Type: ARKLsMassTimesSetupFn

  The user-supplied mass matrix-times-vector product setup 
  function MassTimesSetup and the user-supplied mass 
  matrix-times-vector product function MTimes together must 
  generate the product M*v for v, where M is the mass matrix, or 
  an approximation to it, and v is a given vector. It should 
  return 0 if successful a positive value for a recoverable error 
  or a negative value for an unrecoverable failure.

  A function MassTimesSetup must have the prototype given below.
  Its parameters are as follows:

  t       is the current value of the independent variable.

  mtimes_data  is a pointer to user data - the same as the 
          parameter passed to the ARKLsSetMassTimes function.

  Returned value:
  The value to be returned by the MassTimesSetup function is a 
  flag indicating whether it was successful.  This value should be
    0   if successful,
    > 0 for a recoverable error (step will be retried),
    < 0 for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*ARKLsMassTimesSetupFn)(realtype t, void *mtimes_data);


/*---------------------------------------------------------------
  Type: ARKLsMassTimesVecFn

  The user-supplied function mtimes is to generate the product
  M*v for given v, where M is the mass matrix, or an 
  approximation to it, and v is a given vector. It should return 
  0 if successful or a negative value for an unrecoverable failure.

  A function mtimes must have the prototype given below. Its
  parameters are as follows:

  v        is the N_Vector to be multiplied by M.

  Mv       is the output N_Vector containing M*v.

  t        is the current value of the independent variable.

  mtimes_data  is a pointer to user data - the same as the 
           parameter passed to the ARKLsSetMassTimes function.
  ---------------------------------------------------------------*/
typedef int (*ARKLsMassTimesVecFn)(N_Vector v, N_Vector Mv, 
                                   realtype t, void *mtimes_data);


/*---------------------------------------------------------------
  Type: ARKLsMassPrecSetupFn

  The user-supplied mass matrix preconditioner setup function 
  MPrecSetup and the user-supplied mass matrix preconditioner solve 
  function PrecSolve together must define left and right 
  preconditoner matrices P1 and P2 (either of which may be 
  trivial), such that the product P1*P2 is an approximation to 
  the mass matrix M.  The solution of systems P z = r, with P = P1 
  or P2, is to be carried out by the PrecSolve function, and 
  MPrecSetup is to do any necessary setup operations.

  The user-supplied preconditioner setup function MPrecSetup
  is to evaluate and preprocess any mass-matrix-related data
  needed by the preconditioner solve function PrecSolve.

  A function MPrecSetup must have the prototype given below.
  Its parameters are as follows:

  t       is the current value of the independent variable.

  user_data  is a pointer to user data - the same as the user_data
          parameter passed to the ARKodeSetUserData function.

  Returned value:
  The value to be returned by the MPrecSetup function is a flag
  indicating whether it was successful.  This value should be
    0   if successful,
    < 0 for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*ARKLsMassPrecSetupFn)(realtype t, void *user_data);


/*---------------------------------------------------------------
  Type: ARKLsMassPrecSolveFn

  The user-supplied mass matrix preconditioner solve function 
  MPrecSolve is to solve a linear system P z = r in which the 
  matrix P is one of the preconditioner matrices P1 or P2, 
  depending on the type of preconditioning chosen.

  A function MPrecSolve must have the prototype given below.
  Its parameters are as follows:

  t      is the current value of the independent variable.

  r      is the right-hand side vector of the linear system.

  z      is the output vector computed by MPrecSolve.

  delta  is an input tolerance for use by PSolve if it uses
         an iterative method in its solution.  In that case,
         the residual vector Res = r - P z of the system
         should be made less than delta in weighted L2 norm,
         i.e., sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta.
         Note: the error weight vector ewt can be obtained
         through a call to the routine ARKodeGetErrWeights.

  lr     is an input flag indicating whether MPrecSolve is to use
         the left preconditioner P1 or right preconditioner
         P2: lr = 1 means use P1, and lr = 2 means use P2.

  user_data  is a pointer to user data - the same as the user_data
         parameter passed to the ARKodeSetUserData function.

  Returned value:
  The value to be returned by the MPrecSolve function is a flag
  indicating whether it was successful.  This value should be
    0 if successful,
    negative for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*ARKLsMassPrecSolveFn)(realtype t, N_Vector r, 
                                    N_Vector z, realtype delta, 
                                    int lr, void *user_data);
  
  
#ifdef __cplusplus
}
#endif

#endif
