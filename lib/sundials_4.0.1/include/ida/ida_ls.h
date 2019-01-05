/*----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *         Alan Hindmarsh, Radu Serban and Aaron Collier @ LLNL
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
 * This is the header file for IDA's linear solver interface.
 *-----------------------------------------------------------------*/

#ifndef _IDALS_H
#define _IDALS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  IDALS Constants
  ===============================================================*/

#define IDALS_SUCCESS           0
#define IDALS_MEM_NULL         -1 
#define IDALS_LMEM_NULL        -2 
#define IDALS_ILL_INPUT        -3
#define IDALS_MEM_FAIL         -4
#define IDALS_PMEM_NULL        -5
#define IDALS_JACFUNC_UNRECVR  -6
#define IDALS_JACFUNC_RECVR    -7
#define IDALS_SUNMAT_FAIL      -8
#define IDALS_SUNLS_FAIL       -9

  
/*===============================================================
  IDALS user-supplied function prototypes
  ===============================================================*/
  
/*-----------------------------------------------------------------
  Type: IDALsJacFn

  A Jacobian approximation function Jac must be of type IDALsJacFn.
  Its parameters are:                     
                                                                 
  t   is the current value of the independent variable.
                                                                 
  y   is the current value of the dependent variable vector,     
      namely the predicted value of y(t).                     
                                                                 
  yp  is the current value of the derivative vector y',          
      namely the predicted value of y'(t).                    
                                                                 
  r   is the residual vector F(tt,yy,yp).                     
                                                                 
  c_j is the scalar in the system Jacobian, proportional to 
      the inverse of the step size h.
                                                                 
  user_data is a pointer to user Jacobian data - the same as the    
      user_data parameter passed to IDASetUserData.                     
                                                                 
  Jac is the SUNMatrix to be loaded by an IDALsJacFn routine 
      with an approximation to the system Jacobian matrix                                  
             J = dF/dy + c_j *dF/dy'
      at the given point (t,y,y'), where the ODE system is given 
      by F(t,y,y') = 0.   

      Note that Jac is NOT preset to zero!
                                                                 
  tmp1, tmp2, tmp3 are pointers to memory allocated for          
      N_Vectors which can be used by an IDALsJacFn routine 
      as temporary storage or work space.                     
                                                                 
  A IDALsJacFn should return                                
      0 if successful,                                           
      a positive int if a recoverable error occurred, or         
      a negative int if a nonrecoverable error occurred.         
  In the case of a recoverable error return, the integrator will 
  attempt to recover by reducing the stepsize (which changes c_j).
 
  NOTE: See the relevant SUNMatrix implementation header files
      and documentation for mechanisms to inquire about matrix 
      dimensions, and for efficient ways to set matrix entries.
                                                                
  NOTE: If the user's Jacobian routine needs other quantities,   
      they are accessible as follows: hcur (the current stepsize)
      and ewt (the error weight vector) are accessible through   
      IDAGetCurrentStep and IDAGetErrWeights, respectively, but this
      requires including in user_data a pointer to the solver memory.
      The unit roundoff is available as UNIT_ROUNDOFF defined in
      sundials_types.h.
  -----------------------------------------------------------------*/
typedef int (*IDALsJacFn)(realtype t, realtype c_j, N_Vector y,
                          N_Vector yp, N_Vector r, SUNMatrix Jac,
                          void *user_data, N_Vector tmp1,
                          N_Vector tmp2, N_Vector tmp3);

  
/*---------------------------------------------------------------
  Type: IDALsPrecSetupFn

  The optional user-supplied functions PrecSetup and PrecSolve
  together must define the left preconditoner matrix P
  approximating the system Jacobian matrix
     J = dF/dy + c_j*dF/dy'
  (where the DAE system is F(t,y,y') = 0), and solve the linear
  systems P z = r.   PrecSetup is to do any necessary setup
  operations, and PrecSolve is to compute the solution of
  P z = r.
 
  The preconditioner setup function PrecSetup is to evaluate and
  preprocess any Jacobian-related data needed by the
  preconditioner solve function PrecSolve.  This might include
  forming a crude approximate Jacobian, and performing an LU
  factorization on it.  This function will not be called in
  advance of every call to PrecSolve, but instead will be called
  only as often as necessary to achieve convergence within the
  Newton iteration.  If the PrecSolve function needs no
  preparation, the PrecSetup function can be NULL when passed to
  IDALsSetPreconditioner().
 
  Each call to the PrecSetup function is preceded by a call to
  the system function res with the same (t,y,y') arguments.
  Thus the PrecSetup function can use any auxiliary data that is
  computed and saved by the res function and made accessible
  to PrecSetup.
 
  A preconditioner setup function PrecSetup must have the
  prototype given below.  Its parameters are as follows:
 
  tt  is the current value of the independent variable t.
 
  yy  is the current value of the dependent variable vector,
      namely the predicted value of y(t).
 
  yp  is the current value of the derivative vector y',
      namely the predicted value of y'(t).
 
  rr  is the current value of the residual vector F(t,y,y').
 
  c_j is the scalar in the system Jacobian, proportional to 1/hh.
 
  user_data is a pointer to user data, the same as the user_data
      parameter passed to IDASetUserData.
 
  NOTE: If the user's preconditioner needs other quantities,
      they are accessible as follows: hcur (the current stepsize)
      and ewt (the error weight vector) are accessible through
      IDAGetCurrentStep and IDAGetErrWeights, respectively (see
      ida.h). The unit roundoff is available as
      UNIT_ROUNDOFF defined in sundials_types.h
 
  The IDALsPrecSetupFn should return
      0 if successful,
      a positive int if a recoverable error occurred, or
      a negative int if a nonrecoverable error occurred.
  In the case of a recoverable error return, the integrator will
  attempt to recover by reducing the stepsize (which changes cj).
  ---------------------------------------------------------------*/
typedef int (*IDALsPrecSetupFn)(realtype tt, N_Vector yy,
                                N_Vector yp, N_Vector rr,
                                realtype c_j, void *user_data);

  
/*---------------------------------------------------------------
  Type: IDALsPrecSolveFn

  The optional user-supplied function PrecSolve must compute a
  solution to the linear system P z = r, where P is the left
  preconditioner defined by the user.  If no preconditioning
  is desired, pass NULL for PrecSolve to 
  IDALsSetPreconditioner() (or do not call that routine at 
  all).
 
  A preconditioner solve function PrecSolve must have the
  prototype given below.  Its parameters are as follows:
 
  tt is the current value of the independent variable t.
 
  yy is the current value of the dependent variable vector y.
 
  yp is the current value of the derivative vector y'.
 
  rr is the current value of the residual vector F(t,y,y').
 
  rvec is the input right-hand side vector r.
 
  zvec is the computed solution vector z.
 
  c_j is the scalar in the system Jacobian, proportional to 1/hh.
 
  delta is an input tolerance for use by PrecSolve if it uses an
      iterative method in its solution.   In that case, the
      the residual vector r - P z of the system should be
      made less than delta in weighted L2 norm, i.e.,
             sqrt [ Sum (Res[i]*ewt[i])^2 ] < delta .
      Note: the error weight vector ewt can be obtained
      through a call to the routine IDAGetErrWeights.
 
  user_data is a pointer to user data, the same as the user_data
      parameter passed to IDASetUserData.
 
  The IDALsPrecSolveFn should return
      0 if successful,
      a positive int if a recoverable error occurred, or
      a negative int if a nonrecoverable error occurred.
  Following a recoverable error, the integrator will attempt to
  recover by updating the preconditioner and/or reducing the
  stepsize.
  ---------------------------------------------------------------*/
typedef int (*IDALsPrecSolveFn)(realtype tt, N_Vector yy,
                                N_Vector yp, N_Vector rr,
                                N_Vector rvec, N_Vector zvec,
                                realtype c_j, realtype delta,
                                void *user_data);

  
/*---------------------------------------------------------------
 Type: IDALsJacTimesSetupFn

 The user-supplied Jacobian-times-vector product setup function 
 JacTimesSetup and the user-supplied Jacobian-times-vector 
 product function JTimes together must generate the product
 J*v for v, where J is the Jacobian matrix
     J = dF/dy + c_j*dF/dy'
 or an approximation to it, and v is a given vector. 

 Each call to the JacTimesSetup function is preceded by a call 
 to the residual res with the same (t,y) arguments.  Thus the 
 JacTimesSetup function can use any auxiliary data that is 
 computed and saved by the res function and made accessible to 
 JacTimesSetup.

 A function JacTimesSetup must have the prototype given below.
 Its parameters are as follows:

 t       is the current value of the independent variable.

 y       is the current value of the dependent variable vector,
          namely the predicted value of y(t).

 fy      is the vector f(t,y).

 user_data  is a pointer to user data - the same as the user_data
         parameter passed to the IDASetUserData function.

 Returned value:
 The value to be returned by the JacTimesSetup function is a flag
 indicating whether it was successful.  This value should be
   0   if successful,
   > 0 for a recoverable error (step will be retried),
   < 0 for an unrecoverable error (integration is halted).
  ---------------------------------------------------------------*/
typedef int (*IDALsJacTimesSetupFn)(realtype tt, N_Vector yy,
                                    N_Vector yp, N_Vector rr,
                                    realtype c_j, void *user_data);


/*---------------------------------------------------------------
  Type: IDALsJacTimesVecFn

  The user-supplied function jtimes is to generate the product
  J*v for given v, where J is the Jacobian matrix
     J = dF/dy + c_j*dF/dy'
   or an approximation to it, and v is a given vector.
  It should return 0 if successful and a nonzero int otherwise.
 
  A function jtimes must have the prototype given below. Its
  parameters are as follows:
 
    tt   is the current value of the independent variable.
 
    yy   is the current value of the dependent variable vector,
         namely the predicted value of y(t).
 
    yp   is the current value of the derivative vector y',
         namely the predicted value of y'(t).
 
    rr   is the current value of the residual vector F(t,y,y').
 
    v    is the N_Vector to be multiplied by J.
 
    Jv   is the output N_Vector containing J*v.
 
    c_j  is the scalar in the system Jacobian, proportional
         to 1/hh.
 
    user_data is a pointer to user data, the same as the
         pointer passed to IDASetUserData.
 
    tmp1, tmp2 are two N_Vectors which can be used by Jtimes for
          work space.
  ---------------------------------------------------------------*/
typedef int (*IDALsJacTimesVecFn)(realtype tt, N_Vector yy,
                                  N_Vector yp, N_Vector rr,
                                  N_Vector v, N_Vector Jv,
                                  realtype c_j, void *user_data,
                                  N_Vector tmp1, N_Vector tmp2);


/*===============================================================
  IDALS Exported functions
  ===============================================================*/

/*---------------------------------------------------------------
  Required inputs for the IDALS linear solver interface:

  IDASetLinearSolver specifies the iterative SUNLinearSolver 
  object that IDA should use.  

  The return value is one of:
     IDALS_SUCCESS   if successful
     IDALS_MEM_NULL  if the IDA memory was NULL
     IDALS_ILL_INPUT if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDASetLinearSolver(void *ida_mem, 
                                       SUNLinearSolver LS,
                                       SUNMatrix A);

  
/*---------------------------------------------------------------
  Optional inputs to the IDALS linear solver -- ALL of these 
  must be called AFTER the corresponding linear solver 
  object has been attached to IDA).
  -----------------------------------------------------------------
  IDASetJacFn specifies the dense Jacobian approximation
    routine to be used for a direct linear solver.  By default, 
    a difference quotient approximation is used for dense and 
    band; no default exists for sparse or user-supplied matrix 
    types (so this must be user-supplied).
 
  IDASetPreconditioner specifies the PrecSetup and PrecSolve
    functions.  Default is NULL for both arguments.

  IDASetJacTimes specifies the jtsetup and jtimes functions.
    Default is to use an internal finite difference approximation 
    routine for the jtimes (with no setup).

  IDASetEpsLin specifies the factor in the linear iteration
    convergence test constant.  Default is 0.05

  IDASetIncrementFactor specifies a factor in the increments
    to yy used in the difference quotient approximations to
    matrix-vector products Jv.  Default is 1.0
                                                                 
  The return value of these routines is one of:
     IDALS_SUCCESS   if successful
     IDALS_MEM_NULL  if the ida memory was NULL
     IDALS_LMEM_NULL if the linear solver memory was NULL
     IDALS_ILL_INPUT if an input has an illegal value
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDASetJacFn(void *ida_mem, IDALsJacFn jac);
SUNDIALS_EXPORT int IDASetPreconditioner(void *ida_mem,
                                         IDALsPrecSetupFn pset, 
                                         IDALsPrecSolveFn psolve);
SUNDIALS_EXPORT int IDASetJacTimes(void *ida_mem,
                                   IDALsJacTimesSetupFn jtsetup,
                                   IDALsJacTimesVecFn jtimes);
SUNDIALS_EXPORT int IDASetEpsLin(void *ida_mem, realtype eplifac);
SUNDIALS_EXPORT int IDASetIncrementFactor(void *ida_mem,
                                          realtype dqincfac);


/*---------------------------------------------------------------
  Optional outputs from the IDALS linear solver interface:
  ---------------------------------------------------------------
  IDAGetLinWorkSpace returns the real and integer workspace used 
    by IDALS.                                                  

  IDAGetNumJacEvals returns the number of calls made to the
    Jacobian evaluation routine jac.

  IDAGetNumPrecEvals returns the number of preconditioner   
    evaluations, i.e. the number of calls made to PrecSetup    
    with jok==SUNFALSE.                                           

  IDAGetNumPrecSolves returns the number of calls made to   
    PrecSolve.                                                 

  IDAGetNumLinIters returns the number of linear iterations.

  IDAGetNumLinConvFails returns the number of linear           
    convergence failures.                                      

  IDAGetNumJTSetupEvals returns the number of calls to jtsetup

  IDAGetNumJtimesEvals returns the number of calls to jtimes

  IDAGetNumLinResEvals returns the number of calls to the user 
    res routine due to finite difference Jacobian or 
    Jacobian-times-vector evaluation.                                                

  IDAGetLastLinFlag returns the last error flag set by any of
    the IDALS interface functions.
                                                                 
  The return value of these routines is one of:
     IDALS_SUCCESS   if successful
     IDALS_MEM_NULL  if the ida memory was NULL
     IDALS_LMEM_NULL if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDAGetLinWorkSpace(void *ida_mem,
                                       long int *lenrwLS,
                                       long int *leniwLS);
SUNDIALS_EXPORT int IDAGetNumJacEvals(void *ida_mem,
                                      long int *njevals);
SUNDIALS_EXPORT int IDAGetNumPrecEvals(void *ida_mem,
                                       long int *npevals);
SUNDIALS_EXPORT int IDAGetNumPrecSolves(void *ida_mem,
                                        long int *npsolves);
SUNDIALS_EXPORT int IDAGetNumLinIters(void *ida_mem,
                                      long int *nliters);
SUNDIALS_EXPORT int IDAGetNumLinConvFails(void *ida_mem,
                                          long int *nlcfails);
SUNDIALS_EXPORT int IDAGetNumJTSetupEvals(void *ida_mem,
                                          long int *njtsetups);
SUNDIALS_EXPORT int IDAGetNumJtimesEvals(void *ida_mem,
                                         long int *njvevals);
SUNDIALS_EXPORT int IDAGetNumLinResEvals(void *ida_mem,
                                         long int *nrevalsLS); 
SUNDIALS_EXPORT int IDAGetLastLinFlag(void *ida_mem,
                                      long int *flag);

/*---------------------------------------------------------------
  The following function returns the name of the constant 
  associated with an IDALS return flag
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT char *IDAGetLinReturnFlagName(long int flag);

#ifdef __cplusplus
}
#endif

#endif
