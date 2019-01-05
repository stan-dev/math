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
 * Header file for the Scaled and Preconditioned Iterative Linear 
 * Solvers interface in IDAS.
 *-----------------------------------------------------------------*/

#ifndef _IDASLS_H
#define _IDASLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*===============================================================
  IDASLS Constants
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

/* Return values for the adjoint module */
#define IDALS_NO_ADJ          -101
#define IDALS_LMEMB_NULL      -102

/*===============================================================
  PART I - forward problems
  ===============================================================*/

/*---------------------------------------------------------------
  IDASLS user-supplied function prototypes
  ---------------------------------------------------------------*/
  
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

  
/*-----------------------------------------------------------------
  Type : IDALsPrecSetupFn

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
  preparation, the PrecSetup function can be NULL.
 
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
  -----------------------------------------------------------------*/
typedef int (*IDALsPrecSetupFn)(realtype tt, N_Vector yy,
                                N_Vector yp, N_Vector rr,
                                realtype c_j, void *user_data);

/*-----------------------------------------------------------------
  Type : IDALsPrecSolveFn

  The optional user-supplied function PrecSolve must compute a
  solution to the linear system P z = r, where P is the left
  preconditioner defined by the user.  If no preconditioning
  is desired, pass NULL for PrecSolve to IDASp*.
 
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
  -----------------------------------------------------------------*/
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


/*-----------------------------------------------------------------
  Type : IDALsJacTimesVecFn

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
  -----------------------------------------------------------------*/
typedef int (*IDALsJacTimesVecFn)(realtype tt, N_Vector yy, 
                                  N_Vector yp, N_Vector rr,
                                  N_Vector v, N_Vector Jv,
                                  realtype c_j, void *user_data,
                                  N_Vector tmp1, N_Vector tmp2);


/*---------------------------------------------------------------
  IDASLS Exported functions
  ---------------------------------------------------------------*/


/*---------------------------------------------------------------
  Required inputs for the IDASLS linear solver interface:

  IDASetLinearSolver specifies the SUNLinearSolver object that 
  IDAS should use.  

  The return value is one of:
     IDALS_SUCCESS   if successful
     IDALS_MEM_NULL  if the IDA memory was NULL
     IDALS_ILL_INPUT if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDASetLinearSolver(void *ida_mem, 
                                       SUNLinearSolver LS,
                                       SUNMatrix A);

  
/*-----------------------------------------------------------------
  Optional inputs to the IDASLS linear solver -- ALL of these 
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
                                                                 
  The return value of IDALsSet* is one of:
     IDALS_SUCCESS   if successful
     IDALS_MEM_NULL  if the IDAS memory was NULL
     IDALS_LMEM_NULL if the linear solver memory was NULL
     IDALS_ILL_INPUT if an input has an illegal value
  -----------------------------------------------------------------*/
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

/*-----------------------------------------------------------------
  Optional outputs from the IDASLS linear solver interface:
  ---------------------------------------------------------------
  IDAGetLinWorkSpace returns the real and integer workspace used 
    by IDASLS.                                                  

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
    the IDASLS interface functions.
                                                                 
  The return value of IDALsGet* is one of:
     IDALS_SUCCESS   if successful
     IDALS_MEM_NULL  if the IDAS memory was NULL
     IDALS_LMEM_NULL if the linear solver memory was NULL
  -----------------------------------------------------------------*/
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


  
/*=================================================================
  PART II - backward problems
  =================================================================*/

/*-----------------------------------------------------------------
  Type: IDALsJacFnB

  A Jacobian approximation function JacB for the adjoint (backward) 
  problem must have the prototype given below. 
  -----------------------------------------------------------------*/
typedef int (*IDALsJacFnB)(realtype tt, realtype c_jB, N_Vector yy, 
                           N_Vector yp, N_Vector yyB, N_Vector ypB,
                           N_Vector rrB, SUNMatrix JacB,
                           void *user_dataB, N_Vector tmp1B,
                           N_Vector tmp2B, N_Vector tmp3B);


/*-----------------------------------------------------------------
  Type: IDALsJacFnBS

  A Jacobian approximation function JacBS for the adjoint (backward) 
  problem, sensitivity-dependent case, must have the prototype given 
  below. 
  -----------------------------------------------------------------*/
typedef int (*IDALsJacFnBS)(realtype tt, realtype c_jB, N_Vector yy, 
                            N_Vector yp, N_Vector *yS, N_Vector *ypS,
                            N_Vector yyB, N_Vector ypB, N_Vector rrB,
                            SUNMatrix JacB, void *user_dataB, 
                            N_Vector tmp1B, N_Vector tmp2B,
                            N_Vector tmp3B);


/*-----------------------------------------------------------------
  Type : IDALsPrecSetupFnB

  A function PrecSetupB for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*IDALsPrecSetupFnB)(realtype tt, N_Vector yy, 
                                 N_Vector yp, N_Vector yyB, 
                                 N_Vector ypB, N_Vector rrB, 
                                 realtype c_jB, void *user_dataB);

/*-----------------------------------------------------------------
  Type : IDALsPrecSetupFnBS

  A function PrecSetupBS for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*IDALsPrecSetupFnBS)(realtype tt, N_Vector yy, 
                                  N_Vector yp, N_Vector *yyS, 
                                  N_Vector *ypS, N_Vector yyB, 
                                  N_Vector ypB, N_Vector rrB,
                                  realtype c_jB, void *user_dataB);

/*-----------------------------------------------------------------
  Type : IDALsPrecSolveFnB

  A function PrecSolveB for the adjoint (backward) problem  must 
  have the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*IDALsPrecSolveFnB)(realtype tt, N_Vector yy,  
                                 N_Vector yp, N_Vector yyB, 
                                 N_Vector ypB, N_Vector rrB, 
                                 N_Vector rvecB, N_Vector zvecB,
                                 realtype c_jB, realtype deltaB,
                                 void *user_dataB);

/*-----------------------------------------------------------------
  Type : IDALsPrecSolveFnBS

  A function PrecSolveBS for the adjoint (backward) problem  must 
  have the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*IDALsPrecSolveFnBS)(realtype tt, N_Vector yy, 
                                  N_Vector yp, N_Vector *yyS, 
                                  N_Vector *ypS, N_Vector yyB, 
                                  N_Vector ypB, N_Vector rrB,
                                  N_Vector rvecB, N_Vector zvecB,
                                  realtype c_jB, realtype deltaB,
                                  void *user_dataB);

/*-----------------------------------------------------------------
  Type : IDALsJacTimesSetupFnB

  A function jtsetupB for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*IDALsJacTimesSetupFnB)(realtype t, N_Vector yy,
                                     N_Vector yp, N_Vector yyB,
                                     N_Vector ypB, N_Vector rrB,
                                     realtype c_jB, void *user_dataB);

/*-----------------------------------------------------------------
  Type : IDALsJacTimesSetupFnBS

  A function jtsetupBS for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*IDALsJacTimesSetupFnBS)(realtype t, N_Vector yy, 
                                      N_Vector yp, N_Vector *yyS,
                                      N_Vector *ypS, N_Vector yyB,
                                      N_Vector ypB, N_Vector rrB,
                                      realtype c_jB, void *user_dataB);

/*-----------------------------------------------------------------
  Type : IDALsJacTimesVecFnB

  A function jtimesB for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*IDALsJacTimesVecFnB)(realtype t, N_Vector yy, 
                                   N_Vector yp, N_Vector yyB, 
                                   N_Vector ypB, N_Vector rrB,
                                   N_Vector vB, N_Vector JvB, 
                                   realtype c_jB, void *user_dataB, 
                                   N_Vector tmp1B, N_Vector tmp2B);

/*-----------------------------------------------------------------
  Type : IDALsJacTimesVecFnBS

  A function jtimesBS for the adjoint (backward) problem must have 
  the prototype given below.
  -----------------------------------------------------------------*/
typedef int (*IDALsJacTimesVecFnBS)(realtype t, N_Vector yy, 
                                    N_Vector yp, N_Vector *yyS, 
                                    N_Vector *ypS, N_Vector yyB, 
                                    N_Vector ypB, N_Vector rrB,
                                    N_Vector vB, N_Vector JvB, 
                                    realtype c_jB, void *user_dataB, 
                                    N_Vector tmp1B, N_Vector tmp2B);

  
/*---------------------------------------------------------------
  Required input for the IDASLS linear solver interface:

  IDALsSetLinearSolverB specifies the SUNLinearSolver 
  object that IDAS should use for the backwards integration.  The 
  'which' argument is the int returned by IDACreateB.

  The return value is one of:
     IDALS_SUCCESS   if successful
     IDALS_MEM_NULL  if the IDA memory was NULL
     IDALS_ILL_INPUT if the linear solver memory was NULL
  ---------------------------------------------------------------*/
SUNDIALS_EXPORT int IDASetLinearSolverB(void *ida_mem,
                                        int which,
                                        SUNLinearSolver LS,
                                        SUNMatrix A);

 
/*-----------------------------------------------------------------
  Optional Set Routines
  -----------------------------------------------------------------*/

/*-----------------------------------------------------------------
  Each IDASet***B or IDASet***BS function below links the
  main IDAS integrator with the corresponding IDASet***
  optional input function for the backward integration.
  The 'which' argument is the int returned by IDACreateB.
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT int IDASetJacFnB(void *ida_mem, int which,
                                 IDALsJacFnB jacB);
SUNDIALS_EXPORT int IDASetJacFnBS(void *ida_mem, int which,
                                  IDALsJacFnBS jacBS);
SUNDIALS_EXPORT int IDASetEpsLinB(void *ida_mem, int which,
                                  realtype eplifacB);
SUNDIALS_EXPORT int IDASetIncrementFactorB(void *ida_mem, int which, 
                                           realtype dqincfacB);
SUNDIALS_EXPORT int IDASetPreconditionerB(void *ida_mem, int which,
                                          IDALsPrecSetupFnB psetB,
                                          IDALsPrecSolveFnB psolveB);
SUNDIALS_EXPORT int IDASetPreconditionerBS(void *ida_mem, int which,
                                           IDALsPrecSetupFnBS psetBS,
                                           IDALsPrecSolveFnBS psolveBS);
SUNDIALS_EXPORT int IDASetJacTimesB(void *ida_mem, int which,
                                    IDALsJacTimesSetupFnB jtsetupB,
                                    IDALsJacTimesVecFnB jtimesB);
SUNDIALS_EXPORT int IDASetJacTimesBS(void *ida_mem, int which,
                                     IDALsJacTimesSetupFnBS jtsetupBS,
                                     IDALsJacTimesVecFnBS jtimesBS);

 
#ifdef __cplusplus
}
#endif

#endif
