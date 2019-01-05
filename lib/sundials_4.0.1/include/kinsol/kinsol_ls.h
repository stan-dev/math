/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Scott Cohen, Alan Hindmarsh, Radu Serban, 
 *                  and Aaron Collier @ LLNL
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
 * This is the header file for KINSOL's linear solver interface.
 *-----------------------------------------------------------------*/

#ifndef _KINLS_H
#define _KINLS_H

#include <sundials/sundials_direct.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*==================================================================
  KINLS return values
  ==================================================================*/

#define KINLS_SUCCESS      0

#define KINLS_MEM_NULL    -1
#define KINLS_LMEM_NULL   -2
#define KINLS_ILL_INPUT   -3
#define KINLS_MEM_FAIL    -4
#define KINLS_PMEM_NULL   -5
#define KINLS_JACFUNC_ERR -6
#define KINLS_SUNMAT_FAIL -7
#define KINLS_SUNLS_FAIL  -8

  
/*===============================================================
  KINLS user-supplied function prototypes
  ===============================================================*/

/*-----------------------------------------------------------------
  Type: KINLsJacFn
 
  A Jacobian approximation function Jac must be of type KINLsJacFn.
  Its parameters are:
 
  u        - current iterate (unscaled) [input]
 
  fu       - vector (type N_Vector) containing result of nonlinear
             system function evaluated at current iterate:
             fu = F(u) [input]
 
  J        - SUNMatrix that will be loaded by a KINLsJacFn with
             an approximation to the Jacobian matrix 
             J = (dF_i/dy_j).
 
  user_data  - pointer to user data - the same as the user_data
               parameter passed to KINSetFdata.
 
  tmp1, tmp2 - available scratch vectors (volatile storage)
 
  A KINLsJacFn should return 0 if successful or a non-zero value
  otherwise.
 
  NOTE: See the relevant SUNMatrix implementation header files and
      documentation for mechanisms to inquire about matrix 
      dimensions, and for efficient ways to set matrix entries.
  -----------------------------------------------------------------*/
typedef int (*KINLsJacFn)(N_Vector u, N_Vector fu, SUNMatrix J,
                          void *user_data, N_Vector tmp1, N_Vector tmp2);
  
/*------------------------------------------------------------------
  Type : KINLsPrecSetupFn

  The user-supplied preconditioner setup subroutine should compute
  the right-preconditioner matrix P (stored in memory block
  referenced by pdata pointer) used to form the scaled
  preconditioned linear system:
  
    (Df*J(uu)*(P^-1)*(Du^-1)) * (Du*P*x) = Df*(-F(uu))
  
  where Du and Df denote the diagonal scaling matrices whose
  diagonal elements are stored in the vectors uscale and fscale,
  repsectively.
 
  The preconditioner setup routine (referenced by iterative linear
  solver modules via pset (type KINLsPrecSetupFn)) will not be
  called prior to every call made to the psolve function, but will
  instead be called only as often as necessary to achieve
  convergence of the Newton iteration.
 
  Note: If the psolve routine requires no preparation, then a
  preconditioner setup function need not be given.
 
    uu  current iterate (unscaled) [input]
 
    uscale  vector (type N_Vector) containing diagonal elements
            of scaling matrix for vector uu [input]
 
    fval  vector (type N_Vector) containing result of nonliear
          system function evaluated at current iterate:
          fval = F(uu) [input]
 
    fscale  vector (type N_Vector) containing diagonal elements
            of scaling matrix for fval [input]
 
    user_data  pointer to user-allocated data memory block
 
  If successful, the function should return 0 (zero). If an error
  occurs, then the routine should return a non-zero integer value.
  -----------------------------------------------------------------*/
typedef int (*KINLsPrecSetupFn)(N_Vector uu, N_Vector uscale,
                                N_Vector fval, N_Vector fscale,
                                void *user_data);


/*------------------------------------------------------------------
  Type : KINLsPrecSolveFn

  The user-supplied preconditioner solve subroutine (referenced
  by iterative linear solver modules via psolve (type
  KINLsPrecSolveFn)) should solve a (scaled) preconditioned
  linear system of the generic form P*z = r, where P denotes the
  right-preconditioner matrix computed by the pset routine.
  
   uu  current iterate (unscaled) [input]
 
   uscale  vector (type N_Vector) containing diagonal elements
           of scaling matrix for vector uu [input]
 
   fval  vector (type N_Vector) containing result of nonliear
         system function evaluated at current iterate:
         fval = F(uu) [input]
 
   fscale  vector (type N_Vector) containing diagonal elements
           of scaling matrix for fval [input]
 
   vv  vector initially set to the right-hand side vector r, but
       which upon return contains a solution of the linear system
       P*z = r [input/output]
 
   user_data  pointer to user-allocated data memory block
 
  If successful, the function should return 0 (zero). If a
  recoverable error occurs, then the subroutine should return
  a positive integer value (in this case, KINSOL attempts to
  correct by calling the preconditioner setup function if the 
  preconditioner information is out of date). If an unrecoverable 
  error occurs, then the preconditioner solve function should return 
  a negative integer value.
  ------------------------------------------------------------------*/
typedef int (*KINLsPrecSolveFn)(N_Vector uu, N_Vector uscale, 
                                N_Vector fval, N_Vector fscale, 
                                N_Vector vv, void *user_data);


/*------------------------------------------------------------------
  Type : KINLsJacTimesVecFn
  
  The (optional) user-supplied matrix-vector product subroutine
  (referenced internally via jtimes (type KINLsJacTimesVecFn))
  is used to compute Jv = J(uu)*v (system Jacobian applied to a
  given vector). If a user-defined routine is not given, then the
  private routine is used.
 
   v  unscaled variant of vector to be multiplied by J(uu) [input]
 
   Jv  vector containing result of matrix-vector product J(uu)*v
       [output]
 
   uu  current iterate (unscaled) [input]
 
   new_uu  flag (reset by user) indicating if the iterate uu
           has been updated in the interim - Jacobian needs
           to be updated/reevaluated, if appropriate, unless
           new_uu = SUNFALSE [input/output]
 
   user_data  pointer to user data, the same as the user_data
              parameter passed to the KINSetUserData function.
 
  If successful, the function should return 0 (zero). If an error
  occurs, then the routine should return a non-zero integer value.
  ------------------------------------------------------------------*/
typedef int (*KINLsJacTimesVecFn)(N_Vector v, N_Vector Jv, N_Vector uu, 
                                  booleantype *new_uu, void *J_data);


  
/*==================================================================
  KINLS Exported functions
  ==================================================================*/

/*------------------------------------------------------------------
  Required inputs to the KINLS linear solver interface
 
  KINSetLinearSolver specifies the SUNLinearSolver object that  
  KINSOL should use.  This is required if KINSOL is solving a 
  problem with the Newton or Picard nonlinear solvers (i.e. not the 
  fixed-point or accelerated fixed-point solvers).
 
  The return value is one of:
     KINLS_SUCCESS   if successful
     KINLS_MEM_NULL  if the KINSOL memory was NULL
     KINLS_MEM_FAIL  if a memory allocation request failed
     KINLS_ILL_INPUT if the linear solver memory was NULL
  ------------------------------------------------------------------*/
SUNDIALS_EXPORT int KINSetLinearSolver(void *kinmem, SUNLinearSolver LS,
                                       SUNMatrix A);

/*------------------------------------------------------------------
  Optional inputs to the KINLS linear solver interface
  ------------------------------------------------------------------

  KINSetJacFn specifies the Jacobian approximation routine to be 
    used for a direct linear solver.  By default, a difference 
    quotient approximation is used when constructing J.  By 
    default, a difference quotient approximation is used for 
    dense/band SUNMatrix objects; for all other matrix types passed 
    to KINSetLinearSolver, this routine must be user-supplied).
 
  KINSetPreconditioner specifies the psetup and psolve functions
    for user-supplied preconditioning.  Here, the psetup routine 
    should compute a preconditioner matrix for the given linear 
    system, and the psolve routine is used to apply the 
    preconditioner to the linear system

  KINSetJacTimesVecFn specifies the user-supplied subroutine for 
    computing the matrix-vector product J(u)*v, where J denotes 
    the system Jacobian.

  The return value is one of:
     KINLS_SUCCESS   if successful
     KINLS_MEM_NULL  if the KINSOL memory was NULL
     KINLS_LMEM_NULL if the linear solver memory was NULL
     KINLS_ILL_INPUT if the input has an illegal value
  ------------------------------------------------------------------*/
SUNDIALS_EXPORT int KINSetJacFn(void *kinmem, KINLsJacFn jac);
SUNDIALS_EXPORT int KINSetPreconditioner(void *kinmem,
                                         KINLsPrecSetupFn psetup,
                                         KINLsPrecSolveFn psolve);
SUNDIALS_EXPORT int KINSetJacTimesVecFn(void *kinmem,
                                        KINLsJacTimesVecFn jtv);

/*------------------------------------------------------------------
  Optional outputs from the KINLS linear solver interface
  ------------------------------------------------------------------

  KINGetLinWorkSpace returns the real and integer workspace used by
    the KINLS linear solver.  Here, the integer workspace size is 
    the total number of long int-sized blocks of memory allocated 
    for vector storage.  Similarly, the real workspace size is the 
    total number of realtype-sized blocks of memory allocated for 
    vector storage.

  KINGetNumJacEvals returns the number of calls made to the
    Jacobian evaluation routine.

  KINGetNumLinFuncEvals returns the number of calls to the user's F
    routine due to finite difference Jacobian or Jacobian-vector 
    product approximation.

  KINGetNumPrecEvals returns the total number of preconditioner 
    evaluations (number of calls made to the user-defined psetup
    routine)

  KINGetNumPrecSolves returns the total number of times the
    preconditioner was applied to linear system (number of calls 
    made to the user-supplied psolve function)

  KINGetNumLinIters returns the total number of linear iterations
    performed

  KINGetNumLinConvFails returns the total number of linear solver
    convergence failures

  KINGetNumJtimesEvals returns the total number of times the 
    matrix-vector product J(u)*v was computed (number of calls 
    made to the jtimes subroutine)

  KINGetLastLinFlag returns the last error flag set by any of the 
    KINLS interface functions.

  The return value of these routines is one of:
     KINLS_SUCCESS   if successful
     KINLS_MEM_NULL  if the KINSOL memory was NULL
     KINLS_LMEM_NULL if the linear solver memory was NULL
  ------------------------------------------------------------------*/
SUNDIALS_EXPORT int KINGetLinWorkSpace(void *kinmem,
                                       long int *lenrwLS,
                                       long int *leniwLS);
SUNDIALS_EXPORT int KINGetNumJacEvals(void *kinmem,
                                      long int *njevals);
SUNDIALS_EXPORT int KINGetNumLinFuncEvals(void *kinmem,
                                          long int *nfevals);
SUNDIALS_EXPORT int KINGetNumPrecEvals(void *kinmem,
                                       long int *npevals);
SUNDIALS_EXPORT int KINGetNumPrecSolves(void *kinmem,
                                        long int *npsolves);
SUNDIALS_EXPORT int KINGetNumLinIters(void *kinmem,
                                      long int *nliters);
SUNDIALS_EXPORT int KINGetNumLinConvFails(void *kinmem,
                                          long int *nlcfails);
SUNDIALS_EXPORT int KINGetNumJtimesEvals(void *kinmem,
                                         long int *njvevals);
SUNDIALS_EXPORT int KINGetNumLinFuncEvals(void *kinmem,
                                          long int *nfevals); 
SUNDIALS_EXPORT int KINGetLastLinFlag(void *kinmem,
                                      long int *flag);

/*-----------------------------------------------------------------
  The following function returns the name of the constant
  associated with a KINLS return flag
  -----------------------------------------------------------------*/
SUNDIALS_EXPORT char *KINGetLinReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
