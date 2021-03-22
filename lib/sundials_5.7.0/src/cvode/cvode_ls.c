/* ----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Alan C. Hindmarsh and Radu Serban @ LLNL
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * Implementation file for CVode's linear solver interface.
 * ---------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cvode_impl.h"
#include "cvode_ls_impl.h"
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_band.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>

/* Private constants */
#define MIN_INC_MULT RCONST(1000.0)
#define MAX_DQITERS  3  /* max. number of attempts to recover in DQ J*v */
#define ZERO         RCONST(0.0)
#define PT25         RCONST(0.25)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/*=================================================================
  PRIVATE FUNCTION PROTOTYPES
  =================================================================*/

static int cvLsLinSys(realtype t, N_Vector y, N_Vector fy, SUNMatrix A,
                      booleantype jok, booleantype *jcur, realtype gamma,
                      void *user_data, N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3);

/*===============================================================
  CVLS Exported functions -- Required
  ===============================================================*/

/*---------------------------------------------------------------
  CVodeSetLinearSolver specifies the linear solver
  ---------------------------------------------------------------*/
int CVodeSetLinearSolver(void *cvode_mem, SUNLinearSolver LS,
                         SUNMatrix A)
{
  CVodeMem    cv_mem;
  CVLsMem     cvls_mem;
  int         retval, LSType;
  booleantype iterative;    /* is the solver iterative?    */
  booleantype matrixbased;  /* is a matrix structure used? */

  /* Return immediately if either cvode_mem or LS inputs are NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVLS_MEM_NULL, "CVLS",
                   "CVodeSetLinearSolver", MSG_LS_CVMEM_NULL);
    return(CVLS_MEM_NULL);
  }
  if (LS == NULL) {
    cvProcessError(NULL, CVLS_ILL_INPUT, "CVLS",
                   "CVodeSetLinearSolver",
                   "LS must be non-NULL");
    return(CVLS_ILL_INPUT);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if solver is compatible with LS interface */
  if ( (LS->ops->gettype == NULL) || (LS->ops->solve == NULL) ) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS",
                   "CVodeSetLinearSolver",
                   "LS object is missing a required operation");
    return(CVLS_ILL_INPUT);
  }

  /* Retrieve the LS type */
  LSType = SUNLinSolGetType(LS);

  /* Set flags based on LS type */
  iterative   = (LSType != SUNLINEARSOLVER_DIRECT);
  matrixbased = (LSType != SUNLINEARSOLVER_ITERATIVE);

  /* Test if vector is compatible with LS interface */
  if ( (cv_mem->cv_tempv->ops->nvconst == NULL) ||
       (cv_mem->cv_tempv->ops->nvwrmsnorm == NULL) ) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS",
                   "CVodeSetLinearSolver", MSG_LS_BAD_NVECTOR);
    return(CVLS_ILL_INPUT);
  }

  /* Check for compatible LS type, matrix and "atimes" support */
  if (iterative) {

    if (cv_mem->cv_tempv->ops->nvgetlength == NULL) {
      cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS",
                     "CVodeSetLinearSolver", MSG_LS_BAD_NVECTOR);
      return(CVLS_ILL_INPUT);
    }

    if (!matrixbased && LS->ops->setatimes == NULL) {
      cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "CVodeSetLinearSolver",
                     "Incompatible inputs: iterative LS must support ATimes routine");
      return(CVLS_ILL_INPUT);
    }

    if (matrixbased && A == NULL) {
      cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "CVodeSetLinearSolver",
                     "Incompatible inputs: matrix-iterative LS requires non-NULL matrix");
      return(CVLS_ILL_INPUT);
    }

  } else if (A == NULL) {

    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "CVodeSetLinearSolver",
                   "Incompatible inputs: direct LS requires non-NULL matrix");
    return(CVLS_ILL_INPUT);

  }

  /* free any existing system solver attached to CVode */
  if (cv_mem->cv_lfree)  cv_mem->cv_lfree(cv_mem);

  /* Set four main system linear solver function fields in cv_mem */
  cv_mem->cv_linit  = cvLsInitialize;
  cv_mem->cv_lsetup = cvLsSetup;
  cv_mem->cv_lsolve = cvLsSolve;
  cv_mem->cv_lfree  = cvLsFree;

  /* Allocate memory for CVLsMemRec */
  cvls_mem = NULL;
  cvls_mem = (CVLsMem) malloc(sizeof(struct CVLsMemRec));
  if (cvls_mem == NULL) {
    cvProcessError(cv_mem, CVLS_MEM_FAIL, "CVLS",
                    "CVodeSetLinearSolver", MSG_LS_MEM_FAIL);
    return(CVLS_MEM_FAIL);
  }
  memset(cvls_mem, 0, sizeof(struct CVLsMemRec));

  /* set SUNLinearSolver pointer */
  cvls_mem->LS = LS;

  /* Linear solver type information */
  cvls_mem->iterative   = iterative;
  cvls_mem->matrixbased = matrixbased;

  /* Set defaults for Jacobian-related fields */
  if (A != NULL) {
    cvls_mem->jacDQ  = SUNTRUE;
    cvls_mem->jac    = cvLsDQJac;
    cvls_mem->J_data = cv_mem;
  } else {
    cvls_mem->jacDQ  = SUNFALSE;
    cvls_mem->jac    = NULL;
    cvls_mem->J_data = NULL;
  }

  cvls_mem->jtimesDQ = SUNTRUE;
  cvls_mem->jtsetup  = NULL;
  cvls_mem->jtimes   = cvLsDQJtimes;
  cvls_mem->jt_f     = cv_mem->cv_f;
  cvls_mem->jt_data  = cv_mem;

  cvls_mem->user_linsys = SUNFALSE;
  cvls_mem->linsys      = cvLsLinSys;
  cvls_mem->A_data      = cv_mem;

  /* Set defaults for preconditioner-related fields */
  cvls_mem->pset   = NULL;
  cvls_mem->psolve = NULL;
  cvls_mem->pfree  = NULL;
  cvls_mem->P_data = cv_mem->cv_user_data;

  /* Initialize counters */
  cvLsInitializeCounters(cvls_mem);

  /* Set default values for the rest of the LS parameters */
  cvls_mem->msbj      = CVLS_MSBJ;
  cvls_mem->jbad      = SUNTRUE;
  cvls_mem->eplifac   = CVLS_EPLIN;
  cvls_mem->last_flag = CVLS_SUCCESS;

  /* If LS supports ATimes, attach CVLs routine */
  if (LS->ops->setatimes) {
    retval = SUNLinSolSetATimes(LS, cv_mem, cvLsATimes);
    if (retval != SUNLS_SUCCESS) {
      cvProcessError(cv_mem, CVLS_SUNLS_FAIL, "CVLS",
                     "CVodeSetLinearSolver",
                     "Error in calling SUNLinSolSetATimes");
      free(cvls_mem); cvls_mem = NULL;
      return(CVLS_SUNLS_FAIL);
    }
  }

  /* If LS supports preconditioning, initialize pset/psol to NULL */
  if (LS->ops->setpreconditioner) {
    retval = SUNLinSolSetPreconditioner(LS, cv_mem, NULL, NULL);
    if (retval != SUNLS_SUCCESS) {
      cvProcessError(cv_mem, CVLS_SUNLS_FAIL, "CVLS",
                     "CVodeSetLinearSolver",
                     "Error in calling SUNLinSolSetPreconditioner");
      free(cvls_mem); cvls_mem = NULL;
      return(CVLS_SUNLS_FAIL);
    }
  }

  /* When using a SUNMatrix object, store pointer to A and initialize savedJ */
  if (A != NULL) {
    cvls_mem->A = A;
    cvls_mem->savedJ = NULL; /* allocated in cvLsInitialize */
  }

  /* Allocate memory for ytemp and x */
  cvls_mem->ytemp = N_VClone(cv_mem->cv_tempv);
  if (cvls_mem->ytemp == NULL) {
    cvProcessError(cv_mem, CVLS_MEM_FAIL, "CVLS",
                    "CVodeSetLinearSolver", MSG_LS_MEM_FAIL);
    free(cvls_mem); cvls_mem = NULL;
    return(CVLS_MEM_FAIL);
  }

  cvls_mem->x = N_VClone(cv_mem->cv_tempv);
  if (cvls_mem->x == NULL) {
    cvProcessError(cv_mem, CVLS_MEM_FAIL, "CVLS",
                    "CVodeSetLinearSolver", MSG_LS_MEM_FAIL);
    N_VDestroy(cvls_mem->ytemp);
    free(cvls_mem); cvls_mem = NULL;
    return(CVLS_MEM_FAIL);
  }

  /* For iterative LS, compute default norm conversion factor */
  if (iterative)
    cvls_mem->nrmfac = SUNRsqrt( N_VGetLength(cvls_mem->ytemp) );

  /* Check if soltuion scaling should be enabled */
  if (matrixbased && cv_mem->cv_lmm == CV_BDF)
    cvls_mem->scalesol = SUNTRUE;
  else
    cvls_mem->scalesol = SUNFALSE;

  /* Attach linear solver memory to integrator memory */
  cv_mem->cv_lmem = cvls_mem;

  return(CVLS_SUCCESS);
}


/*===============================================================
  Optional input/output routines
  ===============================================================*/


/* CVodeSetJacFn specifies the Jacobian function. */
int CVodeSetJacFn(void *cvode_mem, CVLsJacFn jac)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetJacFn",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* return with failure if jac cannot be used */
  if ((jac != NULL) && (cvls_mem->A == NULL)) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "CVodeSetJacFn",
                   "Jacobian routine cannot be supplied for NULL SUNMatrix");
    return(CVLS_ILL_INPUT);
  }

  /* set the Jacobian routine pointer, and update relevant flags */
  if (jac != NULL) {
    cvls_mem->jacDQ  = SUNFALSE;
    cvls_mem->jac    = jac;
    cvls_mem->J_data = cv_mem->cv_user_data;
  } else {
    cvls_mem->jacDQ  = SUNTRUE;
    cvls_mem->jac    = cvLsDQJac;
    cvls_mem->J_data = cv_mem;
  }

  /* ensure the internal linear system function is used */
  cvls_mem->user_linsys = SUNFALSE;
  cvls_mem->linsys      = cvLsLinSys;
  cvls_mem->A_data      = cv_mem;

  return(CVLS_SUCCESS);
}


/* CVodeSetEpsLin specifies the nonlinear -> linear tolerance scale factor */
int CVodeSetEpsLin(void *cvode_mem, realtype eplifac)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetEpsLin",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* Check for legal eplifac */
  if(eplifac < ZERO) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS",
                   "CVodeSetEpsLin", MSG_LS_BAD_EPLIN);
    return(CVLS_ILL_INPUT);
  }

  cvls_mem->eplifac = (eplifac == ZERO) ? CVLS_EPLIN : eplifac;

  return(CVLS_SUCCESS);
}


/* CVodeSetLSNormFactor sets or computes the factor to use when converting from
   the integrator tolerance to the linear solver tolerance (WRMS to L2 norm). */
int CVodeSetLSNormFactor(void *cvode_mem, realtype nrmfac)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetLSNormFactor",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS) return(retval);

  if (nrmfac > ZERO) {
    /* user-provided factor */
    cvls_mem->nrmfac = nrmfac;
  } else if (nrmfac < ZERO) {
    /* compute factor for WRMS norm with dot product */
    N_VConst(ONE, cvls_mem->ytemp);
    cvls_mem->nrmfac = SUNRsqrt(N_VDotProd(cvls_mem->ytemp, cvls_mem->ytemp));
  } else {
    /* compute default factor for WRMS norm from vector legnth */
    cvls_mem->nrmfac = SUNRsqrt(N_VGetLength(cvls_mem->ytemp));
  }

  return(CVLS_SUCCESS);
}


/* CVodeSetJacEvalFrequency specifies the frequency for recomputing the Jacobian
   matrix and/or preconditioner */
int CVodeSetJacEvalFrequency(void *cvode_mem, long int msbj)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; store input and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetJacEvalFrequency",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* Check for legal msbj */
  if(msbj < 0) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "CVodeSetJacEvalFrequency",
                   "A negative evaluation frequency was provided.");
    return(CVLS_ILL_INPUT);
  }

  cvls_mem->msbj = (msbj == 0) ? CVLS_MSBJ : msbj;

  return(CVLS_SUCCESS);
}

/* Deprecated */
int CVodeSetMaxStepsBetweenJac(void *cvode_mem, long int msbj)
{
  return(CVodeSetJacEvalFrequency(cvode_mem, msbj));
}


/* CVodeSetLinearSolutionScaling enables or disables scaling the
   linear solver solution to account for changes in gamma. */
int CVodeSetLinearSolutionScaling(void *cvode_mem, booleantype onoff)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; store input and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetLinearSolutionScaling",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS) return(retval);

  /* check for valid solver and method type */
  if (!(cvls_mem->matrixbased) || cv_mem->cv_lmm != CV_BDF)
    return(CVLS_ILL_INPUT);

  /* set solution scaling flag */
  cvls_mem->scalesol = onoff;

  return(CVLS_SUCCESS);
}


/* CVodeSetPreconditioner specifies the user-supplied preconditioner
   setup and solve routines */
int CVodeSetPreconditioner(void *cvode_mem, CVLsPrecSetupFn psetup,
                           CVLsPrecSolveFn psolve)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  PSetupFn cvls_psetup;
  PSolveFn cvls_psolve;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetPreconditioner",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* store function pointers for user-supplied routines in CVLs interface */
  cvls_mem->pset   = psetup;
  cvls_mem->psolve = psolve;

  /* issue error if LS object does not allow user-supplied preconditioning */
  if (cvls_mem->LS->ops->setpreconditioner == NULL) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS",
                   "CVodeSetPreconditioner",
                   "SUNLinearSolver object does not support user-supplied preconditioning");
    return(CVLS_ILL_INPUT);
  }

  /* notify iterative linear solver to call CVLs interface routines */
  cvls_psetup = (psetup == NULL) ? NULL : cvLsPSetup;
  cvls_psolve = (psolve == NULL) ? NULL : cvLsPSolve;
  retval = SUNLinSolSetPreconditioner(cvls_mem->LS, cv_mem,
                                      cvls_psetup, cvls_psolve);
  if (retval != SUNLS_SUCCESS) {
    cvProcessError(cv_mem, CVLS_SUNLS_FAIL, "CVLS",
                   "CVLsSetPreconditioner",
                   "Error in calling SUNLinSolSetPreconditioner");
    return(CVLS_SUNLS_FAIL);
  }

  return(CVLS_SUCCESS);
}


/* CVodeSetJacTimes specifies the user-supplied Jacobian-vector product
   setup and multiply routines */
int CVodeSetJacTimes(void *cvode_mem, CVLsJacTimesSetupFn jtsetup,
                     CVLsJacTimesVecFn jtimes)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetJacTimes",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* issue error if LS object does not allow user-supplied ATimes */
  if (cvls_mem->LS->ops->setatimes == NULL) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS",
                    "CVodeSetJacTimes",
                    "SUNLinearSolver object does not support user-supplied ATimes routine");
    return(CVLS_ILL_INPUT);
  }

  /* store function pointers for user-supplied routines in CVLs
     interface (NULL jtimes implies use of DQ default) */
  if (jtimes != NULL) {
    cvls_mem->jtimesDQ = SUNFALSE;
    cvls_mem->jtsetup  = jtsetup;
    cvls_mem->jtimes   = jtimes;
    cvls_mem->jt_data  = cv_mem->cv_user_data;
  } else {
    cvls_mem->jtimesDQ = SUNTRUE;
    cvls_mem->jtsetup  = NULL;
    cvls_mem->jtimes   = cvLsDQJtimes;
    cvls_mem->jt_f     = cv_mem->cv_f;
    cvls_mem->jt_data  = cv_mem;
  }

  return(CVLS_SUCCESS);
}


/* CVodeSetJacTimesRhsFn specifies an alternative user-supplied ODE right-hand
   side function to use in the internal finite difference Jacobian-vector
   product */
int CVodeSetJacTimesRhsFn(void *cvode_mem, CVRhsFn jtimesRhsFn)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetJacTimesRhsFn",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS) return(retval);

  /* check if using internal finite difference approximation */
  if (!(cvls_mem->jtimesDQ)) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "CVodeSetJacTimesRhsFn",
                   "Internal finite-difference Jacobian-vector product is disabled.");
    return(CVLS_ILL_INPUT);
  }

  /* store function pointers for RHS function (NULL implies use ODE RHS) */
  if (jtimesRhsFn != NULL)
    cvls_mem->jt_f = jtimesRhsFn;
  else
    cvls_mem->jt_f = cv_mem->cv_f;

  return(CVLS_SUCCESS);
}


/* CVodeSetLinSysFn specifies the linear system setup function. */
int CVodeSetLinSysFn(void *cvode_mem, CVLsLinSysFn linsys)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeSetLinSysFn",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS) return(retval);

  /* return with failure if linsys cannot be used */
  if ((linsys != NULL) && (cvls_mem->A == NULL)) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "CVodeSetLinSysFn",
                   "Linear system setup routine cannot be supplied for NULL SUNMatrix");
    return(CVLS_ILL_INPUT);
  }

  /* set the linear system routine pointer, and update relevant flags */
  if (linsys != NULL) {
    cvls_mem->user_linsys = SUNTRUE;
    cvls_mem->linsys      = linsys;
    cvls_mem->A_data      = cv_mem->cv_user_data;
  } else {
    cvls_mem->user_linsys = SUNFALSE;
    cvls_mem->linsys      = cvLsLinSys;
    cvls_mem->A_data      = cv_mem;
  }

  return(CVLS_SUCCESS);
}


/* CVodeGetLinWorkSpace returns the length of workspace allocated
   for the CVLS linear solver interface */
int CVodeGetLinWorkSpace(void *cvode_mem, long int *lenrwLS,
                         long int *leniwLS)
{
  CVodeMem     cv_mem;
  CVLsMem      cvls_mem;
  sunindextype lrw1, liw1;
  long int     lrw, liw;
  int          retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetLinWorkSpace",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* start with fixed sizes plus vector/matrix pointers */
  *lenrwLS = 2;
  *leniwLS = 30;

  /* add NVector sizes */
  if (cv_mem->cv_tempv->ops->nvspace) {
    N_VSpace(cv_mem->cv_tempv, &lrw1, &liw1);
    *lenrwLS += 2*lrw1;
    *leniwLS += 2*liw1;
  }

  /* add SUNMatrix size (only account for the one owned by Ls interface) */
  if (cvls_mem->savedJ)
    if (cvls_mem->savedJ->ops->space) {
      retval = SUNMatSpace(cvls_mem->savedJ, &lrw, &liw);
      if (retval == 0) {
        *lenrwLS += lrw;
        *leniwLS += liw;
      }
    }

  /* add LS sizes */
  if (cvls_mem->LS->ops->space) {
    retval = SUNLinSolSpace(cvls_mem->LS, &lrw, &liw);
    if (retval == 0) {
      *lenrwLS += lrw;
      *leniwLS += liw;
    }
  }

  return(CVLS_SUCCESS);
}


/* CVodeGetNumJacEvals returns the number of Jacobian evaluations */
int CVodeGetNumJacEvals(void *cvode_mem, long int *njevals)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetNumJacEvals",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *njevals = cvls_mem->nje;
  return(CVLS_SUCCESS);
}


/* CVodeGetNumLinRhsEvals returns the number of calls to the ODE
   function needed for the DQ Jacobian approximation or J*v product
   approximation */
int CVodeGetNumLinRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetNumLinRhsEvals",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *nfevalsLS = cvls_mem->nfeDQ;
  return(CVLS_SUCCESS);
}


/* CVodeGetNumPrecEvals returns the number of calls to the
   user- or CVode-supplied preconditioner setup routine */
int CVodeGetNumPrecEvals(void *cvode_mem, long int *npevals)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetNumPrecEvals",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *npevals = cvls_mem->npe;
  return(CVLS_SUCCESS);
}


/* CVodeGetNumPrecSolves returns the number of calls to the
   user- or CVode-supplied preconditioner solve routine */
int CVodeGetNumPrecSolves(void *cvode_mem, long int *npsolves)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetNumPrecSolves",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *npsolves = cvls_mem->nps;
  return(CVLS_SUCCESS);
}


/* CVodeGetNumLinIters returns the number of linear iterations
   (if accessible from the LS object) */
int CVodeGetNumLinIters(void *cvode_mem, long int *nliters)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetNumLinIters",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *nliters = cvls_mem->nli;
  return(CVLS_SUCCESS);
}


/* CVodeGetNumLinConvFails returns the number of linear solver
   convergence failures (as reported by the LS object) */
int CVodeGetNumLinConvFails(void *cvode_mem, long int *nlcfails)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetNumLinConvFails",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *nlcfails = cvls_mem->ncfl;
  return(CVLS_SUCCESS);
}


/* CVodeGetNumJTSetupEvals returns the number of calls to the
   user-supplied Jacobian-vector product setup routine */
int CVodeGetNumJTSetupEvals(void *cvode_mem, long int *njtsetups)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetNumJTSetupEvals",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *njtsetups = cvls_mem->njtsetup;
  return(CVLS_SUCCESS);
}


/* CVodeGetNumJtimesEvals returns the number of calls to the
   Jacobian-vector product multiply routine */
int CVodeGetNumJtimesEvals(void *cvode_mem, long int *njvevals)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetNumJtimesEvals",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *njvevals = cvls_mem->njtimes;
  return(CVLS_SUCCESS);
}


/* CVodeGetLinSolveStats returns statistics related to the linear solve. */
int CVodeGetLinSolveStats(void* cvode_mem, long int* njevals, long int* nfevalsLS,
                          long int* nliters, long int* nlcfails, long int* npevals,
                          long int* npsolves, long int* njtsetups, long int* njtimes)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetLinSolveStats",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  *njevals   = cvls_mem->nje;
  *nfevalsLS = cvls_mem->nfeDQ;
  *nliters   = cvls_mem->nli;
  *nlcfails  = cvls_mem->ncfl;
  *npevals   = cvls_mem->npe;
  *npsolves  = cvls_mem->nps;
  *njtsetups = cvls_mem->njtsetup;
  *njtimes   = cvls_mem->njtimes;

  return(CVLS_SUCCESS);
}


/* CVodeGetLastLinFlag returns the last flag set in a CVLS function */
int CVodeGetLastLinFlag(void *cvode_mem, long int *flag)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure; set output value and return */
  retval = cvLs_AccessLMem(cvode_mem, "CVodeGetLastLinFlag",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);
  *flag = cvls_mem->last_flag;
  return(CVLS_SUCCESS);
}


/* CVodeGetLinReturnFlagName translates from the integer error code
   returned by an CVLs routine to the corresponding string
   equivalent for that flag */
char *CVodeGetLinReturnFlagName(long int flag)
{
  char *name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVLS_SUCCESS:
    sprintf(name,"CVLS_SUCCESS");
    break;
  case CVLS_MEM_NULL:
    sprintf(name,"CVLS_MEM_NULL");
    break;
  case CVLS_LMEM_NULL:
    sprintf(name,"CVLS_LMEM_NULL");
    break;
  case CVLS_ILL_INPUT:
    sprintf(name,"CVLS_ILL_INPUT");
    break;
  case CVLS_MEM_FAIL:
    sprintf(name,"CVLS_MEM_FAIL");
    break;
  case CVLS_PMEM_NULL:
    sprintf(name,"CVLS_PMEM_NULL");
    break;
  case CVLS_JACFUNC_UNRECVR:
    sprintf(name,"CVLS_JACFUNC_UNRECVR");
    break;
  case CVLS_JACFUNC_RECVR:
    sprintf(name,"CVLS_JACFUNC_RECVR");
    break;
  case CVLS_SUNMAT_FAIL:
    sprintf(name,"CVLS_SUNMAT_FAIL");
    break;
  case CVLS_SUNLS_FAIL:
    sprintf(name,"CVLS_SUNLS_FAIL");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


/*=================================================================
  CVLS private functions
  =================================================================*/

/*-----------------------------------------------------------------
  cvLsATimes

  This routine generates the matrix-vector product z = Mv, where
  M = I - gamma*J. The product J*v is obtained by calling the jtimes
  routine. It is then scaled by -gamma and added to v to obtain M*v.
  The return value is the same as the value returned by jtimes --
  0 if successful, nonzero otherwise.
  -----------------------------------------------------------------*/
int cvLsATimes(void *cvode_mem, N_Vector v, N_Vector z)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "cvLsATimes",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* call Jacobian-times-vector product routine
     (either user-supplied or internal DQ) */
  retval = cvls_mem->jtimes(v, z, cv_mem->cv_tn,
                            cvls_mem->ycur,
                            cvls_mem->fcur,
                            cvls_mem->jt_data,
                            cvls_mem->ytemp);
  cvls_mem->njtimes++;
  if (retval != 0) return(retval);

  /* add contribution from identity matrix */
  N_VLinearSum(ONE, v, -cv_mem->cv_gamma, z, z);

  return(0);
}


/*---------------------------------------------------------------
  cvLsPSetup:

  This routine interfaces between the generic iterative linear
  solvers and the user's psetup routine.  It passes to psetup all
  required state information from cvode_mem.  Its return value
  is the same as that returned by psetup. Note that the generic
  iterative linear solvers guarantee that cvLsPSetup will only
  be called in the case that the user's psetup routine is non-NULL.
  ---------------------------------------------------------------*/
int cvLsPSetup(void *cvode_mem)
{
  int      retval;
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "cvLsPSetup",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* Call user pset routine to update preconditioner and possibly
     reset jcur (pass !jbad as update suggestion) */
  retval = cvls_mem->pset(cv_mem->cv_tn, cvls_mem->ycur,
                          cvls_mem->fcur, !(cvls_mem->jbad),
                          &cv_mem->cv_jcur, cv_mem->cv_gamma,
                          cvls_mem->P_data);
  return(retval);
}


/*-----------------------------------------------------------------
  cvLsPSolve

  This routine interfaces between the generic SUNLinSolSolve
  routine and the user's psolve routine.  It passes to psolve all
  required state information from cvode_mem.  Its return value is
  the same as that returned by psolve. Note that the generic
  SUNLinSol solver guarantees that cvLsPSolve will not be called
  in the case in which preconditioning is not done. This is the
  only case in which the user's psolve routine is allowed to be
  NULL.
  -----------------------------------------------------------------*/
int cvLsPSolve(void *cvode_mem, N_Vector r, N_Vector z, realtype tol, int lr)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "cvLsPSolve",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* call the user-supplied psolve routine, and accumulate count */
  retval = cvls_mem->psolve(cv_mem->cv_tn, cvls_mem->ycur,
                            cvls_mem->fcur, r, z,
                            cv_mem->cv_gamma, tol, lr,
                            cvls_mem->P_data);
  cvls_mem->nps++;
  return(retval);
}


/*-----------------------------------------------------------------
  cvLsDQJac

  This routine is a wrapper for the Dense and Band
  implementations of the difference quotient Jacobian
  approximation routines.
  ---------------------------------------------------------------*/
int cvLsDQJac(realtype t, N_Vector y, N_Vector fy,
              SUNMatrix Jac, void *cvode_mem, N_Vector tmp1,
              N_Vector tmp2, N_Vector tmp3)
{
  CVodeMem cv_mem;
  int      retval;

  /* access CVodeMem structure */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVLS_MEM_NULL, "CVLS",
                   "cvLsDQJac", MSG_LS_CVMEM_NULL);
    return(CVLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* verify that Jac is non-NULL */
  if (Jac == NULL) {
    cvProcessError(cv_mem, CVLS_LMEM_NULL, "CVLS",
                   "cvLsDQJac", MSG_LS_LMEM_NULL);
    return(CVLS_LMEM_NULL);
  }

  /* Verify that N_Vector supports required operations */
  if (cv_mem->cv_tempv->ops->nvcloneempty == NULL ||
      cv_mem->cv_tempv->ops->nvwrmsnorm == NULL ||
      cv_mem->cv_tempv->ops->nvlinearsum == NULL ||
      cv_mem->cv_tempv->ops->nvdestroy == NULL ||
      cv_mem->cv_tempv->ops->nvscale == NULL ||
      cv_mem->cv_tempv->ops->nvgetarraypointer == NULL ||
      cv_mem->cv_tempv->ops->nvsetarraypointer == NULL) {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS",
                   "cvLsDQJac", MSG_LS_BAD_NVECTOR);
    return(CVLS_ILL_INPUT);
  }

  /* Call the matrix-structure-specific DQ approximation routine */
  if (SUNMatGetID(Jac) == SUNMATRIX_DENSE) {
    retval = cvLsDenseDQJac(t, y, fy, Jac, cv_mem, tmp1);
  } else if (SUNMatGetID(Jac) == SUNMATRIX_BAND) {
    retval = cvLsBandDQJac(t, y, fy, Jac, cv_mem, tmp1, tmp2);
  } else {
    cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "cvLsDQJac",
                   "unrecognized matrix type for cvLsDQJac");
    retval = CVLS_ILL_INPUT;
  }
  return(retval);
}


/*-----------------------------------------------------------------
  cvLsDenseDQJac

  This routine generates a dense difference quotient approximation
  to the Jacobian of f(t,y). It assumes that a dense SUNMatrix is
  stored column-wise, and that elements within each column are
  contiguous. The address of the jth column of J is obtained via
  the accessor function SUNDenseMatrix_Column, and this pointer
  is associated with an N_Vector using the N_VSetArrayPointer
  function.  Finally, the actual computation of the jth column of
  the Jacobian is done with a call to N_VLinearSum.
  -----------------------------------------------------------------*/
int cvLsDenseDQJac(realtype t, N_Vector y, N_Vector fy,
                   SUNMatrix Jac, CVodeMem cv_mem, N_Vector tmp1)
{
  realtype fnorm, minInc, inc, inc_inv, yjsaved, srur, conj;
  realtype *y_data, *ewt_data, *cns_data;
  N_Vector ftemp, jthCol;
  sunindextype j, N;
  CVLsMem cvls_mem;
  int retval = 0;

  /* initialize cns_data to avoid compiler warning */
  cns_data = NULL;

  /* access LsMem interface structure */
  cvls_mem = (CVLsMem) cv_mem->cv_lmem;

  /* access matrix dimension */
  N = SUNDenseMatrix_Columns(Jac);

  /* Rename work vector for readibility */
  ftemp = tmp1;

  /* Create an empty vector for matrix column calculations */
  jthCol = N_VCloneEmpty(tmp1);

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(cv_mem->cv_ewt);
  y_data   = N_VGetArrayPointer(y);
  if (cv_mem->cv_constraintsSet)
    cns_data = N_VGetArrayPointer(cv_mem->cv_constraints);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(cv_mem->cv_uround);
  fnorm = N_VWrmsNorm(fy, cv_mem->cv_ewt);
  minInc = (fnorm != ZERO) ?
    (MIN_INC_MULT * SUNRabs(cv_mem->cv_h) * cv_mem->cv_uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */
    N_VSetArrayPointer(SUNDenseMatrix_Column(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = SUNMAX(srur*SUNRabs(yjsaved), minInc/ewt_data[j]);

    /* Adjust sign(inc) if y_j has an inequality constraint. */
    if (cv_mem->cv_constraintsSet) {
      conj = cns_data[j];
      if (SUNRabs(conj) == ONE)      {if ((yjsaved+inc)*conj < ZERO)  inc = -inc;}
      else if (SUNRabs(conj) == TWO) {if ((yjsaved+inc)*conj <= ZERO) inc = -inc;}
    }

    y_data[j] += inc;

    retval = cv_mem->cv_f(t, y, ftemp, cv_mem->cv_user_data);
    cvls_mem->nfeDQ++;
    if (retval != 0) break;

    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

  }

  /* Destroy jthCol vector */
  N_VSetArrayPointer(NULL, jthCol);  /* SHOULDN'T BE NEEDED */
  N_VDestroy(jthCol);

  return(retval);
}


/*-----------------------------------------------------------------
  cvLsBandDQJac

  This routine generates a banded difference quotient approximation
  to the Jacobian of f(t,y).  It assumes that a band SUNMatrix is
  stored column-wise, and that elements within each column are
  contiguous. This makes it possible to get the address of a column
  of J via the accessor function SUNBandMatrix_Column, and to write
  a simple for loop to set each of the elements of a column in
  succession.
  -----------------------------------------------------------------*/
int cvLsBandDQJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix Jac,
                  CVodeMem cv_mem, N_Vector tmp1, N_Vector tmp2)
{
  N_Vector ftemp, ytemp;
  realtype fnorm, minInc, inc, inc_inv, srur, conj;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data;
  realtype *y_data, *ytemp_data, *cns_data;
  sunindextype group, i, j, width, ngroups, i1, i2;
  sunindextype N, mupper, mlower;
  CVLsMem cvls_mem;
  int retval = 0;

  /* initialize cns_data to avoid compiler warning */
  cns_data = NULL;

  /* access LsMem interface structure */
  cvls_mem = (CVLsMem) cv_mem->cv_lmem;

  /* access matrix dimensions */
  N = SUNBandMatrix_Columns(Jac);
  mupper = SUNBandMatrix_UpperBandwidth(Jac);
  mlower = SUNBandMatrix_LowerBandwidth(Jac);

  /* Rename work vectors for use as temporary values of y and f */
  ftemp = tmp1;
  ytemp = tmp2;

  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp */
  ewt_data   = N_VGetArrayPointer(cv_mem->cv_ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);
  if (cv_mem->cv_constraintsSet)
    cns_data = N_VGetArrayPointer(cv_mem->cv_constraints);

  /* Load ytemp with y = predicted y vector */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f */
  srur = SUNRsqrt(cv_mem->cv_uround);
  fnorm = N_VWrmsNorm(fy, cv_mem->cv_ewt);
  minInc = (fnorm != ZERO) ?
    (MIN_INC_MULT * SUNRabs(cv_mem->cv_h) * cv_mem->cv_uround * N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mlower + mupper + 1;
  ngroups = SUNMIN(width, N);

  /* Loop over column groups. */
  for (group=1; group <= ngroups; group++) {

    /* Increment all y_j in group */
    for(j=group-1; j < N; j+=width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);

      /* Adjust sign(inc) if yj has an inequality constraint. */
      if (cv_mem->cv_constraintsSet) {
        conj = cns_data[j];
        if (SUNRabs(conj) == ONE)      {if ((ytemp_data[j]+inc)*conj < ZERO)  inc = -inc;}
        else if (SUNRabs(conj) == TWO) {if ((ytemp_data[j]+inc)*conj <= ZERO) inc = -inc;}
      }

      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y */
    retval = cv_mem->cv_f(cv_mem->cv_tn, ytemp, ftemp, cv_mem->cv_user_data);
    cvls_mem->nfeDQ++;
    if (retval != 0) break;

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < N; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = SUNBandMatrix_Column(Jac, j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);

      /* Adjust sign(inc) as before. */
      if (cv_mem->cv_constraintsSet) {
        conj = cns_data[j];
        if (SUNRabs(conj) == ONE)      {if ((ytemp_data[j]+inc)*conj < ZERO)  inc = -inc;}
        else if (SUNRabs(conj) == TWO) {if ((ytemp_data[j]+inc)*conj <= ZERO) inc = -inc;}
      }

      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mupper);
      i2 = SUNMIN(j+mlower, N-1);
      for (i=i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j,i,j) = inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(retval);
}


/*-----------------------------------------------------------------
  cvLsDQJtimes

  This routine generates a difference quotient approximation to
  the Jacobian times vector f_y(t,y) * v. The approximation is
  Jv = [f(y + v*sig) - f(y)]/sig, where sig = 1 / ||v||_WRMS,
  i.e. the WRMS norm of v*sig is 1.
  -----------------------------------------------------------------*/
int cvLsDQJtimes(N_Vector v, N_Vector Jv, realtype t,
                 N_Vector y, N_Vector fy, void *cvode_mem,
                 N_Vector work)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  realtype sig, siginv;
  int      iter, retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "cvLsDQJtimes",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* Initialize perturbation to 1/||v|| */
  sig = ONE/N_VWrmsNorm(v, cv_mem->cv_ewt);

  for (iter=0; iter<MAX_DQITERS; iter++) {

    /* Set work = y + sig*v */
    N_VLinearSum(sig, v, ONE, y, work);

    /* Set Jv = f(tn, y+sig*v) */
    retval = cvls_mem->jt_f(t, work, Jv, cv_mem->cv_user_data);
    cvls_mem->nfeDQ++;
    if (retval == 0) break;
    if (retval < 0)  return(-1);

    /* If f failed recoverably, shrink sig and retry */
    sig *= PT25;
  }

  /* If retval still isn't 0, return with a recoverable failure */
  if (retval > 0) return(+1);

  /* Replace Jv by (Jv - fy)/sig */
  siginv = ONE/sig;
  N_VLinearSum(siginv, Jv, -siginv, fy, Jv);

  return(0);
}


/*-----------------------------------------------------------------
  cvLsLinSys

  Setup the linear system A = I - gamma J
  -----------------------------------------------------------------*/
static int cvLsLinSys(realtype t, N_Vector y, N_Vector fy, SUNMatrix A,
                      booleantype jok, booleantype *jcur, realtype gamma,
                      void *cvode_mem, N_Vector vtemp1, N_Vector vtemp2,
                      N_Vector vtemp3)
{
  CVodeMem cv_mem;
  CVLsMem  cvls_mem;
  int      retval;

  /* access CVLsMem structure */
  retval = cvLs_AccessLMem(cvode_mem, "cvLsLinSys",
                           &cv_mem, &cvls_mem);
  if (retval != CVLS_SUCCESS)  return(retval);

  /* Check if Jacobian needs to be updated */
  if (jok) {

    /* Use saved copy of J */
    *jcur = SUNFALSE;

    /* Overwrite linear system matrix with saved J */
    retval = SUNMatCopy(cvls_mem->savedJ, A);
    if (retval) {
      cvProcessError(cv_mem, CVLS_SUNMAT_FAIL, "CVLS",
                     "cvLsSetup",  MSG_LS_SUNMAT_FAILED);
      cvls_mem->last_flag = CVLS_SUNMAT_FAIL;
      return(cvls_mem->last_flag);
    }

  } else {

    /* Call jac() routine to update J */
    *jcur = SUNTRUE;

    /* Clear the linear system matrix if necessary */
    if (SUNLinSolGetType(cvls_mem->LS) == SUNLINEARSOLVER_DIRECT) {
      retval = SUNMatZero(A);
      if (retval) {
        cvProcessError(cv_mem, CVLS_SUNMAT_FAIL, "CVLS",
                       "cvLsSetup",  MSG_LS_SUNMAT_FAILED);
        cvls_mem->last_flag = CVLS_SUNMAT_FAIL;
        return(cvls_mem->last_flag);
      }
    }

    /* Compute new Jacobian matrix */
    retval = cvls_mem->jac(t, y, fy, A, cvls_mem->J_data,
                           vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      cvProcessError(cv_mem, CVLS_JACFUNC_UNRECVR, "CVLS",
                     "cvLsSetup",  MSG_LS_JACFUNC_FAILED);
      cvls_mem->last_flag = CVLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      cvls_mem->last_flag = CVLS_JACFUNC_RECVR;
      return(1);
    }

    /* Update saved copy of the Jacobian matrix */
    retval = SUNMatCopy(A, cvls_mem->savedJ);
    if (retval) {
      cvProcessError(cv_mem, CVLS_SUNMAT_FAIL, "CVLS",
                     "cvLsSetup",  MSG_LS_SUNMAT_FAILED);
      cvls_mem->last_flag = CVLS_SUNMAT_FAIL;
      return(cvls_mem->last_flag);
    }

  }

  /* Perform linear combination A = I - gamma*J */
  retval = SUNMatScaleAddI(-gamma, A);
  if (retval) {
    cvProcessError(cv_mem, CVLS_SUNMAT_FAIL, "CVLS",
                   "cvLsSetup",  MSG_LS_SUNMAT_FAILED);
    cvls_mem->last_flag = CVLS_SUNMAT_FAIL;
    return(cvls_mem->last_flag);
  }

  return(CVLS_SUCCESS);
}


/*-----------------------------------------------------------------
  cvLsInitialize

  This routine performs remaining initializations specific
  to the iterative linear solver interface (and solver itself)
  -----------------------------------------------------------------*/
int cvLsInitialize(CVodeMem cv_mem)
{
  CVLsMem cvls_mem;
  int     retval;

  /* access CVLsMem structure */
  if (cv_mem->cv_lmem==NULL) {
    cvProcessError(cv_mem, CVLS_LMEM_NULL, "CVLS",
                   "cvLsInitialize", MSG_LS_LMEM_NULL);
    return(CVLS_LMEM_NULL);
  }
  cvls_mem = (CVLsMem) cv_mem->cv_lmem;

  /* Test for valid combinations of matrix & Jacobian routines: */
  if (cvls_mem->A != NULL) {

    /* Matrix-based case */

    if (cvls_mem->user_linsys) {

      /* User-supplied linear system function, reset A_data (just in case) */
      cvls_mem->A_data = cv_mem->cv_user_data;

    } else {

      /* Internal linear system function, reset pointers (just in case) */
      cvls_mem->linsys = cvLsLinSys;
      cvls_mem->A_data = cv_mem;

      /* Check if an internal or user-supplied Jacobian function is used */
      if (cvls_mem->jacDQ) {

        /* Internal difference quotient Jacobian. Check that A is dense or band,
           otherwise return an error */
        retval = 0;
        if (cvls_mem->A->ops->getid) {

          if ( (SUNMatGetID(cvls_mem->A) == SUNMATRIX_DENSE) ||
               (SUNMatGetID(cvls_mem->A) == SUNMATRIX_BAND) ) {
            cvls_mem->jac    = cvLsDQJac;
            cvls_mem->J_data = cv_mem;
          } else {
            retval++;
          }

        } else {
          retval++;
        }
        if (retval) {
          cvProcessError(cv_mem, CVLS_ILL_INPUT, "CVLS", "cvLsInitialize",
                         "No Jacobian constructor available for SUNMatrix type");
          cvls_mem->last_flag = CVLS_ILL_INPUT;
          return(CVLS_ILL_INPUT);
        }

      } else {

        /* User-supplied Jacobian, reset J_data pointer (just in case) */
        cvls_mem->J_data = cv_mem->cv_user_data;

      }

      /* Allocate internally saved Jacobian if not already done */
      if (cvls_mem->savedJ == NULL) {
        cvls_mem->savedJ = SUNMatClone(cvls_mem->A);
        if (cvls_mem->savedJ == NULL) {
          cvProcessError(cv_mem, CVLS_MEM_FAIL, "CVLS", "cvLsInitialize",
                         MSG_LS_MEM_FAIL);
          cvls_mem->last_flag = CVLS_MEM_FAIL;
          return(CVLS_MEM_FAIL);
        }
      }

    } /* end matrix-based case */

  } else {

    /* Matrix-free case: ensure 'jac' and `linsys` function pointers are NULL */
    cvls_mem->jacDQ  = SUNFALSE;
    cvls_mem->jac    = NULL;
    cvls_mem->J_data = NULL;

    cvls_mem->user_linsys = SUNFALSE;
    cvls_mem->linsys      = NULL;
    cvls_mem->A_data      = NULL;

  }

  /* reset counters */
  cvLsInitializeCounters(cvls_mem);

  /* Set Jacobian-vector product related fields, based on jtimesDQ */
  if (cvls_mem->jtimesDQ) {
    cvls_mem->jtsetup = NULL;
    cvls_mem->jtimes  = cvLsDQJtimes;
    cvls_mem->jt_data = cv_mem;
  } else {
    cvls_mem->jt_data = cv_mem->cv_user_data;
  }

  /* if A is NULL and psetup is not present, then cvLsSetup does
     not need to be called, so set the lsetup function to NULL */
  if ( (cvls_mem->A == NULL) && (cvls_mem->pset == NULL) )
    cv_mem->cv_lsetup = NULL;

  /* Call LS initialize routine, and return result */
  cvls_mem->last_flag = SUNLinSolInitialize(cvls_mem->LS);
  return(cvls_mem->last_flag);
}


/*-----------------------------------------------------------------
  cvLsSetup

  This conditionally calls the LS 'setup' routine.

  When using a SUNMatrix object, this determines whether
  to update a Jacobian matrix (or use a stored version), based
  on heuristics regarding previous convergence issues, the number
  of time steps since it was last updated, etc.; it then creates
  the system matrix from this, the 'gamma' factor and the
  identity matrix, A = I-gamma*J.

  This routine then calls the LS 'setup' routine with A.
  -----------------------------------------------------------------*/
int cvLsSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
              N_Vector fpred, booleantype *jcurPtr,
              N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  CVLsMem  cvls_mem;
  realtype dgamma;
  int      retval;

  /* access CVLsMem structure */
  if (cv_mem->cv_lmem==NULL) {
    cvProcessError(cv_mem, CVLS_LMEM_NULL, "CVLS",
                   "cvLsSetup", MSG_LS_LMEM_NULL);
    return(CVLS_LMEM_NULL);
  }
  cvls_mem = (CVLsMem) cv_mem->cv_lmem;

  /* Set CVLs N_Vector pointers to current solution and rhs */
  cvls_mem->ycur = ypred;
  cvls_mem->fcur = fpred;

  /* Use nst, gamma/gammap, and convfail to set J/P eval. flag jok */
  dgamma = SUNRabs((cv_mem->cv_gamma/cv_mem->cv_gammap) - ONE);
  cvls_mem->jbad = (cv_mem->cv_nst == 0) ||
    (cv_mem->cv_nst >= cvls_mem->nstlj + cvls_mem->msbj) ||
    ((convfail == CV_FAIL_BAD_J) && (dgamma < CVLS_DGMAX)) ||
    (convfail == CV_FAIL_OTHER);

  /* Setup the linear system if necessary */
  if (cvls_mem->A != NULL) {

    /* Update J if appropriate and evaluate A = I - gamma J */
    retval = cvls_mem->linsys(cv_mem->cv_tn, ypred, fpred, cvls_mem->A,
                              !(cvls_mem->jbad), jcurPtr, cv_mem->cv_gamma,
                              cvls_mem->A_data, vtemp1, vtemp2, vtemp3);

    /* Update J eval count and step when J was last updated */
    if (*jcurPtr) {
      cvls_mem->nje++;
      cvls_mem->nstlj = cv_mem->cv_nst;
    }

    /* Check linsys() return value and return if necessary */
    if (retval != CVLS_SUCCESS) {
      if (cvls_mem->user_linsys) {
        if (retval < 0) {
          cvProcessError(cv_mem, CVLS_JACFUNC_UNRECVR, "CVLS",
                         "cvLsSetup",  MSG_LS_JACFUNC_FAILED);
          cvls_mem->last_flag = CVLS_JACFUNC_UNRECVR;
          return(-1);
        } else {
          cvls_mem->last_flag = CVLS_JACFUNC_RECVR;
          return(1);
        }
      } else {
        return(retval);
      }
    }

  } else {

    /* Matrix-free case, set jcur to jbad */
    *jcurPtr = cvls_mem->jbad;

  }

  /* Call LS setup routine -- the LS may call cvLsPSetup, who will
     pass the heuristic suggestions above to the user code(s) */
  cvls_mem->last_flag = SUNLinSolSetup(cvls_mem->LS, cvls_mem->A);

  /* If Matrix-free, update heuristics flags */
  if (cvls_mem->A == NULL) {

    /* If user set jcur to SUNTRUE, increment npe and save nst value */
    if (*jcurPtr) {
      cvls_mem->npe++;
      cvls_mem->nstlj = cv_mem->cv_nst;
    }

    /* Update jcur flag if we suggested an update */
    if (cvls_mem->jbad) *jcurPtr = SUNTRUE;
  }

  return(cvls_mem->last_flag);
}


/*-----------------------------------------------------------------
  cvLsSolve

  This routine interfaces between CVode and the generic
  SUNLinearSolver object LS, by setting the appropriate tolerance
  and scaling vectors, calling the solver, and accumulating
  statistics from the solve for use/reporting by CVode.
  -----------------------------------------------------------------*/
int cvLsSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
              N_Vector ynow, N_Vector fnow)
{
  CVLsMem  cvls_mem;
  realtype bnorm, deltar, delta, w_mean;
  int      curiter, nli_inc, retval;

  /* access CVLsMem structure */
  if (cv_mem->cv_lmem==NULL) {
    cvProcessError(cv_mem, CVLS_LMEM_NULL, "CVLS",
                   "cvLsSolve", MSG_LS_LMEM_NULL);
    return(CVLS_LMEM_NULL);
  }
  cvls_mem = (CVLsMem) cv_mem->cv_lmem;

  /* get current nonlinear solver iteration */
  retval = SUNNonlinSolGetCurIter(cv_mem->NLS, &curiter);

  /* If the linear solver is iterative:
     test norm(b), if small, return x = 0 or x = b;
     set linear solver tolerance (in left/right scaled 2-norm) */
  if (cvls_mem->iterative) {
    deltar = cvls_mem->eplifac * cv_mem->cv_tq[4];
    bnorm = N_VWrmsNorm(b, weight);
    if (bnorm <= deltar) {
      if (curiter > 0) N_VConst(ZERO, b);
      cvls_mem->last_flag = CVLS_SUCCESS;
      return(cvls_mem->last_flag);
    }
    /* Adjust tolerance for 2-norm */
    delta = deltar * cvls_mem->nrmfac;
  } else {
    delta = ZERO;
  }

  /* Set vectors ycur and fcur for use by the Atimes and Psolve
     interface routines */
  cvls_mem->ycur = ynow;
  cvls_mem->fcur = fnow;

  /* Set scaling vectors for LS to use (if applicable) */
  if (cvls_mem->LS->ops->setscalingvectors) {
    retval = SUNLinSolSetScalingVectors(cvls_mem->LS,
                                        weight,
                                        weight);
    if (retval != SUNLS_SUCCESS) {
      cvProcessError(cv_mem, CVLS_SUNLS_FAIL, "CVLS", "cvLsSolve",
                     "Error in calling SUNLinSolSetScalingVectors");
      cvls_mem->last_flag = CVLS_SUNLS_FAIL;
      return(cvls_mem->last_flag);
    }

  /* If solver is iterative and does not support scaling vectors, update the
     tolerance in an attempt to account for weight vector.  We make the
     following assumptions:
       1. w_i = w_mean, for i=0,...,n-1 (i.e. the weights are homogeneous)
       2. the linear solver uses a basic 2-norm to measure convergence
     Hence (using the notation from sunlinsol_spgmr.h, with S = diag(w)),
           || bbar - Abar xbar ||_2 < tol
       <=> || S b - S A x ||_2 < tol
       <=> || S (b - A x) ||_2 < tol
       <=> \sum_{i=0}^{n-1} (w_i (b - A x)_i)^2 < tol^2
       <=> w_mean^2 \sum_{i=0}^{n-1} (b - A x_i)^2 < tol^2
       <=> \sum_{i=0}^{n-1} (b - A x_i)^2 < tol^2 / w_mean^2
       <=> || b - A x ||_2 < tol / w_mean
     So we compute w_mean = ||w||_RMS = ||w||_2 and scale the desired tolerance accordingly. */
  } else if (cvls_mem->iterative) {

    N_VConst(ONE, cvls_mem->x);
    w_mean = N_VWrmsNorm(weight, cvls_mem->x);
    delta /= w_mean;

  }

  /* Set initial guess x = 0 to LS */
  N_VConst(ZERO, cvls_mem->x);

  /* If a user-provided jtsetup routine is supplied, call that here */
  if (cvls_mem->jtsetup) {
    cvls_mem->last_flag = cvls_mem->jtsetup(cv_mem->cv_tn, ynow, fnow,
                                            cvls_mem->jt_data);
    cvls_mem->njtsetup++;
    if (cvls_mem->last_flag != 0) {
      cvProcessError(cv_mem, retval, "CVLS",
                     "cvLsSolve", MSG_LS_JTSETUP_FAILED);
      return(cvls_mem->last_flag);
    }
  }

  /* Call solver, and copy x to b */
  retval = SUNLinSolSolve(cvls_mem->LS, cvls_mem->A, cvls_mem->x, b, delta);
  N_VScale(ONE, cvls_mem->x, b);

  /* If using a direct or matrix-iterative solver, BDF method, and gamma has changed,
     scale the correction to account for change in gamma */
  if (cvls_mem->scalesol && cv_mem->cv_gamrat != ONE)
    N_VScale(TWO/(ONE + cv_mem->cv_gamrat), b, b);

  /* Retrieve statistics from iterative linear solvers */
  nli_inc = 0;
  if (cvls_mem->iterative && cvls_mem->LS->ops->numiters)
    nli_inc = SUNLinSolNumIters(cvls_mem->LS);

  /* Increment counters nli and ncfl */
  cvls_mem->nli += nli_inc;
  if (retval != SUNLS_SUCCESS) cvls_mem->ncfl++;

  /* Interpret solver return value  */
  cvls_mem->last_flag = retval;

  switch(retval) {

  case SUNLS_SUCCESS:
    return(0);
    break;
  case SUNLS_RES_REDUCED:
    /* allow reduction but not solution on first Newton iteration,
       otherwise return with a recoverable failure */
    if (curiter == 0) return(0);
    else              return(1);
    break;
  case SUNLS_CONV_FAIL:
  case SUNLS_ATIMES_FAIL_REC:
  case SUNLS_PSOLVE_FAIL_REC:
  case SUNLS_PACKAGE_FAIL_REC:
  case SUNLS_QRFACT_FAIL:
  case SUNLS_LUFACT_FAIL:
    return(1);
    break;
  case SUNLS_MEM_NULL:
  case SUNLS_ILL_INPUT:
  case SUNLS_MEM_FAIL:
  case SUNLS_GS_FAIL:
  case SUNLS_QRSOL_FAIL:
    return(-1);
    break;
  case SUNLS_PACKAGE_FAIL_UNREC:
    cvProcessError(cv_mem, SUNLS_PACKAGE_FAIL_UNREC, "CVLS",
                   "cvLsSolve",
                    "Failure in SUNLinSol external package");
    return(-1);
    break;
  case SUNLS_ATIMES_FAIL_UNREC:
    cvProcessError(cv_mem, SUNLS_ATIMES_FAIL_UNREC, "CVLS",
                   "cvLsSolve", MSG_LS_JTIMES_FAILED);
    return(-1);
    break;
  case SUNLS_PSOLVE_FAIL_UNREC:
    cvProcessError(cv_mem, SUNLS_PSOLVE_FAIL_UNREC, "CVLS",
                   "cvLsSolve", MSG_LS_PSOLVE_FAILED);
    return(-1);
    break;
  }

  return(0);
}


/*-----------------------------------------------------------------
  cvLsFree

  This routine frees memory associates with the CVLs system
  solver interface.
  -----------------------------------------------------------------*/
int cvLsFree(CVodeMem cv_mem)
{
  CVLsMem cvls_mem;

  /* Return immediately if CVodeMem or CVLsMem  are NULL */
  if (cv_mem == NULL)  return (CVLS_SUCCESS);
  if (cv_mem->cv_lmem == NULL)  return(CVLS_SUCCESS);
  cvls_mem = (CVLsMem) cv_mem->cv_lmem;

  /* Free N_Vector memory */
  if (cvls_mem->ytemp) {
    N_VDestroy(cvls_mem->ytemp);
    cvls_mem->ytemp = NULL;
  }
  if (cvls_mem->x) {
    N_VDestroy(cvls_mem->x);
    cvls_mem->x = NULL;
  }

  /* Free savedJ memory */
  if (cvls_mem->savedJ) {
    SUNMatDestroy(cvls_mem->savedJ);
    cvls_mem->savedJ = NULL;
  }

  /* Nullify other N_Vector pointers */
  cvls_mem->ycur = NULL;
  cvls_mem->fcur = NULL;

  /* Nullify other SUNMatrix pointer */
  cvls_mem->A = NULL;

  /* Free preconditioner memory (if applicable) */
  if (cvls_mem->pfree)  cvls_mem->pfree(cv_mem);

  /* free CVLs interface structure */
  free(cv_mem->cv_lmem);

  return(CVLS_SUCCESS);
}


/*-----------------------------------------------------------------
  cvLsInitializeCounters

  This routine resets all counters from an CVLsMem structure.
  -----------------------------------------------------------------*/
int cvLsInitializeCounters(CVLsMem cvls_mem)
{
  cvls_mem->nje      = 0;
  cvls_mem->nfeDQ    = 0;
  cvls_mem->nstlj    = 0;
  cvls_mem->npe      = 0;
  cvls_mem->nli      = 0;
  cvls_mem->nps      = 0;
  cvls_mem->ncfl     = 0;
  cvls_mem->njtsetup = 0;
  cvls_mem->njtimes  = 0;
  return(0);
}


/*---------------------------------------------------------------
  cvLs_AccessLMem

  This routine unpacks the cv_mem and ls_mem structures from
  void* pointer.  If either is missing it returns CVLS_MEM_NULL
  or CVLS_LMEM_NULL.
  ---------------------------------------------------------------*/
int cvLs_AccessLMem(void* cvode_mem, const char *fname,
                    CVodeMem *cv_mem, CVLsMem *cvls_mem)
{
  if (cvode_mem==NULL) {
    cvProcessError(NULL, CVLS_MEM_NULL, "CVLS",
                   fname, MSG_LS_CVMEM_NULL);
    return(CVLS_MEM_NULL);
  }
  *cv_mem = (CVodeMem) cvode_mem;
  if ((*cv_mem)->cv_lmem==NULL) {
    cvProcessError(*cv_mem, CVLS_LMEM_NULL, "CVLS",
                   fname, MSG_LS_LMEM_NULL);
    return(CVLS_LMEM_NULL);
  }
  *cvls_mem = (CVLsMem) (*cv_mem)->cv_lmem;
  return(CVLS_SUCCESS);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
