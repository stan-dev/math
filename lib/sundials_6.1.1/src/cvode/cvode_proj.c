/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * Based on CPODES by Radu Serban @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * Implementation file for projections in CVODE.
 * ---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>

#include "sundials/sundials_math.h"
#include "cvode_impl.h"

/* Private constants */
#define ZERO  RCONST(0.0)  /* real 0.0 */
#define ONE   RCONST(1.0)  /* real 1.0 */

#define ONEPSM RCONST(1.000001)

/* Private utility function prototypes */
static int cvProjCreate(CVodeProjMem *proj_mem);
static int cvProjSetDefaults(CVodeProjMem proj_mem);
static int cvAccessProjMem(void* cvode_mem, const char *fname,
                           CVodeMem *cv_mem, CVodeProjMem *proj_mem);


/* ===========================================================================
 * Exported Functions - projection initialization
 * ===========================================================================*/

/* -----------------------------------------------------------------------------
 * CVodeSetProjFn sets a user defined projection function
 * ---------------------------------------------------------------------------*/
int CVodeSetProjFn(void *cvode_mem, CVProjFn pfun)
{
  int          retval;
  CVodeMem     cv_mem;
  CVodeProjMem proj_mem;

  /* Check the CVODE memory pointer */
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeSetProjFn",
                   MSG_CV_MEM_NULL);
    return(CV_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Check if the projection function is NULL */
  if (pfun == NULL)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeSetProjFn",
                   "The projection function is NULL.");
    return(CV_ILL_INPUT);
  }

  /* Check for compatible method */
  if (cv_mem->cv_lmm != CV_BDF)
  {
    cvProcessError(cv_mem, CV_ILL_INPUT, "CVODE", "CVodeSetProjFn",
                   "Projection is only supported with BDF methods.");
    return(CV_ILL_INPUT);
  }

  /* Create the projection memory (if necessary) */
  retval = cvProjCreate(&(cv_mem->proj_mem));
  if (retval != CV_SUCCESS)
  {
    cvProcessError(cv_mem, CV_MEM_FAIL, "CVODE", "CVodeSetProjFn",
                   MSG_CV_MEM_FAIL);
    return(CV_MEM_FAIL);
  }

  /* Shortcut to projection memory */
  proj_mem = cv_mem->proj_mem;

  /* User-defined projection */
  proj_mem->internal_proj = SUNFALSE;

  /* Set the projection function */
  proj_mem->pfun = pfun;

  /* Enable projection */
  cv_mem->proj_enabled = SUNTRUE;

  return(CV_SUCCESS);
}


/* ===========================================================================
 * Exported Functions - projection set function
 * ===========================================================================*/


int CVodeSetProjErrEst(void *cvode_mem, booleantype onoff)
{
  int          retval;
  CVodeMem     cv_mem;
  CVodeProjMem proj_mem;

  /* Access memory structures */
  retval = cvAccessProjMem(cvode_mem, "CVodeSetProjErrEst",
                           &cv_mem, &proj_mem);
  if (retval != CV_SUCCESS) return(retval);

  /* Set projection error flag */
  proj_mem->err_proj = onoff;

  return(CV_SUCCESS);
}


int CVodeSetProjFrequency(void *cvode_mem, long int freq)
{
  int          retval;
  CVodeMem     cv_mem;
  CVodeProjMem proj_mem;

  /* Access memory structures */
  retval = cvAccessProjMem(cvode_mem, "CVodeSetProjFrequency",
                           &cv_mem, &proj_mem);
  if (retval != CV_SUCCESS) return(retval);

  /* Set projection frequency */
  if (freq < 0)
  {
    /* Restore default */
    proj_mem->freq  = 1;
    cv_mem->proj_enabled = SUNTRUE;
  }
  else if (freq == 0)
  {
    /* Disable projection */
    proj_mem->freq = 0;
    cv_mem->proj_enabled = SUNFALSE;
  }
  else
  {
    /* Enable projection at given frequency */
    proj_mem->freq = freq;
    cv_mem->proj_enabled = SUNTRUE;
  }

  return(CV_SUCCESS);
}


int CVodeSetMaxNumProjFails(void *cvode_mem, int max_fails)
{
  int          retval;
  CVodeMem     cv_mem;
  CVodeProjMem proj_mem;

  /* Access memory structures */
  retval = cvAccessProjMem(cvode_mem, "CVodeSetMaxNumProjFails",
                           &cv_mem, &proj_mem);
  if (retval != CV_SUCCESS) return(retval);

  /* Set maximum number of projection failures in a step attempt */
  if (max_fails < 1)
  {
    /* Restore default */
    proj_mem->max_fails = PROJ_MAX_FAILS;
  }
  else
  {
    /* Update max number of fails */
    proj_mem->max_fails = max_fails;
  }

  return(CV_SUCCESS);
}


int CVodeSetEpsProj(void *cvode_mem, realtype eps)
{
  int          retval;
  CVodeMem     cv_mem;
  CVodeProjMem proj_mem;

  /* Access memory structures */
  retval = cvAccessProjMem(cvode_mem, "CVodeSetEpsProj",
                           &cv_mem, &proj_mem);
  if (retval != CV_SUCCESS) return(retval);

  /* Set the projection tolerance */
  if (eps <= ZERO)
  {
    /* Restore default */
    proj_mem->eps_proj = PROJ_EPS;
  }
  else
  {
    /* Update projection tolerance */
    proj_mem->eps_proj = eps;
  }

  return(CV_SUCCESS);
}


int CVodeSetProjFailEta(void *cvode_mem, realtype eta)
{
  int          retval;
  CVodeMem     cv_mem;
  CVodeProjMem proj_mem;

  /* Access memory structures */
  retval = cvAccessProjMem(cvode_mem, "CVodeSetProjFailEta",
                           &cv_mem, &proj_mem);
  if (retval != CV_SUCCESS) return(retval);

  /* Set the step size reduction factor for a projection failure */
  if ((eta <= ZERO) || (eta > ONE))
  {
    /* Restore detault */
    proj_mem->eta_pfail = PROJ_FAIL_ETA;
  }
  else
  {
    /* Udpate the eta value */
    proj_mem->eta_pfail = PROJ_FAIL_ETA;
  }

  return(CV_SUCCESS);
}


/* ===========================================================================
 * Exported Functions - projection get functions
 * ===========================================================================*/


int CVodeGetNumProjEvals(void *cvode_mem, long int *nproj)
{
  int          retval;
  CVodeMem     cv_mem;
  CVodeProjMem proj_mem;

  /* Access memory structures */
  retval = cvAccessProjMem(cvode_mem, "CVodeGetNumProjectionEvals",
                           &cv_mem, &proj_mem);
  if (retval != CV_SUCCESS) return(retval);

  /* Get number of projection evaluations */
  *nproj = proj_mem->nproj;

  return(CV_SUCCESS);
}


int CVodeGetNumProjFails(void *cvode_mem, long int *npfails)
{
  int          retval;
  CVodeMem     cv_mem;
  CVodeProjMem proj_mem;

  /* Access memory structures */
  retval = cvAccessProjMem(cvode_mem, "CVodeGetNumProjFails",
                           &cv_mem, &proj_mem);
  if (retval != CV_SUCCESS) return(retval);

  /* Get number of projection fails */
  *npfails = proj_mem->npfails;

  return(CV_SUCCESS);
}


/* ===========================================================================
 * Internal Functions
 * ===========================================================================*/


/*
 * cvProjection
 *
 * For user supplied projection function, use ftemp as temporary storage
 * for the current error estimate (acor) and use tempv to store the
 * accumulated corection due to projection, acorP (tempv is not touched
 * until it is potentially used in cvCompleteStep).
 */

int cvDoProjection(CVodeMem cv_mem, int *nflagPtr, realtype saved_t,
                   int *npfailPtr)
{
  int          retval;
  N_Vector     errP;
  N_Vector     acorP;
  CVodeProjMem proj_mem;

  /* Access projection memory */
  if (cv_mem->proj_mem == NULL) {
    cvProcessError(cv_mem, CV_PROJ_MEM_NULL, "CVODE",
                   "cvDoProjection", MSG_CV_PROJ_MEM_NULL);
    return(CV_PROJ_MEM_NULL);
  }
  proj_mem = cv_mem->proj_mem;

  /* Initialize return flag to success */
  retval = CV_SUCCESS;

  /* Use tempv to store acorP and, if projecting the error, ftemp to store
     errP (recall that in this case we did not allocate vectors to for
     acorP and errP). */
  acorP = cv_mem->cv_tempv;
  if (proj_mem->err_proj)
    errP = cv_mem->cv_ftemp;
  else
    errP = NULL;

  /* Copy acor into errP (if projecting the error) */
  if (proj_mem->err_proj) N_VScale(ONE, cv_mem->cv_acor, errP);

  /* Call the user projection function */
  retval = proj_mem->pfun(cv_mem->cv_tn, cv_mem->cv_y, acorP,
                          proj_mem->eps_proj, errP, cv_mem->cv_user_data);
  proj_mem->nproj++;

  /* This is not the first projection anymore */
  proj_mem->first_proj = SUNFALSE;

  /* Check the return value */
  if (retval == CV_SUCCESS)
  {
    /* Recompute acnrm to be used in error test (if projecting the error) */
    if (proj_mem->err_proj)
      cv_mem->cv_acnrm = N_VWrmsNorm(errP, cv_mem->cv_ewt);

    /* The projection was successful, return now */
    cv_mem->proj_applied = SUNTRUE;
    return(CV_SUCCESS);
  }

  /* The projection failed, update the return value */
  if (retval < 0) retval = CV_PROJFUNC_FAIL;
  if (retval > 0) retval = PROJFUNC_RECVR;

  /* Increment cumulative failure count and restore zn */
  proj_mem->npfails++;
  cvRestore(cv_mem, saved_t);

  /* Return if failed unrecoverably */
  if (retval == CV_PROJFUNC_FAIL) return(CV_PROJFUNC_FAIL);

  /* Recoverable failure, increment failure count for this step attempt */
  (*npfailPtr)++;
  cv_mem->cv_etamax = ONE;

  /* Check for maximum number of failures or |h| = hmin */
  if ((SUNRabs(cv_mem->cv_h) <= cv_mem->cv_hmin * ONEPSM) ||
      (*npfailPtr == proj_mem->max_fails))
  {
    if (retval == PROJFUNC_RECVR) return(CV_REPTD_PROJFUNC_ERR);
  }

  /* Reduce step size; return to reattempt the step */
  cv_mem->cv_eta = SUNMAX(proj_mem->eta_pfail,
                          cv_mem->cv_hmin / SUNRabs(cv_mem->cv_h));
  *nflagPtr = PREV_PROJ_FAIL;
  cvRescale(cv_mem);

  return(PREDICT_AGAIN);
}


int cvProjInit(CVodeProjMem proj_mem)
{
  /* check if projection memory exists */
  if (proj_mem == NULL) return(CV_PROJ_MEM_NULL);

  /* reset flags and counters */
  proj_mem->first_proj = SUNTRUE;
  proj_mem->nstlprj    = 0;
  proj_mem->nproj      = 0;
  proj_mem->npfails    = 0;

  return(CV_SUCCESS);
}


int cvProjFree(CVodeProjMem *proj_mem)
{
  if (*proj_mem == NULL) return(CV_SUCCESS);

  free(*proj_mem);
  *proj_mem = NULL;

  return(CV_SUCCESS);
}


/* ===========================================================================
 * Utility Functions
 * ===========================================================================*/

static int cvProjCreate(CVodeProjMem *proj_mem)
{
  int retval;

  /* Allocate projection memory if necessary, otherwise return success */
  if (*proj_mem == NULL)
  {
    *proj_mem = (CVodeProjMem) malloc(sizeof(struct CVodeProjMemRec));
    if (*proj_mem == NULL) return(CV_MEM_FAIL);

    /* Zero out proj_mem */
    memset(*proj_mem, 0, sizeof(struct CVodeProjMemRec));

    /* Initialize projection variables */
    retval = cvProjSetDefaults(*proj_mem);
    if (retval != CV_SUCCESS) return(retval);
  }

  return(CV_SUCCESS);
}


static int cvProjSetDefaults(CVodeProjMem proj_mem)
{
  if (proj_mem == NULL) return(CV_MEM_FAIL);

  proj_mem->internal_proj = SUNTRUE;
  proj_mem->err_proj      = SUNTRUE;
  proj_mem->first_proj    = SUNTRUE;

  proj_mem->freq    = 1;
  proj_mem->nstlprj = 0;

  proj_mem->max_fails = PROJ_MAX_FAILS;

  proj_mem->pfun = NULL;

  proj_mem->eps_proj  = PROJ_EPS;
  proj_mem->eta_pfail = PROJ_FAIL_ETA;

  proj_mem->nproj   = 0;
  proj_mem->npfails = 0;

  return(CV_SUCCESS);
}


static int cvAccessProjMem(void* cvode_mem, const char *fname,
                           CVodeMem *cv_mem, CVodeProjMem *proj_mem)
{
  /* Access cvode memory */
  if (cvode_mem == NULL)
  {
    cvProcessError(NULL, CV_MEM_NULL, "CVODE",
                   fname, MSG_CV_MEM_NULL);
    return(CV_MEM_NULL);
  }
  *cv_mem = (CVodeMem) cvode_mem;

  /* Access projection memory */
  if ((*cv_mem)->proj_mem == NULL)
  {
    cvProcessError(*cv_mem, CV_PROJ_MEM_NULL, "CVODE",
                   fname, MSG_CV_PROJ_MEM_NULL);
    return(CV_PROJ_MEM_NULL);
  }
  *proj_mem = (CVodeProjMem) (*cv_mem)->proj_mem;

  return(CV_SUCCESS);
}
