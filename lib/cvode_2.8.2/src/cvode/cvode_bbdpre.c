/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Michael Wittman, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This file contains implementations of routines for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks, for use with CVODE, a CVSPILS linear
 * solver, and the parallel implementation of NVECTOR.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "cvode_bbdpre_impl.h"
#include "cvode_spils_impl.h"

#include <cvode/cvode_sptfqmr.h>
#include <cvode/cvode_spbcgs.h>
#include <cvode/cvode_spgmr.h>

#include <sundials/sundials_math.h>

#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)

/* Prototypes of functions CVBBDPrecSetup and CVBBDPrecSolve */

static int CVBBDPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                          booleantype jok, booleantype *jcurPtr, 
                          realtype gamma, void *bbd_data, 
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int CVBBDPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                          N_Vector r, N_Vector z, 
                          realtype gamma, realtype delta,
                          int lr, void *bbd_data, N_Vector tmp);

/* Prototype for CVBBDPrecFree */
static void CVBBDPrecFree(CVodeMem cv_mem);


/* Prototype for difference quotient Jacobian calculation routine */

static int CVBBDDQJac(CVBBDPrecData pdata, realtype t, 
                      N_Vector y, N_Vector gy, 
                      N_Vector ytemp, N_Vector gtemp);

/* Redability replacements */

#define uround   (cv_mem->cv_uround)
#define vec_tmpl (cv_mem->cv_tempv)

/*
 * -----------------------------------------------------------------
 * User-Callable Functions: initialization, reinit and free
 * -----------------------------------------------------------------
 */

int CVBBDPrecInit(void *cvode_mem, long int Nlocal, 
                   long int mudq, long int mldq,
                   long int mukeep, long int mlkeep, 
                   realtype dqrely, 
                   CVLocalFn gloc, CVCommFn cfn)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBBDPrecData pdata;
  long int muk, mlk, storage_mu;
  int flag;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBBDPRE", "CVBBDPrecInit", MSGBBD_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBBDPRE", "CVBBDPrecInit", MSGBBD_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Test if the NVECTOR package is compatible with the BLOCK BAND preconditioner */
  if(vec_tmpl->ops->nvgetarraypointer == NULL) {
    cvProcessError(cv_mem, CVSPILS_ILL_INPUT, "CVBBDPRE", "CVBBDPrecInit", MSGBBD_BAD_NVECTOR);
    return(CVSPILS_ILL_INPUT);
  }

  /* Allocate data memory */
  pdata = NULL;
  pdata = (CVBBDPrecData) malloc(sizeof *pdata);  
  if (pdata == NULL) {
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBBDPRE", "CVBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Set pointers to gloc and cfn; load half-bandwidths */
  pdata->cvode_mem = cvode_mem;
  pdata->gloc = gloc;
  pdata->cfn = cfn;
  pdata->mudq = SUNMIN(Nlocal-1, SUNMAX(0,mudq));
  pdata->mldq = SUNMIN(Nlocal-1, SUNMAX(0,mldq));
  muk = SUNMIN(Nlocal-1, SUNMAX(0,mukeep));
  mlk = SUNMIN(Nlocal-1, SUNMAX(0,mlkeep));
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* Allocate memory for saved Jacobian */
  pdata->savedJ = NewBandMat(Nlocal, muk, mlk, muk);
  if (pdata->savedJ == NULL) { 
    free(pdata); pdata = NULL; 
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBBDPRE", "CVBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CVSPILS_MEM_FAIL); 
  }

  /* Allocate memory for preconditioner matrix */
  storage_mu = SUNMIN(Nlocal-1, muk + mlk);
  pdata->savedP = NULL;
  pdata->savedP = NewBandMat(Nlocal, muk, mlk, storage_mu);
  if (pdata->savedP == NULL) {
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBBDPRE", "CVBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }
  /* Allocate memory for lpivots */
  pdata->lpivots = NULL;
  pdata->lpivots = NewLintArray(Nlocal);
  if (pdata->lpivots == NULL) {
    DestroyMat(pdata->savedP);
    DestroyMat(pdata->savedJ);
    free(pdata); pdata = NULL;
    cvProcessError(cv_mem, CVSPILS_MEM_FAIL, "CVBBDPRE", "CVBBDPrecInit", MSGBBD_MEM_FAIL);
    return(CVSPILS_MEM_FAIL);
  }

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : SUNRsqrt(uround);

  /* Store Nlocal to be used in CVBBDPrecSetup */
  pdata->n_local = Nlocal;

  /* Set work space sizes and initialize nge */
  pdata->rpwsize = Nlocal*(muk + 2*mlk + storage_mu + 2);
  pdata->ipwsize = Nlocal;
  pdata->nge = 0;

  /* Overwrite the P_data field in the SPILS memory */
  cvspils_mem->s_P_data = pdata;

  /* Attach the pfree function */
  cvspils_mem->s_pfree = CVBBDPrecFree;

  /* Attach preconditioner solve and setup functions */
  flag = CVSpilsSetPreconditioner(cvode_mem, CVBBDPrecSetup, CVBBDPrecSolve);

  return(flag);
}


int CVBBDPrecReInit(void *cvode_mem, 
                    long int mudq, long int mldq, 
                    realtype dqrely)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBBDPrecData pdata;
  long int Nlocal;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBBDPRE", "CVBBDPrecReInit", MSGBBD_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if one of the SPILS linear solvers has been attached */
  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBBDPRE", "CVBBDPrecReInit", MSGBBD_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  /* Test if the preconditioner data is non-NULL */
  if (cvspils_mem->s_P_data == NULL) {
    cvProcessError(cv_mem, CVSPILS_PMEM_NULL, "CVBBDPRE", "CVBBDPrecReInit", MSGBBD_PMEM_NULL);
    return(CVSPILS_PMEM_NULL);
  } 
  pdata = (CVBBDPrecData) cvspils_mem->s_P_data;

  /* Load half-bandwidths */
  Nlocal = pdata->n_local;
  pdata->mudq = SUNMIN(Nlocal-1, SUNMAX(0,mudq));
  pdata->mldq = SUNMIN(Nlocal-1, SUNMAX(0,mldq));

  /* Set pdata->dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? dqrely : SUNRsqrt(uround);

  /* Re-initialize nge */
  pdata->nge = 0;

  return(CVSPILS_SUCCESS);
}

int CVBBDPrecGetWorkSpace(void *cvode_mem, long int *lenrwBBDP, long int *leniwBBDP)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBBDPrecData pdata;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBBDPRE", "CVBBDPrecGetWorkSpace", MSGBBD_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBBDPRE", "CVBBDPrecGetWorkSpace", MSGBBD_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  if (cvspils_mem->s_P_data == NULL) {
    cvProcessError(cv_mem, CVSPILS_PMEM_NULL, "CVBBDPRE", "CVBBDPrecGetWorkSpace", MSGBBD_PMEM_NULL);
    return(CVSPILS_PMEM_NULL);
  } 
  pdata = (CVBBDPrecData) cvspils_mem->s_P_data;

  *lenrwBBDP = pdata->rpwsize;
  *leniwBBDP = pdata->ipwsize;

  return(CVSPILS_SUCCESS);
}

int CVBBDPrecGetNumGfnEvals(void *cvode_mem, long int *ngevalsBBDP)
{
  CVodeMem cv_mem;
  CVSpilsMem cvspils_mem;
  CVBBDPrecData pdata;

  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSPILS_MEM_NULL, "CVBBDPRE", "CVBBDPrecGetNumGfnEvals", MSGBBD_MEM_NULL);
    return(CVSPILS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (cv_mem->cv_lmem == NULL) {
    cvProcessError(cv_mem, CVSPILS_LMEM_NULL, "CVBBDPRE", "CVBBDPrecGetNumGfnEvals", MSGBBD_LMEM_NULL);
    return(CVSPILS_LMEM_NULL);
  }
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;

  if (cvspils_mem->s_P_data == NULL) {
    cvProcessError(cv_mem, CVSPILS_PMEM_NULL, "CVBBDPRE", "CVBBDPrecGetNumGfnEvals", MSGBBD_PMEM_NULL);
    return(CVSPILS_PMEM_NULL);
  } 
  pdata = (CVBBDPrecData) cvspils_mem->s_P_data;

  *ngevalsBBDP = pdata->nge;

  return(CVSPILS_SUCCESS);
}

/* Readability Replacements */

#define Nlocal  (pdata->n_local)
#define mudq    (pdata->mudq)
#define mldq    (pdata->mldq)
#define mukeep  (pdata->mukeep)
#define mlkeep  (pdata->mlkeep)
#define dqrely  (pdata->dqrely)
#define gloc    (pdata->gloc)
#define cfn     (pdata->cfn)
#define savedJ  (pdata->savedJ)
#define savedP  (pdata->savedP)
#define lpivots (pdata->lpivots)
#define nge     (pdata->nge)

/*
 * -----------------------------------------------------------------
 * Function : CVBBDPrecSetup                                      
 * -----------------------------------------------------------------
 * CVBBDPrecSetup generates and factors a banded block of the
 * preconditioner matrix on each processor, via calls to the
 * user-supplied gloc and cfn functions. It uses difference
 * quotient approximations to the Jacobian elements.
 *
 * CVBBDPrecSetup calculates a new J,if necessary, then calculates
 * P = I - gamma*J, and does an LU factorization of P.
 *
 * The parameters of CVBBDPrecSetup used here are as follows:
 *
 * t       is the current value of the independent variable.
 *
 * y       is the current value of the dependent variable vector,
 *         namely the predicted value of y(t).
 *
 * fy      is the vector f(t,y).
 *
 * jok     is an input flag indicating whether Jacobian-related
 *         data needs to be recomputed, as follows:
 *           jok == FALSE means recompute Jacobian-related data
 *                  from scratch.
 *           jok == TRUE  means that Jacobian data from the
 *                  previous CVBBDPrecon call can be reused
 *                  (with the current value of gamma).
 *         A CVBBDPrecon call with jok == TRUE should only occur
 *         after a call with jok == FALSE.
 *
 * jcurPtr is a pointer to an output integer flag which is
 *         set by CVBBDPrecon as follows:
 *           *jcurPtr = TRUE if Jacobian data was recomputed.
 *           *jcurPtr = FALSE if Jacobian data was not recomputed,
 *                      but saved data was reused.
 *
 * gamma   is the scalar appearing in the Newton matrix.
 *
 * bbd_data is a pointer to the preconditioner data set by
 *          CVBBDPrecInit
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated
 *           for NVectors which are be used by CVBBDPrecSetup
 *           as temporary storage or work space.
 *
 * Return value:
 * The value returned by this CVBBDPrecSetup function is the int
 *   0  if successful,
 *   1  for a recoverable error (step will be retried).
 * -----------------------------------------------------------------
 */

static int CVBBDPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                          booleantype jok, booleantype *jcurPtr, 
                          realtype gamma, void *bbd_data, 
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  long int ier;
  CVBBDPrecData pdata;
  CVodeMem cv_mem;
  int retval;

  pdata = (CVBBDPrecData) bbd_data;

  cv_mem = (CVodeMem) pdata->cvode_mem;

  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    BandCopy(savedJ, savedP, mukeep, mlkeep);

  } else {

    /* Otherwise call CVBBDDQJac for new J value */
    *jcurPtr = TRUE;
    SetToZero(savedJ);

    retval = CVBBDDQJac(pdata, t, y, tmp1, tmp2, tmp3);
    if (retval < 0) {
      cvProcessError(cv_mem, -1, "CVBBDPRE", "CVBBDPrecSetup", MSGBBD_FUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    BandCopy(savedJ, savedP, mukeep, mlkeep);

  }
  
  /* Scale and add I to get P = I - gamma*J */
  BandScale(-gamma, savedP);
  AddIdentity(savedP);
 
  /* Do LU factorization of P in place */
  ier = BandGBTRF(savedP, lpivots);
 
  /* Return 0 if the LU was complete; otherwise return 1 */
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * Function : CVBBDPrecSolve
 * -----------------------------------------------------------------
 * CVBBDPrecSolve solves a linear system P z = r, with the
 * band-block-diagonal preconditioner matrix P generated and
 * factored by CVBBDPrecSetup.
 *
 * The parameters of CVBBDPrecSolve used here are as follows:
 *
 * r is the right-hand side vector of the linear system.
 *
 * bbd_data is a pointer to the preconditioner data set by
 *   CVBBDPrecInit.
 *
 * z is the output vector computed by CVBBDPrecSolve.
 *
 * The value returned by the CVBBDPrecSolve function is always 0,
 * indicating success.
 * -----------------------------------------------------------------
 */

static int CVBBDPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                          N_Vector r, N_Vector z, 
                          realtype gamma, realtype delta,
                          int lr, void *bbd_data, N_Vector tmp)
{
  CVBBDPrecData pdata;
  realtype *zd;

  pdata = (CVBBDPrecData) bbd_data;

  /* Copy r to z, then do backsolve and return */
  N_VScale(ONE, r, z);
  
  zd = N_VGetArrayPointer(z);

  BandGBTRS(savedP, lpivots, zd);

  return(0);
}


static void CVBBDPrecFree(CVodeMem cv_mem)
{
  CVSpilsMem cvspils_mem;
  CVBBDPrecData pdata;
  
  if (cv_mem->cv_lmem == NULL) return;
  cvspils_mem = (CVSpilsMem) cv_mem->cv_lmem;
  
  if (cvspils_mem->s_P_data == NULL) return;
  pdata = (CVBBDPrecData) cvspils_mem->s_P_data;

  DestroyMat(savedJ);
  DestroyMat(savedP);
  DestroyArray(lpivots);

  free(pdata);
  pdata = NULL;
}


#define ewt       (cv_mem->cv_ewt)
#define h         (cv_mem->cv_h)
#define user_data (cv_mem->cv_user_data)

/*
 * -----------------------------------------------------------------
 * Function : CVBBDDQJac
 * -----------------------------------------------------------------
 * This routine generates a banded difference quotient approximation
 * to the local block of the Jacobian of g(t,y). It assumes that a
 * band matrix of type DlsMat is stored columnwise, and that elements
 * within each column are contiguous. All matrix elements are generated
 * as difference quotients, by way of calls to the user routine gloc.
 * By virtue of the band structure, the number of these calls is
 * bandwidth + 1, where bandwidth = mldq + mudq + 1.
 * But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 * This routine also assumes that the local elements of a vector are
 * stored contiguously.
 * -----------------------------------------------------------------
 */

static int CVBBDDQJac(CVBBDPrecData pdata, realtype t, 
                      N_Vector y, N_Vector gy, 
                      N_Vector ytemp, N_Vector gtemp)
{
  CVodeMem cv_mem;
  realtype gnorm, minInc, inc, inc_inv;
  long int group, i, j, width, ngroups, i1, i2;
  realtype *y_data, *ewt_data, *gy_data, *gtemp_data, *ytemp_data, *col_j;
  int retval;

  cv_mem = (CVodeMem) pdata->cvode_mem;

  /* Load ytemp with y = predicted solution vector */
  N_VScale(ONE, y, ytemp);

  /* Call cfn and gloc to get base value of g(t,y) */
  if (cfn != NULL) {
    retval = cfn(Nlocal, t, y, user_data);
    if (retval != 0) return(retval);
  }

  retval = gloc(Nlocal, t, ytemp, gy, user_data);
  nge++;
  if (retval != 0) return(retval);

  /* Obtain pointers to the data for various vectors */
  y_data     =  N_VGetArrayPointer(y);
  gy_data    =  N_VGetArrayPointer(gy);
  ewt_data   =  N_VGetArrayPointer(ewt);
  ytemp_data =  N_VGetArrayPointer(ytemp);
  gtemp_data =  N_VGetArrayPointer(gtemp);

  /* Set minimum increment based on uround and norm of g */
  gnorm = N_VWrmsNorm(gy, ewt);
  minInc = (gnorm != ZERO) ?
           (MIN_INC_MULT * SUNRabs(h) * uround * Nlocal * gnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = mldq + mudq + 1;
  ngroups = SUNMIN(width, Nlocal);

  /* Loop over groups */  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < Nlocal; j+=width) {
      inc = SUNMAX(dqrely*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate g with incremented y */
    retval = gloc(Nlocal, t, ytemp, gtemp, user_data);
    nge++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < Nlocal; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = BAND_COL(savedJ,j);
      inc = SUNMAX(dqrely*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-mukeep);
      i2 = SUNMIN(j+mlkeep, Nlocal-1);
      for (i=i1; i <= i2; i++)
        BAND_COL_ELEM(col_j,i,j) =
          inc_inv * (gtemp_data[i] - gy_data[i]);
    }
  }

  return(0);
}
