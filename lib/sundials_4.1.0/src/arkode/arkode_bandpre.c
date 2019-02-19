/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * Based off of cvode_bandpre.c by Scott D. Cohen, 
 *      Alan C. Hindmarsh, Radu Serban, and Aaron Collier @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * This file contains implementations of the banded difference
 * quotient Jacobian-based preconditioner and solver routines for
 * use with the ARKLS linear solver interface.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_bandpre_impl.h"
#include "arkode_ls_impl.h"
#include <sundials/sundials_math.h>

#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)


/* Prototypes of ARKBandPrecSetup and ARKBandPrecSolve */
static int ARKBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                            booleantype jok, booleantype *jcurPtr, 
                            realtype gamma, void *bp_data);
static int ARKBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                            N_Vector r, N_Vector z, 
                            realtype gamma, realtype delta,
                            int lr, void *bp_data);

/* Prototype for ARKBandPrecFree */
static int ARKBandPrecFree(ARKodeMem ark_mem);

/* Prototype for difference quotient Jacobian calculation routine */
static int ARKBandPDQJac(ARKBandPrecData pdata,
                         realtype t, N_Vector y, N_Vector fy, 
                         N_Vector ftemp, N_Vector ytemp);


/*---------------------------------------------------------------
 Initialization, Free, and Get Functions
 NOTE: The band linear solver assumes a serial implementation
       of the NVECTOR package. Therefore, ARKBandPrecInit will
       first test for a compatible N_Vector internal 
       representation by checking that the function 
       N_VGetArrayPointer exists.
---------------------------------------------------------------*/
int ARKBandPrecInit(void *arkode_mem, sunindextype N, 
                    sunindextype mu, sunindextype ml)
{
  ARKodeMem       ark_mem;
  ARKLsMem        arkls_mem;
  ARKBandPrecData pdata;
  sunindextype    mup, mlp, storagemu;
  int             retval;
  
  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "ARKBandPrecInit",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Test compatibility of NVECTOR package with the BAND preconditioner */
  if(ark_mem->tempv1->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSG_BP_BAD_NVECTOR);
    return(ARKLS_ILL_INPUT);
  }

  /* Allocate data memory */
  pdata = NULL;
  pdata = (ARKBandPrecData) malloc(sizeof *pdata);
  if (pdata == NULL) {
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSG_BP_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }

  /* Load pointers and bandwidths into pdata block. */
  pdata->arkode_mem = arkode_mem;
  pdata->N = N;
  pdata->mu = mup = SUNMIN(N-1, SUNMAX(0,mu));
  pdata->ml = mlp = SUNMIN(N-1, SUNMAX(0,ml));

  /* Initialize nfeBP counter */
  pdata->nfeBP = 0;

  /* Allocate memory for saved banded Jacobian approximation. */
  pdata->savedJ = NULL;
  pdata->savedJ = SUNBandMatrixStorage(N, mup, mlp, mup);
  if (pdata->savedJ == NULL) {
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSG_BP_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }

  /* Allocate memory for banded preconditioner. */
  storagemu = SUNMIN(N-1, mup+mlp);
  pdata->savedP = NULL;
  pdata->savedP = SUNBandMatrixStorage(N, mup, mlp, storagemu);
  if (pdata->savedP == NULL) {
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSG_BP_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }

  /* Allocate memory for banded linear solver */
  pdata->LS = NULL;
  pdata->LS = SUNLinSol_Band(ark_mem->tempv1, pdata->savedP);
  if (pdata->LS == NULL) {
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSG_BP_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }
  
  /* allocate memory for temporary N_Vectors */
  pdata->tmp1 = NULL;
  pdata->tmp1 = N_VClone(ark_mem->tempv1);
  if (pdata->tmp1 == NULL) {
    SUNLinSolFree(pdata->LS);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSG_BP_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }
  pdata->tmp2 = NULL;
  pdata->tmp2 = N_VClone(ark_mem->tempv1);
  if (pdata->tmp2 == NULL) {
    SUNLinSolFree(pdata->LS);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    N_VDestroy(pdata->tmp1);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSG_BP_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }

  /* initialize band linear solver object */
  retval = SUNLinSolInitialize(pdata->LS);
  if (retval != SUNLS_SUCCESS) {
    SUNLinSolFree(pdata->LS);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    N_VDestroy(pdata->tmp1);
    N_VDestroy(pdata->tmp2);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKBANDPRE", 
                    "ARKBandPrecInit", MSG_BP_SUNLS_FAIL);
    return(ARKLS_SUNLS_FAIL);
  }
  
  /* make sure s_P_data is free from any previous allocations */
  if (arkls_mem->pfree)
    arkls_mem->pfree(ark_mem);

  /* Point to the new P_data field in the LS memory */
  arkls_mem->P_data = pdata;

  /* Attach the pfree function */
  arkls_mem->pfree = ARKBandPrecFree;

  /* Attach preconditioner solve and setup functions */
  retval = arkLSSetPreconditioner(arkode_mem, 
                                  ARKBandPrecSetup, 
                                  ARKBandPrecSolve);
  return(retval);
}


int ARKBandPrecGetWorkSpace(void *arkode_mem, long int *lenrwBP, 
                            long int *leniwBP)
{
  ARKodeMem       ark_mem;
  ARKLsMem        arkls_mem;
  ARKBandPrecData pdata;
  sunindextype    lrw1, liw1;
  long int        lrw, liw;
  int             retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "ARKBandPrecGetWorkSpace",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return immediately if ARKBandPrecData is NULL */
  if (arkls_mem->P_data == NULL) {
    arkProcessError(ark_mem, ARKLS_PMEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecGetWorkSpace", MSG_BP_PMEM_NULL);
    return(ARKLS_PMEM_NULL);
  } 
  pdata = (ARKBandPrecData) arkls_mem->P_data;

  /* sum space requirements for all objects in pdata */
  *leniwBP = 4;
  *lenrwBP = 0;
  if (ark_mem->tempv1->ops->nvspace) {
    N_VSpace(ark_mem->tempv1, &lrw1, &liw1);
    *leniwBP += 2*liw1;
    *lenrwBP += 2*lrw1;
  }
  if (pdata->savedJ->ops->space) {
    retval = SUNMatSpace(pdata->savedJ, &lrw, &liw);
    if (retval == 0) {
      *leniwBP += liw;
      *lenrwBP += lrw;
    }
  }
  if (pdata->savedP->ops->space) {
    retval = SUNMatSpace(pdata->savedP, &lrw, &liw);
    if (retval == 0) {
      *leniwBP += liw;
      *lenrwBP += lrw;
    }
  }
  if (pdata->LS->ops->space) {
    retval = SUNLinSolSpace(pdata->LS, &lrw, &liw);
    if (retval == SUNLS_SUCCESS) {
      *leniwBP += liw;
      *lenrwBP += lrw;
    }
  }

  return(ARKLS_SUCCESS);
}


int ARKBandPrecGetNumRhsEvals(void *arkode_mem, long int *nfevalsBP)
{
  ARKodeMem       ark_mem;
  ARKLsMem        arkls_mem;
  ARKBandPrecData pdata;
  int             retval;

  /* access ARKLsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "ARKBandPrecGetNumRhsEvals",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return immediately if ARKBandPrecData is NULL */
  if (arkls_mem->P_data == NULL) {
    arkProcessError(ark_mem, ARKLS_PMEM_NULL, "ARKBANDPRE", 
                    "ARKBandPrecGetNumRhsEvals", MSG_BP_PMEM_NULL);
    return(ARKLS_PMEM_NULL);
  } 
  pdata = (ARKBandPrecData) arkls_mem->P_data;

  /* set output */
  *nfevalsBP = pdata->nfeBP;

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKBandPrecSetup:

 Together ARKBandPrecSetup and ARKBandPrecSolve use a banded
 difference quotient Jacobian to create a preconditioner.
 ARKBandPrecSetup calculates a new J, if necessary, then
 calculates P = I - gamma*J, and does an LU factorization of P.

 The parameters of ARKBandPrecSetup are as follows:

 t       is the current value of the independent variable.

 y       is the current value of the dependent variable vector,
         namely the predicted value of y(t).

 fy      is the vector f(t,y).

 jok     is an input flag indicating whether Jacobian-related
         data needs to be recomputed, as follows:
           jok == SUNFALSE means recompute Jacobian-related data
                  from scratch.
           jok == SUNTRUE means that Jacobian data from the
                  previous PrecSetup call will be reused
                  (with the current value of gamma).
         A ARKBandPrecSetup call with jok == SUNTRUE should only
         occur after a call with jok == SUNFALSE.

 *jcurPtr is a pointer to an output integer flag which is
          set by ARKBandPrecond as follows:
            *jcurPtr = SUNTRUE if Jacobian data was recomputed.
            *jcurPtr = SUNFALSE if Jacobian data was not recomputed,
                       but saved data was reused.

 gamma   is the scalar appearing in the Newton matrix.

 bp_data is a pointer to preconditoner data (set by ARKBandPrecInit)

 The value to be returned by the ARKBandPrecSetup function is
   0  if successful, or
   1  if the band factorization failed.
---------------------------------------------------------------*/
static int ARKBandPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bp_data)
{
  ARKBandPrecData pdata;
  ARKodeMem ark_mem;
  int retval;
  sunindextype ier;

  /* Assume matrix and lpivots have already been allocated. */
  pdata = (ARKBandPrecData) bp_data;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  if (jok) {

    /* If jok = SUNTRUE, use saved copy of J. */
    *jcurPtr = SUNFALSE;
    retval = SUNMatCopy(pdata->savedJ, pdata->savedP);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                      "ARKBandPrecSetup", MSG_BP_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

  } else {

    /* If jok = SUNFALSE, call ARKBandPDQJac for new J value. */
    *jcurPtr = SUNTRUE;
    retval = SUNMatZero(pdata->savedJ);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                      "ARKBandPrecSetup", MSG_BP_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    retval = ARKBandPDQJac(pdata, t, y, fy, 
                           pdata->tmp1, pdata->tmp2);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                      "ARKBandPrecSetup", MSG_BP_RHSFUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    retval = SUNMatCopy(pdata->savedJ, pdata->savedP);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                      "ARKBandPrecSetup", MSG_BP_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

  }
  
  /* Scale and add identity to get savedP = I - gamma*J. */
  retval = SUNMatScaleAddI(-gamma, pdata->savedP);
  if (retval) {
    arkProcessError(ark_mem, -1, "ARKBANDPRE", 
                    "ARKBandPrecSetup", MSG_BP_SUNMAT_FAIL);
    return(-1);
  }

  /* Do LU factorization of matrix and return error flag */
  ier = SUNLinSolSetup_Band(pdata->LS, pdata->savedP);
  return(ier);
}


/*---------------------------------------------------------------
 ARKBandPrecSolve:

 ARKBandPrecSolve solves a linear system P z = r, where P is the
 matrix computed by ARKBandPrecond.

 The parameters of ARKBandPrecSolve used here are as follows:

 r is the right-hand side vector of the linear system.

 bp_data is a pointer to preconditoner data (set by ARKBandPrecInit)

 z is the output vector computed by ARKBandPrecSolve.

 The value returned by the ARKBandPrecSolve function is always 0,
 indicating success.
---------------------------------------------------------------*/ 
static int ARKBandPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                            N_Vector r, N_Vector z, 
                            realtype gamma, realtype delta,
                            int lr, void *bp_data)
{
  ARKBandPrecData pdata;
  int retval;

  /* Assume matrix and linear solver have already been allocated. */
  pdata = (ARKBandPrecData) bp_data;

  /* Call banded solver object to do the work */
  retval = SUNLinSolSolve(pdata->LS, pdata->savedP, z, r, ZERO);
  return(retval);
}


/*---------------------------------------------------------------
 ARKBandPrecFree:

 Frees data associated with the ARKBand preconditioner.
---------------------------------------------------------------*/ 
static int ARKBandPrecFree(ARKodeMem ark_mem)
{
  ARKLsMem        arkls_mem;
  void*           ark_step_lmem;
  ARKBandPrecData pdata;

  /* Return immediately if ARKodeMem, ARKLsMem or ARKBandPrecData are NULL */
  if (ark_mem == NULL) return(0);
  ark_step_lmem = ark_mem->step_getlinmem((void*) ark_mem);
  if (ark_step_lmem == NULL) return(0);
  arkls_mem = (ARKLsMem) ark_step_lmem;
  if (arkls_mem->P_data == NULL) return(0);
  pdata = (ARKBandPrecData) arkls_mem->P_data;

  SUNLinSolFree(pdata->LS);
  SUNMatDestroy(pdata->savedP);
  SUNMatDestroy(pdata->savedJ);
  N_VDestroy(pdata->tmp1);
  N_VDestroy(pdata->tmp2);

  free(pdata);
  pdata = NULL;

  return(0);
}


/*---------------------------------------------------------------
 ARKBandPDQJac:

 This routine generates a banded difference quotient approximation to
 the Jacobian of f(t,y). It assumes that a band matrix of type
 DlsMat is stored column-wise, and that elements within each column
 are contiguous. This makes it possible to get the address of a column
 of J via the macro BAND_COL and to write a simple for loop to set
 each of the elements of a column in succession.
---------------------------------------------------------------*/
static int ARKBandPDQJac(ARKBandPrecData pdata, 
                         realtype t, N_Vector y, N_Vector fy, 
                         N_Vector ftemp, N_Vector ytemp)
{
  ARKodeMem ark_mem;
  ARKRhsFn fi;
  realtype fnorm, minInc, inc, inc_inv, srur;
  sunindextype group, i, j, width, ngroups, i1, i2;
  realtype *col_j, *ewt_data, *fy_data, *ftemp_data, *y_data, *ytemp_data;
  int retval;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  /* Access implicit RHS function */
  fi = NULL;
  fi = ark_mem->step_getimplicitrhs((void*) ark_mem);
  if (fi == NULL)  return(-1);
  
  /* Obtain pointers to the data for ewt, fy, ftemp, y, ytemp. */
  ewt_data   = N_VGetArrayPointer(ark_mem->ewt);
  fy_data    = N_VGetArrayPointer(fy);
  ftemp_data = N_VGetArrayPointer(ftemp);
  y_data     = N_VGetArrayPointer(y);
  ytemp_data = N_VGetArrayPointer(ytemp);

  /* Load ytemp with y = predicted y vector. */
  N_VScale(ONE, y, ytemp);

  /* Set minimum increment based on uround and norm of f. */
  srur = SUNRsqrt(ark_mem->uround);
  fnorm = N_VWrmsNorm(fy, ark_mem->rwt);
  minInc = (fnorm != ZERO) ?
    (MIN_INC_MULT * SUNRabs(ark_mem->h) * 
     ark_mem->uround * pdata->N * fnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing. */
  width = pdata->ml + pdata->mu + 1;
  ngroups = SUNMIN(width, pdata->N);
  
  for (group = 1; group <= ngroups; group++) {
    
    /* Increment all y_j in group. */
    for(j = group-1; j < pdata->N; j += width) {
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate f with incremented y. */
    retval = fi(t, ytemp, ftemp, ark_mem->user_data);
    pdata->nfeBP++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients. */
    for (j = group-1; j < pdata->N; j += width) {
      ytemp_data[j] = y_data[j];
      col_j = SUNBandMatrix_Column(pdata->savedJ,j);
      inc = SUNMAX(srur*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-pdata->mu);
      i2 = SUNMIN(j+pdata->ml, pdata->N-1);
      for (i=i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j,i,j) =
          inc_inv * (ftemp_data[i] - fy_data[i]);
    }
  }

  return(0);
}


/*---------------------------------------------------------------
   EOF
---------------------------------------------------------------*/ 
