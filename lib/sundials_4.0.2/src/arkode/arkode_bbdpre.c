/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * This file contains implementations of routines for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks, for use with ARKode, the ARKLS 
 * linear solver interface, and the MPI-parallel implementation 
 * of NVECTOR.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include "arkode_bbdpre_impl.h"
#include "arkode_ls_impl.h"
#include <sundials/sundials_math.h>
#include <nvector/nvector_serial.h>


#define MIN_INC_MULT RCONST(1000.0)
#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)


/* Prototypes of functions ARKBBDPrecSetup and ARKBBDPrecSolve */
static int ARKBBDPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bbd_data);
static int ARKBBDPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bbd_data);

/* Prototype for ARKBBDPrecFree */
static int ARKBBDPrecFree(ARKodeMem ark_mem);

/* Prototype for difference quotient Jacobian calculation routine */
static int ARKBBDDQJac(ARKBBDPrecData pdata, realtype t, 
                       N_Vector y, N_Vector gy, 
                       N_Vector ytemp, N_Vector gtemp);


/*---------------------------------------------------------------
 User-Callable Functions: initialization, reinit and free
---------------------------------------------------------------*/
int ARKBBDPrecInit(void *arkode_mem, sunindextype Nlocal, 
                   sunindextype mudq, sunindextype mldq,
                   sunindextype mukeep, sunindextype mlkeep, 
                   realtype dqrely, 
                   ARKLocalFn gloc, ARKCommFn cfn)
{
  ARKodeMem      ark_mem;
  ARKLsMem       arkls_mem;
  ARKBBDPrecData pdata;
  sunindextype   muk, mlk, storage_mu, lrw1, liw1;
  long int       lrw, liw;
  int            retval;

  /* access ARKMilsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "ARKBBDPrecInit",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Test compatibility of NVECTOR package with the BBD preconditioner */
  if(ark_mem->tempv1->ops->nvgetarraypointer == NULL) {
    arkProcessError(ark_mem, ARKLS_ILL_INPUT, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_BAD_NVECTOR);
    return(ARKLS_ILL_INPUT);
  }

  /* Allocate data memory */
  pdata = NULL;
  pdata = (ARKBBDPrecData) malloc(sizeof *pdata);  
  if (pdata == NULL) {
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }

  /* Set pointers to gloc and cfn; load half-bandwidths */
  pdata->arkode_mem = arkode_mem;
  pdata->gloc = gloc;
  pdata->cfn = cfn;
  pdata->mudq = SUNMIN(Nlocal-1, SUNMAX(0,mudq));
  pdata->mldq = SUNMIN(Nlocal-1, SUNMAX(0,mldq));
  muk = SUNMIN(Nlocal-1, SUNMAX(0,mukeep));
  mlk = SUNMIN(Nlocal-1, SUNMAX(0,mlkeep));
  pdata->mukeep = muk;
  pdata->mlkeep = mlk;

  /* Allocate memory for saved Jacobian */
  pdata->savedJ = SUNBandMatrixStorage(Nlocal, muk, mlk, muk);
  if (pdata->savedJ == NULL) { 
    free(pdata); pdata = NULL; 
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL); 
  }

  /* Allocate memory for preconditioner matrix */
  storage_mu = SUNMIN(Nlocal-1, muk + mlk);
  pdata->savedP = NULL;
  pdata->savedP = SUNBandMatrixStorage(Nlocal, muk, mlk, storage_mu);
  if (pdata->savedP == NULL) {
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }

  /* Allocate memory for temporary N_Vectors */
  pdata->zlocal = NULL;
  pdata->zlocal = N_VNewEmpty_Serial(Nlocal);
  if (pdata->zlocal == NULL) {
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }
  pdata->rlocal = NULL;
  pdata->rlocal = N_VNewEmpty_Serial(Nlocal);
  if (pdata->rlocal == NULL) {
    N_VDestroy(pdata->zlocal);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }
  pdata->tmp1 = NULL;
  pdata->tmp1 = N_VClone(ark_mem->tempv1);
  if (pdata->tmp1 == NULL) {
    N_VDestroy(pdata->zlocal);
    N_VDestroy(pdata->rlocal);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }
  pdata->tmp2 = NULL;
  pdata->tmp2 = N_VClone(ark_mem->tempv1);
  if (pdata->tmp2 == NULL) {
    N_VDestroy(pdata->tmp1);
    N_VDestroy(pdata->zlocal);
    N_VDestroy(pdata->rlocal);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }
  pdata->tmp3 = NULL;
  pdata->tmp3 = N_VClone(ark_mem->tempv1);
  if (pdata->tmp3 == NULL) {
    N_VDestroy(pdata->tmp1);
    N_VDestroy(pdata->tmp2);
    N_VDestroy(pdata->zlocal);
    N_VDestroy(pdata->rlocal);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }

  /* Allocate memory for banded linear solver */
  pdata->LS = NULL;
  pdata->LS = SUNLinSol_Band(pdata->rlocal, pdata->savedP);
  if (pdata->LS == NULL) {
    N_VDestroy(pdata->tmp1);
    N_VDestroy(pdata->tmp2);
    N_VDestroy(pdata->tmp3);
    N_VDestroy(pdata->zlocal);
    N_VDestroy(pdata->rlocal);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_MEM_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_MEM_FAIL);
    return(ARKLS_MEM_FAIL);
  }

  /* initialize band linear solver object */
  retval = SUNLinSolInitialize(pdata->LS);
  if (pdata->LS == NULL) {
    N_VDestroy(pdata->tmp1);
    N_VDestroy(pdata->tmp2);
    N_VDestroy(pdata->tmp3);
    N_VDestroy(pdata->zlocal);
    N_VDestroy(pdata->rlocal);
    SUNMatDestroy(pdata->savedP);
    SUNMatDestroy(pdata->savedJ);
    SUNLinSolFree(pdata->LS);
    free(pdata); pdata = NULL;
    arkProcessError(ark_mem, ARKLS_SUNLS_FAIL, "ARKBBDPRE", 
                    "ARKBBDPrecInit", MSG_BBD_SUNLS_FAIL);
    return(ARKLS_SUNLS_FAIL);
  }
  
  /* Set dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? 
    dqrely : SUNRsqrt(ark_mem->uround);

  /* Store Nlocal to be used in ARKBBDPrecSetup */
  pdata->n_local = Nlocal;

  /* Set work space sizes and initialize nge */
  pdata->rpwsize = 0;
  pdata->ipwsize = 0;
  if (ark_mem->tempv1->ops->nvspace) {
    N_VSpace(ark_mem->tempv1, &lrw1, &liw1);
    pdata->rpwsize += 3*lrw1;
    pdata->ipwsize += 3*liw1;
  }
  if (pdata->rlocal->ops->nvspace) {
    N_VSpace(pdata->rlocal, &lrw1, &liw1);
    pdata->rpwsize += 2*lrw1;
    pdata->ipwsize += 2*liw1;
  }
  if (pdata->savedJ->ops->space) {
    retval = SUNMatSpace(pdata->savedJ, &lrw, &liw);
    pdata->rpwsize += lrw;
    pdata->ipwsize += liw;
  }
  if (pdata->savedP->ops->space) {
    retval = SUNMatSpace(pdata->savedP, &lrw, &liw);
    pdata->rpwsize += lrw;
    pdata->ipwsize += liw;
  }
  if (pdata->LS->ops->space) {
    retval = SUNLinSolSpace(pdata->LS, &lrw, &liw);
    pdata->rpwsize += lrw;
    pdata->ipwsize += liw;
  }
  pdata->nge = 0;

  /* make sure P_data is free from any previous allocations */
  if (arkls_mem->pfree) 
    arkls_mem->pfree(ark_mem);

  /* Point to the new P_data field in the LS memory */
  arkls_mem->P_data = pdata;

  /* Attach the pfree function */
  arkls_mem->pfree = ARKBBDPrecFree;

  /* Attach preconditioner solve and setup functions */
  retval = arkLSSetPreconditioner(arkode_mem, 
                                  ARKBBDPrecSetup, 
                                  ARKBBDPrecSolve);

  return(retval);
}


/*-------------------------------------------------------------*/
int ARKBBDPrecReInit(void *arkode_mem, sunindextype mudq, 
                     sunindextype mldq, realtype dqrely)
{
  ARKodeMem      ark_mem;
  ARKLsMem       arkls_mem;
  ARKBBDPrecData pdata;
  sunindextype   Nlocal;
  int            retval;
  
  /* access ARKMilsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "ARKBBDPrecReInit",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return immediately ARKBBDPrecData is NULL */
  if (arkls_mem->P_data == NULL) {
    arkProcessError(ark_mem, ARKLS_PMEM_NULL, "ARKBBDPRE", 
                    "ARKBBDPrecReInit", MSG_BBD_PMEM_NULL);
    return(ARKLS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkls_mem->P_data;

  /* Load half-bandwidths */
  Nlocal = pdata->n_local;
  pdata->mudq = SUNMIN(Nlocal-1, SUNMAX(0,mudq));
  pdata->mldq = SUNMIN(Nlocal-1, SUNMAX(0,mldq));

  /* Set dqrely based on input dqrely (0 implies default). */
  pdata->dqrely = (dqrely > ZERO) ? 
    dqrely : SUNRsqrt(ark_mem->uround);

  /* Re-initialize nge */
  pdata->nge = 0;

  return(ARKLS_SUCCESS);
}


/*-------------------------------------------------------------*/
int ARKBBDPrecGetWorkSpace(void *arkode_mem, 
                           long int *lenrwBBDP, 
                           long int *leniwBBDP)
{
  ARKodeMem      ark_mem;
  ARKLsMem       arkls_mem;
  ARKBBDPrecData pdata;
  int            retval;
  
  /* access ARKMilsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "ARKBBDPrecGetWorkSpace",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return immediately ARKBBDPrecData is NULL */
  if (arkls_mem->P_data == NULL) {
    arkProcessError(ark_mem, ARKLS_PMEM_NULL, "ARKBBDPRE", 
                    "ARKBBDPrecGetWorkSpace", MSG_BBD_PMEM_NULL);
    return(ARKLS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkls_mem->P_data;

  /* set outputs */
  *lenrwBBDP = pdata->rpwsize;
  *leniwBBDP = pdata->ipwsize;

  return(ARKLS_SUCCESS);
}


/*-------------------------------------------------------------*/
int ARKBBDPrecGetNumGfnEvals(void *arkode_mem, 
                             long int *ngevalsBBDP)
{
  ARKodeMem      ark_mem;
  ARKLsMem       arkls_mem;
  ARKBBDPrecData pdata;
  int            retval;
  
  /* access ARKMilsMem structure */
  retval = arkLs_AccessLMem(arkode_mem, "ARKBBDPrecGetNumGfnEvals",
                            &ark_mem, &arkls_mem);
  if (retval != ARK_SUCCESS)  return(retval);

  /* Return immediately if ARKBBDPrecData is NULL */
  if (arkls_mem->P_data == NULL) {
    arkProcessError(ark_mem, ARKLS_PMEM_NULL, "ARKBBDPRE", 
                    "ARKBBDPrecGetNumGfnEvals", MSG_BBD_PMEM_NULL);
    return(ARKLS_PMEM_NULL);
  } 
  pdata = (ARKBBDPrecData) arkls_mem->P_data;

  /* set output */
  *ngevalsBBDP = pdata->nge;

  return(ARKLS_SUCCESS);
}


/*---------------------------------------------------------------
 ARKBBDPrecSetup:

 ARKBBDPrecSetup generates and factors a banded block of the
 preconditioner matrix on each processor, via calls to the
 user-supplied gloc and cfn functions. It uses difference
 quotient approximations to the Jacobian elements.

 ARKBBDPrecSetup calculates a new J, if necessary, then 
 calculates P = M - gamma*J, and does an LU factorization of P.

 The parameters of ARKBBDPrecSetup used here are as follows:

 t       is the current value of the independent variable.

 y       is the current value of the dependent variable vector,
         namely the predicted value of y(t).

 fy      is the vector f(t,y).

 jok     is an input flag indicating whether Jacobian-related
         data needs to be recomputed, as follows:
           jok == SUNFALSE means recompute Jacobian-related data
                  from scratch.
           jok == SUNTRUE  means that Jacobian data from the
                  previous ARKBBDPrecon call can be reused
                  (with the current value of gamma).
         A ARKBBDPrecon call with jok == SUNTRUE should only occur
         after a call with jok == SUNFALSE.

 jcurPtr is a pointer to an output integer flag which is
         set by ARKBBDPrecon as follows:
           *jcurPtr = SUNTRUE if Jacobian data was recomputed.
           *jcurPtr = SUNFALSE if Jacobian data was not recomputed,
                      but saved data was reused.

 gamma   is the scalar appearing in the Newton matrix.

 bbd_data is a pointer to the preconditioner data set by
          ARKBBDPrecInit

 Return value:
 The value returned by this ARKBBDPrecSetup function is the int
   0  if successful,
   1  for a recoverable error (step will be retried).
---------------------------------------------------------------*/
static int ARKBBDPrecSetup(realtype t, N_Vector y, N_Vector fy, 
                           booleantype jok, booleantype *jcurPtr, 
                           realtype gamma, void *bbd_data)
{
  sunindextype ier;
  ARKBBDPrecData pdata;
  ARKodeMem ark_mem;
  int retval;

  pdata = (ARKBBDPrecData) bbd_data;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  /* If jok = SUNTRUE, use saved copy of J */
  if (jok) {
    *jcurPtr = SUNFALSE;
    retval = SUNMatCopy(pdata->savedJ, pdata->savedP);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBBDPRE", 
                      "ARKBBDPrecSetup", MSG_BBD_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

  /* Otherwise call ARKBBDDQJac for new J value */
  } else {
    
    *jcurPtr = SUNTRUE;
    retval = SUNMatZero(pdata->savedJ);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBBDPRE", 
                      "ARKBBDPrecSetup", MSG_BBD_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    retval = ARKBBDDQJac(pdata, t, y, pdata->tmp1, 
                         pdata->tmp2, pdata->tmp3);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBBDPRE", "ARKBBDPrecSetup", 
                      MSG_BBD_FUNC_FAILED);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

    retval = SUNMatCopy(pdata->savedJ, pdata->savedP);
    if (retval < 0) {
      arkProcessError(ark_mem, -1, "ARKBBDPRE", 
                      "ARKBBDPrecSetup", MSG_BBD_SUNMAT_FAIL);
      return(-1);
    }
    if (retval > 0) {
      return(1);
    }

  }
  
  /* Scale and add I to get P = I - gamma*J */
  retval = SUNMatScaleAddI(-gamma, pdata->savedP);
  if (retval) {
    arkProcessError(ark_mem, -1, "ARKBBDPRE", 
                    "ARKBBDPrecSetup", MSG_BBD_SUNMAT_FAIL);
    return(-1);
  }
 
  /* Do LU factorization of matrix and return error flag */
  ier = SUNLinSolSetup_Band(pdata->LS, pdata->savedP);
  return(ier);
}


/*---------------------------------------------------------------
 ARKBBDPrecSolve:

 ARKBBDPrecSolve solves a linear system P z = r, with the
 band-block-diagonal preconditioner matrix P generated and
 factored by ARKBBDPrecSetup.

 The parameters of ARKBBDPrecSolve used here are as follows:

 r is the right-hand side vector of the linear system.

 bbd_data is a pointer to the preconditioner data set by
   ARKBBDPrecInit.

 z is the output vector computed by ARKBBDPrecSolve.

 The value returned by the ARKBBDPrecSolve function is the same 
 as the value returned from the linear solver object.
---------------------------------------------------------------*/
static int ARKBBDPrecSolve(realtype t, N_Vector y, N_Vector fy, 
                           N_Vector r, N_Vector z, 
                           realtype gamma, realtype delta,
                           int lr, void *bbd_data)
{
  int retval;
  ARKBBDPrecData pdata;

  pdata = (ARKBBDPrecData) bbd_data;

  /* Attach local data arrays for r and z to rlocal and zlocal */
  N_VSetArrayPointer(N_VGetArrayPointer(r), pdata->rlocal);
  N_VSetArrayPointer(N_VGetArrayPointer(z), pdata->zlocal);
  
  /* Call banded solver object to do the work */
  retval = SUNLinSolSolve(pdata->LS, pdata->savedP, pdata->zlocal, 
                          pdata->rlocal, ZERO);

  /* Detach local data arrays from rlocal and zlocal */
  N_VSetArrayPointer(NULL, pdata->rlocal);
  N_VSetArrayPointer(NULL, pdata->zlocal);

  return(retval);
}


/*-------------------------------------------------------------*/
static int ARKBBDPrecFree(ARKodeMem ark_mem)
{
  ARKLsMem       arkls_mem;
  void*          ark_step_lmem;
  ARKBBDPrecData pdata;
  
  /* Return immediately if ARKodeMem, ARKLsMem or ARKBandPrecData are NULL */
  if (ark_mem == NULL) return(0);
  ark_step_lmem = ark_mem->step_getlinmem((void*) ark_mem);
  if (ark_step_lmem == NULL) return(0);
  arkls_mem = (ARKLsMem) ark_step_lmem;
  if (arkls_mem->P_data == NULL) return(0);
  pdata = (ARKBBDPrecData) arkls_mem->P_data;

  SUNLinSolFree(pdata->LS);
  N_VDestroy(pdata->tmp1);
  N_VDestroy(pdata->tmp2);
  N_VDestroy(pdata->tmp3);
  N_VDestroy(pdata->zlocal);
  N_VDestroy(pdata->rlocal);
  SUNMatDestroy(pdata->savedP);
  SUNMatDestroy(pdata->savedJ);

  free(pdata);
  pdata = NULL;

  return(0);
}


/*---------------------------------------------------------------
 ARKBBDDQJac:

 This routine generates a banded difference quotient approximation
 to the local block of the Jacobian of g(t,y). It assumes that a
 band matrix of type SUNMatrix is stored columnwise, and that 
 elements within each column are contiguous. All matrix elements 
 are generated as difference quotients, by way of calls to the 
 user routine gloc.  By virtue of the band structure, the number 
 of these calls is bandwidth + 1, where bandwidth = mldq + mudq + 1.
 But the band matrix kept has bandwidth = mlkeep + mukeep + 1.
 This routine also assumes that the local elements of a vector are
 stored contiguously.
---------------------------------------------------------------*/
static int ARKBBDDQJac(ARKBBDPrecData pdata, realtype t, 
                       N_Vector y, N_Vector gy, 
                       N_Vector ytemp, N_Vector gtemp)
{
  ARKodeMem ark_mem;
  realtype gnorm, minInc, inc, inc_inv;
  sunindextype group, i, j, width, ngroups, i1, i2;
  realtype *y_data, *ewt_data, *gy_data, *gtemp_data, *ytemp_data, *col_j;
  int retval;

  ark_mem = (ARKodeMem) pdata->arkode_mem;

  /* Load ytemp with y = predicted solution vector */
  N_VScale(ONE, y, ytemp);

  /* Call cfn and gloc to get base value of g(t,y) */
  if (pdata->cfn != NULL) {
    retval = pdata->cfn(pdata->n_local, t, y, ark_mem->user_data);
    if (retval != 0) return(retval);
  }

  retval = pdata->gloc(pdata->n_local, t, ytemp, gy, 
                       ark_mem->user_data);
  pdata->nge++;
  if (retval != 0) return(retval);

  /* Obtain pointers to the data for various vectors */
  y_data     =  N_VGetArrayPointer(y);
  gy_data    =  N_VGetArrayPointer(gy);
  ewt_data   =  N_VGetArrayPointer(ark_mem->ewt);
  ytemp_data =  N_VGetArrayPointer(ytemp);
  gtemp_data =  N_VGetArrayPointer(gtemp);

  /* Set minimum increment based on uround and norm of g */
  gnorm = N_VWrmsNorm(gy, ark_mem->rwt);
  minInc = (gnorm != ZERO) ? 
    (MIN_INC_MULT * SUNRabs(ark_mem->h) *
     ark_mem->uround * pdata->n_local * gnorm) : ONE;

  /* Set bandwidth and number of column groups for band differencing */
  width = pdata->mldq + pdata->mudq + 1;
  ngroups = SUNMIN(width, pdata->n_local);

  /* Loop over groups */  
  for (group=1; group <= ngroups; group++) {
    
    /* Increment all y_j in group */
    for(j=group-1; j < pdata->n_local; j+=width) {
      inc = SUNMAX(pdata->dqrely*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      ytemp_data[j] += inc;
    }

    /* Evaluate g with incremented y */
    retval = pdata->gloc(pdata->n_local, t, ytemp, gtemp, 
                         ark_mem->user_data);
    pdata->nge++;
    if (retval != 0) return(retval);

    /* Restore ytemp, then form and load difference quotients */
    for (j=group-1; j < pdata->n_local; j+=width) {
      ytemp_data[j] = y_data[j];
      col_j = SUNBandMatrix_Column(pdata->savedJ,j);
      inc = SUNMAX(pdata->dqrely*SUNRabs(y_data[j]), minInc/ewt_data[j]);
      inc_inv = ONE/inc;
      i1 = SUNMAX(0, j-pdata->mukeep);
      i2 = SUNMIN(j+pdata->mlkeep, pdata->n_local-1);
      for (i=i1; i <= i2; i++)
        SM_COLUMN_ELEMENT_B(col_j,i,j) =
          inc_inv * (gtemp_data[i] - gy_data[i]);
    }
  }

  return(0);
}



/*---------------------------------------------------------------
    EOF
---------------------------------------------------------------*/
