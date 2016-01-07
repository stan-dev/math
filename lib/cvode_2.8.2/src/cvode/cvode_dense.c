/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
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
 * This is the impleentation file for the CVDENSE linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvode/cvode_dense.h>
#include "cvode_direct_impl.h"
#include "cvode_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* CVDENSE linit, lsetup, lsolve, and lfree routines */
 
static int cvDenseInit(CVodeMem cv_mem);

static int cvDenseSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, 
                        N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int cvDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ycur, N_Vector fcur);

static void cvDenseFree(CVodeMem cv_mem);

/* Readability Replacements */

#define lmm       (cv_mem->cv_lmm)
#define f         (cv_mem->cv_f)
#define nst       (cv_mem->cv_nst)
#define tn        (cv_mem->cv_tn)
#define h         (cv_mem->cv_h)
#define gamma     (cv_mem->cv_gamma)
#define gammap    (cv_mem->cv_gammap)
#define gamrat    (cv_mem->cv_gamrat)
#define ewt       (cv_mem->cv_ewt)
#define linit     (cv_mem->cv_linit)
#define lsetup    (cv_mem->cv_lsetup)
#define lsolve    (cv_mem->cv_lsolve)
#define lfree     (cv_mem->cv_lfree)
#define lmem      (cv_mem->cv_lmem)
#define vec_tmpl     (cv_mem->cv_tempv)
#define setupNonNull (cv_mem->cv_setupNonNull)

#define mtype     (cvdls_mem->d_type)
#define n         (cvdls_mem->d_n)
#define jacDQ     (cvdls_mem->d_jacDQ)
#define jac       (cvdls_mem->d_djac)
#define M         (cvdls_mem->d_M)
#define lpivots   (cvdls_mem->d_lpivots)
#define savedJ    (cvdls_mem->d_savedJ)
#define nstlj     (cvdls_mem->d_nstlj)
#define nje       (cvdls_mem->d_nje)
#define nfeDQ     (cvdls_mem->d_nfeDQ)
#define J_data    (cvdls_mem->d_J_data)
#define last_flag (cvdls_mem->d_last_flag)
                  
/*
 * -----------------------------------------------------------------
 * CVDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module.  CVDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be cvDenseInit, cvDenseSetup, cvDenseSolve, and cvDenseFree,
 * respectively.  It allocates memory for a structure of type
 * CVDlsMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure.  It sets setupNonNull in (*cvode_mem) to
 * TRUE, and the d_jac field to the default cvDlsDenseDQJac.
 * Finally, it allocates memory for M, savedJ, and lpivots.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int CVDense(void *cvode_mem, long int N)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVDENSE", "CVDense", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVDENSE", "CVDense", MSGD_BAD_NVECTOR);
    return(CVDLS_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  linit  = cvDenseInit;
  lsetup = cvDenseSetup;
  lsolve = cvDenseSolve;
  lfree  = cvDenseFree;

  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(struct CVDlsMemRec));
  if (cvdls_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVDENSE", "CVDense", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Initialize Jacobian-related data */
  jacDQ = TRUE;
  jac = NULL;
  J_data = NULL;

  last_flag = CVDLS_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, savedJ, and pivot array */

  M = NULL;
  M = NewDenseMat(N, N);
  if (M == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVDENSE", "CVDense", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  savedJ = NULL;
  savedJ = NewDenseMat(N, N);
  if (savedJ == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVDENSE", "CVDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  lpivots = NULL;
  lpivots = NewLintArray(N);
  if (lpivots == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVDENSE", "CVDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    DestroyMat(savedJ);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cvdls_mem;

  return(CVDLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * cvDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cvDenseInit(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (jacDQ) {
    jac = cvDlsDenseDQJac;
    J_data = cv_mem;
  } else {
    J_data = cv_mem->cv_user_data;
  }

  last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the dense LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int cvDenseSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, booleantype *jcurPtr, 
                        N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  realtype dgamma;
  long int ier;
  CVDlsMem cvdls_mem;
  int retval;

  cvdls_mem = (CVDlsMem) lmem;
 
  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
 
  dgamma = SUNRabs((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + CVD_MSBJ) ||
         ((convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
         (convfail == CV_FAIL_OTHER);
  jok = !jbad;
 
  if (jok) {

    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    DenseCopy(savedJ, M);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    SetToZero(M);

    retval = jac(n, tn, ypred, fpred, M, J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      cvProcessError(cv_mem, CVDLS_JACFUNC_UNRECVR, "CVDENSE", "cvDenseSetup", MSGD_JACFUNC_FAILED);
      last_flag = CVDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      last_flag = CVDLS_JACFUNC_RECVR;
      return(1);
    }

    DenseCopy(M, savedJ);

  }
  
  /* Scale and add I to get M = I - gamma*J */
  DenseScale(-gamma, M);
  AddIdentity(M);

  /* Do LU factorization of M */
  ier = DenseGETRF(M, lpivots); 

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int cvDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ycur, N_Vector fcur)
{
  CVDlsMem cvdls_mem;
  realtype *bd;

  cvdls_mem = (CVDlsMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(M, lpivots, bd);

  /* If CV_BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }
  
  last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void cvDenseFree(CVodeMem cv_mem)
{
  CVDlsMem  cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  DestroyMat(M);
  DestroyMat(savedJ);
  DestroyArray(lpivots);
  free(cvdls_mem);
  cv_mem->cv_lmem = NULL;
}

