
/*
 * -----------------------------------------------------------------
 * $Revision: 4435 $
 * $Date: 2015-03-23 18:26:14 -0700 (Mon, 23 Mar 2015) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Carol S. Woodward @ LLNL
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
 * This is the implementation file for the CVSUPERLUMT linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "cvode_sparse_impl.h"
#include "cvode/cvode_superlumt.h"
#include "sundials/sundials_superlumt_impl.h"
#include "sundials/sundials_math.h"

/* Constants */

#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* CVSUPERLUMT linit, lsetup, lsolve, and lfree routines */
 
static int cvSuperLUMTInit(CVodeMem cv_mem);

static int cvSuperLUMTSetup(CVodeMem cv_mem, int convfail, N_Vector ypred, 
			    N_Vector fpred, booleantype *jcurPtr,
			    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int cvSuperLUMTSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
			    N_Vector ycur, N_Vector fcur);

static void cvSuperLUMTFree(CVodeMem cv_mem);

/*
 * -----------------------------------------------------------------
 * CVSuperLUMT
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the CVODE / SuperLUMT linear solver module.  
 * CVSUPERLUMT first calls the existing lfree routine if this is not NULL.
 * Then it sets the cv_linit, cv_lsetup, cv_lsolve, and
 * cv_lfree fields in (*cv_mem) to be cvSuperLUMTInit, cvSuperLUMTSetup,
 * cvSuperLUMTSolve, and cvSuperLUMTFree, respectively.
 * It allocates memory for a structure of type CVsluMemRec and sets
 * the cv_lmem field in (*cvode_mem) to the address of this structure.
 * It sets setupNonNull in (*cvode_mem) to TRUE.
 * Finally, it allocates memory for SuperLUMT.
 * The return value is CVSLS_SUCCESS = 0, CVSLS_LMEM_FAIL = -1,
 * or CVSLS_ILL_INPUT = -2.
 *
 * NOTE: The SuperLUMT linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVSuperLUMT will first 
 *       test for a compatible N_Vector internal representation
 *       by checking that the function N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */

int CVSuperLUMT(void *cvode_mem, int num_threads, int n, int nnz)
{
  CVodeMem cv_mem;
  CVSlsMem cvsls_mem;
  SLUMTData slumt_data;
  int *perm_c, *perm_r;
  int nrhs, panel_size, relax;
  double *bd;
  SuperMatrix *B;

  /* Return immediately if cv_mem is NULL. */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSLS", "cvSuperLUMT", 
		    MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the Direct solver */
  if (cv_mem->cv_tempv->ops->nvgetarraypointer == NULL) {
    cvProcessError(cv_mem, CVSLS_ILL_INPUT, "CVSLS", "cvSuperLUMT", 
		    MSGSP_BAD_NVECTOR);
    return(CVSLS_ILL_INPUT);
  }

  if (cv_mem->cv_lfree != NULL) cv_mem->cv_lfree(cv_mem);

  /* Set five main function fields in cv_mem. */
  cv_mem->cv_linit  = cvSuperLUMTInit;
  cv_mem->cv_lsetup = cvSuperLUMTSetup;
  cv_mem->cv_lsolve = cvSuperLUMTSolve;
  cv_mem->cv_lfree  = cvSuperLUMTFree;

  /* Get memory for CVSlsMemRec. */
  cvsls_mem = (CVSlsMem) malloc(sizeof(struct CVSlsMemRec));
  if (cvsls_mem == NULL) {
    cvProcessError(cv_mem, CVSLS_MEM_FAIL, "CVSLS", "cvSuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(CVSLS_MEM_FAIL);
  }

  /* Get memory for SLUMT_data. */
  slumt_data = (SLUMTData)malloc(sizeof(struct SLUMTDataRec));
  if (slumt_data == NULL) {
    cvProcessError(cv_mem, CVSLS_MEM_FAIL, "CVSLS", "cvSuperLUMT", 
		    MSGSP_MEM_FAIL);
    return(CVSLS_MEM_FAIL);
  }

  cv_mem->cv_setupNonNull = TRUE;

  /* Set default Jacobian routine and Jacobian data */
  cvsls_mem->s_jaceval = NULL;
  cvsls_mem->s_jacdata = cv_mem->cv_user_data;

  /* Allocate memory for the sparse Jacobian */
  cvsls_mem->s_JacMat = NewSparseMat(n, n, nnz);
  if (cvsls_mem->s_JacMat == NULL) {
    cvProcessError(cv_mem, CVSLS_MEM_FAIL, "CVSLS", "cvSuperLUMT", 
		    MSGSP_MEM_FAIL);
    free(cvsls_mem);
    return(CVSLS_MEM_FAIL);
  }

  /* Allocate memory for saved sparse Jacobian */
  cvsls_mem->s_savedJ = NewSparseMat(n, n, nnz);
  if (cvsls_mem->s_savedJ == NULL) {
    cvProcessError(cv_mem, CVSLS_MEM_FAIL, "CVSLS", "cvSuperLUMT", 
		    MSGSP_MEM_FAIL);
    DestroySparseMat(cvsls_mem->s_JacMat);
    free(cvsls_mem);
    return(CVSLS_MEM_FAIL);
  }

  /* Set up memory for the permutations */
  perm_r = (int *)malloc(n*sizeof(int));
  if (perm_r == NULL) {
    cvProcessError(cv_mem, CVSLS_MEM_FAIL, "CVSLS", "cvSuperLUMT", 
		   MSGSP_MEM_FAIL);
    return(CVSLS_MEM_FAIL);
  }
  perm_c = (int *)malloc(n*sizeof(int));
  if (perm_c == NULL) {
    cvProcessError(cv_mem, CVSLS_MEM_FAIL, "CVSLS", "cvSuperLUMT", 
		   MSGSP_MEM_FAIL);
    free(perm_r);
    return(CVSLS_MEM_FAIL);
  }
  slumt_data->perm_r = perm_r;
  slumt_data->perm_c = perm_c;

  /* Set default parameters for SuperLU */
  slumt_data->num_threads = num_threads;
  slumt_data->diag_pivot_thresh = 1.0;

  /* Allocate structures for SuperLU */
  slumt_data->Gstat = (Gstat_t *)malloc(sizeof(Gstat_t));
  slumt_data->s_A = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  slumt_data->s_AC = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  slumt_data->s_L = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  slumt_data->s_U = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  slumt_data->s_A->Store  = NULL;
  slumt_data->s_AC->Store = NULL;
  slumt_data->s_L->Store  = NULL;
  slumt_data->s_U->Store  = NULL;
  slumt_data->superlumt_options = (superlumt_options_t *)malloc(sizeof(superlumt_options_t));

  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  StatAlloc(cvsls_mem->s_JacMat->N, num_threads, panel_size, relax, 
	    slumt_data->Gstat);
  
  /* Create RHS matrix */
  nrhs = 1;
  bd = NULL;
  B = (SuperMatrix *)malloc(sizeof(SuperMatrix));
  B->Store = NULL;
  dCreate_Dense_Matrix(B, n, nrhs, bd, n, 
		       SLU_DN, SLU_D, SLU_GE);
  slumt_data->s_B = B;

  /* Set ordering to COLAMD as the cvode default use.
     Users can set a different value with CVSuperLUMTSetOrdering,
     and the user-set value is loaded before any call to factorize the
     matrix in cvSuperLUMTSetup.  */
  slumt_data->s_ordering = 3;

  /* Attach linear solver memory to the integrator memory */
  cvsls_mem->s_solver_data = (void *) slumt_data;
  cv_mem->cv_lmem = cvsls_mem;

  cvsls_mem->s_last_flag = CVSLS_SUCCESS;

  return(CVSLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVSuperLUMT interface functions
 * -----------------------------------------------------------------
 */

/*
  This routine does remaining initializations specific to the CVSuperLUMT
  linear solver module.  
  It returns 0 if successful.
*/

static int cvSuperLUMTInit(CVodeMem cv_mem)
{
  int num_threads, n;
  CVSlsMem cvsls_mem;
  SLUMTData slumt_data;

  cvsls_mem = (CVSlsMem)cv_mem->cv_lmem;
  slumt_data = (SLUMTData) cvsls_mem->s_solver_data;

  cvsls_mem->s_nje = 0;
  cvsls_mem->s_first_factorize = 1;
  cvsls_mem->s_nstlj = 0;

  /* ------------------------------------------------------------
     Allocate storage and initialize statistics variables. 
     ------------------------------------------------------------*/
  n = cvsls_mem->s_JacMat->N;
  num_threads = slumt_data->num_threads;

  StatInit(n, num_threads, slumt_data->Gstat);

  cvsls_mem->s_last_flag = 0;
  return(0);
}

/*
  This routine does the setup operations for the CVSuperLUMT linear 
  solver module.  It calls the Jacobian evaluation routine,
  updates counters, and calls the LU factorization routine.
  The return value is either
     CVSLS_SUCCESS = 0  if successful,
     +1  if the jac routine failed recoverably or the
         LU factorization failed, or
     -1  if the jac routine failed unrecoverably.
*/

static int cvSuperLUMTSetup(CVodeMem cv_mem, int convfail, N_Vector ypred, 
			   N_Vector fpred, booleantype *jcurPtr,
			   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  booleantype jbad, jok;
  int retval, info;
  int nprocs, panel_size, relax, permc_spec, lwork;
  int *perm_r, *perm_c;
  long int nst, nstlj;
  realtype tn, gamma, gammap, dgamma;
  double diag_pivot_thresh, drop_tol;
  fact_t fact;
  trans_t trans;
  yes_no_t refact, usepr;
  CVSlsMem cvsls_mem;
  CVSlsSparseJacFn jaceval;
  SuperMatrix *A, *AC, *L, *U;
  Gstat_t *Gstat;
  superlumt_options_t *superlumt_options;
  SLUMTData slumt_data;
  SlsMat JacMat, savedJ;
  void *jacdata;
  void *work;
  
  cvsls_mem = (CVSlsMem) (cv_mem->cv_lmem);
  tn = cv_mem->cv_tn; 
  gamma = cv_mem->cv_gamma;
  gammap = cv_mem->cv_gammap;
  nst = cv_mem->cv_nst;

  slumt_data = (SLUMTData) cvsls_mem->s_solver_data;

  jaceval = cvsls_mem->s_jaceval;
  jacdata = cvsls_mem->s_jacdata;
  JacMat = cvsls_mem->s_JacMat;
  savedJ = cvsls_mem->s_savedJ;
  nstlj = cvsls_mem->s_nstlj;

  superlumt_options = slumt_data->superlumt_options;
  A = slumt_data->s_A;
  AC = slumt_data->s_AC;
  L = slumt_data->s_L;
  U = slumt_data->s_U;
  Gstat = slumt_data->Gstat;
  perm_r = slumt_data->perm_r;
  perm_c = slumt_data->perm_c;
  nprocs = slumt_data->num_threads;
  diag_pivot_thresh = slumt_data->diag_pivot_thresh;

  /* Set option values for SuperLU_MT */
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  fact = EQUILIBRATE;
  trans = NOTRANS;
  usepr = NO;
  drop_tol = 0.0;
  lwork = 0;
  work = NULL;

  /* Check that Jacobian eval routine is set */
  if (jaceval == NULL) {
    cvProcessError(cv_mem, CVSLS_JAC_NOSET, "CVSLS", "cvSuperLUMTSetup", 
		    MSGSP_JAC_NOSET);
    free(cvsls_mem); cvsls_mem = NULL;
    return(CVSLS_JAC_NOSET);
  }

  /* Determine whether Jacobian needs to be recalculated */
  dgamma = SUNRabs((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + CVS_MSBJ) ||
         ((convfail == CV_FAIL_BAD_J) && (dgamma < CVS_DGMAX)) ||
         (convfail == CV_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    CopySparseMat(savedJ, JacMat);
  } else {
    /* If jok = FALSE, call jac routine for new J value */
    cvsls_mem->s_nje++;
    cvsls_mem->s_nstlj = nst;
    *jcurPtr = TRUE;
    SlsSetToZero(JacMat);
    retval = jaceval(tn, ypred, fpred, JacMat, jacdata, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      cvProcessError(cv_mem, CVSLS_JACFUNC_UNRECVR, "CVSLS", "cvSuperLUMTSetup", MSGSP_JACFUNC_FAILED);
      cvsls_mem->s_last_flag = CVSLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      cvsls_mem->s_last_flag = CVSLS_JACFUNC_RECVR;
      return(1);
    }

    CopySparseMat(JacMat, savedJ);
  }

  /* Scale and add I to get M = I - gamma*J */
  ScaleSparseMat(-gamma, JacMat);
  AddIdentitySparseMat(JacMat);

  if (A->Store) {
    SUPERLU_FREE(A->Store);
  }
  dCreate_CompCol_Matrix(A, JacMat->M, JacMat->N, JacMat->NNZ, 
			 JacMat->data, JacMat->rowvals, JacMat->colptrs, 
			 SLU_NC, SLU_D, SLU_GE);

  if (cvsls_mem->s_first_factorize) {
    /* ------------------------------------------------------------
       Get column permutation vector perm_c[], according to permc_spec:
       permc_spec = 3: approximate minimum degree for unsymmetric matrices
       ------------------------------------------------------------*/ 
    permc_spec = slumt_data->s_ordering;
    get_perm_c(permc_spec, A, perm_c);

    refact= NO;
    cvsls_mem->s_first_factorize = 0;
  }
  else {
    /* ------------------------------------------------------------
       Re-initialize statistics variables 
       ------------------------------------------------------------*/
    StatInit(JacMat->N, nprocs, Gstat);
    Destroy_CompCol_Permuted(AC);
    refact= YES;
  }

  /* ------------------------------------------------------------
     Initialize the option structure superlumt_options using the
     user-input parameters;  Subsequent calls will re-initialize
     options.
     Apply perm_c to the columns of original A to form AC.
     ------------------------------------------------------------*/
  pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
	       diag_pivot_thresh, usepr, drop_tol, perm_c, perm_r,
	       work, lwork, A, AC, superlumt_options, Gstat);
  /* ------------------------------------------------------------
     Compute the LU factorization of A.
     The following routine will create nprocs threads.
     ------------------------------------------------------------*/
  pdgstrf(superlumt_options, AC, perm_r, L, U, Gstat, &info);
    
  if (info != 0) {
    cvsls_mem->s_last_flag = info;
    return(+1);
  }

  cvsls_mem->s_last_flag = CVSLS_SUCCESS;

  return(0);
}

/*
  This routine handles the solve operation for the CVSuperLUMT linear
  solver module.  It calls the SuperLU_MT solve routine,
  then returns CVSLS_SUCCESS = 0.
*/

static int cvSuperLUMTSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
			    N_Vector ycur, N_Vector fcur)
{
  int info, trans, lmm;
  int *perm_r, *perm_c;
  realtype gamrat;
  CVSlsMem cvsls_mem;
  SuperMatrix *L, *U, *B;
  Gstat_t *Gstat;
  DNformat *Bstore;
  SLUMTData slumt_data;
  realtype *bd;
  
  gamrat = cv_mem->cv_gamrat;
  lmm = cv_mem->cv_lmm;

  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;

  slumt_data = (SLUMTData) cvsls_mem->s_solver_data;

  L = slumt_data->s_L;
  U = slumt_data->s_U;
  perm_r = slumt_data->perm_r;
  perm_c = slumt_data->perm_c;
  Gstat = slumt_data->Gstat;
  B = slumt_data->s_B;
   
  bd = N_VGetArrayPointer(b);
  Bstore = (DNformat *) (B->Store);
  Bstore->nzval = bd;

  /* Call SuperLUMT to solve the linear system using L and U */
  trans = NOTRANS;
  dgstrs(trans, L, U, perm_r, perm_c, B, Gstat, &info);

  /* Scale the correction to account for change in gamma. */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }

  Bstore->nzval = NULL;

  cvsls_mem->s_last_flag = CVSLS_SUCCESS;
  return(CVSLS_SUCCESS);
}

/*
  This routine frees memory specific to the CVSuperLUMT linear solver.
*/

static void cvSuperLUMTFree(CVodeMem cv_mem)
{
  CVSlsMem cvsls_mem;
  SLUMTData slumt_data;
  
  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;

  slumt_data = (SLUMTData) cvsls_mem->s_solver_data;

  pxgstrf_finalize(slumt_data->superlumt_options, slumt_data->s_AC);

  free(slumt_data->perm_r);
  free(slumt_data->perm_c);
  free(slumt_data->superlumt_options);
  Destroy_SuperNode_SCP( (slumt_data->s_L) );
  Destroy_CompCol_NCP( (slumt_data->s_U) );
  StatFree( (slumt_data->Gstat) );
  free(slumt_data->Gstat);
  
  Destroy_SuperMatrix_Store(slumt_data->s_B);
  SUPERLU_FREE(slumt_data->s_A->Store);
  if (cvsls_mem->s_JacMat) {
    DestroySparseMat(cvsls_mem->s_JacMat);
    cvsls_mem->s_JacMat = NULL;
  }
  if (cvsls_mem->s_savedJ) {
    DestroySparseMat(cvsls_mem->s_savedJ);
    cvsls_mem->s_savedJ = NULL;
  }

  free(slumt_data->s_B);
  free(slumt_data->s_A);
  free(slumt_data->s_AC);
  free(slumt_data->s_L);
  free(slumt_data->s_U);

  free(slumt_data); 
  slumt_data = NULL;
  free(cv_mem->cv_lmem); 

  return;
}

/* 
 * -----------------------------------------------------------------
 * Optional Input Specification Functions
 * -----------------------------------------------------------------
 *
 * CVSuperLUMTSetOrdering sets the ordering used by SuperLUMT for reducing fill.
 * Options are: 
 * 0 for natural ordering
 * 1 for minimal degree ordering on A'*A
 * 2 for minimal degree ordering on A'+A
 * 3 for approximate minimal degree ordering for unsymmetric matrices
 * The default used in SUNDIALS is 3 for COLAMD.
 * -----------------------------------------------------------------
 */

int CVSuperLUMTSetOrdering(void *cv_mem_v, int ordering_choice)
{
  CVodeMem cv_mem;
  CVSlsMem cvsls_mem;
  SLUMTData slumt_data;

 /* Return immediately if cv_mem is NULL */
  if (cv_mem_v == NULL) {
    cvProcessError(NULL, CVSLS_MEM_NULL, "CVSLS", "CVSuperLUMTSetOrdering",
		    MSGSP_CVMEM_NULL);
    return(CVSLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cv_mem_v;

 /* Return if ordering choice argument is not valid */
  if ( (ordering_choice != 0) && (ordering_choice != 1) && 
       (ordering_choice != 2) && (ordering_choice != 3) ) {
    cvProcessError(NULL, CVSLS_ILL_INPUT, "CVSLS", "CVSuperLUMTSetOrdering",
		    MSGSP_ILL_INPUT);
    return(CVSLS_ILL_INPUT);
  }

  cvsls_mem = (CVSlsMem) cv_mem->cv_lmem;
  slumt_data = (SLUMTData) cvsls_mem->s_solver_data;

  slumt_data->s_ordering = ordering_choice;

  return(CVSLS_SUCCESS);
}
