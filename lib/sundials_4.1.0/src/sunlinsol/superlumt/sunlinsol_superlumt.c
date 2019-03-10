/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 * Based on codes <solver>_superlumt.c, written by 
 *     Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the SuperLUMT implementation of 
 * the SUNLINSOL package.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sunlinsol/sunlinsol_superlumt.h>
#include <sundials/sundials_math.h>

#define ZERO      RCONST(0.0)
#define ONE       RCONST(1.0)
#define TWO       RCONST(2.0)

/* Private function prototypes */
sunindextype GlobalVectorLength_SuperLUMT(N_Vector y);

/*
 * -----------------------------------------------------------------
 * SuperLUMT solver structure accessibility macros: 
 * -----------------------------------------------------------------
 */

#define SLUMT_CONTENT(S)    ( (SUNLinearSolverContent_SuperLUMT)(S->content) )
#define LASTFLAG(S)         ( SLUMT_CONTENT(S)->last_flag )
#define FIRSTFACTORIZE(S)   ( SLUMT_CONTENT(S)->first_factorize )
#define SM_A(S)             ( SLUMT_CONTENT(S)->A )
#define SM_AC(S)            ( SLUMT_CONTENT(S)->AC )
#define SM_L(S)             ( SLUMT_CONTENT(S)->L )
#define SM_U(S)             ( SLUMT_CONTENT(S)->U )
#define SM_B(S)             ( SLUMT_CONTENT(S)->B )
#define GSTAT(S)            ( SLUMT_CONTENT(S)->Gstat )
#define PERMR(S)            ( SLUMT_CONTENT(S)->perm_r )
#define PERMC(S)            ( SLUMT_CONTENT(S)->perm_c )
#define SIZE(S)             ( SLUMT_CONTENT(S)->N )
#define NUMTHREADS(S)       ( SLUMT_CONTENT(S)->num_threads )
#define DIAGPIVOTTHRESH(S)  ( SLUMT_CONTENT(S)->diag_pivot_thresh )
#define ORDERING(S)         ( SLUMT_CONTENT(S)->ordering )
#define OPTIONS(S)          ( SLUMT_CONTENT(S)->options )

/*
 * -----------------------------------------------------------------
 * deprecated wrapper functions
 * -----------------------------------------------------------------
 */

SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A, int num_threads)
{ return(SUNLinSol_SuperLUMT(y, A, num_threads)); }

int SUNSuperLUMTSetOrdering(SUNLinearSolver S, int ordering_choice)
{ return(SUNLinSol_SuperLUMTSetOrdering(S, ordering_choice)); }


/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new SuperLUMT linear solver
 */

SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y, SUNMatrix A, int num_threads)
{
  SUNLinearSolver S;
  SUNLinearSolver_Ops ops;
  SUNLinearSolverContent_SuperLUMT content;
  sunindextype MatrixRows, VecLength;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE)
    return(NULL);
  if (SUNSparseMatrix_Rows(A) != SUNSparseMatrix_Columns(A))
    return(NULL);
  MatrixRows = SUNSparseMatrix_Rows(A);
  if ( (N_VGetVectorID(y) != SUNDIALS_NVEC_SERIAL) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_OPENMP) &&
       (N_VGetVectorID(y) != SUNDIALS_NVEC_PTHREADS) )
    return(NULL);

  /* optimally this function would be replaced with a generic N_Vector routine */
  VecLength = GlobalVectorLength_SuperLUMT(y);
  if (MatrixRows != VecLength)
    return(NULL);
  
  /* Create linear solver */
  S = NULL;
  S = (SUNLinearSolver) malloc(sizeof *S);
  if (S == NULL) return(NULL);
  
  /* Create linear solver operation structure */
  ops = NULL;
  ops = (SUNLinearSolver_Ops) malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  if (ops == NULL) { free(S); return(NULL); }

  /* Attach operations */
  ops->gettype           = SUNLinSolGetType_SuperLUMT;
  ops->initialize        = SUNLinSolInitialize_SuperLUMT;
  ops->setup             = SUNLinSolSetup_SuperLUMT;
  ops->solve             = SUNLinSolSolve_SuperLUMT;
  ops->lastflag          = SUNLinSolLastFlag_SuperLUMT;
  ops->space             = SUNLinSolSpace_SuperLUMT;
  ops->free              = SUNLinSolFree_SuperLUMT;
  ops->setatimes         = NULL;
  ops->setpreconditioner = NULL;
  ops->setscalingvectors = NULL;
  ops->numiters          = NULL;
  ops->resnorm           = NULL;
  ops->resid             = NULL;

  /* Create content */
  content = NULL;
  content = (SUNLinearSolverContent_SuperLUMT)
    malloc(sizeof(struct _SUNLinearSolverContent_SuperLUMT));
  if (content == NULL) { free(ops); free(S); return(NULL); }

  /* Fill content */
  content->N = MatrixRows;
  content->last_flag = 0;
  content->num_threads = num_threads;
  content->diag_pivot_thresh = ONE;
  content->ordering = SUNSLUMT_ORDERING_DEFAULT;

  content->perm_r = NULL;
  content->perm_r = (sunindextype *) malloc(MatrixRows*sizeof(sunindextype));
  if (content->perm_r == NULL) {
    free(content); free(ops); free(S); return(NULL); }

  content->perm_c = NULL;
  content->perm_c = (sunindextype *) malloc(MatrixRows*sizeof(sunindextype));
  if (content->perm_c == NULL) {
    free(content->perm_r); free(content); free(ops); free(S); return(NULL); }

  content->Gstat = (Gstat_t *) malloc(sizeof(Gstat_t));
  if (content->Gstat == NULL) {
    free(content->perm_c); free(content->perm_r); free(content); free(ops);
    free(S); return(NULL); }

  content->A = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->A == NULL) {
    free(content->Gstat); free(content->perm_c); free(content->perm_r);
    free(content); free(ops); free(S); return(NULL); }
  content->A->Store = NULL;

  content->AC = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->AC == NULL) {
    free(content->A); free(content->Gstat); free(content->perm_c);
    free(content->perm_r); free(content); free(ops); free(S); return(NULL); }
  content->AC->Store = NULL;

  content->L = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->L == NULL) {
    free(content->AC); free(content->A); free(content->Gstat); free(content->perm_c);
    free(content->perm_r); free(content); free(ops); free(S); return(NULL); }
  content->L->Store = NULL;

  content->U = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->U == NULL) {
    free(content->L); free(content->AC); free(content->A); free(content->Gstat);
    free(content->perm_c); free(content->perm_r); free(content); free(ops); free(S);
    return(NULL); }
  content->U->Store = NULL;

  content->B = (SuperMatrix *) malloc(sizeof(SuperMatrix));
  if (content->B == NULL) {
    free(content->U); free(content->L); free(content->AC); free(content->A);
    free(content->Gstat); free(content->perm_c); free(content->perm_r); free(content);
    free(ops); free(S); return(NULL); }
  content->B->Store = NULL;
  xCreate_Dense_Matrix(content->B, MatrixRows, 1, NULL, MatrixRows, SLU_DN, SLU_D, SLU_GE);
  
  content->options = (superlumt_options_t *) malloc(sizeof(superlumt_options_t));
  if (content->options == NULL) {
    free(content->B); free(content->U); free(content->L); free(content->AC);
    free(content->A); free(content->Gstat); free(content->perm_c); free(content->perm_r);
    free(content); free(ops); free(S); return(NULL); }
  StatAlloc(MatrixRows, num_threads, sp_ienv(1), sp_ienv(2), content->Gstat);

  /* Attach content and ops */
  S->content = content;
  S->ops     = ops;

  return(S);
}


/* ----------------------------------------------------------------------------
 * Function to set the ordering type for a SuperLUMT linear solver
 */

int SUNLinSol_SuperLUMTSetOrdering(SUNLinearSolver S, int ordering_choice)
{
  /* Check for legal ordering_choice */ 
  if ((ordering_choice < 0) || (ordering_choice > 3))
    return(SUNLS_ILL_INPUT);

  /* Check for non-NULL SUNLinearSolver */
  if (S == NULL) return(SUNLS_MEM_NULL);

  /* Set ordering_choice */
  ORDERING(S) = ordering_choice;

  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_SuperLUMT(SUNLinearSolver S)
{
  return(SUNLINEARSOLVER_DIRECT);
}


int SUNLinSolInitialize_SuperLUMT(SUNLinearSolver S)
{
  /* force a first factorization */
  FIRSTFACTORIZE(S) = 1;

  /* Initialize statistics variables */
  StatInit(SIZE(S), NUMTHREADS(S), GSTAT(S));
  
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSetup_SuperLUMT(SUNLinearSolver S, SUNMatrix A)
{
  int_t retval;
  int panel_size, relax, lwork;
  double drop_tol;
  fact_t fact;
  trans_t trans;
  yes_no_t refact, usepr;
  void *work;
  
  /* Set option values for SuperLU_MT */
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  fact = EQUILIBRATE;
  trans = (SUNSparseMatrix_SparseType(A) == CSC_MAT) ? NOTRANS : TRANS;
  usepr = NO;
  drop_tol = ZERO;
  lwork = 0;
  work = NULL;

  /* free and reallocate sparse matrix */
  if (SM_A(S)->Store)
    SUPERLU_FREE(SM_A(S)->Store);
  xCreate_CompCol_Matrix(SM_A(S), SUNSparseMatrix_Rows(A), 
			 SUNSparseMatrix_Columns(A),
                         SUNSparseMatrix_NNZ(A), 
                         SUNSparseMatrix_Data(A), 
                         (int_t*) SUNSparseMatrix_IndexValues(A), 
                         (int_t*) SUNSparseMatrix_IndexPointers(A), 
			 SLU_NC, SLU_D, SLU_GE);

  /* On first decomposition, set up reusable pieces */ 
  if (FIRSTFACTORIZE(S)) {

    /* Get column permutation vector perm_c[], according to ordering */
    get_perm_c(ORDERING(S), SM_A(S), (int_t *) PERMC(S));
    refact = NO;
    FIRSTFACTORIZE(S) = 0;

  } else {

    /* Re-initialize statistics variables */
    StatInit(SIZE(S), NUMTHREADS(S), GSTAT(S)); 
    Destroy_CompCol_Permuted(SM_AC(S));
    refact = YES;
    
  }

  /* Initialize the option structure using the user-input parameters. 
     Subsequent calls will re-initialize options.  Apply perm_c to 
     columns of original A to form AC */
  pxgstrf_init(NUMTHREADS(S), fact, trans, refact, panel_size, relax, 
	       DIAGPIVOTTHRESH(S), usepr, drop_tol, (int_t *) PERMC(S), (int_t *) PERMR(S), 
               work, lwork, SM_A(S), SM_AC(S), OPTIONS(S), GSTAT(S));

  /* Compute the LU factorization of A.
     The following routine will create num_threads threads. */
  pxgstrf(OPTIONS(S), SM_AC(S), (int_t *) PERMR(S), SM_L(S), SM_U(S),
          GSTAT(S), &retval);
  if (retval != 0) {
    LASTFLAG(S) = (retval < 0) ? 
      SUNLS_PACKAGE_FAIL_UNREC : SUNLS_PACKAGE_FAIL_REC;
    return(LASTFLAG(S));
  }
  
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


int SUNLinSolSolve_SuperLUMT(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                             N_Vector b, realtype tol)
{
  int_t retval;
  realtype *xdata;
  DNformat *Bstore;
  trans_t trans;
  
  /* copy b into x */
  N_VScale(ONE, b, x);

  /* access x data array */
  xdata = N_VGetArrayPointer(x);
  if (xdata == NULL) {
    LASTFLAG(S) = SUNLS_MEM_FAIL;
    return(LASTFLAG(S));
  }

  Bstore = (DNformat *) (SM_B(S)->Store);
  Bstore->nzval = xdata;
  
  /* Call SuperLUMT to solve the linear system using L and U */
  trans = (SUNSparseMatrix_SparseType(A) == CSC_MAT) ? NOTRANS : TRANS;
  xgstrs(trans, SM_L(S), SM_U(S), (int_t *) PERMR(S), (int_t *) PERMC(S), SM_B(S), GSTAT(S), &retval);
  if (retval != 0) {
    LASTFLAG(S) = SUNLS_PACKAGE_FAIL_UNREC;
    return(LASTFLAG(S));
  }
  
  LASTFLAG(S) = SUNLS_SUCCESS;
  return(LASTFLAG(S));
}


long int SUNLinSolLastFlag_SuperLUMT(SUNLinearSolver S)
{
  /* return the stored 'last_flag' value */
  return(LASTFLAG(S));
}


int SUNLinSolSpace_SuperLUMT(SUNLinearSolver S, 
                             long int *lenrwLS, 
                             long int *leniwLS)
{
  /* since the SuperLU_MT structures are opaque objects, we 
     omit those from these results */
  *leniwLS = 5 + 2*SIZE(S);
  *lenrwLS = 1;
  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_SuperLUMT(SUNLinearSolver S)
{
  /* return with success if already freed */
  if (S == NULL)
    return(SUNLS_SUCCESS);
  
  /* delete items from the contents structure (if it exists) */
  if (S->content) {
    pxgstrf_finalize(OPTIONS(S), SM_AC(S));
    free(PERMR(S));
    free(PERMC(S));
    free(OPTIONS(S));
    Destroy_SuperNode_SCP(SM_L(S));
    Destroy_CompCol_NCP(SM_U(S));
    StatFree(GSTAT(S));
    free(GSTAT(S));
  
    Destroy_SuperMatrix_Store(SM_B(S));
    SUPERLU_FREE(SM_A(S)->Store);

    free(SM_B(S));
    free(SM_A(S));
    free(SM_AC(S));
    free(SM_L(S));
    free(SM_U(S));

    free(S->content);  
    S->content = NULL;
  }
  
  /* delete generic structures */
  if (S->ops) {
    free(S->ops);  
    S->ops = NULL;
  }
  free(S); S = NULL;
  return(SUNLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

/* Inefficient kludge for determining the number of entries in a N_Vector 
   object (replace if such a routine is ever added to the N_Vector API).

   Returns "-1" on an error. */
sunindextype GlobalVectorLength_SuperLUMT(N_Vector y)
{
  realtype len;
  N_Vector tmp = NULL;
  tmp = N_VClone(y);
  if (tmp == NULL)  return(-1);
  N_VConst(ONE, tmp);
  len = N_VDotProd(tmp, tmp);
  N_VDestroy(tmp);
  return( (sunindextype) len );
}
