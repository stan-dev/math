/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
 * Based on code sundials_dense.c by: Scott D. Cohen,
 *     Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for the dense implementation of
 * the SUNMATRIX package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)


/* Private function prototypes */
static booleantype SMCompatible_Dense(SUNMatrix A, SUNMatrix B);
static booleantype SMCompatible2_Dense(SUNMatrix A, N_Vector x, N_Vector y);


/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new dense matrix
 */

SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N)
{
  SUNMatrix A;
  SUNMatrixContent_Dense content;
  sunindextype j;

  /* return with NULL matrix on illegal dimension input */
  if ( (M <= 0) || (N <= 0) ) return(NULL);

  /* Create an empty matrix object */
  A = NULL;
  A = SUNMatNewEmpty();
  if (A == NULL) return(NULL);

  /* Attach operations */
  A->ops->getid     = SUNMatGetID_Dense;
  A->ops->clone     = SUNMatClone_Dense;
  A->ops->destroy   = SUNMatDestroy_Dense;
  A->ops->zero      = SUNMatZero_Dense;
  A->ops->copy      = SUNMatCopy_Dense;
  A->ops->scaleadd  = SUNMatScaleAdd_Dense;
  A->ops->scaleaddi = SUNMatScaleAddI_Dense;
  A->ops->matvec    = SUNMatMatvec_Dense;
  A->ops->space     = SUNMatSpace_Dense;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_Dense) malloc(sizeof *content);
  if (content == NULL) { SUNMatDestroy(A); return(NULL); }

  /* Attach content */
  A->content = content;

  /* Fill content */
  content->M     = M;
  content->N     = N;
  content->ldata = M*N;
  content->data  = NULL;
  content->cols  = NULL;

  /* Allocate content */
  content->data = (realtype *) calloc(M * N, sizeof(realtype));
  if (content->data == NULL) { SUNMatDestroy(A); return(NULL); }

  content->cols = (realtype **) malloc(N * sizeof(realtype *));
  if (content->cols == NULL) { SUNMatDestroy(A); return(NULL); }
  for (j=0; j<N; j++) content->cols[j] = content->data + j * M;

  return(A);
}


/* ----------------------------------------------------------------------------
 * Function to print the dense matrix 
 */
 
void SUNDenseMatrix_Print(SUNMatrix A, FILE* outfile)
{
  sunindextype i, j;
  
  /* should not be called unless A is a dense matrix; 
     otherwise return immediately */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
    return;

  /* perform operation */
  STAN_SUNDIALS_FPRINTF(outfile,"\n");
  for (i=0; i<SM_ROWS_D(A); i++) {
    for (j=0; j<SM_COLUMNS_D(A); j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      STAN_SUNDIALS_FPRINTF(outfile,"%12Lg  ", SM_ELEMENT_D(A,i,j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      STAN_SUNDIALS_FPRINTF(outfile,"%12g  ", SM_ELEMENT_D(A,i,j));
#else
      STAN_SUNDIALS_FPRINTF(outfile,"%12g  ", SM_ELEMENT_D(A,i,j));
#endif
    }
    STAN_SUNDIALS_FPRINTF(outfile,"\n");
  }
  STAN_SUNDIALS_FPRINTF(outfile,"\n");
  return;
}


/* ----------------------------------------------------------------------------
 * Functions to access the contents of the dense matrix structure
 */

sunindextype SUNDenseMatrix_Rows(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_DENSE)
    return SM_ROWS_D(A);
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNDenseMatrix_Columns(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_DENSE)
    return SM_COLUMNS_D(A);
  else
    return SUNMAT_ILL_INPUT;
}

sunindextype SUNDenseMatrix_LData(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_DENSE)
    return SM_LDATA_D(A);
  else
    return SUNMAT_ILL_INPUT;
}

realtype* SUNDenseMatrix_Data(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_DENSE)
    return SM_DATA_D(A);
  else
    return NULL;
}

realtype** SUNDenseMatrix_Cols(SUNMatrix A)
{
  if (SUNMatGetID(A) == SUNMATRIX_DENSE)
    return SM_COLS_D(A);
  else
    return NULL;
}

realtype* SUNDenseMatrix_Column(SUNMatrix A, sunindextype j)
{
  if (SUNMatGetID(A) == SUNMATRIX_DENSE)
    return SM_COLUMN_D(A,j);
  else
    return NULL;
}


/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_Dense(SUNMatrix A)
{
  return SUNMATRIX_DENSE;
}

SUNMatrix SUNMatClone_Dense(SUNMatrix A)
{
  SUNMatrix B = SUNDenseMatrix(SM_ROWS_D(A), SM_COLUMNS_D(A));
  return(B);
}

void SUNMatDestroy_Dense(SUNMatrix A)
{
  if (A == NULL) return;

  /* free content */
  if (A->content != NULL) {
    /* free data array */
    if (SM_DATA_D(A) != NULL) {
      free(SM_DATA_D(A));
      SM_DATA_D(A) = NULL;
    }
    /* free column pointers */
    if (SM_CONTENT_D(A)->cols != NULL) {
      free(SM_CONTENT_D(A)->cols);
      SM_CONTENT_D(A)->cols = NULL;
    }
    /* free content struct */
    free(A->content);
    A->content = NULL;
  }

  /* free ops and matrix */
  if (A->ops) { free(A->ops); A->ops = NULL; }
  free(A); A = NULL;

  return;
}

int SUNMatZero_Dense(SUNMatrix A)
{
  sunindextype i;
  realtype *Adata;

  /* Perform operation */
  Adata = SM_DATA_D(A);
  for (i=0; i<SM_LDATA_D(A); i++)
    Adata[i] = ZERO;
  return SUNMAT_SUCCESS;
}

int SUNMatCopy_Dense(SUNMatrix A, SUNMatrix B)
{
  sunindextype i, j;

  /* Verify that A and B are compatible */
  if (!SMCompatible_Dense(A, B))
    return SUNMAT_ILL_INPUT;

  /* Perform operation */
  for (j=0; j<SM_COLUMNS_D(A); j++)
    for (i=0; i<SM_ROWS_D(A); i++)
      SM_ELEMENT_D(B,i,j) = SM_ELEMENT_D(A,i,j);
  return SUNMAT_SUCCESS;
}

int SUNMatScaleAddI_Dense(realtype c, SUNMatrix A)
{
  sunindextype i, j;

  /* Perform operation */
  for (j=0; j<SM_COLUMNS_D(A); j++)
    for (i=0; i<SM_ROWS_D(A); i++) {
      SM_ELEMENT_D(A,i,j) *= c;
      if (i == j) 
        SM_ELEMENT_D(A,i,j) += ONE;
    }
  return SUNMAT_SUCCESS;
}

int SUNMatScaleAdd_Dense(realtype c, SUNMatrix A, SUNMatrix B)
{
  sunindextype i, j;

  /* Verify that A and B are compatible */
  if (!SMCompatible_Dense(A, B))
    return SUNMAT_ILL_INPUT;

  /* Perform operation */
  for (j=0; j<SM_COLUMNS_D(A); j++)
    for (i=0; i<SM_ROWS_D(A); i++)
      SM_ELEMENT_D(A,i,j) = c*SM_ELEMENT_D(A,i,j) + SM_ELEMENT_D(B,i,j);
  return SUNMAT_SUCCESS;
}

int SUNMatMatvec_Dense(SUNMatrix A, N_Vector x, N_Vector y)
{
  sunindextype i, j;
  realtype *col_j, *xd, *yd;
  
  /* Verify that A, x and y are compatible */
  if (!SMCompatible2_Dense(A, x, y))
    return SUNMAT_ILL_INPUT;

  /* access vector data (return if failure) */
  xd = N_VGetArrayPointer(x);
  yd = N_VGetArrayPointer(y);
  if ((xd == NULL) || (yd == NULL) || (xd == yd))
    return SUNMAT_MEM_FAIL;

  /* Perform operation */
  for (i=0; i<SM_ROWS_D(A); i++)
    yd[i] = ZERO;
  for(j=0; j<SM_COLUMNS_D(A); j++) {
    col_j = SM_COLUMN_D(A,j);
    for (i=0; i<SM_ROWS_D(A); i++)
      yd[i] += col_j[i]*xd[j];
  }
  return SUNMAT_SUCCESS;
}

int SUNMatSpace_Dense(SUNMatrix A, long int *lenrw, long int *leniw)
{
  *lenrw = SM_LDATA_D(A);
  *leniw = 3 + SM_COLUMNS_D(A);
  return SUNMAT_SUCCESS;
}


/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static booleantype SMCompatible_Dense(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be SUNMATRIX_DENSE */
  if (SUNMatGetID(A) != SUNMATRIX_DENSE)
    return SUNFALSE;
  if (SUNMatGetID(B) != SUNMATRIX_DENSE)
    return SUNFALSE;

  /* both matrices must have the same shape */
  if (SM_ROWS_D(A) != SM_ROWS_D(B))
    return SUNFALSE;
  if (SM_COLUMNS_D(A) != SM_COLUMNS_D(B))
    return SUNFALSE;

  return SUNTRUE;
}


static booleantype SMCompatible2_Dense(SUNMatrix A, N_Vector x, N_Vector y)
{
  /*   vectors must be one of {SERIAL, OPENMP, PTHREADS} */ 
  if ( (N_VGetVectorID(x) != SUNDIALS_NVEC_SERIAL) &&
       (N_VGetVectorID(x) != SUNDIALS_NVEC_OPENMP) &&
       (N_VGetVectorID(x) != SUNDIALS_NVEC_PTHREADS) )
    return SUNFALSE;

  /* Optimally we would verify that the dimensions of A, x and y agree, 
   but since there is no generic 'length' routine for N_Vectors we cannot */

  return SUNTRUE;
}

