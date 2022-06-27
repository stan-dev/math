/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David J. Gardner, Carol S. Woodward, and
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a generic SUNMATRIX package.
 * It contains the implementation of the SUNMatrix operations listed
 * in sundials_matrix.h
 * -----------------------------------------------------------------*/

#include <stdlib.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include "sundials_context_impl.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(SUNMatrix A)
{
  return(A->sunctx->profiler);
}
#endif

/* -----------------------------------------------------------------
 * Create a new empty SUNMatrix object
 * ----------------------------------------------------------------- */

SUNMatrix SUNMatNewEmpty(SUNContext sunctx)
{
  SUNMatrix     A;
  SUNMatrix_Ops ops;

  /* a context is required */
  if (sunctx == NULL) return(NULL);

  /* create matrix object */
  A = NULL;
  A = (SUNMatrix) malloc(sizeof *A);
  if (A == NULL) return(NULL);

  /* create matrix ops structure */
  ops = NULL;
  ops = (SUNMatrix_Ops) malloc(sizeof *ops);
  if (ops == NULL) { free(A); return(NULL); }

  /* initialize operations to NULL */
  ops->getid       = NULL;
  ops->clone       = NULL;
  ops->destroy     = NULL;
  ops->zero        = NULL;
  ops->copy        = NULL;
  ops->scaleadd    = NULL;
  ops->scaleaddi   = NULL;
  ops->matvecsetup = NULL;
  ops->matvec      = NULL;
  ops->space       = NULL;

  /* attach ops and initialize content to NULL */
  A->ops     = ops;
  A->content = NULL;
  A->sunctx  = sunctx;

  return(A);
}


/* -----------------------------------------------------------------
 * Free a generic SUNMatrix (assumes content is already empty)
 * ----------------------------------------------------------------- */

void SUNMatFreeEmpty(SUNMatrix A)
{
  if (A == NULL)  return;

  /* free non-NULL ops structure */
  if (A->ops)  free(A->ops);
  A->ops = NULL;

  /* free overall SUNMatrix object and return */
  free(A);
  return;
}


/* -----------------------------------------------------------------
 * Copy a matrix 'ops' structure
 * -----------------------------------------------------------------*/

int SUNMatCopyOps(SUNMatrix A, SUNMatrix B)
{
  /* Check that ops structures exist */
  if (A == NULL || B == NULL) return(-1);
  if (A->ops == NULL || B->ops == NULL) return(-1);

  /* Copy ops from A to B */
  B->ops->getid       = A->ops->getid;
  B->ops->clone       = A->ops->clone;
  B->ops->destroy     = A->ops->destroy;
  B->ops->zero        = A->ops->zero;
  B->ops->copy        = A->ops->copy;
  B->ops->scaleadd    = A->ops->scaleadd;
  B->ops->scaleaddi   = A->ops->scaleaddi;
  B->ops->matvecsetup = A->ops->matvecsetup;
  B->ops->matvec      = A->ops->matvec;
  B->ops->space       = A->ops->space;

  return(0);
}


/* -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * ----------------------------------------------------------------- */

SUNMatrix_ID SUNMatGetID(SUNMatrix A)
{
  SUNMatrix_ID id;
  id = A->ops->getid(A);
  return(id);
}

SUNMatrix SUNMatClone(SUNMatrix A)
{
  SUNMatrix B = NULL;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(A));
  B = A->ops->clone(A);
  B->sunctx = A->sunctx;
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(A));
  return(B);
}

void SUNMatDestroy(SUNMatrix A)
{
  if (A == NULL) return;

  /* if the destroy operation exists use it */
  if (A->ops)
    if (A->ops->destroy) { A->ops->destroy(A); return; }

  /* if we reach this point, either ops == NULL or destroy == NULL,
     try to cleanup by freeing the content, ops, and matrix */
  if (A->content) { free(A->content); A->content = NULL; }
  if (A->ops) { free(A->ops); A->ops = NULL; }
  free(A); A = NULL;

  return;
}

int SUNMatZero(SUNMatrix A)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(A));
  ier = A->ops->zero(A);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(A));
  return(ier);
}

int SUNMatCopy(SUNMatrix A, SUNMatrix B)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(A));
  ier = A->ops->copy(A, B);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(A));
  return(ier);
}

int SUNMatScaleAdd(realtype c, SUNMatrix A, SUNMatrix B)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(A));
  ier = A->ops->scaleadd(c, A, B);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(A));
  return(ier);
}

int SUNMatScaleAddI(realtype c, SUNMatrix A)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(A));
  ier = A->ops->scaleaddi(c, A);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(A));
  return(ier);
}

int SUNMatMatvecSetup(SUNMatrix A)
{
  int ier = 0;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(A));
  if (A->ops->matvecsetup)
    ier = A->ops->matvecsetup(A);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(A));
  return(ier);
}

int SUNMatMatvec(SUNMatrix A, N_Vector x, N_Vector y)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(A));
  ier = A->ops->matvec(A, x, y);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(A));
  return(ier);
}

int SUNMatSpace(SUNMatrix A, long int *lenrw, long int *leniw)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(A));
  ier = A->ops->space(A, lenrw, leniw);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(A));
  return(ier);
}
