/* --------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * --------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * --------------------------------------------------------------------------
 * This is the implementation file for the SUNDIALS + XBraid interface.
 * -------------------------------------------------------------------------- */

#include "sundials/sundials_xbraid.h"
#include "sundials/sundials_math.h"

#define ONE RCONST(1.0)


/* -------------------------
 * Create and free utilities
 * ------------------------- */


/* Create an empty SUNBraidApp instance */
int SUNBraidApp_NewEmpty(braid_App *app)
{
  SUNBraidOps ops;

  /* Create XBraid interface object */
  *app = NULL;
  *app = (braid_App) malloc(sizeof(struct _braid_App_struct));
  if (*app == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create operations structure */
  ops = NULL;
  ops = (SUNBraidOps) malloc(sizeof(struct _SUNBraidOps));
  if (ops == NULL)
  {
    free(*app);
    *app = NULL;
    return SUNBRAID_ALLOCFAIL;
  }

  /* Initialize operations to NULL */
  ops->getvectmpl = NULL;

  /* Attach operations and initialize content to NULL */
  (*app)->ops     = ops;
  (*app)->content = NULL;

  return SUNBRAID_SUCCESS;
}


/* Free and empty SUNBraidApp instance */
int SUNBraidApp_FreeEmpty(braid_App *app)
{
  if (*app == NULL) return SUNBRAID_SUCCESS;

  if ((*app)->ops) free((*app)->ops);
  (*app)->ops = NULL;

  free(*app);
  *app = NULL;

  return SUNBRAID_SUCCESS;
}


/* ----------------------
 * Generic app operations
 * ---------------------- */


/* Get a template vector from the integrator */
int SUNBraidApp_GetVecTmpl(braid_App app, N_Vector *y)
{
  if (app->ops->getvectmpl == NULL) return SUNBRAID_OPNULL;
  return app->ops->getvectmpl(app, y);
}


/* -------------------------
 * SUNBraid Vector Functions
 * ------------------------- */


/* Create a new vector wrapper */
int SUNBraidVector_New(N_Vector y, SUNBraidVector *u)
{
  /* Check for valid N_Vector */
  if (y == NULL) return SUNBRAID_ILLINPUT;

  /* Create new vector wrapper */
  *u = NULL;
  *u = (SUNBraidVector) malloc(sizeof(struct _braid_Vector_struct));
  if (*u == NULL) return SUNBRAID_ALLOCFAIL;

  /* Attach N_Vector */
  (*u)->y = y;

  return SUNBRAID_SUCCESS;
}


/* Get the wrapped NVector */
int SUNBraidVector_GetNVector(SUNBraidVector u, N_Vector *y)
{
  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;
  if (u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Extract NVector */
  *y = u->y;

  return SUNBRAID_SUCCESS;
}


/* Create clone of an existing vector */
int SUNBraidVector_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr)
{
  int      flag;
  N_Vector vy;

  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;
  if (u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Clone input NVector */
  vy = N_VClone(u->y);
  if (vy == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create new vector wrapper */
  flag = SUNBraidVector_New(vy, v_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Copy data from u to v */
  N_VScale(ONE, u->y, vy);

  return SUNBRAID_SUCCESS;
}


/* Free vector */
int SUNBraidVector_Free(braid_App app, braid_Vector u)
{
  /* Check for valid input */
  if (u == NULL) return SUNBRAID_SUCCESS;

  /* Destroy N_Vector */
  if (u->y != NULL)
  {
    N_VDestroy(u->y);
    u->y = NULL;
  }

  /* Destroy SUNBraidVector wrapper */
  free(u);
  u = NULL;

  return SUNBRAID_SUCCESS;
}


/* Compute alpha x + beta y -> y */
int SUNBraidVector_Sum(braid_App app, braid_Real alpha, braid_Vector x,
                       braid_Real beta, braid_Vector y)
{
  /* Check for valid wrappers */
  if (x == NULL || y == NULL) return SUNBRAID_ILLINPUT;
  if (x->y == NULL || y->y == NULL) return SUNBRAID_MEMFAIL;

  /* Compute linear sum */
  N_VLinearSum(alpha, x->y, beta, y->y, y->y);
  return SUNBRAID_SUCCESS;
}


/* Compute L2 norm */
int SUNBraidVector_SpatialNorm(braid_App app, braid_Vector u,
                               braid_Real *norm_ptr)
{
  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;
  if (u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Compute L2 norm */
  *norm_ptr = SUNRsqrt(N_VDotProd(u->y, u->y));

  return SUNBRAID_SUCCESS;
}


/* Compute message buffer size */
int SUNBraidVector_BufSize(braid_App app, braid_Int *size_ptr,
                           braid_BufferStatus bstatus)
{
  int      flag;  /* return flag     */
  N_Vector ytmpl; /* template vector */

  /* Get template vector */
  flag = SUNBraidApp_GetVecTmpl(app, &ytmpl);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Get buffer size */
  flag = N_VBufSize(ytmpl, size_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  return SUNBRAID_SUCCESS;
}


/* Pack message buffer */
int SUNBraidVector_BufPack(braid_App app, braid_Vector u, void *buffer,
                           braid_BufferStatus bstatus)
{
  int flag; /* return flag */

  /* Check for valid wrapper */
  if (u == NULL) return SUNBRAID_ILLINPUT;
  if (u->y == NULL) return SUNBRAID_MEMFAIL;

  /* Fill buffer */
  flag = N_VBufPack(u->y, buffer);
  if (flag != SUNBRAID_SUCCESS) return flag;

  return SUNBRAID_SUCCESS;
}


/* Unpack message buffer */
int SUNBraidVector_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr,
                             braid_BufferStatus bstatus)
{
  int      flag;  /* return flag     */
  N_Vector ytmpl; /* template vector */
  N_Vector y;     /* new NVector     */

  /* Check for valid input */
  if (buffer == NULL) return SUNBRAID_ILLINPUT;

  /* Get template vector */
  flag = SUNBraidApp_GetVecTmpl(app, &ytmpl);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Create new NVector */
  y = N_VClone(ytmpl);
  if (y == NULL) return SUNBRAID_ALLOCFAIL;

  /* Create new XBraid vector */
  flag = SUNBraidVector_New(y, u_ptr);
  if (flag != SUNBRAID_SUCCESS) return flag;

  /* Unpack buffer */
  flag = N_VBufUnpack(y, buffer);
  if (flag != SUNBRAID_SUCCESS) return flag;

  return SUNBRAID_SUCCESS;
}
