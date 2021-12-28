/* -----------------------------------------------------------------
 * Programmer(s): Radu Serban, Aaron Collier, and
 *                David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a generic NVECTOR package.
 * It contains the implementation of the N_Vector operations listed
 * in nvector.h.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_nvector.h>
#include "sundials_context_impl.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING)
static SUNProfiler getSUNProfiler(N_Vector v)
{
  return(v->sunctx->profiler);
}
#endif

/* -----------------------------------------------------------------
 * Methods that are not ops (i.e., non-virtual and not overridable).
 * -----------------------------------------------------------------*/

/* Create an empty NVector object */
N_Vector N_VNewEmpty(SUNContext sunctx)
{
  N_Vector     v;
  N_Vector_Ops ops;

  if (sunctx == NULL) return(NULL);

  /* create vector object */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* create vector ops structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof *ops);
  if (ops == NULL) { free(v); return(NULL); }

  /* initialize operations to NULL */

  /*
   * REQUIRED operations.
   *
   * These must be implemented by derivations of the generic N_Vector.
   */

  /* constructors, destructors, and utility operations */
  ops->nvgetvectorid           = NULL;
  ops->nvclone                 = NULL;
  ops->nvcloneempty            = NULL;
  ops->nvdestroy               = NULL;
  ops->nvspace                 = NULL;
  ops->nvgetarraypointer       = NULL;
  ops->nvgetdevicearraypointer = NULL;
  ops->nvsetarraypointer       = NULL;
  ops->nvgetcommunicator       = NULL;
  ops->nvgetlength             = NULL;

  /* standard vector operations */
  ops->nvlinearsum    = NULL;
  ops->nvconst        = NULL;
  ops->nvprod         = NULL;
  ops->nvdiv          = NULL;
  ops->nvscale        = NULL;
  ops->nvabs          = NULL;
  ops->nvinv          = NULL;
  ops->nvaddconst     = NULL;
  ops->nvdotprod      = NULL;
  ops->nvmaxnorm      = NULL;
  ops->nvwrmsnorm     = NULL;
  ops->nvwrmsnormmask = NULL;
  ops->nvmin          = NULL;
  ops->nvwl2norm      = NULL;
  ops->nvl1norm       = NULL;
  ops->nvcompare      = NULL;
  ops->nvinvtest      = NULL;
  ops->nvconstrmask   = NULL;
  ops->nvminquotient  = NULL;

  /*
   * OPTIONAL operations.
   *
   * These operations provide default implementations that may be overriden.
   */

  /* fused vector operations (optional) */
  ops->nvlinearcombination = NULL;
  ops->nvscaleaddmulti     = NULL;
  ops->nvdotprodmulti      = NULL;

  /* vector array operations (optional) */
  ops->nvlinearsumvectorarray         = NULL;
  ops->nvscalevectorarray             = NULL;
  ops->nvconstvectorarray             = NULL;
  ops->nvwrmsnormvectorarray          = NULL;
  ops->nvwrmsnormmaskvectorarray      = NULL;
  ops->nvscaleaddmultivectorarray     = NULL;
  ops->nvlinearcombinationvectorarray = NULL;

  /*
   * OPTIONAL operations with no default implementation.
   */

  /* local reduction operations (optional) */
  ops->nvdotprodlocal      = NULL;
  ops->nvmaxnormlocal      = NULL;
  ops->nvminlocal          = NULL;
  ops->nvl1normlocal       = NULL;
  ops->nvinvtestlocal      = NULL;
  ops->nvconstrmasklocal   = NULL;
  ops->nvminquotientlocal  = NULL;
  ops->nvwsqrsumlocal      = NULL;
  ops->nvwsqrsummasklocal  = NULL;

  /* single buffer reduction operations */
  ops->nvdotprodmultilocal     = NULL;
  ops->nvdotprodmultiallreduce = NULL;

  /* XBraid interface operations */
  ops->nvbufsize   = NULL;
  ops->nvbufpack   = NULL;
  ops->nvbufunpack = NULL;

  /* debugging functions (called when SUNDIALS_DEBUG_PRINTVEC is defined) */
  ops->nvprint     = NULL;
  ops->nvprintfile = NULL;

  /* attach ops */
  v->ops = ops;

  /* initialize content and sunctx to NULL */
  v->content = NULL;
  v->sunctx  = sunctx;

  return(v);
}

/* Free a generic N_Vector (assumes content is already empty) */
void N_VFreeEmpty(N_Vector v)
{
  if (v == NULL)  return;

  /* free non-NULL ops structure */
  if (v->ops)  free(v->ops);
  v->ops = NULL;

  /* free overall N_Vector object and return */
  free(v);
  return;
}

/* Copy a vector 'ops' structure */
int N_VCopyOps(N_Vector w, N_Vector v)
{
  /* Check that ops structures exist */
  if (w == NULL || v == NULL) return(-1);
  if (w->ops == NULL || v->ops == NULL) return(-1);

  /* Copy ops from w to v */

  /*
   * REQUIRED operations.
   *
   * These must be implemented by derivations of the generic N_Vector.
   */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = w->ops->nvgetvectorid;
  v->ops->nvclone                 = w->ops->nvclone;
  v->ops->nvcloneempty            = w->ops->nvcloneempty;
  v->ops->nvdestroy               = w->ops->nvdestroy;
  v->ops->nvspace                 = w->ops->nvspace;
  v->ops->nvgetarraypointer       = w->ops->nvgetarraypointer;
  v->ops->nvgetdevicearraypointer = w->ops->nvgetdevicearraypointer;
  v->ops->nvsetarraypointer       = w->ops->nvsetarraypointer;
  v->ops->nvgetcommunicator       = w->ops->nvgetcommunicator;
  v->ops->nvgetlength             = w->ops->nvgetlength;

  /* standard vector operations */
  v->ops->nvlinearsum    = w->ops->nvlinearsum;
  v->ops->nvconst        = w->ops->nvconst;
  v->ops->nvprod         = w->ops->nvprod;
  v->ops->nvdiv          = w->ops->nvdiv;
  v->ops->nvscale        = w->ops->nvscale;
  v->ops->nvabs          = w->ops->nvabs;
  v->ops->nvinv          = w->ops->nvinv;
  v->ops->nvaddconst     = w->ops->nvaddconst;
  v->ops->nvdotprod      = w->ops->nvdotprod;
  v->ops->nvmaxnorm      = w->ops->nvmaxnorm;
  v->ops->nvwrmsnorm     = w->ops->nvwrmsnorm;
  v->ops->nvwrmsnormmask = w->ops->nvwrmsnormmask;
  v->ops->nvmin          = w->ops->nvmin;
  v->ops->nvwl2norm      = w->ops->nvwl2norm;
  v->ops->nvl1norm       = w->ops->nvl1norm;
  v->ops->nvcompare      = w->ops->nvcompare;
  v->ops->nvinvtest      = w->ops->nvinvtest;
  v->ops->nvconstrmask   = w->ops->nvconstrmask;
  v->ops->nvminquotient  = w->ops->nvminquotient;

  /*
   * OPTIONAL operations.
   *
   * These operations provide default implementations that may be overriden.
   */

  /* fused vector operations */
  v->ops->nvlinearcombination = w->ops->nvlinearcombination;
  v->ops->nvscaleaddmulti     = w->ops->nvscaleaddmulti;
  v->ops->nvdotprodmulti      = w->ops->nvdotprodmulti;

  /* vector array operations */
  v->ops->nvlinearsumvectorarray         = w->ops->nvlinearsumvectorarray;
  v->ops->nvscalevectorarray             = w->ops->nvscalevectorarray;
  v->ops->nvconstvectorarray             = w->ops->nvconstvectorarray;
  v->ops->nvwrmsnormvectorarray          = w->ops->nvwrmsnormvectorarray;
  v->ops->nvwrmsnormmaskvectorarray      = w->ops->nvwrmsnormmaskvectorarray;
  v->ops->nvscaleaddmultivectorarray     = w->ops->nvscaleaddmultivectorarray;
  v->ops->nvlinearcombinationvectorarray = w->ops->nvlinearcombinationvectorarray;

  /*
   * OPTIONAL operations with no default implementation.
   */

  /* local reduction operations */
  v->ops->nvdotprodlocal      = w->ops->nvdotprodlocal;
  v->ops->nvmaxnormlocal      = w->ops->nvmaxnormlocal;
  v->ops->nvminlocal          = w->ops->nvminlocal;
  v->ops->nvl1normlocal       = w->ops->nvl1normlocal;
  v->ops->nvinvtestlocal      = w->ops->nvinvtestlocal;
  v->ops->nvconstrmasklocal   = w->ops->nvconstrmasklocal;
  v->ops->nvminquotientlocal  = w->ops->nvminquotientlocal;
  v->ops->nvwsqrsumlocal      = w->ops->nvwsqrsumlocal;
  v->ops->nvwsqrsummasklocal  = w->ops->nvwsqrsummasklocal;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal     = w->ops->nvdotprodmultilocal;
  v->ops->nvdotprodmultiallreduce = w->ops->nvdotprodmultiallreduce;

  /* XBraid interface operations */
  v->ops->nvbufsize   = w->ops->nvbufsize;
  v->ops->nvbufpack   = w->ops->nvbufpack;
  v->ops->nvbufunpack = w->ops->nvbufunpack;

  /* debugging functions (called when SUNDIALS_DEBUG_PRINTVEC is defined) */
  v->ops->nvprint     = w->ops->nvprint;
  v->ops->nvprintfile = w->ops->nvprintfile;

  return(0);
}

/* -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * -----------------------------------------------------------------*/

N_Vector_ID N_VGetVectorID(N_Vector w)
{
  return(w->ops->nvgetvectorid(w));
}

N_Vector N_VClone(N_Vector w)
{
  N_Vector result = NULL;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(w));
  result = w->ops->nvclone(w);
  result->sunctx = w->sunctx;
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(w));
  return result;
}

N_Vector N_VCloneEmpty(N_Vector w)
{
  N_Vector result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(w));
  result = w->ops->nvcloneempty(w);
  result->sunctx = w->sunctx;
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(w));
  return result;
}

void N_VDestroy(N_Vector v)
{
  if (v == NULL) return;

  /* if the destroy operation exists use it */
  if (v->ops)
    if (v->ops->nvdestroy) { v->ops->nvdestroy(v); return; }

  /* if we reach this point, either ops == NULL or nvdestroy == NULL,
     try to cleanup by freeing the content, ops, and vector */
  if (v->content) { free(v->content); v->content = NULL; }
  if (v->ops) { free(v->ops); v->ops = NULL; }
  free(v); v = NULL;

  return;
}

void N_VSpace(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  v->ops->nvspace(v, lrw, liw);
  return;
}

realtype *N_VGetArrayPointer(N_Vector v)
{
  return((realtype *) v->ops->nvgetarraypointer(v));
}

realtype *N_VGetDeviceArrayPointer(N_Vector v)
{
  if (v->ops->nvgetdevicearraypointer)
    return((realtype *) v->ops->nvgetdevicearraypointer(v));
  else
    return(NULL);
}

void N_VSetArrayPointer(realtype *v_data, N_Vector v)
{
  v->ops->nvsetarraypointer(v_data, v);
  return;
}

void *N_VGetCommunicator(N_Vector v)
{
  if (v->ops->nvgetcommunicator)
    return(v->ops->nvgetcommunicator(v));
  else
    return(NULL);
}

sunindextype N_VGetLength(N_Vector v)
{
  return((sunindextype) v->ops->nvgetlength(v));
}

/* -----------------------------------------------------------------
 * standard vector operations
 * -----------------------------------------------------------------*/

void N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  z->ops->nvlinearsum(a, x, b, y, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return;
}

void N_VConst(realtype c, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(z));
  z->ops->nvconst(c, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(z));
  return;
}

void N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  z->ops->nvprod(x, y, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return;
}

void N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  z->ops->nvdiv(x, y, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return;
}

void N_VScale(realtype c, N_Vector x, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  z->ops->nvscale(c, x, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return;
}

void N_VAbs(N_Vector x, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  z->ops->nvabs(x, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return;
}

void N_VInv(N_Vector x, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  z->ops->nvinv(x, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return;
}

void N_VAddConst(N_Vector x, realtype b, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  z->ops->nvaddconst(x, b, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return;
}

realtype N_VDotProd(N_Vector x, N_Vector y)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) y->ops->nvdotprod(x, y));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VMaxNorm(N_Vector x)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvmaxnorm(x));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VWrmsNorm(N_Vector x, N_Vector w)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvwrmsnorm(x, w));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VWrmsNormMask(N_Vector x, N_Vector w, N_Vector id)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvwrmsnormmask(x, w, id));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VMin(N_Vector x)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvmin(x));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VWL2Norm(N_Vector x, N_Vector w)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvwl2norm(x, w));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VL1Norm(N_Vector x)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvl1norm(x));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

void N_VCompare(realtype c, N_Vector x, N_Vector z)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  z->ops->nvcompare(c, x, z);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return;
}

booleantype N_VInvTest(N_Vector x, N_Vector z)
{
  booleantype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((booleantype) z->ops->nvinvtest(x, z));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

booleantype N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
  booleantype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(c));
  result = ((booleantype) x->ops->nvconstrmask(c, x, m));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(c));
  return(result);
}

realtype N_VMinQuotient(N_Vector num, N_Vector denom)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(num));
  result = ((realtype) num->ops->nvminquotient(num, denom));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(num));
  return(result);
}



/* -----------------------------------------------------------------
 * OPTIONAL fused vector operations
 * -----------------------------------------------------------------*/

int N_VLinearCombination(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int i, ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(X[0]));

  if (z->ops->nvlinearcombination != NULL) {

    ier = z->ops->nvlinearcombination(nvec, c, X, z);

  } else {

    z->ops->nvscale(c[0], X[0], z);
    for (i=1; i<nvec; i++) {
      z->ops->nvlinearsum(c[i], X[i], RCONST(1.0), z, z);
    }
    ier = 0;

  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(X[0]));
  return(ier);
}

int N_VScaleAddMulti(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  int i, ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));

  if (x->ops->nvscaleaddmulti != NULL) {

    ier = x->ops->nvscaleaddmulti(nvec, a, x, Y, Z);

  } else {

    for (i=0; i<nvec; i++) {
      x->ops->nvlinearsum(a[i], x, RCONST(1.0), Y[i], Z[i]);
    }
    ier = 0;

  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(ier);
}

int N_VDotProdMulti(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  int i, ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));

  if (x->ops->nvdotprodmulti != NULL) {

    ier = x->ops->nvdotprodmulti(nvec, x, Y, dotprods);

  } else {

    for (i=0; i<nvec; i++) {
      dotprods[i] = x->ops->nvdotprod(x, Y[i]);
    }
    ier = 0;

  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(ier);
}

/* -----------------------------------------------------------------
 * OPTIONAL vector array operations
 * -----------------------------------------------------------------*/

int N_VLinearSumVectorArray(int nvec, realtype a, N_Vector* X,
                            realtype b, N_Vector* Y, N_Vector* Z)
{
  int i, ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(X[0]));

  if (Z[0]->ops->nvlinearsumvectorarray != NULL) {

    ier = Z[0]->ops->nvlinearsumvectorarray(nvec, a, X, b, Y, Z);

  } else {

    for (i=0; i<nvec; i++) {
      Z[0]->ops->nvlinearsum(a, X[i], b, Y[i], Z[i]);
    }
    ier = 0;

  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(X[0]));
  return(ier);
}

int N_VScaleVectorArray(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int i, ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(X[0]));

  if (Z[0]->ops->nvscalevectorarray != NULL) {

    ier = Z[0]->ops->nvscalevectorarray(nvec, c, X, Z);

  } else {

    for (i=0; i<nvec; i++) {
      Z[0]->ops->nvscale(c[i], X[i], Z[i]);
    }
    ier = 0;

  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(X[0]));
  return(ier);
}

int N_VConstVectorArray(int nvec, realtype c, N_Vector* Z)
{
  int i, ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(Z[0]));

  if (Z[0]->ops->nvconstvectorarray != NULL) {

    ier = Z[0]->ops->nvconstvectorarray(nvec, c, Z);

  } else {

    for (i=0; i<nvec; i++) {
      Z[0]->ops->nvconst(c, Z[i]);
    }
    ier = 0;

  }

  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(Z[0]));
  return(ier);
}

int N_VWrmsNormVectorArray(int nvec, N_Vector* X, N_Vector* W, realtype* nrm)
{
  int i, ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(X[0]));

  if (X[0]->ops->nvwrmsnormvectorarray != NULL) {

    ier = X[0]->ops->nvwrmsnormvectorarray(nvec, X, W, nrm);

  } else {

    for (i=0; i<nvec; i++) {
      nrm[i] = X[0]->ops->nvwrmsnorm(X[i], W[i]);
    }
    ier = 0;

  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(X[0]));
  return(ier);
}

int N_VWrmsNormMaskVectorArray(int nvec, N_Vector* X, N_Vector* W, N_Vector id,
                               realtype* nrm)
{
  int i, ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(X[0]));

  if (id->ops->nvwrmsnormmaskvectorarray != NULL) {

    ier = id->ops->nvwrmsnormmaskvectorarray(nvec, X, W, id, nrm);

  } else {

    for (i=0; i<nvec; i++) {
      nrm[i] = id->ops->nvwrmsnormmask(X[i], W[i], id);
    }
    ier = 0;

  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(X[0]));
  return(ier);
}

int N_VScaleAddMultiVectorArray(int nvec, int nsum, realtype* a, N_Vector* X,
                                 N_Vector** Y, N_Vector** Z)
{
  int i, j, ier;
  N_Vector* YY = NULL;
  N_Vector* ZZ = NULL;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(X[0]));

  if (X[0]->ops->nvscaleaddmultivectorarray != NULL) {

    ier = X[0]->ops->nvscaleaddmultivectorarray(nvec, nsum, a, X, Y, Z);

  } else if (X[0]->ops->nvscaleaddmulti != NULL ) {

    /* allocate arrays of vectors */
    YY = (N_Vector*) malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector*) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nvec; i++) {

      for (j=0; j<nsum; j++) {
        YY[j] = Y[j][i];
        ZZ[j] = Z[j][i];
      }

      ier = X[0]->ops->nvscaleaddmulti(nsum, a, X[i], YY, ZZ);
      if (ier != 0) break;
    }

    /* free array of vectors */
    free(YY);
    free(ZZ);

  } else {

    for (i=0; i<nvec; i++) {
      for (j=0; j<nsum; j++) {
        X[0]->ops->nvlinearsum(a[j], X[i], RCONST(1.0), Y[j][i], Z[j][i]);
      }
    }
    ier = 0;
  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(X[0]));
  return(ier);
}

int N_VLinearCombinationVectorArray(int nvec, int nsum, realtype* c,
                                    N_Vector** X, N_Vector* Z)
{
  int i, j, ier;
  N_Vector* Y = NULL;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(X[0][0]));

  if (Z[0]->ops->nvlinearcombinationvectorarray != NULL) {

    ier = Z[0]->ops->nvlinearcombinationvectorarray(nvec, nsum, c, X, Z);

  } else if (Z[0]->ops->nvlinearcombination != NULL ) {

    /* allocate array of vectors */
    Y = (N_Vector* ) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nvec; i++) {

      for (j=0; j<nsum; j++) {
        Y[j] = X[j][i];
      }

      ier = Z[0]->ops->nvlinearcombination(nsum, c, Y, Z[i]);
      if (ier != 0) break;
    }

    /* free array of vectors */
    free(Y);

  } else {

    for (i=0; i<nvec; i++) {
      Z[0]->ops->nvscale(c[0], X[0][i], Z[i]);
      for (j=1; j<nsum; j++) {
        Z[0]->ops->nvlinearsum(c[j], X[j][i], RCONST(1.0), Z[i], Z[i]);
      }
    }
    ier = 0;
  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(X[0][0]));
  return(ier);
}

/* -----------------------------------------------------------------
 * OPTIONAL local reduction kernels (no parallel communication)
 * -----------------------------------------------------------------*/

realtype N_VDotProdLocal(N_Vector x, N_Vector y)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) y->ops->nvdotprodlocal(x, y));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VMaxNormLocal(N_Vector x)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvmaxnormlocal(x));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VMinLocal(N_Vector x)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvminlocal(x));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VL1NormLocal(N_Vector x)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvl1normlocal(x));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VWSqrSumLocal(N_Vector x, N_Vector w)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvwsqrsumlocal(x,w));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VWSqrSumMaskLocal(N_Vector x, N_Vector w, N_Vector id)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((realtype) x->ops->nvwsqrsummasklocal(x,w,id));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

booleantype N_VInvTestLocal(N_Vector x, N_Vector z)
{
  booleantype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((booleantype) z->ops->nvinvtestlocal(x,z));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

booleantype N_VConstrMaskLocal(N_Vector c, N_Vector x, N_Vector m)
{
  booleantype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  result = ((booleantype) x->ops->nvconstrmasklocal(c,x,m));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(result);
}

realtype N_VMinQuotientLocal(N_Vector num, N_Vector denom)
{
  realtype result;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(num));
  result = ((realtype) num->ops->nvminquotientlocal(num,denom));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(num));
  return(result);
}

/* -------------------------------------------
 * OPTIONAL single buffer reduction operations
 * -------------------------------------------*/

int N_VDotProdMultiLocal(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  int i;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));

  if (x->ops->nvdotprodmultilocal)
    return((int) x->ops->nvdotprodmultilocal(nvec, x, Y, dotprods));

  if (x->ops->nvdotprodlocal) {
    for (i = 0; i < nvec; i++) {
      dotprods[i] = x->ops->nvdotprodlocal(x, Y[i]);
    }
    return(0);
  }

  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(-1);
}

int N_VDotProdMultiAllReduce(int nvec, N_Vector x, realtype* sum)
{
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  if (x->ops->nvdotprodmultiallreduce)
    return(x->ops->nvdotprodmultiallreduce(nvec, x, sum));
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(-1);
}

/* ------------------------------------
 * OPTIONAL XBraid interface operations
 * ------------------------------------*/

int N_VBufSize(N_Vector x, sunindextype *size)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  if (x->ops->nvbufsize == NULL)
    ier = -1;
  else
    ier = x->ops->nvbufsize(x, size);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(ier);
}

int N_VBufPack(N_Vector x, void *buf)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  if (x->ops->nvbufpack == NULL)
    ier = -1;
  else
    ier = x->ops->nvbufpack(x, buf);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(ier);
}

int N_VBufUnpack(N_Vector x, void *buf)
{
  int ier;
  SUNDIALS_MARK_FUNCTION_BEGIN(getSUNProfiler(x));
  if (x->ops->nvbufunpack == NULL)
    ier = -1;
  else
    ier = x->ops->nvbufunpack(x, buf);
  SUNDIALS_MARK_FUNCTION_END(getSUNProfiler(x));
  return(ier);
}

/* -----------------------------------------------------------------
 * Additional functions exported by the generic NVECTOR:
 *   N_VNewVectorArray
 *   N_VCloneEmptyVectorArray
 *   N_VCloneVectorArray
 *   N_VDestroyVectorArray
 * -----------------------------------------------------------------*/
N_Vector* N_VNewVectorArray(int count)
{
  N_Vector* vs = NULL;

  if (count <= 0) return(NULL);

  vs = (N_Vector* ) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  return(vs);
}

N_Vector* N_VCloneEmptyVectorArray(int count, N_Vector w)
{
  N_Vector* vs = NULL;
  int j;

  if (count <= 0) return(NULL);

  vs = (N_Vector* ) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = N_VCloneEmpty(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

N_Vector* N_VCloneVectorArray(int count, N_Vector w)
{
  N_Vector* vs = NULL;
  int j;

  if (count <= 0) return(NULL);

  vs = (N_Vector* ) malloc(count * sizeof(N_Vector));
  if (vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = N_VClone(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

void N_VDestroyVectorArray(N_Vector* vs, int count)
{
  int j;

  if (vs == NULL) return;

  for (j = 0; j < count; j++) {
    N_VDestroy(vs[j]);
    vs[j] = NULL;
  }

  free(vs); vs = NULL;

  return;
}

/* These function are really only for users of the Fortran interface */
N_Vector N_VGetVecAtIndexVectorArray(N_Vector* vs, int index)
{
  if (vs==NULL)       return NULL;
  else if (index < 0) return NULL;
  else                return vs[index];
}

void N_VSetVecAtIndexVectorArray(N_Vector* vs, int index, N_Vector w)
{
  if (vs==NULL)       return;
  else if (index < 0) return;
  else                vs[index] = w;
}


/* -----------------------------------------------------------------
 * Debugging functions
 * ----------------------------------------------------------------- */

void N_VPrint(N_Vector v)
{
  if (v == NULL) {
   STAN_SUNDIALS_PRINTF("NULL Vector\n");
    return;
  }
  if (v->ops->nvprint == NULL) {
   STAN_SUNDIALS_PRINTF("NULL Print Op\n");
    return;
  }
  v->ops->nvprint(v);
  return;
}


void N_VPrintFile(N_Vector v, FILE* outfile)
{
  if (v == NULL) {
    STAN_SUNDIALS_FPRINTF(outfile, "NULL Vector\n");
    return;
  }
  if (v->ops->nvprintfile == NULL) {
    STAN_SUNDIALS_FPRINTF(outfile, "NULL PrintFile Op\n");
    return;
  }
  v->ops->nvprintfile(v, outfile);
  return;
}
