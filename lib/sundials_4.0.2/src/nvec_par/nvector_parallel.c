/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
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
 * This is the implementation file for a parallel MPI implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_mpi.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Private functions for special cases of vector operations */
static void VCopy_Parallel(N_Vector x, N_Vector z);                              /* z=x       */
static void VSum_Parallel(N_Vector x, N_Vector y, N_Vector z);                   /* z=x+y     */
static void VDiff_Parallel(N_Vector x, N_Vector y, N_Vector z);                  /* z=x-y     */
static void VNeg_Parallel(N_Vector x, N_Vector z);                               /* z=-x      */
static void VScaleSum_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z);  /* z=c(x+y)  */
static void VScaleDiff_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z); /* z=c(x-y)  */
static void VLin1_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax+y    */
static void VLin2_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax-y    */
static void Vaxpy_Parallel(realtype a, N_Vector x, N_Vector y);                  /* y <- ax+y */
static void VScaleBy_Parallel(realtype a, N_Vector x);                           /* x <- ax   */

/* Private functions for special cases of vector array operations */
static int VSumVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z);                   /* Z=X+Y     */
static int VDiffVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z);                  /* Z=X-Y     */
static int VScaleSumVectorArray_Parallel(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z);  /* Z=c(X+Y)  */
static int VScaleDiffVectorArray_Parallel(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z); /* Z=c(X-Y)  */
static int VLin1VectorArray_Parallel(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z);      /* Z=aX+Y    */
static int VLin2VectorArray_Parallel(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z);      /* Z=aX-Y    */
static int VaxpyVectorArray_Parallel(int nvec, realtype a, N_Vector* X, N_Vector* Y);                   /* Y <- aX+Y */

/* Error Message */
#define BAD_N1 "N_VNew_Parallel -- Sum of local vector lengths differs from "
#define BAD_N2 "input global length. \n\n"
#define BAD_N   BAD_N1 BAD_N2

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */

N_Vector_ID N_VGetVectorID_Parallel(N_Vector v)
{
  return SUNDIALS_NVEC_PARALLEL;
}

/* ----------------------------------------------------------------
 * Function to create a new parallel vector with empty data array
 */

N_Vector N_VNewEmpty_Parallel(MPI_Comm comm,
                              sunindextype local_length,
                              sunindextype global_length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Parallel content;
  sunindextype n, Nsum;

  /* Compute global length as sum of local lengths */
  n = local_length;
  MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  if (Nsum != global_length) {
    STAN_SUNDIALS_FPRINTF(stderr, BAD_N);
    return(NULL);
  }

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = N_VGetVectorID_Parallel;
  ops->nvclone           = N_VClone_Parallel;
  ops->nvcloneempty      = N_VCloneEmpty_Parallel;
  ops->nvdestroy         = N_VDestroy_Parallel;
  ops->nvspace           = N_VSpace_Parallel;
  ops->nvgetarraypointer = N_VGetArrayPointer_Parallel;
  ops->nvsetarraypointer = N_VSetArrayPointer_Parallel;

  /* standard vector operations */
  ops->nvlinearsum    = N_VLinearSum_Parallel;
  ops->nvconst        = N_VConst_Parallel;
  ops->nvprod         = N_VProd_Parallel;
  ops->nvdiv          = N_VDiv_Parallel;
  ops->nvscale        = N_VScale_Parallel;
  ops->nvabs          = N_VAbs_Parallel;
  ops->nvinv          = N_VInv_Parallel;
  ops->nvaddconst     = N_VAddConst_Parallel;
  ops->nvdotprod      = N_VDotProd_Parallel;
  ops->nvmaxnorm      = N_VMaxNorm_Parallel;
  ops->nvwrmsnormmask = N_VWrmsNormMask_Parallel;
  ops->nvwrmsnorm     = N_VWrmsNorm_Parallel;
  ops->nvmin          = N_VMin_Parallel;
  ops->nvwl2norm      = N_VWL2Norm_Parallel;
  ops->nvl1norm       = N_VL1Norm_Parallel;
  ops->nvcompare      = N_VCompare_Parallel;
  ops->nvinvtest      = N_VInvTest_Parallel;
  ops->nvconstrmask   = N_VConstrMask_Parallel;
  ops->nvminquotient  = N_VMinQuotient_Parallel;

  /* fused vector operations (optional, NULL means disabled by default) */
  ops->nvlinearcombination = NULL;
  ops->nvscaleaddmulti     = NULL;
  ops->nvdotprodmulti      = NULL;

  /* vector array operations (optional, NULL means disabled by default) */
  ops->nvlinearsumvectorarray         = NULL;
  ops->nvscalevectorarray             = NULL;
  ops->nvconstvectorarray             = NULL;
  ops->nvwrmsnormvectorarray          = NULL;
  ops->nvwrmsnormmaskvectorarray      = NULL;
  ops->nvscaleaddmultivectorarray     = NULL;
  ops->nvlinearcombinationvectorarray = NULL;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_data      = SUNFALSE;
  content->data          = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/* ----------------------------------------------------------------
 * Function to create a new parallel vector
 */

N_Vector N_VNew_Parallel(MPI_Comm comm,
                         sunindextype local_length,
                         sunindextype global_length)
{
  N_Vector v;
  realtype *data;

  v = NULL;
  v = N_VNewEmpty_Parallel(comm, local_length, global_length);
  if (v == NULL) return(NULL);

  /* Create data */
  if(local_length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(local_length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_Parallel(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_P(v) = SUNTRUE;
    NV_DATA_P(v)     = data;

  }

  return(v);
}

/* ----------------------------------------------------------------
 * Function to create a parallel N_Vector with user data component
 */

N_Vector N_VMake_Parallel(MPI_Comm comm,
                          sunindextype local_length,
                          sunindextype global_length,
                          realtype *v_data)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Parallel(comm, local_length, global_length);
  if (v == NULL) return(NULL);

  if (local_length > 0) {
    /* Attach data */
    NV_OWN_DATA_P(v) = SUNFALSE;
    NV_DATA_P(v)     = v_data;
  }

  return(v);
}

/* ----------------------------------------------------------------
 * Function to create an array of new parallel vectors.
 */

N_Vector *N_VCloneVectorArray_Parallel(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_Parallel(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Parallel(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to create an array of new parallel vectors with empty
 * (NULL) data array.
 */

N_Vector *N_VCloneVectorArrayEmpty_Parallel(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_Parallel(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Parallel(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Parallel
 */

void N_VDestroyVectorArray_Parallel(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Parallel(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ----------------------------------------------------------------
 * Function to return global vector length
 */

sunindextype N_VGetLength_Parallel(N_Vector v)
{
  return NV_GLOBLENGTH_P(v);
}

/* ----------------------------------------------------------------
 * Function to return local vector length
 */

sunindextype N_VGetLocalLength_Parallel(N_Vector v)
{
  return NV_LOCLENGTH_P(v);
}

/* ----------------------------------------------------------------
 * Function to print the local data in a parallel vector to stdout
 */

void N_VPrint_Parallel(N_Vector x)
{
  N_VPrintFile_Parallel(x, stdout);
}

/* ----------------------------------------------------------------
 * Function to print the local data in a parallel vector to outfile
 */

void N_VPrintFile_Parallel(N_Vector x, FILE* outfile)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    STAN_SUNDIALS_FPRINTF(outfile, "%Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    STAN_SUNDIALS_FPRINTF(outfile, "%g\n", xd[i]);
#else
    STAN_SUNDIALS_FPRINTF(outfile, "%g\n", xd[i]);
#endif
  }
  STAN_SUNDIALS_FPRINTF(outfile, "\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Parallel(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Parallel content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = w->ops->nvgetvectorid;
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;

  /* standard vector operations */
  ops->nvlinearsum    = w->ops->nvlinearsum;
  ops->nvconst        = w->ops->nvconst;
  ops->nvprod         = w->ops->nvprod;
  ops->nvdiv          = w->ops->nvdiv;
  ops->nvscale        = w->ops->nvscale;
  ops->nvabs          = w->ops->nvabs;
  ops->nvinv          = w->ops->nvinv;
  ops->nvaddconst     = w->ops->nvaddconst;
  ops->nvdotprod      = w->ops->nvdotprod;
  ops->nvmaxnorm      = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm     = w->ops->nvwrmsnorm;
  ops->nvmin          = w->ops->nvmin;
  ops->nvwl2norm      = w->ops->nvwl2norm;
  ops->nvl1norm       = w->ops->nvl1norm;
  ops->nvcompare      = w->ops->nvcompare;
  ops->nvinvtest      = w->ops->nvinvtest;
  ops->nvconstrmask   = w->ops->nvconstrmask;
  ops->nvminquotient  = w->ops->nvminquotient;

  /* fused vector operations */
  ops->nvlinearcombination = w->ops->nvlinearcombination;
  ops->nvscaleaddmulti     = w->ops->nvscaleaddmulti;
  ops->nvdotprodmulti      = w->ops->nvdotprodmulti;

  /* vector array operations */
  ops->nvlinearsumvectorarray         = w->ops->nvlinearsumvectorarray;
  ops->nvscalevectorarray             = w->ops->nvscalevectorarray;
  ops->nvconstvectorarray             = w->ops->nvconstvectorarray;
  ops->nvwrmsnormvectorarray          = w->ops->nvwrmsnormvectorarray;
  ops->nvwrmsnormmaskvectorarray      = w->ops->nvwrmsnormmaskvectorarray;
  ops->nvscaleaddmultivectorarray     = w->ops->nvscaleaddmultivectorarray;
  ops->nvlinearcombinationvectorarray = w->ops->nvlinearcombinationvectorarray;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = NV_LOCLENGTH_P(w);
  content->global_length = NV_GLOBLENGTH_P(w);
  content->comm          = NV_COMM_P(w);
  content->own_data      = SUNFALSE;
  content->data          = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

N_Vector N_VClone_Parallel(N_Vector w)
{
  N_Vector v;
  realtype *data;
  sunindextype local_length;

  v = NULL;
  v = N_VCloneEmpty_Parallel(w);
  if (v == NULL) return(NULL);

  local_length  = NV_LOCLENGTH_P(w);

  /* Create data */
  if(local_length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(local_length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_Parallel(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_P(v) = SUNTRUE;
    NV_DATA_P(v)     = data;
  }

  return(v);
}

void N_VDestroy_Parallel(N_Vector v)
{
  if ((NV_OWN_DATA_P(v) == SUNTRUE) && (NV_DATA_P(v) != NULL)) {
    free(NV_DATA_P(v));
    NV_DATA_P(v) = NULL;
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VSpace_Parallel(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  MPI_Comm comm;
  int npes;

  comm = NV_COMM_P(v);
  MPI_Comm_size(comm, &npes);

  *lrw = NV_GLOBLENGTH_P(v);
  *liw = 2*npes;

  return;
}

realtype *N_VGetArrayPointer_Parallel(N_Vector v)
{
  return((realtype *) NV_DATA_P(v));
}

void N_VSetArrayPointer_Parallel(realtype *v_data, N_Vector v)
{
  if (NV_LOCLENGTH_P(v) > 0) NV_DATA_P(v) = v_data;

  return;
}

void N_VLinearSum_Parallel(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Parallel(a, x, y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy_Parallel(b, y, x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_Parallel(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Parallel(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Parallel(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Parallel(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_Parallel(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_Parallel(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+(b*yd[i]);

  return;
}

void N_VConst_Parallel(realtype c, N_Vector z)
{
  sunindextype i, N;
  realtype *zd;

  zd = NULL;

  N  = NV_LOCLENGTH_P(z);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) zd[i] = c;

  return;
}

void N_VProd_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]*yd[i];

  return;
}

void N_VDiv_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]/yd[i];

  return;
}

void N_VScale_Parallel(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VScaleBy_Parallel(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_Parallel(x, z);
  } else if (c == -ONE) {
    VNeg_Parallel(x, z);
  } else {
    N  = NV_LOCLENGTH_P(x);
    xd = NV_DATA_P(x);
    zd = NV_DATA_P(z);
    for (i = 0; i < N; i++)
      zd[i] = c*xd[i];
  }

  return;
}

void N_VAbs_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = SUNRabs(xd[i]);

  return;
}

void N_VInv_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = ONE/xd[i];

  return;
}

void N_VAddConst_Parallel(N_Vector x, realtype b, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) zd[i] = xd[i]+b;

  return;
}

realtype N_VDotProd_Parallel(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  realtype sum, *xd, *yd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = yd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  comm = NV_COMM_P(x);

  for (i = 0; i < N; i++) sum += xd[i]*yd[i];

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(gsum);
}

realtype N_VMaxNorm_Parallel(N_Vector x)
{
  sunindextype i, N;
  realtype max, *xd, gmax;
  MPI_Comm comm;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  comm = NV_COMM_P(x);

  max = ZERO;

  for (i = 0; i < N; i++) {
    if (SUNRabs(xd[i]) > max) max = SUNRabs(xd[i]);
  }

  gmax = SUNMPI_Allreduce_scalar(max, 2, comm);

  return(gmax);
}

realtype N_VWrmsNorm_Parallel(N_Vector x, N_Vector w)
{
  sunindextype i, N, N_global;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  N        = NV_LOCLENGTH_P(x);
  N_global = NV_GLOBLENGTH_P(x);
  xd       = NV_DATA_P(x);
  wd       = NV_DATA_P(w);
  comm     = NV_COMM_P(x);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(SUNRsqrt(gsum/N_global));
}

realtype N_VWrmsNormMask_Parallel(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N, N_global;
  realtype sum, prodi, *xd, *wd, *idd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = idd = NULL;

  N        = NV_LOCLENGTH_P(x);
  N_global = NV_GLOBLENGTH_P(x);
  xd       = NV_DATA_P(x);
  wd       = NV_DATA_P(w);
  idd      = NV_DATA_P(id);
  comm = NV_COMM_P(x);

  for (i = 0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i]*wd[i];
      sum += SUNSQR(prodi);
    }
  }

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(SUNRsqrt(gsum/N_global));
}

realtype N_VMin_Parallel(N_Vector x)
{
  sunindextype i, N;
  realtype min, *xd, gmin;
  MPI_Comm comm;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  comm = NV_COMM_P(x);

  min = BIG_REAL;

  if (N > 0) {

    xd = NV_DATA_P(x);

    min = xd[0];

    for (i = 1; i < N; i++) {
      if (xd[i] < min) min = xd[i];
    }

  }

  gmin = SUNMPI_Allreduce_scalar(min, 3, comm);

  return(gmin);
}

realtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  wd = NV_DATA_P(w);
  comm = NV_COMM_P(x);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(SUNRsqrt(gsum));
}

realtype N_VL1Norm_Parallel(N_Vector x)
{
  sunindextype i, N;
  realtype sum, gsum, *xd;
  MPI_Comm comm;

  sum = ZERO;
  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  comm = NV_COMM_P(x);

  for (i = 0; i<N; i++)
    sum += SUNRabs(xd[i]);

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(gsum);
}

void N_VCompare_Parallel(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) {
    zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO;
  }

  return;
}

booleantype N_VInvTest_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd, val, gval;
  MPI_Comm comm;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);
  comm = NV_COMM_P(x);

  val = ONE;
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO)
      val = ZERO;
    else
      zd[i] = ONE/xd[i];
  }

  gval = SUNMPI_Allreduce_scalar(val, 3, comm);

  if (gval == ZERO)
    return(SUNFALSE);
  else
    return(SUNTRUE);
}

booleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  realtype temp;
  realtype *cd, *xd, *md;
  booleantype test;
  MPI_Comm comm;

  cd = xd = md = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  cd = NV_DATA_P(c);
  md = NV_DATA_P(m);
  comm = NV_COMM_P(x);

  temp = ZERO;

  for (i = 0; i < N; i++) {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO)
      continue;

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cd[i]) > ONEPT5 && xd[i]*cd[i] <= ZERO) ||
           (SUNRabs(cd[i]) > HALF   && xd[i]*cd[i] <  ZERO);
    if (test) {
      temp = md[i] = ONE;
    }
  }

  /* Find max temp across all MPI ranks */
  temp = SUNMPI_Allreduce_scalar(temp, 2, comm);

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

realtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  sunindextype i, N;
  realtype *nd, *dd, min;
  MPI_Comm comm;

  nd = dd = NULL;

  N  = NV_LOCLENGTH_P(num);
  nd = NV_DATA_P(num);
  dd = NV_DATA_P(denom);
  comm = NV_COMM_P(num);

  notEvenOnce = SUNTRUE;
  min = BIG_REAL;

  for (i = 0; i < N; i++) {
    if (dd[i] == ZERO) continue;
    else {
      if (!notEvenOnce) min = SUNMIN(min, nd[i]/dd[i]);
      else {
        min = nd[i]/dd[i];
        notEvenOnce = SUNFALSE;
      }
    }
  }

  return(SUNMPI_Allreduce_scalar(min, 3, comm));
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination_Parallel(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Parallel(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_Parallel(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_LOCLENGTH_P(z);
  zd = NV_DATA_P(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE)) {
    for (i=1; i<nvec; i++) {
      xd = NV_DATA_P(X[i]);
      for (j=0; j<N; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z) {
    for (j=0; j<N; j++) {
      zd[j] *= c[0];
    }
    for (i=1; i<nvec; i++) {
      xd = NV_DATA_P(X[i]);
      for (j=0; j<N; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = NV_DATA_P(X[0]);
  for (j=0; j<N; j++) {
    zd[j] = c[0] * xd[j];
  }
  for (i=1; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    for (j=0; j<N; j++) {
      zd[j] += c[i] * xd[j];
    }
  }
  return(0);
}


int N_VScaleAddMulti_Parallel(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Parallel(a[0], x, ONE, Y[0], Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      yd = NV_DATA_P(Y[i]);
      for (j=0; j<N; j++) {
        yd[j] += a[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = a[i] * xd[j] + yd[j];
    }
  }
  return(0);
}


int N_VDotProdMulti_Parallel(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VDotProd */
  if (nvec == 1) {
    dotprods[0] = N_VDotProd_Parallel(x, Y[0]);
    return(0);
  }

  /* get vector length, data array, and communicator */
  N    = NV_LOCLENGTH_P(x);
  xd   = NV_DATA_P(x);
  comm = NV_COMM_P(x);

  /* compute multiple dot products */
  for (i=0; i<nvec; i++) {
    yd = NV_DATA_P(Y[i]);
    dotprods[i] = ZERO;
    for (j=0; j<N; j++) {
      dotprods[i] += xd[j] * yd[j];
    }
  }
  SUNMPI_Allreduce(dotprods, nvec, 1, comm);

  return(0);
}


/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Parallel(int nvec,
                                   realtype a, N_Vector* X,
                                   realtype b, N_Vector* Y,
                                   N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;
  realtype     c;
  N_Vector*    V1;
  N_Vector*    V2;
  booleantype  test;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Parallel(a, X[0], b, Y[0], Z[0]);
    return(0);
  }

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (Z == Y))
    return(VaxpyVectorArray_Parallel(nvec, a, X, Y));

  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (Z == X))
    return(VaxpyVectorArray_Parallel(nvec, b, Y, X));

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE))
    return(VSumVectorArray_Parallel(nvec, X, Y, Z));

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VDiffVectorArray_Parallel(nvec, V2, V1, Z));
  }

  /* Cases:                                                  */
  /*   (1) a == 1.0, b == other or 0.0,                      */
  /*   (2) a == other or 0.0, b == 1.0                       */
  /* if a or b is 0.0, then user should have called N_VScale */
  if ((test = (a == ONE)) || (b == ONE)) {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VLin1VectorArray_Parallel(nvec, c, V1, V2, Z));
  }

  /* Cases:                     */
  /*   (1) a == -1.0, b != 1.0, */
  /*   (2) a != 1.0, b == -1.0  */
  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VLin2VectorArray_Parallel(nvec, c, V1, V2, Z));
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
    return(VScaleSumVectorArray_Parallel(nvec, a, X, Y, Z));

  /* Case: a == -b */
  if (a == -b)
    return(VScaleDiffVectorArray_Parallel(nvec, a, X, Y, Z));

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */
  
  /* get vector length */
  N = NV_LOCLENGTH_P(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = a * xd[j] + b * yd[j];
    }
  }

  return(0);
}


int N_VScaleVectorArray_Parallel(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Parallel(c[0], X[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LOCLENGTH_P(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_P(X[i]);
      for (j=0; j<N; j++) {
        xd[j] *= c[i];
      }
    }
    return(0);
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = c[i] * xd[j];
    }
  }
  return(0);
}


int N_VConstVectorArray_Parallel(int nvec, realtype c, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VConst */
  if (nvec == 1) {
    N_VConst_Parallel(c, Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LOCLENGTH_P(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i=0; i<nvec; i++) {
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = c;
    }
  }

  return(0);
}


int N_VWrmsNormVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* W, realtype* nrm)
{
  int          i;
  sunindextype j, Nl, Ng;
  realtype*    wd=NULL;
  realtype*    xd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNorm_Parallel(X[0], W[0]);
    return(0);
  }

  /* get vector lengths and communicator */
  Nl   = NV_LOCLENGTH_P(X[0]);
  Ng   = NV_GLOBLENGTH_P(X[0]);
  comm = NV_COMM_P(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    wd = NV_DATA_P(W[i]);
    nrm[i] = ZERO;
    for (j=0; j<Nl; j++) {
      nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
  }
  SUNMPI_Allreduce(nrm, nvec, 1, comm);

  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/Ng);

  return(0);
}


int N_VWrmsNormMaskVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* W,
                                        N_Vector id, realtype* nrm)
{
  int          i;
  sunindextype j, Nl, Ng;
  realtype*    wd=NULL;
  realtype*    xd=NULL;
  realtype*    idd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNormMask_Parallel(X[0], W[0], id);
    return(0);
  }

  /* get vector lengths, communicator, and mask data */
  Nl   = NV_LOCLENGTH_P(X[0]);
  Ng   = NV_GLOBLENGTH_P(X[0]);
  comm = NV_COMM_P(X[0]);
  idd  = NV_DATA_P(id);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    wd = NV_DATA_P(W[i]);
    nrm[i] = ZERO;
    for (j=0; j<Nl; j++) {
      if (idd[j] > ZERO)
        nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
  }
  SUNMPI_Allreduce(nrm, nvec, 1, comm);

  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/Ng);

  return(0);
}


int N_VScaleAddMultiVectorArray_Parallel(int nvec, int nsum, realtype* a,
                                          N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  int          i, j;
  sunindextype k, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  int          retval;
  N_Vector*    YY;
  N_Vector*    ZZ;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VLinearSum */
    if (nsum == 1) {
      N_VLinearSum_Parallel(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      return(0);
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector *) malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (j=0; j<nsum; j++) {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    retval = N_VScaleAddMulti_Parallel(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);
    return(retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1) {
    retval = N_VLinearSumVectorArray_Parallel(nvec, a[0], X, ONE, Y[0], Z[0]);
    return(retval);
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N  = NV_LOCLENGTH_P(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_P(X[i]);
      for (j=0; j<nsum; j++){
        yd = NV_DATA_P(Y[j][i]);
        for (k=0; k<N; k++) {
          yd[k] += a[j] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    for (j=0; j<nsum; j++) {
      yd = NV_DATA_P(Y[j][i]);
      zd = NV_DATA_P(Z[j][i]);
      for (k=0; k<N; k++) {
        zd[k] = a[j] * xd[k] + yd[k];
      }
    }
  }
  return(0);
}


int N_VLinearCombinationVectorArray_Parallel(int nvec, int nsum,
                                             realtype* c,
                                             N_Vector** X,
                                             N_Vector* Z)
{
  int          i; /* vector arrays index in summation [0,nsum) */
  int          j; /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

  realtype*    ctmp;
  N_Vector*    Y;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VScale */
    if (nsum == 1) {
      N_VScale_Parallel(c[0], X[0][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearSum */
    if (nsum == 2) {
      N_VLinearSum_Parallel(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nsum; i++) {
      Y[i] = X[i][0];
    }

    N_VLinearCombination_Parallel(nsum, c, Y, Z[0]);

    free(Y);
    return(0);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VScaleVectorArray */
  if (nsum == 1) {

    ctmp = (realtype*) malloc(nvec * sizeof(realtype));

    for (j=0; j<nvec; j++) {
      ctmp[j] = c[0];
    }

    N_VScaleVectorArray_Parallel(nvec, ctmp, X[0], Z);

    free(ctmp);
    return(0);
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2) {
    N_VLinearSumVectorArray_Parallel(nvec, c[0], X[0], c[1], X[1], Z);
    return(0);
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LOCLENGTH_P(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE)) {
    for (j=0; j<nvec; j++) {
      zd = NV_DATA_P(Z[j]);
      for (i=1; i<nsum; i++) {
        xd = NV_DATA_P(X[i][j]);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z) {
    for (j=0; j<nvec; j++) {
      zd = NV_DATA_P(Z[j]);
      for (k=0; k<N; k++) {
        zd[k] *= c[0];
      }
      for (i=1; i<nsum; i++) {
        xd = NV_DATA_P(X[i][j]);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j=0; j<nvec; j++) {
    xd = NV_DATA_P(X[0][j]);
    zd = NV_DATA_P(Z[j]);
    for (k=0; k<N; k++) {
      zd[k] = c[0] * xd[k];
    }
    for (i=1; i<nsum; i++) {
      xd = NV_DATA_P(X[i][j]);
      for (k=0; k<N; k++) {
        zd[k] += c[i] * xd[k];
      }
    }
  }
  return(0);
}


/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static void VCopy_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i];

  return;
}

static void VSum_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]+yd[i];

  return;
}

static void VDiff_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]-yd[i];

  return;
}

static void VNeg_Parallel(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = -xd[i];

  return;
}

static void VScaleSum_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]+yd[i]);

  return;
}

static void VScaleDiff_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]-yd[i]);

  return;
}

static void VLin1_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+yd[i];

  return;
}

static void VLin2_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])-yd[i];

  return;
}

static void Vaxpy_Parallel(realtype a, N_Vector x, N_Vector y)
{
  sunindextype i, N;
  realtype *xd, *yd;

  xd = yd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);

  if (a == ONE) {
    for (i = 0; i < N; i++)
      yd[i] += xd[i];
    return;
  }

  if (a == -ONE) {
    for (i = 0; i < N; i++)
      yd[i] -= xd[i];
    return;
  }

  for (i = 0; i < N; i++)
    yd[i] += a*xd[i];

  return;
}

static void VScaleBy_Parallel(realtype a, N_Vector x)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i = 0; i < N; i++)
    xd[i] *= a;

  return;
}


/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static int VSumVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = xd[j] + yd[j];
  }

  return(0);
}

static int VDiffVectorArray_Parallel(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = xd[j] - yd[j];
  }

  return(0);
}

static int VScaleSumVectorArray_Parallel(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = c * (xd[j] + yd[j]);
  }

  return(0);
}

static int VScaleDiffVectorArray_Parallel(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = c * (xd[j] - yd[j]);
  }

  return(0);
}

static int VLin1VectorArray_Parallel(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = (a * xd[j]) + yd[j];
  }

  return(0);
}

static int VLin2VectorArray_Parallel(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  N = NV_LOCLENGTH_P(X[0]);

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    zd = NV_DATA_P(Z[i]);
    for (j=0; j<N; j++)
      zd[j] = (a * xd[j]) - yd[j];
  }

  return(0);
}

static int VaxpyVectorArray_Parallel(int nvec, realtype a, N_Vector* X, N_Vector* Y)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;

  N = NV_LOCLENGTH_P(X[0]);

  if (a == ONE) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_P(X[i]);
      yd = NV_DATA_P(Y[i]);
      for (j=0; j<N; j++)
        yd[j] += xd[j];
    }

    return(0);
  }

  if (a == -ONE) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_P(X[i]);
      yd = NV_DATA_P(Y[i]);
      for (j=0; j<N; j++)
        yd[j] -= xd[j];
    }

    return(0);
  }    

  for (i=0; i<nvec; i++) {
    xd = NV_DATA_P(X[i]);
    yd = NV_DATA_P(Y[i]);
    for (j=0; j<N; j++)
      yd[j] += a * xd[j];
  }

  return(0);
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Parallel;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Parallel;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Parallel;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Parallel;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Parallel;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Parallel;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Parallel;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Parallel;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Parallel;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Parallel;
  } else {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Parallel;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Parallel;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_Parallel;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Parallel;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Parallel;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Parallel;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Parallel;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Parallel;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Parallel;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Parallel(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Parallel;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}
