/*
 * -----------------------------------------------------------------
 * $Revision: 4396 $
 * $Date: 2015-02-26 16:59:39 -0800 (Thu, 26 Feb 2015) $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
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
 * This is the implementation file for a parallel MPI implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_parallel.h>
#include <sundials/sundials_math.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Error Message */

#define BAD_N1 "N_VNew_Parallel -- Sum of local vector lengths differs from "
#define BAD_N2 "input global length. \n\n"
#define BAD_N   BAD_N1 BAD_N2

/* Private function prototypes */

/* Reduction operations add/max/min over the processor group */
static realtype VAllReduce_Parallel(realtype d, int op, MPI_Comm comm);
/* z=x */
static void VCopy_Parallel(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_Parallel(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_Parallel(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_Parallel(N_Vector x, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_Parallel(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_Parallel(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_Parallel(realtype a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_Parallel(realtype a, N_Vector x);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Function to create a new parallel vector with empty data array
 */

N_Vector N_VNewEmpty_Parallel(MPI_Comm comm, 
                              long int local_length,
                              long int global_length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Parallel content;
  long int n, Nsum;

  /* Compute global length as sum of local lengths */
  n = local_length;
  MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  if (Nsum != global_length) {
    printf(BAD_N);
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

  ops->nvclone           = N_VClone_Parallel;
  ops->nvcloneempty      = N_VCloneEmpty_Parallel;
  ops->nvdestroy         = N_VDestroy_Parallel;
  ops->nvspace           = N_VSpace_Parallel;
  ops->nvgetarraypointer = N_VGetArrayPointer_Parallel;
  ops->nvsetarraypointer = N_VSetArrayPointer_Parallel;
  ops->nvlinearsum       = N_VLinearSum_Parallel;
  ops->nvconst           = N_VConst_Parallel;
  ops->nvprod            = N_VProd_Parallel;
  ops->nvdiv             = N_VDiv_Parallel;
  ops->nvscale           = N_VScale_Parallel;
  ops->nvabs             = N_VAbs_Parallel;
  ops->nvinv             = N_VInv_Parallel;
  ops->nvaddconst        = N_VAddConst_Parallel;
  ops->nvdotprod         = N_VDotProd_Parallel;
  ops->nvmaxnorm         = N_VMaxNorm_Parallel;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Parallel;
  ops->nvwrmsnorm        = N_VWrmsNorm_Parallel;
  ops->nvmin             = N_VMin_Parallel;
  ops->nvwl2norm         = N_VWL2Norm_Parallel;
  ops->nvl1norm          = N_VL1Norm_Parallel;
  ops->nvcompare         = N_VCompare_Parallel;
  ops->nvinvtest         = N_VInvTest_Parallel;
  ops->nvconstrmask      = N_VConstrMask_Parallel;
  ops->nvminquotient     = N_VMinQuotient_Parallel;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_data      = FALSE;
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
                         long int local_length,
                         long int global_length)
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
    NV_OWN_DATA_P(v) = TRUE;
    NV_DATA_P(v)     = data; 

  }

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create a parallel N_Vector with user data component 
 */

N_Vector N_VMake_Parallel(MPI_Comm comm, 
                          long int local_length,
                          long int global_length,
                          realtype *v_data)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Parallel(comm, local_length, global_length);
  if (v == NULL) return(NULL);

  if (local_length > 0) {
    /* Attach data */
    NV_OWN_DATA_P(v) = FALSE;
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
 * Function to print a parallel vector 
 */

void N_VPrint_Parallel(N_Vector x)
{
  long int i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%g\n", xd[i]);
#else
    printf("%g\n", xd[i]);
#endif
  }
  printf("\n");

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
  
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;  
  ops->nvprod            = w->ops->nvprod;   
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale; 
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;    
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

  /* Create content */  
  content = NULL;
  content = (N_VectorContent_Parallel) malloc(sizeof(struct _N_VectorContent_Parallel));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = NV_LOCLENGTH_P(w);
  content->global_length = NV_GLOBLENGTH_P(w);
  content->comm          = NV_COMM_P(w);
  content->own_data      = FALSE;
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
  long int local_length;

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
    NV_OWN_DATA_P(v) = TRUE;
    NV_DATA_P(v)     = data;
  }

  return(v);
}

void N_VDestroy_Parallel(N_Vector v)
{
  if ((NV_OWN_DATA_P(v) == TRUE) && (NV_DATA_P(v) != NULL)) {
    free(NV_DATA_P(v));
    NV_DATA_P(v) = NULL;
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VSpace_Parallel(N_Vector v, long int *lrw, long int *liw)
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
  long int i, N;
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
  long int i, N;
  realtype *zd;

  zd = NULL;

  N  = NV_LOCLENGTH_P(z);
  zd = NV_DATA_P(z);

  for (i = 0; i < N; i++) zd[i] = c;

  return;
}

void N_VProd_Parallel(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
  realtype sum, *xd, *yd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = yd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  yd = NV_DATA_P(y);
  comm = NV_COMM_P(x);

  for (i = 0; i < N; i++) sum += xd[i]*yd[i];

  gsum = VAllReduce_Parallel(sum, 1, comm);

  return(gsum);
}

realtype N_VMaxNorm_Parallel(N_Vector x)
{
  long int i, N;
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
   
  gmax = VAllReduce_Parallel(max, 2, comm);

  return(gmax);
}

realtype N_VWrmsNorm_Parallel(N_Vector x, N_Vector w)
{
  long int i, N, N_global;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  N        = NV_LOCLENGTH_P(x);
  N_global = NV_GLOBLENGTH_P(x);
  xd       = NV_DATA_P(x);
  wd       = NV_DATA_P(w);
  comm = NV_COMM_P(x);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  gsum = VAllReduce_Parallel(sum, 1, comm);

  return(SUNRsqrt(gsum/N_global));
}

realtype N_VWrmsNormMask_Parallel(N_Vector x, N_Vector w, N_Vector id)
{
  long int i, N, N_global;
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

  gsum = VAllReduce_Parallel(sum, 1, comm);

  return(SUNRsqrt(gsum/N_global));
}

realtype N_VMin_Parallel(N_Vector x)
{
  long int i, N;
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

  gmin = VAllReduce_Parallel(min, 3, comm);

  return(gmin);
}

realtype N_VWL2Norm_Parallel(N_Vector x, N_Vector w)
{
  long int i, N;
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

  gsum = VAllReduce_Parallel(sum, 1, comm);

  return(SUNRsqrt(gsum));
}

realtype N_VL1Norm_Parallel(N_Vector x)
{
  long int i, N;
  realtype sum, gsum, *xd;
  MPI_Comm comm;

  sum = ZERO;
  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  comm = NV_COMM_P(x);

  for (i = 0; i<N; i++) 
    sum += SUNRabs(xd[i]);

  gsum = VAllReduce_Parallel(sum, 1, comm);

  return(gsum);
}

void N_VCompare_Parallel(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
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
  long int i, N;
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

  gval = VAllReduce_Parallel(val, 3, comm);

  if (gval == ZERO)
    return(FALSE);
  else
    return(TRUE);
}

booleantype N_VConstrMask_Parallel(N_Vector c, N_Vector x, N_Vector m)
{
  long int i, N;
  realtype temp;
  realtype *cd, *xd, *md;
  MPI_Comm comm;

  cd = xd = md = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);
  cd = NV_DATA_P(c);
  md = NV_DATA_P(m);
  comm = NV_COMM_P(x);

  temp = ONE;

  for (i = 0; i < N; i++) {
    md[i] = ZERO;
    if (cd[i] == ZERO) continue;
    if (cd[i] > ONEPT5 || cd[i] < -ONEPT5) {
      if (xd[i]*cd[i] <= ZERO) { temp = ZERO; md[i] = ONE; }
      continue;
    }
    if (cd[i] > HALF || cd[i] < -HALF) {
      if (xd[i]*cd[i] < ZERO ) { temp = ZERO; md[i] = ONE; }
    }
  }

  temp = VAllReduce_Parallel(temp, 3, comm);

  if (temp == ONE) return(TRUE);
  else return(FALSE);
}

realtype N_VMinQuotient_Parallel(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  long int i, N;
  realtype *nd, *dd, min;
  MPI_Comm comm;

  nd = dd = NULL;

  N  = NV_LOCLENGTH_P(num);
  nd = NV_DATA_P(num);
  dd = NV_DATA_P(denom);
  comm = NV_COMM_P(num);

  notEvenOnce = TRUE;
  min = BIG_REAL;

  for (i = 0; i < N; i++) {
    if (dd[i] == ZERO) continue;
    else {
      if (!notEvenOnce) min = SUNMIN(min, nd[i]/dd[i]);
      else {
        min = nd[i]/dd[i];
        notEvenOnce = FALSE;
      }
    }
  }

  return(VAllReduce_Parallel(min, 3, comm));
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static realtype VAllReduce_Parallel(realtype d, int op, MPI_Comm comm)
{
  /* 
   * This function does a global reduction.  The operation is
   *   sum if op = 1,
   *   max if op = 2,
   *   min if op = 3.
   * The operation is over all processors in the communicator 
   */

  realtype out;

  switch (op) {
   case 1: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);
}

static void VCopy_Parallel(N_Vector x, N_Vector z)
{
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
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
  long int i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LOCLENGTH_P(x);
  xd = NV_DATA_P(x);

  for (i = 0; i < N; i++)
    xd[i] *= a;

  return;
}
