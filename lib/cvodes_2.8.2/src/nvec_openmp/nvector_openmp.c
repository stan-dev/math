/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner and Carol S. Woodward @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR 
 *                   Serial module by Scott D. Cohen, Alan C. 
 *                   Hindmarsh, Radu Serban, and Aaron Collier 
 *                   @ LLNL
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
 * This is the implementation file for an OpenMP implementation
 * of the NVECTOR module.
 * -----------------------------------------------------------------
 */

#include <omp.h>

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_openmp.h>
#include <sundials/sundials_math.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Private function prototypes */
/* z=x */
static void VCopy_OpenMP(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_OpenMP(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_OpenMP(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_OpenMP(N_Vector x, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_OpenMP(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_OpenMP(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_OpenMP(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_OpenMP(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_OpenMP(realtype a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_OpenMP(realtype a, N_Vector x);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new empty vector 
 */

N_Vector N_VNewEmpty_OpenMP(long int length, int num_threads)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_OpenMP content;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvclone           = N_VClone_OpenMP;
  ops->nvcloneempty      = N_VCloneEmpty_OpenMP;
  ops->nvdestroy         = N_VDestroy_OpenMP;
  ops->nvspace           = N_VSpace_OpenMP;
  ops->nvgetarraypointer = N_VGetArrayPointer_OpenMP;
  ops->nvsetarraypointer = N_VSetArrayPointer_OpenMP;
  ops->nvlinearsum       = N_VLinearSum_OpenMP;
  ops->nvconst           = N_VConst_OpenMP;
  ops->nvprod            = N_VProd_OpenMP;
  ops->nvdiv             = N_VDiv_OpenMP;
  ops->nvscale           = N_VScale_OpenMP;
  ops->nvabs             = N_VAbs_OpenMP;
  ops->nvinv             = N_VInv_OpenMP;
  ops->nvaddconst        = N_VAddConst_OpenMP;
  ops->nvdotprod         = N_VDotProd_OpenMP;
  ops->nvmaxnorm         = N_VMaxNorm_OpenMP;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_OpenMP;
  ops->nvwrmsnorm        = N_VWrmsNorm_OpenMP;
  ops->nvmin             = N_VMin_OpenMP;
  ops->nvwl2norm         = N_VWL2Norm_OpenMP;
  ops->nvl1norm          = N_VL1Norm_OpenMP;
  ops->nvcompare         = N_VCompare_OpenMP;
  ops->nvinvtest         = N_VInvTest_OpenMP;
  ops->nvconstrmask      = N_VConstrMask_OpenMP;
  ops->nvminquotient     = N_VMinQuotient_OpenMP;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_OpenMP) malloc(sizeof(struct _N_VectorContent_OpenMP));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  content->length   = length;
  content->num_threads = num_threads;
  content->own_data = FALSE;
  content->data     = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new vector 
 */

N_Vector N_VNew_OpenMP(long int length, int num_threads)
{
  N_Vector v;
  realtype *data;

  v = NULL;
  v = N_VNewEmpty_OpenMP(length, num_threads);
  if (v == NULL) return(NULL);

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_OpenMP(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_OMP(v) = TRUE;
    NV_DATA_OMP(v)     = data;

  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a vector with user data component 
 */

N_Vector N_VMake_OpenMP(long int length, realtype *v_data, int num_threads)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_OpenMP(length, num_threads);
  if (v == NULL) return(NULL);

  if (length > 0) {
    /* Attach data */
    NV_OWN_DATA_OMP(v) = FALSE;
    NV_DATA_OMP(v)     = v_data;
  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new vectors. 
 */

N_Vector *N_VCloneVectorArray_OpenMP(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_OpenMP(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_OpenMP(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new vectors with NULL data array. 
 */

N_Vector *N_VCloneVectorArrayEmpty_OpenMP(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_OpenMP(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_OpenMP(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_OpenMP
 */

void N_VDestroyVectorArray_OpenMP(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_OpenMP(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Function to print a vector 
 */
 
void N_VPrint_OpenMP(N_Vector x)
{
  long int i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%11.8Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%11.8g\n", xd[i]);
#else
    printf("%11.8g\n", xd[i]);
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

/* ----------------------------------------------------------------------------
 * Create new vector from existing vector without attaching data
 */

N_Vector N_VCloneEmpty_OpenMP(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_OpenMP content;

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
  content = (N_VectorContent_OpenMP) malloc(sizeof(struct _N_VectorContent_OpenMP));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  content->length   = NV_LENGTH_OMP(w);
  content->num_threads   = NV_NUM_THREADS_OMP(w);
  content->own_data = FALSE;
  content->data     = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}


/* ----------------------------------------------------------------------------
 * Create new vector from existing vector and attach data
 */

N_Vector N_VClone_OpenMP(N_Vector w)
{
  N_Vector v;
  realtype *data;
  long int length;

  v = NULL;
  v = N_VCloneEmpty_OpenMP(w);
  if (v == NULL) return(NULL);

  length = NV_LENGTH_OMP(w);

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_OpenMP(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_OMP(v) = TRUE;
    NV_DATA_OMP(v)     = data;

  }

  return(v);
}


/* ----------------------------------------------------------------------------
 * Destroy vector and free vector memory
 */

void N_VDestroy_OpenMP(N_Vector v)
{
  if (NV_OWN_DATA_OMP(v) == TRUE) {
    free(NV_DATA_OMP(v));
    NV_DATA_OMP(v) = NULL;
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}


/* ----------------------------------------------------------------------------
 * Get storage requirement for N_Vector
 */

void N_VSpace_OpenMP(N_Vector v, long int *lrw, long int *liw)
{
  *lrw = NV_LENGTH_OMP(v);
  *liw = 1;

  return;
}


/* ----------------------------------------------------------------------------
 * Get vector data pointer
 */

realtype *N_VGetArrayPointer_OpenMP(N_Vector v)
{
  return((realtype *) NV_DATA_OMP(v));
}


/* ----------------------------------------------------------------------------
 * Set vector data pointer
 */

void N_VSetArrayPointer_OpenMP(realtype *v_data, N_Vector v)
{
  if (NV_LENGTH_OMP(v) > 0) NV_DATA_OMP(v) = v_data;

  return;
}


/* ----------------------------------------------------------------------------
 * Compute linear combination z[i] = a*x[i]+b*y[i]
 */

void N_VLinearSum_OpenMP(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy_OpenMP(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy_OpenMP(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_OpenMP(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_OpenMP(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_OpenMP(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_OpenMP(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_OpenMP(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_OpenMP(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
  
  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,a,b,xd,yd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+(b*yd[i]);

  return;
}


/* ----------------------------------------------------------------------------
 * Assigns constant value to all vector elements, z[i] = c
 */

void N_VConst_OpenMP(realtype c, N_Vector z)
{
  long int i, N;
  realtype *zd;

  zd = NULL;

  N  = NV_LENGTH_OMP(z);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,c,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(z))
  for (i = 0; i < N; i++) zd[i] = c;

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise product z[i] = x[i]*y[i]
 */

void N_VProd_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = xd[i]*yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise division z[i] = x[i]/y[i]
 */

void N_VDiv_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = xd[i]/yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute scaler multiplication z[i] = c*x[i]
 */

void N_VScale_OpenMP(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  if (z == x) {  /* BLAS usage: scale x <- cx */
    VScaleBy_OpenMP(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_OpenMP(x, z);
  } else if (c == -ONE) {
    VNeg_OpenMP(x, z);
  } else {
    N  = NV_LENGTH_OMP(x);
    xd = NV_DATA_OMP(x);
    zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,c,xd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
    for (i = 0; i < N; i++) 
      zd[i] = c*xd[i];
  }

  return;
}


/* ----------------------------------------------------------------------------
 * Compute absolute value of vector components z[i] = SUNRabs(x[i])
 */

void N_VAbs_OpenMP(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = SUNRabs(xd[i]);

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = 1 / x[i]
 */

void N_VInv_OpenMP(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,xd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = ONE/xd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise addition of a scaler to a vector z[i] = x[i] + b
 */

void N_VAddConst_OpenMP(N_Vector x, realtype b, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,b,xd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) 
    zd[i] = xd[i]+b;

  return;
}


/* ----------------------------------------------------------------------------
 * Computes the dot product of two vectors, a = sum(x[i]*y[i])
 */

realtype N_VDotProd_OpenMP(N_Vector x, N_Vector y)
{
  long int i, N;
  realtype sum, *xd, *yd;

  sum = ZERO;
  xd = yd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);

#pragma omp parallel for default(none) private(i) shared(N,xd,yd) \
  reduction(+:sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) {
    sum += xd[i]*yd[i];
  }
  
  return(sum);
}


/* ----------------------------------------------------------------------------
 * Computes max norm of a vector 
 */

realtype N_VMaxNorm_OpenMP(N_Vector x)
{
  long int i, N;
  realtype tmax, max, *xd;

  max = ZERO;
  xd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

#pragma omp parallel default(none) private(i,tmax) shared(N,max,xd) \
   num_threads(NV_NUM_THREADS_OMP(x))
  {
    tmax = ZERO;
#pragma omp for schedule(static)
    for (i = 0; i < N; i++) {
      if (SUNRabs(xd[i]) > tmax) tmax = SUNRabs(xd[i]);
    }
#pragma omp critical 
    {
      if (tmax > max)
	max = tmax;
    }
  }
  return(max);
}


/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a vector 
 */

realtype N_VWrmsNorm_OpenMP(N_Vector x, N_Vector w)
{
  long int i, N;
  realtype sum, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  wd = NV_DATA_OMP(w);

#pragma omp parallel for default(none) private(i) shared(N,xd,wd) \
  reduction(+:sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) {
    sum += SUNSQR(xd[i]*wd[i]);
  }

  return(SUNRsqrt(sum/N));
}


/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a masked vector 
 */

realtype N_VWrmsNormMask_OpenMP(N_Vector x, N_Vector w, N_Vector id)
{
  long int i, N;
  realtype sum, *xd, *wd, *idd;

  sum = ZERO;
  xd = wd = idd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd  = NV_DATA_OMP(x);
  wd  = NV_DATA_OMP(w);
  idd = NV_DATA_OMP(id);

#pragma omp parallel for default(none) private(i) shared(N,xd,wd,idd) \
  reduction(+:sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) {
    if (idd[i] > ZERO) {
      sum += SUNSQR(xd[i]*wd[i]);
    }
  }

  return(SUNRsqrt(sum / N));
}


/* ----------------------------------------------------------------------------
 * Finds the minimun component of a vector 
 */

realtype N_VMin_OpenMP(N_Vector x)
{
  long int i, N;
  realtype min, *xd;
  realtype tmin;

  xd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

  min = xd[0];

#pragma omp parallel default(none) private(i,tmin) shared(N,min,xd) \
            num_threads(NV_NUM_THREADS_OMP(x))
  {
    tmin = xd[0];
#pragma omp for schedule(static)
    for (i = 1; i < N; i++) {
      if (xd[i] < tmin) tmin = xd[i];
    }
    if (tmin < min) {
#pragma omp critical
      {
	if (tmin < min) min = tmin;
      }
    }
  }

  return(min);
}


/* ----------------------------------------------------------------------------
 * Computes weighted L2 norm of a vector
 */

realtype N_VWL2Norm_OpenMP(N_Vector x, N_Vector w)
{
  long int i, N;
  realtype sum, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  wd = NV_DATA_OMP(w);

#pragma omp parallel for default(none) private(i) shared(N,xd,wd) \
  reduction(+:sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) {
    sum += SUNSQR(xd[i]*wd[i]);
  }

  return(SUNRsqrt(sum));
}


/* ----------------------------------------------------------------------------
 * Computes L1 norm of a vector
 */

realtype N_VL1Norm_OpenMP(N_Vector x)
{
  long int i, N;
  realtype sum, *xd;

  sum = ZERO;
  xd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

#pragma omp parallel for default(none) private(i) shared(N,xd) \
  reduction(+:sum) schedule(static) num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i<N; i++)  
    sum += SUNRabs(xd[i]);

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Compare vector component values to a scaler   
 */

void N_VCompare_OpenMP(realtype c, N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,c,xd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) {
    zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO;
  }

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = ONE/x[i] and checks if x[i] == ZERO
 */

booleantype N_VInvTest_OpenMP(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd, val;

  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

  val = ZERO;

#pragma omp parallel for default(none) private(i) shared(N,val,xd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO) 
      val = ONE;
    else
      zd[i] = ONE/xd[i];
  }

  if (val > ZERO)
    return (FALSE);
  else
    return (TRUE);
}


/* ----------------------------------------------------------------------------
 * Compute constraint mask of a vector 
 */

booleantype N_VConstrMask_OpenMP(N_Vector c, N_Vector x, N_Vector m)
{
  long int i, N;
  realtype temp;
  realtype *cd, *xd, *md;

  cd = xd = md = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  cd = NV_DATA_OMP(c);
  md = NV_DATA_OMP(m);

  temp = ONE;

#pragma omp parallel for default(none) private(i) shared(N,xd,cd,md,temp) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++) {
    md[i] = ZERO;
    if (cd[i] == ZERO) continue;
    if (cd[i] > ONEPT5 || cd[i] < -ONEPT5) {
      if ( xd[i]*cd[i] <= ZERO) { temp = ZERO; md[i] = ONE; }
      continue;
    }
    if ( cd[i] > HALF || cd[i] < -HALF) {
      if (xd[i]*cd[i] < ZERO ) { temp = ZERO; md[i] = ONE; }
    }
  }

  if (temp == ONE) return (TRUE);
  else return(FALSE);
}


/* ----------------------------------------------------------------------------
 * Compute minimum componentwise quotient 
 */

realtype N_VMinQuotient_OpenMP(N_Vector num, N_Vector denom)
{
  long int i, N;
  realtype *nd, *dd, min, tmin, val;

  nd = dd = NULL;

  N  = NV_LENGTH_OMP(num);
  nd = NV_DATA_OMP(num);
  dd = NV_DATA_OMP(denom);

  min = BIG_REAL;

#pragma omp parallel default(none) private(i,tmin,val) shared(N,min,nd,dd) \
   num_threads(NV_NUM_THREADS_OMP(num))
  {
    tmin = BIG_REAL;
#pragma omp for schedule(static)
    for (i = 0; i < N; i++) {
      if (dd[i] != ZERO) {
	val = nd[i]/dd[i];
	if (val < tmin) tmin = val;
      }
    }
    if (tmin < min) {
#pragma omp critical
      {
	if (tmin < min) min = tmin;
      }
    }
  }

  return(min);
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */


/* ----------------------------------------------------------------------------
 * Copy vector components into a second vector   
 */

static void VCopy_OpenMP(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,xd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = xd[i]; 

  return;
}


/* ----------------------------------------------------------------------------
 * Compute vector sum    
 */

static void VSum_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static) \
 num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = xd[i]+yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute vector difference   
 */

static void VDiff_OpenMP(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,xd,yd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = xd[i]-yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute the negative of a vector   
 */

static void VNeg_OpenMP(N_Vector x, N_Vector z)
{
  long int i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  zd = NV_DATA_OMP(z);
  
#pragma omp parallel for default(none) private(i) shared(N,xd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = -xd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector sum    
 */

static void VScaleSum_OpenMP(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,c,xd,yd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]+yd[i]);

  return;
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector difference
 */

static void VScaleDiff_OpenMP(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,c,xd,yd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]-yd[i]);

  return;
}


/* ----------------------------------------------------------------------------
 * Compute vector sum z[i] = a*x[i]+y[i]    
 */

static void VLin1_OpenMP(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,a,xd,yd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute vector difference z[i] = a*x[i]-y[i]    
 */

static void VLin2_OpenMP(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);
  zd = NV_DATA_OMP(z);

#pragma omp parallel for default(none) private(i) shared(N,a,xd,yd,zd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])-yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute special cases of linear sum
 */

static void Vaxpy_OpenMP(realtype a, N_Vector x, N_Vector y)
{
  long int i, N;
  realtype *xd, *yd;

  xd = yd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);
  yd = NV_DATA_OMP(y);

  if (a == ONE) {
#pragma omp parallel for default(none) private(i) shared(N,xd,yd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
    for (i = 0; i < N; i++)
      yd[i] += xd[i];
    return;
  }

  if (a == -ONE) {
#pragma omp parallel for default(none) private(i) shared(N,xd,yd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
    for (i = 0; i < N; i++)
      yd[i] -= xd[i];
    return;
  }    

#pragma omp parallel for default(none) private(i) shared(N,a,xd,yd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    yd[i] += a*xd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector x[i] = a*x[i]    
 */

static void VScaleBy_OpenMP(realtype a, N_Vector x)
{
  long int i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LENGTH_OMP(x);
  xd = NV_DATA_OMP(x);

#pragma omp parallel for default(none) private(i) shared(N,a,xd) schedule(static) \
   num_threads(NV_NUM_THREADS_OMP(x))
  for (i = 0; i < N; i++)
    xd[i] *= a;

  return;
}
