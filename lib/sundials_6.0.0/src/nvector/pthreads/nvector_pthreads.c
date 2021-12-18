/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR
 *                   Serial module by Scott D. Cohen, Alan C.
 *                   Hindmarsh, Radu Serban, and Aaron Collier
 *                   @ LLNL
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
 * This is the implementation file for a POSIX Threads (Pthreads)
 * implementation of the NVECTOR package using a LOCAL array of
 * structures to pass data to threads.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_pthreads.h>
#include <sundials/sundials_math.h>
#include <math.h> /* define NAN */

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Private functions for special cases of vector operations */
static void VCopy_Pthreads(N_Vector x, N_Vector z);                              /* z=x       */
static void VSum_Pthreads(N_Vector x, N_Vector y, N_Vector z);                   /* z=x+y     */
static void VDiff_Pthreads(N_Vector x, N_Vector y, N_Vector z);                  /* z=x-y     */
static void VNeg_Pthreads(N_Vector x, N_Vector z);                               /* z=-x      */
static void VScaleSum_Pthreads(realtype c, N_Vector x, N_Vector y, N_Vector z);  /* z=c(x+y)  */
static void VScaleDiff_Pthreads(realtype c, N_Vector x, N_Vector y, N_Vector z); /* z=c(x-y)  */
static void VLin1_Pthreads(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax+y    */
static void VLin2_Pthreads(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax-y    */
static void Vaxpy_Pthreads(realtype a, N_Vector x, N_Vector y);                  /* y <- ax+y */
static void VScaleBy_Pthreads(realtype a, N_Vector x);                           /* x <- ax   */

/* Private functions for special cases of vector array operations */
static int VSumVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z);                   /* Z=X+Y     */
static int VDiffVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z);                  /* Z=X-Y     */
static int VScaleSumVectorArray_Pthreads(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z);  /* Z=c(X+Y)  */
static int VScaleDiffVectorArray_Pthreads(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z); /* Z=c(X-Y)  */
static int VLin1VectorArray_Pthreads(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z);      /* Z=aX+Y    */
static int VLin2VectorArray_Pthreads(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z);      /* Z=aX-Y    */
static int VaxpyVectorArray_Pthreads(int nvec, realtype a, N_Vector* X, N_Vector* Y);                    /* Y <- aX+Y */

/* Pthread companion functions for vector operations */
static void *N_VLinearSum_PT(void *thread_data);
static void *N_VConst_PT(void *thread_data);
static void *N_VProd_PT(void *thread_data);
static void *N_VDiv_PT(void *thread_data);
static void *N_VScale_PT(void *thread_data);
static void *N_VAbs_PT(void *thread_data);
static void *N_VInv_PT(void *thread_data);
static void *N_VAddConst_PT(void *thread_data);
static void *N_VCompare_PT(void *thread_data);
static void *N_VDotProd_PT(void *thread_data);
static void *N_VMaxNorm_PT(void *thread_data);
static void *N_VWSqrSum_PT(void *thread_data);
static void *N_VMin_PT(void *thread_data);
static void *N_VWL2Norm_PT(void *thread_data);
static void *N_VL1Norm_PT(void *thread_data);
static void *N_VInvTest_PT(void *thread_data);
static void *N_VWSqrSumMask_PT(void *thread_data);
static void *N_VConstrMask_PT(void *thread_data);
static void *N_VMinQuotient_PT(void *thread_data);

/* Pthread companion functions special cases of vector operations */
static void *VCopy_PT(void *thread_data);
static void *VSum_PT(void *thread_data);
static void *VDiff_PT(void *thread_data);
static void *VNeg_PT(void *thread_data);
static void *VScaleSum_PT(void *thread_data);
static void *VScaleDiff_PT(void *thread_data);
static void *VLin1_PT(void *thread_data);
static void *VLin2_PT(void *thread_data);
static void *VScaleBy_PT(void *thread_data);
static void *Vaxpy_PT(void *thread_data);

/* Pthread companion functions for fused vector operations */
static void *N_VLinearCombination_PT(void *thread_data);
static void *N_VScaleAddMulti_PT(void *thread_data);
static void *N_VDotProdMulti_PT(void *thread_data);

/* Pthread companion functions for vector array operations */
static void *N_VLinearSumVectorArray_PT(void *thread_data);
static void *N_VScaleVectorArray_PT(void *thread_data);
static void *N_VConstVectorArray_PT(void *thread_data);
static void *N_VWrmsNormVectorArray_PT(void *thread_data);
static void *N_VWrmsNormMaskVectorArray_PT(void *thread_data);
static void *N_VScaleAddMultiVectorArray_PT(void *thread_data);
static void *N_VLinearCombinationVectorArray_PT(void *thread_data);

/* Pthread companion functions special cases of vector array operations */
static void *VSumVectorArray_PT(void *thread_data);
static void *VDiffVectorArray_PT(void *thread_data);
static void *VScaleSumVectorArray_PT(void *thread_data);
static void *VScaleDiffVectorArray_PT(void *thread_data);
static void *VLin1VectorArray_PT(void *thread_data);
static void *VLin2VectorArray_PT(void *thread_data);
static void *VaxpyVectorArray_PT(void *thread_data);

/* Pthread companion functions for XBraid interface operations */
static void *VBufPack_PT(void *thread_data);
static void *VBufUnpack_PT(void *thread_data);

/* Function to determine loop values for threads */
static void N_VSplitLoop(int myid, int *nthreads, sunindextype *N,
                         sunindextype *start, sunindextype *end);

/* Function to initialize thread data */
static void N_VInitThreadData(Pthreads_Data *thread_data);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Pthreads(N_Vector v)
{
  return SUNDIALS_NVEC_PTHREADS;
}

/* ----------------------------------------------------------------------------
 * Function to create a new empty vector
 */

N_Vector N_VNewEmpty_Pthreads(sunindextype length, int num_threads,
                              SUNContext sunctx)
{
  N_Vector v;
  N_VectorContent_Pthreads content;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Pthreads;
  v->ops->nvclone           = N_VClone_Pthreads;
  v->ops->nvcloneempty      = N_VCloneEmpty_Pthreads;
  v->ops->nvdestroy         = N_VDestroy_Pthreads;
  v->ops->nvspace           = N_VSpace_Pthreads;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_Pthreads;
  v->ops->nvsetarraypointer = N_VSetArrayPointer_Pthreads;
  v->ops->nvgetlength       = N_VGetLength_Pthreads;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Pthreads;
  v->ops->nvconst        = N_VConst_Pthreads;
  v->ops->nvprod         = N_VProd_Pthreads;
  v->ops->nvdiv          = N_VDiv_Pthreads;
  v->ops->nvscale        = N_VScale_Pthreads;
  v->ops->nvabs          = N_VAbs_Pthreads;
  v->ops->nvinv          = N_VInv_Pthreads;
  v->ops->nvaddconst     = N_VAddConst_Pthreads;
  v->ops->nvdotprod      = N_VDotProd_Pthreads;
  v->ops->nvmaxnorm      = N_VMaxNorm_Pthreads;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Pthreads;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Pthreads;
  v->ops->nvmin          = N_VMin_Pthreads;
  v->ops->nvwl2norm      = N_VWL2Norm_Pthreads;
  v->ops->nvl1norm       = N_VL1Norm_Pthreads;
  v->ops->nvcompare      = N_VCompare_Pthreads;
  v->ops->nvinvtest      = N_VInvTest_Pthreads;
  v->ops->nvconstrmask   = N_VConstrMask_Pthreads;
  v->ops->nvminquotient  = N_VMinQuotient_Pthreads;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProd_Pthreads;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Pthreads;
  v->ops->nvminlocal         = N_VMin_Pthreads;
  v->ops->nvl1normlocal      = N_VL1Norm_Pthreads;
  v->ops->nvinvtestlocal     = N_VInvTest_Pthreads;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Pthreads;
  v->ops->nvminquotientlocal = N_VMinQuotient_Pthreads;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Pthreads;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Pthreads;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = N_VDotProdMulti_Pthreads;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_Pthreads;
  v->ops->nvbufpack   = N_VBufPack_Pthreads;
  v->ops->nvbufunpack = N_VBufUnpack_Pthreads;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Pthreads) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length      = length;
  content->num_threads = num_threads;
  content->own_data    = SUNFALSE;
  content->data        = NULL;

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new vector
 */

N_Vector N_VNew_Pthreads(sunindextype length, int num_threads,
                         SUNContext sunctx)
{
  N_Vector v;
  realtype *data;

  v = NULL;
  v = N_VNewEmpty_Pthreads(length, num_threads, sunctx);
  if (v == NULL) return(NULL);

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_Pthreads(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_PT(v) = SUNTRUE;
    NV_DATA_PT(v)     = data;

  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a vector with user data component
 */

N_Vector N_VMake_Pthreads(sunindextype length, int num_threads,
                          realtype *v_data, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Pthreads(length, num_threads, sunctx);
  if (v == NULL) return(NULL);

  if (length > 0) {
    /* Attach data */
    NV_OWN_DATA_PT(v) = SUNFALSE;
    NV_DATA_PT(v)     = v_data;
  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new vectors.
 */

N_Vector* N_VCloneVectorArray_Pthreads(int count, N_Vector w)
{
  return(N_VCloneVectorArray(count, w));
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new vectors with NULL data array.
 */

N_Vector* N_VCloneVectorArrayEmpty_Pthreads(int count, N_Vector w)
{
  return(N_VCloneEmptyVectorArray(count, w));
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Pthreads
 */

void N_VDestroyVectorArray_Pthreads(N_Vector* vs, int count)
{
  N_VDestroyVectorArray(vs, count);
  return;
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_Pthreads(N_Vector v)
{
  return NV_LENGTH_PT(v);
}

/* ----------------------------------------------------------------------------
 * Function to print a vector to stdout
 */

void N_VPrint_Pthreads(N_Vector x)
{
  N_VPrintFile_Pthreads(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print a vector to outfile
 */

void N_VPrintFile_Pthreads(N_Vector x, FILE *outfile)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LENGTH_PT(x);
  xd = NV_DATA_PT(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    STAN_SUNDIALS_FPRINTF(outfile, "%11.8Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    STAN_SUNDIALS_FPRINTF(outfile, "%11.8g\n", xd[i]);
#else
    STAN_SUNDIALS_FPRINTF(outfile, "%11.8g\n", xd[i]);
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

/* ----------------------------------------------------------------------------
 * Create new vector from existing vector without attaching data
 */

N_Vector N_VCloneEmpty_Pthreads(N_Vector w)
{
  N_Vector v;
  N_VectorContent_Pthreads content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Pthreads) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length      = NV_LENGTH_PT(w);
  content->num_threads = NV_NUM_THREADS_PT(w);
  content->own_data    = SUNFALSE;
  content->data        = NULL;

  return(v);
}


/* ----------------------------------------------------------------------------
 * Create new vector from existing vector and attach data
 */

N_Vector N_VClone_Pthreads(N_Vector w)
{
  N_Vector v;
  realtype *data;
  sunindextype length;

  v = NULL;
  v = N_VCloneEmpty_Pthreads(w);
  if (v == NULL) return(NULL);

  length = NV_LENGTH_PT(w);

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_Pthreads(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_PT(v) = SUNTRUE;
    NV_DATA_PT(v)     = data;

  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Destroy vector and free vector memory
 */

void N_VDestroy_Pthreads(N_Vector v)
{
  if (v == NULL) return;

  /* free content */
  if (v->content != NULL) {
    if (NV_OWN_DATA_PT(v) && NV_DATA_PT(v) != NULL) {
      free(NV_DATA_PT(v));
      NV_DATA_PT(v) = NULL;
    }
    free(v->content);
    v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL) { free(v->ops); v->ops = NULL; }
  free(v); v = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Get storage requirement for vector
 */

void N_VSpace_Pthreads(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NV_LENGTH_PT(v);
  *liw = 1;

  return;
}


/* ----------------------------------------------------------------------------
 * Get vector data pointer
 */

realtype *N_VGetArrayPointer_Pthreads(N_Vector v)
{
  return((realtype *) NV_DATA_PT(v));
}


/* ----------------------------------------------------------------------------
 * Set vector data pointer
 */

void N_VSetArrayPointer_Pthreads(realtype *v_data, N_Vector v)
{
  if (NV_LENGTH_PT(v) > 0) NV_DATA_PT(v) = v_data;

  return;
}


/* ----------------------------------------------------------------------------
 * Compute linear sum z[i] = a*x[i]+b*y[i]
 */

void N_VLinearSum_Pthreads(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  realtype c;
  N_Vector v1, v2;
  booleantype test;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Pthreads(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy_Pthreads(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_Pthreads(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Pthreads(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Pthreads(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Pthreads(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_Pthreads(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_Pthreads(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = (pthread_t *) malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].c2 = b;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VLinearSum_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VLinearSum
 */

static void *N_VLinearSum_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype a, b;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  a  = my_data->c1;
  b  = my_data->c2;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute linear sum */
  for (i = start; i < end; i++){
    zd[i] = (a*xd[i])+(b*yd[i]);
  }

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Assigns constant value to all vector elements, z[i] = c
 */

void N_VConst_Pthreads(realtype c, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(z);
  nthreads     = NV_NUM_THREADS_PT(z);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = c;
    thread_data[i].v1 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VConst_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VConst
 */

static void *N_VConst_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype c;
  realtype *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  c  = my_data->c1;
  zd = my_data->v1;

  start = my_data->start;
  end   = my_data->end;

  /* assign constant values */
  for (i = start; i < end; i++)
    zd[i] = c;

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute componentwise product z[i] = x[i]*y[i]
 */

void N_VProd_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VProd_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and exit */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VProd
 */

static void *N_VProd_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute componentwise product */
  for (i = start; i < end; i++)
    zd[i] = xd[i]*yd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute componentwise division z[i] = x[i]/y[i]
 */

void N_VDiv_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VDiv_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VDiv
 */

static void *N_VDiv_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute componentwise division */
  for (i = start; i < end; i++)
    zd[i] = xd[i]/yd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute scaler multiplication z[i] = c*x[i]
 */

void N_VScale_Pthreads(realtype c, N_Vector x, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  if (z == x) {  /* BLAS usage: scale x <- cx */
    VScaleBy_Pthreads(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_Pthreads(x, z);
  } else if (c == -ONE) {
    VNeg_Pthreads(x, z);
  } else {
    /* allocate threads and thread data structs */
    N            = NV_LENGTH_PT(x);
    nthreads     = NV_NUM_THREADS_PT(x);
    threads      = malloc(nthreads*sizeof(pthread_t));
    thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

    /* set thread attributes */
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (i=0; i<nthreads; i++) {
      /* initialize thread data */
      N_VInitThreadData(&thread_data[i]);

      /* compute start and end loop index for thread */
      N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

      /* pack thread data */
      thread_data[i].c1 = c;
      thread_data[i].v1 = NV_DATA_PT(x);
      thread_data[i].v2 = NV_DATA_PT(z);

      /* create threads and call pthread companion function */
      pthread_create(&threads[i], &attr, N_VScale_PT, (void *) &thread_data[i]);
    }

    /* wait for all threads to finish */
    for (i=0; i<nthreads; i++) {
      pthread_join(threads[i], NULL);
    }

    /* clean up */
    pthread_attr_destroy(&attr);
    free(threads);
    free(thread_data);
  }

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VScale
 */

static void *N_VScale_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype c;
  realtype *xd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  c  = my_data->c1;
  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute scaler multiplication */
  for (i = start; i < end; i++)
    zd[i] = c*xd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute absolute value of vector components z[i] = SUNRabs(x[i])
 */

void N_VAbs_Pthreads(N_Vector x, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VAbs_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VAbs
 */

static void *N_VAbs_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute absolute value of components */
  for (i = start; i < end; i++)
    zd[i] = SUNRabs(xd[i]);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = 1 / x[i]
 */

void N_VInv_Pthreads(N_Vector x, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VInv_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VInv
 */

static void *N_VInv_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute componentwise inverse */
  for (i = start; i < end; i++)
    zd[i] = ONE/xd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute componentwise addition of a scaler to a vector z[i] = x[i] + b
 */

void N_VAddConst_Pthreads(N_Vector x, realtype b, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = b;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VAddConst_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VAddConst
 */

static void *N_VAddConst_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype b;
  realtype *xd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  b  = my_data->c1;
  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute componentwise constant addition */
  for (i = start; i < end; i++)
    zd[i] = xd[i] + b;

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Computes the dot product of two vectors, a = sum(x[i]*y[i])
 */

realtype N_VDotProd_Pthreads(N_Vector x, N_Vector y)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;
  realtype        sum = ZERO;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VDotProd_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VDotProd
 */

static void *N_VDotProd_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *yd;
  realtype local_sum, *global_sum;
  Pthreads_Data *my_data;
  pthread_mutex_t *global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  yd = my_data->v2;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute dot product */
  local_sum = ZERO;
  for (i = start; i < end; i++)
    local_sum += xd[i] * yd[i];

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Computes max norm of the vector
 */

realtype N_VMaxNorm_Pthreads(N_Vector x)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;
  realtype        max = ZERO;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].global_val   = &max;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VMaxNorm_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(max);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VMaxNorm
 */

static void *N_VMaxNorm_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd;
  realtype local_max, *global_max;
  Pthreads_Data *my_data;
  pthread_mutex_t *global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;

  global_max   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* find local max */
  local_max = ZERO;
  for (i = start; i < end; i++)
    if (SUNRabs(xd[i]) > local_max) local_max = SUNRabs(xd[i]);

  /* update global max */
  pthread_mutex_lock(global_mutex);
  if (local_max > *global_max) {
    *global_max = local_max;
  }
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a vector
 */

realtype N_VWrmsNorm_Pthreads(N_Vector x, N_Vector w)
{
  return(SUNRsqrt(N_VWSqrSumLocal_Pthreads(x, w)/(NV_LENGTH_PT(x))));
}


/* ----------------------------------------------------------------------------
 * Computes weighted square sum of a vector
 */

realtype N_VWSqrSumLocal_Pthreads(N_Vector x, N_Vector w)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;
  realtype        sum = ZERO;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(w);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VWSqrSum_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWrmsNorm
 */

static void *N_VWSqrSum_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *wd;
  realtype local_sum, *global_sum;
  Pthreads_Data *my_data;
  pthread_mutex_t *global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  wd = my_data->v2;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute wrms norm */
  local_sum = ZERO;
  for (i = start; i < end; i++)
    local_sum += SUNSQR(xd[i] * wd[i]);

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a masked vector
 */

realtype N_VWrmsNormMask_Pthreads(N_Vector x, N_Vector w, N_Vector id)
{
  return(SUNRsqrt(N_VWSqrSumMaskLocal_Pthreads(x, w, id)/(NV_LENGTH_PT(x))));
}


/* ----------------------------------------------------------------------------
 * Computes weighted square sum of a masked vector
 */

realtype N_VWSqrSumMaskLocal_Pthreads(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;
  realtype        sum = ZERO;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(w);
    thread_data[i].v3 = NV_DATA_PT(id);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VWSqrSumMask_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWSqrSumMask
 */

static void *N_VWSqrSumMask_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *wd, *idd;
  realtype local_sum, *global_sum;
  Pthreads_Data *my_data;
  pthread_mutex_t *global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd  = my_data->v1;
  wd  = my_data->v2;
  idd = my_data->v3;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute wrms norm with mask */
  local_sum = ZERO;
  for (i = start; i < end; i++) {
    if (idd[i] > ZERO)
      local_sum += SUNSQR(xd[i]*wd[i]);
  }

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Finds the minimun component of a vector
 */

realtype N_VMin_Pthreads(N_Vector x)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;
  realtype        min;

  /* initialize global min */
  min = NV_Ith_PT(x,0);

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].global_val   = &min;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VMin_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(min);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VMin
 */

static void *N_VMin_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd;
  realtype local_min, *global_min;
  Pthreads_Data *my_data;
  pthread_mutex_t *global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;

  global_min   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* find local min */
  local_min = *global_min;
  for (i = start; i < end; i++) {
    if (xd[i] < local_min)
      local_min = xd[i];
  }

  /* update global min */
  pthread_mutex_lock(global_mutex);
  if (local_min < *global_min)
    *global_min = local_min;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Computes weighted L2 norm of a vector
 */

realtype N_VWL2Norm_Pthreads(N_Vector x, N_Vector w)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;
  realtype        sum = ZERO;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(w);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VWL2Norm_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(SUNRsqrt(sum));
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWL2Norm
 */

static void *N_VWL2Norm_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *wd;
  realtype local_sum, *global_sum;
  Pthreads_Data *my_data;
  pthread_mutex_t *global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  wd = my_data->v2;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute WL2 norm */
  local_sum = ZERO;
  for (i = start; i < end; i++)
    local_sum += SUNSQR(xd[i]*wd[i]);

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Computes L1 norm of a vector
 */

realtype N_VL1Norm_Pthreads(N_Vector x)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;
  realtype        sum = ZERO;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VL1Norm_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VL1Norm
 */

static void *N_VL1Norm_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd;
  realtype local_sum, *global_sum;
  Pthreads_Data *my_data;
  pthread_mutex_t *global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;

  global_sum   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute L1 norm */
  local_sum = ZERO;
  for (i = start; i < end; i++)
    local_sum += SUNRabs(xd[i]);

  /* update global sum */
  pthread_mutex_lock(global_mutex);
  *global_sum += local_sum;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compare vector component values to a scaler
 */

void N_VCompare_Pthreads(realtype c, N_Vector x, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1  = c;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VCompare_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VCompare
 */

static void *N_VCompare_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype c;
  realtype *xd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  c  = my_data->c1;
  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compare component to scaler */
  for (i = start; i < end; i++)
    zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO;

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = ONE/x[i] and check if x[i] == ZERO
 */

booleantype N_VInvTest_Pthreads(N_Vector x, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  realtype val = ZERO;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);
    thread_data[i].global_val = &val;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VInvTest_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  if (val > ZERO)
    return (SUNFALSE);
  else
    return (SUNTRUE);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VInvTest
 */

static void *N_VInvTest_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *zd;
  realtype local_val, *global_val;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  global_val = my_data->global_val;

  start = my_data->start;
  end   = my_data->end;

  /* compute inverse with check for divide by ZERO */
  local_val = ZERO;
  for (i = start; i < end; i++) {
    if (xd[i] == ZERO)
      local_val = ONE;
    else
      zd[i] = ONE/xd[i];
  }

  /* update global val */
  if (local_val > ZERO) {
    *global_val = local_val;
  }

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute constraint mask of a vector
 */

booleantype N_VConstrMask_Pthreads(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  realtype val = ZERO;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(c);
    thread_data[i].v2 = NV_DATA_PT(x);
    thread_data[i].v3 = NV_DATA_PT(m);
    thread_data[i].global_val = &val;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VConstrMask_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  if (val > ZERO)
    return(SUNFALSE);
  else
    return(SUNTRUE);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VConstrMask
 */

static void *N_VConstrMask_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *cd, *xd, *md;
  realtype local_val, *global_val;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  cd = my_data->v1;
  xd = my_data->v2;
  md = my_data->v3;

  global_val = my_data->global_val;

  start = my_data->start;
  end   = my_data->end;

  /* compute constraint mask */
  local_val = ZERO;
  for (i = start; i < end; i++) {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO)
      continue;

    /* Check if a set constraint has been violated */
    if ((SUNRabs(cd[i]) > ONEPT5 && xd[i]*cd[i] <= ZERO) ||
        (SUNRabs(cd[i]) > HALF   && xd[i]*cd[i] <  ZERO)) {
      local_val = md[i] = ONE;
    }
  }

  /* update global val */
  if (local_val > ZERO) {
    *global_val = local_val;
  }

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute minimum componentwise quotient
 */

realtype N_VMinQuotient_Pthreads(N_Vector num, N_Vector denom)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;
  realtype        min = BIG_REAL;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(num);
  nthreads    = NV_NUM_THREADS_PT(num);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(num);
    thread_data[i].v2 = NV_DATA_PT(denom);
    thread_data[i].global_val   = &min;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VMinQuotient_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(min);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VConstrMask
 */

static void *N_VMinQuotient_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *nd, *dd;
  realtype local_min, *global_min;
  Pthreads_Data *my_data;
  pthread_mutex_t *global_mutex;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  nd = my_data->v1;
  dd = my_data->v2;

  global_min   = my_data->global_val;
  global_mutex = my_data->global_mutex;

  start = my_data->start;
  end   = my_data->end;

  /* compute minimum quotient */
  local_min = BIG_REAL;
  for (i = start; i < end; i++) {
    if (dd[i] == ZERO)
      continue;
    local_min = SUNMIN(local_min, nd[i]/dd[i]);
  }

  /* update global min */
  pthread_mutex_lock(global_mutex);
  if (local_min < *global_min)
    *global_min = local_min;
  pthread_mutex_unlock(global_mutex);

  /* exit */
  pthread_exit(NULL);
}


/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */


/* -----------------------------------------------------------------------------
 * Compute the linear combination z = c[i]*X[i]
 */

int N_VLinearCombination_Pthreads(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Pthreads(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_Pthreads(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  /* get vector length and data array */
  N           = NV_LENGTH_PT(z);
  nthreads    = NV_NUM_THREADS_PT(z);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].cvals = c;
    thread_data[i].Y1    = X;
    thread_data[i].x1    = z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VLinearCombination_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VLinearCombination
 */

static void *N_VLinearCombination_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype* c=NULL;
  realtype* xd=NULL;
  realtype* zd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;

  c  = my_data->cvals;
  zd = NV_DATA_PT(my_data->x1);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((my_data->Y1[0] == my_data->x1) && (c[0] == ONE)) {
    for (i=1; i<my_data->nvec; i++) {
      xd = NV_DATA_PT(my_data->Y1[i]);
      for (j=start; j<end; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    pthread_exit(NULL);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (my_data->Y1[0] == my_data->x1) {
    for (j=start; j<end; j++) {
      zd[j] *= c[0];
    }
    for (i=1; i<my_data->nvec; i++) {
      xd = NV_DATA_PT(my_data->Y1[i]);
      for (j=start; j<end; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    pthread_exit(NULL);
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = NV_DATA_PT(my_data->Y1[0]);
  for (j=start; j<end; j++) {
    zd[j] = c[0] * xd[j];
  }
  for (i=1; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    for (j=start; j<end; j++) {
      zd[j] += c[i] * xd[j];
    }
  }
  pthread_exit(NULL);
}


/* -----------------------------------------------------------------------------
 * Compute multiple linear sums Z[i] = Y[i] + a*x
 */

int N_VScaleAddMulti_Pthreads(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Pthreads(a[0], x, ONE, Y[0], Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].cvals = a;
    thread_data[i].x1    = x;
    thread_data[i].Y1    = Y;
    thread_data[i].Y2    = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VScaleAddMulti_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VScaleAddMulti
 */

static void *N_VScaleAddMulti_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype* a=NULL;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;

  a  = my_data->cvals;
  xd = NV_DATA_PT(my_data->x1);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (my_data->Y1 == my_data->Y2) {
    for (i=0; i<my_data->nvec; i++) {
      yd = NV_DATA_PT(my_data->Y1[i]);
      for (j=start; j<end; j++) {
        yd[j] += a[i] * xd[j];
      }
    }
    pthread_exit(NULL);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<my_data->nvec; i++) {
    yd = NV_DATA_PT(my_data->Y1[i]);
    zd = NV_DATA_PT(my_data->Y2[i]);
    for (j=start; j<end; j++) {
      zd[j] = a[i] * xd[j] + yd[j];
    }
  }
  pthread_exit(NULL);
}


/* -----------------------------------------------------------------------------
 * Compute the dot product of a vector with multiple vectors, a[i] = sum(x*Y[i])
 */

int N_VDotProdMulti_Pthreads(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VDotProd */
  if (nvec == 1) {
    dotprods[0] = N_VDotProd_Pthreads(x, Y[0]);
    return(0);
  }

  /* initialize output array */
  for (i=0; i<nvec; i++)
    dotprods[i] = ZERO;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].x1    = x;
    thread_data[i].Y1    = Y;
    thread_data[i].cvals = dotprods;

    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VDotProdMulti_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VDotProdMulti
 */

static void *N_VDotProdMulti_PT(void *thread_data)
{
  Pthreads_Data   *my_data;
  sunindextype    j, start, end;
  pthread_mutex_t *lock;

  int       i;
  realtype  sum;
  realtype* dotprods=NULL;
  realtype* xd=NULL;
  realtype* yd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;
  lock  = my_data->global_mutex;

  xd       = NV_DATA_PT(my_data->x1);
  dotprods = my_data->cvals;

  /* compute multiple dot products */
  for (i=0; i<my_data->nvec; i++) {
    yd  = NV_DATA_PT(my_data->Y1[i]);
    sum = ZERO;
    for (j=start; j<end; j++) {
      sum += xd[j] * yd[j];
    }
    /* update global sum */
    pthread_mutex_lock(lock);
    dotprods[i] += sum;
    pthread_mutex_unlock(lock);
  }

  /* exit */
  pthread_exit(NULL);
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */


/* -----------------------------------------------------------------------------
 * Compute multiple linear sums Z[i] = a*X[i] + b*Y[i]
 */

int N_VLinearSumVectorArray_Pthreads(int nvec, realtype a, N_Vector* X,
                                     realtype b, N_Vector* Y,  N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  realtype    c;
  N_Vector*  V1;
  N_Vector*  V2;
  booleantype test;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Pthreads(a, X[0], b, Y[0], Z[0]);
    return(0);
  }

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (Z == Y))
    return(VaxpyVectorArray_Pthreads(nvec, a, X, Y));

  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (Z == X))
    return(VaxpyVectorArray_Pthreads(nvec, b, Y, X));

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE))
    return(VSumVectorArray_Pthreads(nvec, X, Y, Z));

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VDiffVectorArray_Pthreads(nvec, V2, V1, Z));
  }

  /* Cases:                                                  */
  /*   (1) a == 1.0, b == other or 0.0,                      */
  /*   (2) a == other or 0.0, b == 1.0                       */
  /* if a or b is 0.0, then user should have called N_VScale */
  if ((test = (a == ONE)) || (b == ONE)) {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VLin1VectorArray_Pthreads(nvec, c, V1, V2, Z));
  }

  /* Cases:                     */
  /*   (1) a == -1.0, b != 1.0, */
  /*   (2) a != 1.0, b == -1.0  */
  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VLin2VectorArray_Pthreads(nvec, c, V1, V2, Z));
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
    return(VScaleSumVectorArray_Pthreads(nvec, a, X, Y, Z));

  /* Case: a == -b */
  if (a == -b)
    return(VScaleDiffVectorArray_Pthreads(nvec, a, X, Y, Z));

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length and data array */
  N           = NV_LENGTH_PT(Z[0]);
  nthreads    = NV_NUM_THREADS_PT(Z[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec = nvec;
    thread_data[i].c1   = a;
    thread_data[i].c2   = b;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VLinearSumVectorArray_PT,
                   (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VLinearSumVectorArray
 */

static void *N_VLinearSumVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype  a, b;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;

  a = my_data->c1;
  b = my_data->c2;

  /* compute linear sum for each vector pair in vector arrays */
  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j=start; j<end; j++) {
      zd[j] = a * xd[j] + b * yd[j];
    }
  }

  /* exit */
  pthread_exit(NULL);
}


/* -----------------------------------------------------------------------------
 * Scale multiple vectors Z[i] = c[i]*X[i]
 */

int N_VScaleVectorArray_Pthreads(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Pthreads(c[0], X[0], Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N           = NV_LENGTH_PT(Z[0]);
  nthreads    = NV_NUM_THREADS_PT(Z[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].cvals = c;
    thread_data[i].Y1    = X;
    thread_data[i].Y2    = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VScaleVectorArray_PT,
                   (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VScaleVectorArray
 */

static void *N_VScaleVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype* c;
  realtype* xd=NULL;
  realtype* zd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;

  c = my_data->cvals;

  /*
   * X[i] *= c[i]
   */
  if (my_data->Y1 == my_data->Y2) {
    for (i=0; i<my_data->nvec; i++) {
      xd = NV_DATA_PT(my_data->Y1[i]);
      for (j=start; j<end; j++) {
        xd[j] *= c[i];
      }
    }
    pthread_exit(NULL);
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    zd = NV_DATA_PT(my_data->Y2[i]);
    for (j=start; j<end; j++) {
      zd[j] = c[i] * xd[j];
    }
  }
  pthread_exit(NULL);
}


/* -----------------------------------------------------------------------------
 * Set multiple vectors to a constant value Z[i] = c
 */

int N_VConstVectorArray_Pthreads(int nvec, realtype c, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VConst */
  if (nvec == 1) {
    N_VConst_Pthreads(c, Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N           = NV_LENGTH_PT(Z[0]);
  nthreads    = NV_NUM_THREADS_PT(Z[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec = nvec;
    thread_data[i].c1   = c;
    thread_data[i].Y1   = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VConstVectorArray_PT,
                   (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VConstVectorArray
 */

static void *N_VConstVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype* zd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;

  /* set each vector in the vector array to a constant */
  for (i=0; i<my_data->nvec; i++) {
    zd = NV_DATA_PT(my_data->Y1[i]);
    for (j=start; j<end; j++) {
      zd[j] = my_data->c1;
    }
  }

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute the weighted root mean square norm of multiple vectors
 */

int N_VWrmsNormVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* W, realtype* nrm)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNorm_Pthreads(X[0], W[0]);
    return(0);
  }

  /* initialize output array */
  for (i=0; i<nvec; i++)
    nrm[i] = ZERO;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].Y1    = X;
    thread_data[i].Y2    = W;
    thread_data[i].cvals = nrm;

    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VWrmsNormVectorArray_PT,
                   (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* finalize wrms calculation */
  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/N);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(0);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWrmsNormVectorArray
 */

static void *N_VWrmsNormVectorArray_PT(void *thread_data)
{
  Pthreads_Data   *my_data;
  sunindextype    j, start, end;
  pthread_mutex_t *lock;

  int       i;
  realtype  sum;
  realtype* nrm=NULL;
  realtype* xd=NULL;
  realtype* wd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;
  lock  = my_data->global_mutex;

  nrm = my_data->cvals;

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    wd = NV_DATA_PT(my_data->Y2[i]);
    sum = ZERO;
    for (j=start; j<end; j++) {
      sum += SUNSQR(xd[j] * wd[j]);
    }
    /* update global sum */
    pthread_mutex_lock(lock);
    nrm[i] += sum;
    pthread_mutex_unlock(lock);
  }

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute the weighted root mean square norm of multiple vectors
 */

int N_VWrmsNormMaskVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* W,
                                        N_Vector id, realtype* nrm)
{
  sunindextype    N;
  int             i, nthreads;
  pthread_t       *threads;
  Pthreads_Data   *thread_data;
  pthread_attr_t  attr;
  pthread_mutex_t global_mutex;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNormMask_Pthreads(X[0], W[0], id);
    return(0);
  }

  /* initialize output array */
  for (i=0; i<nvec; i++)
    nrm[i] = ZERO;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* lock for reduction */
  pthread_mutex_init(&global_mutex, NULL);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].Y1    = X;
    thread_data[i].Y2    = W;
    thread_data[i].x1    = id;
    thread_data[i].cvals = nrm;

    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VWrmsNormMaskVectorArray_PT,
                   (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* finalize wrms calculation */
  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/N);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&global_mutex);
  free(threads);
  free(thread_data);

  return(0);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWrmsNormVectorArray
 */

static void *N_VWrmsNormMaskVectorArray_PT(void *thread_data)
{
  Pthreads_Data   *my_data;
  sunindextype    j, start, end;
  pthread_mutex_t *lock;

  int       i;
  realtype  sum;
  realtype* nrm=NULL;
  realtype* xd=NULL;
  realtype* wd=NULL;
  realtype* idd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;
  lock  = my_data->global_mutex;

  nrm = my_data->cvals;
  idd = NV_DATA_PT(my_data->x1);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    wd = NV_DATA_PT(my_data->Y2[i]);
    sum = ZERO;
    for (j=start; j<end; j++) {
      if (idd[j] > ZERO)
        sum += SUNSQR(xd[j] * wd[j]);
    }
    /* update global sum */
    pthread_mutex_lock(lock);
    nrm[i] += sum;
    pthread_mutex_unlock(lock);
  }

  /* exit */
  pthread_exit(NULL);
}


/* -----------------------------------------------------------------------------
 * Scale and add a vector to multiple vectors Z = Y + a*X
 */

int N_VScaleAddMultiVectorArray_Pthreads(int nvec, int nsum, realtype* a,
                                          N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  sunindextype   N;
  int            i, j, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  int          retval;
  N_Vector*   YY;
  N_Vector*   ZZ;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VLinearSum */
    if (nsum == 1) {
      N_VLinearSum_Pthreads(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      return(0);
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector*) malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector*) malloc(nsum * sizeof(N_Vector));

    for (j=0; j<nsum; j++) {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    retval = N_VScaleAddMulti_Pthreads(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);
    return(retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1) {
    retval = N_VLinearSumVectorArray_Pthreads(nvec, a[0], X, ONE, Y[0], Z[0]);
    return(retval);
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length and data array */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].nsum  = nsum;
    thread_data[i].cvals = a;
    thread_data[i].Y1    = X;
    thread_data[i].ZZ1   = Y;
    thread_data[i].ZZ2   = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VScaleAddMultiVectorArray_PT,
                   (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VScaleAddMultiVectorArray
 */

static void *N_VScaleAddMultiVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  k, start, end;

  int       i, j;
  realtype* a=NULL;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;

  a = my_data->cvals;

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (my_data->ZZ1 == my_data->ZZ2) {
    for (i=0; i<my_data->nvec; i++) {
      xd = NV_DATA_PT(my_data->Y1[i]);
      for (j=0; j<my_data->nsum; j++){
        yd = NV_DATA_PT(my_data->ZZ1[j][i]);
        for (k=start; k<end; k++) {
          yd[k] += a[j] * xd[k];
        }
      }
    }
    pthread_exit(NULL);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    for (j=0; j<my_data->nsum; j++){
      yd = NV_DATA_PT(my_data->ZZ1[j][i]);
      zd = NV_DATA_PT(my_data->ZZ2[j][i]);
      for (k=start; k<end; k++) {
        zd[k] = a[j] * xd[k] + yd[k];
      }
    }
  }
  pthread_exit(NULL);
}


/* -----------------------------------------------------------------------------
 * Compute a linear combination for multiple vectors
 */

int N_VLinearCombinationVectorArray_Pthreads(int nvec, int nsum, realtype* c,
                                           N_Vector** X, N_Vector* Z)
{
  sunindextype   N;
  int            i, j, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  int          retval;
  realtype*    ctmp;
  N_Vector*   Y;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VScale */
    if (nsum == 1) {
      N_VScale_Pthreads(c[0], X[0][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearSum */
    if (nsum == 2) {
      N_VLinearSum_Pthreads(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nsum; i++) {
      Y[i] = X[i][0];
    }

    retval = N_VLinearCombination_Pthreads(nsum, c, Y, Z[0]);

    free(Y);
    return(retval);
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

    retval = N_VScaleVectorArray_Pthreads(nvec, ctmp, X[0], Z);

    free(ctmp);
    return(retval);
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2) {
    retval = N_VLinearSumVectorArray_Pthreads(nvec, c[0], X[0], c[1], X[1], Z);
    return(retval);
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length and data array */
  N           = NV_LENGTH_PT(Z[0]);
  nthreads    = NV_NUM_THREADS_PT(Z[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].nvec  = nvec;
    thread_data[i].nsum  = nsum;
    thread_data[i].cvals = c;
    thread_data[i].ZZ1   = X;
    thread_data[i].Y1    = Z;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VLinearCombinationVectorArray_PT,
                   (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VLinearCombinationVectorArray
 */

static void *N_VLinearCombinationVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  k, start, end;

  int       i; /* vector arrays index in summation [0,nsum) */
  int       j; /* vector index in vector array     [0,nvec) */
  realtype* c=NULL;
  realtype* zd=NULL;
  realtype* xd=NULL;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  start = my_data->start;
  end   = my_data->end;

  c = my_data->cvals;

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((my_data->ZZ1[0] == my_data->Y1) && (c[0] == ONE)) {
    for (j=0; j<my_data->nvec; j++) {
      zd = NV_DATA_PT(my_data->Y1[j]);
      for (i=1; i<my_data->nsum; i++) {
        xd = NV_DATA_PT(my_data->ZZ1[i][j]);
        for (k=start; k<end; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    pthread_exit(NULL);
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (my_data->ZZ1[0] == my_data->Y1) {
    for (j=0; j<my_data->nvec; j++) {
      zd = NV_DATA_PT(my_data->Y1[j]);
      for (k=start; k<end; k++) {
        zd[k] *= c[0];
      }
      for (i=1; i<my_data->nsum; i++) {
        xd = NV_DATA_PT(my_data->ZZ1[i][j]);
        for (k=start; k<end; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    pthread_exit(NULL);
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j=0; j<my_data->nvec; j++) {
    xd = NV_DATA_PT(my_data->ZZ1[0][j]);
    zd = NV_DATA_PT(my_data->Y1[j]);
    for (k=start; k<end; k++) {
      zd[k] = c[0] * xd[k];
    }
    for (i=1; i<my_data->nsum; i++) {
      xd = NV_DATA_PT(my_data->ZZ1[i][j]);
      for (k=start; k<end; k++) {
        zd[k] += c[i] * xd[k];
      }
    }
  }
  pthread_exit(NULL);
}


/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */



/* -----------------------------------------------------------------------------
 * Set buffer size
 */

int N_VBufSize_Pthreads(N_Vector x, sunindextype *size)
{
  if (x == NULL) return(-1);
  *size = NV_LENGTH_PT(x) * ((sunindextype)sizeof(realtype));
  return(0);
}


/* -----------------------------------------------------------------------------
 * Pack butter
 */

int N_VBufPack_Pthreads(N_Vector x, void *buf)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  if (x == NULL || buf == NULL) return(-1);

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = (realtype*)buf;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VBufPack_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VBufPack
 */


static void *VBufPack_PT(void *thread_data)
{
  sunindextype  i, start, end;
  realtype      *xd, *bd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  bd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* pack the buffer */
  for (i = start; i < end; i++)
    bd[i] = xd[i];

  /* exit */
  pthread_exit(NULL);
}


/* -----------------------------------------------------------------------------
 * Unpack butter
 */

int N_VBufUnpack_Pthreads(N_Vector x, void *buf)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  if (x == NULL || buf == NULL) return(-1);

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(x);
  nthreads    = NV_NUM_THREADS_PT(x);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = (realtype*)buf;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VBufUnpack_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}


/* -----------------------------------------------------------------------------
 * Pthread companion function to N_VBufUnpack
 */


static void *VBufUnpack_PT(void *thread_data)
{
  sunindextype  i, start, end;
  realtype      *xd, *bd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  bd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* unpack the buffer */
  for (i = start; i < end; i++)
    xd[i] = bd[i];

  /* exit */
  pthread_exit(NULL);
}


/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector operations
 * -----------------------------------------------------------------
 */


/* ----------------------------------------------------------------------------
 * Copy vector components into second vector
 */

static void VCopy_Pthreads(N_Vector x, N_Vector z)
{
  sunindextype      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VCopy_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to VCopy
 */

static void *VCopy_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* copy vector components */
  for (i = start; i < end; i++)
    zd[i] = xd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute vector sum
 */

static void VSum_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VSum_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to VSum
 */

static void *VSum_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute vector sum */
  for (i = start; i < end; i++)
    zd[i] = xd[i] + yd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute vector difference
 */

static void VDiff_Pthreads(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype  N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VDiff_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to VDiff
 */

static void *VDiff_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute vector difference */
  for (i = start; i < end; i++)
    zd[i] = xd[i] - yd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute the negative of a vector
 */

static void VNeg_Pthreads(N_Vector x, N_Vector z)
{
  sunindextype  N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VNeg_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to VNeg
 */

static void *VNeg_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype *xd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  xd = my_data->v1;
  zd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute negative of vector */
  for (i = start; i < end; i++)
    zd[i] = -xd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector sum
 */

static void VScaleSum_Pthreads(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype  N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = c;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VScaleSum_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to VScaleSum
 */

static void *VScaleSum_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype c;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  c  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute scaled vector sum */
  for (i = start; i < end; i++)
    zd[i] = c*(xd[i] + yd[i]);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector difference
 */

static void VScaleDiff_Pthreads(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype  N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = c;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VScaleDiff_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to VScaleDiff
 */

static void *VScaleDiff_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype c;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  c  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute scaled vector difference */
  for (i = start; i < end; i++)
    zd[i] = c*(xd[i] - yd[i]);

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute vector sum z[i] = a*x[i]+y[i]
 */

static void VLin1_Pthreads(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype  N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VLin1_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to VLin1
 */

static void *VLin1_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype a;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  a  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute vector sum */
  for (i = start; i < end; i++)
    zd[i] = (a*xd[i]) + yd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute vector difference z[i] = a*x[i]-y[i]
 */

static void VLin2_Pthreads(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype  N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);
    thread_data[i].v3 = NV_DATA_PT(z);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VLin2_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to VLin2
 */

static void *VLin2_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype a;
  realtype *xd, *yd, *zd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  a  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;
  zd = my_data->v3;

  start = my_data->start;
  end   = my_data->end;

  /* compute vector difference */
  for (i = start; i < end; i++)
    zd[i] = (a*xd[i]) - yd[i];

  /* exit */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute special cases of linear sum
 */

static void Vaxpy_Pthreads(realtype a, N_Vector x, N_Vector y)
{
  sunindextype  N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].v1 = NV_DATA_PT(x);
    thread_data[i].v2 = NV_DATA_PT(y);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, Vaxpy_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}

/* ----------------------------------------------------------------------------
 * Pthread companion function to Vaxpy
 */

static void *Vaxpy_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype a;
  realtype *xd, *yd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  a  = my_data->c1;
  xd = my_data->v1;
  yd = my_data->v2;

  start = my_data->start;
  end   = my_data->end;

  /* compute axpy */
  if (a == ONE) {
    for (i = start; i < end; i++)
      yd[i] += xd[i];

    /* exit */
    pthread_exit(NULL);
  }

  if (a == -ONE) {
    for (i = start; i < end; i++)
      yd[i] -= xd[i];

    /* exit */
    pthread_exit(NULL);
  }

  for (i = start; i < end; i++)
    yd[i] += a*xd[i];

  /* return */
  pthread_exit(NULL);
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector
 */

static void VScaleBy_Pthreads(realtype a, N_Vector x)
{
  sunindextype  N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(x);
  nthreads     = NV_NUM_THREADS_PT(x);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i<nthreads; i++) {
    /* initialize thread data */
    N_VInitThreadData(&thread_data[i]);

    /* compute start and end loop index for thread */
    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    /* pack thread data */
    thread_data[i].c1 = a;
    thread_data[i].v1 = NV_DATA_PT(x);

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, VScaleBy_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++) {
    pthread_join(threads[i], NULL);
  }

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return;
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to VScaleBy
 */

static void *VScaleBy_PT(void *thread_data)
{
  sunindextype i, start, end;
  realtype a;
  realtype *xd;
  Pthreads_Data *my_data;

  /* extract thread data */
  my_data = (Pthreads_Data *) thread_data;

  a  = my_data->c1;
  xd = my_data->v1;

  start = my_data->start;
  end   = my_data->end;

  /* compute scaled vector */
  for (i = start; i < end; i++)
    xd[i] *= a;

  /* exit */
  pthread_exit(NULL);
}


/*
 * -----------------------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------------------
 */

static int VSumVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i=0; i<nthreads; i++) {
    N_VInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VSumVectorArray_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++)
    pthread_join(threads[i], NULL);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}

static void *VSumVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  my_data = (Pthreads_Data *) thread_data;
  start   = my_data->start;
  end     = my_data->end;

  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j=start; j<end; j++)
      zd[j] = xd[j] + yd[j];
  }

  pthread_exit(NULL);
}


static int VDiffVectorArray_Pthreads(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i=0; i<nthreads; i++) {
    N_VInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VDiffVectorArray_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++)
    pthread_join(threads[i], NULL);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}

static void *VDiffVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  my_data = (Pthreads_Data *) thread_data;
  start   = my_data->start;
  end     = my_data->end;

  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j=start; j<end; j++)
      zd[j] = xd[j] - yd[j];
  }

  pthread_exit(NULL);
}


static int VScaleSumVectorArray_Pthreads(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N            = NV_LENGTH_PT(X[0]);
  nthreads     = NV_NUM_THREADS_PT(X[0]);
  threads      = malloc(nthreads*sizeof(pthread_t));
  thread_data  = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i=0; i<nthreads; i++) {
    N_VInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = c;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VScaleSumVectorArray_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++)
    pthread_join(threads[i], NULL);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}

static void *VScaleSumVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype  c;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  my_data = (Pthreads_Data *) thread_data;
  start   = my_data->start;
  end     = my_data->end;
  c       = my_data->c1;

  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j=start; j<end; j++)
      zd[j] = c * (xd[j] + yd[j]);
  }

  pthread_exit(NULL);
}


static int VScaleDiffVectorArray_Pthreads(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i=0; i<nthreads; i++) {
    N_VInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = c;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VScaleDiffVectorArray_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++)
    pthread_join(threads[i], NULL);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}

static void *VScaleDiffVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype  c;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  my_data = (Pthreads_Data *) thread_data;
  start   = my_data->start;
  end     = my_data->end;
  c       = my_data->c1;

  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j=start; j<end; j++)
      zd[j] = c * (xd[j] - yd[j]);
  }

  pthread_exit(NULL);
}


static int VLin1VectorArray_Pthreads(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i=0; i<nthreads; i++) {
    N_VInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = a;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VLin1VectorArray_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++)
    pthread_join(threads[i], NULL);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}

static void *VLin1VectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype  a;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  my_data = (Pthreads_Data *) thread_data;
  start   = my_data->start;
  end     = my_data->end;
  a       = my_data->c1;

  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j=start; j<end; j++)
      zd[j] = (a * xd[j]) + yd[j];
  }

  pthread_exit(NULL);
}


static int VLin2VectorArray_Pthreads(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i=0; i<nthreads; i++) {
    N_VInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = a;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;
    thread_data[i].Y3   = Z;

    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VLin2VectorArray_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++)
    pthread_join(threads[i], NULL);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}

static void *VLin2VectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype  a;
  realtype* xd=NULL;
  realtype* yd=NULL;
  realtype* zd=NULL;

  my_data = (Pthreads_Data *) thread_data;
  start   = my_data->start;
  end     = my_data->end;
  a       = my_data->c1;

  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    zd = NV_DATA_PT(my_data->Y3[i]);
    for (j=start; j<end; j++)
      zd[j] = (a * xd[j]) - yd[j];
  }

  pthread_exit(NULL);
}


static int VaxpyVectorArray_Pthreads(int nvec, realtype a, N_Vector* X, N_Vector* Y)
{
  sunindextype   N;
  int            i, nthreads;
  pthread_t      *threads;
  Pthreads_Data  *thread_data;
  pthread_attr_t attr;

  /* allocate threads and thread data structs */
  N           = NV_LENGTH_PT(X[0]);
  nthreads    = NV_NUM_THREADS_PT(X[0]);
  threads     = malloc(nthreads*sizeof(pthread_t));
  thread_data = (Pthreads_Data *) malloc(nthreads*sizeof(struct _Pthreads_Data));

  /* set thread attributes */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* pack thread data, distribute loop indices, and create threads/call kernel */
  for (i=0; i<nthreads; i++) {
    N_VInitThreadData(&thread_data[i]);

    thread_data[i].nvec = nvec;
    thread_data[i].c1   = a;
    thread_data[i].Y1   = X;
    thread_data[i].Y2   = Y;

    N_VSplitLoop(i, &nthreads, &N, &thread_data[i].start, &thread_data[i].end);

    pthread_create(&threads[i], &attr, VaxpyVectorArray_PT, (void *) &thread_data[i]);
  }

  /* wait for all threads to finish */
  for (i=0; i<nthreads; i++)
    pthread_join(threads[i], NULL);

  /* clean up and return */
  pthread_attr_destroy(&attr);
  free(threads);
  free(thread_data);

  return(0);
}

static void *VaxpyVectorArray_PT(void *thread_data)
{
  Pthreads_Data *my_data;
  sunindextype  j, start, end;

  int       i;
  realtype  a;
  realtype* xd=NULL;
  realtype* yd=NULL;

  my_data = (Pthreads_Data *) thread_data;
  start   = my_data->start;
  end     = my_data->end;
  a       = my_data->c1;

  if (a == ONE) {
    for (i=0; i<my_data->nvec; i++) {
      xd = NV_DATA_PT(my_data->Y1[i]);
      yd = NV_DATA_PT(my_data->Y2[i]);
      for (j=start; j<end; j++)
        yd[j] += xd[j];
    }
    pthread_exit(NULL);
  }

  if (a == -ONE) {
    for (i=0; i<my_data->nvec; i++) {
      xd = NV_DATA_PT(my_data->Y1[i]);
      yd = NV_DATA_PT(my_data->Y2[i]);
      for (j=start; j<end; j++)
        yd[j] -= xd[j];
    }
    pthread_exit(NULL);
  }

  for (i=0; i<my_data->nvec; i++) {
    xd = NV_DATA_PT(my_data->Y1[i]);
    yd = NV_DATA_PT(my_data->Y2[i]);
    for (j=start; j<end; j++)
      yd[j] += a * xd[j];
  }
  pthread_exit(NULL);
}


/*
 * -----------------------------------------------------------------
 * private utility functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Determine loop indices for a thread
 */

static void N_VSplitLoop(int myid, int *nthreads, sunindextype *N,
			 sunindextype *start, sunindextype *end)
{
  sunindextype q, r; /* quotient and remainder */

  /* work per thread and leftover work */
  q = *N / *nthreads;
  r = *N % *nthreads;

  /* assign work */
  if (myid < r) {
    *start = myid * q + myid;
    *end   = *start + q + 1;
  } else {
    *start = myid * q + r;
    *end   = *start + q;
  }
}


/* ----------------------------------------------------------------------------
 * Initialize values of local thread data struct
 */

static void N_VInitThreadData(Pthreads_Data *thread_data)
{
  thread_data->start = -1;
  thread_data->end   = -1;

#if __STDC_VERSION__ >= 199901L
  thread_data->c1 = NAN;
  thread_data->c2 = NAN;
#else
  thread_data->c1 = ZERO;
  thread_data->c2 = ZERO;
#endif

  thread_data->v1 = NULL;
  thread_data->v2 = NULL;
  thread_data->v3 = NULL;
  thread_data->global_val = NULL;
  thread_data->global_mutex = NULL;

  thread_data->nvec = ZERO;
  thread_data->nsum = ZERO;
  thread_data->cvals = NULL;
  thread_data->Y1 = NULL;
  thread_data->Y2 = NULL;
  thread_data->Y3 = NULL;
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Pthreads;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Pthreads;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Pthreads;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Pthreads;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Pthreads;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Pthreads;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Pthreads;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Pthreads;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Pthreads;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Pthreads;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = N_VDotProdMulti_Pthreads;
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
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Pthreads;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Pthreads;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf) {
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Pthreads;
    v->ops->nvdotprodmultilocal = N_VDotProdMulti_Pthreads;
  } else {
    v->ops->nvdotprodmulti      = NULL;
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Pthreads;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Pthreads;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Pthreads;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Pthreads;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Pthreads;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Pthreads;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Pthreads(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Pthreads;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}
