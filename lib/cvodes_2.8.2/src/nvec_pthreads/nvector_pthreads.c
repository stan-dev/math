/*
 * -----------------------------------------------------------------
 * $Revision: 4272 $
 * $Date: 2014-12-02 11:19:41 -0800 (Tue, 02 Dec 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the implementation file for a POSIX Threads (Pthreads)
 * implementation of the NVECTOR package using a LOCAL array of 
 * structures to pass data to threads.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_pthreads.h>
#include <sundials/sundials_math.h>
#include <math.h> /* define NAN */

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Private function prototypes */
/* z=x */
static void VCopy_Pthreads(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_Pthreads(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_Pthreads(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_Pthreads(N_Vector x, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_Pthreads(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_Pthreads(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_Pthreads(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_Pthreads(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_Pthreads(realtype a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_Pthreads(realtype a, N_Vector x);

/* Pthread companion function prototypes */
static void *N_VLinearSum_PT(void *thread_data);
static void *N_VConst_PT(void *thread_data);
static void *N_VProd_PT(void *thread_data);
static void *N_VDiv_PT(void *thread_data);
static void *N_VScale_PT(void *thread_data);
static void *N_VAbs_PT(void *thread_data);
static void *N_VInv_PT(void *thread_data);
static void *N_VAddConst_PT(void *thread_data);
static void *N_VCompare_PT(void *thread_data);
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
static void *N_VDotProd_PT(void *thread_data);
static void *N_VMaxNorm_PT(void *thread_data);
static void *N_VWrmsNorm_PT(void *thread_data);
static void *N_VMin_PT(void *thread_data);
static void *N_VWL2Norm_PT(void *thread_data);
static void *N_VL1Norm_PT(void *thread_data);
static void *N_VInvTest_PT(void *thread_data);
static void *N_VWrmsNormMask_PT(void *thread_data);
static void *N_VConstrMask_PT(void *thread_data);
static void *N_VMinQuotient_PT(void *thread_data);

/* Function to determine loop values for threads */
static void N_VSplitLoop(int myid, int *nthreads, long int *N, 
			 long int *start, long int *end);

/* Function to initialize thread data */
static void N_VInitThreadData(Pthreads_Data *thread_data);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new empty vector 
 */

N_Vector N_VNewEmpty_Pthreads(long int length, int num_threads)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Pthreads content;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvclone           = N_VClone_Pthreads;
  ops->nvcloneempty      = N_VCloneEmpty_Pthreads;
  ops->nvdestroy         = N_VDestroy_Pthreads;
  ops->nvspace           = N_VSpace_Pthreads;
  ops->nvgetarraypointer = N_VGetArrayPointer_Pthreads;
  ops->nvsetarraypointer = N_VSetArrayPointer_Pthreads;
  ops->nvlinearsum       = N_VLinearSum_Pthreads;
  ops->nvconst           = N_VConst_Pthreads;
  ops->nvprod            = N_VProd_Pthreads;
  ops->nvdiv             = N_VDiv_Pthreads;
  ops->nvscale           = N_VScale_Pthreads;
  ops->nvabs             = N_VAbs_Pthreads;
  ops->nvinv             = N_VInv_Pthreads;
  ops->nvaddconst        = N_VAddConst_Pthreads;
  ops->nvdotprod         = N_VDotProd_Pthreads;
  ops->nvmaxnorm         = N_VMaxNorm_Pthreads;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Pthreads;
  ops->nvwrmsnorm        = N_VWrmsNorm_Pthreads;
  ops->nvmin             = N_VMin_Pthreads;
  ops->nvwl2norm         = N_VWL2Norm_Pthreads;
  ops->nvl1norm          = N_VL1Norm_Pthreads;
  ops->nvcompare         = N_VCompare_Pthreads;
  ops->nvinvtest         = N_VInvTest_Pthreads;
  ops->nvconstrmask      = N_VConstrMask_Pthreads;
  ops->nvminquotient     = N_VMinQuotient_Pthreads;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Pthreads) malloc(sizeof(struct _N_VectorContent_Pthreads));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  content->length      = length;
  content->num_threads = num_threads;
  content->own_data    = FALSE;
  content->data        = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new vector 
 */

N_Vector N_VNew_Pthreads(long int length, int num_threads)
{
  N_Vector v;
  realtype *data;

  v = NULL;
  v = N_VNewEmpty_Pthreads(length, num_threads);
  if (v == NULL) return(NULL);

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (realtype *) malloc(length * sizeof(realtype));
    if(data == NULL) { N_VDestroy_Pthreads(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_PT(v) = TRUE;
    NV_DATA_PT(v)     = data;

  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a vector with user data component 
 */

N_Vector N_VMake_Pthreads(long int length, int num_threads, realtype *v_data)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Pthreads(length, num_threads);
  if (v == NULL) return(NULL);

  if (length > 0) {
    /* Attach data */
    NV_OWN_DATA_PT(v) = FALSE;
    NV_DATA_PT(v)     = v_data;
  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new vectors. 
 */

N_Vector *N_VCloneVectorArray_Pthreads(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_Pthreads(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Pthreads(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new vectors with NULL data array. 
 */

N_Vector *N_VCloneVectorArrayEmpty_Pthreads(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_Pthreads(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Pthreads(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Pthreads
 */

void N_VDestroyVectorArray_Pthreads(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Pthreads(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Function to print a vector 
 */
 
void N_VPrint_Pthreads(N_Vector x)
{
  long int i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LENGTH_PT(x);
  xd = NV_DATA_PT(x);

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

N_Vector N_VCloneEmpty_Pthreads(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Pthreads content;

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
  content = (N_VectorContent_Pthreads) malloc(sizeof(struct _N_VectorContent_Pthreads));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  content->length      = NV_LENGTH_PT(w);
  content->num_threads = NV_NUM_THREADS_PT(w);
  content->own_data    = FALSE;
  content->data        = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}


/* ----------------------------------------------------------------------------
 * Create new vector from existing vector and attach data
 */

N_Vector N_VClone_Pthreads(N_Vector w)
{
  N_Vector v;
  realtype *data;
  long int length;

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
    NV_OWN_DATA_PT(v) = TRUE;
    NV_DATA_PT(v)     = data;

  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Destroy vector and free vector memory
 */

void N_VDestroy_Pthreads(N_Vector v)
{
  if (NV_OWN_DATA_PT(v) == TRUE) {
    free(NV_DATA_PT(v));
    NV_DATA_PT(v) = NULL;
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Get storage requirement for vector
 */

void N_VSpace_Pthreads(N_Vector v, long int *lrw, long int *liw)
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
 * Compute linear combination z[i] = a*x[i]+b*y[i]
 */

void N_VLinearSum_Pthreads(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  realtype c;
  N_Vector v1, v2;
  booleantype test;

  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;

  pthread_attr_t attr;

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
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;

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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  if (z == x) {  /* BLAS usage: scale x <- cx */
    VScaleBy_Pthreads(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_Pthreads(x, z);
  } 
  else if (c == -ONE) {
    VNeg_Pthreads(x, z);
  } 
  else {      
      long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      sum = ZERO;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;  

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
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      max = ZERO;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;  

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
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      sum = ZERO;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;  

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
    thread_data[i].global_val   = &sum;
    thread_data[i].global_mutex = &global_mutex;

    /* create threads and call pthread companion function */
    pthread_create(&threads[i], &attr, N_VWrmsNorm_PT, (void *) &thread_data[i]);
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
  
  return(SUNRsqrt(sum/N));
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWrmsNorm
 */

static void *N_VWrmsNorm_PT(void *thread_data)
{
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      sum = ZERO;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;  

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
    pthread_create(&threads[i], &attr, N_VWrmsNormMask_PT, (void *) &thread_data[i]);
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
  
  return(SUNRsqrt(sum/N));
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VWrmsNormMask
 */

static void *N_VWrmsNormMask_PT(void *thread_data)
{
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      min;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;  

  /* initialize global min */
  min = NV_Ith_PT(x,0);

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
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      sum = ZERO;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;  

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
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      sum = ZERO;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;  

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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      val = ZERO;
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
    return (FALSE);
  else
    return (TRUE);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VInvTest
 */

static void *N_VInvTest_PT(void *thread_data)
{
  long int i, start, end;
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      val = ZERO;
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
    return(FALSE);
  else 
    return(TRUE);
}


/* ----------------------------------------------------------------------------
 * Pthread companion function to N_VConstrMask
 */

static void *N_VConstrMask_PT(void *thread_data)
{
  long int i, start, end;
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
    
    /* c[i] = 0, do nothing */
    if (cd[i] == ZERO) 
      continue;

    /* c[i] = +/- 2, check x[i] > or < 0 */
    if (cd[i] > ONEPT5 || cd[i] < -ONEPT5) {
      if (xd[i]*cd[i] <= ZERO) {
	local_val = ONE; 
	md[i] = ONE;
      }
      continue;
    }

    /* c[i] = +/- 1, check x[i] >= or <= 0 */
    if (cd[i] > HALF || cd[i] < -HALF) {
      if (xd[i]*cd[i] < ZERO) {
	local_val = ONE; 
	md[i] = ONE;
      }
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
  long int      N;
  int           i, nthreads;
  pthread_t     *threads;
  Pthreads_Data *thread_data;
  realtype      min = BIG_REAL;
  pthread_attr_t attr;
  pthread_mutex_t global_mutex;  

  /* allocate threads and thread data structs */ 
  N            = NV_LENGTH_PT(num);
  nthreads     = NV_NUM_THREADS_PT(num);
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
  long int i, start, end;
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
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */


/* ----------------------------------------------------------------------------
 * Copy vector components into second vector   
 */

static void VCopy_Pthreads(N_Vector x, N_Vector z)
{
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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
  long int      N;
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
  long int i, start, end;
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


/* ----------------------------------------------------------------------------
 * Determine loop indices for a thread
 */

static void N_VSplitLoop(int myid, int *nthreads, long int *N, 
			 long int *start, long int *end)
{
  long int q, r; /* quotient and remainder */

  /* work per thread and leftover work */
  q = *N / *nthreads;
  r = *N % *nthreads;
  
  /* assign work */
  if (myid < r) {
    *start = myid * q + myid;
    *end   = *start + q + 1;
  }
  else {
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
}
