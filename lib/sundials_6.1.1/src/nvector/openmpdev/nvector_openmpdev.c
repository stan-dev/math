/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner and Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: This NVECTOR module is based on the NVECTOR
 *                   Serial module by Scott D. Cohen, Alan C.
 *                   Hindmarsh, Radu Serban, and Aaron Collier
 *                   @ LLNL
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
 * This is the implementation file for an OpenMP DEV implementation
 * of the NVECTOR module.
 * -----------------------------------------------------------------*/

#include <omp.h>

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_openmpdev.h>
#include <sundials/sundials_math.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Private functions for special cases of vector operations */
static void VCopy_OpenMPDEV(N_Vector x, N_Vector z);                              /* z=x       */
static void VSum_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z);                   /* z=x+y     */
static void VDiff_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z);                  /* z=x-y     */
static void VNeg_OpenMPDEV(N_Vector x, N_Vector z);                               /* z=-x      */
static void VScaleSum_OpenMPDEV(realtype c, N_Vector x, N_Vector y, N_Vector z);  /* z=c(x+y)  */
static void VScaleDiff_OpenMPDEV(realtype c, N_Vector x, N_Vector y, N_Vector z); /* z=c(x-y)  */
static void VLin1_OpenMPDEV(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax+y    */
static void VLin2_OpenMPDEV(realtype a, N_Vector x, N_Vector y, N_Vector z);      /* z=ax-y    */
static void Vaxpy_OpenMPDEV(realtype a, N_Vector x, N_Vector y);                  /* y <- ax+y */
static void VScaleBy_OpenMPDEV(realtype a, N_Vector x);                           /* x <- ax   */

/* Private functions for special cases of vector array operations */
static int VSumVectorArray_OpenMPDEV(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z);                   /* Z=X+Y     */
static int VDiffVectorArray_OpenMPDEV(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z);                  /* Z=X-Y     */
static int VScaleSumVectorArray_OpenMPDEV(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z);  /* Z=c(X+Y)  */
static int VScaleDiffVectorArray_OpenMPDEV(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z); /* Z=c(X-Y)  */
static int VLin1VectorArray_OpenMPDEV(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z);      /* Z=aX+Y    */
static int VLin2VectorArray_OpenMPDEV(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z);      /* Z=aX-Y    */
static int VaxpyVectorArray_OpenMPDEV(int nvec, realtype a, N_Vector* X, N_Vector* Y);                   /* Y <- aX+Y */

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_OpenMPDEV(N_Vector v)
{
  return SUNDIALS_NVEC_OPENMPDEV;
}

/* ----------------------------------------------------------------------------
 * Function to create a new empty vector
 */

N_Vector N_VNewEmpty_OpenMPDEV(sunindextype length, SUNContext sunctx)
{
  N_Vector v;
  N_VectorContent_OpenMPDEV content;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = N_VGetVectorID_OpenMPDEV;
  v->ops->nvclone                 = N_VClone_OpenMPDEV;
  v->ops->nvcloneempty            = N_VCloneEmpty_OpenMPDEV;
  v->ops->nvdestroy               = N_VDestroy_OpenMPDEV;
  v->ops->nvspace                 = N_VSpace_OpenMPDEV;
  v->ops->nvgetlength             = N_VGetLength_OpenMPDEV;
  v->ops->nvgetarraypointer       = N_VGetHostArrayPointer_OpenMPDEV;
  v->ops->nvgetdevicearraypointer = N_VGetDeviceArrayPointer_OpenMPDEV;
  v->ops->nvprint                 = N_VPrint_OpenMPDEV;
  v->ops->nvprintfile             = N_VPrintFile_OpenMPDEV;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_OpenMPDEV;
  v->ops->nvconst        = N_VConst_OpenMPDEV;
  v->ops->nvprod         = N_VProd_OpenMPDEV;
  v->ops->nvdiv          = N_VDiv_OpenMPDEV;
  v->ops->nvscale        = N_VScale_OpenMPDEV;
  v->ops->nvabs          = N_VAbs_OpenMPDEV;
  v->ops->nvinv          = N_VInv_OpenMPDEV;
  v->ops->nvaddconst     = N_VAddConst_OpenMPDEV;
  v->ops->nvdotprod      = N_VDotProd_OpenMPDEV;
  v->ops->nvmaxnorm      = N_VMaxNorm_OpenMPDEV;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_OpenMPDEV;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_OpenMPDEV;
  v->ops->nvmin          = N_VMin_OpenMPDEV;
  v->ops->nvwl2norm      = N_VWL2Norm_OpenMPDEV;
  v->ops->nvl1norm       = N_VL1Norm_OpenMPDEV;
  v->ops->nvcompare      = N_VCompare_OpenMPDEV;
  v->ops->nvinvtest      = N_VInvTest_OpenMPDEV;
  v->ops->nvconstrmask   = N_VConstrMask_OpenMPDEV;
  v->ops->nvminquotient  = N_VMinQuotient_OpenMPDEV;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProd_OpenMPDEV;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_OpenMPDEV;
  v->ops->nvminlocal         = N_VMin_OpenMPDEV;
  v->ops->nvl1normlocal      = N_VL1Norm_OpenMPDEV;
  v->ops->nvinvtestlocal     = N_VInvTest_OpenMPDEV;
  v->ops->nvconstrmasklocal  = N_VConstrMask_OpenMPDEV;
  v->ops->nvminquotientlocal = N_VMinQuotient_OpenMPDEV;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_OpenMPDEV;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_OpenMPDEV;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = N_VDotProdMulti_OpenMPDEV;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_OpenMPDEV) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length    = length;
  content->own_data  = SUNFALSE;
  content->host_data = NULL;
  content->dev_data  = NULL;

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new vector
 */

N_Vector N_VNew_OpenMPDEV(sunindextype length)
{
  N_Vector v;
  realtype *data;
  realtype *dev_data;
  int dev;

  v = NULL;
  v = N_VNewEmpty_OpenMPDEV(length);
  if (v == NULL) return(NULL);

  /* Create data */
  if (length > 0) {

    /* Update ownership */
    NV_OWN_DATA_OMPDEV(v) = SUNTRUE;

    /* Allocate memory on host */
    data = NULL;
    data = (realtype *) malloc(length * sizeof(realtype));
    if (data == NULL) { N_VDestroy(v); return(NULL); }

    /* Allocate memory on device */
    dev = omp_get_default_device();
    dev_data = omp_target_alloc(length * sizeof(realtype), dev);
    if (dev_data == NULL) { N_VDestroy(v); return(NULL); }

    /* Attach data */
    NV_DATA_HOST_OMPDEV(v) = data;
    NV_DATA_DEV_OMPDEV(v)  = dev_data;

  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a vector with user data component
 */

N_Vector N_VMake_OpenMPDEV(sunindextype length, realtype *h_vdata,
                           realtype *d_vdata)
{
  N_Vector v;
  int dev, host;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_OpenMPDEV(length);
  if (v == NULL) return(NULL);

  if (length > 0) {
    /* Get device and host identifiers */
    dev  = omp_get_default_device();
    host = omp_get_initial_device();

    /* Attach data */
    NV_OWN_DATA_OMPDEV(v)  = SUNFALSE;
    NV_DATA_HOST_OMPDEV(v) = h_vdata;
    NV_DATA_DEV_OMPDEV(v)  = d_vdata;
  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new vectors.
 */

N_Vector *N_VCloneVectorArray_OpenMPDEV(int count, N_Vector w)
{
  return(N_VCloneVectorArray(count, w));
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new vectors with NULL data array.
 */

N_Vector *N_VCloneVectorArrayEmpty_OpenMPDEV(int count, N_Vector w)
{
  return(N_VCloneEmptyVectorArray(count, w));
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_OpenMPDEV
 */

void N_VDestroyVectorArray_OpenMPDEV(N_Vector *vs, int count)
{
  N_VDestroyVectorArray(vs, count);
  return;
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_OpenMPDEV(N_Vector v)
{
  return NV_LENGTH_OMPDEV(v);
}

/* ----------------------------------------------------------------------------
 * Function to return a pointer to the data array on the host.
 */
realtype *N_VGetHostArrayPointer_OpenMPDEV(N_Vector v)
{
  return((realtype *) NV_DATA_HOST_OMPDEV(v));
}

/* ----------------------------------------------------------------------------
 * Function to return a pointer to the data array on the device.
 */
realtype *N_VGetDeviceArrayPointer_OpenMPDEV(N_Vector v)
{
  return((realtype *) NV_DATA_DEV_OMPDEV(v));
}

/* ----------------------------------------------------------------------------
 * Function to print a vector to stdout
 */

void N_VPrint_OpenMPDEV(N_Vector x)
{
  N_VPrintFile_OpenMPDEV(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print a vector to outfile
 */

void N_VPrintFile_OpenMPDEV(N_Vector x, FILE *outfile)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LENGTH_OMPDEV(x);
  xd = NV_DATA_HOST_OMPDEV(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%11.8Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%11.8g\n", xd[i]);
#else
    fprintf(outfile, "%11.8g\n", xd[i]);
#endif
  }
  fprintf(outfile, "\n");

  return;
}

/* ----------------------------------------------------------------------------
 * Function to copy host array into device array
 */

void N_VCopyToDevice_OpenMPDEV(N_Vector x)
{
  int dev, host;
  sunindextype length;
  realtype *host_ptr;
  realtype *dev_ptr;

  /* Get array information */
  length   = NV_LENGTH_OMPDEV(x);
  host_ptr = NV_DATA_HOST_OMPDEV(x);
  dev_ptr  = NV_DATA_DEV_OMPDEV(x);

  /* Get device and host identifiers */
  dev  = omp_get_default_device();
  host = omp_get_initial_device();

  /* Copy array from host to device */
  omp_target_memcpy(dev_ptr, host_ptr, sizeof(realtype) * length, 0, 0, dev, host);

  return;
}

/* ----------------------------------------------------------------------------
 * Function to copy device array into host array
 */

void N_VCopyFromDevice_OpenMPDEV(N_Vector x)
{
  int dev, host;
  sunindextype length;
  realtype *host_ptr;
  realtype *dev_ptr;

  /* Get array information */
  length   = NV_LENGTH_OMPDEV(x);
  host_ptr = NV_DATA_HOST_OMPDEV(x);
  dev_ptr  = NV_DATA_DEV_OMPDEV(x);

  /* Get device and host identifiers */
  dev  = omp_get_default_device();
  host = omp_get_initial_device();

  /* Copy array from device to host */
  omp_target_memcpy(host_ptr, dev_ptr, sizeof(realtype) * length, 0, 0, host, dev);

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

N_Vector N_VCloneEmpty_OpenMPDEV(N_Vector w)
{
  N_Vector v;
  N_VectorContent_OpenMPDEV content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Create content */
  content = NULL;
  content = (N_VectorContent_OpenMPDEV) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length    = NV_LENGTH_OMPDEV(w);
  content->own_data  = SUNFALSE;
  content->host_data = NULL;
  content->dev_data  = NULL;

  return(v);
}


/* ----------------------------------------------------------------------------
 * Create new vector from existing vector and attach data
 */

N_Vector N_VClone_OpenMPDEV(N_Vector w)
{
  N_Vector v;
  realtype *data;
  realtype *dev_data;
  sunindextype length;
  int dev;

  v = NULL;
  v = N_VCloneEmpty_OpenMPDEV(w);
  if (v == NULL) return(NULL);

  length = NV_LENGTH_OMPDEV(w);

  /* Create data */
  if (length > 0) {

    /* Update ownership flag */
    NV_OWN_DATA_OMPDEV(v) = SUNTRUE;

    /* Allocate memory on host */
    data = NULL;
    data = (realtype *) malloc(length * sizeof(realtype));
    if (data == NULL) { N_VDestroy(v); return(NULL); }

    /* Allocate memory on device */
    dev = omp_get_default_device();
    dev_data = omp_target_alloc(length * sizeof(realtype), dev);
    if (dev_data == NULL) { N_VDestroy(v); return(NULL); }

    /* Attach data */
    NV_DATA_HOST_OMPDEV(v)= data;
    NV_DATA_DEV_OMPDEV(v) = dev_data;

  }

  return(v);
}


/* ----------------------------------------------------------------------------
 * Destroy vector and free vector memory
 */

void N_VDestroy_OpenMPDEV(N_Vector v)
{
  int dev;

  if (v == NULL) return;

  /* free content */
  if (v->content != NULL) {
    /* free data arrays if they are owned by the vector */
    if (NV_OWN_DATA_OMPDEV(v)) {

      if (NV_DATA_HOST_OMPDEV(v) != NULL) {
        free(NV_DATA_HOST_OMPDEV(v));
        NV_DATA_HOST_OMPDEV(v) = NULL;
      }

      if (NV_DATA_DEV_OMPDEV(v) != NULL) {
        dev = omp_get_default_device();
        omp_target_free(NV_DATA_DEV_OMPDEV(v), dev);
        NV_DATA_DEV_OMPDEV(v) = NULL;
      }
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
 * Get storage requirement for N_Vector
 */

void N_VSpace_OpenMPDEV(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NV_LENGTH_OMPDEV(v);
  *liw = 1;

  return;
}

/* ----------------------------------------------------------------------------
 * Compute linear combination z[i] = a*x[i]+b*y[i]
 */

void N_VLinearSum_OpenMPDEV(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype c, *xd_dev, *yd_dev, *zd_dev;
  N_Vector v1, v2;
  booleantype test;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy_OpenMPDEV(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy_OpenMPDEV(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_OpenMPDEV(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_OpenMPDEV(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_OpenMPDEV(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_OpenMPDEV(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_OpenMPDEV(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_OpenMPDEV(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = (a*xd_dev[i])+(b*yd_dev[i]);

  return;
}


/* ----------------------------------------------------------------------------
 * Assigns constant value to all vector elements, z[i] = c
 */

void N_VConst_OpenMPDEV(realtype c, N_Vector z)
{
  sunindextype i, N;
  realtype *zd_dev;
  int dev;

  zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(z);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev  = omp_get_default_device();

#pragma omp target is_device_ptr(zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
    for (i = 0; i < N; i++) zd_dev[i] = c;
  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise product z[i] = x[i]*y[i]
 */

void N_VProd_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev, *zd_dev;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev  = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = xd_dev[i]*yd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise division z[i] = x[i]/y[i]
 */

void N_VDiv_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev, *zd_dev;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev  = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = xd_dev[i]/yd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute scaler multiplication z[i] = c*x[i]
 */

void N_VScale_OpenMPDEV(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *zd_dev;
  int dev;

  xd_dev = zd_dev = NULL;

  if (z == x) {  /* BLAS usage: scale x <- cx */
    VScaleBy_OpenMPDEV(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_OpenMPDEV(x, z);
  } else if (c == -ONE) {
    VNeg_OpenMPDEV(x, z);
  } else {
    N  = NV_LENGTH_OMPDEV(x);
    xd_dev = NV_DATA_DEV_OMPDEV(x);
    zd_dev = NV_DATA_DEV_OMPDEV(z);

    /* get default device identifier */
    dev  = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
    for (i = 0; i < N; i++)
      zd_dev[i] = c*xd_dev[i];
  }

  return;
}


/* ----------------------------------------------------------------------------
 * Compute absolute value of vector components z[i] = SUNRabs(x[i])
 */

void N_VAbs_OpenMPDEV(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *zd_dev;
  int dev;

  xd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev  = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = SUNRabs(xd_dev[i]);

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = 1 / x[i]
 */

void N_VInv_OpenMPDEV(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *zd_dev;
  int dev;

  xd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = ONE/xd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise addition of a scaler to a vector z[i] = x[i] + b
 */

void N_VAddConst_OpenMPDEV(N_Vector x, realtype b, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *zd_dev;
  int dev;

  xd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = xd_dev[i]+b;

  return;
}


/* ----------------------------------------------------------------------------
 * Computes the dot product of two vectors, a = sum(x[i]*y[i])
 */

realtype N_VDotProd_OpenMPDEV(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  realtype sum, *xd_dev, *yd_dev;
  int dev;

  xd_dev = yd_dev = NULL;

  sum = ZERO;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target map(tofrom:sum) is_device_ptr(xd_dev, yd_dev) device(dev)
#pragma omp teams distribute parallel for reduction(+:sum) schedule(static, 1)
  for (i = 0; i < N; i++) {
    sum += xd_dev[i]*yd_dev[i];
  }

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Computes max norm of a vector
 */

realtype N_VMaxNorm_OpenMPDEV(N_Vector x)
{
  sunindextype i, N;
  realtype max, *xd_dev;
  int dev;

  max = ZERO;
  xd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target map(tofrom:max) is_device_ptr(xd_dev) device(dev)
#pragma omp teams distribute parallel for reduction(max:max) schedule(static, 1)
    for (i = 0; i < N; i++) {
      max = SUNMAX(SUNRabs(xd_dev[i]), max);
    }

  return(max);
}


/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a vector
 */

realtype N_VWrmsNorm_OpenMPDEV(N_Vector x, N_Vector w)
{
  return(SUNRsqrt(N_VWSqrSumLocal_OpenMPDEV(x, w)/(NV_LENGTH_OMPDEV(x))));
}


/* ----------------------------------------------------------------------------
 * Computes weighted root mean square norm of a masked vector
 */

realtype N_VWrmsNormMask_OpenMPDEV(N_Vector x, N_Vector w, N_Vector id)
{
  return(SUNRsqrt(N_VWSqrSumMaskLocal_OpenMPDEV(x, w, id) / (NV_LENGTH_OMPDEV(x))));
}


/* ----------------------------------------------------------------------------
 * Computes weighted square sum of a vector
 */

realtype N_VWSqrSumLocal_OpenMPDEV(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum, *xd_dev, *wd_dev;
  int dev;

  sum = ZERO;
  xd_dev = wd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  wd_dev = NV_DATA_DEV_OMPDEV(w);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target map(tofrom:sum) is_device_ptr(xd_dev, wd_dev) device(dev)
#pragma omp teams distribute parallel for reduction(+:sum) schedule(static, 1)
  for (i = 0; i < N; i++) {
    sum += SUNSQR(xd_dev[i]*wd_dev[i]);
  }

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Computes weighted square sum of a masked vector
 */

realtype N_VWSqrSumMaskLocal_OpenMPDEV(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  realtype sum, *xd_dev, *wd_dev, *idd_dev;
  int dev;

  sum = ZERO;
  xd_dev = wd_dev = idd_dev = NULL;

  N       = NV_LENGTH_OMPDEV(x);
  xd_dev  = NV_DATA_DEV_OMPDEV(x);
  wd_dev  = NV_DATA_DEV_OMPDEV(w);
  idd_dev = NV_DATA_DEV_OMPDEV(id);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target map(tofrom:sum) is_device_ptr(xd_dev, wd_dev, idd_dev) device(dev)
#pragma omp teams distribute parallel for reduction(+:sum) schedule(static, 1)
  for (i = 0; i < N; i++) {
    if (idd_dev[i] > ZERO) {
      sum += SUNSQR(xd_dev[i]*wd_dev[i]);
    }
  }

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Finds the minimun component of a vector
 */

realtype N_VMin_OpenMPDEV(N_Vector x)
{
  sunindextype i, N;
  realtype min, *xd_dev;
  int dev;

  xd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target map(from:min) is_device_ptr(xd_dev) device(dev)
#pragma omp teams num_teams(1)
  {
    min = xd_dev[0];
#pragma omp distribute parallel for reduction(min:min) schedule(static, 1)
    for (i = 1; i < N; i++) {
      min = SUNMIN(xd_dev[i], min);
    }
  }

  return(min);
}


/* ----------------------------------------------------------------------------
 * Computes weighted L2 norm of a vector
 */

realtype N_VWL2Norm_OpenMPDEV(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum, *xd_dev, *wd_dev;
  int dev;

  sum = ZERO;
  xd_dev = wd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  wd_dev = NV_DATA_DEV_OMPDEV(w);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target map(tofrom:sum) is_device_ptr(xd_dev, wd_dev) device(dev)
#pragma omp teams distribute parallel for reduction(+:sum) schedule(static, 1)
  for (i = 0; i < N; i++) {
    sum += SUNSQR(xd_dev[i]*wd_dev[i]);
  }

  return(SUNRsqrt(sum));
}


/* ----------------------------------------------------------------------------
 * Computes L1 norm of a vector
 */

realtype N_VL1Norm_OpenMPDEV(N_Vector x)
{
  sunindextype i, N;
  realtype sum, *xd_dev;
  int dev;

  sum = ZERO;
  xd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target map(tofrom:sum) is_device_ptr(xd_dev) device(dev)
#pragma omp teams distribute parallel for reduction(+:sum) schedule(static, 1)
  for (i = 0; i<N; i++)
    sum += SUNRabs(xd_dev[i]);

  return(sum);
}


/* ----------------------------------------------------------------------------
 * Compare vector component values to a scaler
 */

void N_VCompare_OpenMPDEV(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *zd_dev;
  int dev;

  xd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = (SUNRabs(xd_dev[i]) >= c) ? ONE : ZERO;

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise inverse z[i] = ONE/x[i] and checks if x[i] == ZERO
 */

booleantype N_VInvTest_OpenMPDEV(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *zd_dev, val;
  int dev;

  xd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

  val = ZERO;

#pragma omp target map(tofrom:val) is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for reduction(max:val) schedule(static, 1)
  for (i = 0; i < N; i++) {
    if (xd_dev[i] == ZERO)
      val = ONE;
    else
      zd_dev[i] = ONE/xd_dev[i];
  }

  if (val > ZERO)
    return (SUNFALSE);
  else
    return (SUNTRUE);
}


/* ----------------------------------------------------------------------------
 * Compute constraint mask of a vector
 */

booleantype N_VConstrMask_OpenMPDEV(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  realtype temp;
  realtype *cd_dev, *xd_dev, *md_dev;
  int dev;

  cd_dev = xd_dev = md_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  cd_dev = NV_DATA_DEV_OMPDEV(c);
  md_dev = NV_DATA_DEV_OMPDEV(m);

  /* get default device identifier */
  dev = omp_get_default_device();

  temp = ONE;

#pragma omp target map(tofrom:temp) is_device_ptr(xd_dev, cd_dev, md_dev) device(dev)
#pragma omp teams distribute parallel for reduction(min:temp) schedule(static, 1)
  for (i = 0; i < N; i++) {
    md_dev[i] = ZERO;
    if (cd_dev[i] == ZERO) continue;
    if (cd_dev[i] > ONEPT5 || cd_dev[i] < -ONEPT5) {
      if ( xd_dev[i]*cd_dev[i] <= ZERO) { temp = ZERO; md_dev[i] = ONE; }
      continue;
    }
    if ( cd_dev[i] > HALF || cd_dev[i] < -HALF) {
      if (xd_dev[i]*cd_dev[i] < ZERO ) { temp = ZERO; md_dev[i] = ONE; }
    }
  }

  if (temp == ONE) return (SUNTRUE);
  else return(SUNFALSE);
}


/* ----------------------------------------------------------------------------
 * Compute minimum componentwise quotient
 */

realtype N_VMinQuotient_OpenMPDEV(N_Vector num, N_Vector denom)
{
  sunindextype i, N;
  realtype *nd_dev, *dd_dev, min;
  int dev;

  nd_dev = dd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(num);
  nd_dev = NV_DATA_DEV_OMPDEV(num);
  dd_dev = NV_DATA_DEV_OMPDEV(denom);

  /* get default device identifier */
  dev = omp_get_default_device();

  min = BIG_REAL;

#pragma omp target map(tofrom:min) is_device_ptr(nd_dev, dd_dev) device(dev)
#pragma omp teams distribute parallel for reduction(min:min) schedule(static, 1)
  for (i = 0; i < N; i++)
    if (dd_dev[i] != ZERO)  min = SUNMIN(nd_dev[i]/dd_dev[i], min);

  return(min);
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination_OpenMPDEV(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int          i, dev;
  realtype     to_add; /* temporary variable to hold sum being added in atomic operation */
  sunindextype j, N;
  realtype*    zd_dev=NULL;
  realtype*    xd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_OpenMPDEV(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_OpenMPDEV(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  /* get vector length and data array */
  N      = NV_LENGTH_OMPDEV(z);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store X dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE)) {
#pragma omp target map(to:N,nvec,c[:nvec],xd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,zd_dev) device(dev)
#pragma omp teams distribute
    {
      for (i=1; i<nvec; i++) {
        xd_dev = xd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
        for (j=0; j<N; j++) {
          to_add = c[i] * xd_dev[j];
#pragma omp atomic
          zd_dev[j] += to_add;
        }
      }
    }
    free(xd_dev_ptrs);
    return(0);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z) {
#pragma omp target map(to:N,nvec,c[:nvec],xd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,zd_dev)
    {
#pragma omp teams distribute parallel for schedule(static,1)
      for (j=0; j<N; j++)
        zd_dev[j] *= c[0];
    }

#pragma omp target map(to:N,nvec,c[:nvec],xd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,zd_dev)
#pragma omp teams distribute
    {
      for (i=1; i<nvec; i++) {
        xd_dev = xd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
        for (j=0; j<N; j++) {
          to_add = c[i] * xd_dev[j];
#pragma omp atomic
          zd_dev[j] += to_add;
        }
      }
    }
    free(xd_dev_ptrs);
    return(0);
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd_dev = NV_DATA_DEV_OMPDEV(X[0]);
#pragma omp target map(to:N,c[:nvec]) \
  is_device_ptr(xd_dev, zd_dev) device(dev)
  {
#pragma omp teams distribute parallel for schedule(static, 1)
    for (j=0; j<N; j++) {
      zd_dev[j] = c[0] * xd_dev[j];
    }
  }

#pragma omp target map(to:N,nvec,c[:nvec],xd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=1; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++) {
        to_add = c[i] * xd_dev[j];
#pragma omp atomic
        zd_dev[j] += to_add;
      }
    }
  }
  free(xd_dev_ptrs);
  return(0);
}

int N_VScaleAddMulti_OpenMPDEV(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_OpenMPDEV(a[0], x, ONE, Y[0], Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
#pragma omp target map(to:N,nvec,a[:nvec],yd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev, yd_dev) device(dev)
#pragma omp teams distribute
    {
      for (i=0; i<nvec; i++) {
        yd_dev = yd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
        for (j=0; j<N; j++)
          yd_dev[j] += a[i] * xd_dev[j];
      }
    }
    free(yd_dev_ptrs);
    return(0);
  }

  /* Allocate and store dev pointers to copy to device */
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
#pragma omp target map(to:N,nvec,a[:nvec],yd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      yd_dev = yd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = a[i] * xd_dev[j] + yd_dev[j];
    }
  }
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

int N_VDotProdMulti_OpenMPDEV(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  int          i, dev;
  sunindextype j, N;
  realtype     sum;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype**   yd_dev_ptrs=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VDotProd */
  if (nvec == 1) {
    dotprods[0] = N_VDotProd_OpenMPDEV(x, Y[0]);
    return(0);
  }

  /* get vector length and data array */
  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* initialize dot products */
  for (i=0; i<nvec; i++) {
    dotprods[i] = ZERO;
  }

  /* Allocate and store dev pointers to copy to device */
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);

  /* compute multiple dot products */
#pragma omp target map(to:N,nvec,yd_dev_ptrs[:nvec]) map(tofrom:dotprods[:nvec]) \
  is_device_ptr(xd_dev,yd_dev) device(dev)
#pragma omp teams distribute
  for (i=0; i<nvec; i++) {
    yd_dev = yd_dev_ptrs[i];
    sum = ZERO;
#pragma omp parallel for reduction(+:sum) schedule(static, 1)
    for (j=0; j<N; j++)
      sum += xd_dev[j] * yd_dev[j];
    dotprods[i] += sum;
  }

  free(yd_dev_ptrs);
  return(0);
}


/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

int N_VLinearSumVectorArray_OpenMPDEV(int nvec,
                                     realtype a, N_Vector* X,
                                     realtype b, N_Vector* Y,
                                     N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  N_Vector*    V1;
  N_Vector*    V2;
  booleantype  test;
  realtype     c;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_OpenMPDEV(a, X[0], b, Y[0], Z[0]);
    return(0);
  }

  /* BLAS usage: axpy y <- ax+y */
  if ((b == ONE) && (Z == Y))
    return(VaxpyVectorArray_OpenMPDEV(nvec, a, X, Y));

  /* BLAS usage: axpy x <- by+x */
  if ((a == ONE) && (Z == X))
    return(VaxpyVectorArray_OpenMPDEV(nvec, b, Y, X));

  /* Case: a == b == 1.0 */
  if ((a == ONE) && (b == ONE))
    return(VSumVectorArray_OpenMPDEV(nvec, X, Y, Z));

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VDiffVectorArray_OpenMPDEV(nvec, V2, V1, Z));
  }

  /* Cases:                                                  */
  /*   (1) a == 1.0, b == other or 0.0,                      */
  /*   (2) a == other or 0.0, b == 1.0                       */
  /* if a or b is 0.0, then user should have called N_VScale */
  if ((test = (a == ONE)) || (b == ONE)) {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VLin1VectorArray_OpenMPDEV(nvec, c, V1, V2, Z));
  }

  /* Cases:                     */
  /*   (1) a == -1.0, b != 1.0, */
  /*   (2) a != 1.0, b == -1.0  */
  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    return(VLin2VectorArray_OpenMPDEV(nvec, c, V1, V2, Z));
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
    return(VScaleSumVectorArray_OpenMPDEV(nvec, a, X, Y, Z));

  /* Case: a == -b */
  if (a == -b)
    return(VScaleDiffVectorArray_OpenMPDEV(nvec, a, X, Y, Z));

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length */
  N = NV_LENGTH_OMPDEV(Z[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

  /* compute linear sum for each vector pair in vector arrays */
#pragma omp target map(to:N,nvec,a,b,xd_dev_ptrs[:nvec], yd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      yd_dev = yd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = a * xd_dev[j] + b * yd_dev[j];
    }
  }

  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

int N_VScaleVectorArray_OpenMPDEV(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_OpenMPDEV(c[0], X[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LENGTH_OMPDEV(Z[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++) {
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  }

  /*
   * X[i] *= c[i]
   */
  if (X == Z) {
#pragma omp target map(to:N,nvec,c[:nvec],xd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev) device(dev)
#pragma omp teams distribute
    {
      for (i=0; i<nvec; i++) {
        xd_dev = xd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
        for (j=0; j<N; j++)
          xd_dev[j] *= c[i];
      }
    }
    free(xd_dev_ptrs);
    return(0);
  }

  /* Allocate and store dev pointers to copy to device */
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

  /*
   * Z[i] = c[i] * X[i]
   */
#pragma omp target map(to:N,nvec,c[:nvec],xd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = c[i] * xd_dev[j];
    }
  }
  free(xd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

int N_VConstVectorArray_OpenMPDEV(int nvec, realtype c, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    zd_dev=NULL;
  realtype**   zd_dev_ptrs=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VConst */
  if (nvec == 1) {
    N_VConst_OpenMPDEV(c, Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LENGTH_OMPDEV(Z[0]);

  /* get device */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

  /* set each vector in the vector array to a constant */
#pragma omp target map(to:N,nvec,zd_dev_ptrs[:nvec]) \
  is_device_ptr(zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = c;
    }
  }

  free(zd_dev_ptrs);
  return(0);
}

int N_VWrmsNormVectorArray_OpenMPDEV(int nvec, N_Vector* X, N_Vector* W, realtype* nrm)
{
  int          i, dev;
  sunindextype j, N;
  realtype     sum;
  realtype*    wd_dev=NULL;
  realtype*    xd_dev=NULL;
  realtype**   wd_dev_ptrs=NULL;
  realtype**   xd_dev_ptrs=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNorm_OpenMPDEV(X[0], W[0]);
    return(0);
  }

  /* get vector length */
  N  = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* initialize norms */
  for (i=0; i<nvec; i++)
    nrm[i] = ZERO;

  /* Allocate and store dev pointers to copy to device */
  wd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    wd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(W[i]);
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);

  /* compute the WRMS norm for each vector in the vector array */
#pragma omp target map(to:N,nvec,xd_dev_ptrs[:nvec],wd_dev_ptrs[:nvec]) map(tofrom:nrm[:nvec]) \
  is_device_ptr(xd_dev, wd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      wd_dev = wd_dev_ptrs[i];
      sum = ZERO;
#pragma omp parallel for reduction(+:sum) schedule(static, 1)
      {
        for (j=0; j<N; j++)
          sum += SUNSQR(xd_dev[j] * wd_dev[j]);
      }
      nrm[i] = SUNRsqrt(sum/N);
    }
  }

  free(wd_dev_ptrs);
  free(xd_dev_ptrs);
  return(0);
}


int N_VWrmsNormMaskVectorArray_OpenMPDEV(int nvec, N_Vector* X, N_Vector* W,
                                        N_Vector id, realtype* nrm)
{
  int          i, dev;
  sunindextype j, N;
  realtype     sum;
  realtype*    wd_dev=NULL;
  realtype*    xd_dev=NULL;
  realtype*    idd_dev=NULL;
  realtype**   wd_dev_ptrs=NULL;
  realtype**   xd_dev_ptrs=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNormMask_OpenMPDEV(X[0], W[0], id);
    return(0);
  }

  /* get vector length and mask data array */
  N   = NV_LENGTH_OMPDEV(X[0]);
  idd_dev = NV_DATA_DEV_OMPDEV(id);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* initialize norms */
  for (i=0; i<nvec; i++)
    nrm[i] = ZERO;

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  wd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    wd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(W[i]);

  /* compute the WRMS norm for each vector in the vector array */
#pragma omp target map(to:N,nvec,xd_dev_ptrs[:nvec],wd_dev_ptrs[:nvec]) map(tofrom:nrm[:nvec]) \
  is_device_ptr(idd_dev,xd_dev,wd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      wd_dev = wd_dev_ptrs[i];
      sum = ZERO;
#pragma omp parallel for reduction(+:sum) schedule(static, 1)
      {
        for (j=0; j<N; j++) {
          if (idd_dev[j] > ZERO)
            sum += SUNSQR(xd_dev[j] * wd_dev[j]);
        }
      }
      nrm[i] = SUNRsqrt(sum/N);
    }
  }

  free(xd_dev_ptrs);
  free(wd_dev_ptrs);
  return(0);
}

int N_VScaleAddMultiVectorArray_OpenMPDEV(int nvec, int nsum, realtype* a,
                                          N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  int          i, j, dev;
  sunindextype k, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

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
      N_VLinearSum_OpenMPDEV(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      return(0);
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector *) malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (j=0; j<nsum; j++) {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    retval = N_VScaleAddMulti_OpenMPDEV(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);
    return(retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1) {
    retval = N_VLinearSumVectorArray_OpenMPDEV(nvec, a[0], X, ONE, Y[0], Z[0]);
    return(retval);
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N  = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * nsum * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++) {
    for (j=0; j<nsum; j++)
      yd_dev_ptrs[i * nsum + j] = NV_DATA_DEV_OMPDEV(Y[j][i]);
  }

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
#pragma omp target map(to:N,nvec,nsum,a[:nsum],xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec*nsum]) \
  is_device_ptr(xd_dev, yd_dev) device(dev)
#pragma omp teams distribute
    {
      for (i=0; i<nvec; i++) {
        xd_dev = xd_dev_ptrs[i];
        for (j=0; j<nsum; j++) {
          yd_dev = yd_dev_ptrs[i*nsum+j];
#pragma omp parallel for schedule(static, 1)
          for (k=0; k<N; k++)
            yd_dev[k] += a[j] * xd_dev[k];
        }
      }
    }
    free(xd_dev_ptrs);
    free(yd_dev_ptrs);
    return(0);
  }

  /* Allocate and store dev pointers to copy to device */
  zd_dev_ptrs = (realtype**) malloc(nvec * nsum * sizeof(realtype*));
  for (i=0; i<nvec; i++) {
    for (j=0; j<nsum; j++)
      zd_dev_ptrs[i * nsum + j] = NV_DATA_DEV_OMPDEV(Z[j][i]);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
#pragma omp target map(to:N,nvec,nsum,a[:nsum],xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec*nsum],zd_dev_ptrs[:nvec*nsum]) \
  is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      for (j=0; j<nsum; j++) {
        yd_dev = yd_dev_ptrs[i*nsum+j];
        zd_dev = zd_dev_ptrs[i*nsum+j];
#pragma omp parallel for schedule(static, 1)
        for (k=0; k<N; k++)
          zd_dev[k] = a[j] * xd_dev[k] + yd_dev[k];
      }
    }
  }

  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

int N_VLinearCombinationVectorArray_OpenMPDEV(int nvec, int nsum,
                                             realtype* c,
                                             N_Vector** X,
                                             N_Vector* Z)
{
  int          i; /* vector arrays index in summation [0,nsum) */
  int          j; /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  realtype*    zd_dev=NULL;
  realtype*    xd_dev=NULL;
  realtype**   zd_dev_ptrs=NULL;
  realtype**   xd_dev_ptrs=NULL;
  int dev;

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
      N_VScale_OpenMPDEV(c[0], X[0][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearSum */
    if (nsum == 2) {
      N_VLinearSum_OpenMPDEV(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nsum; i++) {
      Y[i] = X[i][0];
    }

    N_VLinearCombination_OpenMPDEV(nsum, c, Y, Z[0]);

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

    N_VScaleVectorArray_OpenMPDEV(nvec, ctmp, X[0], Z);

    free(ctmp);
    return(0);
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2) {
    N_VLinearSumVectorArray_OpenMPDEV(nvec, c[0], X[0], c[1], X[1], Z);
    return(0);
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LENGTH_OMPDEV(Z[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  xd_dev_ptrs = (realtype**) malloc(nvec * nsum * sizeof(realtype*));
  for (j=0; j<nvec; j++)
    zd_dev_ptrs[j] = NV_DATA_DEV_OMPDEV(Z[j]);
  for (j=0; j<nvec; j++) {
    for (i=0; i<nsum; i++)
      xd_dev_ptrs[j * nsum + i] = NV_DATA_DEV_OMPDEV(X[i][j]);
  }

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE)) {
#pragma omp target map(to:N,nvec,c[:nsum],xd_dev_ptrs[:nvec*nsum],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute
    {
      for (j=0; j<nvec; j++) {
        zd_dev = zd_dev_ptrs[j];
        for (i=1; i<nsum; i++) {
          xd_dev = xd_dev_ptrs[j*nsum+i];
#pragma omp parallel for schedule(static, 1)
          for (k=0; k<N; k++)
            zd_dev[k] += c[i] * xd_dev[k];
        }
      }
    }
    free(xd_dev_ptrs);
    free(zd_dev_ptrs);
    return(0);
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z) {
#pragma omp target map(to:N,nvec,c[:nsum],xd_dev_ptrs[:nvec*nsum],zd_dev_ptrs[:nvec]) \
  is_device_ptr(zd_dev) device(dev)
#pragma omp teams distribute
    {
      for (j=0; j<nvec; j++) {
        zd_dev = zd_dev_ptrs[j];
#pragma omp parallel for schedule(static, 1)
        for (k=0; k<N; k++)
          zd_dev[k] *= c[0];

        for (i=1; i<nsum; i++) {
          xd_dev = xd_dev_ptrs[j*nsum+i];
#pragma omp parallel for schedule(static, 1)
          for (k=0; k<N; k++)
            zd_dev[k] += c[i] * xd_dev[k];
        }
      }
    }
    free(xd_dev_ptrs);
    free(zd_dev_ptrs);
    return(0);
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
#pragma omp target map(to:N,nvec,c[:nsum],xd_dev_ptrs[:nvec*nsum],zd_dev_ptrs[:nvec]) \
  is_device_ptr(zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (j=0; j<nvec; j++) {
      /* scale first vector in the sum into the output vector */
      xd_dev = xd_dev_ptrs[j*nsum];
      zd_dev = zd_dev_ptrs[j];
#pragma omp parallel for schedule(static, 1)
      for (k=0; k<N; k++)
        zd_dev[k] = c[0] * xd_dev[k];

      /* scale and sum remaining vectors into the output vector */
      for (i=1; i<nsum; i++) {
        xd_dev = xd_dev_ptrs[j*nsum+i];
#pragma omp parallel for schedule(static, 1)
        for (k=0; k<N; k++)
          zd_dev[k] += c[i] * xd_dev[k];
      }
    }
  }
  free(xd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}


/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */


/* ----------------------------------------------------------------------------
 * Copy vector components into a second vector
 */

static void VCopy_OpenMPDEV(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *zd_dev;
  int dev;

  xd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = xd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute vector sum
 */

static void VSum_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev, *zd_dev;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = xd_dev[i]+yd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute vector difference
 */

static void VDiff_OpenMPDEV(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev, *zd_dev;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = xd_dev[i]-yd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute the negative of a vector
 */

static void VNeg_OpenMPDEV(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *zd_dev;
  int dev;

  xd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = -xd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector sum
 */

static void VScaleSum_OpenMPDEV(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev, *zd_dev;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = c*(xd_dev[i]+yd_dev[i]);

  return;
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector difference
 */

static void VScaleDiff_OpenMPDEV(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev, *zd_dev;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = c*(xd_dev[i]-yd_dev[i]);

  return;
}


/* ----------------------------------------------------------------------------
 * Compute vector sum z[i] = a*x[i]+y[i]
 */

static void VLin1_OpenMPDEV(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev, *zd_dev;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = (a*xd_dev[i])+yd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute vector difference z[i] = a*x[i]-y[i]
 */

static void VLin2_OpenMPDEV(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev, *zd_dev;
  int dev;

  xd_dev = yd_dev = zd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);
  zd_dev = NV_DATA_DEV_OMPDEV(z);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    zd_dev[i] = (a*xd_dev[i])-yd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute special cases of linear sum
 */

static void Vaxpy_OpenMPDEV(realtype a, N_Vector x, N_Vector y)
{
  sunindextype i, N;
  realtype *xd_dev, *yd_dev;
  int dev;

  xd_dev = yd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);
  yd_dev = NV_DATA_DEV_OMPDEV(y);

  /* get default device identifier */
  dev = omp_get_default_device();

  if (a == ONE) {
#pragma omp target is_device_ptr(xd_dev, yd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
    for (i = 0; i < N; i++)
      yd_dev[i] += xd_dev[i];
    return;
  }

  if (a == -ONE) {
#pragma omp target is_device_ptr(xd_dev, yd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
    for (i = 0; i < N; i++)
      yd_dev[i] -= xd_dev[i];
    return;
  }

#pragma omp target is_device_ptr(xd_dev, yd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    yd_dev[i] += a*xd_dev[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute scaled vector x[i] = a*x[i]
 */

static void VScaleBy_OpenMPDEV(realtype a, N_Vector x)
{
  sunindextype i, N;
  realtype *xd_dev;
  int dev;

  xd_dev = NULL;

  N      = NV_LENGTH_OMPDEV(x);
  xd_dev = NV_DATA_DEV_OMPDEV(x);

  /* get default device identifier */
  dev = omp_get_default_device();

#pragma omp target is_device_ptr(xd_dev) device(dev)
#pragma omp teams distribute parallel for schedule(static, 1)
  for (i = 0; i < N; i++)
    xd_dev[i] *= a;

  return;
}


/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static int VSumVectorArray_OpenMPDEV(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  N = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev, yd_dev, zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      yd_dev = yd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = xd_dev[j] + yd_dev[j];
    }
  }

  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

static int VDiffVectorArray_OpenMPDEV(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  N = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,yd_dev,zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      yd_dev = yd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = xd_dev[j] - yd_dev[j];
    }
  }

  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

static int VScaleSumVectorArray_OpenMPDEV(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  N = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,yd_dev,zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      yd_dev = yd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = c * (xd_dev[j] + yd_dev[j]);
    }
  }

  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

static int VScaleDiffVectorArray_OpenMPDEV(int nvec, realtype c, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  N = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev ointer to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,yd_dev,zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      yd_dev = yd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = c * (xd_dev[j] - yd_dev[j]);
    }
  }

  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

static int VLin1VectorArray_OpenMPDEV(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  N = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,yd_dev,zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      yd_dev = yd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = (a * xd_dev[j]) + yd_dev[j];
    }
  }

  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

static int VLin2VectorArray_OpenMPDEV(int nvec, realtype a, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype*    zd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;
  realtype**   zd_dev_ptrs=NULL;

  N = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  zd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);
  for (i=0; i<nvec; i++)
    zd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Z[i]);

#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec],zd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,yd_dev,zd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
      xd_dev = xd_dev_ptrs[i];
      yd_dev = yd_dev_ptrs[i];
      zd_dev = zd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        zd_dev[j] = (a * xd_dev[j]) - yd_dev[j];
    }
  }

  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  free(zd_dev_ptrs);
  return(0);
}

static int VaxpyVectorArray_OpenMPDEV(int nvec, realtype a, N_Vector* X, N_Vector* Y)
{
  int          i, dev;
  sunindextype j, N;
  realtype*    xd_dev=NULL;
  realtype*    yd_dev=NULL;
  realtype**   xd_dev_ptrs=NULL;
  realtype**   yd_dev_ptrs=NULL;

  N = NV_LENGTH_OMPDEV(X[0]);

  /* get default device identifier */
  dev = omp_get_default_device();

  /* Allocate and store dev pointers to copy to device */
  xd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  yd_dev_ptrs = (realtype**) malloc(nvec * sizeof(realtype*));
  for (i=0; i<nvec; i++)
    xd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(X[i]);
  for (i=0; i<nvec; i++)
    yd_dev_ptrs[i] = NV_DATA_DEV_OMPDEV(Y[i]);

  if (a == ONE) {
#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,yd_dev) device(dev)
#pragma omp teams distribute
    {
      for (i=0; i<nvec; i++) {
        xd_dev = xd_dev_ptrs[i];
        yd_dev = yd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
        for (j=0; j<N; j++)
          yd_dev[j] += xd_dev[j];
      }
    }
    free(xd_dev_ptrs);
    free(yd_dev_ptrs);
    return(0);
  }

  if (a == -ONE) {
#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,yd_dev) device(dev)
#pragma omp teams distribute
    {
      for (i=0; i<nvec; i++) {
        xd_dev = xd_dev_ptrs[i];
        yd_dev = yd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
        for (j=0; j<N; j++)
          yd_dev[j] -= xd_dev[j];
      }
    }
    free(xd_dev_ptrs);
    free(yd_dev_ptrs);
    return(0);
  }

#pragma omp target map(to:N,xd_dev_ptrs[:nvec],yd_dev_ptrs[:nvec]) \
  is_device_ptr(xd_dev,yd_dev) device(dev)
#pragma omp teams distribute
  {
    for (i=0; i<nvec; i++) {
        xd_dev = xd_dev_ptrs[i];
        yd_dev = yd_dev_ptrs[i];
#pragma omp parallel for schedule(static, 1)
      for (j=0; j<N; j++)
        yd_dev[j] += a * xd_dev[j];
    }
  }
  free(xd_dev_ptrs);
  free(yd_dev_ptrs);
  return(0);
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_OpenMPDEV;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_OpenMPDEV;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_OpenMPDEV;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_OpenMPDEV;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_OpenMPDEV;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_OpenMPDEV;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_OpenMPDEV;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_OpenMPDEV;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_OpenMPDEV;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_OpenMPDEV;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = N_VDotProdMultiLocal_OpenMPDEV;
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
    /* disable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_OpenMPDEV;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_OpenMPDEV;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf) {
    v->ops->nvdotprodmulti      = N_VDotProdMulti_OpenMPDEV;
    v->ops->nvdotprodmultilocal = N_VDotProdMulti_OpenMPDEV;
  } else {
    v->ops->nvdotprodmulti      = NULL;
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_OpenMPDEV;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_OpenMPDEV;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_OpenMPDEV;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_OpenMPDEV;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_OpenMPDEV;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_OpenMPDEV;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_OpenMPDEV(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_OpenMPDEV;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}
