/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Cody J. Balos, Daniel McGreer @ LLNL
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
 * This is the implementation file for a RAJA implementation
 * of the NVECTOR package. This will support CUDA and HIP
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <RAJA/RAJA.hpp>
#include <nvector/nvector_raja.h>

#include "sundials_debug.h"

// RAJA defines
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
#include <sunmemory/sunmemory_cuda.h>
#include "sundials_cuda.h"
#define RAJA_NODE_TYPE RAJA::cuda_exec< 256 >
#define RAJA_REDUCE_TYPE RAJA::cuda_reduce
#define SUNDIALS_GPU_PREFIX(val) cuda ## val
#define SUNDIALS_GPU_VERIFY SUNDIALS_CUDA_VERIFY
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
#include <sunmemory/sunmemory_hip.h>
#include "sundials_hip.h"
#define RAJA_NODE_TYPE RAJA::hip_exec< 512 >
#define RAJA_REDUCE_TYPE RAJA::hip_reduce
#define SUNDIALS_GPU_PREFIX(val) hip ## val
#define SUNDIALS_GPU_VERIFY SUNDIALS_HIP_VERIFY
#endif

#define RAJA_LAMBDA [=] __device__

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

// Static constants
static constexpr sunindextype zeroIdx = 0;

// Helpful macros
#define NVEC_RAJA_CONTENT(x) ((N_VectorContent_Raja)(x->content))
#define NVEC_RAJA_PRIVATE(x) ((N_PrivateVectorContent_Raja)(NVEC_RAJA_CONTENT(x)->priv))
#define NVEC_RAJA_MEMSIZE(x) (NVEC_RAJA_CONTENT(x)->length * sizeof(realtype))
#define NVEC_RAJA_MEMHELP(x) (NVEC_RAJA_CONTENT(x)->mem_helper)
#define NVEC_RAJA_HDATAp(x)  ((realtype*) NVEC_RAJA_CONTENT(x)->host_data->ptr)
#define NVEC_RAJA_DDATAp(x)  ((realtype*) NVEC_RAJA_CONTENT(x)->device_data->ptr)

/*
 * Private structure definition
 */

struct _N_PrivateVectorContent_Raja
{
  booleantype use_managed_mem; /* indicates if the data pointers and buffer pointers are managed memory */
};

typedef struct _N_PrivateVectorContent_Raja *N_PrivateVectorContent_Raja;


/*
 * Utility functions
 */

static int AllocateData(N_Vector v);
static void CreateArrayOfPointersOnDevice(realtype*** d_ptrs, SUNMemory* d_ref,
                                          int nvec, N_Vector *V);
static void Create2DArrayOfPointersOnDevice(realtype*** d_ptrs, SUNMemory* d_ref,
                                            int nvec, int nsum, N_Vector **V);

N_Vector N_VNewEmpty_Raja()
{
  N_Vector v;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty();
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = N_VGetVectorID_Raja;
  v->ops->nvclone                 = N_VClone_Raja;
  v->ops->nvcloneempty            = N_VCloneEmpty_Raja;
  v->ops->nvdestroy               = N_VDestroy_Raja;
  v->ops->nvspace                 = N_VSpace_Raja;
  v->ops->nvgetlength             = N_VGetLength_Raja;
  v->ops->nvgetarraypointer       = N_VGetHostArrayPointer_Raja;
  v->ops->nvgetdevicearraypointer = N_VGetDeviceArrayPointer_Raja;
  v->ops->nvsetarraypointer       = N_VSetHostArrayPointer_Raja;


  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Raja;
  v->ops->nvconst        = N_VConst_Raja;
  v->ops->nvprod         = N_VProd_Raja;
  v->ops->nvdiv          = N_VDiv_Raja;
  v->ops->nvscale        = N_VScale_Raja;
  v->ops->nvabs          = N_VAbs_Raja;
  v->ops->nvinv          = N_VInv_Raja;
  v->ops->nvaddconst     = N_VAddConst_Raja;
  v->ops->nvdotprod      = N_VDotProd_Raja;
  v->ops->nvmaxnorm      = N_VMaxNorm_Raja;
  v->ops->nvmin          = N_VMin_Raja;
  v->ops->nvl1norm       = N_VL1Norm_Raja;
  v->ops->nvinvtest      = N_VInvTest_Raja;
  v->ops->nvconstrmask   = N_VConstrMask_Raja;
  v->ops->nvminquotient  = N_VMinQuotient_Raja;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Raja;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Raja;
  v->ops->nvwl2norm      = N_VWL2Norm_Raja;
  v->ops->nvcompare      = N_VCompare_Raja;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Raja;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Raja;
  v->ops->nvdotprodlocal     = N_VDotProd_Raja;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Raja;
  v->ops->nvminlocal         = N_VMin_Raja;
  v->ops->nvl1normlocal      = N_VL1Norm_Raja;
  v->ops->nvinvtestlocal     = N_VInvTest_Raja;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Raja;
  v->ops->nvminquotientlocal = N_VMinQuotient_Raja;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_Raja;
  v->ops->nvbufpack   = N_VBufPack_Raja;
  v->ops->nvbufunpack = N_VBufUnpack_Raja;

  /* print operation for debugging */
  v->ops->nvprint            = N_VPrint_Raja;
  v->ops->nvprintfile        = N_VPrintFile_Raja;

  v->content = (N_VectorContent_Raja) malloc(sizeof(_N_VectorContent_Raja));
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  NVEC_RAJA_CONTENT(v)->priv = malloc(sizeof(_N_PrivateVectorContent_Raja));
  if (NVEC_RAJA_CONTENT(v)->priv == NULL)
  {
    N_VDestroy(v);
    return NULL;
  }

  NVEC_RAJA_CONTENT(v)->length          = 0;
  NVEC_RAJA_CONTENT(v)->mem_helper      = NULL;
  NVEC_RAJA_CONTENT(v)->own_helper      = SUNFALSE;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNFALSE;

  return(v);
}

N_Vector N_VNew_Raja(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda();
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Hip();
#endif
  NVEC_RAJA_CONTENT(v)->own_helper      = SUNTRUE;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNFALSE;

  if (NVEC_RAJA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Raja: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Raja: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return(v);
}

N_Vector N_VNewWithMemHelp_Raja(sunindextype length, booleantype use_managed_mem, SUNMemoryHelper helper)
{
  N_Vector v;

  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Raja: helper is NULL\n");
    return(NULL);
  }

  if (!SUNMemoryHelper_ImplementsRequiredOps(helper))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Raja: helper doesn't implement all required ops\n");
    return(NULL);
  }

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->mem_helper      = helper;
  NVEC_RAJA_CONTENT(v)->own_helper      = SUNFALSE;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = use_managed_mem;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Raja: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VNewManaged_Raja(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda();
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Hip();
#endif
  NVEC_RAJA_CONTENT(v)->own_helper      = SUNTRUE;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNTRUE;

  if (NVEC_RAJA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Raja: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Raja: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return(v);
}

N_Vector N_VMake_Raja(sunindextype length, realtype *h_vdata, realtype *d_vdata)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->host_data       = SUNMemoryHelper_Wrap(h_vdata, SUNMEMTYPE_HOST);
  NVEC_RAJA_CONTENT(v)->device_data     = SUNMemoryHelper_Wrap(d_vdata, SUNMEMTYPE_DEVICE);
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda();
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Hip();
#endif
  NVEC_RAJA_CONTENT(v)->own_helper      = SUNTRUE;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNFALSE;

  if (NVEC_RAJA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Raja: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }


  if (NVEC_RAJA_CONTENT(v)->device_data == NULL ||
      NVEC_RAJA_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Raja: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VMakeManaged_Raja(sunindextype length, realtype *vdata)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->host_data       = SUNMemoryHelper_Wrap(vdata, SUNMEMTYPE_UVM);
  NVEC_RAJA_CONTENT(v)->device_data     = SUNMemoryHelper_Alias(NVEC_RAJA_CONTENT(v)->host_data);
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda();
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Hip();
#endif
  NVEC_RAJA_CONTENT(v)->own_helper      = SUNTRUE;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = SUNTRUE;

  if (NVEC_RAJA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Raja: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (NVEC_RAJA_CONTENT(v)->device_data == NULL ||
      NVEC_RAJA_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Raja: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

/* -----------------------------------------------------------------
 * Function to return the global length of the vector.
 * This is defined as an inline function in nvector_raja.h, so
 * we just mark it as extern here.
 */
extern sunindextype N_VGetLength_Raja(N_Vector v);

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data.
 * This is defined as an inline function in nvector_raja.h, so
 * we just mark it as extern here.
 */

extern realtype *N_VGetHostArrayPointer_Raja(N_Vector x);

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data.
 * This is defined as an inline function in nvector_raja.h, so
 * we just mark it as extern here.
 */

extern realtype *N_VGetDeviceArrayPointer_Raja(N_Vector x);


/* ----------------------------------------------------------------------------
 * Set pointer to the raw host data. Does not free the existing pointer.
 */

void N_VSetHostArrayPointer_Raja(realtype* h_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Raja(v))
  {
    if (NVEC_RAJA_CONTENT(v)->host_data)
    {
      NVEC_RAJA_CONTENT(v)->host_data->ptr = (void*) h_vdata;
      NVEC_RAJA_CONTENT(v)->device_data->ptr = (void*) h_vdata;
    }
    else
    {
      NVEC_RAJA_CONTENT(v)->host_data = SUNMemoryHelper_Wrap((void*) h_vdata, SUNMEMTYPE_UVM);
      NVEC_RAJA_CONTENT(v)->device_data = SUNMemoryHelper_Alias(NVEC_RAJA_CONTENT(v)->host_data);
    }
  }
  else
  {
    if (NVEC_RAJA_CONTENT(v)->host_data)
    {
      NVEC_RAJA_CONTENT(v)->host_data->ptr = (void*) h_vdata;
    }
    else
    {
      NVEC_RAJA_CONTENT(v)->host_data = SUNMemoryHelper_Wrap((void*) h_vdata, SUNMEMTYPE_HOST);
    }
  }
}

/* ----------------------------------------------------------------------------
 * Set pointer to the raw device data
 */

void N_VSetDeviceArrayPointer_Raja(realtype* d_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Raja(v))
  {
    if (NVEC_RAJA_CONTENT(v)->device_data)
    {
      NVEC_RAJA_CONTENT(v)->device_data->ptr = (void*) d_vdata;
      NVEC_RAJA_CONTENT(v)->host_data->ptr = (void*) d_vdata;
    }
    else
    {
      NVEC_RAJA_CONTENT(v)->device_data = SUNMemoryHelper_Wrap((void*) d_vdata, SUNMEMTYPE_UVM);
      NVEC_RAJA_CONTENT(v)->host_data = SUNMemoryHelper_Alias(NVEC_RAJA_CONTENT(v)->device_data);
    }
  }
  else
  {
    if (NVEC_RAJA_CONTENT(v)->device_data)
    {
      NVEC_RAJA_CONTENT(v)->device_data->ptr = (void*) d_vdata;
    }
    else
    {
      NVEC_RAJA_CONTENT(v)->device_data = SUNMemoryHelper_Wrap((void*) d_vdata, SUNMEMTYPE_DEVICE);
    }
  }
}

/* ----------------------------------------------------------------------------
 * Return a flag indicating if the memory for the vector data is managed
 */
booleantype N_VIsManagedMemory_Raja(N_Vector x)
{
  return NVEC_RAJA_PRIVATE(x)->use_managed_mem;
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Raja(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(x),
                                        NVEC_RAJA_CONTENT(x)->device_data,
                                        NVEC_RAJA_CONTENT(x)->host_data,
                                        NVEC_RAJA_MEMSIZE(x),
                                        0);

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyToDevice_Raja: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but on the default stream currently */
  SUNDIALS_GPU_VERIFY(SUNDIALS_GPU_PREFIX(StreamSynchronize)(0));
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Raja(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(x),
                                        NVEC_RAJA_CONTENT(x)->host_data,
                                        NVEC_RAJA_CONTENT(x)->device_data,
                                        NVEC_RAJA_MEMSIZE(x),
                                        0);

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyFromDevice_Raja: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  SUNDIALS_GPU_VERIFY(SUNDIALS_GPU_PREFIX(StreamSynchronize)(0));
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

void N_VPrint_Raja(N_Vector X)
{
  N_VPrintFile_Raja(X, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void N_VPrintFile_Raja(N_Vector X, FILE *outfile)
{
  sunindextype i;

  for (i = 0; i < NVEC_RAJA_CONTENT(X)->length; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", NVEC_RAJA_HDATAp(X)[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", NVEC_RAJA_HDATAp(X)[i]);
#else
    fprintf(outfile, "%11.8g\n", NVEC_RAJA_HDATAp(X)[i]);
#endif
  }
  fprintf(outfile, "\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Raja(N_Vector w)
{
  N_Vector v;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Set content */
  NVEC_RAJA_CONTENT(v)->length          = NVEC_RAJA_CONTENT(w)->length;
  NVEC_RAJA_CONTENT(v)->host_data       = NULL;
  NVEC_RAJA_CONTENT(v)->device_data     = NULL;
  NVEC_RAJA_PRIVATE(v)->use_managed_mem = NVEC_RAJA_PRIVATE(w)->use_managed_mem;


  return(v);
}

N_Vector N_VClone_Raja(N_Vector w)
{
  N_Vector v;
  v = NULL;
  v = N_VCloneEmpty_Raja(w);
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->mem_helper = SUNMemoryHelper_Clone(NVEC_RAJA_MEMHELP(w));
  NVEC_RAJA_CONTENT(v)->own_helper = SUNTRUE;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Raja: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

return(v);

}


void N_VDestroy_Raja(N_Vector v)
{
  N_VectorContent_Raja vc;
  N_PrivateVectorContent_Raja vcp;

  if (v == NULL) return;

  /* free ops structure */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }

  /* extract content */
  vc = NVEC_RAJA_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

  /* free private content */
  vcp = (N_PrivateVectorContent_Raja) vc->priv;
  if (vcp != NULL)
  {
    /* free items in private content */
    free(vcp);
    vc->priv = NULL;
  }

  /* free items in content */
  if (NVEC_RAJA_MEMHELP(v))
  {
    SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(v), vc->host_data);
    vc->host_data = NULL;
    SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(v), vc->device_data);
    vc->device_data = NULL;
    if (vc->own_helper) SUNMemoryHelper_Destroy(vc->mem_helper);
    vc->mem_helper = NULL;
  }
  else
  {
    SUNDIALS_DEBUG_PRINT("WARNING in N_VDestroy_Raja: mem_helper was NULL when trying to dealloc data, this could result in a memory leak\n");
  }

  /* free content struct */
  free(vc);

  /* free vector */
  free(v);

  return;
}

void N_VSpace_Raja(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NVEC_RAJA_CONTENT(X)->length;
  *liw = 2;
}

void N_VConst_Raja(realtype c, N_Vector Z)
{
  const sunindextype N = NVEC_RAJA_CONTENT(Z)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N), RAJA_LAMBDA(sunindextype i) {
     zdata[i] = c;
  });
}

void N_VLinearSum_Raja(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const realtype *ydata = NVEC_RAJA_DDATAp(Y);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = a*xdata[i] + b*ydata[i];
    }
  );
}

void N_VProd_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const realtype *ydata = NVEC_RAJA_DDATAp(Y);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] * ydata[i];
    }
  );
}

void N_VDiv_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const realtype *ydata = NVEC_RAJA_DDATAp(Y);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] / ydata[i];
    }
  );
}

void N_VScale_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = c * xdata[i];
    }
  );
}

void N_VAbs_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = abs(xdata[i]);
    }
  );
}

void N_VInv_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = ONE / xdata[i];
    }
  );
}

void N_VAddConst_Raja(N_Vector X, realtype b, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] + b;
    }
  );
}

realtype N_VDotProd_Raja(N_Vector X, N_Vector Y)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const realtype *ydata = NVEC_RAJA_DDATAp(Y);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += xdata[i] * ydata[i] ;
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VMaxNorm_Raja(N_Vector X)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceMax< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result.max(abs(xdata[i]));
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWSqrSumLocal_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const realtype *wdata = NVEC_RAJA_DDATAp(W);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWrmsNorm_Raja(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Raja(X, W);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  return std::sqrt(sum/N);
}

realtype N_VWSqrSumMaskLocal_Raja(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const realtype *wdata = NVEC_RAJA_DDATAp(W);
  const realtype *iddata = NVEC_RAJA_DDATAp(ID);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (iddata[i] > ZERO)
        gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWrmsNormMask_Raja(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype sum = N_VWSqrSumMaskLocal_Raja(X, W, ID);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  return std::sqrt(sum/N);
}

realtype N_VMin_Raja(N_Vector X)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceMin< RAJA_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result.min(xdata[i]);
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VWL2Norm_Raja(N_Vector X, N_Vector W)
{
  return std::sqrt(N_VWSqrSumLocal_Raja(X, W));
}

realtype N_VL1Norm_Raja(N_Vector X)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (abs(xdata[i]));
    }
  );

  return (static_cast<realtype>(gpu_result));
}

void N_VCompare_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
    }
  );
}

booleantype N_VInvTest_Raja(N_Vector x, N_Vector z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(x);
  const sunindextype N = NVEC_RAJA_CONTENT(x)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(z);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(ZERO);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (xdata[i] == ZERO) {
        gpu_result += ONE;
      } else {
        zdata[i] = ONE/xdata[i];
      }
    }
  );
  realtype minimum = static_cast<realtype>(gpu_result);
  return (minimum < HALF);
}

booleantype N_VConstrMask_Raja(N_Vector c, N_Vector x, N_Vector m)
{
  const realtype *cdata = NVEC_RAJA_DDATAp(c);
  const realtype *xdata = NVEC_RAJA_DDATAp(x);
  const sunindextype N = NVEC_RAJA_CONTENT(x)->length;
  realtype *mdata = NVEC_RAJA_DDATAp(m);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(ZERO);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      bool test = (abs(cdata[i]) > ONEPT5 && cdata[i]*xdata[i] <= ZERO) ||
                  (abs(cdata[i]) > HALF   && cdata[i]*xdata[i] <  ZERO);
      mdata[i] = test ? ONE : ZERO;
      gpu_result += mdata[i];
    }
  );

  realtype sum = static_cast<realtype>(gpu_result);
  return(sum < HALF);
}

realtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom)
{
  const realtype *ndata = NVEC_RAJA_DDATAp(num);
  const realtype *ddata = NVEC_RAJA_DDATAp(denom);
  const sunindextype N = NVEC_RAJA_CONTENT(num)->length;

  RAJA::ReduceMin< RAJA_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (ddata[i] != ZERO)
        gpu_result.min(ndata[i]/ddata[i]);
    }
  );
  return (static_cast<realtype>(gpu_result));
}


/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearCombination_Raja(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int retval;

  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(X[0]);
  sunindextype N = NVEC_RAJA_CONTENT(z)->length;
  realtype* d_zd = NVEC_RAJA_DDATAp(z);

  // Create device c array for device
  SUNMemory h_c, d_c;
  h_c = SUNMemoryHelper_Wrap(c, SUNMEMTYPE_HOST);
  retval = SUNMemoryHelper_Alloc(h, &d_c, sizeof(realtype)*nvec, SUNMEMTYPE_DEVICE);
  if (retval) return(-1);

  // Copy c array to device
  retval = SUNMemoryHelper_Copy(h, d_c, h_c, sizeof(realtype)*nvec);
  if (retval) return(-1);

  SUNMemory d_X;
  realtype** d_Xd;
  CreateArrayOfPointersOnDevice(&d_Xd, &d_X, nvec, X);

  // Shortcut to the arrays to work on
  realtype* d_cd  = (realtype*) d_c->ptr;
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      d_zd[i] = d_cd[0] * d_Xd[0][i];
      for (int j=1; j<nvec; j++)
        d_zd[i] += d_cd[j] * d_Xd[j][i];
    }
  );

  SUNMemoryHelper_Dealloc(h, h_c);
  SUNMemoryHelper_Dealloc(h, d_c);
  SUNMemoryHelper_Dealloc(h, d_X);

  return(0);
}


int N_VScaleAddMulti_Raja(int nvec, realtype* c, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  int retval;

  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(x);
  sunindextype N = NVEC_RAJA_CONTENT(x)->length;
  realtype* d_xd = NVEC_RAJA_DDATAp(x);

  // Create c array for device
  SUNMemory h_c, d_c;
  h_c = SUNMemoryHelper_Wrap(c, SUNMEMTYPE_HOST);
  retval = SUNMemoryHelper_Alloc(h, &d_c, sizeof(realtype)*nvec, SUNMEMTYPE_DEVICE);
  if (retval) return(-1);

  // Copy c array to device
  retval = SUNMemoryHelper_Copy(h, d_c, h_c, sizeof(realtype)*nvec);
  if (retval) return(-1);

  SUNMemory d_Y, d_Z;
  realtype **d_Yd, **d_Zd;
  CreateArrayOfPointersOnDevice(&d_Yd, &d_Y, nvec, Y);
  CreateArrayOfPointersOnDevice(&d_Zd, &d_Z, nvec, Z);

  // Perform operation
  realtype* d_cd  = (realtype*) d_c->ptr;
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
     RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = d_cd[j] * d_xd[i] + d_Yd[j][i];
    }
  );

  SUNMemoryHelper_Dealloc(h, h_c);
  SUNMemoryHelper_Dealloc(h, d_c);
  SUNMemoryHelper_Dealloc(h, d_Y);
  SUNMemoryHelper_Dealloc(h, d_Z);

  return(0);
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Raja(int nvec,
                                 realtype a, N_Vector* X,
                                 realtype b, N_Vector* Y,
                                 N_Vector* Z)
{
  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(Z[0]);
  sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  SUNMemory d_X, d_Y, d_Z;
  realtype **d_Xd, **d_Yd, **d_Zd;
  CreateArrayOfPointersOnDevice(&d_Xd, &d_X, nvec, X);
  CreateArrayOfPointersOnDevice(&d_Yd, &d_Y, nvec, Y);
  CreateArrayOfPointersOnDevice(&d_Zd, &d_Z, nvec, Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = a * d_Xd[j][i] + b * d_Yd[j][i];
    }
  );

  SUNMemoryHelper_Dealloc(h, d_X);
  SUNMemoryHelper_Dealloc(h, d_Y);
  SUNMemoryHelper_Dealloc(h, d_Z);

  return(0);
}


int N_VScaleVectorArray_Raja(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int retval;

  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(Z[0]);
  sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Create c array for device
  SUNMemory h_c, d_c;
  h_c = SUNMemoryHelper_Wrap(c, SUNMEMTYPE_HOST);
  retval = SUNMemoryHelper_Alloc(h, &d_c, sizeof(realtype)*nvec, SUNMEMTYPE_DEVICE);
  if (retval) return(-1);

  // Copy c array to device
  retval = SUNMemoryHelper_Copy(h, d_c, h_c, sizeof(realtype)*nvec);
  if (retval) return(-1);

  SUNMemory d_X, d_Z;
  realtype **d_Xd, **d_Zd;
  CreateArrayOfPointersOnDevice(&d_Xd, &d_X, nvec, X);
  CreateArrayOfPointersOnDevice(&d_Zd, &d_Z, nvec, Z);

  realtype* d_cd  = (realtype*) d_c->ptr;
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = d_cd[j] * d_Xd[j][i];
    }
  );

  SUNMemoryHelper_Dealloc(h, h_c);
  SUNMemoryHelper_Dealloc(h, d_c);
  SUNMemoryHelper_Dealloc(h, d_X);
  SUNMemoryHelper_Dealloc(h, d_Z);

  return(0);
}


int N_VConstVectorArray_Raja(int nvec, realtype c, N_Vector* Z)
{
  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(Z[0]);
  sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  SUNMemory d_Z;
  realtype** d_Zd;
  CreateArrayOfPointersOnDevice(&d_Zd, &d_Z, nvec, Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = c;
    }
  );

  SUNMemoryHelper_Dealloc(h, d_Z);

  return(0);
}


int N_VScaleAddMultiVectorArray_Raja(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  int retval;

  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(X[0]);
  sunindextype N = NVEC_RAJA_CONTENT(X[0])->length;

  // Create c array for device
  SUNMemory h_c, d_c;
  h_c = SUNMemoryHelper_Wrap(c, SUNMEMTYPE_HOST);
  retval = SUNMemoryHelper_Alloc(h, &d_c, sizeof(realtype)*nsum, SUNMEMTYPE_DEVICE);
  if (retval) return(-1);

  // Copy c array to device
  retval = SUNMemoryHelper_Copy(h, d_c, h_c, sizeof(realtype)*nsum);
  if (retval) return(-1);

  SUNMemory d_X, d_Y, d_Z;
  realtype **d_Xd, **d_Yd, **d_Zd;
  CreateArrayOfPointersOnDevice(&d_Xd, &d_X, nvec, X);
  Create2DArrayOfPointersOnDevice(&d_Yd, &d_Y, nvec, nsum, Y);
  Create2DArrayOfPointersOnDevice(&d_Zd, &d_Z, nvec, nsum, Z);

  realtype* d_cd = (realtype*) d_c->ptr;
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        for (int k=0; k<nsum; k++)
          d_Zd[j*nsum+k][i] = d_cd[k] * d_Xd[j][i] + d_Yd[j*nsum+k][i];
    }
  );

  SUNMemoryHelper_Dealloc(h, h_c);
  SUNMemoryHelper_Dealloc(h, d_c);
  SUNMemoryHelper_Dealloc(h, d_X);
  SUNMemoryHelper_Dealloc(h, d_Y);
  SUNMemoryHelper_Dealloc(h, d_Z);

  return(0);
}


int N_VLinearCombinationVectorArray_Raja(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  int retval;

  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(Z[0]);
  sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Create c array for device
  SUNMemory h_c, d_c;
  h_c = SUNMemoryHelper_Wrap(c, SUNMEMTYPE_HOST);
  retval = SUNMemoryHelper_Alloc(h, &d_c, sizeof(realtype)*nsum, SUNMEMTYPE_DEVICE);
  if (retval) return(-1);

  // Copy c array to device
  retval = SUNMemoryHelper_Copy(h, d_c, h_c, sizeof(realtype)*nsum);
  if (retval) return(-1);

  SUNMemory d_X, d_Z;
  realtype **d_Xd, **d_Zd;
  CreateArrayOfPointersOnDevice(&d_Zd, &d_Z, nvec, Z);
  Create2DArrayOfPointersOnDevice(&d_Xd, &d_X, nvec, nsum, X);

  realtype *d_cd = (realtype*) d_c->ptr;
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++) {
        d_Zd[j][i] = d_cd[0] * d_Xd[j*nsum][i];
        for (int k=1; k<nsum; k++) {
          d_Zd[j][i] += d_cd[k] * d_Xd[j*nsum+k][i];
        }
      }
    }
  );

  SUNMemoryHelper_Dealloc(h, h_c);
  SUNMemoryHelper_Dealloc(h, d_c);
  SUNMemoryHelper_Dealloc(h, d_X);
  SUNMemoryHelper_Dealloc(h, d_Z);

  return(0);
}


/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */


int N_VBufSize_Raja(N_Vector x, sunindextype *size)
{
  if (x == NULL) return(-1);
  *size = (sunindextype)NVEC_RAJA_MEMSIZE(x);
  return(0);
}


int N_VBufPack_Raja(N_Vector x, void *buf)
{
  int copy_fail = 0;
  SUNDIALS_GPU_PREFIX(Error_t) cuerr;

  if (x == NULL || buf == NULL) return(-1);

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return(-1);

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(x),
                                        buf_mem,
                                        NVEC_RAJA_CONTENT(x)->device_data,
                                        NVEC_RAJA_MEMSIZE(x),
                                        0);

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = SUNDIALS_GPU_PREFIX(StreamSynchronize)(0);

  SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(x), buf_mem);

  return (!SUNDIALS_GPU_VERIFY(cuerr) || copy_fail ? -1 : 0);
}


int N_VBufUnpack_Raja(N_Vector x, void *buf)
{
  int copy_fail = 0;
  SUNDIALS_GPU_PREFIX(Error_t) cuerr;

  if (x == NULL || buf == NULL) return(-1);

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return(-1);

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(x),
                                        NVEC_RAJA_CONTENT(x)->device_data,
                                        buf_mem,
                                        NVEC_RAJA_MEMSIZE(x),
                                        0);

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = SUNDIALS_GPU_PREFIX(StreamSynchronize)(0);

  SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(x), buf_mem);

  return (!SUNDIALS_GPU_VERIFY(cuerr) || copy_fail ? -1 : 0);
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Raja;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Raja;
    v->ops->nvdotprodmulti      = NULL;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Raja;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Raja;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Raja;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Raja;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Raja;
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

int N_VEnableLinearCombination_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Raja;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Raja;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Raja;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Raja;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Raja;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Raja;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Raja;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}


/*
 * -----------------------------------------------------------------
 * Private utility functions
 * -----------------------------------------------------------------
 */

int AllocateData(N_Vector v)
{
  int alloc_fail = 0;
  N_VectorContent_Raja vc = NVEC_RAJA_CONTENT(v);
  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

  if (N_VGetLength_Raja(v) == 0) return(0);

  if (vcp->use_managed_mem)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v), &(vc->device_data),
                                       NVEC_RAJA_MEMSIZE(v), SUNMEMTYPE_UVM);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed for SUNMEMTYPE_UVM\n");
    }
    vc->host_data = SUNMemoryHelper_Alias(vc->device_data);
  }
  else
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v), &(vc->host_data),
                                       NVEC_RAJA_MEMSIZE(v), SUNMEMTYPE_HOST);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
    }

    alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v), &(vc->device_data),
                                       NVEC_RAJA_MEMSIZE(v), SUNMEMTYPE_DEVICE);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  return(alloc_fail ? -1 : 0);
}


void CreateArrayOfPointersOnDevice(realtype*** d_ptrs, SUNMemory* d_ref,
                                   int nvec, N_Vector *V)
{
  size_t bytes = sizeof(realtype*)*nvec;
  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(V[0]);

  // Default return values
  *d_ref  = nullptr;
  *d_ptrs = nullptr;

  // Create space for host and device pointers
  SUNMemory h_mem;
  SUNMemoryHelper_Alloc(h, &h_mem, bytes, SUNMEMTYPE_HOST);

  SUNMemory d_mem;
  SUNMemoryHelper_Alloc(h, &d_mem, bytes, SUNMEMTYPE_DEVICE);

  // Fill the host memory with the pointers
  realtype** h_array = (realtype**) h_mem->ptr;
  for (int j=0; j<nvec; j++) {
    h_array[j] = NVEC_RAJA_DDATAp(V[j]);
  }

  // Copy the host memory to the device
  SUNMemoryHelper_Copy(h, d_mem, h_mem, bytes);

  // Return the device SUNMemory, and the raw pointer array
  *d_ref  = d_mem;
  *d_ptrs = (realtype**) d_mem->ptr;

  // Free the host SUNMemory
  SUNMemoryHelper_Dealloc(h, h_mem);
}

void Create2DArrayOfPointersOnDevice(realtype*** d_ptrs, SUNMemory* d_ref,
                                     int nvec, int nsum, N_Vector **V)
{
  size_t bytes = sizeof(realtype*)*nsum*nvec;
  SUNMemoryHelper h = NVEC_RAJA_MEMHELP(V[0][0]);

  // Default return values
  *d_ref  = nullptr;
  *d_ptrs = nullptr;

  // Create space for host and device pointers
  SUNMemory h_mem;
  SUNMemoryHelper_Alloc(h, &h_mem, bytes, SUNMEMTYPE_HOST);

  SUNMemory d_mem;
  SUNMemoryHelper_Alloc(h, &d_mem, bytes, SUNMEMTYPE_DEVICE);

  // Fill the host memory with the pointers
  realtype** h_array = (realtype**) h_mem->ptr;
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_array[j*nsum+k] = NVEC_RAJA_DDATAp(V[k][j]);

  // Copy the host memory to the device
  SUNMemoryHelper_Copy(h, d_mem, h_mem, bytes);

  // Return the device SUNMemory, and the raw pointer array
  *d_ref  = d_mem;
  *d_ptrs = (realtype**) d_mem->ptr;

  // Free the host SUNMemory
  SUNMemoryHelper_Dealloc(h, h_mem);
}

} // extern "C"
