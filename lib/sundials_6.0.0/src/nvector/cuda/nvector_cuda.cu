/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, and Cody J. Balos @ LLNL
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
 * This is the implementation file for a CUDA implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <iostream>

#include <nvector/nvector_cuda.h>
#include "VectorArrayKernels.cuh"
#include "VectorKernels.cuh"

#include "sundials_cuda.h"
#include "sundials_debug.h"

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)

extern "C" {

using namespace sundials;
using namespace sundials::cuda;
using namespace sundials::nvector_cuda;

/*
 * Macro definitions
 */

// Macros to access vector content
#define NVEC_CUDA_CONTENT(x)  ((N_VectorContent_Cuda)(x->content))
#define NVEC_CUDA_MEMSIZE(x)  (NVEC_CUDA_CONTENT(x)->length * sizeof(realtype))
#define NVEC_CUDA_MEMHELP(x)  (NVEC_CUDA_CONTENT(x)->mem_helper)
#define NVEC_CUDA_HDATAp(x)   ((realtype*) NVEC_CUDA_CONTENT(x)->host_data->ptr)
#define NVEC_CUDA_DDATAp(x)   ((realtype*) NVEC_CUDA_CONTENT(x)->device_data->ptr)
#define NVEC_CUDA_STREAM(x)   (NVEC_CUDA_CONTENT(x)->stream_exec_policy->stream())

// Macros to access vector private content
#define NVEC_CUDA_PRIVATE(x)  ((N_PrivateVectorContent_Cuda)(NVEC_CUDA_CONTENT(x)->priv))
#define NVEC_CUDA_HBUFFERp(x) ((realtype*) NVEC_CUDA_PRIVATE(x)->reduce_buffer_host->ptr)
#define NVEC_CUDA_DBUFFERp(x) ((realtype*) NVEC_CUDA_PRIVATE(x)->reduce_buffer_dev->ptr)

/*
 * Private structure definition
 */

struct _N_PrivateVectorContent_Cuda
{
  booleantype use_managed_mem; /* do data pointers use managed memory */

  // reduction workspace
  SUNMemory reduce_buffer_dev;   // device memory for reductions
  SUNMemory reduce_buffer_host;  // host memory for reductions
  size_t    reduce_buffer_bytes; // current size of reduction buffers

  // fused op workspace
  SUNMemory fused_buffer_dev;    // device memory for fused ops
  SUNMemory fused_buffer_host;   // host memory for fused ops
  size_t    fused_buffer_bytes;  // current size of the buffers
  size_t    fused_buffer_offset; // current offset into the buffer
};

typedef struct _N_PrivateVectorContent_Cuda *N_PrivateVectorContent_Cuda;

/*
 * Private function definitions
 */

// Allocate vector data
static int AllocateData(N_Vector v);

// Reduction buffer functions
static int InitializeReductionBuffer(N_Vector v, const realtype* value,
                                     size_t n = 1);
static void FreeReductionBuffer(N_Vector v);
static int CopyReductionBufferFromDevice(N_Vector v, size_t n = 1);

// Fused operation buffer functions
static int FusedBuffer_Init(N_Vector v, int nreal, int nptr);
static int FusedBuffer_CopyRealArray(N_Vector v, realtype *r_data, int nval,
                                     realtype **shortcut);
static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector *X, int nvec,
                                      realtype ***shortcut);
static int FusedBuffer_CopyPtrArray2D(N_Vector v, N_Vector **X, int nvec,
                                      int nsum, realtype ***shortcut);
static int FusedBuffer_CopyToDevice(N_Vector v);
static int FusedBuffer_Free(N_Vector v);

// Kernel launch parameters
static int GetKernelParameters(N_Vector v, booleantype reduction, size_t& grid, size_t& block,
                               size_t& shMemSize, cudaStream_t& stream, size_t n = 0);
static void PostKernelLaunch();

N_Vector N_VNewEmpty_Cuda(SUNContext sunctx)
{
  N_Vector v;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = N_VGetVectorID_Cuda;
  v->ops->nvclone                 = N_VClone_Cuda;
  v->ops->nvcloneempty            = N_VCloneEmpty_Cuda;
  v->ops->nvdestroy               = N_VDestroy_Cuda;
  v->ops->nvspace                 = N_VSpace_Cuda;
  v->ops->nvgetlength             = N_VGetLength_Cuda;
  v->ops->nvgetarraypointer       = N_VGetHostArrayPointer_Cuda;
  v->ops->nvgetdevicearraypointer = N_VGetDeviceArrayPointer_Cuda;
  v->ops->nvsetarraypointer       = N_VSetHostArrayPointer_Cuda;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Cuda;
  v->ops->nvconst        = N_VConst_Cuda;
  v->ops->nvprod         = N_VProd_Cuda;
  v->ops->nvdiv          = N_VDiv_Cuda;
  v->ops->nvscale        = N_VScale_Cuda;
  v->ops->nvabs          = N_VAbs_Cuda;
  v->ops->nvinv          = N_VInv_Cuda;
  v->ops->nvaddconst     = N_VAddConst_Cuda;
  v->ops->nvdotprod      = N_VDotProd_Cuda;
  v->ops->nvmaxnorm      = N_VMaxNorm_Cuda;
  v->ops->nvmin          = N_VMin_Cuda;
  v->ops->nvl1norm       = N_VL1Norm_Cuda;
  v->ops->nvinvtest      = N_VInvTest_Cuda;
  v->ops->nvconstrmask   = N_VConstrMask_Cuda;
  v->ops->nvminquotient  = N_VMinQuotient_Cuda;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Cuda;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Cuda;
  v->ops->nvwl2norm      = N_VWL2Norm_Cuda;
  v->ops->nvcompare      = N_VCompare_Cuda;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProd_Cuda;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Cuda;
  v->ops->nvminlocal         = N_VMin_Cuda;
  v->ops->nvl1normlocal      = N_VL1Norm_Cuda;
  v->ops->nvinvtestlocal     = N_VInvTest_Cuda;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Cuda;
  v->ops->nvminquotientlocal = N_VMinQuotient_Cuda;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Cuda;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Cuda;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = N_VDotProdMulti_Cuda;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_Cuda;
  v->ops->nvbufpack   = N_VBufPack_Cuda;
  v->ops->nvbufunpack = N_VBufUnpack_Cuda;

  /* print operation for debugging */
  v->ops->nvprint     = N_VPrint_Cuda;
  v->ops->nvprintfile = N_VPrintFile_Cuda;

  /* Create content */

  v->content = (N_VectorContent_Cuda) malloc(sizeof(_N_VectorContent_Cuda));
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return(NULL);
  }

  NVEC_CUDA_CONTENT(v)->priv = malloc(sizeof(_N_PrivateVectorContent_Cuda));
  if (NVEC_CUDA_CONTENT(v)->priv == NULL)
  {
    N_VDestroy(v);
    return(NULL);
  }

  // Initialize content
  NVEC_CUDA_CONTENT(v)->length             = 0;
  NVEC_CUDA_CONTENT(v)->host_data          = NULL;
  NVEC_CUDA_CONTENT(v)->device_data        = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = NULL;
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = NULL;
  NVEC_CUDA_CONTENT(v)->mem_helper         = NULL;
  NVEC_CUDA_CONTENT(v)->own_helper         = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->own_exec           = SUNTRUE;

  // Initialize private content
  NVEC_CUDA_PRIVATE(v)->use_managed_mem      = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev    = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host   = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_bytes  = 0;
  NVEC_CUDA_PRIVATE(v)->fused_buffer_dev     = NULL;
  NVEC_CUDA_PRIVATE(v)->fused_buffer_host    = NULL;
  NVEC_CUDA_PRIVATE(v)->fused_buffer_bytes   = 0;
  NVEC_CUDA_PRIVATE(v)->fused_buffer_offset  = 0;

  return(v);
}

N_Vector N_VNew_Cuda(sunindextype length, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length             = length;
  NVEC_CUDA_CONTENT(v)->mem_helper         = SUNMemoryHelper_Cuda(sunctx);
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem    = SUNFALSE;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VNewWithMemHelp_Cuda(sunindextype length, booleantype use_managed_mem, SUNMemoryHelper helper, SUNContext sunctx)
{
  N_Vector v;

  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Cuda: helper is NULL\n");
    return(NULL);
  }

  if (!SUNMemoryHelper_ImplementsRequiredOps(helper))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Cuda: helper doesn't implement all required ops\n");
    return(NULL);
  }

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length             = length;
  NVEC_CUDA_CONTENT(v)->mem_helper         = helper;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->own_helper         = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem    = use_managed_mem;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VNewManaged_Cuda(sunindextype length, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length             = length;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->mem_helper         = SUNMemoryHelper_Cuda(sunctx);
  NVEC_CUDA_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem    = SUNTRUE;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VMake_Cuda(sunindextype length, realtype *h_vdata, realtype *d_vdata, SUNContext sunctx)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length             = length;
  NVEC_CUDA_CONTENT(v)->host_data          = SUNMemoryHelper_Wrap(h_vdata, SUNMEMTYPE_HOST);
  NVEC_CUDA_CONTENT(v)->device_data        = SUNMemoryHelper_Wrap(d_vdata, SUNMEMTYPE_DEVICE);
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->mem_helper         = SUNMemoryHelper_Cuda(sunctx);
  NVEC_CUDA_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem    = SUNFALSE;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (NVEC_CUDA_CONTENT(v)->device_data == NULL ||
      NVEC_CUDA_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Cuda: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VMakeManaged_Cuda(sunindextype length, realtype *vdata, SUNContext sunctx)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda(sunctx);
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length             = length;
  NVEC_CUDA_CONTENT(v)->host_data          = SUNMemoryHelper_Wrap(vdata, SUNMEMTYPE_UVM);
  NVEC_CUDA_CONTENT(v)->device_data        = SUNMemoryHelper_Alias(NVEC_CUDA_CONTENT(v)->host_data);
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->mem_helper         = SUNMemoryHelper_Cuda(sunctx);
  NVEC_CUDA_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem    = SUNTRUE;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (NVEC_CUDA_CONTENT(v)->device_data == NULL ||
      NVEC_CUDA_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Cuda: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Set pointer to the raw host data. Does not free the existing pointer.
 */

void N_VSetHostArrayPointer_Cuda(realtype* h_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Cuda(v))
  {
    if (NVEC_CUDA_CONTENT(v)->host_data)
    {
      NVEC_CUDA_CONTENT(v)->host_data->ptr = (void*) h_vdata;
      NVEC_CUDA_CONTENT(v)->device_data->ptr = (void*) h_vdata;
    }
    else
    {
      NVEC_CUDA_CONTENT(v)->host_data = SUNMemoryHelper_Wrap((void*) h_vdata, SUNMEMTYPE_UVM);
      NVEC_CUDA_CONTENT(v)->device_data = SUNMemoryHelper_Alias(NVEC_CUDA_CONTENT(v)->host_data);
    }
  }
  else
  {
    if (NVEC_CUDA_CONTENT(v)->host_data)
    {
      NVEC_CUDA_CONTENT(v)->host_data->ptr = (void*) h_vdata;
    }
    else
    {
      NVEC_CUDA_CONTENT(v)->host_data = SUNMemoryHelper_Wrap((void*) h_vdata, SUNMEMTYPE_HOST);
    }
  }
}

/* ----------------------------------------------------------------------------
 * Set pointer to the raw device data
 */

void N_VSetDeviceArrayPointer_Cuda(realtype* d_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Cuda(v))
  {
    if (NVEC_CUDA_CONTENT(v)->device_data)
    {
      NVEC_CUDA_CONTENT(v)->device_data->ptr = (void*) d_vdata;
      NVEC_CUDA_CONTENT(v)->host_data->ptr = (void*) d_vdata;
    }
    else
    {
      NVEC_CUDA_CONTENT(v)->device_data = SUNMemoryHelper_Wrap((void*) d_vdata, SUNMEMTYPE_UVM);
      NVEC_CUDA_CONTENT(v)->host_data = SUNMemoryHelper_Alias(NVEC_CUDA_CONTENT(v)->device_data);
    }
  }
  else
  {
    if (NVEC_CUDA_CONTENT(v)->device_data)
    {
      NVEC_CUDA_CONTENT(v)->device_data->ptr = (void*) d_vdata;
    }
    else
    {
      NVEC_CUDA_CONTENT(v)->device_data = SUNMemoryHelper_Wrap((void*) d_vdata, SUNMEMTYPE_DEVICE);
    }
  }
}

/* ----------------------------------------------------------------------------
 * Return a flag indicating if the memory for the vector data is managed
 */
booleantype N_VIsManagedMemory_Cuda(N_Vector x)
{
  return NVEC_CUDA_PRIVATE(x)->use_managed_mem;
}

int N_VSetKernelExecPolicy_Cuda(N_Vector x,
                                SUNCudaExecPolicy* stream_exec_policy,
                                SUNCudaExecPolicy* reduce_exec_policy)
{
  if (x == NULL || stream_exec_policy == NULL || reduce_exec_policy == NULL)
    return(-1);

  if (NVEC_CUDA_CONTENT(x)->own_exec)
  {
    delete NVEC_CUDA_CONTENT(x)->stream_exec_policy;
    delete NVEC_CUDA_CONTENT(x)->reduce_exec_policy;
  }

  NVEC_CUDA_CONTENT(x)->stream_exec_policy = stream_exec_policy;
  NVEC_CUDA_CONTENT(x)->reduce_exec_policy = reduce_exec_policy;
  NVEC_CUDA_CONTENT(x)->own_exec = SUNFALSE;

  return(0);
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Cuda(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(x),
                                        NVEC_CUDA_CONTENT(x)->device_data,
                                        NVEC_CUDA_CONTENT(x)->host_data,
                                        NVEC_CUDA_MEMSIZE(x),
                                        (void*) NVEC_CUDA_STREAM(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyToDevice_Cuda: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  SUNDIALS_CUDA_VERIFY(cudaStreamSynchronize(*NVEC_CUDA_STREAM(x)));
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Cuda(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(x),
                                        NVEC_CUDA_CONTENT(x)->host_data,
                                        NVEC_CUDA_CONTENT(x)->device_data,
                                        NVEC_CUDA_MEMSIZE(x),
                                        (void*) NVEC_CUDA_STREAM(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyFromDevice_Cuda: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  SUNDIALS_CUDA_VERIFY(cudaStreamSynchronize(*NVEC_CUDA_STREAM(x)));
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to stdout
 */

void N_VPrint_Cuda(N_Vector x)
{
  N_VPrintFile_Cuda(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to outfile
 */

void N_VPrintFile_Cuda(N_Vector x, FILE *outfile)
{
  sunindextype i;

#ifdef SUNDIALS_DEBUG_PRINTVEC
  N_VCopyFromDevice_Cuda(x);
#endif

  for (i = 0; i < NVEC_CUDA_CONTENT(x)->length; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Le\n", NVEC_CUDA_HDATAp(x)[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16e\n", NVEC_CUDA_HDATAp(x)[i]);
#else
    fprintf(outfile, "%11.8e\n", NVEC_CUDA_HDATAp(x)[i]);
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

N_Vector N_VCloneEmpty_Cuda(N_Vector w)
{
  N_Vector v;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty_Cuda(w->sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Set content */
  NVEC_CUDA_CONTENT(v)->length          = NVEC_CUDA_CONTENT(w)->length;
  NVEC_CUDA_CONTENT(v)->own_exec        = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem = NVEC_CUDA_PRIVATE(w)->use_managed_mem;

  return(v);
}

N_Vector N_VClone_Cuda(N_Vector w)
{
  N_Vector v;

  v = NULL;
  v = N_VCloneEmpty_Cuda(w);
  if (v == NULL) return(NULL);

  NVEC_CUDA_MEMHELP(v) = SUNMemoryHelper_Clone(NVEC_CUDA_MEMHELP(w));
  NVEC_CUDA_CONTENT(v)->own_helper = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy = NVEC_CUDA_CONTENT(w)->stream_exec_policy->clone();
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy = NVEC_CUDA_CONTENT(w)->reduce_exec_policy->clone();

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Cuda: SUNMemoryHelper_Clone returned NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

void N_VDestroy_Cuda(N_Vector v)
{
  N_VectorContent_Cuda vc;
  N_PrivateVectorContent_Cuda vcp;

  if (v == NULL) return;

  /* free ops structure */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }

  /* extract content */
  vc = NVEC_CUDA_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

  /* free private content */
  vcp = (N_PrivateVectorContent_Cuda) vc->priv;
  if (vcp != NULL)
  {
    /* free items in private content */
    FreeReductionBuffer(v);
    FusedBuffer_Free(v);
    free(vcp);
    vc->priv = NULL;
  }

  /* free items in content */
  if (vc->own_exec)
  {
    delete vc->stream_exec_policy;
    vc->stream_exec_policy = NULL;
    delete vc->reduce_exec_policy;
    vc->reduce_exec_policy = NULL;
  }

  if (NVEC_CUDA_MEMHELP(v))
  {
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vc->host_data, nullptr);
    vc->host_data = NULL;
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vc->device_data, nullptr);
    vc->device_data = NULL;
    if (vc->own_helper) SUNMemoryHelper_Destroy(vc->mem_helper);
    vc->mem_helper = NULL;
  }

  /* free content struct */
  free(vc);

  /* free vector */
  free(v);

  return;
}

void N_VSpace_Cuda(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NVEC_CUDA_CONTENT(X)->length;
  *liw = 2;
}

void N_VConst_Cuda(realtype a, N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConst_Cuda: GetKernelParameters returned nonzero\n");
  }

  setConstKernel<<<grid, block, shMemSize, stream>>>
  (
    a,
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VLinearSum_Cuda(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSum_Cuda: GetKernelParameters returned nonzero\n");
  }

  linearSumKernel<<<grid, block, shMemSize, stream>>>
  (
    a,
    NVEC_CUDA_DDATAp(X),
    b,
    NVEC_CUDA_DDATAp(Y),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VProd_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VProd_Cuda: GetKernelParameters returned nonzero\n");
  }


  prodKernel<<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Y),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VDiv_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDiv_Cuda: GetKernelParameters returned nonzero\n");
  }

  divKernel<<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Y),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VScale_Cuda(realtype a, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScale_Cuda: GetKernelParameters returned nonzero\n");
  }

  scaleKernel<<<grid, block, shMemSize, stream>>>
  (
    a,
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VAbs_Cuda(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VAbs_Cuda: GetKernelParameters returned nonzero\n");
  }

  absKernel<<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VInv_Cuda(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInv_Cuda: GetKernelParameters returned nonzero\n");
  }

  invKernel<<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VAddConst_Cuda(N_Vector X, realtype b, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VAddConst_Cuda: GetKernelParameters returned nonzero\n");
  }

  addConstKernel<<<grid, block, shMemSize, stream>>>
  (
    b,
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

realtype N_VDotProd_Cuda(N_Vector X, N_Vector Y)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = ZERO;

  if (InitializeReductionBuffer(X, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProd_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProd_Cuda: GetKernelParameters returned nonzero\n");
  }

  dotProdKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Y),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VMaxNorm_Cuda(N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = ZERO;

  if (InitializeReductionBuffer(X, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMaxNorm_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMaxNorm_Cuda: GetKernelParameters returned nonzero\n");
  }

  maxNormKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWSqrSumLocal_Cuda(N_Vector X, N_Vector W)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = ZERO;

  if (InitializeReductionBuffer(X, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumLocal_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumLocal_Cuda: GetKernelParameters returned nonzero\n");
  }

  wL2NormSquareKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(W),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWrmsNorm_Cuda(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Cuda(X, W);
  return std::sqrt(sum/NVEC_CUDA_CONTENT(X)->length);
}

realtype N_VWSqrSumMaskLocal_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = ZERO;

  if (InitializeReductionBuffer(X, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumMaskLocal_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumMaskLocal_Cuda: GetKernelParameters returned nonzero\n");
  }

  wL2NormSquareMaskKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(W),
    NVEC_CUDA_DDATAp(Id),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWrmsNormMask_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  const realtype sum = N_VWSqrSumMaskLocal_Cuda(X, W, Id);
  return std::sqrt(sum/NVEC_CUDA_CONTENT(X)->length);
}

realtype N_VMin_Cuda(N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = std::numeric_limits<realtype>::max();

  if (InitializeReductionBuffer(X, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMin_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMin_Cuda: GetKernelParameters returned nonzero\n");
  }

  findMinKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    gpu_result,
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWL2Norm_Cuda(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Cuda(X, W);
  return std::sqrt(sum);
}

realtype N_VL1Norm_Cuda(N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = ZERO;

  if (InitializeReductionBuffer(X, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VL1Norm_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VL1Norm_Cuda: GetKernelParameters returned nonzero\n");
  }

  L1NormKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

void N_VCompare_Cuda(realtype c, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCompare_Cuda: GetKernelParameters returned nonzero\n");
  }

  compareKernel<<<grid, block, shMemSize, stream>>>
  (
    c,
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();
}

booleantype N_VInvTest_Cuda(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = ZERO;

  if (InitializeReductionBuffer(X, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInvTest_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInvTest_Cuda: GetKernelParameters returned nonzero\n");
  }

  invTestKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return (gpu_result < HALF);
}

booleantype N_VConstrMask_Cuda(N_Vector C, N_Vector X, N_Vector M)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = ZERO;

  if (InitializeReductionBuffer(X, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstrMask_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(X, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstrMask_Cuda: GetKernelParameters returned nonzero\n");
  }

  constrMaskKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(C),
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DDATAp(M),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return (gpu_result < HALF);
}

realtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom)
{
  // Starting value for min reduction
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  realtype gpu_result = std::numeric_limits<realtype>::max();;

  if (InitializeReductionBuffer(num, &gpu_result))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMinQuotient_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(num, true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMinQuotient_Cuda: GetKernelParameters returned nonzero\n");
  }

  minQuotientKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    gpu_result,
    NVEC_CUDA_DDATAp(num),
    NVEC_CUDA_DDATAp(denom),
    NVEC_CUDA_DBUFFERp(num),
    NVEC_CUDA_CONTENT(num)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(num);
  gpu_result = NVEC_CUDA_HBUFFERp(num)[0];

  return gpu_result;
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */


int N_VLinearCombination_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  // Fused op workspace shortcuts
  realtype*  cdata = NULL;
  realtype** xdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(z, nvec, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(z, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Cuda: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(z, X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(z))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Set kernel parameters and launch
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }

  linearCombinationKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    cdata,
    xdata,
    NVEC_CUDA_DDATAp(z),
    NVEC_CUDA_CONTENT(z)->length
  );
  PostKernelLaunch();

  return 0;
}


int N_VScaleAddMulti_Cuda(int nvec, realtype* c, N_Vector x, N_Vector* Y,
                          N_Vector* Z)
{
  // Shortcuts to the fused op workspace
  realtype*  cdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(x, nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(x, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Cuda: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(x, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }

  scaleAddMultiKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    cdata,
    NVEC_CUDA_DDATAp(x),
    ydata,
    zdata,
    NVEC_CUDA_CONTENT(x)->length
  );
  PostKernelLaunch();

  return 0;
}


int N_VDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y, realtype* dots)
{
  // Fused op workspace shortcuts
  realtype** ydata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(x, 0, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProdMulti_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProdMulti_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProdMulti_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Setup the reduction buffer
  for (int i = 0; i < nvec; ++i)
  {
    dots[i] = ZERO;
  }

  if (InitializeReductionBuffer(x, dots, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProd_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(x, false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProdMulti_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }
  grid = nvec;

  dotProdMultiKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    nvec,
    NVEC_CUDA_DDATAp(x),
    ydata,
    NVEC_CUDA_DBUFFERp(x),
    NVEC_CUDA_CONTENT(x)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(x, nvec);
  for (int i = 0; i < nvec; ++i)
  {
    dots[i] = NVEC_CUDA_HBUFFERp(x)[i];
  }

  return 0;
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */


int N_VLinearSumVectorArray_Cuda(int nvec,
                                 realtype a, N_Vector* X,
                                 realtype b, N_Vector* Y,
                                 N_Vector* Z)
{
  // Shortcuts to the fused op workspace
  realtype** xdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], 0, 3 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinaerSumVectorArray_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }

  linearSumVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    a,
    xdata,
    b,
    ydata,
    zdata,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  return 0;
}


int N_VScaleVectorArray_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  // Shortcuts to the fused op workspace arrays
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(Z[0], c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Cuda: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }

  scaleVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    cdata,
    xdata,
    zdata,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  return 0;
}


int N_VConstVectorArray_Cuda(int nvec, realtype c, N_Vector* Z)
{
  // Shortcuts to the fused op workspace arrays
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], 0, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstVectorArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstVectorArray_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstVectorArray_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }

  constVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    c,
    zdata,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  return 0;
}


int N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                realtype* norms)
{
  // Fused op workspace shortcuts
  realtype** xdata = NULL;
  realtype** wdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(W[0], 0, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(W[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(W[0], W, nvec, &wdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(W[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Setup the reduction buffer
  for (int i = 0; i < nvec; ++i)
  {
    norms[i] = ZERO;
  }

  if (InitializeReductionBuffer(W[0], norms, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(W[0], true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }
  grid = nvec;

  wL2NormSquareVectorArrayKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    nvec,
    xdata,
    wdata,
    NVEC_CUDA_DBUFFERp(W[0]),
    NVEC_CUDA_CONTENT(W[0])->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(W[0], nvec);
  for (int i = 0; i < nvec; ++i)
  {
    norms[i] = std::sqrt(NVEC_CUDA_HBUFFERp(W[0])[i] /
                         NVEC_CUDA_CONTENT(W[0])->length);
  }

  return 0;
}


int N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                    N_Vector id, realtype* norms)
{
  // Fused op workspace shortcuts
  realtype** xdata = NULL;
  realtype** wdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(W[0], 0, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(W[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(W[0], W, nvec, &wdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(W[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Setup the reduction buffer
  for (int i = 0; i < nvec; ++i)
  {
    norms[i] = ZERO;
  }

  if (InitializeReductionBuffer(W[0], norms, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormVectorArray_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(W[0], true, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWrmsNormMaskVectorArray_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }
  grid = nvec;

  wL2NormSquareMaskVectorArrayKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    nvec,
    xdata,
    wdata,
    NVEC_CUDA_DDATAp(id),
    NVEC_CUDA_DBUFFERp(W[0]),
    NVEC_CUDA_CONTENT(W[0])->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(W[0], nvec);
  for (int i = 0; i < nvec; ++i)
  {
    norms[i] = std::sqrt(NVEC_CUDA_HBUFFERp(W[0])[i] /
                         NVEC_CUDA_CONTENT(W[0])->length);
  }

  return 0;
}


int N_VScaleAddMultiVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  // Shortcuts to the fused op workspace
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(X[0], nsum, nvec + 2 * nvec * nsum))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(X[0], c, nsum, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiArray_Cuda: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(X[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray2D(X[0], Y, nvec, nsum, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiVectorArray_Cuda: FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray2D(X[0], Z, nvec, nsum, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiVectorArray_Cuda: FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(X[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiVectorArray_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }

  scaleAddMultiVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    nsum,
    cdata,
    xdata,
    ydata,
    zdata,
    NVEC_CUDA_CONTENT(X[0])->length
  );
  PostKernelLaunch();

  return 0;
}


int N_VLinearCombinationVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  // Shortcuts to the fused op workspace arrays
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], nsum, nvec + nvec * nsum))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Cuda: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(Z[0], c, nsum, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Cuda: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray2D(Z[0], X, nvec, nsum, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Cuda: FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Cuda: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Cuda: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Cuda: GetKernelParameters returned nonzero\n");
    return -1;
  }

  linearCombinationVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    nsum,
    cdata,
    xdata,
    zdata,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  return 0;
}


/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */


int N_VBufSize_Cuda(N_Vector x, sunindextype *size)
{
  if (x == NULL) return(-1);
  *size = (sunindextype)NVEC_CUDA_MEMSIZE(x);
  return(0);
}


int N_VBufPack_Cuda(N_Vector x, void *buf)
{
  int copy_fail = 0;
  cudaError_t cuerr;

  if (x == NULL || buf == NULL) return(-1);

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return(-1);

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(x),
                                        buf_mem,
                                        NVEC_CUDA_CONTENT(x)->device_data,
                                        NVEC_CUDA_MEMSIZE(x),
                                        (void*) NVEC_CUDA_STREAM(x));

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = cudaStreamSynchronize(*NVEC_CUDA_STREAM(x));

  SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(x), buf_mem, nullptr);

  return (!SUNDIALS_CUDA_VERIFY(cuerr) || copy_fail ? -1 : 0);
}


int N_VBufUnpack_Cuda(N_Vector x, void *buf)
{
  int copy_fail = 0;
  cudaError_t cuerr;

  if (x == NULL || buf == NULL) return(-1);

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return(-1);

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(x),
                                        NVEC_CUDA_CONTENT(x)->device_data,
                                        buf_mem,
                                        NVEC_CUDA_MEMSIZE(x),
                                        (void*) NVEC_CUDA_STREAM(x));

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = cudaStreamSynchronize(*NVEC_CUDA_STREAM(x));

  SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(x), buf_mem, nullptr);

  return (!SUNDIALS_CUDA_VERIFY(cuerr) || copy_fail ? -1 : 0);
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */


int N_VEnableFusedOps_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Cuda;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Cuda;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Cuda;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Cuda;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Cuda;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Cuda;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Cuda;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Cuda;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Cuda;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Cuda;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = N_VDotProdMulti_Cuda;
  }
  else
  {
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

int N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvlinearcombination = tf ? N_VLinearCombination_Cuda : NULL;
  return 0;
}


int N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvscaleaddmulti = tf ? N_VScaleAddMulti_Cuda : NULL;
  return 0;
}


int N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvdotprodmulti      = tf ? N_VDotProdMulti_Cuda : NULL;
  v->ops->nvdotprodmultilocal = tf ? N_VDotProdMulti_Cuda : NULL;
  return 0;
}


int N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvlinearsumvectorarray = tf ? N_VLinearSumVectorArray_Cuda : NULL;
  return 0;
}


int N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvscalevectorarray = tf ? N_VScaleVectorArray_Cuda : NULL;
  return 0;
}


int N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvconstvectorarray = tf ? N_VConstVectorArray_Cuda : NULL;
  return 0;
}


int N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvwrmsnormvectorarray = tf ? N_VWrmsNormVectorArray_Cuda : NULL;
  return 0;
}


int N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvwrmsnormmaskvectorarray = tf ?
    N_VWrmsNormMaskVectorArray_Cuda : NULL;
  return 0;
}


int N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvscaleaddmultivectorarray = tf ?
    N_VScaleAddMultiVectorArray_Cuda : NULL;
  return 0;
}


int N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvlinearcombinationvectorarray = tf ?
    N_VLinearCombinationVectorArray_Cuda : NULL;
  return 0;
}


/*
 * Private helper functions.
 */


int AllocateData(N_Vector v)
{
  int alloc_fail = 0;
  N_VectorContent_Cuda vc = NVEC_CUDA_CONTENT(v);
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  if (N_VGetLength_Cuda(v) == 0) return(0);

  if (vcp->use_managed_mem)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v), &(vc->device_data),
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_UVM,
                                       nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed for SUNMEMTYPE_UVM\n");
    }
    vc->host_data = SUNMemoryHelper_Alias(vc->device_data);
  }
  else
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v), &(vc->host_data),
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_HOST,
                                       nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
    }

    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v), &(vc->device_data),
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_DEVICE,
                                       nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  return(alloc_fail ? -1 : 0);
}

/*
 * Initializes the internal buffer used for reductions.
 * If the buffer is already allocated, it will only be reallocated
 * if it is no longer large enough. This may occur if the length
 * of the vector is increased. The buffer is initialized to the
 * value given.
 */
int InitializeReductionBuffer(N_Vector v, const realtype* value, size_t n)
{
  int         alloc_fail = 0;
  int         copy_fail  = 0;
  booleantype alloc_mem  = SUNFALSE;
  size_t      bytes      = n * sizeof(realtype);

  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Wrap the initial value as SUNMemory object
  SUNMemory value_mem = SUNMemoryHelper_Wrap((void*) value, SUNMEMTYPE_HOST);

  // Check if the existing reduction memory is not large enough
  if (vcp->reduce_buffer_bytes < bytes)
  {
    FreeReductionBuffer(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
    // Allocate pinned memory on the host
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->reduce_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED, nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("WARNING in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      // If pinned alloc failed, allocate plain host memory
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                         &(vcp->reduce_buffer_host), bytes,
                                         SUNMEMTYPE_HOST, nullptr);
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
      }
    }

    // Allocate device memory
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->reduce_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE, nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  if (!alloc_fail)
  {
    // Store the size of the reduction memory buffer
    vcp->reduce_buffer_bytes = bytes;

    // Initialize the memory with the value
    copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(v),
                                          vcp->reduce_buffer_dev, value_mem,
                                          bytes, (void*) NVEC_CUDA_STREAM(v));

    if (copy_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_CopyAsync failed\n");
    }
  }

  // Deallocate the wrapper
  SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), value_mem, nullptr);

  return((alloc_fail || copy_fail) ? -1 : 0);
}

/* Free the reduction buffer
 */
void FreeReductionBuffer(N_Vector v)
{
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  if (vcp == NULL) return;

  // Free device mem
  if (vcp->reduce_buffer_dev != NULL)
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vcp->reduce_buffer_dev,
                            nullptr);
  vcp->reduce_buffer_dev  = NULL;

  // Free host mem
  if (vcp->reduce_buffer_host != NULL)
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vcp->reduce_buffer_host,
                            nullptr);
  vcp->reduce_buffer_host = NULL;

  // Reset allocated memory size
  vcp->reduce_buffer_bytes = 0;
}

/* Copy the reduction buffer from the device to the host.
 */
int CopyReductionBufferFromDevice(N_Vector v, size_t n)
{
  int copy_fail;
  cudaError_t cuerr;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(v),
                                        NVEC_CUDA_PRIVATE(v)->reduce_buffer_host,
                                        NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev,
                                        n * sizeof(realtype),
                                        (void*) NVEC_CUDA_STREAM(v));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in CopyReductionBufferFromDevice: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = cudaStreamSynchronize(*NVEC_CUDA_STREAM(v));
  return (!SUNDIALS_CUDA_VERIFY(cuerr) || copy_fail ? -1 : 0);
}


static int FusedBuffer_Init(N_Vector v, int nreal, int nptr)
{
  int         alloc_fail = 0;
  booleantype alloc_mem  = SUNFALSE;

  // pad buffer with single precision data
#if defined(SUNDIALS_SINGLE_PRECISION)
  size_t bytes = nreal * 2 * sizeof(realtype) + nptr * sizeof(realtype*);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  size_t bytes = nreal * sizeof(realtype) + nptr * sizeof(realtype*);
#else
#error Incompatible precision for CUDA
#endif

  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check if the existing memory is not large enough
  if (vcp->fused_buffer_bytes < bytes)
  {
    FusedBuffer_Free(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
    // Allocate pinned memory on the host
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->fused_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED, nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("WARNING in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      // If pinned alloc failed, allocate plain host memory
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                         &(vcp->fused_buffer_host), bytes,
                                         SUNMEMTYPE_HOST, nullptr);
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
        return -1;
      }
    }

    // Allocate device memory
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->fused_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE, nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
      return -1;
    }

    // Store the size of the fused op buffer
    vcp->fused_buffer_bytes = bytes;
  }

  // Reset the buffer offset
  vcp->fused_buffer_offset = 0;

  return 0;
}


static int FusedBuffer_CopyRealArray(N_Vector v, realtype *rdata, int nval,
                                     realtype **shortcut)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check buffer space and fill the host buffer
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyRealArray: Buffer offset is exceedes the buffer size\n");
    return -1;
  }

  realtype* h_buffer = (realtype*) ((char*)(vcp->fused_buffer_host->ptr) +
                                    vcp->fused_buffer_offset);

  for (int j = 0; j < nval; j++)
  {
    h_buffer[j] = rdata[j];
  }

  // Set shortcut to the device buffer and update offset
  *shortcut = (realtype*) ((char*)(vcp->fused_buffer_dev->ptr) +
                           vcp->fused_buffer_offset);

  // accounting for buffer padding
#if defined(SUNDIALS_SINGLE_PRECISION)
  vcp->fused_buffer_offset += nval * 2 * sizeof(realtype);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  vcp->fused_buffer_offset += nval * sizeof(realtype);
#else
#error Incompatible precision for CUDA
#endif

  return 0;
}


static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector *X, int nvec,
                                      realtype ***shortcut)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check buffer space and fill the host buffer
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyPtrArray1D: Buffer offset is exceedes the buffer size\n");
    return -1;
  }

  realtype** h_buffer = (realtype**) ((char*)(vcp->fused_buffer_host->ptr) +
                                      vcp->fused_buffer_offset);

  for (int j = 0; j < nvec; j++)
  {
    h_buffer[j] = NVEC_CUDA_DDATAp(X[j]);
  }

  // Set shortcut to the device buffer and update offset
  *shortcut = (realtype**) ((char*)(vcp->fused_buffer_dev->ptr) +
                            vcp->fused_buffer_offset);

  vcp->fused_buffer_offset += nvec * sizeof(realtype*);

  return 0;
}


static int FusedBuffer_CopyPtrArray2D(N_Vector v, N_Vector **X, int nvec,
                                      int nsum, realtype ***shortcut)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Check buffer space and fill the host buffer
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyPtrArray2D: Buffer offset is exceedes the buffer size\n");
    return -1;
  }

  realtype** h_buffer = (realtype**) ((char*)(vcp->fused_buffer_host->ptr) +
                                      vcp->fused_buffer_offset);

  for (int j = 0; j < nvec; j++)
  {
    for (int k = 0; k < nsum; k++)
    {
      h_buffer[j * nsum + k] = NVEC_CUDA_DDATAp(X[k][j]);
    }
  }

  // Set shortcut to the device buffer and update offset
  *shortcut = (realtype**) ((char*)(vcp->fused_buffer_dev->ptr) +
                            vcp->fused_buffer_offset);

  // Update the offset
  vcp->fused_buffer_offset += nvec * nsum * sizeof(realtype*);

  return 0;
}


static int FusedBuffer_CopyToDevice(N_Vector v)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  // Copy the fused buffer to the device
  int copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(v),
                                            vcp->fused_buffer_dev,
                                            vcp->fused_buffer_host,
                                            vcp->fused_buffer_offset,
                                            (void*) NVEC_CUDA_STREAM(v));
  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyToDevice: SUNMemoryHelper_CopyAsync failed\n");
    return -1;
  }

  // Synchronize with respect to the host, but only in this stream
  SUNDIALS_CUDA_VERIFY(cudaStreamSynchronize(*NVEC_CUDA_STREAM(v)));

  return 0;
}


static int FusedBuffer_Free(N_Vector v)
{
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  if (vcp == NULL) return 0;

  if (vcp->fused_buffer_host)
  {
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v),
                            vcp->fused_buffer_host, nullptr);
    vcp->fused_buffer_host = NULL;
  }

  if (vcp->fused_buffer_dev)
  {
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v),
                            vcp->fused_buffer_dev, nullptr);
    vcp->fused_buffer_dev = NULL;
  }

  vcp->fused_buffer_bytes  = 0;
  vcp->fused_buffer_offset = 0;

  return 0;
}


/* Get the kernel launch parameters based on the kernel type (reduction or not),
 * using the appropriate kernel execution policy.
 */
static int GetKernelParameters(N_Vector v, booleantype reduction, size_t& grid,
                               size_t& block, size_t& shMemSize,
                               cudaStream_t& stream, size_t n)
{
  n = (n == 0) ? NVEC_CUDA_CONTENT(v)->length : n;
  if (reduction)
  {
    SUNCudaExecPolicy* reduce_exec_policy = NVEC_CUDA_CONTENT(v)->reduce_exec_policy;
    grid      = reduce_exec_policy->gridSize(n);
    block     = reduce_exec_policy->blockSize();
    shMemSize = 0;
    stream    = *(reduce_exec_policy->stream());
    if (block % CUDA_WARP_SIZE)
    {
#ifdef SUNDIALS_DEBUG
      throw std::runtime_error("the block size must be a multiple must be of CUDA warp size");
#endif
      return(-1);
    }
  }
  else
  {
    SUNCudaExecPolicy* stream_exec_policy = NVEC_CUDA_CONTENT(v)->stream_exec_policy;
    grid      = stream_exec_policy->gridSize(n);
    block     = stream_exec_policy->blockSize();
    shMemSize = 0;
    stream    = *(stream_exec_policy->stream());
  }

  if (grid == 0)
  {
#ifdef SUNDIALS_DEBUG
    throw std::runtime_error("the grid size must be > 0");
#endif
    return(-1);
  }
  if (block == 0)
  {
#ifdef SUNDIALS_DEBUG
    throw std::runtime_error("the block size must be > 0");
#endif
    return(-1);
  }

  return(0);
}

/* Should be called after a kernel launch.
 * If SUNDIALS_DEBUG_CUDA_LASTERROR is not defined, then the function does nothing.
 * If it is defined, the function will synchronize and check the last CUDA error.
 */
void PostKernelLaunch()
{
#ifdef SUNDIALS_DEBUG_CUDA_LASTERROR
  cudaDeviceSynchronize();
  SUNDIALS_CUDA_VERIFY(cudaGetLastError());
#endif
}


} // extern "C"
