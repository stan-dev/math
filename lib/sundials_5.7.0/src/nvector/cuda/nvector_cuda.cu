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

#include <nvector/nvector_cuda.h>
#include "VectorArrayKernels.cuh"
#include "VectorKernels.cuh"

#include "sundials_cuda.h"
#include "sundials_debug.h"

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)

extern "C" {

using namespace sundials;
using namespace sundials::nvector_cuda;

/*
 * Macro definitions
 */

#define NVEC_CUDA_CONTENT(x)  ((N_VectorContent_Cuda)(x->content))
#define NVEC_CUDA_PRIVATE(x)  ((N_PrivateVectorContent_Cuda)(NVEC_CUDA_CONTENT(x)->priv))
#define NVEC_CUDA_MEMSIZE(x)  (NVEC_CUDA_CONTENT(x)->length * sizeof(realtype))
#define NVEC_CUDA_MEMHELP(x)  (NVEC_CUDA_CONTENT(x)->mem_helper)
#define NVEC_CUDA_HDATAp(x)   ((realtype*) NVEC_CUDA_CONTENT(x)->host_data->ptr)
#define NVEC_CUDA_DDATAp(x)   ((realtype*) NVEC_CUDA_CONTENT(x)->device_data->ptr)
#define NVEC_CUDA_HBUFFERp(x) ((realtype*) NVEC_CUDA_PRIVATE(x)->reduce_buffer_host->ptr)
#define NVEC_CUDA_DBUFFERp(x) ((realtype*) NVEC_CUDA_PRIVATE(x)->reduce_buffer_dev->ptr)
#define NVEC_CUDA_STREAM(x)   (NVEC_CUDA_CONTENT(x)->stream_exec_policy->stream())


/*
 * Private structure definition
 */

struct _N_PrivateVectorContent_Cuda
{
  booleantype     use_managed_mem;               /* indicates if the data pointers and buffer pointers are managed memory */
  size_t          reduce_buffer_allocated_bytes; /* current size of the reduction buffer */
  SUNMemory       reduce_buffer_dev;             /* device buffer used for reductions */
  SUNMemory       reduce_buffer_host;            /* host buffer used for reductions */
};

typedef struct _N_PrivateVectorContent_Cuda *N_PrivateVectorContent_Cuda;

/*
 * Private function definitions
 */

static int AllocateData(N_Vector v);
static int InitializeReductionBuffer(N_Vector v, const realtype value);
static void FreeReductionBuffer(N_Vector v);
static int CopyReductionBufferFromDevice(N_Vector v, size_t n = 1);
static int GetKernelParameters(N_Vector v, booleantype reduction, size_t& grid, size_t& block,
                               size_t& shMemSize, cudaStream_t& stream, size_t n = 0);
static void PostKernelLaunch();

/*
 * Private functions needed for N_VMakeWithManagedAllocator_Cuda
 * backwards compatibility.
 */

/* DEPRECATION NOTICE: The 4 functions below can be removed once
   N_VMakeWithManagedAllocator_Cuda (deprecated) is removed in the
   next major release. The UserAllocHelper struct can also be removed. */

/* Struct that we use to pack up the user
   provided alloc and free functions. */
typedef struct _UserAllocHelper
{
  void*  (*userallocfn)(size_t);
  void   (*userfreefn)(void*);
} UserAllocHelper;

static int UserAlloc(SUNMemoryHelper helper, SUNMemory* memptr,
                     size_t memsize, SUNMemoryType mem_type)
{
  UserAllocHelper* ua = (UserAllocHelper*) helper->content;
  SUNMemory mem = SUNMemoryNewEmpty();

  mem->type = SUNMEMTYPE_UVM;
  mem->ptr  = ua->userallocfn(memsize);
  mem->own  = SUNTRUE;
  if (mem->ptr == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in UserAlloc: user provided alloc failed\n");
    free(mem);
    return(-1);
  }

  *memptr = mem;
  return(0);
}

static int UserDealloc(SUNMemoryHelper helper, SUNMemory mem)
{
  UserAllocHelper* ua = (UserAllocHelper*) helper->content;
  if (mem->own)
  {
    ua->userfreefn(mem->ptr);
    mem->ptr = NULL;
  }
  free(mem);
  return(0);
}

static SUNMemoryHelper HelperClone(SUNMemoryHelper helper)
{
  UserAllocHelper* uaclone;
  UserAllocHelper* ua = (UserAllocHelper*) helper->content;
  SUNMemoryHelper hclone = SUNMemoryHelper_NewEmpty();

  SUNMemoryHelper_CopyOps(helper, hclone);

  uaclone = (UserAllocHelper*) malloc(sizeof(UserAllocHelper));
  uaclone->userallocfn = ua->userallocfn;
  uaclone->userfreefn  = ua->userfreefn;

  hclone->content = uaclone;

  return(hclone);
}

static int HelperDestroy(SUNMemoryHelper helper)
{
  free(helper->content);
  helper->content = NULL;
  free(helper->ops);
  free(helper);
  return(0);
}

N_Vector N_VNewEmpty_Cuda()
{
  N_Vector v;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty();
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

  NVEC_CUDA_CONTENT(v)->length                        = 0;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = NULL;
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = NULL;
  NVEC_CUDA_CONTENT(v)->mem_helper                    = NULL;
  NVEC_CUDA_CONTENT(v)->own_helper                    = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  return(v);
}

N_Vector N_VNew_Cuda(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Cuda();
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaBlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

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

N_Vector N_VNewWithMemHelp_Cuda(sunindextype length, booleantype use_managed_mem, SUNMemoryHelper helper)
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
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->mem_helper                    = helper;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaBlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->own_helper                    = SUNFALSE;
  NVEC_CUDA_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = use_managed_mem;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Cuda: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VNewManaged_Cuda(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaBlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Cuda();
  NVEC_CUDA_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

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

N_Vector N_VMake_Cuda(sunindextype length, realtype *h_vdata, realtype *d_vdata)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->host_data                     = SUNMemoryHelper_Wrap(h_vdata, SUNMEMTYPE_HOST);
  NVEC_CUDA_CONTENT(v)->device_data                   = SUNMemoryHelper_Wrap(d_vdata, SUNMEMTYPE_DEVICE);
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaBlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Cuda();
  NVEC_CUDA_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

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

N_Vector N_VMakeManaged_Cuda(sunindextype length, realtype *vdata)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->host_data                     = SUNMemoryHelper_Wrap(vdata, SUNMEMTYPE_UVM);
  NVEC_CUDA_CONTENT(v)->device_data                   = SUNMemoryHelper_Alias(NVEC_CUDA_CONTENT(v)->host_data);
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaBlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Cuda();
  NVEC_CUDA_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

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

N_Vector N_VMakeWithManagedAllocator_Cuda(sunindextype length,
                                          void* (*allocfn)(size_t),
                                          void (*freefn)(void*))
{
  UserAllocHelper* ua;
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  NVEC_CUDA_CONTENT(v)->length                        = length;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->stream_exec_policy            = new CudaThreadDirectExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->reduce_exec_policy            = new CudaBlockReduceExecPolicy(256);
  NVEC_CUDA_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Cuda();
  NVEC_CUDA_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_CUDA_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  if (NVEC_CUDA_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeWithManagedAllocator_Cuda: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  ua = (UserAllocHelper*) malloc(sizeof(UserAllocHelper));
  ua->userallocfn                    = allocfn;
  ua->userfreefn                     = freefn;
  NVEC_CUDA_MEMHELP(v)->content      = (void*) ua;
  NVEC_CUDA_MEMHELP(v)->ops->alloc   = UserAlloc;
  NVEC_CUDA_MEMHELP(v)->ops->dealloc = UserDealloc;
  NVEC_CUDA_MEMHELP(v)->ops->clone   = HelperClone;
  NVEC_CUDA_MEMHELP(v)->ops->destroy = HelperDestroy;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeWithManagedAllocator_Cuda: AllocateData returned nonzero\n");
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

/*
 * ----------------------------------------------------------------------------
 * DEPRECATED: will be removed in SUNDIALS v6.
 * Sets the cudaStream_t to use for execution of the CUDA kernels.
 */
void N_VSetCudaStream_Cuda(N_Vector x, cudaStream_t *stream)
{
  const CudaExecPolicy* xs = NVEC_CUDA_CONTENT(x)->stream_exec_policy;
  const CudaExecPolicy* xr = NVEC_CUDA_CONTENT(x)->reduce_exec_policy;
  CudaThreadDirectExecPolicy* s =
    new CudaThreadDirectExecPolicy(xs->blockSize(), *stream);
  CudaBlockReduceExecPolicy* r =
    new CudaBlockReduceExecPolicy(xr->blockSize(), xr->gridSize(), *stream);
  N_VSetKernelExecPolicy_Cuda(x, s, r);
  NVEC_CUDA_CONTENT(x)->own_exec = SUNTRUE;
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

  for (i = 0; i < NVEC_CUDA_CONTENT(x)->length; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", NVEC_CUDA_HDATAp(x)[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", NVEC_CUDA_HDATAp(x)[i]);
#else
    fprintf(outfile, "%11.8g\n", NVEC_CUDA_HDATAp(x)[i]);
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
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Set content */
  NVEC_CUDA_CONTENT(v)->length                        = NVEC_CUDA_CONTENT(w)->length;
  NVEC_CUDA_CONTENT(v)->host_data                     = NULL;
  NVEC_CUDA_CONTENT(v)->device_data                   = NULL;
  NVEC_CUDA_CONTENT(v)->mem_helper                    = NULL;
  NVEC_CUDA_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_CUDA_PRIVATE(v)->use_managed_mem               = NVEC_CUDA_PRIVATE(w)->use_managed_mem;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_CUDA_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

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
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vc->host_data);
    vc->host_data = NULL;
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vc->device_data);
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

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProd_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
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
  realtype gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VMaxNorm_Cuda(N_Vector X)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMaxNorm_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  maxNormKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWSqrSumLocal_Cuda(N_Vector X, N_Vector W)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumLocal_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
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
  realtype gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

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

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumMaskLocal_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
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
  realtype gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWrmsNormMask_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  const realtype sum = N_VWSqrSumMaskLocal_Cuda(X, W, Id);
  return std::sqrt(sum/NVEC_CUDA_CONTENT(X)->length);
}

realtype N_VMin_Cuda(N_Vector X)
{
  const realtype maxVal = std::numeric_limits<realtype>::max();

  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (InitializeReductionBuffer(X, maxVal))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMin_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  findMinKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    maxVal,
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

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

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VL1Norm_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  L1NormKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    NVEC_CUDA_DDATAp(X),
    NVEC_CUDA_DBUFFERp(X),
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return gpu_result;
}

void N_VCompare_Cuda(realtype c, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
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

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInvTest_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
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
  realtype gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return (gpu_result < HALF);
}

booleantype N_VConstrMask_Cuda(N_Vector C, N_Vector X, N_Vector M)
{
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstrMask_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
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
  realtype gpu_result = NVEC_CUDA_HBUFFERp(X)[0];

  return (gpu_result < HALF);
}

realtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom)
{
  // Starting value for min reduction
  const realtype maxVal = std::numeric_limits<realtype>::max();
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (InitializeReductionBuffer(num, maxVal))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMinQuotient_Cuda: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(num, true, grid, block, shMemSize, stream);
  minQuotientKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    maxVal,
    NVEC_CUDA_DDATAp(num),
    NVEC_CUDA_DDATAp(denom),
    NVEC_CUDA_DBUFFERp(num),
    NVEC_CUDA_CONTENT(num)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(num);
  realtype gpu_result = NVEC_CUDA_HBUFFERp(num)[0];

  return gpu_result;
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_DDATAp(X[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters and launch
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X[0], false, grid, block, shMemSize, stream)) return(-1);
  linearCombinationKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_c,
    d_Xd,
    NVEC_CUDA_DDATAp(Z),
    NVEC_CUDA_CONTENT(Z)->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;

  // Free device arrays
  err = cudaFree(d_c);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}

int N_VScaleAddMulti_Cuda(int nvec, realtype* c, N_Vector X, N_Vector* Y,
                          N_Vector* Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_CUDA_DDATAp(Y[i]);

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream)) return(-1);
  scaleAddMultiKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_c,
    NVEC_CUDA_DDATAp(X),
    d_Yd,
    d_Zd,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_c);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Yd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}

int N_VDotProdMulti_Cuda(int nvec, N_Vector X, N_Vector* Y, realtype* dots)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_CUDA_DDATAp(Y[i]);

  // Copy array of device pointers to device from host
  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream)) return(-1);
  grid = nvec;

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = cudaMalloc((void**) &d_buff, grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemsetAsync(d_buff, 0, grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  dotProdMultiKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    nvec,
    NVEC_CUDA_DDATAp(X),
    d_Yd,
    d_buff,
    NVEC_CUDA_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Copy GPU result to the cpu.
  err = cudaMemcpy(dots, d_buff, grid*sizeof(realtype), cudaMemcpyDeviceToHost);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Free host array
  delete[] h_Yd;

  // Free device arrays
  err = cudaFree(d_Yd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_buff);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Cuda(int nvec, realtype a, N_Vector* X, realtype b,
                                 N_Vector* Y, N_Vector* Z)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_DDATAp(X[i]);

  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_CUDA_DDATAp(Y[i]);

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream)) return(-1);
  linearSumVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    a,
    d_Xd,
    b,
    d_Yd,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Yd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}

int N_VScaleVectorArray_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_DDATAp(X[i]);

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream)) return(-1);
  scaleVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_c,
    d_Xd,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_c);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}

int N_VConstVectorArray_Cuda(int nvec, realtype c, N_Vector* Z)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream)) return(-1);
  constVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    c,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}

int N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                realtype* norms)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_DDATAp(X[i]);
  realtype** h_Wd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = NVEC_CUDA_DDATAp(W[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Wd;
  err = cudaMalloc((void**) &d_Wd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Wd, h_Wd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X[0], true, grid, block, shMemSize, stream)) return(-1);
  grid = nvec;

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = cudaMalloc((void**) &d_buff, grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemsetAsync(d_buff, 0, grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  wL2NormSquareVectorArrayKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_Xd,
    d_Wd,
    d_buff,
    NVEC_CUDA_CONTENT(X[0])->length
  );
  PostKernelLaunch();

  // Copy GPU result to the cpu.
  err = cudaMemcpy(norms, d_buff, grid*sizeof(realtype), cudaMemcpyDeviceToHost);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Finish computation
  for (int k=0; k<nvec; ++k)
    norms[k] = std::sqrt(norms[k]/NVEC_CUDA_CONTENT(X[0])->length);

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Wd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_buff);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}

int N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                    N_Vector id, realtype* norms)
{
  cudaError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_DDATAp(X[i]);

  realtype** h_Wd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = NVEC_CUDA_DDATAp(W[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Wd;
  err = cudaMalloc((void**) &d_Wd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Wd, h_Wd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(X[0], true, grid, block, shMemSize, stream)) return(-1);
  grid = nvec;

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = cudaMalloc((void**) &d_buff, grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemsetAsync(d_buff, 0, grid*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  wL2NormSquareMaskVectorArrayKernel<realtype, sunindextype><<<grid, block, shMemSize, stream>>>
  (
    nvec,
    d_Xd,
    d_Wd,
    NVEC_CUDA_DDATAp(id),
    d_buff,
    NVEC_CUDA_CONTENT(X[0])->length
  );
  PostKernelLaunch();

  // Copy GPU result to the cpu.
  err = cudaMemcpy(norms, d_buff, grid*sizeof(realtype), cudaMemcpyDeviceToHost);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Finish computation
  for (int k=0; k<nvec; ++k)
    norms[k] = std::sqrt(norms[k]/NVEC_CUDA_CONTENT(X[0])->length);

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Wd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_buff);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}

int N_VScaleAddMultiVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_CUDA_DDATAp(X[i]);

  realtype** h_Yd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Yd[j*nsum+i] = NVEC_CUDA_DDATAp(Y[i][j]);

  realtype** h_Zd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Zd[j*nsum+i] = NVEC_CUDA_DDATAp(Z[i][j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Yd, h_Yd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Zd, h_Zd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0][0], false, grid, block, shMemSize, stream)) return(-1);
  scaleAddMultiVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    nsum,
    d_c,
    d_Xd,
    d_Yd,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0][0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_c);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Yd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return(0);
}

int N_VLinearCombinationVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  cudaError_t err;

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Xd[j*nsum+i] = NVEC_CUDA_DDATAp(X[i][j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_CUDA_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Xd, h_Xd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  cudaStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream)) return(-1);
  linearCombinationVectorArrayKernel<<<grid, block, shMemSize, stream>>>
  (
    nvec,
    nsum,
    d_c,
    d_Xd,
    d_Zd,
    NVEC_CUDA_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_c);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Xd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);
  err = cudaFree(d_Zd);
  if (!SUNDIALS_CUDA_VERIFY(err)) return(-1);

  return cudaGetLastError();
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

  SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(x), buf_mem);

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

  SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(x), buf_mem);

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
  }

  /* return success */
  return(0);
}

int N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Cuda;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Cuda;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_Cuda;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Cuda;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Cuda;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Cuda;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Cuda;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Cuda;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Cuda;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Cuda;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
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
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_UVM);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed for SUNMEMTYPE_UVM\n");
    }
    vc->host_data = SUNMemoryHelper_Alias(vc->device_data);
  }
  else
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v), &(vc->host_data),
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_HOST);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
    }

    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v), &(vc->device_data),
                                       NVEC_CUDA_MEMSIZE(v), SUNMEMTYPE_DEVICE);
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
int InitializeReductionBuffer(N_Vector v, const realtype value)
{
  int alloc_fail = 0, copy_fail = 0;
  size_t bytes = sizeof(realtype);
  booleantype need_to_allocate = SUNFALSE;
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);
  SUNMemory value_mem = SUNMemoryHelper_Wrap((void*) &value, SUNMEMTYPE_HOST);

  /* we allocate if the existing reduction buffer is not large enough */
  if (vcp->reduce_buffer_allocated_bytes < bytes)
  {
    FreeReductionBuffer(v);
    need_to_allocate = SUNTRUE;
  }

  if (need_to_allocate)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->reduce_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("WARNING in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      /* try to allocate just plain host memory instead */
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                         &(vcp->reduce_buffer_host), bytes,
                                         SUNMEMTYPE_HOST);
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
      }
    }
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_CUDA_MEMHELP(v),
                                       &(vcp->reduce_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  if (!alloc_fail)
  {
    /* store the size of the buffer */
    vcp->reduce_buffer_allocated_bytes = bytes;

    /* initialize the memory with the value */
    copy_fail = SUNMemoryHelper_CopyAsync(NVEC_CUDA_MEMHELP(v),
                                          vcp->reduce_buffer_dev, value_mem,
                                          bytes, (void*) NVEC_CUDA_STREAM(v));

    if (copy_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_CopyAsync failed\n");
    }
  }

  SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), value_mem);
  return((alloc_fail || copy_fail) ? -1 : 0);
}

/* Free the reduction buffer
 */
void FreeReductionBuffer(N_Vector v)
{
  N_PrivateVectorContent_Cuda vcp = NVEC_CUDA_PRIVATE(v);

  if (vcp == NULL) return;

  if (vcp->reduce_buffer_dev != NULL)
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vcp->reduce_buffer_dev);
  vcp->reduce_buffer_dev  = NULL;
  if (vcp->reduce_buffer_host != NULL)
    SUNMemoryHelper_Dealloc(NVEC_CUDA_MEMHELP(v), vcp->reduce_buffer_host);
  vcp->reduce_buffer_host = NULL;
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
                                        n*sizeof(realtype),
                                        (void*) NVEC_CUDA_STREAM(v));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in CopyReductionBufferFromDevice: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = cudaStreamSynchronize(*NVEC_CUDA_STREAM(v));
  return (!SUNDIALS_CUDA_VERIFY(cuerr) || copy_fail ? -1 : 0);
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
