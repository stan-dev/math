/* ---------------------------------------------------------------------------
 * Programmer(s): Slaven Peles, Cody J. Balos, Daniel McGreer, and
 *                David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * This is the implementation file for a RAJA implementation of the NVECTOR
 * class with support for CUDA, HIP, and SYCL backends.
 * ---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <RAJA/RAJA.hpp>
#include <nvector/nvector_raja.h>

#include "sundials_debug.h"

// RAJA defines
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
#include <sunmemory/sunmemory_cuda.h>
#include "sundials_cuda.h"
#define SUNDIALS_RAJA_EXEC_STREAM RAJA::cuda_exec< 256 >
#define SUNDIALS_RAJA_EXEC_REDUCE RAJA::cuda_exec< 256 >
#define SUNDIALS_RAJA_REDUCE RAJA::cuda_reduce
#define SUNDIALS_GPU_PREFIX(val) cuda ## val
#define SUNDIALS_GPU_VERIFY SUNDIALS_CUDA_VERIFY
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
#include <sunmemory/sunmemory_hip.h>
#include "sundials_hip.h"
#define SUNDIALS_RAJA_EXEC_STREAM RAJA::hip_exec< 512 >
#define SUNDIALS_RAJA_EXEC_REDUCE RAJA::hip_exec< 512 >
#define SUNDIALS_RAJA_REDUCE RAJA::hip_reduce
#define SUNDIALS_GPU_PREFIX(val) hip ## val
#define SUNDIALS_GPU_VERIFY SUNDIALS_HIP_VERIFY
#elif defined(SUNDIALS_RAJA_BACKENDS_SYCL)
#include <sunmemory/sunmemory_sycl.h>
#include <CL/sycl.hpp>
#define SUNDIALS_RAJA_EXEC_STREAM RAJA::sycl_exec< 256 >
#define SUNDIALS_RAJA_EXEC_REDUCE RAJA::sycl_exec_nontrivial< 256 >
#define SUNDIALS_RAJA_REDUCE RAJA::sycl_reduce
#else
#error "Unknown RAJA backend"
#endif

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

// Static constants
static constexpr sunindextype zeroIdx = 0;

// Helpful macros
#define NVEC_RAJA_CONTENT(x) ((N_VectorContent_Raja)(x->content))
#define NVEC_RAJA_MEMSIZE(x) (NVEC_RAJA_CONTENT(x)->length * sizeof(realtype))
#define NVEC_RAJA_MEMHELP(x) (NVEC_RAJA_CONTENT(x)->mem_helper)
#define NVEC_RAJA_HDATAp(x)  ((realtype*) NVEC_RAJA_CONTENT(x)->host_data->ptr)
#define NVEC_RAJA_DDATAp(x)  ((realtype*) NVEC_RAJA_CONTENT(x)->device_data->ptr)

// Macros to access vector private content
#define NVEC_RAJA_PRIVATE(x)  ((N_PrivateVectorContent_Raja)(NVEC_RAJA_CONTENT(x)->priv))
#define NVEC_RAJA_HBUFFERp(x) ((realtype*) NVEC_RAJA_PRIVATE(x)->reduce_buffer_host->ptr)
#define NVEC_RAJA_DBUFFERp(x) ((realtype*) NVEC_RAJA_PRIVATE(x)->reduce_buffer_dev->ptr)

/*
 * Private structure definition
 */

struct _N_PrivateVectorContent_Raja
{
  booleantype use_managed_mem; /* do data pointers use managed memory */

  /* fused op workspace */
  SUNMemory fused_buffer_dev;    /* device memory for fused ops    */
  SUNMemory fused_buffer_host;   /* host memory for fused ops      */
  size_t    fused_buffer_bytes;  /* current size of the buffers    */
  size_t    fused_buffer_offset; /* current offset into the buffer */
};

typedef struct _N_PrivateVectorContent_Raja *N_PrivateVectorContent_Raja;


/*
 * Utility functions
 */


// Allocate vector data
static int AllocateData(N_Vector v);

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


N_Vector N_VNewEmpty_Raja(SUNContext sunctx)
{
  N_Vector v;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
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

  NVEC_RAJA_CONTENT(v)->length      = 0;
  NVEC_RAJA_CONTENT(v)->mem_helper  = NULL;
  NVEC_RAJA_CONTENT(v)->own_helper  = SUNFALSE;
  NVEC_RAJA_CONTENT(v)->host_data   = NULL;
  NVEC_RAJA_CONTENT(v)->device_data = NULL;

  NVEC_RAJA_PRIVATE(v)->use_managed_mem      = SUNFALSE;
  NVEC_RAJA_PRIVATE(v)->fused_buffer_dev     = NULL;
  NVEC_RAJA_PRIVATE(v)->fused_buffer_host    = NULL;
  NVEC_RAJA_PRIVATE(v)->fused_buffer_bytes   = 0;
  NVEC_RAJA_PRIVATE(v)->fused_buffer_offset  = 0;

  return(v);
}

N_Vector N_VNew_Raja(sunindextype length, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja(sunctx);
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda(sunctx);
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Hip(sunctx);
#elif defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Sycl(sunctx);
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

N_Vector N_VNewWithMemHelp_Raja(sunindextype length,
                                booleantype use_managed_mem,
                                SUNMemoryHelper helper,
                                SUNContext sunctx)
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
  v = N_VNewEmpty_Raja(sunctx);
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

N_Vector N_VNewManaged_Raja(sunindextype length, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja(sunctx);
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda(sunctx);
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Hip(sunctx);
#elif defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Sycl(sunctx);
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

N_Vector N_VMake_Raja(sunindextype length, realtype *h_vdata, realtype *d_vdata,
                      SUNContext sunctx)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Raja(sunctx);
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->host_data       = SUNMemoryHelper_Wrap(h_vdata, SUNMEMTYPE_HOST);
  NVEC_RAJA_CONTENT(v)->device_data     = SUNMemoryHelper_Wrap(d_vdata, SUNMEMTYPE_DEVICE);
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda(sunctx);
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Hip(sunctx);
#elif defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Sycl(sunctx);
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

N_Vector N_VMakeManaged_Raja(sunindextype length, realtype *vdata,
                             SUNContext sunctx)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Raja(sunctx);
  if (v == NULL) return(NULL);

  NVEC_RAJA_CONTENT(v)->length          = length;
  NVEC_RAJA_CONTENT(v)->host_data       = SUNMemoryHelper_Wrap(vdata, SUNMEMTYPE_UVM);
  NVEC_RAJA_CONTENT(v)->device_data     = SUNMemoryHelper_Alias(NVEC_RAJA_CONTENT(v)->host_data);
#if defined(SUNDIALS_RAJA_BACKENDS_CUDA)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Cuda(sunctx);
#elif defined(SUNDIALS_RAJA_BACKENDS_HIP)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Hip(sunctx);
#elif defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  NVEC_RAJA_CONTENT(v)->mem_helper      = SUNMemoryHelper_Sycl(sunctx);
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

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
  void* queue = nullptr;
#endif

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(x),
                                        NVEC_RAJA_CONTENT(x)->device_data,
                                        NVEC_RAJA_CONTENT(x)->host_data,
                                        NVEC_RAJA_MEMSIZE(x),
                                        queue);

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyToDevice_Raja: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but on the default stream currently */
#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  ::sycl::queue* q = ::RAJA::sycl::detail::getQueue();
  q->wait_and_throw();
#else
  SUNDIALS_GPU_VERIFY(SUNDIALS_GPU_PREFIX(StreamSynchronize)(0));
#endif
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Raja(N_Vector x)
{
  int copy_fail;

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
  void* queue = nullptr;
#endif

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(x),
                                        NVEC_RAJA_CONTENT(x)->host_data,
                                        NVEC_RAJA_CONTENT(x)->device_data,
                                        NVEC_RAJA_MEMSIZE(x),
                                        queue);

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyFromDevice_Raja: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  ::sycl::queue* q = ::RAJA::sycl::detail::getQueue();
  q->wait_and_throw();
#else
  SUNDIALS_GPU_VERIFY(SUNDIALS_GPU_PREFIX(StreamSynchronize)(0));
#endif
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
  v = N_VNewEmpty_Raja(w->sunctx);
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
    FusedBuffer_Free(v);
    free(vcp);
    vc->priv = NULL;
  }

  /* free items in content */
  if (NVEC_RAJA_MEMHELP(v))
  {
#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
    void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
    void* queue = nullptr;
#endif

    SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(v), vc->host_data, queue);
    vc->host_data = NULL;
    SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(v), vc->device_data, queue);
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

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N), [=] RAJA_DEVICE (sunindextype i) {
     zdata[i] = c;
  });
}

void N_VLinearSum_Raja(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const realtype *ydata = NVEC_RAJA_DDATAp(Y);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      zdata[i] = xdata[i] / ydata[i];
    }
  );
}

void N_VScale_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      zdata[i] = c * xdata[i];
    }
  );
}

void N_VAbs_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      zdata[i] = abs(xdata[i]);
    }
  );
}

void N_VInv_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      zdata[i] = ONE / xdata[i];
    }
  );
}

void N_VAddConst_Raja(N_Vector X, realtype b, N_Vector Z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(Z);

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      zdata[i] = xdata[i] + b;
    }
  );
}

realtype N_VDotProd_Raja(N_Vector X, N_Vector Y)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const realtype *ydata = NVEC_RAJA_DDATAp(Y);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceSum< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(0.0);
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      gpu_result += xdata[i] * ydata[i] ;
    }
  );

  return (static_cast<realtype>(gpu_result));
}

realtype N_VMaxNorm_Raja(N_Vector X)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(X);
  const sunindextype N = NVEC_RAJA_CONTENT(X)->length;

  RAJA::ReduceMax< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(0.0);
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::ReduceSum< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(0.0);
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::ReduceSum< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(0.0);
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::ReduceMin< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::ReduceSum< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(0.0);
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
    }
  );
}

booleantype N_VInvTest_Raja(N_Vector x, N_Vector z)
{
  const realtype *xdata = NVEC_RAJA_DDATAp(x);
  const sunindextype N = NVEC_RAJA_CONTENT(x)->length;
  realtype *zdata = NVEC_RAJA_DDATAp(z);

  RAJA::ReduceSum< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(ZERO);
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::ReduceSum< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(ZERO);
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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

  RAJA::ReduceMin< SUNDIALS_RAJA_REDUCE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< SUNDIALS_RAJA_EXEC_REDUCE >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
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
  const sunindextype N     = NVEC_RAJA_CONTENT(z)->length;
  realtype*          zdata = NVEC_RAJA_DDATAp(z);

  // Fused op workspace shortcuts
  realtype*  cdata = NULL;
  realtype** xdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(z, nvec, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Raja: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(z, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Raja: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(z, X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Raja: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(z))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Raja: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      zdata[i] = cdata[0] * xdata[0][i];
      for (int j=1; j<nvec; j++)
        zdata[i] += cdata[j] * xdata[j][i];
    }
  );

  return 0;
}


int N_VScaleAddMulti_Raja(int nvec, realtype* c, N_Vector x, N_Vector* Y,
                          N_Vector* Z)
{
  const sunindextype N      = NVEC_RAJA_CONTENT(x)->length;
  const realtype     *xdata = NVEC_RAJA_DDATAp(x);

  // Shortcuts to the fused op workspace
  realtype*  cdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(x, nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Raja: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(x, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Raja: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Raja: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Raja: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Raja: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
     [=] RAJA_DEVICE (sunindextype i) {
      for (int j=0; j<nvec; j++)
        zdata[j][i] = cdata[j] * xdata[i] + ydata[j][i];
    }
  );

  return 0;
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
  const sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Shortcuts to the fused op workspace
  realtype** xdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], 0, 3 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Sycl: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSumVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinaerSumVectorArray_Sycl: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      for (int j=0; j<nvec; j++)
        zdata[j][i] = a * xdata[j][i] + b * ydata[j][i];
    }
  );

  return 0;
}


int N_VScaleVectorArray_Raja(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  const sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Shortcuts to the fused op workspace arrays
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Sycl: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(Z[0], c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Sycl: FusedBuffer_CopyReadArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Sycl: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      for (int j=0; j<nvec; j++)
        zdata[j][i] = cdata[j] * xdata[j][i];
    }
  );

  return 0;
}


int N_VConstVectorArray_Raja(int nvec, realtype c, N_Vector* Z)
{
  const sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Shortcuts to the fused op workspace arrays
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], 0, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstVectorArray_Sycl: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstVectorArray_Sycl: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      for (int j=0; j<nvec; j++)
        zdata[j][i] = c;
    }
  );

  return 0;
}


int N_VScaleAddMultiVectorArray_Raja(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  const sunindextype N = NVEC_RAJA_CONTENT(X[0])->length;

  // Shortcuts to the fused op workspace
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(X[0], nsum, nvec + 2 * nvec * nsum))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiArray_Sycl: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(X[0], c, nsum, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiArray_Sycl: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(X[0], X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray2D(X[0], Y, nvec, nsum, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiVectorArray_Sycl: FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray2D(X[0], Z, nvec, nsum, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiVectorArray_Sycl: FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(X[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Sycl: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      for (int j=0; j<nvec; j++)
        for (int k=0; k<nsum; k++)
          zdata[j * nsum + k][i] =
            cdata[k] * xdata[j][i] + ydata[j * nsum + k][i];
    }
  );

  return 0;
}


int N_VLinearCombinationVectorArray_Raja(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  const sunindextype N = NVEC_RAJA_CONTENT(Z[0])->length;

  // Shortcuts to the fused op workspace arrays
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** zdata = NULL;

  // Setup the fused op workspace
  if (FusedBuffer_Init(Z[0], nsum, nvec + nvec * nsum))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Sycl: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(Z[0], c, nsum, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Sycl: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray2D(Z[0], X, nvec, nsum, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Sycl: FusedBuffer_CopyPtrArray2D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(Z[0], Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(Z[0]))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Sycl: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  RAJA::forall< SUNDIALS_RAJA_EXEC_STREAM >(RAJA::RangeSegment(zeroIdx, N),
    [=] RAJA_DEVICE (sunindextype i) {
      for (int j=0; j<nvec; j++) {
        zdata[j][i] = cdata[0] * xdata[j * nsum][i];
        for (int k=1; k<nsum; k++) {
          zdata[j][i] += cdata[k] * xdata[j * nsum + k][i];
        }
      }
    }
  );

  return 0;
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
#if !defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  SUNDIALS_GPU_PREFIX(Error_t) cuerr;
#endif

  if (x == NULL || buf == NULL) return(-1);

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return(-1);

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
  void* queue = nullptr;
#endif

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(x),
                                        buf_mem,
                                        NVEC_RAJA_CONTENT(x)->device_data,
                                        NVEC_RAJA_MEMSIZE(x),
                                        queue);

  /* we synchronize with respect to the host, but only in this stream */
#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  ::sycl::queue* q = ::RAJA::sycl::detail::getQueue();
  q->wait_and_throw();
#else
  cuerr = SUNDIALS_GPU_PREFIX(StreamSynchronize)(0);
#endif

  SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(x), buf_mem, queue);

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  return (copy_fail ? -1 : 0);
#else
  return (!SUNDIALS_GPU_VERIFY(cuerr) || copy_fail ? -1 : 0);
#endif
}


int N_VBufUnpack_Raja(N_Vector x, void *buf)
{
  int copy_fail = 0;
#if !defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  SUNDIALS_GPU_PREFIX(Error_t) cuerr;
#endif

  if (x == NULL || buf == NULL) return(-1);

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return(-1);

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
  void* queue = nullptr;
#endif

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(x),
                                        NVEC_RAJA_CONTENT(x)->device_data,
                                        buf_mem,
                                        NVEC_RAJA_MEMSIZE(x),
                                        queue);

  /* we synchronize with respect to the host, but only in this stream */
#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  ::sycl::queue* q = ::RAJA::sycl::detail::getQueue();
  q->wait_and_throw();
#else
  cuerr = SUNDIALS_GPU_PREFIX(StreamSynchronize)(0);
#endif

  SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(x), buf_mem, queue);

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  return (copy_fail ? -1 : 0);
#else
  return (!SUNDIALS_GPU_VERIFY(cuerr) || copy_fail ? -1 : 0);
#endif
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

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
  void* queue = nullptr;
#endif

  if (vcp->use_managed_mem)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v), &(vc->device_data),
                                       NVEC_RAJA_MEMSIZE(v), SUNMEMTYPE_UVM,
                                       queue);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed for SUNMEMTYPE_UVM\n");
    }
    vc->host_data = SUNMemoryHelper_Alias(vc->device_data);
  }
  else
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v), &(vc->host_data),
                                       NVEC_RAJA_MEMSIZE(v), SUNMEMTYPE_HOST,
                                       queue);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
    }

    alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v), &(vc->device_data),
                                       NVEC_RAJA_MEMSIZE(v), SUNMEMTYPE_DEVICE,
                                       queue);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  return(alloc_fail ? -1 : 0);
}


static int FusedBuffer_Init(N_Vector v, int nreal, int nptr)
{
  int         alloc_fail = 0;
  booleantype alloc_mem  = SUNFALSE;
  size_t      bytes      = nreal * sizeof(realtype) + nptr * sizeof(realtype*);

  // Get the vector private memory structure
  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

  // Check if the existing memory is not large enough
  if (vcp->fused_buffer_bytes < bytes)
  {
    FusedBuffer_Free(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
    void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
    void* queue = nullptr;
#endif

    // allocate pinned memory on the host
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v),
                                       &(vcp->fused_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED, queue);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("WARNING in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      // if pinned alloc failed, allocate plain host memory
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v),
                                         &(vcp->fused_buffer_host), bytes,
                                         SUNMEMTYPE_HOST, queue);
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
        return -1;
      }
    }

    // allocate device memory
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_RAJA_MEMHELP(v),
                                       &(vcp->fused_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE, queue);
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
  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

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

  vcp->fused_buffer_offset += nval * sizeof(realtype);

  return 0;
}


static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector *X, int nvec,
                                      realtype ***shortcut)
{
  // Get the vector private memory structure
  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

  // Check buffer space and fill the host buffer
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyPtrArray1D: Buffer offset is exceedes the buffer size\n");    return -1;
    return -1;
  }

  realtype** h_buffer = (realtype**) ((char*)(vcp->fused_buffer_host->ptr) +
                                      vcp->fused_buffer_offset);

  for (int j = 0; j < nvec; j++)
  {
    h_buffer[j] = NVEC_RAJA_DDATAp(X[j]);
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
  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

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
      h_buffer[j * nsum + k] = NVEC_RAJA_DDATAp(X[k][j]);
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
  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
  void* queue = nullptr;
#endif

  // Copy the fused buffer to the device
  int copy_fail = SUNMemoryHelper_CopyAsync(NVEC_RAJA_MEMHELP(v),
                                            vcp->fused_buffer_dev,
                                            vcp->fused_buffer_host,
                                            vcp->fused_buffer_offset,
                                            queue);
  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyToDevice: SUNMemoryHelper_CopyAsync failed\n");
    return -1;
  }

  return 0;
}


static int FusedBuffer_Free(N_Vector v)
{
  N_PrivateVectorContent_Raja vcp = NVEC_RAJA_PRIVATE(v);

  if (vcp == NULL) return 0;

#if defined(SUNDIALS_RAJA_BACKENDS_SYCL)
  void* queue = static_cast<void*>(::RAJA::sycl::detail::getQueue());
#else
  void* queue = nullptr;
#endif

  if (vcp->fused_buffer_host)
  {
    SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(v),
                            vcp->fused_buffer_host, queue);
    vcp->fused_buffer_host = NULL;
  }

  if (vcp->fused_buffer_dev)
  {
    SUNMemoryHelper_Dealloc(NVEC_RAJA_MEMHELP(v),
                            vcp->fused_buffer_dev, queue);
    vcp->fused_buffer_dev = NULL;
  }

  vcp->fused_buffer_bytes  = 0;
  vcp->fused_buffer_offset = 0;

  return 0;
}


} // extern "C"
