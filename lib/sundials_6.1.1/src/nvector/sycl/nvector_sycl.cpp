/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the implementation file for a SYCL implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <CL/sycl.hpp>

/* SUNDIALS public headers */
#include <nvector/nvector_sycl.h>
#include <sunmemory/sunmemory_sycl.h>

/* SUNDIALS private headers */
#include "sundials_debug.h"
#include "sundials_sycl.h"

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

using namespace sundials;
using namespace sundials::sycl;


/* --------------------------------------------------------------------------
 * Helpful macros
 * -------------------------------------------------------------------------- */


/* Macros to access vector content */
#define NVEC_SYCL_CONTENT(x)  ((N_VectorContent_Sycl)(x->content))
#define NVEC_SYCL_LENGTH(x)   (NVEC_SYCL_CONTENT(x)->length)
#define NVEC_SYCL_MEMHELP(x)  (NVEC_SYCL_CONTENT(x)->mem_helper)
#define NVEC_SYCL_MEMSIZE(x)  (NVEC_SYCL_CONTENT(x)->length * sizeof(realtype))
#define NVEC_SYCL_HDATAp(x)   ((realtype*) NVEC_SYCL_CONTENT(x)->host_data->ptr)
#define NVEC_SYCL_DDATAp(x)   ((realtype*) NVEC_SYCL_CONTENT(x)->device_data->ptr)
#define NVEC_SYCL_QUEUE(x)    (NVEC_SYCL_CONTENT(x)->queue)

/* Macros to access vector private content */
#define NVEC_SYCL_PRIVATE(x)      ((N_PrivateVectorContent_Sycl)(NVEC_SYCL_CONTENT(x)->priv))
#define NVEC_SYCL_HBUFFERp(x)     ((realtype*) NVEC_SYCL_PRIVATE(x)->reduce_buffer_host->ptr)
#define NVEC_SYCL_DBUFFERp(x)     ((realtype*) NVEC_SYCL_PRIVATE(x)->reduce_buffer_dev->ptr)


/* --------------------------------------------------------------------------
 * Private structure definition
 * -------------------------------------------------------------------------- */


struct _N_PrivateVectorContent_Sycl
{
  booleantype use_managed_mem; /* do data pointers use managed memory */

  /* reduction workspace */
  SUNMemory reduce_buffer_dev;   /* device memory for reductions      */
  SUNMemory reduce_buffer_host;  /* host memory for reductions        */
  size_t    reduce_buffer_bytes; /* current size of reduction buffers */

  /* fused op workspace */
  SUNMemory fused_buffer_dev;    /* device memory for fused ops    */
  SUNMemory fused_buffer_host;   /* host memory for fused ops      */
  size_t    fused_buffer_bytes;  /* current size of the buffers    */
  size_t    fused_buffer_offset; /* current offset into the buffer */
};

typedef struct _N_PrivateVectorContent_Sycl *N_PrivateVectorContent_Sycl;


/* --------------------------------------------------------------------------
 * Utility functions
 * -------------------------------------------------------------------------- */


/* Allocate vector data */
static int AllocateData(N_Vector v);

/* Reduction buffer functions */
static int InitializeReductionBuffer(N_Vector v, const realtype value,
                                     size_t n = 1);
static void FreeReductionBuffer(N_Vector v);
static int CopyReductionBufferFromDevice(N_Vector v, size_t n = 1);

/* Fused operation buffer functions */
static int FusedBuffer_Init(N_Vector v, int nreal, int nptr);
static int FusedBuffer_CopyRealArray(N_Vector v, realtype *r_data, int nval,
                                     realtype **shortcut);
static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector *X, int nvec,
                                      realtype ***shortcut);
static int FusedBuffer_CopyPtrArray2D(N_Vector v, N_Vector **X, int nvec,
                                      int nsum, realtype ***shortcut);
static int FusedBuffer_CopyToDevice(N_Vector v);
static int FusedBuffer_Free(N_Vector v);

/* Kernel launch parameters */
static int GetKernelParameters(N_Vector v, booleantype reduction,
                               size_t& nthreads_total,
                               size_t& nthreads_per_block);


/* --------------------------------------------------------------------------
 * Constructors
 * -------------------------------------------------------------------------- */


N_Vector N_VNewEmpty_Sycl(SUNContext sunctx)
{
  /* Create an empty vector object */
  N_Vector v = N_VNewEmpty(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewEmpty_Sycl: N_VNewEmpty returned NULL\n");
    return NULL;
  }

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = N_VGetVectorID_Sycl;
  v->ops->nvclone                 = N_VClone_Sycl;
  v->ops->nvcloneempty            = N_VCloneEmpty_Sycl;
  v->ops->nvdestroy               = N_VDestroy_Sycl;
  v->ops->nvspace                 = N_VSpace_Sycl;
  v->ops->nvgetlength             = N_VGetLength_Sycl;
  v->ops->nvgetarraypointer       = N_VGetHostArrayPointer_Sycl;
  v->ops->nvgetdevicearraypointer = N_VGetDeviceArrayPointer_Sycl;
  v->ops->nvsetarraypointer       = N_VSetHostArrayPointer_Sycl;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Sycl;
  v->ops->nvconst        = N_VConst_Sycl;
  v->ops->nvprod         = N_VProd_Sycl;
  v->ops->nvdiv          = N_VDiv_Sycl;
  v->ops->nvscale        = N_VScale_Sycl;
  v->ops->nvabs          = N_VAbs_Sycl;
  v->ops->nvinv          = N_VInv_Sycl;
  v->ops->nvaddconst     = N_VAddConst_Sycl;
  v->ops->nvdotprod      = N_VDotProd_Sycl;
  v->ops->nvmaxnorm      = N_VMaxNorm_Sycl;
  v->ops->nvmin          = N_VMin_Sycl;
  v->ops->nvl1norm       = N_VL1Norm_Sycl;
  v->ops->nvinvtest      = N_VInvTest_Sycl;
  v->ops->nvconstrmask   = N_VConstrMask_Sycl;
  v->ops->nvminquotient  = N_VMinQuotient_Sycl;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Sycl;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Sycl;
  v->ops->nvwl2norm      = N_VWL2Norm_Sycl;
  v->ops->nvcompare      = N_VCompare_Sycl;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Sycl;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Sycl;
  v->ops->nvdotprodlocal     = N_VDotProd_Sycl;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Sycl;
  v->ops->nvminlocal         = N_VMin_Sycl;
  v->ops->nvl1normlocal      = N_VL1Norm_Sycl;
  v->ops->nvinvtestlocal     = N_VInvTest_Sycl;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Sycl;
  v->ops->nvminquotientlocal = N_VMinQuotient_Sycl;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_Sycl;
  v->ops->nvbufpack   = N_VBufPack_Sycl;
  v->ops->nvbufunpack = N_VBufUnpack_Sycl;

  /* print operation for debugging */
  v->ops->nvprint     = N_VPrint_Sycl;
  v->ops->nvprintfile = N_VPrintFile_Sycl;

  /* Allocate content structure */
  v->content = (N_VectorContent_Sycl) malloc(sizeof(_N_VectorContent_Sycl));
  if (v->content == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewEmpty_Sycl: content allocation failed\n");
    N_VDestroy(v);
    return NULL;
  }

  /* Allocate private content structure */
  NVEC_SYCL_CONTENT(v)->priv = NULL;
  NVEC_SYCL_CONTENT(v)->priv = malloc(sizeof(_N_PrivateVectorContent_Sycl));
  if (NVEC_SYCL_CONTENT(v)->priv == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewEmpty_Sycl: private content allocation failed\n");
    N_VDestroy(v);
    return NULL;
  }

  /* Initialize content */
  NVEC_SYCL_CONTENT(v)->length             = 0;
  NVEC_SYCL_CONTENT(v)->own_exec           = SUNFALSE;
  NVEC_SYCL_CONTENT(v)->own_helper         = SUNFALSE;
  NVEC_SYCL_CONTENT(v)->host_data          = NULL;
  NVEC_SYCL_CONTENT(v)->device_data        = NULL;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy = NULL;
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy = NULL;
  NVEC_SYCL_CONTENT(v)->mem_helper         = NULL;
  NVEC_SYCL_CONTENT(v)->queue              = NULL;

  /* Initialize private content */
  NVEC_SYCL_PRIVATE(v)->use_managed_mem      = SUNFALSE;
  NVEC_SYCL_PRIVATE(v)->reduce_buffer_dev    = NULL;
  NVEC_SYCL_PRIVATE(v)->reduce_buffer_host   = NULL;
  NVEC_SYCL_PRIVATE(v)->reduce_buffer_bytes  = 0;
  NVEC_SYCL_PRIVATE(v)->fused_buffer_dev     = NULL;
  NVEC_SYCL_PRIVATE(v)->fused_buffer_host    = NULL;
  NVEC_SYCL_PRIVATE(v)->fused_buffer_bytes   = 0;
  NVEC_SYCL_PRIVATE(v)->fused_buffer_offset  = 0;

  return v;
}


N_Vector N_VNew_Sycl(sunindextype length, ::sycl::queue *Q, SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Sycl: queue is not in-order\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length             = length;
  NVEC_SYCL_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->mem_helper         = SUNMemoryHelper_Sycl(sunctx);
  NVEC_SYCL_CONTENT(v)->queue              = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem    = SUNFALSE;

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Sycl: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}


N_Vector N_VNewWithMemHelp_Sycl(sunindextype length,
                                booleantype use_managed_mem,
                                SUNMemoryHelper helper,
                                ::sycl::queue *Q,
                                SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: queue is not in-order\n");
    return NULL;
  }

  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: helper is NULL\n");
    return NULL;
  }

  if (!SUNMemoryHelper_ImplementsRequiredOps(helper))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: helper doesn't implement all required ops\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length             = length;
  NVEC_SYCL_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->own_helper         = SUNFALSE;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->mem_helper         = helper;
  NVEC_SYCL_CONTENT(v)->queue              = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem    = use_managed_mem;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Sycl: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}


N_Vector N_VNewManaged_Sycl(sunindextype length, ::sycl::queue *Q,
                            SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Sycl: queue is not in-order\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length             = length;
  NVEC_SYCL_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->mem_helper         = SUNMemoryHelper_Sycl(sunctx);
  NVEC_SYCL_CONTENT(v)->queue              = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem    = SUNTRUE;

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Sycl: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}


N_Vector N_VMake_Sycl(sunindextype length, realtype *h_vdata, realtype *d_vdata,
                      ::sycl::queue *Q, SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: queue is not in-order\n");
    return NULL;
  }

  if (h_vdata == NULL || d_vdata == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: host or device data is NULL\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length             = length;
  NVEC_SYCL_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->host_data          = SUNMemoryHelper_Wrap(h_vdata, SUNMEMTYPE_HOST);
  NVEC_SYCL_CONTENT(v)->device_data        = SUNMemoryHelper_Wrap(d_vdata, SUNMEMTYPE_DEVICE);
  NVEC_SYCL_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->mem_helper         = SUNMemoryHelper_Sycl(sunctx);
  NVEC_SYCL_CONTENT(v)->queue              = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem    = SUNFALSE;

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (NVEC_SYCL_CONTENT(v)->device_data == NULL ||
      NVEC_SYCL_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Sycl: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}


N_Vector N_VMakeManaged_Sycl(sunindextype length, realtype *vdata,
                             ::sycl::queue *Q, SUNContext sunctx)
{
  /* Check inputs */
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Sycl: queue is NULL\n");
    return NULL;
  }

  if (!(Q->is_in_order()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Sycl: queue is not in-order\n");
    return NULL;
  }

  if (vdata == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Sycl: host or device data is NULL\n");
    return NULL;
  }

  /* Create vector with empty content */
  N_Vector v = N_VNewEmpty_Sycl(sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Sycl: N_VNewEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Fill content */
  NVEC_SYCL_CONTENT(v)->length             = length;
  NVEC_SYCL_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->host_data          = SUNMemoryHelper_Wrap(vdata, SUNMEMTYPE_UVM);
  NVEC_SYCL_CONTENT(v)->device_data        = SUNMemoryHelper_Alias(NVEC_SYCL_CONTENT(v)->host_data);
  NVEC_SYCL_CONTENT(v)->stream_exec_policy = new ThreadDirectExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy = new BlockReduceExecPolicy(SYCL_BLOCKDIM(Q));
  NVEC_SYCL_CONTENT(v)->mem_helper         = SUNMemoryHelper_Sycl(sunctx);
  NVEC_SYCL_CONTENT(v)->queue              = Q;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem    = SUNTRUE;

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (NVEC_SYCL_CONTENT(v)->device_data == NULL ||
      NVEC_SYCL_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Sycl: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}


/* --------------------------------------------------------------------------
 * Vector get, set, and utility functions
 * -------------------------------------------------------------------------- */


/* Function to return the global length of the vector. This is defined as an
 * inline function in nvector_sycl.h, so we just mark it as extern here. */
extern sunindextype N_VGetLength_Sycl(N_Vector v);


/* Return pointer to the raw host data. This is defined as an inline function in
 * nvector_sycl.h, so we just mark it as extern here. */
extern realtype *N_VGetHostArrayPointer_Sycl(N_Vector x);


/* Return pointer to the raw device data. This is defined as an inline function
 * in nvector_sycl.h, so we just mark it as extern here. */
extern realtype *N_VGetDeviceArrayPointer_Sycl(N_Vector x);


/* Set pointer to the raw host data. Does not free the existing pointer. */
void N_VSetHostArrayPointer_Sycl(realtype* h_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Sycl(v))
  {
    if (NVEC_SYCL_CONTENT(v)->host_data)
    {
      NVEC_SYCL_CONTENT(v)->host_data->ptr = (void*) h_vdata;
      NVEC_SYCL_CONTENT(v)->device_data->ptr = (void*) h_vdata;
    }
    else
    {
      NVEC_SYCL_CONTENT(v)->host_data = SUNMemoryHelper_Wrap((void*) h_vdata, SUNMEMTYPE_UVM);
      NVEC_SYCL_CONTENT(v)->device_data = SUNMemoryHelper_Alias(NVEC_SYCL_CONTENT(v)->host_data);
    }
  }
  else
  {
    if (NVEC_SYCL_CONTENT(v)->host_data)
    {
      NVEC_SYCL_CONTENT(v)->host_data->ptr = (void*) h_vdata;
    }
    else
    {
      NVEC_SYCL_CONTENT(v)->host_data = SUNMemoryHelper_Wrap((void*) h_vdata, SUNMEMTYPE_HOST);
    }
  }
}


/* Set pointer to the raw device data */
void N_VSetDeviceArrayPointer_Sycl(realtype* d_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Sycl(v))
  {
    if (NVEC_SYCL_CONTENT(v)->device_data)
    {
      NVEC_SYCL_CONTENT(v)->device_data->ptr = (void*) d_vdata;
      NVEC_SYCL_CONTENT(v)->host_data->ptr = (void*) d_vdata;
    }
    else
    {
      NVEC_SYCL_CONTENT(v)->device_data = SUNMemoryHelper_Wrap((void*) d_vdata, SUNMEMTYPE_UVM);
      NVEC_SYCL_CONTENT(v)->host_data = SUNMemoryHelper_Alias(NVEC_SYCL_CONTENT(v)->device_data);
    }
  }
  else
  {
    if (NVEC_SYCL_CONTENT(v)->device_data)
    {
      NVEC_SYCL_CONTENT(v)->device_data->ptr = (void*) d_vdata;
    }
    else
    {
      NVEC_SYCL_CONTENT(v)->device_data = SUNMemoryHelper_Wrap((void*) d_vdata, SUNMEMTYPE_DEVICE);
    }
  }
}


/* Return a flag indicating if the memory for the vector data is managed */
booleantype N_VIsManagedMemory_Sycl(N_Vector x)
{
  return NVEC_SYCL_PRIVATE(x)->use_managed_mem;
}


int N_VSetKernelExecPolicy_Sycl(N_Vector x,
                                SUNSyclExecPolicy* stream_exec_policy,
                                SUNSyclExecPolicy* reduce_exec_policy)
{
  if (x == NULL || stream_exec_policy == NULL || reduce_exec_policy == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VSetKernelExecPolicy_Sycl: An input is NULL\n");
    return -1;
  }

  if (NVEC_SYCL_CONTENT(x)->own_exec)
  {
    delete NVEC_SYCL_CONTENT(x)->stream_exec_policy;
    delete NVEC_SYCL_CONTENT(x)->reduce_exec_policy;
  }

  NVEC_SYCL_CONTENT(x)->stream_exec_policy = stream_exec_policy;
  NVEC_SYCL_CONTENT(x)->reduce_exec_policy = reduce_exec_policy;
  NVEC_SYCL_CONTENT(x)->own_exec = SUNFALSE;

  return 0;
}


/* Copy vector data to the device */
void N_VCopyToDevice_Sycl(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_Copy(NVEC_SYCL_MEMHELP(x),
                                   NVEC_SYCL_CONTENT(x)->device_data,
                                   NVEC_SYCL_CONTENT(x)->host_data,
                                   NVEC_SYCL_MEMSIZE(x),
                                   NVEC_SYCL_QUEUE(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyToDevice_Sycl: SUNMemoryHelper_Copy returned nonzero\n");
  }

  /* synchronize with the host */
  NVEC_SYCL_QUEUE(x)->wait_and_throw();
}


/* Copy vector data from the device to the host */
void N_VCopyFromDevice_Sycl(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_Copy(NVEC_SYCL_MEMHELP(x),
                                   NVEC_SYCL_CONTENT(x)->host_data,
                                   NVEC_SYCL_CONTENT(x)->device_data,
                                   NVEC_SYCL_MEMSIZE(x),
                                   NVEC_SYCL_QUEUE(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyFromDevice_Sycl: SUNMemoryHelper_Copy returned nonzero\n");
  }

  /* synchronize with the host */
  NVEC_SYCL_QUEUE(x)->wait_and_throw();
}


/* Function to print the a serial vector to stdout */
void N_VPrint_Sycl(N_Vector X)
{
  N_VPrintFile_Sycl(X, stdout);
}


/* Function to print the a serial vector to outfile */
void N_VPrintFile_Sycl(N_Vector X, FILE *outfile)
{
  sunindextype i;

  for (i = 0; i < NVEC_SYCL_CONTENT(X)->length; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", NVEC_SYCL_HDATAp(X)[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", NVEC_SYCL_HDATAp(X)[i]);
#else
    fprintf(outfile, "%11.8g\n", NVEC_SYCL_HDATAp(X)[i]);
#endif
  }
  fprintf(outfile, "\n");

  return;
}


/* --------------------------------------------------------------------------
 * Vector operations
 * -------------------------------------------------------------------------- */


N_Vector N_VCloneEmpty_Sycl(N_Vector w)
{
  /* Check input */
  if (w == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCloneEmpty_Sycl: input vector is NULL\n");
    return NULL;
  }

  /* Create vector */
  N_Vector v = N_VNewEmpty_Sycl(w->sunctx);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCloneEmpty_Sycl: N_VNewEmpty returned NULL\n");
    return NULL;
  }

  /* Attach operations */
  if (N_VCopyOps(w, v))
  {
    N_VDestroy(v);
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCloneEmpty_Sycl: Error in N_VCopyOps\n");
    return NULL;
  }

  /* Copy some content */
  NVEC_SYCL_CONTENT(v)->length          = NVEC_SYCL_CONTENT(w)->length;
  NVEC_SYCL_CONTENT(v)->queue           = NVEC_SYCL_CONTENT(w)->queue;
  NVEC_SYCL_PRIVATE(v)->use_managed_mem = NVEC_SYCL_PRIVATE(w)->use_managed_mem;

  return v;
}


N_Vector N_VClone_Sycl(N_Vector w)
{
  /* Check inputs */
  if (w == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Sycl: vector is NULL\n");
    return NULL;
  }

  /* Create an empty clone vector */
  N_Vector v = N_VCloneEmpty_Sycl(w);
  if (v == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Sycl: N_VCloneEmpty_Sycl returned NULL\n");
    return NULL;
  }

  /* Allocate content */
  NVEC_SYCL_CONTENT(v)->own_exec           = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->own_helper         = SUNTRUE;
  NVEC_SYCL_CONTENT(v)->stream_exec_policy = NVEC_SYCL_CONTENT(w)->stream_exec_policy->clone();
  NVEC_SYCL_CONTENT(v)->reduce_exec_policy = NVEC_SYCL_CONTENT(w)->reduce_exec_policy->clone();
  NVEC_SYCL_CONTENT(v)->mem_helper         = SUNMemoryHelper_Clone(NVEC_SYCL_MEMHELP(w));

  if (NVEC_SYCL_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Sycl: memory helper is NULL\n");
    N_VDestroy(v);
    return NULL;
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Sycl: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return NULL;
  }

  return v;
}


void N_VDestroy_Sycl(N_Vector v)
{
  N_VectorContent_Sycl vc;
  N_PrivateVectorContent_Sycl vcp;

  if (v == NULL) return;

  /* free ops structure */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }

  /* extract content */
  vc = NVEC_SYCL_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

  /* free private content */
  vcp = (N_PrivateVectorContent_Sycl) vc->priv;
  if (vcp != NULL)
  {
    /* free items in private content */
    FreeReductionBuffer(v);
    FusedBuffer_Free(v);
    free(vcp);
    vc->priv = NULL;
  }

  /* free items in content */
  if (NVEC_SYCL_MEMHELP(v))
  {
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vc->host_data,
                            NVEC_SYCL_QUEUE(v));
    vc->host_data = NULL;
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vc->device_data,
                            NVEC_SYCL_QUEUE(v));
    vc->device_data = NULL;
    if (vc->own_helper) SUNMemoryHelper_Destroy(vc->mem_helper);
    vc->mem_helper = NULL;
  }
  else
  {
    SUNDIALS_DEBUG_PRINT("WARNING in N_VDestroy_Sycl: mem_helper was NULL when trying to dealloc data, this could result in a memory leak\n");
  }

  /* free content struct and vector */
  free(vc);
  free(v);

  return;
}


void N_VSpace_Sycl(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NVEC_SYCL_CONTENT(X)->length;
  *liw = 2;
}


void N_VConst_Sycl(realtype c, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConst_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = c;
           });
}


void N_VLinearSum_Sycl(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  const realtype     *ydata = NVEC_SYCL_DDATAp(y);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearSum_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = (a * xdata[i]) + (b * ydata[i]);
           });
}


void N_VProd_Sycl(N_Vector x, N_Vector y, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  const realtype     *ydata = NVEC_SYCL_DDATAp(y);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VProd_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = xdata[i] * ydata[i];
           });
}


void N_VDiv_Sycl(N_Vector x, N_Vector y, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  const realtype     *ydata = NVEC_SYCL_DDATAp(y);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDiv_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = xdata[i] / ydata[i];
           });
}


void N_VScale_Sycl(realtype c, N_Vector x, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScale_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = c * xdata[i];
           });
}


void N_VAbs_Sycl(N_Vector x, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VAbs_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = abs(xdata[i]);
           });
}


void N_VInv_Sycl(N_Vector x, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInv_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = ONE / xdata[i];
           });
}


void N_VAddConst_Sycl(N_Vector x, realtype b, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VAddConst_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = xdata[i] + b;
           });
}


realtype N_VDotProd_Sycl(N_Vector x, N_Vector y)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(x);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  const realtype     *ydata = NVEC_SYCL_DDATAp(y);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(x);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProd_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProd_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  sum, ::sycl::plus<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    sum += xdata[i] * ydata[i];
                  });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProd_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}


realtype N_VMaxNorm_Sycl(N_Vector x)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(x);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(x);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMaxNorm_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMaxNorm_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *max = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  max, ::sycl::maximum<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    max.combine(abs(xdata[i]));
                  });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMaxNorm_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}


realtype N_VWSqrSumLocal_Sycl(N_Vector x, N_Vector w)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(x);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  const realtype     *wdata = NVEC_SYCL_DDATAp(w);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(x);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumLocal_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumLocal_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  sum, ::sycl::plus<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    sum += xdata[i] * wdata[i] * xdata[i] * wdata[i];
                  });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumLocal_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}


realtype N_VWrmsNorm_Sycl(N_Vector x, N_Vector w)
{
  const sunindextype N   = NVEC_SYCL_LENGTH(x);
  const realtype     sum = N_VWSqrSumLocal_Sycl(x, w);
  return std::sqrt(sum/N);
}


realtype N_VWSqrSumMaskLocal_Sycl(N_Vector x, N_Vector w, N_Vector id)
{
  const sunindextype N       = NVEC_SYCL_LENGTH(x);
  const realtype     *xdata  = NVEC_SYCL_DDATAp(x);
  const realtype     *wdata  = NVEC_SYCL_DDATAp(w);
  const realtype     *iddata = NVEC_SYCL_DDATAp(id);
  ::sycl::queue      *Q      = NVEC_SYCL_QUEUE(x);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumMaskLocal_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumMaskLocal_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  sum, ::sycl::plus<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    if (iddata[i] > ZERO)
                      sum += xdata[i] * wdata[i] * xdata[i] * wdata[i];
                  });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumMaskLocal_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}


realtype N_VWrmsNormMask_Sycl(N_Vector x, N_Vector w, N_Vector id)
{
  const sunindextype N   = NVEC_SYCL_LENGTH(x);
  const realtype     sum = N_VWSqrSumMaskLocal_Sycl(x, w, id);
  return std::sqrt(sum/N);
}


realtype N_VMin_Sycl(N_Vector x)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(x);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(x);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, std::numeric_limits<realtype>::max()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMin_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMin_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *min = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  min, ::sycl::minimum<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    min.combine(xdata[i]);
                  });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMin_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}


realtype N_VWL2Norm_Sycl(N_Vector x, N_Vector w)
{
  return std::sqrt(N_VWSqrSumLocal_Sycl(x, w));
}


realtype N_VL1Norm_Sycl(N_Vector x)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(x);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(x);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VL1Norm_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VL1Norm_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  sum, ::sycl::plus<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    sum += abs(xdata[i]);
                  });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VL1Norm_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(x)[0];
}


void N_VCompare_Sycl(realtype c, N_Vector x, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCompare_Sycl: GetKernelParameters returned nonzero\n");
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
           });
}


booleantype N_VInvTest_Sycl(N_Vector x, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInvTest_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInvTest_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  sum, ::sycl::plus<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    if (xdata[i] == ZERO)
                    {
                      sum += ONE;
                    }
                    else
                    {
                      zdata[i] = ONE / xdata[i];
                    }
                  });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInvTest_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return (NVEC_SYCL_HBUFFERp(x)[0] < HALF);
}


booleantype N_VConstrMask_Sycl(N_Vector c, N_Vector x, N_Vector m)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(x);
  const realtype     *cdata = NVEC_SYCL_DDATAp(c);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  realtype           *mdata = NVEC_SYCL_DDATAp(m);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(x);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(x, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstrMask_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(x, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstrMask_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *sum = NVEC_SYCL_DBUFFERp(x);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  sum, ::sycl::plus<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    bool test =
                      (abs(cdata[i]) > ONEPT5 && cdata[i] * xdata[i] <= ZERO) ||
                      (abs(cdata[i]) > HALF   && cdata[i] * xdata[i] <  ZERO);
                    mdata[i] = test ? ONE : ZERO;
                    sum += mdata[i];
                  });

  if (CopyReductionBufferFromDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstrMask_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return (NVEC_SYCL_HBUFFERp(x)[0] < HALF);
}


realtype N_VMinQuotient_Sycl(N_Vector num, N_Vector denom)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(num);
  const realtype     *ndata = NVEC_SYCL_DDATAp(num);
  const realtype     *ddata = NVEC_SYCL_DDATAp(denom);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(num);
  size_t             nthreads_total, nthreads_per_block;

  if (InitializeReductionBuffer(num, std::numeric_limits<realtype>::max()))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMinQuotient_Sycl: InitializeReductionBuffer returned nonzero\n");
  }

  if (GetKernelParameters(num, SUNTRUE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMinQuotient_Sycl: GetKernelParameters returned nonzero\n");
  }

  /* Shortcut to the reduction buffer */
  realtype *min = NVEC_SYCL_DBUFFERp(num);

  SYCL_FOR_REDUCE(Q, nthreads_total, nthreads_per_block, item,
                  min, ::sycl::minimum<realtype>(),
                  GRID_STRIDE_XLOOP(item, i, N)
                  {
                    if (ddata[i] != ZERO)
                      min.combine(ndata[i] / ddata[i]);
                  });

  if (CopyReductionBufferFromDevice(num))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMinQuotient_Sycl: CopyReductionBufferFromDevice returned nonzero\n");
  }

  return NVEC_SYCL_HBUFFERp(num)[0];
}


/* --------------------------------------------------------------------------
 * fused vector operations
 * -------------------------------------------------------------------------- */


int N_VLinearCombination_Sycl(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(z);
  realtype           *zdata = NVEC_SYCL_DDATAp(z);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(z);
  size_t             nthreads_total, nthreads_per_block;

  /* Fused op workspace shortcuts */
  realtype*  cdata = NULL;
  realtype** xdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(z, nvec, nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Sycl: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(z, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Sycl: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(z, X, nvec, &xdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(z))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Sycl: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  if (GetKernelParameters(z, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombination_Sycl: GetKernelParameters returned nonzero\n");
    return -1;
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             zdata[i] = cdata[0] * xdata[0][i];
             for (int j = 1; j < nvec; j++)
             {
               zdata[i] += cdata[j] * xdata[j][i];
             }
           });

  return 0;
}


int N_VScaleAddMulti_Sycl(int nvec, realtype* c, N_Vector x, N_Vector* Y,
                          N_Vector* Z)
{
  const sunindextype N      = NVEC_SYCL_LENGTH(x);
  const realtype     *xdata = NVEC_SYCL_DDATAp(x);
  ::sycl::queue      *Q     = NVEC_SYCL_QUEUE(x);
  size_t             nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace */
  realtype*  cdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  /* Setup the fused op workspace */
  if (FusedBuffer_Init(x, nvec, 2 * nvec))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Sycl: FusedBuffer_Init returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyRealArray(x, c, nvec, &cdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Sycl: FusedBuffer_CopyRealArray returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Y, nvec, &ydata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyPtrArray1D(x, Z, nvec, &zdata))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  if (FusedBuffer_CopyToDevice(x))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Sycl: FusedBuffer_CopyToDevice returned nonzero\n");
    return -1;
  }

  if (GetKernelParameters(x, SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMulti_Sycl: GetKernelParameters returned nonzero\n");
    return -1;
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             for (int j = 0; j < nvec; j++)
             {
               zdata[j][i] = cdata[j] * xdata[i] + ydata[j][i];
             }
           });

  return 0;
}


/* --------------------------------------------------------------------------
 * vector array operations
 * -------------------------------------------------------------------------- */


int N_VLinearSumVectorArray_Sycl(int nvec,
                                 realtype a, N_Vector* X,
                                 realtype b, N_Vector* Y,
                                 N_Vector* Z)
{
  const sunindextype N  = NVEC_SYCL_LENGTH(Z[0]);
  ::sycl::queue      *Q = NVEC_SYCL_QUEUE(Z[0]);
  size_t             nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace */
  realtype** xdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  /* Setup the fused op workspace */
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

  if (GetKernelParameters(Z[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinaerSumVectorArray_Sycl: GetKernelParameters returned nonzero\n");
    return -1;
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             for (int j = 0; j < nvec; j++)
             {
               zdata[j][i] = a * xdata[j][i] + b * ydata[j][i];
             }
           });

  return 0;
}


int N_VScaleVectorArray_Sycl(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  const sunindextype N  = NVEC_SYCL_LENGTH(Z[0]);
  ::sycl::queue      *Q = NVEC_SYCL_QUEUE(Z[0]);
  size_t             nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace arrays */
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** zdata = NULL;

  /* Setup the fused op workspace */
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

  if (GetKernelParameters(Z[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleVectorArray_Sycl: FusedBuffer_CopyPtrArray1D returned nonzero\n");
    return -1;
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             for (int j = 0; j < nvec; j++)
             {
               zdata[j][i] = cdata[j] * xdata[j][i];
             }
           });

  return 0;
}


int N_VConstVectorArray_Sycl(int nvec, realtype c, N_Vector* Z)
{
  const sunindextype N  = NVEC_SYCL_LENGTH(Z[0]);
  ::sycl::queue      *Q = NVEC_SYCL_QUEUE(Z[0]);
  size_t             nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace arrays */
  realtype** zdata = NULL;

  /* Setup the fused op workspace */
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

  if (GetKernelParameters(Z[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstVectorArray_Sycl: GetKernelParameters returned nonzero\n");
    return -1;
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             for (int j = 0; j < nvec; j++)
             {
               zdata[j][i] = c;
             }
           });

  return 0;
}


int N_VScaleAddMultiVectorArray_Sycl(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  const sunindextype N  = NVEC_SYCL_LENGTH(X[0]);
  ::sycl::queue      *Q = NVEC_SYCL_QUEUE(X[0]);
  size_t             nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace */
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** ydata = NULL;
  realtype** zdata = NULL;

  /* Setup the fused op workspace */
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

  if (GetKernelParameters(X[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VScaleAddMultiVectorArray_Sycl: GetKernelParameters returned nonzero\n");
    return -1;
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             for (int j = 0; j < nvec; j++)
             {
               for (int k = 0; k < nsum; k++)
               {
                 zdata[j * nsum + k][i] =
                   cdata[k] * xdata[j][i] + ydata[j * nsum + k][i];
               }
             }
           });

  return 0;
}


int N_VLinearCombinationVectorArray_Sycl(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  const sunindextype N  = NVEC_SYCL_LENGTH(Z[0]);
  ::sycl::queue      *Q = NVEC_SYCL_QUEUE(Z[0]);
  size_t             nthreads_total, nthreads_per_block;

  /* Shortcuts to the fused op workspace arrays */
  realtype*  cdata = NULL;
  realtype** xdata = NULL;
  realtype** zdata = NULL;

  /* Setup the fused op workspace */
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

  if (GetKernelParameters(Z[0], SUNFALSE, nthreads_total, nthreads_per_block))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VLinearCombinationVectorArray_Sycl: GetKernelParameters returned nonzero\n");
    return -1;
  }

  SYCL_FOR(Q, nthreads_total, nthreads_per_block, item,
           GRID_STRIDE_XLOOP(item, i, N)
           {
             for (int j = 0; j < nvec; j++)
             {
               zdata[j][i] = cdata[0] * xdata[j * nsum][i];
               for (int k = 1; k < nsum; k++)
               {
                 zdata[j][i] += cdata[k] * xdata[j * nsum + k][i];
               }
             }
           });

  return 0;
}


/* --------------------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -------------------------------------------------------------------------- */


int N_VBufSize_Sycl(N_Vector x, sunindextype *size)
{
  if (x == NULL) return -1;
  *size = (sunindextype)NVEC_SYCL_MEMSIZE(x);
  return 0;
}


int N_VBufPack_Sycl(N_Vector x, void *buf)
{
  int copy_fail = 0;

  if (x == NULL || buf == NULL) return -1;

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return -1;

  copy_fail = SUNMemoryHelper_Copy(NVEC_SYCL_MEMHELP(x),
                                   buf_mem,
                                   NVEC_SYCL_CONTENT(x)->device_data,
                                   NVEC_SYCL_MEMSIZE(x),
                                   NVEC_SYCL_QUEUE(x));

  /* synchronize with the host */
  NVEC_SYCL_QUEUE(x)->wait_and_throw();

  SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(x), buf_mem, NVEC_SYCL_QUEUE(x));

  return (copy_fail ? -1 : 0);
}


int N_VBufUnpack_Sycl(N_Vector x, void *buf)
{
  int copy_fail = 0;

  if (x == NULL || buf == NULL) return -1;

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return -1;

  copy_fail = SUNMemoryHelper_Copy(NVEC_SYCL_MEMHELP(x),
                                   NVEC_SYCL_CONTENT(x)->device_data,
                                   buf_mem,
                                   NVEC_SYCL_MEMSIZE(x),
                                   NVEC_SYCL_QUEUE(x));

  /* synchronize with the host */
  NVEC_SYCL_QUEUE(x)->wait_and_throw();

  SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(x), buf_mem, NVEC_SYCL_QUEUE(x));

  return (copy_fail ? -1 : 0);
}


/* --------------------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -------------------------------------------------------------------------- */


int N_VEnableFusedOps_Sycl(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return -1;

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return -1;

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Sycl;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Sycl;
    v->ops->nvdotprodmulti      = NULL;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Sycl;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Sycl;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Sycl;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Sycl;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Sycl;
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
  return 0;
}


int N_VEnableLinearCombination_Sycl(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvlinearcombination = tf ? N_VLinearCombination_Sycl : NULL;
  return 0;
}


int N_VEnableScaleAddMulti_Sycl(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvscaleaddmulti = tf ? N_VScaleAddMulti_Sycl : NULL;
  return 0;
}


int N_VEnableLinearSumVectorArray_Sycl(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvlinearsumvectorarray = tf ? N_VLinearSumVectorArray_Sycl : NULL;
  return 0;
}


int N_VEnableScaleVectorArray_Sycl(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvscalevectorarray = tf ? N_VScaleVectorArray_Sycl : NULL;
  return 0;
}


int N_VEnableConstVectorArray_Sycl(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvconstvectorarray = tf ? N_VConstVectorArray_Sycl : NULL;
  return 0;
}


int N_VEnableScaleAddMultiVectorArray_Sycl(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvscaleaddmultivectorarray = tf ?
    N_VScaleAddMultiVectorArray_Sycl : NULL;
  return 0;
}


int N_VEnableLinearCombinationVectorArray_Sycl(N_Vector v, booleantype tf)
{
  if (v == NULL) return -1;
  if (v->ops == NULL) return -1;
  v->ops->nvlinearcombinationvectorarray = tf ?
    N_VLinearCombinationVectorArray_Sycl : NULL;
  return 0;
}


/* --------------------------------------------------------------------------
 * Private utility functions
 * -------------------------------------------------------------------------- */


static int AllocateData(N_Vector v)
{
  int                         alloc_fail = 0;
  N_VectorContent_Sycl        vc         = NVEC_SYCL_CONTENT(v);
  N_PrivateVectorContent_Sycl vcp        = NVEC_SYCL_PRIVATE(v);

  if (N_VGetLength_Sycl(v) == 0) return 0;

  if (vcp->use_managed_mem)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v), &(vc->device_data),
                                       NVEC_SYCL_MEMSIZE(v), SUNMEMTYPE_UVM,
                                       NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed for SUNMEMTYPE_UVM\n");
    }
    vc->host_data = SUNMemoryHelper_Alias(vc->device_data);
  }
  else
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v), &(vc->host_data),
                                       NVEC_SYCL_MEMSIZE(v), SUNMEMTYPE_HOST,
                                       NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
    }

    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v), &(vc->device_data),
                                       NVEC_SYCL_MEMSIZE(v), SUNMEMTYPE_DEVICE,
                                       NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  return (alloc_fail ? -1 : 0);
}


/* Allocate and initializes the internal memory used for reductions */
static int InitializeReductionBuffer(N_Vector v, const realtype value, size_t n)
{
  int         alloc_fail = 0;
  int         copy_fail  = 0;
  booleantype alloc_mem  = SUNFALSE;
  size_t      bytes      = n * sizeof(realtype);

  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Wrap the initial value as SUNMemory object */
  SUNMemory value_mem = SUNMemoryHelper_Wrap((void*) &value, SUNMEMTYPE_HOST);

  /* check if the existing reduction memory is not large enough */
  if (vcp->reduce_buffer_bytes < bytes)
  {
    FreeReductionBuffer(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
    /* allocate pinned memory on the host */
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                       &(vcp->reduce_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED, NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("WARNING in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      /* if pinned alloc failed, allocate plain host memory */
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                         &(vcp->reduce_buffer_host), bytes,
                                         SUNMEMTYPE_HOST, NVEC_SYCL_QUEUE(v));
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
      }
    }

    /* allocate device memory */
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                       &(vcp->reduce_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE, NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
    }
  }

  if (!alloc_fail)
  {
    /* store the size of the reduction memory */
    vcp->reduce_buffer_bytes = bytes;

    /* initialize the memory with the value */
    copy_fail = SUNMemoryHelper_CopyAsync(NVEC_SYCL_MEMHELP(v),
                                          vcp->reduce_buffer_dev, value_mem,
                                          bytes, (void*) NVEC_SYCL_QUEUE(v));

    /* wait for copy to finish (possible bug in minimum reduction object) */
    NVEC_SYCL_QUEUE(v)->wait_and_throw();

    if (copy_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_CopyAsync failed\n");
    }
  }

  /* deallocate the wrapper */
  SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), value_mem, NVEC_SYCL_QUEUE(v));

  return ((alloc_fail || copy_fail) ? -1 : 0);
}


/* Free the reduction memory */
static void FreeReductionBuffer(N_Vector v)
{
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  if (vcp == NULL) return;

  /* free device mem */
  if (vcp->reduce_buffer_dev != NULL)
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vcp->reduce_buffer_dev,
                            NVEC_SYCL_QUEUE(v));
  vcp->reduce_buffer_dev = NULL;

  /* free host mem */
  if (vcp->reduce_buffer_host != NULL)
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v), vcp->reduce_buffer_host,
                            NVEC_SYCL_QUEUE(v));
  vcp->reduce_buffer_host = NULL;

  /* reset allocated memory size */
  vcp->reduce_buffer_bytes = 0;
}


/* Copy the reduction memory from the device to the host. */
static int CopyReductionBufferFromDevice(N_Vector v, size_t n)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_SYCL_MEMHELP(v),
                                        NVEC_SYCL_PRIVATE(v)->reduce_buffer_host,
                                        NVEC_SYCL_PRIVATE(v)->reduce_buffer_dev,
                                        n * sizeof(realtype),
                                        (void*) NVEC_SYCL_QUEUE(v));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in CopyReductionBufferFromDevice: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* synchronize with respect to the host */
  NVEC_SYCL_QUEUE(v)->wait_and_throw();

  return (copy_fail ? -1 : 0);
}


static int FusedBuffer_Init(N_Vector v, int nreal, int nptr)
{
  int         alloc_fail = 0;
  booleantype alloc_mem  = SUNFALSE;
  size_t      bytes      = nreal * sizeof(realtype) + nptr * sizeof(realtype*);

  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Check if the existing memory is not large enough */
  if (vcp->fused_buffer_bytes < bytes)
  {
    FusedBuffer_Free(v);
    alloc_mem = SUNTRUE;
  }

  if (alloc_mem)
  {
    /* allocate pinned memory on the host */
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                       &(vcp->fused_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED, NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("WARNING in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      /* if pinned alloc failed, allocate plain host memory */
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                         &(vcp->fused_buffer_host), bytes,
                                         SUNMEMTYPE_HOST, NVEC_SYCL_QUEUE(v));
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
        return -1;
      }
    }

    /* allocate device memory */
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_SYCL_MEMHELP(v),
                                       &(vcp->fused_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE, NVEC_SYCL_QUEUE(v));
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_Init: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_DEVICE\n");
      return -1;
    }

    /* Store the size of the fused op buffer */
    vcp->fused_buffer_bytes = bytes;
  }

  /* Reset the buffer offset */
  vcp->fused_buffer_offset = 0;

  return 0;
}


static int FusedBuffer_CopyRealArray(N_Vector v, realtype *rdata, int nval,
                                     realtype **shortcut)
{
  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Check buffer space and fill the host buffer */
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

  /* Set shortcut to the device buffer and update offset*/
  *shortcut = (realtype*) ((char*)(vcp->fused_buffer_dev->ptr) +
                           vcp->fused_buffer_offset);

  vcp->fused_buffer_offset += nval * sizeof(realtype);

  return 0;
}


static int FusedBuffer_CopyPtrArray1D(N_Vector v, N_Vector *X, int nvec,
                                      realtype ***shortcut)
{
  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Check buffer space and fill the host buffer */
  if (vcp->fused_buffer_offset >= vcp->fused_buffer_bytes)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyPtrArray1D: Buffer offset is exceedes the buffer size\n");    return -1;
    return -1;
  }

  realtype** h_buffer = (realtype**) ((char*)(vcp->fused_buffer_host->ptr) +
                                      vcp->fused_buffer_offset);

  for (int j = 0; j < nvec; j++)
  {
    h_buffer[j] = NVEC_SYCL_DDATAp(X[j]);
  }

  /* Set shortcut to the device buffer and update offset*/
  *shortcut = (realtype**) ((char*)(vcp->fused_buffer_dev->ptr) +
                            vcp->fused_buffer_offset);

  vcp->fused_buffer_offset += nvec * sizeof(realtype*);

  return 0;
}


static int FusedBuffer_CopyPtrArray2D(N_Vector v, N_Vector **X, int nvec,
                                      int nsum, realtype ***shortcut)
{
  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Check buffer space and fill the host buffer */
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
      h_buffer[j * nsum + k] = NVEC_SYCL_DDATAp(X[k][j]);
    }
  }

  /* Set shortcut to the device buffer and update offset*/
  *shortcut = (realtype**) ((char*)(vcp->fused_buffer_dev->ptr) +
                            vcp->fused_buffer_offset);

  /* Update the offset */
  vcp->fused_buffer_offset += nvec * nsum * sizeof(realtype*);

  return 0;
}


static int FusedBuffer_CopyToDevice(N_Vector v)
{
  /* Get the vector private memory structure */
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  /* Copy the fused buffer to the device */
  int copy_fail = SUNMemoryHelper_CopyAsync(NVEC_SYCL_MEMHELP(v),
                                            vcp->fused_buffer_dev,
                                            vcp->fused_buffer_host,
                                            vcp->fused_buffer_offset,
                                            (void*) NVEC_SYCL_QUEUE(v));
  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in FusedBuffer_CopyToDevice: SUNMemoryHelper_CopyAsync failed\n");
    return -1;
  }

  return 0;
}


static int FusedBuffer_Free(N_Vector v)
{
  N_PrivateVectorContent_Sycl vcp = NVEC_SYCL_PRIVATE(v);

  if (vcp == NULL) return 0;

  if (vcp->fused_buffer_host)
  {
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v),
                            vcp->fused_buffer_host, NVEC_SYCL_QUEUE(v));
    vcp->fused_buffer_host = NULL;
  }

  if (vcp->fused_buffer_dev)
  {
    SUNMemoryHelper_Dealloc(NVEC_SYCL_MEMHELP(v),
                            vcp->fused_buffer_dev, NVEC_SYCL_QUEUE(v));
    vcp->fused_buffer_dev = NULL;
  }

  vcp->fused_buffer_bytes  = 0;
  vcp->fused_buffer_offset = 0;

  return 0;
}


/* Get the kernel launch parameters based on the kernel type (reduction or not),
 * using the appropriate kernel execution policy. */
static int GetKernelParameters(N_Vector v, booleantype reduction,
                               size_t& nthreads_total,
                               size_t& nthreads_per_block)
{
  /* Get the execution policy */
  SUNSyclExecPolicy* exec_policy = NULL;
  exec_policy = reduction ?
    NVEC_SYCL_CONTENT(v)->reduce_exec_policy :
    NVEC_SYCL_CONTENT(v)->stream_exec_policy;

  if (exec_policy == NULL)
  {
    SUNDIALS_DEBUG_ERROR("The execution policy is NULL\n");
    return -1;
  }

  /* Get the number of threads per block and total number threads */
  nthreads_per_block = exec_policy->blockSize();
  nthreads_total     = nthreads_per_block *
                       exec_policy->gridSize(NVEC_SYCL_LENGTH(v));

  if (nthreads_per_block == 0)
  {
    SUNDIALS_DEBUG_ERROR("The number of threads per block must be > 0\n");
    return -1;
  }

  if (nthreads_total == 0)
  {
    SUNDIALS_DEBUG_ERROR("the total number of threads must be > 0\n");
    return -1;
  }

  return 0;
}


} /* extern "C" */
