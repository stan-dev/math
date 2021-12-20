/* -----------------------------------------------------------------
 * Programmer(s): Daniel McGreer, and Cody J. Balos @ LLNL
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
 * This is the implementation file for a HIP implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>

#include <nvector/nvector_hip.h>
#include "VectorArrayKernels.hip.hpp"
#include "VectorKernels.hip.hpp"
#include "sundials_hip.h"
#include "sundials_debug.h"

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)

extern "C" {

using namespace sundials;
using namespace sundials::hip;
using namespace sundials::nvector_hip;

/*
 * Macro definitions
 */

#define NVEC_HIP_CONTENT(x)  ((N_VectorContent_Hip)(x->content))
#define NVEC_HIP_PRIVATE(x)  ((N_PrivateVectorContent_Hip)(NVEC_HIP_CONTENT(x)->priv))
#define NVEC_HIP_MEMSIZE(x)  (NVEC_HIP_CONTENT(x)->length * sizeof(realtype))
#define NVEC_HIP_MEMHELP(x)  (NVEC_HIP_CONTENT(x)->mem_helper)
#define NVEC_HIP_HDATAp(x)   ((realtype*) NVEC_HIP_CONTENT(x)->host_data->ptr)
#define NVEC_HIP_DDATAp(x)   ((realtype*) NVEC_HIP_CONTENT(x)->device_data->ptr)
#define NVEC_HIP_HBUFFERp(x) ((realtype*) NVEC_HIP_PRIVATE(x)->reduce_buffer_host->ptr)
#define NVEC_HIP_DBUFFERp(x) ((realtype*) NVEC_HIP_PRIVATE(x)->reduce_buffer_dev->ptr)
#define NVEC_HIP_STREAM(x)   (NVEC_HIP_CONTENT(x)->stream_exec_policy->stream())


/*
 * Private structure definition
 */

struct _N_PrivateVectorContent_Hip
{
  booleantype     use_managed_mem;               /* indicates if the data pointers and buffer pointers are managed memory */
  size_t          reduce_buffer_allocated_bytes; /* current size of the reduction buffer */
  SUNMemory       reduce_buffer_dev;             /* device buffer used for reductions */
  SUNMemory       reduce_buffer_host;            /* host buffer used for reductions */
};

typedef struct _N_PrivateVectorContent_Hip *N_PrivateVectorContent_Hip;

/*
 * Private function definitions
 */

static int AllocateData(N_Vector v);
static int InitializeReductionBuffer(N_Vector v, const realtype value, size_t n = 1);
static void FreeReductionBuffer(N_Vector v);
static int CopyReductionBufferFromDevice(N_Vector v, size_t n = 1);
static int GetKernelParameters(N_Vector v, booleantype reduction, size_t& grid, size_t& block,
                               size_t& shMemSize, hipStream_t& stream, size_t n = 0);
static void PostKernelLaunch();

/*
 * Defaults
 */

static ThreadDirectExecPolicy NVEC_HIP_DEFAULT_STREAM_POLICY(512);
static BlockReduceExecPolicy NVEC_HIP_DEFAULT_REDUCE_POLICY(512);


/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Hip(N_Vector v)
{
  return SUNDIALS_NVEC_HIP;
}

N_Vector N_VNewEmpty_Hip(SUNContext sunctx)
{
  N_Vector v;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid           = N_VGetVectorID_Hip;
  v->ops->nvclone                 = N_VClone_Hip;
  v->ops->nvcloneempty            = N_VCloneEmpty_Hip;
  v->ops->nvdestroy               = N_VDestroy_Hip;
  v->ops->nvspace                 = N_VSpace_Hip;
  v->ops->nvgetlength             = N_VGetLength_Hip;
  v->ops->nvgetarraypointer       = N_VGetHostArrayPointer_Hip;
  v->ops->nvgetdevicearraypointer = N_VGetDeviceArrayPointer_Hip;
  v->ops->nvsetarraypointer       = N_VSetHostArrayPointer_Hip;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Hip;
  v->ops->nvconst        = N_VConst_Hip;
  v->ops->nvprod         = N_VProd_Hip;
  v->ops->nvdiv          = N_VDiv_Hip;
  v->ops->nvscale        = N_VScale_Hip;
  v->ops->nvabs          = N_VAbs_Hip;
  v->ops->nvinv          = N_VInv_Hip;
  v->ops->nvaddconst     = N_VAddConst_Hip;
  v->ops->nvdotprod      = N_VDotProd_Hip;
  v->ops->nvmaxnorm      = N_VMaxNorm_Hip;
  v->ops->nvmin          = N_VMin_Hip;
  v->ops->nvl1norm       = N_VL1Norm_Hip;
  v->ops->nvinvtest      = N_VInvTest_Hip;
  v->ops->nvconstrmask   = N_VConstrMask_Hip;
  v->ops->nvminquotient  = N_VMinQuotient_Hip;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Hip;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Hip;
  v->ops->nvwl2norm      = N_VWL2Norm_Hip;
  v->ops->nvcompare      = N_VCompare_Hip;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProd_Hip;
  v->ops->nvmaxnormlocal     = N_VMaxNorm_Hip;
  v->ops->nvminlocal         = N_VMin_Hip;
  v->ops->nvl1normlocal      = N_VL1Norm_Hip;
  v->ops->nvinvtestlocal     = N_VInvTest_Hip;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Hip;
  v->ops->nvminquotientlocal = N_VMinQuotient_Hip;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Hip;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Hip;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal = N_VDotProdMulti_Hip;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_Hip;
  v->ops->nvbufpack   = N_VBufPack_Hip;
  v->ops->nvbufunpack = N_VBufUnpack_Hip;

  /* print operation for debugging */
  v->ops->nvprint     = N_VPrint_Hip;
  v->ops->nvprintfile = N_VPrintFile_Hip;

  /* Create content */

  v->content = (N_VectorContent_Hip) malloc(sizeof(_N_VectorContent_Hip));
  if (v->content == NULL)
  {
    N_VDestroy(v);
    return(NULL);
  }

  NVEC_HIP_CONTENT(v)->priv = malloc(sizeof(_N_PrivateVectorContent_Hip));
  if (NVEC_HIP_CONTENT(v)->priv == NULL)
  {
    N_VDestroy(v);
    return(NULL);
  }

  NVEC_HIP_CONTENT(v)->length                        = 0;
  NVEC_HIP_CONTENT(v)->host_data                     = NULL;
  NVEC_HIP_CONTENT(v)->device_data                   = NULL;
  NVEC_HIP_CONTENT(v)->stream_exec_policy            = NULL;
  NVEC_HIP_CONTENT(v)->reduce_exec_policy            = NULL;
  NVEC_HIP_CONTENT(v)->mem_helper                    = NULL;
  NVEC_HIP_CONTENT(v)->own_helper                    = SUNFALSE;
  NVEC_HIP_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  return(v);
}

N_Vector N_VNew_Hip(sunindextype length, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Hip(sunctx);
  if (v == NULL) return(NULL);

  NVEC_HIP_CONTENT(v)->length                        = length;
  NVEC_HIP_CONTENT(v)->host_data                     = NULL;
  NVEC_HIP_CONTENT(v)->device_data                   = NULL;
  NVEC_HIP_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Hip(sunctx);
  NVEC_HIP_CONTENT(v)->stream_exec_policy            = NVEC_HIP_DEFAULT_STREAM_POLICY.clone();
  NVEC_HIP_CONTENT(v)->reduce_exec_policy            = NVEC_HIP_DEFAULT_REDUCE_POLICY.clone();
  NVEC_HIP_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_HIP_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  if (NVEC_HIP_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Hip: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNew_Hip: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VNewWithMemHelp_Hip(sunindextype length, booleantype use_managed_mem, SUNMemoryHelper helper, SUNContext sunctx)
{
  N_Vector v;

  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Hip: helper is NULL\n");
    return(NULL);
  }

  if (!SUNMemoryHelper_ImplementsRequiredOps(helper))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Hip: helper doesn't implement all required ops\n");
    return(NULL);
  }

  v = NULL;
  v = N_VNewEmpty_Hip(sunctx);
  if (v == NULL) return(NULL);

  NVEC_HIP_CONTENT(v)->length                        = length;
  NVEC_HIP_CONTENT(v)->host_data                     = NULL;
  NVEC_HIP_CONTENT(v)->device_data                   = NULL;
  NVEC_HIP_CONTENT(v)->mem_helper                    = helper;
  NVEC_HIP_CONTENT(v)->stream_exec_policy            = NVEC_HIP_DEFAULT_STREAM_POLICY.clone();
  NVEC_HIP_CONTENT(v)->reduce_exec_policy            = NVEC_HIP_DEFAULT_REDUCE_POLICY.clone();
  NVEC_HIP_CONTENT(v)->own_helper                    = SUNFALSE;
  NVEC_HIP_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->use_managed_mem               = use_managed_mem;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewWithMemHelp_Hip: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VNewManaged_Hip(sunindextype length, SUNContext sunctx)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Hip(sunctx);
  if (v == NULL) return(NULL);

  NVEC_HIP_CONTENT(v)->length                        = length;
  NVEC_HIP_CONTENT(v)->host_data                     = NULL;
  NVEC_HIP_CONTENT(v)->device_data                   = NULL;
  NVEC_HIP_CONTENT(v)->stream_exec_policy            = NVEC_HIP_DEFAULT_STREAM_POLICY.clone();
  NVEC_HIP_CONTENT(v)->reduce_exec_policy            = NVEC_HIP_DEFAULT_REDUCE_POLICY.clone();
  NVEC_HIP_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Hip(sunctx);
  NVEC_HIP_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_HIP_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->use_managed_mem               = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  if (NVEC_HIP_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Hip: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VNewManaged_Hip: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VMake_Hip(sunindextype length, realtype *h_vdata, realtype *d_vdata, SUNContext sunctx)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Hip(sunctx);
  if (v == NULL) return(NULL);

  NVEC_HIP_CONTENT(v)->length                        = length;
  NVEC_HIP_CONTENT(v)->host_data                     = SUNMemoryHelper_Wrap(h_vdata, SUNMEMTYPE_HOST);
  NVEC_HIP_CONTENT(v)->device_data                   = SUNMemoryHelper_Wrap(d_vdata, SUNMEMTYPE_DEVICE);
  NVEC_HIP_CONTENT(v)->stream_exec_policy            = NVEC_HIP_DEFAULT_STREAM_POLICY.clone();
  NVEC_HIP_CONTENT(v)->reduce_exec_policy            = NVEC_HIP_DEFAULT_REDUCE_POLICY.clone();
  NVEC_HIP_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Hip(sunctx);
  NVEC_HIP_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_HIP_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->use_managed_mem               = SUNFALSE;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  if (NVEC_HIP_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Hip: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (NVEC_HIP_CONTENT(v)->device_data == NULL ||
      NVEC_HIP_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMake_Hip: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

N_Vector N_VMakeManaged_Hip(sunindextype length, realtype *vdata, SUNContext sunctx)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Hip(sunctx);
  if (v == NULL) return(NULL);

  NVEC_HIP_CONTENT(v)->length                        = length;
  NVEC_HIP_CONTENT(v)->host_data                     = SUNMemoryHelper_Wrap(vdata, SUNMEMTYPE_UVM);
  NVEC_HIP_CONTENT(v)->device_data                   = SUNMemoryHelper_Alias(NVEC_HIP_CONTENT(v)->host_data);
  NVEC_HIP_CONTENT(v)->stream_exec_policy            = NVEC_HIP_DEFAULT_STREAM_POLICY.clone();
  NVEC_HIP_CONTENT(v)->reduce_exec_policy            = NVEC_HIP_DEFAULT_REDUCE_POLICY.clone();
  NVEC_HIP_CONTENT(v)->mem_helper                    = SUNMemoryHelper_Hip(sunctx);
  NVEC_HIP_CONTENT(v)->own_helper                    = SUNTRUE;
  NVEC_HIP_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->use_managed_mem               = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  if (NVEC_HIP_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Hip: memory helper is NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (NVEC_HIP_CONTENT(v)->device_data == NULL ||
      NVEC_HIP_CONTENT(v)->host_data == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMakeManaged_Hip: SUNMemoryHelper_Wrap returned NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

void N_VSetHostArrayPointer_Hip(realtype* h_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Hip(v))
  {
    if (NVEC_HIP_CONTENT(v)->host_data)
    {
      NVEC_HIP_CONTENT(v)->host_data->ptr = (void*) h_vdata;
      NVEC_HIP_CONTENT(v)->device_data->ptr = (void*) h_vdata;
    }
    else
    {
      NVEC_HIP_CONTENT(v)->host_data = SUNMemoryHelper_Wrap((void*) h_vdata, SUNMEMTYPE_UVM);
      NVEC_HIP_CONTENT(v)->device_data = SUNMemoryHelper_Alias(NVEC_HIP_CONTENT(v)->host_data);
    }
  }
  else
  {
    if (NVEC_HIP_CONTENT(v)->host_data)
    {
      NVEC_HIP_CONTENT(v)->host_data->ptr = (void*) h_vdata;
    }
    else
    {
      NVEC_HIP_CONTENT(v)->host_data = SUNMemoryHelper_Wrap((void*) h_vdata, SUNMEMTYPE_HOST);
    }
  }
}

void N_VSetDeviceArrayPointer_Hip(realtype* d_vdata, N_Vector v)
{
  if (N_VIsManagedMemory_Hip(v))
  {
    if (NVEC_HIP_CONTENT(v)->device_data)
    {
      NVEC_HIP_CONTENT(v)->device_data->ptr = (void*) d_vdata;
      NVEC_HIP_CONTENT(v)->host_data->ptr = (void*) d_vdata;
    }
    else
    {
      NVEC_HIP_CONTENT(v)->device_data = SUNMemoryHelper_Wrap((void*) d_vdata, SUNMEMTYPE_UVM);
      NVEC_HIP_CONTENT(v)->host_data = SUNMemoryHelper_Alias(NVEC_HIP_CONTENT(v)->device_data);
    }
  }
  else
  {
    if (NVEC_HIP_CONTENT(v)->device_data)
    {
      NVEC_HIP_CONTENT(v)->device_data->ptr = (void*) d_vdata;
    }
    else
    {
      NVEC_HIP_CONTENT(v)->device_data = SUNMemoryHelper_Wrap((void*) d_vdata, SUNMEMTYPE_DEVICE);
    }
  }
}

booleantype N_VIsManagedMemory_Hip(N_Vector x)
{
  return NVEC_HIP_PRIVATE(x)->use_managed_mem;
}

int N_VSetKernelExecPolicy_Hip(N_Vector x,
                               SUNHipExecPolicy* stream_exec_policy,
                               SUNHipExecPolicy* reduce_exec_policy)
{
  if (x == NULL || stream_exec_policy == NULL || reduce_exec_policy == NULL)
    return(-1);

  if (NVEC_HIP_CONTENT(x)->own_exec)
  {
    delete NVEC_HIP_CONTENT(x)->stream_exec_policy;
    delete NVEC_HIP_CONTENT(x)->reduce_exec_policy;
  }

  NVEC_HIP_CONTENT(x)->stream_exec_policy = stream_exec_policy;
  NVEC_HIP_CONTENT(x)->reduce_exec_policy = reduce_exec_policy;
  NVEC_HIP_CONTENT(x)->own_exec = SUNFALSE;

  return(0);
}

void N_VCopyToDevice_Hip(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_HIP_MEMHELP(x),
                                        NVEC_HIP_CONTENT(x)->device_data,
                                        NVEC_HIP_CONTENT(x)->host_data,
                                        NVEC_HIP_MEMSIZE(x),
                                        (void*) NVEC_HIP_STREAM(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyToDevice_Hip: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  SUNDIALS_HIP_VERIFY(hipStreamSynchronize(*NVEC_HIP_STREAM(x)));
}

void N_VCopyFromDevice_Hip(N_Vector x)
{
  int copy_fail;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_HIP_MEMHELP(x),
                                        NVEC_HIP_CONTENT(x)->host_data,
                                        NVEC_HIP_CONTENT(x)->device_data,
                                        NVEC_HIP_MEMSIZE(x),
                                        (void*) NVEC_HIP_STREAM(x));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VCopyFromDevice_Hip: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  SUNDIALS_HIP_VERIFY(hipStreamSynchronize(*NVEC_HIP_STREAM(x)));
}

void N_VPrint_Hip(N_Vector x)
{
  N_VPrintFile_Hip(x, stdout);
}

void N_VPrintFile_Hip(N_Vector x, FILE *outfile)
{
  sunindextype i;

#ifdef SUNDIALS_DEBUG_PRINTVEC
  N_VCopyFromDevice_Hip(x);
#endif

  for (i = 0; i < NVEC_HIP_CONTENT(x)->length; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Le\n", NVEC_HIP_HDATAp(x)[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16e\n", NVEC_HIP_HDATAp(x)[i]);
#else
    fprintf(outfile, "%11.8e\n", NVEC_HIP_HDATAp(x)[i]);
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

N_Vector N_VCloneEmpty_Hip(N_Vector w)
{
  N_Vector v;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty_Hip(w->sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Set content */
  NVEC_HIP_CONTENT(v)->length                        = NVEC_HIP_CONTENT(w)->length;
  NVEC_HIP_CONTENT(v)->host_data                     = NULL;
  NVEC_HIP_CONTENT(v)->device_data                   = NULL;
  NVEC_HIP_CONTENT(v)->mem_helper                    = NULL;
  NVEC_HIP_CONTENT(v)->own_exec                      = SUNTRUE;
  NVEC_HIP_PRIVATE(v)->use_managed_mem               = NVEC_HIP_PRIVATE(w)->use_managed_mem;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_dev             = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_host            = NULL;
  NVEC_HIP_PRIVATE(v)->reduce_buffer_allocated_bytes = 0;

  return(v);
}

N_Vector N_VClone_Hip(N_Vector w)
{
  N_Vector v;

  v = NULL;
  v = N_VCloneEmpty_Hip(w);
  if (v == NULL) return(NULL);

  NVEC_HIP_MEMHELP(v) = SUNMemoryHelper_Clone(NVEC_HIP_MEMHELP(w));
  NVEC_HIP_CONTENT(v)->own_helper = SUNTRUE;
  NVEC_HIP_CONTENT(v)->stream_exec_policy = NVEC_HIP_CONTENT(w)->stream_exec_policy->clone();
  NVEC_HIP_CONTENT(v)->reduce_exec_policy = NVEC_HIP_CONTENT(w)->reduce_exec_policy->clone();

  if (NVEC_HIP_MEMHELP(v) == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Hip: SUNMemoryHelper_Clone returned NULL\n");
    N_VDestroy(v);
    return(NULL);
  }

  if (AllocateData(v))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VClone_Hip: AllocateData returned nonzero\n");
    N_VDestroy(v);
    return(NULL);
  }

  return(v);
}

void N_VDestroy_Hip(N_Vector v)
{
  N_VectorContent_Hip vc;
  N_PrivateVectorContent_Hip vcp;

  if (v == NULL) return;

  /* free ops structure */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }

  /* extract content */
  vc = NVEC_HIP_CONTENT(v);
  if (vc == NULL)
  {
    free(v);
    v = NULL;
    return;
  }

  /* free private content */
  vcp = (N_PrivateVectorContent_Hip) vc->priv;
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

  if (NVEC_HIP_MEMHELP(v))
  {
    SUNMemoryHelper_Dealloc(NVEC_HIP_MEMHELP(v), vc->host_data, nullptr);
    vc->host_data = NULL;
    SUNMemoryHelper_Dealloc(NVEC_HIP_MEMHELP(v), vc->device_data, nullptr);
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

void N_VSpace_Hip(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  *lrw = NVEC_HIP_CONTENT(X)->length;
  *liw = 2;
}

void N_VConst_Hip(realtype a, N_Vector X)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(setConstKernel, grid, block, shMemSize, stream,
    a,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VLinearSum_Hip(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(linearSumKernel, grid, block, shMemSize, stream,
    a,
    NVEC_HIP_DDATAp(X),
    b,
    NVEC_HIP_DDATAp(Y),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VProd_Hip(N_Vector X, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(prodKernel, grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Y),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VDiv_Hip(N_Vector X, N_Vector Y, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(divKernel, grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Y),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VScale_Hip(realtype a, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(scaleKernel, grid, block, shMemSize, stream,
    a,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VAbs_Hip(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(absKernel, grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VInv_Hip(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(invKernel, grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

void N_VAddConst_Hip(N_Vector X, realtype b, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(addConstKernel, grid, block, shMemSize, stream,
    b,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

realtype N_VDotProd_Hip(N_Vector X, N_Vector Y)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VDotProd_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(dotProdKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Y),
    NVEC_HIP_DBUFFERp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_HIP_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VMaxNorm_Hip(N_Vector X)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMaxNorm_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(maxNormKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DBUFFERp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Finish reduction on CPU if there are less than two blocks of data left.
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_HIP_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWSqrSumLocal_Hip(N_Vector X, N_Vector W)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumLocal_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(wL2NormSquareKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(W),
    NVEC_HIP_DBUFFERp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_HIP_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWrmsNorm_Hip(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Hip(X, W);
  return std::sqrt(sum/NVEC_HIP_CONTENT(X)->length);
}

realtype N_VWSqrSumMaskLocal_Hip(N_Vector X, N_Vector W, N_Vector Id)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VWSqrSumMaskLocal_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(wL2NormSquareMaskKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(W),
    NVEC_HIP_DDATAp(Id),
    NVEC_HIP_DBUFFERp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_HIP_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWrmsNormMask_Hip(N_Vector X, N_Vector W, N_Vector Id)
{
  const realtype sum = N_VWSqrSumMaskLocal_Hip(X, W, Id);
  return std::sqrt(sum/NVEC_HIP_CONTENT(X)->length);
}

realtype N_VMin_Hip(N_Vector X)
{
  const realtype maxVal = std::numeric_limits<realtype>::max();

  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(X, maxVal))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMin_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(findMinKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    maxVal,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DBUFFERp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_HIP_HBUFFERp(X)[0];

  return gpu_result;
}

realtype N_VWL2Norm_Hip(N_Vector X, N_Vector W)
{
  const realtype sum = N_VWSqrSumLocal_Hip(X, W);
  return std::sqrt(sum);
}

realtype N_VL1Norm_Hip(N_Vector X)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VL1Norm_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(L1NormKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DBUFFERp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_HIP_HBUFFERp(X)[0];

  return gpu_result;
}

void N_VCompare_Hip(realtype c, N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  GetKernelParameters(X, false, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(compareKernel, grid, block, shMemSize, stream,
    c,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();
}

booleantype N_VInvTest_Hip(N_Vector X, N_Vector Z)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VInvTest_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(invTestKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_DBUFFERp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_HIP_HBUFFERp(X)[0];

  return (gpu_result < HALF);
}

booleantype N_VConstrMask_Hip(N_Vector C, N_Vector X, N_Vector M)
{
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(X, ZERO))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VConstrMask_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(X, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(constrMaskKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    NVEC_HIP_DDATAp(C),
    NVEC_HIP_DDATAp(X),
    NVEC_HIP_DDATAp(M),
    NVEC_HIP_DBUFFERp(X),
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(X);
  realtype gpu_result = NVEC_HIP_HBUFFERp(X)[0];

  return (gpu_result < HALF);
}

realtype N_VMinQuotient_Hip(N_Vector num, N_Vector denom)
{
  // Starting value for min reduction
  const realtype maxVal = std::numeric_limits<realtype>::max();
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (InitializeReductionBuffer(num, maxVal))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in N_VMinQuotient_Hip: InitializeReductionBuffer returned nonzero\n");
  }

  GetKernelParameters(num, true, grid, block, shMemSize, stream);
  hipLaunchKernelGGL(HIP_KERNEL_NAME(minQuotientKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    maxVal,
    NVEC_HIP_DDATAp(num),
    NVEC_HIP_DDATAp(denom),
    NVEC_HIP_DBUFFERp(num),
    NVEC_HIP_CONTENT(num)->length
  );
  PostKernelLaunch();

  // Get result from the GPU
  CopyReductionBufferFromDevice(num);
  realtype gpu_result = NVEC_HIP_HBUFFERp(num)[0];

  return gpu_result;
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination_Hip(int nvec, realtype* c, N_Vector* X, N_Vector Z)
{
  hipError_t err;

  // Copy c array to device
  realtype* d_c;
  err = hipMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_c, c, nvec*sizeof(realtype), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_HIP_DDATAp(X[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = hipMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters and launch
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(X[0], false, grid, block, shMemSize, stream)) return(-1);
  hipLaunchKernelGGL(linearCombinationKernel, grid, block, shMemSize, stream,
    nvec,
    d_c,
    d_Xd,
    NVEC_HIP_DDATAp(Z),
    NVEC_HIP_CONTENT(Z)->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;

  // Free device arrays
  err = hipFree(d_c);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Xd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}

int N_VScaleAddMulti_Hip(int nvec, realtype* c, N_Vector X, N_Vector* Y,
                          N_Vector* Z)
{
  hipError_t err;

  // Copy c array to device
  realtype* d_c;
  err = hipMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_c, c, nvec*sizeof(realtype), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_HIP_DDATAp(Y[i]);

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_HIP_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Yd;
  err = hipMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = hipMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream)) return(-1);
  hipLaunchKernelGGL(scaleAddMultiKernel, grid, block, shMemSize, stream,
    nvec,
    d_c,
    NVEC_HIP_DDATAp(X),
    d_Yd,
    d_Zd,
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = hipFree(d_c);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Yd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Zd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}

int N_VDotProdMulti_Hip(int nvec, N_Vector X, N_Vector* Y, realtype* dots)
{
  hipError_t err;

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_HIP_DDATAp(Y[i]);

  // Copy array of device pointers to device from host
  realtype** d_Yd;
  err = hipMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(X, false, grid, block, shMemSize, stream)) return(-1);
  grid = nvec;

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = hipMalloc((void**) &d_buff, grid*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemsetAsync(d_buff, 0, grid*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  hipLaunchKernelGGL(HIP_KERNEL_NAME(dotProdMultiKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    nvec,
    NVEC_HIP_DDATAp(X),
    d_Yd,
    d_buff,
    NVEC_HIP_CONTENT(X)->length
  );
  PostKernelLaunch();

  // Copy GPU result to the cpu.
  err = hipMemcpy(dots, d_buff, grid*sizeof(realtype), hipMemcpyDeviceToHost);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Free host array
  delete[] h_Yd;

  // Free device arrays
  err = hipFree(d_Yd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_buff);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Hip(int nvec, realtype a, N_Vector* X, realtype b,
                                 N_Vector* Y, N_Vector* Z)
{
  hipError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_HIP_DDATAp(X[i]);

  realtype** h_Yd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Yd[i] = NVEC_HIP_DDATAp(Y[i]);

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_HIP_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = hipMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Yd;
  err = hipMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = hipMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream)) return(-1);
  hipLaunchKernelGGL(linearSumVectorArrayKernel, grid, block, shMemSize, stream,
    nvec,
    a,
    d_Xd,
    b,
    d_Yd,
    d_Zd,
    NVEC_HIP_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = hipFree(d_Xd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Yd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Zd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}

int N_VScaleVectorArray_Hip(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  hipError_t err;

  // Copy c array to device
  realtype* d_c;
  err = hipMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_c, c, nvec*sizeof(realtype), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_HIP_DDATAp(X[i]);

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_HIP_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = hipMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = hipMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream)) return(-1);
  hipLaunchKernelGGL(scaleVectorArrayKernel, grid, block, shMemSize, stream,
    nvec,
    d_c,
    d_Xd,
    d_Zd,
    NVEC_HIP_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = hipFree(d_c);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Xd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Zd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}

int N_VConstVectorArray_Hip(int nvec, realtype c, N_Vector* Z)
{
  hipError_t err;

  // Create array of device pointers on host
  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_HIP_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Zd;
  err = hipMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream)) return(-1);
  hipLaunchKernelGGL(constVectorArrayKernel, grid, block, shMemSize, stream,
    nvec,
    c,
    d_Zd,
    NVEC_HIP_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Zd;

  // Free device arrays
  err = hipFree(d_Zd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}

int N_VWrmsNormVectorArray_Hip(int nvec, N_Vector* X, N_Vector* W,
                                realtype* norms)
{
  hipError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_HIP_DDATAp(X[i]);
  realtype** h_Wd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = NVEC_HIP_DDATAp(W[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = hipMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Wd;
  err = hipMalloc((void**) &d_Wd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Wd, h_Wd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(X[0], true, grid, block, shMemSize, stream)) return(-1);
  grid = nvec;

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = hipMalloc((void**) &d_buff, grid*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemsetAsync(d_buff, 0, grid*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  hipLaunchKernelGGL(HIP_KERNEL_NAME(wL2NormSquareVectorArrayKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    nvec,
    d_Xd,
    d_Wd,
    d_buff,
    NVEC_HIP_CONTENT(X[0])->length
  );
  PostKernelLaunch();

  // Copy GPU result to the cpu.
  err = hipMemcpy(norms, d_buff, grid*sizeof(realtype), hipMemcpyDeviceToHost);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Finish computation
  for (int k=0; k<nvec; ++k)
    norms[k] = std::sqrt(norms[k]/NVEC_HIP_CONTENT(X[0])->length);

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;

  // Free device arrays
  err = hipFree(d_Xd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Wd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_buff);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}

int N_VWrmsNormMaskVectorArray_Hip(int nvec, N_Vector* X, N_Vector* W,
                                    N_Vector id, realtype* norms)
{
  hipError_t err;

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_HIP_DDATAp(X[i]);

  realtype** h_Wd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Wd[i] = NVEC_HIP_DDATAp(W[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = hipMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Wd;
  err = hipMalloc((void**) &d_Wd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Wd, h_Wd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(X[0], true, grid, block, shMemSize, stream)) return(-1);
  grid = nvec;

  // Allocate reduction buffer on device
  realtype* d_buff;
  err = hipMalloc((void**) &d_buff, grid*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemsetAsync(d_buff, 0, grid*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  hipLaunchKernelGGL(HIP_KERNEL_NAME(wL2NormSquareMaskVectorArrayKernel<realtype, sunindextype>), grid, block, shMemSize, stream,
    nvec,
    d_Xd,
    d_Wd,
    NVEC_HIP_DDATAp(id),
    d_buff,
    NVEC_HIP_CONTENT(X[0])->length
  );
  PostKernelLaunch();

  // Copy GPU result to the cpu.
  err = hipMemcpy(norms, d_buff, grid*sizeof(realtype), hipMemcpyDeviceToHost);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Finish computation
  for (int k=0; k<nvec; ++k)
    norms[k] = std::sqrt(norms[k]/NVEC_HIP_CONTENT(X[0])->length);

  // Free host array
  delete[] h_Xd;
  delete[] h_Wd;

  // Free device arrays
  err = hipFree(d_Xd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Wd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_buff);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}

int N_VScaleAddMultiVectorArray_Hip(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  hipError_t err;

  // Copy c array to device
  realtype* d_c;
  err = hipMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_c, c, nsum*sizeof(realtype), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Xd[i] = NVEC_HIP_DDATAp(X[i]);

  realtype** h_Yd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Yd[j*nsum+i] = NVEC_HIP_DDATAp(Y[i][j]);

  realtype** h_Zd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Zd[j*nsum+i] = NVEC_HIP_DDATAp(Z[i][j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = hipMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Yd;
  err = hipMalloc((void**) &d_Yd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Yd, h_Yd, nsum*nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = hipMalloc((void**) &d_Zd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Zd, h_Zd, nsum*nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(Z[0][0], false, grid, block, shMemSize, stream)) return(-1);
  hipLaunchKernelGGL(scaleAddMultiVectorArrayKernel, grid, block, shMemSize, stream,
    nvec,
    nsum,
    d_c,
    d_Xd,
    d_Yd,
    d_Zd,
    NVEC_HIP_CONTENT(Z[0][0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = hipFree(d_c);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Xd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Yd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Zd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return(0);
}

int N_VLinearCombinationVectorArray_Hip(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  hipError_t err;

  // Copy c array to device
  realtype* d_c;
  err = hipMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_c, c, nsum*sizeof(realtype), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int i=0; i<nsum; i++)
      h_Xd[j*nsum+i] = NVEC_HIP_DDATAp(X[i][j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int i=0; i<nvec; i++)
    h_Zd[i] = NVEC_HIP_DDATAp(Z[i]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = hipMalloc((void**) &d_Xd, nsum*nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Xd, h_Xd, nsum*nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  realtype** d_Zd;
  err = hipMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), hipMemcpyHostToDevice);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  // Set kernel parameters
  size_t grid, block, shMemSize;
  hipStream_t stream;

  if (GetKernelParameters(Z[0], false, grid, block, shMemSize, stream)) return(-1);
  hipLaunchKernelGGL(linearCombinationVectorArrayKernel, grid, block, shMemSize, stream,
    nvec,
    nsum,
    d_c,
    d_Xd,
    d_Zd,
    NVEC_HIP_CONTENT(Z[0])->length
  );
  PostKernelLaunch();

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = hipFree(d_c);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Xd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);
  err = hipFree(d_Zd);
  if (!SUNDIALS_HIP_VERIFY(err)) return(-1);

  return hipGetLastError();
}


/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */


int N_VBufSize_Hip(N_Vector x, sunindextype *size)
{
  if (x == NULL) return(-1);
  *size = (sunindextype)NVEC_HIP_MEMSIZE(x);
  return(0);
}


int N_VBufPack_Hip(N_Vector x, void *buf)
{
  int copy_fail = 0;
  hipError_t cuerr;

  if (x == NULL || buf == NULL) return(-1);

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return(-1);

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_HIP_MEMHELP(x),
                                        buf_mem,
                                        NVEC_HIP_CONTENT(x)->device_data,
                                        NVEC_HIP_MEMSIZE(x),
                                        (void*) NVEC_HIP_STREAM(x));

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = hipStreamSynchronize(*NVEC_HIP_STREAM(x));

  SUNMemoryHelper_Dealloc(NVEC_HIP_MEMHELP(x), buf_mem, nullptr);

  return (!SUNDIALS_HIP_VERIFY(cuerr) || copy_fail ? -1 : 0);
}


int N_VBufUnpack_Hip(N_Vector x, void *buf)
{
  int copy_fail = 0;
  hipError_t cuerr;

  if (x == NULL || buf == NULL) return(-1);

  SUNMemory buf_mem = SUNMemoryHelper_Wrap(buf, SUNMEMTYPE_HOST);
  if (buf_mem == NULL) return(-1);

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_HIP_MEMHELP(x),
                                        NVEC_HIP_CONTENT(x)->device_data,
                                        buf_mem,
                                        NVEC_HIP_MEMSIZE(x),
                                        (void*) NVEC_HIP_STREAM(x));

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = hipStreamSynchronize(*NVEC_HIP_STREAM(x));

  SUNMemoryHelper_Dealloc(NVEC_HIP_MEMHELP(x), buf_mem, nullptr);

  return (!SUNDIALS_HIP_VERIFY(cuerr) || copy_fail ? -1 : 0);
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf)
  {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Hip;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Hip;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Hip;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Hip;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Hip;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Hip;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Hip;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Hip;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Hip;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Hip;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = N_VDotProdMulti_Hip;
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

int N_VEnableLinearCombination_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Hip;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Hip;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf) {
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Hip;
    v->ops->nvdotprodmultilocal = N_VDotProdMulti_Hip;
  } else {
    v->ops->nvdotprodmulti      = NULL;
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Hip;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Hip;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Hip;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Hip;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

/* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Hip;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Hip;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Hip(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Hip;
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
  N_VectorContent_Hip vc = NVEC_HIP_CONTENT(v);
  N_PrivateVectorContent_Hip vcp = NVEC_HIP_PRIVATE(v);

  if (N_VGetLength_Hip(v) == 0) return(0);

  if (vcp->use_managed_mem)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_HIP_MEMHELP(v), &(vc->device_data),
                                       NVEC_HIP_MEMSIZE(v), SUNMEMTYPE_UVM,
                                       nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed for SUNMEMTYPE_UVM\n");
    }
    vc->host_data = SUNMemoryHelper_Alias(vc->device_data);
  }
  else
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_HIP_MEMHELP(v), &(vc->host_data),
                                       NVEC_HIP_MEMSIZE(v), SUNMEMTYPE_HOST,
                                       nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in AllocateData: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
    }

    alloc_fail = SUNMemoryHelper_Alloc(NVEC_HIP_MEMHELP(v), &(vc->device_data),
                                       NVEC_HIP_MEMSIZE(v), SUNMEMTYPE_DEVICE,
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
int InitializeReductionBuffer(N_Vector v, const realtype value, size_t n)
{
  int alloc_fail = 0, copy_fail = 0;
  size_t bytes = sizeof(realtype);
  booleantype need_to_allocate = SUNFALSE;
  N_PrivateVectorContent_Hip vcp = NVEC_HIP_PRIVATE(v);
  SUNMemory value_mem = SUNMemoryHelper_Wrap((void*) &value, SUNMEMTYPE_HOST);

  /* we allocate if the existing reduction buffer is not large enough */
  if (vcp->reduce_buffer_allocated_bytes < bytes)
  {
    FreeReductionBuffer(v);
    need_to_allocate = SUNTRUE;
  }

  if (need_to_allocate)
  {
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_HIP_MEMHELP(v),
                                       &(vcp->reduce_buffer_host), bytes,
                                       SUNMEMTYPE_PINNED, nullptr);
    if (alloc_fail)
    {
      SUNDIALS_DEBUG_PRINT("WARNING in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_PINNED, using SUNMEMTYPE_HOST instead\n");

      /* try to allocate just plain host memory instead */
      alloc_fail = SUNMemoryHelper_Alloc(NVEC_HIP_MEMHELP(v),
                                         &(vcp->reduce_buffer_host), bytes,
                                         SUNMEMTYPE_HOST, nullptr);
      if (alloc_fail)
      {
        SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_Alloc failed to alloc SUNMEMTYPE_HOST\n");
      }
    }
    alloc_fail = SUNMemoryHelper_Alloc(NVEC_HIP_MEMHELP(v),
                                       &(vcp->reduce_buffer_dev), bytes,
                                       SUNMEMTYPE_DEVICE, nullptr);
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
    copy_fail = SUNMemoryHelper_CopyAsync(NVEC_HIP_MEMHELP(v),
                                          vcp->reduce_buffer_dev, value_mem,
                                          bytes, (void*) NVEC_HIP_STREAM(v));

    if (copy_fail)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in InitializeReductionBuffer: SUNMemoryHelper_CopyAsync failed\n");
    }
  }

  SUNMemoryHelper_Dealloc(NVEC_HIP_MEMHELP(v), value_mem, nullptr);
  return((alloc_fail || copy_fail) ? -1 : 0);
}

/* Free the reduction buffer
 */
void FreeReductionBuffer(N_Vector v)
{
  N_PrivateVectorContent_Hip vcp = NVEC_HIP_PRIVATE(v);

  if (vcp == NULL) return;

  if (vcp->reduce_buffer_dev != NULL)
    SUNMemoryHelper_Dealloc(NVEC_HIP_MEMHELP(v), vcp->reduce_buffer_dev,
                            nullptr);
  vcp->reduce_buffer_dev  = NULL;
  if (vcp->reduce_buffer_host != NULL)
    SUNMemoryHelper_Dealloc(NVEC_HIP_MEMHELP(v), vcp->reduce_buffer_host,
                            nullptr);
  vcp->reduce_buffer_host = NULL;
}

/* Copy the reduction buffer from the device to the host.
 */
int CopyReductionBufferFromDevice(N_Vector v, size_t n)
{
  int copy_fail;
  hipError_t cuerr;

  copy_fail = SUNMemoryHelper_CopyAsync(NVEC_HIP_MEMHELP(v),
                                        NVEC_HIP_PRIVATE(v)->reduce_buffer_host,
                                        NVEC_HIP_PRIVATE(v)->reduce_buffer_dev,
                                        n*sizeof(realtype),
                                        (void*) NVEC_HIP_STREAM(v));

  if (copy_fail)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in CopyReductionBufferFromDevice: SUNMemoryHelper_CopyAsync returned nonzero\n");
  }

  /* we synchronize with respect to the host, but only in this stream */
  cuerr = hipStreamSynchronize(*NVEC_HIP_STREAM(v));
  return (!SUNDIALS_HIP_VERIFY(cuerr) || copy_fail ? -1 : 0);
}

/* Get the kernel launch parameters based on the kernel type (reduction or not),
 * using the appropriate kernel execution policy.
 */
int GetKernelParameters(N_Vector v, booleantype reduction, size_t& grid, size_t& block,
                        size_t& shMemSize, hipStream_t& stream, size_t n)
{
  n = (n == 0) ? NVEC_HIP_CONTENT(v)->length : n;
  if (reduction)
  {
    SUNHipExecPolicy* reduce_exec_policy = NVEC_HIP_CONTENT(v)->reduce_exec_policy;
    grid      = reduce_exec_policy->gridSize(n);
    block     = reduce_exec_policy->blockSize();
    shMemSize = 0;
    stream    = *(reduce_exec_policy->stream());
    if (block % warpSize)
    {
#ifdef SUNDIALS_DEBUG
      throw std::runtime_error("the block size must be a multiple must be of HIP warp size");
#endif
      return(-1);
    }
  }
  else
  {
    SUNHipExecPolicy* stream_exec_policy = NVEC_HIP_CONTENT(v)->stream_exec_policy;
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
 * If SUNDIALS_DEBUG_HIP_LASTERROR is not defined, then the function does nothing.
 * If it is defined, the function will synchronize and check the last HIP error.
 */
void PostKernelLaunch()
{
#ifdef SUNDIALS_DEBUG_HIP_LASTERROR
  hipDeviceSynchronize();
  SUNDIALS_HIP_VERIFY(hipGetLastError());
#endif
}


} // extern "C"
