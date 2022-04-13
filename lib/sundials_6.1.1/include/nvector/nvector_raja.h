/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, Cody J. Balos, Daniel McGreer @ LLNL
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
 * This is the header file for the RAJA implementation of the
 * NVECTOR module.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_RAJA_H
#define _NVECTOR_RAJA_H

#include <stdio.h>

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_memory.h>
#include <sundials/sundials_config.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * RAJA implementation of N_Vector
 * -----------------------------------------------------------------
 */

/* RAJA implementation of the N_Vector 'content' structure
   contains the length of the vector, pointers to host and device
   arrays of 'realtype' components, a flag indicating ownership of
   the data, and a private data pointer  */

struct _N_VectorContent_Raja {
  sunindextype    length;
  booleantype     own_helper;
  SUNMemory       host_data;
  SUNMemory       device_data;
  SUNMemoryHelper mem_helper;
  void*           priv; /* 'private' data */
};

typedef struct _N_VectorContent_Raja *N_VectorContent_Raja;

/*
 * -----------------------------------------------------------------
 * NVECTOR_RAJA implementation specific functions
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Raja(SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNew_Raja(sunindextype length, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewManaged_Raja(sunindextype length,
                                            SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VNewWithMemHelp_Raja(sunindextype length,
                                                booleantype use_managed_mem,
                                                SUNMemoryHelper helper,
                                                SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMake_Raja(sunindextype length, realtype *h_vdata,
                                      realtype *d_vdata, SUNContext sunctx);
SUNDIALS_EXPORT N_Vector N_VMakeManaged_Raja(sunindextype length,
                                             realtype *vdata,
                                             SUNContext sunctx);
SUNDIALS_EXPORT void N_VSetHostArrayPointer_Raja(realtype* h_vdata, N_Vector v);
SUNDIALS_EXPORT void N_VSetDeviceArrayPointer_Raja(realtype* d_vdata,
                                                   N_Vector v);
SUNDIALS_EXPORT booleantype N_VIsManagedMemory_Raja(N_Vector x);
SUNDIALS_EXPORT void N_VCopyToDevice_Raja(N_Vector v);
SUNDIALS_EXPORT void N_VCopyFromDevice_Raja(N_Vector v);

SUNDIALS_STATIC_INLINE
sunindextype N_VGetLength_Raja(N_Vector x)
{
  N_VectorContent_Raja content = (N_VectorContent_Raja)x->content;
  return content->length;
}

SUNDIALS_STATIC_INLINE
realtype *N_VGetHostArrayPointer_Raja(N_Vector x)
{
  N_VectorContent_Raja content = (N_VectorContent_Raja)x->content;
  return(content->host_data == NULL ? NULL : (realtype*)content->host_data->ptr);
}

SUNDIALS_STATIC_INLINE
realtype *N_VGetDeviceArrayPointer_Raja(N_Vector x)
{
  N_VectorContent_Raja content = (N_VectorContent_Raja)x->content;
  return(content->device_data == NULL ? NULL : (realtype*)content->device_data->ptr);
}


/*
 * -----------------------------------------------------------------
 * NVECTOR API functions
 * -----------------------------------------------------------------
 */

SUNDIALS_STATIC_INLINE
N_Vector_ID N_VGetVectorID_Raja(N_Vector v)
{
  return SUNDIALS_NVEC_RAJA;
}

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Raja(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Raja(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Raja(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Raja(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT void N_VSetArrayPointer_Raja(realtype *v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Raja(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Raja(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Raja(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Raja(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Raja(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Raja(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Raja(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Raja(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Raja(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Raja(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Raja(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Raja(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Raja(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Raja(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Raja(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Raja(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Raja(int nvec, realtype* c, N_Vector* X,
                                              N_Vector z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Raja(int nvec, realtype* c, N_Vector x,
                                          N_Vector* Y, N_Vector* Z);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Raja(int nvec,
                                                 realtype a, N_Vector* X,
                                                 realtype b, N_Vector* Y,
                                                 N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Raja(int nvec, realtype* c, N_Vector* X,
                                             N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Raja(int nvec, realtype c, N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Raja(int nvec, int nsum,
                                                     realtype* a,
                                                     N_Vector* X, N_Vector** Y,
                                                     N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Raja(int nvec, int nsum,
                                                         realtype* c,
                                                         N_Vector** X,
                                                         N_Vector* Z);

/* OPTIONAL local reduction kernels (no parallel communication) */
SUNDIALS_EXPORT realtype N_VWSqrSumLocal_Raja(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWSqrSumMaskLocal_Raja(N_Vector x, N_Vector w, N_Vector id);

/* OPTIONAL XBraid interface operations */
SUNDIALS_EXPORT int N_VBufSize_Raja(N_Vector x, sunindextype *size);
SUNDIALS_EXPORT int N_VBufPack_Raja(N_Vector x, void *buf);
SUNDIALS_EXPORT int N_VBufUnpack_Raja(N_Vector x, void *buf);

/* OPTIONAL operations for debugging */
SUNDIALS_EXPORT void N_VPrint_Raja(N_Vector v);
SUNDIALS_EXPORT void N_VPrintFile_Raja(N_Vector v, FILE *outfile);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Raja(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Raja(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Raja(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Raja(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Raja(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Raja(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Raja(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Raja(N_Vector v, booleantype tf);

#ifdef __cplusplus
}
#endif

#endif
