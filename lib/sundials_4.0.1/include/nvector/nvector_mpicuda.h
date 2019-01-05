/* -----------------------------------------------------------------
 * Programmer(s): Cody Balos @ LLNL
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
 * This is the header file for the MPI+CUDA implementation of the
 * NVECTOR module.
 *
 * Part I contains declarations specific to the MPI+CUDA
 * implementation of the supplied NVECTOR module.
 *
 * Part II contains the prototype for the constructor N_VNew_Cuda
 * as well as implementation-specific prototypes for various useful
 * vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be found
 *     in the header file sundials_nvector.h.
 *
 *   - The definitions of the types 'realtype' and 'sunindextype' can
 *     be found in the header file sundials_types.h, and it may be
 *     changed (at the configuration stage) according to the user's needs.
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype'.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *       N_VLinearSum_Cuda(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_CUDA_H
#define _NVECTOR_CUDA_H

#include <mpi.h>
#include <stdio.h>

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_config.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif



/*
 * -----------------------------------------------------------------
 * PART I: CUDA implementation of N_Vector
 * -----------------------------------------------------------------
 */

/*
 * CUDA implementation of the N_Vector 'content' is in C++ class
 * Vector. The class inherits from structure _N_VectorContent_Cuda
 * to create C <--> C++ interface.
 */

struct _N_VectorContent_Cuda {};

typedef struct _N_VectorContent_Cuda *N_VectorContent_Cuda;




/*
 * -----------------------------------------------------------------
 * PART II: functions exported by nvector_cuda
 *
 * CONSTRUCTORS:
 *    N_VNew_Cuda
 *    N_VNewManaged_Cuda
 *    N_VNewEmpty_Cuda
 *    N_VMake_Cuda
 *    N_VMakeManaged_Cuda
 * DESTRUCTORS:
 *    N_VDestroy_Cuda
 * ENABLE/DISABLE FUSED OPS:
 *    N_VEnableFusedOps_Cuda
 *    N_VEnableLinearCombination_Cuda
 *    N_VEnableScaleAddMulti_Cuda
 *    N_VEnableDotProdMulti_Cuda
 *    N_VEnableLinearSumVectorArray_Cuda
 *    N_VEnableScaleVectorArray_Cuda
 *    N_VEnableConstVectorArray_Cuda
 *    N_VEnableWrmsNormVectorArray_Cuda
 *    N_VEnableWrmsNormMaskVectorArray_Cuda
 *    N_VEnableScaleAddMultiVectorArray_Cuda
 *    N_VEnableLinearCombinationVectorArray_Cuda
 * OTHER:
 *    N_VGetLength_Cuda
 *    N_VGetLocalLength_Cuda
 *    N_VGetHostArrayPointer_Cuda
 *    N_VGetDeviceArrayPointer_Cuda
 *    N_VGetMPIComm_Cuda
 *    N_VIsManagedMemory_Cuda
 *    N_VSetCudaStream_Cuda
 *    N_VPrint_Cuda
 *    N_VPrintFile_Cuda
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNew_Cuda
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a distributed
 * memory CUDA vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_Cuda(MPI_Comm comm,
                                     sunindextype local_length,
                                     sunindextype global_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VNewManaged_Cuda
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a distributed
 * memory CUDA vector with managed memory vector data.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewManaged_Cuda(MPI_Comm comm,
                                            sunindextype local_length,
                                            sunindextype global_length);
/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_Cuda
 * -----------------------------------------------------------------
 * This function creates a new CUDA N_Vector with an empty (NULL)
 * content.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Cuda();

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_Cuda
 * -----------------------------------------------------------------
 * This function creates a CUDA vector with user-supplied data 
 * arrays. Both h_vdata and d_vdata must be non-NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_Cuda(MPI_Comm comm,
                                      sunindextype local_length,
                                      sunindextype global_length,
                                      realtype *h_vdata,
                                      realtype *d_vdata);

/*
 * -----------------------------------------------------------------
 * Function : N_VMakeManaged_Cuda
 * -----------------------------------------------------------------
 * This function creates a CUDA vector with a user-supplied managed
 * memory data array. vdata must be non-NULL.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMakeManaged_Cuda(MPI_Comm comm,
                                             sunindextype local_length,
                                             sunindextype global_length,
                                             realtype *vdata);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetLength_Cuda
 * -----------------------------------------------------------------
 * This function returns the global length of the vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunindextype N_VGetLength_Cuda(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetLocalLength_Cuda
 * -----------------------------------------------------------------
 * This function returns the local length of the vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunindextype N_VGetLocalLength_Cuda(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetMPIComm_Cuda
 * -----------------------------------------------------------------
 * This function returns the MPI communicator for the vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT MPI_Comm N_VGetMPIComm_Cuda(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetHostArrayPointer_Cuda
 * -----------------------------------------------------------------
 * This function returns pointer to the host raw data.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT realtype *N_VGetHostArrayPointer_Cuda(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetDeviceArrayPointer_Cuda
 * -----------------------------------------------------------------
 * This function returns pointer to the device raw data.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT realtype *N_VGetDeviceArrayPointer_Cuda(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VIsManagedMemory_Cuda
 * -----------------------------------------------------------------
 * This function returns a boolean flag indicating if the vector
 * data is managed memory.
 * -----------------------------------------------------------------
 */
SUNDIALS_EXPORT booleantype N_VIsManagedMemory_Cuda(N_Vector x);

/*
 * -----------------------------------------------------------------
 * Function : N_VSetCudaStream_Cuda
 * -----------------------------------------------------------------
 * This function sets the cudaStream_t to use for execution of 
 * the CUDA kernels.
 * -----------------------------------------------------------------
 */
SUNDIALS_EXPORT void N_VSetCudaStream_Cuda(N_Vector x, cudaStream_t *stream);

/*
 * -----------------------------------------------------------------
 * Function : N_VCopyTotDevice_Cuda
 * -----------------------------------------------------------------
 * This function copies host data to the device.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VCopyToDevice_Cuda(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VCopyTotDevice_Cuda
 * -----------------------------------------------------------------
 * This function copies vector data from the device to the host.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VCopyFromDevice_Cuda(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_Cuda
 * -----------------------------------------------------------------
 * This function prints the content of a CUDA vector to stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_Cuda(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrintFile_Cuda
 * -----------------------------------------------------------------
 * This function prints the content of a CUDA vector to outfile.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrintFile_Cuda(N_Vector v, FILE *outfile);

/*
 * -----------------------------------------------------------------
 * CUDA implementations of various useful vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Cuda(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Cuda(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Cuda(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Cuda(N_Vector v, sunindextype *lrw, sunindextype *liw);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Cuda(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Cuda(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Cuda(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Cuda(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Cuda(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Cuda(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Cuda(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Cuda(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Cuda(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Cuda(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Cuda(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Cuda(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Cuda(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Cuda(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Cuda(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Cuda(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Cuda(int nvec, realtype* c, N_Vector* X,
                                              N_Vector Z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Cuda(int nvec, realtype* c, N_Vector X,
                                          N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y,
                                         realtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Cuda(int nvec,
                                                 realtype a, N_Vector* X,
                                                 realtype b, N_Vector* Y,
                                                 N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Cuda(int nvec, realtype* c, N_Vector* X,
                                             N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Cuda(int nvec, realtype c, N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Cuda(int nvec, int nsum,
                                                     realtype* a, N_Vector* X,
                                                     N_Vector** Y, N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Cuda(int nvec, int nsum,
                                                         realtype* c,
                                                         N_Vector** X,
                                                         N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X,
                                                N_Vector* W, realtype* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X,
                                                    N_Vector* W, N_Vector id,
                                                    realtype* nrm);


/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Cuda(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf);

#ifdef __cplusplus
}
#endif

#endif
