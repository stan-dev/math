/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
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
 * This is the header file for the serial implementation of the
 * NVECTOR module.
 *
 * Part I contains declarations specific to the serial
 * implementation of the supplied NVECTOR module.
 *
 * Part II defines accessor macros that allow the user to
 * efficiently use the type N_Vector without making explicit
 * references to the underlying data structure.
 *
 * Part III contains the prototype for the constructor N_VNew_Serial
 * as well as implementation-specific prototypes for various useful
 * vector operations.
 *
 * Notes:
 *
 *   - The definition of the generic N_Vector structure can be found
 *     in the header file sundials_nvector.h.
 *
 *   - The definition of the type 'realtype' can be found in the
 *     header file sundials_types.h, and it may be changed (at the 
 *     configuration stage) according to the user's needs. 
 *     The sundials_types.h file also contains the definition
 *     for the type 'booleantype'.
 *
 *   - N_Vector arguments to arithmetic vector operations need not
 *     be distinct. For example, the following call:
 *
 *       N_VLinearSum_Serial(a,x,b,y,y);
 *
 *     (which stores the result of the operation a*x+b*y in y)
 *     is legal.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_SERIAL_H
#define _NVECTOR_SERIAL_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/*
 * -----------------------------------------------------------------
 * PART I: SERIAL implementation of N_Vector
 * -----------------------------------------------------------------
 */

/* serial implementation of the N_Vector 'content' structure
   contains the length of the vector, a pointer to an array
   of 'realtype' components, and a flag indicating ownership of
   the data */

struct _N_VectorContent_Serial {
  sunindextype length;
  booleantype own_data;
  realtype *data;
};

typedef struct _N_VectorContent_Serial *N_VectorContent_Serial;

/*
 * -----------------------------------------------------------------
 * PART II: macros NV_CONTENT_S, NV_DATA_S, NV_OWN_DATA_S,
 *          NV_LENGTH_S, and NV_Ith_S
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * sunindextype i;
 *
 * (1) NV_CONTENT_S
 *
 *     This routines gives access to the contents of the serial
 *     vector N_Vector.
 *
 *     The assignment v_cont = NV_CONTENT_S(v) sets v_cont to be
 *     a pointer to the serial N_Vector content structure.
 *
 * (2) NV_DATA_S NV_OWN_DATA_S and NV_LENGTH_S
 *
 *     These routines give access to the individual parts of
 *     the content structure of a serial N_Vector.
 *
 *     The assignment v_data = NV_DATA_S(v) sets v_data to be
 *     a pointer to the first component of v. The assignment
 *     NV_DATA_S(v) = data_V sets the component array of v to
 *     be data_v by storing the pointer data_v.
 *
 *     The assignment v_len = NV_LENGTH_S(v) sets v_len to be
 *     the length of v. The call NV_LENGTH_S(v) = len_v sets
 *     the length of v to be len_v.
 *
 * (3) NV_Ith_S
 *
 *     In the following description, the components of an
 *     N_Vector are numbered 0..n-1, where n is the length of v.
 *
 *     The assignment r = NV_Ith_S(v,i) sets r to be the value of
 *     the ith component of v. The assignment NV_Ith_S(v,i) = r
 *     sets the value of the ith component of v to be r.
 *
 * Note: When looping over the components of an N_Vector v, it is
 * more efficient to first obtain the component array via
 * v_data = NV_DATA_S(v) and then access v_data[i] within the
 * loop than it is to use NV_Ith_S(v,i) within the loop.
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_S(v)  ( (N_VectorContent_Serial)(v->content) )

#define NV_LENGTH_S(v)   ( NV_CONTENT_S(v)->length )

#define NV_OWN_DATA_S(v) ( NV_CONTENT_S(v)->own_data )

#define NV_DATA_S(v)     ( NV_CONTENT_S(v)->data )

#define NV_Ith_S(v,i)    ( NV_DATA_S(v)[i] )

/*
 * -----------------------------------------------------------------
 * PART III: functions exported by nvector_serial
 * 
 * CONSTRUCTORS:
 *    N_VNew_Serial
 *    N_VNewEmpty_Serial
 *    N_VMake_Serial
 *    N_VCloneVectorArray_Serial
 *    N_VCloneVectorArrayEmpty_Serial
 * DESTRUCTORS:
 *    N_VDestroy_Serial
 *    N_VDestroyVectorArray_Serial
 * ENABLE/DISABLE FUSED OPS:
 *    N_VEnableFusedOps_Serial
 *    N_VEnableLinearCombination_Serial
 *    N_VEnableScaleAddMulti_Serial
 *    N_VEnableDotProdMulti_Serial
 *    N_VEnableLinearSumVectorArray_Serial
 *    N_VEnableScaleVectorArray_Serial
 *    N_VEnableConstVectorArray_Serial
 *    N_VEnableWrmsNormVectorArray_Serial
 *    N_VEnableWrmsNormMaskVectorArray_Serial
 *    N_VEnableScaleAddMultiVectorArray_Serial
 *    N_VEnableLinearCombinationVectorArray_Serial
 * OTHER:
 *    N_VGetLength_Serial
 *    N_VPrint_Serial
 *    N_VPrintFile_Serial
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : N_VNew_Serial
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a serial vector.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNew_Serial(sunindextype vec_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VNewEmpty_Serial
 * -----------------------------------------------------------------
 * This function creates a new serial N_Vector with an empty (NULL)
 * data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VNewEmpty_Serial(sunindextype vec_length);

/*
 * -----------------------------------------------------------------
 * Function : N_VMake_Serial
 * -----------------------------------------------------------------
 * This function creates and allocates memory for a serial vector
 * with a user-supplied data array.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector N_VMake_Serial(sunindextype vec_length, realtype *v_data);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArray_Serial
 * -----------------------------------------------------------------
 * This function creates an array of 'count' SERIAL vectors by
 * cloning a given vector w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArray_Serial(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VCloneVectorArrayEmpty_Serial
 * -----------------------------------------------------------------
 * This function creates an array of 'count' SERIAL vectors each
 * with an empty (NULL) data array by cloning w.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector *N_VCloneVectorArrayEmpty_Serial(int count, N_Vector w);

/*
 * -----------------------------------------------------------------
 * Function : N_VDestroyVectorArray_Serial
 * -----------------------------------------------------------------
 * This function frees an array of SERIAL vectors created with 
 * N_VCloneVectorArray_Serial or N_VCloneVectorArrayEmpty_Serial.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VDestroyVectorArray_Serial(N_Vector *vs, int count);

/*
 * -----------------------------------------------------------------
 * Function : N_VGetLength_Serial
 * -----------------------------------------------------------------
 * This function returns number of vector elements.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT sunindextype N_VGetLength_Serial(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrint_Serial
 * -----------------------------------------------------------------
 * This function prints the content of a serial vector to stdout.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrint_Serial(N_Vector v);

/*
 * -----------------------------------------------------------------
 * Function : N_VPrintFile_Serial
 * -----------------------------------------------------------------
 * This function prints the content of a serial vector to outfile.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT void N_VPrintFile_Serial(N_Vector v, FILE *outfile);

/*
 * -----------------------------------------------------------------
 * serial implementations of various useful vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_Serial(N_Vector v);
SUNDIALS_EXPORT N_Vector N_VCloneEmpty_Serial(N_Vector w);
SUNDIALS_EXPORT N_Vector N_VClone_Serial(N_Vector w);
SUNDIALS_EXPORT void N_VDestroy_Serial(N_Vector v);
SUNDIALS_EXPORT void N_VSpace_Serial(N_Vector v, sunindextype *lrw, sunindextype *liw);
SUNDIALS_EXPORT realtype *N_VGetArrayPointer_Serial(N_Vector v);
SUNDIALS_EXPORT void N_VSetArrayPointer_Serial(realtype *v_data, N_Vector v);

/* standard vector operations */
SUNDIALS_EXPORT void N_VLinearSum_Serial(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VConst_Serial(realtype c, N_Vector z);
SUNDIALS_EXPORT void N_VProd_Serial(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VDiv_Serial(N_Vector x, N_Vector y, N_Vector z);
SUNDIALS_EXPORT void N_VScale_Serial(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAbs_Serial(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VInv_Serial(N_Vector x, N_Vector z);
SUNDIALS_EXPORT void N_VAddConst_Serial(N_Vector x, realtype b, N_Vector z);
SUNDIALS_EXPORT realtype N_VDotProd_Serial(N_Vector x, N_Vector y);
SUNDIALS_EXPORT realtype N_VMaxNorm_Serial(N_Vector x);
SUNDIALS_EXPORT realtype N_VWrmsNorm_Serial(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VWrmsNormMask_Serial(N_Vector x, N_Vector w, N_Vector id);
SUNDIALS_EXPORT realtype N_VMin_Serial(N_Vector x);
SUNDIALS_EXPORT realtype N_VWL2Norm_Serial(N_Vector x, N_Vector w);
SUNDIALS_EXPORT realtype N_VL1Norm_Serial(N_Vector x);
SUNDIALS_EXPORT void N_VCompare_Serial(realtype c, N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VInvTest_Serial(N_Vector x, N_Vector z);
SUNDIALS_EXPORT booleantype N_VConstrMask_Serial(N_Vector c, N_Vector x, N_Vector m);
SUNDIALS_EXPORT realtype N_VMinQuotient_Serial(N_Vector num, N_Vector denom);

/* fused vector operations */
SUNDIALS_EXPORT int N_VLinearCombination_Serial(int nvec, realtype* c, N_Vector* V,
                                                N_Vector z);
SUNDIALS_EXPORT int N_VScaleAddMulti_Serial(int nvec, realtype* a, N_Vector x,
                                            N_Vector* Y, N_Vector* Z);
SUNDIALS_EXPORT int N_VDotProdMulti_Serial(int nvec, N_Vector x,
                                           N_Vector *Y, realtype* dotprods);

/* vector array operations */
SUNDIALS_EXPORT int N_VLinearSumVectorArray_Serial(int nvec, 
                                                   realtype a, N_Vector* X,
                                                   realtype b, N_Vector* Y,
                                                   N_Vector* Z);
SUNDIALS_EXPORT int N_VScaleVectorArray_Serial(int nvec, realtype* c,
                                               N_Vector* X, N_Vector* Z);
SUNDIALS_EXPORT int N_VConstVectorArray_Serial(int nvecs, realtype c,
                                               N_Vector* Z);
SUNDIALS_EXPORT int N_VWrmsNormVectorArray_Serial(int nvecs, N_Vector* X,
                                                  N_Vector* W, realtype* nrm);
SUNDIALS_EXPORT int N_VWrmsNormMaskVectorArray_Serial(int nvecs, N_Vector* X,
                                                      N_Vector* W, N_Vector id,
                                                      realtype* nrm);
SUNDIALS_EXPORT int N_VScaleAddMultiVectorArray_Serial(int nvec, int nsum,
                                                       realtype* a,
                                                       N_Vector* X,
                                                       N_Vector** Y,
                                                       N_Vector** Z);
SUNDIALS_EXPORT int N_VLinearCombinationVectorArray_Serial(int nvec, int nsum,
                                                           realtype* c,
                                                           N_Vector** X,
                                                           N_Vector* Z);

/*
 * -----------------------------------------------------------------
 * Enable / disable fused vector operations
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int N_VEnableFusedOps_Serial(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearCombination_Serial(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMulti_Serial(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableDotProdMulti_Serial(N_Vector v, booleantype tf);

SUNDIALS_EXPORT int N_VEnableLinearSumVectorArray_Serial(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleVectorArray_Serial(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableConstVectorArray_Serial(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormVectorArray_Serial(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableWrmsNormMaskVectorArray_Serial(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableScaleAddMultiVectorArray_Serial(N_Vector v, booleantype tf);
SUNDIALS_EXPORT int N_VEnableLinearCombinationVectorArray_Serial(N_Vector v, booleantype tf);

#ifdef __cplusplus
}
#endif

#endif
