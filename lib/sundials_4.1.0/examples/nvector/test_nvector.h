/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file contains the prototypes for functions to
 * test an NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#include <math.h>

/* define constants */
#define NEG_TWO  RCONST(-2.0)
#define NEG_ONE  RCONST(-1.0)
#define NEG_HALF RCONST(-0.5)
#define ZERO     RCONST(0.0)
#define HALF     RCONST(0.5)
#define ONE      RCONST(1.0)
#define TWO      RCONST(2.0)

/* NAN and floating point "equality" check, failure update macro */
#if __STDC_VERSION__ >= 199901L
#define FNEQ(a,b) (isnan(a) ? 1 : ( SUNRabs((a)-(b))/SUNRabs(b) > 1.0e-15 ))
#else
#define FNEQ(a,b) (( SUNRabs((a)-(b))/SUNRabs(b) > 1.0e-15 ))
#endif


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  /* Forward declarations for implementation specific utility functions */
  int check_ans(realtype ans, N_Vector X, sunindextype local_length);
  booleantype has_data(N_Vector X);
  void set_element(N_Vector X, sunindextype i, realtype val);
  realtype get_element(N_Vector X, sunindextype i);
  double max_time(N_Vector X, double time);
  void sync_device();

  /* Vector ID test */
  int Test_N_VGetVectorID(N_Vector X, N_Vector_ID ID, int myid);

  /* Vector constructor tests */
  int Test_N_VMake(N_Vector X, sunindextype local_length, int myid);

  /* Vector clone tests */
  int Test_N_VCloneVectorArray(int count, N_Vector W, sunindextype local_length, int myid);
  int Test_N_VCloneEmptyVectorArray(int count, N_Vector W, int myid);
  int Test_N_VCloneEmpty(N_Vector W, int myid);
  int Test_N_VClone(N_Vector W, sunindextype local_length, int myid);

  /* Vector get/set function tests */
  int Test_N_VGetArrayPointer(N_Vector W, sunindextype local_length, int myid);
  int Test_N_VSetArrayPointer(N_Vector W, sunindextype local_length, int myid);

  /* Standard vector operation tests */
  int Test_N_VLinearSum(N_Vector X, N_Vector Y, N_Vector Z, 
                        sunindextype local_length, int myid);
  int Test_N_VConst(N_Vector X, sunindextype local_length, int myid);
  int Test_N_VProd(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length, int myid);
  int Test_N_VDiv(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length, int myid);
  int Test_N_VScale(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
  int Test_N_VAbs(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
  int Test_N_VInv(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
  int Test_N_VAddConst(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
  int Test_N_VDotProd(N_Vector X, N_Vector Y, 
                      sunindextype local_length, sunindextype global_length, int myid);
  int Test_N_VMaxNorm(N_Vector X, sunindextype local_length, int myid); 
  int Test_N_VWrmsNorm(N_Vector X, N_Vector W, sunindextype local_length, int myid);
  int Test_N_VWrmsNormMask(N_Vector X, N_Vector W, N_Vector ID, 
                           sunindextype local_length, sunindextype global_length, int myid);
  int Test_N_VMin(N_Vector X, sunindextype local_length, int myid); 
  int Test_N_VWL2Norm(N_Vector X, N_Vector W, 
                      sunindextype local_length, sunindextype global_length, int myid);
  int Test_N_VL1Norm(N_Vector X, sunindextype local_length, sunindextype global_length, int myid);
  int Test_N_VCompare(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
  int Test_N_VInvTest(N_Vector X, N_Vector Z, sunindextype local_length, int myid);
  int Test_N_VConstrMask(N_Vector C, N_Vector X, N_Vector M, 
                         sunindextype local_length, int myid);
  int Test_N_VMinQuotient(N_Vector NUM, N_Vector DENOM, sunindextype local_length, int myid);

  /* Fused vector operation tests */
  int Test_N_VLinearCombination(N_Vector X, sunindextype local_length, int myid);
  int Test_N_VScaleAddMulti(N_Vector X, sunindextype local_length, int myid);
  int Test_N_VDotProdMulti(N_Vector X, sunindextype local_length,
                           sunindextype global_length, int myid);

  /* Vector array operation tests */
  int Test_N_VLinearSumVectorArray(N_Vector X, sunindextype local_length, int myid);
  int Test_N_VScaleVectorArray(N_Vector X, sunindextype local_length, int myid);
  int Test_N_VConstVectorArray(N_Vector X, sunindextype local_length, int myid);
  int Test_N_VWrmsNormVectorArray(N_Vector X, sunindextype local_length, int myid);
  int Test_N_VWrmsNormMaskVectorArray(N_Vector X, sunindextype local_length,
                                      sunindextype global_length, int myid);
  int Test_N_VScaleAddMultiVectorArray(N_Vector X, sunindextype local_length, int myid);
  int Test_N_VLinearCombinationVectorArray(N_Vector X, sunindextype local_length, int myid);

  void SetTiming(int onoff, int myid);
  
#ifdef __cplusplus
}
#endif
