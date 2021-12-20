/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the header file contains the prototypes for functions to
 * evaluate the performance of an NVECTOR module implementation.
 * -----------------------------------------------------------------*/

#include <math.h>

/* define constatnts */
#define NEG_ONE  RCONST(-1.0)
#define ZERO     RCONST(0.0)
#define ONE      RCONST(1.0)
#define TWO      RCONST(2.0)
#define TEN      RCONST(10.0)

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  /* Forward declarations for implementation specific utility functions */
  void collect_times(N_Vector X, double *times, int ntimes);
  void sync_device(N_Vector x);
  void ClearCache();
  void N_VRand(N_Vector Xvec, sunindextype Xlen, realtype lower, realtype upper);
  void N_VRandZeroOne(N_Vector Xvec, sunindextype Xlen);
  void N_VRandConstraints(N_Vector Xvec, sunindextype Xlen);

  /* Standard vector operation tests */
  int Test_N_VLinearSum(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VConst(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VProd(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VDiv(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VScale(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VAbs(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VInv(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VAddConst(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VDotProd(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VMaxNorm(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VWrmsNorm(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VWrmsNormMask(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VMin(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VWL2Norm(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VL1Norm(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VCompare(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VInvTest(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VConstrMask(N_Vector X, sunindextype local_length, int ntests);
  int Test_N_VMinQuotient(N_Vector X, sunindextype local_length, int ntests);

  /* Fused vector operation tests */
  int Test_N_VLinearCombination(N_Vector X, sunindextype local_length,
                                int nvecs, int ntests);
  int Test_N_VScaleAddMulti(N_Vector X, sunindextype local_length,
                            int nvecs, int ntests);
  int Test_N_VDotProdMulti(N_Vector X, sunindextype local_length,
                           int nvecs, int ntests);

  /* Vector array operation tests */
  int Test_N_VLinearSumVectorArray(N_Vector X, sunindextype local_length,
                                   int nvecs, int ntests);
  int Test_N_VScaleVectorArray(N_Vector X, sunindextype local_length,
                               int nvecs, int ntests);
  int Test_N_VConstVectorArray(N_Vector X, sunindextype local_length,
                               int nvecs, int ntests);
  int Test_N_VWrmsNormVectorArray(N_Vector X, sunindextype local_length,
                                  int nvecs, int ntests);
  int Test_N_VWrmsNormMaskVectorArray(N_Vector X, sunindextype local_length,
                                      int nvecs, int ntests);
  int Test_N_VScaleAddMultiVectorArray(N_Vector X, sunindextype local_length,
                                       int nvecs, int nsums, int ntests);
  int Test_N_VLinearCombinationVectorArray(N_Vector X, sunindextype local_length,
                                           int nvecs, int nsums, int tests);
  /* Turn timing on/off */
  void SetTiming(int onoff, int myid);

  /* Print output table header */
  void PrintTableHeader(int type);

  /* Fill realtype arrays with random data */
  void rand_realtype(realtype *data, sunindextype len, realtype lower, realtype upper);
  void rand_realtype_zero_one(realtype *data, sunindextype len);
  void rand_realtype_constraints(realtype *data, sunindextype len);

#ifdef __cplusplus
}
#endif
