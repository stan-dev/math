/* ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner @ LLNL
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
 * This is the testing routine to check the OpenMP NVECTOR module
 * implementation. 
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_openmp.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"


/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int          fails = 0;     /* counter for test failures */
  sunindextype veclen;        /* vector length             */
  N_Vector     W, X, Y, Z;    /* test vectors              */
  int          num_threads;
  int          print_timing;

  /* check input and set vector length */
  if (argc < 4){
    printf("ERROR: THREE (3) arguments required: <vector length> <number of threads> <print timing>\n");
    return(-1);
  }

  veclen = atol(argv[1]); 
  if (veclen <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    return(-1); 
  }

  num_threads = atoi(argv[2]);
  if (num_threads <= 0) {
    printf("ERROR: numbber of threads must be a positive integer \n");
    return(-1); 
  }

  print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  printf("\nRunning with vector length %ld \n \n", (long int) veclen);
  printf("\nRunning with number of threads %d \n \n", num_threads);

  /* Create vectors */
  W = N_VNewEmpty_OpenMP(veclen, num_threads);
  X = N_VNew_OpenMP(veclen, num_threads);
  Y = N_VNew_OpenMP(veclen, num_threads);
  Z = N_VNew_OpenMP(veclen, num_threads);

  if(N_VGetVectorID(W) == SUNDIALS_NVEC_OPENMP) {
    /*printf("Testing OpenMP variant of N_Vector...\n");*/
  }
  
  /* NVector Tests */
  fails += Test_N_VSetArrayPointer(W, veclen, 0);
  fails += Test_N_VGetArrayPointer(X, veclen, 0);
  fails += Test_N_VLinearSum(X, Y, Z, veclen, 0);
  fails += Test_N_VConst(X, veclen, 0);
  fails += Test_N_VProd(X, Y, Z, veclen, 0);
  fails += Test_N_VDiv(X, Y, Z, veclen, 0);
  fails += Test_N_VScale(X, Z, veclen, 0);
  fails += Test_N_VAbs(X, Z, veclen, 0);
  fails += Test_N_VInv(X, Z, veclen, 0);
  fails += Test_N_VAddConst(X, Z, veclen, 0);
  fails += Test_N_VDotProd(X, Y, veclen, veclen, 0);
  fails += Test_N_VMaxNorm(X, veclen, 0);
  fails += Test_N_VWrmsNorm(X, Y, veclen, 0);
  fails += Test_N_VWrmsNormMask(X, Y, Z, veclen, veclen, 0);
  fails += Test_N_VMin(X, veclen, 0);
  fails += Test_N_VWL2Norm(X, Y, veclen, veclen, 0);
  fails += Test_N_VL1Norm(X, veclen, veclen, 0);
  fails += Test_N_VCompare(X, Z, veclen, 0);
  fails += Test_N_VInvTest(X, Z, veclen, 0);
  fails += Test_N_VConstrMask(X, Y, Z, veclen, 0);
  fails += Test_N_VMinQuotient(X, Y, veclen, 0);
  fails += Test_N_VCloneVectorArray(5, X, veclen, 0);
  fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
  fails += Test_N_VCloneEmpty(X, 0);
  fails += Test_N_VClone(X, veclen, 0);

  /* Free vectors */
  N_VDestroy_OpenMP(W);
  N_VDestroy_OpenMP(X);
  N_VDestroy_OpenMP(Y);
  N_VDestroy_OpenMP(Z);

  /* Print results */
  if (fails) {
    printf("FAIL: NVector module failed %i tests \n \n", fails);
  } else {
    printf("SUCCESS: NVector module passed all tests \n \n");
  }

  return(fails);
}

/* ----------------------------------------------------------------------
 * Check vector
 * --------------------------------------------------------------------*/
int check_ans(realtype ans, N_Vector X, sunindextype local_length)
{
  int      failure = 0;
  sunindextype i;
  realtype *Xdata;
  
  Xdata = N_VGetArrayPointer(X);

  /* check vector data */
  for(i=0; i < local_length; i++){
    failure += FNEQ(Xdata[i], ans);
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(N_Vector X)
{
  realtype *Xdata = N_VGetArrayPointer(X);
  if (Xdata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  NV_Ith_OMP(X,i) = val;
}
 
realtype get_element(N_Vector X, sunindextype i)
{
  return NV_Ith_OMP(X,i);
}

