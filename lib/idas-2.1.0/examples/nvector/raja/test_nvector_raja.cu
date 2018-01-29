/* ----------------------------------------------------------------- 
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the testing routine to check the NVECTOR Raja module
 * implementation. 
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <nvector/nvector_raja.h>
#include <sundials/sundials_math.h>
#include "test_nvector.h"

#include <nvector/raja/Vector.hpp>

/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int          fails = 0;   /* counter for test failures  */
  sunindextype veclen;      /* vector length              */
  N_Vector     W, X, Y, Z;  /* test vectors               */
  int print_timing;
  /*  sunindextype liw, lrw; */


  /* check input and set vector length */
  if (argc < 3){
    printf("ERROR: TWO (2) Inputs required: vector length, print timing \n");
    return(-1);
  }

  veclen = atol(argv[1]); 
  if (veclen <= 0) {
    printf("ERROR: length of vector must be a positive integer \n");
    return(-1); 
  }

  print_timing = atoi(argv[2]);
  SetTiming(print_timing);


  printf("\nRunning with vector length %ld \n\n", (long) veclen);

  /* Create vectors */
  W = N_VNewEmpty_Raja(veclen);
  X = N_VNew_Raja(veclen);

  /* NVector Tests */
  
  /* Raja specific tests */
  
  /* Memory allocation tests */
  fails += Test_N_VCloneEmpty(X, 0);
  fails += Test_N_VClone(X, veclen, 0);
  fails += Test_N_VCloneEmptyVectorArray(5, X, 0);
  fails += Test_N_VCloneVectorArray(5, X, veclen, 0);

  Y = N_VClone_Raja(X);
  Z = N_VClone_Raja(X);

  /* Skipped tests */
  /*   fails += Test_N_VSetArrayPointer(W, veclen, 0); */
  /*   fails += Test_N_VGetArrayPointer(X, veclen, 0); */
  
  /* Vector operation tests */
  fails += Test_N_VConst(X, veclen, 0);
  fails += Test_N_VLinearSum(X, Y, Z, veclen, 0);
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

  /*  N_VSpace_Raja(X, &lrw, &liw);               */
  /*  printf("lrw = %ld, liw = %ld\n", lrw, liw); */
  
  /* Free vectors */
  N_VDestroy_Raja(W);
  N_VDestroy_Raja(X);
  N_VDestroy_Raja(Y);
  N_VDestroy_Raja(Z);

  /* Print result */
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
  sunrajavec::Vector<double, sunindextype>* xv = sunrajavec::extract<realtype, sunindextype>(X);
  realtype *xdata;
  
  xv->copyFromDev();
  
  xdata = xv->host();
  /* check vector data */
  for(i=0; i < local_length; i++){
    failure += FNEQ(xdata[i], ans);
    //printf("%g ?= %g\n", xdata[i], ans);
  }
  return (failure > ZERO) ? (1) : (0);
}

booleantype has_data(N_Vector X)
{
  sunrajavec::Vector<double, sunindextype>* xv = sunrajavec::extract<realtype, sunindextype>(X);

  return (xv == NULL ? SUNFALSE : SUNTRUE);
}

void set_element(N_Vector X, sunindextype i, realtype val)
{
  sunrajavec::Vector<double, sunindextype>* xv = sunrajavec::extract<realtype, sunindextype>(X);
  xv->copyFromDev();
  (xv->host())[i] = val;
  xv->copyToDev();
}

realtype get_element(N_Vector X, sunindextype i)
{
  sunrajavec::Vector<double, sunindextype>* xv = sunrajavec::extract<realtype, sunindextype>(X);
  xv->copyFromDev();
  return (xv->host())[i];
}
