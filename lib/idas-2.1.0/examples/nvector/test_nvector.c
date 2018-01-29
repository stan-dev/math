/* ----------------------------------------------------------------- 
 * Programmer(s): David J. Gardner and Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Acknowledgements: These testing routines are based on an
 *                   NVECTOR testing routine by Daniel R. Reynolds
 *                   @ SMU.
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
 * These test functions are designed to check an NVECTOR module 
 * implementation. 
 *
 * NOTE: Many of these tests rely on the N_VGetArrayPointer routine 
 *       to get a pointer to the data component of an N_Vector. This 
 *       assumes the internal data is stored in a contiguous 
 *       realtype array.
 * -----------------------------------------------------------------*/

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <math.h> /* include isnan */
#include <stdio.h>
#include <stdlib.h>

#include "test_nvector.h"

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif


/* private functions */
static double get_time();

int print_time = 0;

#define PRINT_TIME(format, time) if(print_time) printf(format, time)

/* ----------------------------------------------------------------------
 * N_VCloneVectorArray Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VCloneVectorArray(int count, N_Vector W, sunindextype local_length, int myid)
{
  int      i, failure;
  double   start_time, stop_time;
  N_Vector *vs;
  
  /* clone array of vectors */
  start_time = get_time(); 
  vs = N_VCloneVectorArray(count, W);
  stop_time = get_time(); 
  
  /* check array of vectors */
  if (count <= 0 && vs != NULL) {
    printf(">>> FAILED test -- N_VCloneVectorArray, Proc %d \n", myid);
    printf("    count = %d, expected *vs = NULL \n \n",count);
    return(1);
  } 
  
  /* check vectors in array */
  for(i=0; i<count; i++) {
    if (vs[i] == NULL) {
      printf(">>> FAILED test -- N_VCloneVectorArray, Proc %d \n", myid);
      printf("    Vector[%d] = NULL \n \n",i);
      N_VDestroyVectorArray(vs, count);
      return(1);
    }    
    
    N_VConst(ONE,vs[i]);
    failure = check_ans(ONE, vs[i], local_length);
    if (failure) {
      printf(">>> FAILED test -- N_VCloneVectorArray, Proc %d \n", myid);
      printf("    Vector[%d] failed N_VConst check \n \n",i);
      N_VDestroyVectorArray(vs, count);
      return(1);
    }    
  }

  N_VDestroyVectorArray(vs, count);
  
  if (myid == 0) {
    printf("    PASSED test -- N_VCloneVectorArray \n");
    PRINT_TIME("    N_VCloneVectorArray Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}

/* ----------------------------------------------------------------------
 * N_VCloneVectorArrayEmpty Test
 * --------------------------------------------------------------------*/
int Test_N_VCloneEmptyVectorArray(int count, N_Vector W, int myid)
{
  int      i;
  double   start_time, stop_time;
  N_Vector *vs;

  /* clone empty array */
  start_time = get_time(); 
  vs = N_VCloneEmptyVectorArray(count, W);
  stop_time = get_time(); 
  
  /* check array of vectors */
  if (count <= 0 && vs != NULL) {
    printf(">>> FAILED test -- N_VCloneEmptyVectorArray, Proc %d \n", myid);
    printf("    count = %d, expected *vs = NULL \n \n",count);
    return(1);
  } 

  /* check vectors in array */
  for(i=0; i<count; i++) {
    if (vs[i] == NULL) {
      printf(">>> FAILED test -- N_VCloneEmptyVectorArray, Proc %d \n", myid);
      printf("    Vector[%d] = NULL \n \n",i);
      N_VDestroyVectorArray(vs, count);
      return(1);
    }    

    if (has_data(vs[i])) {
      printf(">>> FAILED test -- N_VCloneEmptyVectorArray, Proc %d \n", myid);
      printf("    Vector[%d] data != NULL \n \n",i);
      N_VDestroyVectorArray(vs, count);
      return(1);
    }    
  }

  N_VDestroyVectorArray(vs, count);
  
  if (myid == 0) {
    printf("    PASSED test -- N_VCloneEmptyVectorArray \n");
    PRINT_TIME("    N_VCloneEmptyVectorArray Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VCloneEmpty Test
 * --------------------------------------------------------------------*/
int Test_N_VCloneEmpty(N_Vector W, int myid)
{
  double   start_time, stop_time;
  N_Vector X;

  /* clone empty vector */
  start_time = get_time();   
  X = N_VCloneEmpty(W);
  stop_time = get_time(); 

  /* check vector */
  if (X == NULL) {
    printf(">>> FAILED test -- N_VCloneEmpty, Proc %d \n", myid);
    printf("    After N_VCloneEmpty, X == NULL \n \n");
    return(1);
  } 

  /* check vector data */
  if (has_data(X)) {
    printf(">>> FAILED test -- N_VCloneEmpty, Proc %d \n", myid);
    printf("    Vector data != NULL \n \n");
    N_VDestroy(X);
    return(1);
  }    

  N_VDestroy(X); 

  if (myid == 0) {
    printf("    PASSED test -- N_VCloneEmpty \n");
    PRINT_TIME("    N_VCloneEmpty Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VClone Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VClone(N_Vector W, sunindextype local_length, int myid)
{
  int      failure;
  double   start_time, stop_time;
  N_Vector X;

  /* clone vector */
  start_time = get_time();   
  X = N_VClone(W);
  stop_time = get_time();   

  /* check cloned vector */
  if (X == NULL) {
    printf(">>> FAILED test -- N_VClone, Proc %d \n", myid);
    printf("    After N_VClone, X == NULL \n \n");
    return(1);
  } 

  /* check cloned vector data */
  if (!has_data(X)) {
    printf(">>> FAILED test -- N_VClone, Proc %d \n", myid);
    printf("    Vector data == NULL \n \n");
    N_VDestroy(X);
    return(1);
  }    

  N_VConst(ONE,X);
  failure = check_ans(ONE, X, local_length);
  if (failure) {
    printf(">>> FAILED test -- N_VClone, Proc %d \n", myid);
    printf("    Failed N_VClone check \n \n");
    N_VDestroy(X);
    return(1);
  }    

  N_VDestroy(X); 

  if (myid == 0) {
    printf("    PASSED test -- N_VClone \n");
    PRINT_TIME("    N_VClone Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * N_VGetArrayPointer Test
 * 
 * For now commenting this out to surpress warning messages (pointer set,
 * but not used). Do we really need to time access to the vector 
 * data pointer? 
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VGetArrayPointer(N_Vector W, sunindextype local_length, int myid)
{
  int      failure = 0;
  double   start_time, stop_time;
  realtype *Wdata;

  /* get vector data, time it and set it to NULL */
  start_time = get_time();   
  Wdata = N_VGetArrayPointer(W);
  stop_time = get_time();
  Wdata++; Wdata=NULL; /* Do something with pointer to surpress warning */
  
  /* check vector data */
  if (!has_data(W)) {
    printf(">>> FAILED test -- N_VGetArrayPointer, Proc %d \n", myid);
    printf("    Vector data == NULL \n \n");
    return(1);
  }    

  N_VConst(NEG_HALF,W);
  failure = check_ans(NEG_HALF, W, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VGetArrayPointer, Proc %d \n", myid);
    printf("    Failed N_VConst check \n \n");
    return(1);
  }

  if (myid == 0) {
    printf("    PASSED test -- N_VGetArrayPointer \n");
    PRINT_TIME("    N_VGetArrayPointer Time: %22.15e \n \n", stop_time - start_time);
  }
  
  return(0);
}


/* ----------------------------------------------------------------------
 * N_VSetArrayPointer Test
 *
 * NOTE: This routine depends on N_VConst to check vector data.
 * --------------------------------------------------------------------*/
int Test_N_VSetArrayPointer(N_Vector W, sunindextype local_length, int myid)
{
  int      failure = 0;
  sunindextype i;  
  double   start_time, stop_time;
  realtype *Wdata;

  /* create vector data */
  Wdata = malloc(local_length * sizeof(realtype));
  for(i=0; i < local_length; i++){
    Wdata[i] = ONE;
  }
  
  /* attach data to vector */
  start_time = get_time();   
  N_VSetArrayPointer(Wdata, W);
  stop_time = get_time();   

  /* check vector data */
  N_VConst(NEG_HALF,W);
  for(i=0; i < local_length; i++){
    failure += FNEQ(Wdata[i], NEG_HALF);
  }

  if (failure) {
    printf(">>> FAILED test -- N_VSetArrayPointer, Proc %d \n", myid);
    printf("    Failed N_VConst check \n \n");
    free(Wdata);
    return(1);
  }

  free(Wdata);

  if (myid == 0) {
    printf("    PASSED test -- N_VSetArrayPointer \n");
    PRINT_TIME("    N_VSetArrayPointer Time: %22.15e \n \n", stop_time - start_time);
  }
  
  return(0);
}


/* ----------------------------------------------------------------------
 * N_VLinearSum Tests
 * --------------------------------------------------------------------*/
int Test_N_VLinearSum(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  
  /* 
   * Case 1a: y = x + y, (Vaxpy Case 1) 
   */

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(NEG_TWO, Y);

  start_time = get_time(); 
  N_VLinearSum(ONE, X, ONE, Y, Y);
  stop_time = get_time(); 
  
  /* Y should be vector of -1 */
  failure = check_ans(NEG_ONE, Y, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 1a, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 1a \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /*
   * Case 1b: y = -x + y, (Vaxpy Case 2)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(TWO, Y);

  start_time = get_time(); 
  N_VLinearSum(NEG_ONE, X, ONE, Y, Y);
  stop_time = get_time(); 
  
  /* Y should be vector of +1 */
  failure = check_ans(ONE, Y, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 1b, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 1b \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /*
   * Case 1c: y = ax + y, (Vaxpy Case 3)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_TWO, Y);

  start_time = get_time(); 
  N_VLinearSum(HALF, X, ONE, Y, Y);
  stop_time = get_time(); 
  
  /* Y should be vector of -1 */
  failure = check_ans(NEG_ONE, Y, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 1c, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 1c \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /*
   * Case 2a: x = x + y, (Vaxpy Case 1)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_ONE, Y);

  start_time = get_time(); 
  N_VLinearSum(ONE, X, ONE, Y, X);
  stop_time = get_time(); 
  
  /* Y should be vector of +1 */
  failure = check_ans(ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 2a, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 2a \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /*
   * Case 2b: x = x - y, (Vaxpy Case 2)
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(TWO, Y);

  start_time = get_time(); 
  N_VLinearSum(ONE, X, NEG_ONE, Y, X);
  stop_time = get_time(); 
  
  /* Y should be vector of -1 */
  failure = check_ans(NEG_ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 2b, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 2b \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /* 
   * Case 2c: x = x + by, (Vaxpy Case 3) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_HALF, Y);

  start_time = get_time(); 
  N_VLinearSum(ONE, X, TWO, Y, X);
  stop_time = get_time(); 
  
  /* X should be vector of +1 */
  failure = check_ans(ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 2c, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 2c \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /* 
   * Case 3: z = x + y, (VSum) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(ONE, X, ONE, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 3, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 3 \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /* 
   * Case 4a: z = x - y, (VDiff) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(ONE, X, NEG_ONE, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 4a, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 4a \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /* 
   * Case 4b: z = -x + y, (VDiff) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(NEG_ONE, X, ONE, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 4b, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 4b \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /* 
   * Case 5a: z = x + by, (VLin1) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(ONE, X, TWO, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 5a, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 5a \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 5b: z = ax + y, (VLin1) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(HALF, X);
  N_VConst(NEG_TWO, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(TWO, X, ONE, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 5b, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 5b \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 6a: z = -x + by, (VLin2) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(NEG_ONE, X, TWO, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 6a, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 6a \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 6b: z = ax - y, (VLin2) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(HALF, X);
  N_VConst(TWO, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(TWO, X, NEG_ONE, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 6b, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 6b \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 7: z = a(x + y), (VScaleSum) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(TWO, X, TWO, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 7, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 7 \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 8: z = a(x - y), (VScaleDiff) 
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(HALF, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(TWO, X, NEG_TWO, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 8, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 8 \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 9: z = ax + by, All Other Cases
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(NEG_TWO, Y);
  N_VConst(ZERO, Z);

  start_time = get_time(); 
  N_VLinearSum(TWO, X, HALF, Y, Z);
  stop_time = get_time(); 
  
  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VLinearSum Case 9, Proc %d \n", myid);
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VLinearSum Case 9 \n");
    PRINT_TIME("    N_VLinearSum Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VConst Test
 * --------------------------------------------------------------------*/
int Test_N_VConst(N_Vector X, sunindextype local_length, int myid)
{
  int      i, fails = 0, failure = 0;
  double   start_time, stop_time;

  /* fill vector data with zeros to prevent passing in the case where
     the input vector is a vector of ones */
  for(i=0; i < local_length; i++){
    set_element(X, i, ZERO);
  }

  start_time = get_time();
  N_VConst(ONE,X);
  stop_time = get_time(); 

  /* X should be vector of +1 */
  failure = check_ans(ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VConst, Proc %d \n", myid);
    PRINT_TIME("    N_VConst Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VConst \n");
    PRINT_TIME("    N_VConst Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VProd Test
 * --------------------------------------------------------------------*/
int Test_N_VProd(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VProd(X, Y, Z);
  stop_time = get_time(); 

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VProd, Proc %d \n", myid);
    PRINT_TIME("    N_VProd Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VProd \n");
    PRINT_TIME("    N_VProd Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VDiv Test
 * --------------------------------------------------------------------*/
int Test_N_VDiv(N_Vector X, N_Vector Y, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(TWO, Y);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VDiv(X, Y, Z);
  stop_time = get_time(); 

  /* Z should be vector of +1/2 */
  failure = check_ans(HALF, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VDiv, Proc %d \n", myid);
    PRINT_TIME("    N_VDiv Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VDiv \n");
    PRINT_TIME("    N_VDiv Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VScale Tests
 * --------------------------------------------------------------------*/
int Test_N_VScale(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;

  /* 
   * Case 1: x = cx, VScaleBy
   */

  /* fill vector data */
  N_VConst(HALF, X);

  start_time = get_time();
  N_VScale(TWO, X, X);
  stop_time = get_time(); 

  /* X should be vector of +1 */
  failure = check_ans(ONE, X, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VScale Case 1, Proc %d \n", myid);
    PRINT_TIME("    N_VScale Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VScale Case 1 \n");
    PRINT_TIME("    N_VScale Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 2: z = x, VCopy
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_ONE, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VScale(ONE, X, Z);
  stop_time = get_time(); 

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VScale Case 2, Proc %d \n", myid);
    PRINT_TIME("    N_VScale Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VScale Case 2 \n");
    PRINT_TIME("    N_VScale Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 3: z = -x, VNeg
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_ONE, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VScale(NEG_ONE, X, Z);
  stop_time = get_time(); 

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VScale Case 3, Proc %d \n", myid);
    PRINT_TIME("    N_VScale Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VScale Case 3 \n");
    PRINT_TIME("    N_VScale Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 4: z = cx, All other cases
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VScale(TWO, X, Z);
  stop_time = get_time(); 

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VScale Case 4, Proc %d \n", myid);
    PRINT_TIME("    N_VScale Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VScale Case 4 \n");
    PRINT_TIME("    N_VScale Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VAbs Test
 * --------------------------------------------------------------------*/
int Test_N_VAbs(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;

  /* fill vector data */
  N_VConst(NEG_ONE, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VAbs(X,Z);
  stop_time = get_time(); 

  /* Z should be vector of +1 */
  failure = check_ans(ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VAbs, Proc %d \n", myid);
    PRINT_TIME("    N_VAbs Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VAbs \n");
    PRINT_TIME("    N_VAbs Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VInv Test
 * --------------------------------------------------------------------*/
int Test_N_VInv(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VInv(X,Z);
  stop_time = get_time(); 

  /* Z should be vector of +1/2 */
  failure = check_ans(HALF, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VInv, Proc %d \n", myid);
    PRINT_TIME("    N_VInv Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VInv \n");
    PRINT_TIME("    N_VInv Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VAddConst Test
 * --------------------------------------------------------------------*/
int Test_N_VAddConst(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;

  /* fill vector data */
  N_VConst(ONE, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  N_VAddConst(X,NEG_TWO,Z);
  stop_time = get_time(); 

  /* Z should be vector of -1 */
  failure = check_ans(NEG_ONE, Z, local_length);

  if (failure) {
    printf(">>> FAILED test -- N_VAddConst, Proc %d \n", myid);
    PRINT_TIME("    N_VAddConst Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VAddConst \n");
    PRINT_TIME("    N_VAddConst Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VDotProd Test
 * --------------------------------------------------------------------*/
int Test_N_VDotProd(N_Vector X, N_Vector Y, 
                    sunindextype local_length, sunindextype global_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  realtype ans;

  /* fill vector data */
  N_VConst(TWO, X);
  N_VConst(HALF, Y);

  start_time = get_time();
  ans = N_VDotProd(X,Y);
  stop_time = get_time(); 

  /* ans should equal global vector length */
  failure = FNEQ(ans, global_length);

  if (failure) {
    printf(">>> FAILED test -- N_VDotProd, Proc %d \n", myid);
    PRINT_TIME("    N_VDotProd Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VDotProd \n");
    PRINT_TIME("    N_VDotProd Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMaxNorm Test
 * --------------------------------------------------------------------*/
int Test_N_VMaxNorm(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  realtype ans;

  /* fill vector data */
  N_VConst(NEG_ONE, X);
  set_element(X, local_length-1, NEG_TWO);
  
  start_time = get_time();
  ans = N_VMaxNorm(X);
  stop_time = get_time(); 

  /* ans should equal 2 */
  failure = (ans < ZERO) ? 1 : FNEQ(ans, TWO);

  if (failure) {
    printf(">>> FAILED test -- N_VMaxNorm, Proc %d \n", myid);
    PRINT_TIME("    N_VMaxNorm Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VMaxNorm \n");
    PRINT_TIME("    N_VMaxNorm Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNorm Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNorm(N_Vector X, N_Vector W, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  realtype ans;

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(HALF, W);

  start_time = get_time();
  ans = N_VWrmsNorm(X, W);
  stop_time = get_time(); 

  /* ans should equal 1/4 */
  failure = (ans < ZERO) ? 1 : FNEQ(ans, HALF*HALF);

  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNorm, Proc %d \n", myid);
    PRINT_TIME("    N_VWrmsNorm Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VWrmsNorm \n");
    PRINT_TIME("    N_VWrmsNorm Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWrmsNormMask Test
 * --------------------------------------------------------------------*/
int Test_N_VWrmsNormMask(N_Vector X, N_Vector W, N_Vector ID, 
			 sunindextype local_length, sunindextype global_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  realtype ans;

  /* 
   * Case 1: use all elements, ID = 1
   */

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(HALF, W);
  N_VConst(ONE, ID);

  start_time = get_time();
  ans = N_VWrmsNormMask(X, W, ID);
  stop_time = get_time(); 
  
  /* ans equals 1/4 (same as wrms norm) */
  failure = (ans < ZERO) ? 1 : FNEQ(ans, HALF*HALF);
    
  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNormMask Case 1, Proc %d \n", myid);
    PRINT_TIME("    N_VWrmsNormMask Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VWrmsNormMask Case 1 \n");
    PRINT_TIME("    N_VWrmsNormMask Time: %22.15e \n \n", stop_time - start_time);
  }    

  /* 
   * Case 2: use no elements, ID = 0
   */
  
  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(HALF, W);
  N_VConst(ZERO, ID);

  start_time = get_time();
  ans = N_VWrmsNormMask(X, W, ID);
  stop_time = get_time(); 
  
  /* ans equals 0 (skips all elements) */
  failure = (ans < ZERO) ? 1 : FNEQ(ans, ZERO);
    
  if (failure) {
    printf(">>> FAILED test -- N_VWrmsNormMask Case 2, Proc %d \n", myid);
    PRINT_TIME("    N_VWrmsNormMask Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VWrmsNormMask Case 2 \n");
    PRINT_TIME("    N_VWrmsNormMask Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMin Test
 * --------------------------------------------------------------------*/
int Test_N_VMin(N_Vector X, sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  realtype ans;

  /* fill vector data */
  N_VConst(TWO, X);
  set_element(X, local_length-1, NEG_ONE);

  start_time = get_time();
  ans = N_VMin(X);
  stop_time = get_time(); 

  /* ans should equal -1 */
  failure = FNEQ(ans, NEG_ONE);

  if (failure) {
    printf(">>> FAILED test -- N_VMin, Proc %d \n", myid);
    PRINT_TIME("    N_VMin Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VMin \n");
    PRINT_TIME("    N_VMin Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VWL2Norm Test
 * --------------------------------------------------------------------*/
int Test_N_VWL2Norm(N_Vector X, N_Vector W, 
                    sunindextype local_length, sunindextype global_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  realtype ans;

  /* fill vector data */
  N_VConst(NEG_HALF, X);
  N_VConst(HALF, W);

  start_time = get_time();
  ans = N_VWL2Norm(X, W);
  stop_time = get_time(); 

  /* ans should equal 1/4 * sqrt(global_length) */
  failure = (ans < ZERO) ? 1 : FNEQ(ans, HALF*HALF*SUNRsqrt((realtype) global_length));

  if (failure) {
    printf(">>> FAILED test -- N_VWL2Norm, Proc %d \n", myid);
    PRINT_TIME("    N_VWL2Norm Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VWL2Norm \n");
    PRINT_TIME("    N_VWL2Norm Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VL1Norm Test
 * --------------------------------------------------------------------*/
int Test_N_VL1Norm(N_Vector X, sunindextype local_length, 
                   sunindextype global_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  realtype ans;

  /* fill vector data */
  N_VConst(NEG_ONE, X);

  start_time = get_time();
  ans = N_VL1Norm(X);
  stop_time = get_time(); 

  /* ans should equal global_length */
  failure = (ans < ZERO) ? 1 : FNEQ(ans, global_length);

  if (failure) {
    printf(">>> FAILED test -- N_VL1Norm, Proc %d \n", myid);
    PRINT_TIME("    N_VL1Norm Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VL1Norm \n");
    PRINT_TIME("    N_VL1Norm Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VCompare
 * --------------------------------------------------------------------*/
int Test_N_VCompare(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int      mask, fails = 0, failure = 0;
  double   start_time, stop_time;
  sunindextype i;

  if (local_length < 3) {
    printf("Error Test_N_VCompare: Local vector length is %ld, length must be >= 3",
           (long int) local_length);
    return(-1);
  }

  /* fill vector data */
  for(i=0; i < local_length; i++){
    set_element(Z, i, NEG_ONE);

    mask = i % 3;
    switch(mask) {

    case 0 :
      /* abs(X[i]) < c */
      set_element(X, i, ZERO);
      break;

    case 1 :
      /* abs(X[i]) = c */
      set_element(X, i, NEG_ONE);
      break;
      
    case 2 :
      /* abs(X[i]) > c */
      set_element(X, i, NEG_TWO);
      break;
    }
  }

  start_time = get_time();
  N_VCompare(ONE, X, Z);
  stop_time = get_time(); 

  /* check return vector */
  for(i=0; i < local_length; i++){
    mask = i % 3;

    switch(mask) {

    case 0 :
      /* Z[i] == 0 */
      if (get_element(Z, i) != ZERO)
        failure = 1;
      break;

    case 1 :
      /* Z[i] == 1 */
      if (get_element(Z, i) != ONE)
        failure = 1;
      break;
      
    case 2 :
      /* Z[i] == 1 */
      if (get_element(Z, i) != ONE)
        failure = 1;
      break;
    }
  }

  if (failure) {
    printf(">>> FAILED test -- N_VCompare, Proc %d \n", myid);
    PRINT_TIME("    N_VCompare Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VCompare \n");
    PRINT_TIME("    N_VCompare Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VInvTest
 * --------------------------------------------------------------------*/
int Test_N_VInvTest(N_Vector X, N_Vector Z, sunindextype local_length, int myid)
{
  int         mask, fails = 0, failure = 0;
  double      start_time, stop_time;
  sunindextype    i;
  booleantype test;

  if (local_length < 2) {
    printf("Error Test_N_VCompare: Local vector length is %ld, length must be >= 2",
           (long int) local_length);
    return(-1);
  }

  /*
   * Case 1: All elements Nonzero, z[i] = 1/x[i], return True
   */

  /* fill vector data */
  N_VConst(HALF, X);
  N_VConst(ZERO, Z);

  start_time = get_time();
  test = N_VInvTest(X, Z);
  stop_time = get_time(); 

  /* Z should be vector of +2 */
  failure = check_ans(TWO, Z, local_length);

  if (failure || !test) {
    printf(">>> FAILED test -- N_VInvTest Case 1, Proc %d \n", myid);
    PRINT_TIME("    N_VInvTest Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VInvTest Case 1 \n");
    PRINT_TIME("    N_VInvTest Time: %22.15e \n \n", stop_time - start_time);
  }    

  /*
   * Case 2: Some elements Zero, z[i] = 1/x[i] for x[i] != 0, return False
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(ZERO, Z);
  for(i=0; i < local_length; i++){
    mask = i % 2;   
    if (mask)
      set_element(X, i, HALF);
    else
      set_element(X, i, ZERO);
  }

  start_time = get_time();
  test = N_VInvTest(X, Z);
  stop_time = get_time();

  /* check return vector */
  for(i=0; i < local_length; i++){
    mask = i % 2;

    if (mask) {
      if (get_element(Z, i) != TWO) 
        failure = 1;
    } else {
      if (get_element(Z, i) != ZERO) 
        failure = 1;
    }
  }

  if (failure || test) {
    printf(">>> FAILED test -- N_VInvTest Case 2, Proc %d \n", myid);
    PRINT_TIME("    N_VInvTest Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VInvTest Case 2 \n");
    PRINT_TIME("    N_VInvTest Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VConstrMask
 * --------------------------------------------------------------------*/
int Test_N_VConstrMask(N_Vector C, N_Vector X, N_Vector M, 
                       sunindextype local_length, int myid)
{
  int         mask, fails = 0, failure = 0;
  double      start_time, stop_time;
  sunindextype    i;
  booleantype test;

  if (local_length < 7) {
    printf("Error Test_N_VCompare: Local vector length is %ld, length must be >= 7",
           (long int) local_length);
    return(-1);
  }

  /*
   * Case 1: Return True
   */

  /* fill vector data */
  for(i=0; i < local_length; i++){
    set_element(M, i, NEG_ONE);

    mask = i % 7;  
    switch(mask) {
    case 0 :
      /* c = -2, test for < 0*/
      set_element(C, i, NEG_TWO);
      set_element(X, i, NEG_TWO);
      break;
      
    case 1 :
      /* c = -1, test for <= 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, NEG_ONE);	
      break;
      
    case 2 :
      /* c = -1, test for == 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, ZERO); 
      break;
      
    case 3 :
      /* c = 0, no test */
      set_element(C, i, ZERO);
      set_element(X, i, HALF);
      break;
      
    case 4 :
      /* c = 1, test for == 0*/
      set_element(C, i, ONE);
      set_element(X, i, ZERO);
      break;
      
    case 5 :
      /* c = 1, test for >= 0*/
      set_element(C, i, ONE);
      set_element(X, i, ONE);
      break;
      
    case 6:
      /* c = 2, test for > 0 */
      set_element(C, i, TWO);
      set_element(X, i, TWO);
      break;
    }
  }

  start_time = get_time(); 
  test = N_VConstrMask(C, X, M);
  stop_time = get_time();

  /* M should be vector of 0 */
  failure = check_ans(ZERO, M, local_length);

  if (failure || !test) {
    printf(">>> FAILED test -- N_VConstrMask Case 1, Proc %d \n", myid);
    PRINT_TIME("    N_VConstrMask Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VConstrMask Case 1 \n");
    PRINT_TIME("    N_VConstrMask Time: %22.15e \n \n", stop_time - start_time);
  }    

  /*
   * Case 2: Return False
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  for(i=0; i < local_length; i++){
    set_element(M, i, NEG_ONE);
    
    mask = i % 5;  
    switch(mask) {
    case 0 :
      /* c = -2, test for < 0*/
      set_element(C, i, NEG_TWO);
      set_element(X, i, TWO);
      break;

    case 1 :
      /* c = -1, test for <= 0 */
      set_element(C, i, NEG_ONE);
      set_element(X, i, ONE);	
      break;
      
    case 2 :
      /* c = 0, no test */
      set_element(C, i, ZERO);
      set_element(X, i, HALF);
      break;

    case 3 :
      /* c = 1, test for >= 0*/
      set_element(C, i, ONE);
      set_element(X, i, NEG_ONE);
      break;

    case 4 :
      /* c = 2, test for > 0 */
      set_element(C, i, TWO);
      set_element(X, i, NEG_TWO);
      break;
    }
  }

  start_time = get_time();  
  test = N_VConstrMask(C, X, M);
  stop_time = get_time();

  /* check mask vector */
  for(i=0; i < local_length; i++){
    mask = i % 5;
    
    if (mask == 2){
      if (get_element(M, i) != ZERO) 
        failure = 1;
    } else {
      if (get_element(M, i) != ONE)
        failure = 1;
    }
  }
  
  if (failure || test) {
    printf(">>> FAILED test -- N_VConstrMask Case 2, Proc %d \n", myid);
    PRINT_TIME("    N_VConstrMask Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VConstrMask Case 2 \n");
    PRINT_TIME("    N_VConstrMask Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ----------------------------------------------------------------------
 * N_VMinQuotient Test
 * --------------------------------------------------------------------*/
int Test_N_VMinQuotient(N_Vector NUM, N_Vector DENOM, 
                        sunindextype local_length, int myid)
{
  int      fails = 0, failure = 0;
  double   start_time, stop_time;
  realtype ans;

  /*
   * Case 1: Pass
   */

  /* fill vector data */
  N_VConst(TWO, NUM);
  N_VConst(TWO, DENOM);
  set_element(NUM, local_length-1, ONE);

  start_time = get_time();
  ans = N_VMinQuotient(NUM, DENOM);
  stop_time = get_time(); 

  /* ans should equal 1/2 */
  failure = FNEQ(ans, HALF);
  
  if (failure) {
    printf(">>> FAILED test -- N_VMinQuotient Case 1, Proc %d \n", myid);
    PRINT_TIME("    N_VMinQuotient Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VMinQuotient Case 1 \n");
    PRINT_TIME("    N_VMinQuotient Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  /*
   * Case 2: Fail
   */

  /* reset failure */
  failure = 0;

  /* fill vector data */
  N_VConst(TWO, NUM);
  N_VConst(ZERO, DENOM);

  start_time = get_time();
  ans = N_VMinQuotient(NUM, DENOM);
  stop_time = get_time(); 
  
  /* ans should equal BIG_REAL */
  failure = FNEQ(ans, BIG_REAL);

  if (failure) {
    printf(">>> FAILED test -- N_VMinQuotient Case 2, Proc %d \n", myid);
    PRINT_TIME("    N_VMinQuotient Time: %22.15e \n \n", stop_time - start_time);
    fails++;
  }
  else if (myid == 0) {
    printf("    PASSED test -- N_VMinQuotient Case 2 \n");
    PRINT_TIME("    N_VMinQuotient Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(fails);
}


/* ======================================================================
 * Private functions
 * ====================================================================*/

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
time_t base_time_tv_sec = 0; /* Base time; makes time values returned
                                by get_time easier to read when
                                printed since they will be zero
                                based.
                              */
#endif

void SetTiming(int onoff)
{
   print_time = onoff;

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;  
  clock_gettime( CLOCK_MONOTONIC_RAW, &spec );
  base_time_tv_sec = spec.tv_sec;
#endif
}

/* ----------------------------------------------------------------------
 * Timer
 * --------------------------------------------------------------------*/
static double get_time()
{
#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
  struct timespec spec;  
  clock_gettime( CLOCK_MONOTONIC_RAW, &spec );
  double time = (double)(spec.tv_sec - base_time_tv_sec) + ((double)(spec.tv_nsec) / 1E9);
#else
  double time = 0;
#endif
  return time;
}


