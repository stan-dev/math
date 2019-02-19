/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                David J. Gardner @ LLNL
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
 * These test functions are designed to check a SUNLinSol module
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include <math.h> /* include isnan */
#include <stdio.h>
#include <stdlib.h>

#include "test_sunlinsol.h"

#if defined( SUNDIALS_HAVE_POSIX_TIMERS) && defined(_POSIX_TIMERS)
#include <time.h>
#include <unistd.h>
#endif


/* private functions */
static double get_time();

int print_time = 0;

#define PRINT_TIME(format, time) if(print_time) printf(format, time)


/* ----------------------------------------------------------------------
 * SUNLinSolGetType Test
 * --------------------------------------------------------------------*/
int Test_SUNLinSolGetType(SUNLinearSolver S, SUNLinearSolver_Type sunid, int myid)
{
  double               start_time, stop_time;
  SUNLinearSolver_Type mysunid;

  start_time = get_time();   
  mysunid = SUNLinSolGetType(S);
  stop_time = get_time();   

  if (sunid != mysunid) {
    printf(">>> FAILED test -- SUNLinSolGetType, Proc %d \n", myid);
    PRINT_TIME("    SUNLinSolGetType Time: %22.15e \n \n", stop_time - start_time);
    return(1);
  } else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolGetType \n");
    PRINT_TIME("    SUNLinSolGetType Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * Test_SUNLinSolLastFlag Test
 * --------------------------------------------------------------------*/
int Test_SUNLinSolLastFlag(SUNLinearSolver S, int myid)
{
  double   start_time, stop_time;
  long int lastflag;

  /* the only way to fail this test is if the function is NULL, 
     which will cause a seg-fault */
  start_time = get_time();   
  lastflag = SUNLinSolLastFlag(S);
  stop_time = get_time();   

  if (myid == 0) {
    printf("    PASSED test -- SUNLinSolLastFlag (%ld) \n", lastflag);
    PRINT_TIME("    SUNLinSolLastFlag Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * Test_SUNLinSolSpace Test
 * --------------------------------------------------------------------*/
int Test_SUNLinSolSpace(SUNLinearSolver S, int myid)
{
  int      failure;
  double   start_time, stop_time;
  long int lenrw, leniw;

  /* call SUNLinSolSpace (failure based on output flag) */
  start_time = get_time();   
  failure = SUNLinSolSpace(S, &lenrw, &leniw);
  stop_time = get_time();   

  if (failure) {
    printf(">>> FAILED test -- SUNLinSolSpace, Proc %d \n", myid);
    PRINT_TIME("    SUNLinSolSpace Time: %22.15e \n \n", stop_time - start_time);
    return(1);
  } else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolSpace, lenrw = %li, leniw = %li\n", lenrw, leniw);
    PRINT_TIME("    SUNLinSolSpace Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolNumIters Test
 * --------------------------------------------------------------------*/
int Test_SUNLinSolNumIters(SUNLinearSolver S, int myid)
{
  int    numiters;
  double start_time, stop_time;

  /* the only way to fail this test is if the function is NULL, 
     which will cause a seg-fault */
  start_time = get_time();   
  numiters = SUNLinSolNumIters(S);
  stop_time = get_time();   

  if (myid == 0) {
    printf("    PASSED test -- SUNLinSolNumIters (%d) \n", numiters);
    PRINT_TIME("    SUNLinSolNumIters Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolResNorm Test
 * --------------------------------------------------------------------*/
int Test_SUNLinSolResNorm(SUNLinearSolver S, int myid)
{
  double start_time, stop_time, resnorm;

  /* this test can fail if the function is NULL, which will cause a seg-fault */
  start_time = get_time();   
  resnorm = (double) SUNLinSolResNorm(S);
  stop_time = get_time();   

  /* this test can also fail if the return value is negative */
  if (resnorm < ZERO){
    printf(">>> FAILED test -- SUNLinSolResNorm returned %g on Proc %d \n", 
           resnorm, myid);
    return(1);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolResNorm\n");
    PRINT_TIME("    SUNLinSolResNorm Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolResid Test
 * --------------------------------------------------------------------*/
int Test_SUNLinSolResid(SUNLinearSolver S, int myid)
{
  double start_time, stop_time;
  N_Vector resid;

  /* this test can fail if the function returns NULL */
  start_time = get_time();
  resid = SUNLinSolResid(S);
  stop_time = get_time();   

  /* this test can also fail if the return value is negative */
  if (resid == NULL){
    printf(">>> FAILED test -- SUNLinSolResid returned NULL N_Vector on Proc %d \n", 
           myid);
    return(1);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolResid\n");
    PRINT_TIME("    SUNLinSolResid Time: %22.15e \n \n", stop_time - start_time);
  }

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolSetATimes Test
 * --------------------------------------------------------------------*/
int Test_SUNLinSolSetATimes(SUNLinearSolver S, void *ATdata,
                            ATimesFn ATimes, int myid)
{
  int     failure;
  double  start_time, stop_time;

  /* try calling SetATimes routine: should pass/fail based on expected input */
  start_time = get_time();
  failure = SUNLinSolSetATimes(S, ATdata, ATimes);
  stop_time = get_time(); 

  if (failure) {
    printf(">>> FAILED test -- SUNLinSolSetATimes returned %d on Proc %d \n", 
           failure, myid);
    return(1);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolSetATimes \n");
    PRINT_TIME("    SUNLinSolSetATimes Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolSetPreconditioner
 * --------------------------------------------------------------------*/
int Test_SUNLinSolSetPreconditioner(SUNLinearSolver S, void *Pdata,
                                    PSetupFn PSetup, PSolveFn PSolve, int myid)
{
  int       failure;
  double    start_time, stop_time;

  /* try calling SetPreconditioner routine: should pass/fail based on expected input */
  start_time = get_time();
  failure = SUNLinSolSetPreconditioner(S, Pdata, PSetup, PSolve);
  stop_time = get_time(); 

  if (failure) {
    printf(">>> FAILED test -- SUNLinSolSetPreconditioner returned %d on Proc %d \n", 
           failure, myid);
    return(1);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolSetPreconditioner \n");
    PRINT_TIME("    SUNLinSolSetPreconditioner Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolSetScalingVectors
 * --------------------------------------------------------------------*/
int Test_SUNLinSolSetScalingVectors(SUNLinearSolver S, N_Vector s1,
                                    N_Vector s2, int myid)
{
  int       failure;
  double    start_time, stop_time;

  /* try calling SetScalingVectors routine: should pass/fail based on expected input */
  start_time = get_time();
  failure = SUNLinSolSetScalingVectors(S, s1, s2);
  stop_time = get_time(); 

  if (failure) {
    printf(">>> FAILED test -- SUNLinSolSetScalingVectors returned %d on Proc %d \n", 
           failure, myid);
    return(1);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolSetScalingVectors \n");
    PRINT_TIME("    SUNLinSolSetScalingVectors Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolInitialize Test
 * --------------------------------------------------------------------*/
int Test_SUNLinSolInitialize(SUNLinearSolver S, int myid)
{
  int       failure;
  double    start_time, stop_time;

  start_time = get_time(); 
  failure = SUNLinSolInitialize(S);
  stop_time = get_time(); 
  
  if (failure) {
    printf(">>> FAILED test -- SUNLinSolInitialize check, Proc %d \n", myid);
    PRINT_TIME("    SUNLinSolInitialize Time: %22.15e \n \n", stop_time - start_time);
    return(1);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolInitialize \n");
    PRINT_TIME("    SUNLinSolInitialize Time: %22.15e \n \n", stop_time - start_time);
  }    
  
  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolSetup Test
 *
 * This test must follow Test_SUNLinSolInitialize
 * --------------------------------------------------------------------*/
int Test_SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A, int myid)
{
  int       failure;
  double    start_time, stop_time;

  start_time = get_time();
  failure = SUNLinSolSetup(S, A);
  stop_time = get_time(); 

  if (failure) {
    printf(">>> FAILED test -- SUNLinSolSetup check, Proc %d \n", myid);
    PRINT_TIME("    SUNLinSolSetup Time: %22.15e \n \n", stop_time - start_time);
    return(1);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolSetup \n");
    PRINT_TIME("    SUNLinSolSetup Time: %22.15e \n \n", stop_time - start_time);
  }    

  return(0);
}


/* ----------------------------------------------------------------------
 * SUNLinSolSolve Test
 *
 * This test must follow Test_SUNLinSolSetup.  Also, x must be the 
 * solution to the linear system A*x = b (for the original A matrix); 
 * while the 'A' that is supplied to this function should have been 
 * 'setup' by the Test_SUNLinSolSetup() function prior to this call.
 * --------------------------------------------------------------------*/
int Test_SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x, 
                        N_Vector b, realtype tol, int myid)
{
  int       failure;
  double    start_time, stop_time;
  N_Vector  y;

  /* clone to create solution vector */
  y = N_VClone(x);
  N_VConst(ZERO, y);

  /* perform solve */
  start_time = get_time();
  failure = SUNLinSolSolve(S, A, y, b, tol);
  stop_time = get_time(); 
  if (failure) {
    printf(">>> FAILED test -- SUNLinSolSolve returned %d on Proc %d \n", 
           failure, myid);
  }
  if (failure < 0) {
    N_VDestroy(y);
    return(1);
  }

  /* Check solution, and copy y into x for return */
  failure = check_vector(x, y, 10.0*tol);
  N_VScale(ONE, y, x);
  if (failure) {
    printf(">>> FAILED test -- SUNLinSolSolve check, Proc %d \n", myid);
    PRINT_TIME("    SUNLinSolSolve Time: %22.15e \n \n", stop_time - start_time);
    N_VDestroy(y);
    return(1);
  }
  else if (myid == 0) {
    printf("    PASSED test -- SUNLinSolSolve \n");
    PRINT_TIME("    SUNLinSolSolve Time: %22.15e \n \n", stop_time - start_time);
  }    

  N_VDestroy(y);
  return(0);
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


