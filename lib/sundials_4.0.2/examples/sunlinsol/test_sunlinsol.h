/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * test SUNLinearSolver module implementations. 
 * -----------------------------------------------------------------
 */

#include <math.h>

/* define constatnts */
#define ZERO     RCONST(0.0)
#define ONE      RCONST(1.0)

/* NAN and floating point "equality" check, failure update macro */
#if __STDC_VERSION__ >= 199901L
#define FNEQ(a,b,tol) (isnan(a) ? 1 : ( SUNRabs((a)-(b)) > tol ))
#else
#define FNEQ(a,b,tol) (( SUNRabs((a)-(b)) > tol ))
#endif


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  /* Forward declarations for implementation specific utility functions */
  int check_vector(N_Vector X, N_Vector Y, realtype tol);
   
  /* Test function declarations */
  int Test_SUNLinSolGetType(SUNLinearSolver S, SUNLinearSolver_Type sunid, int myid);
  int Test_SUNLinSolLastFlag(SUNLinearSolver S, int myid);
  int Test_SUNLinSolSpace(SUNLinearSolver S, int myid);
  int Test_SUNLinSolNumIters(SUNLinearSolver S, int myid);
  int Test_SUNLinSolResNorm(SUNLinearSolver S, int myid);
  int Test_SUNLinSolResid(SUNLinearSolver S, int myid);
  int Test_SUNLinSolSetATimes(SUNLinearSolver S, void *ATdata,
                              ATimesFn ATimes, int myid);
  int Test_SUNLinSolSetPreconditioner(SUNLinearSolver S, void *Pdata,
                                      PSetupFn PSetup, PSolveFn PSolve, int myid);
  int Test_SUNLinSolSetScalingVectors(SUNLinearSolver S, N_Vector s1, 
                                      N_Vector s2, int myid);
  int Test_SUNLinSolInitialize(SUNLinearSolver S, int myid);
  int Test_SUNLinSolSetup(SUNLinearSolver S, SUNMatrix A, int myid);
  int Test_SUNLinSolSolve(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                          N_Vector b, realtype tol, int myid);

  /* Timing function */
  void SetTiming(int onoff);
  
#ifdef __cplusplus
}
#endif
