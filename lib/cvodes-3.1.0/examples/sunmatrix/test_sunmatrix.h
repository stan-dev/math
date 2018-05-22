/*
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 *                Daniel R. Reynolds @ SMU
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
 * This is the header file contains the prototypes for functions to 
 * test SUNMatrix module implementation. 
 * -----------------------------------------------------------------
 */

#include <math.h>

/* define constatnts */
#define NEG_TWO  RCONST(-2.0)
#define NEG_ONE  RCONST(-1.0)
#define NEG_HALF RCONST(-0.5)
#define ZERO     RCONST(0.0)
#define HALF     RCONST(0.5)
#define ONE      RCONST(1.0)
#define TWO      RCONST(2.0)
#define THREE    RCONST(3.0)

/* NAN and floating point "equality" check, failure update macro */
#if __STDC_VERSION__ >= 199901L
#define FNEQ(a,b,tol) (isnan(a) ? 1 : ( SUNRabs((a)-(b))/SUNRabs(b) > tol ))
#else
#define FNEQ(a,b,tol) (( SUNRabs((a)-(b))/SUNRabs(b) > tol ))
#endif


#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  /* Forward declarations for implementation specific utility functions */
  int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol);
  int check_matrix_entry(SUNMatrix A, realtype val, realtype tol);
  int check_vector(N_Vector X, N_Vector Y, realtype tol);
  booleantype has_data(SUNMatrix A);
  booleantype is_square(SUNMatrix A);
   
  /* Test function declarations */
  int Test_SUNMatGetID(SUNMatrix A, SUNMatrix_ID sunid, int myid);
  int Test_SUNMatClone(SUNMatrix A, int myid);
  int Test_SUNMatZero(SUNMatrix A, int myid);
  int Test_SUNMatCopy(SUNMatrix A, int myid);
  int Test_SUNMatScaleAdd(SUNMatrix A, SUNMatrix I, int myid);
  int Test_SUNMatScaleAddI(SUNMatrix A, SUNMatrix I, int myid);
  int Test_SUNMatMatvec(SUNMatrix A, N_Vector x, N_Vector y, int myid);
  int Test_SUNMatSpace(SUNMatrix A, int myid);

  /* Timing function */
  void SetTiming(int onoff);
  
#ifdef __cplusplus
}
#endif
