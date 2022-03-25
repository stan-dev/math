/*
 * -----------------------------------------------------------------
 * Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
 *                Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
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

/* Helpers for printing out test status information */
#define TEST_STATUS(fmt, myrank)                     \
  if (print_all_ranks == 0 && myrank == 0) {         \
    printf(fmt);                                     \
  } else if (print_all_ranks) {                      \
    printf("process %06d: ", myrank);                \
    printf(fmt);                                     \
  }
#define TEST_STATUS2(fmt, retval, myrank)            \
  if (print_all_ranks == 0 && myrank == 0) {         \
    printf(fmt, retval);                             \
  } else if (print_all_ranks) {                      \
    printf("process %06d: ", myrank);                \
    printf(fmt, retval);                             \
  }
#define TEST_STATUS3(fmt, val1, val2, myrank)        \
  if (print_all_ranks == 0 && myrank == 0) {         \
    printf(fmt, val1, val2);                         \
  } else if (print_all_ranks) {                      \
    printf("process %06d: ", myrank);                \
    printf(fmt, val1, val2);                         \
  }

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  /* Forward declarations for implementation specific utility functions */
  int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol);
  int check_matrix_entry(SUNMatrix A, realtype val, realtype tol);
  int check_vector(N_Vector expected, N_Vector computed, realtype tol);
  booleantype has_data(SUNMatrix A);
  booleantype is_square(SUNMatrix A);
  void sync_device(SUNMatrix A);

  /* Test function declarations */
  int Test_SUNMatGetID(SUNMatrix A, SUNMatrix_ID sunid, int myid);
  int Test_SUNMatClone(SUNMatrix A, int myid);
  int Test_SUNMatZero(SUNMatrix A, int myid);
  int Test_SUNMatCopy(SUNMatrix A, int myid);
  int Test_SUNMatScaleAdd(SUNMatrix A, SUNMatrix I, int myid);
  int Test_SUNMatScaleAddI(SUNMatrix A, SUNMatrix I, int myid);
  int Test_SUNMatMatvecSetup(SUNMatrix A, int myid);
  int Test_SUNMatMatvec(SUNMatrix A, N_Vector x, N_Vector y, int myid);
  int Test_SUNMatSpace(SUNMatrix A, int myid);

  /* Timing function */
  void SetTiming(int onoff);

  /* Function to set print_all_ranks */
  void SetPrintAllRanks(int onoff);
#ifdef __cplusplus
}
#endif
