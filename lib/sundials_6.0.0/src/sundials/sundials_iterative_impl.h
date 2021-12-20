/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner and Shelby Lockhart @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the implementation header file for SUNDIALS functions used by
 * different iterative solvers.
 * ---------------------------------------------------------------------------*/

#include <sundials/sundials_iterative.h>

/* -----------------------------------------------------------------------------
 * Type: SUNQRData
 * -----------------------------------------------------------------------------
 * A SUNQRData struct holds temporary workspace vectors and realtype arrays for
 * a SUNQRAddFn. The N_Vectors and realtype arrays it contains are created by
 * the routine calling a SUNQRAdd function.
 * ---------------------------------------------------------------------------*/

typedef struct _SUNQRData *SUNQRData;

struct _SUNQRData
{
  N_Vector vtemp;
  N_Vector vtemp2;
  realtype *temp_array;
};
