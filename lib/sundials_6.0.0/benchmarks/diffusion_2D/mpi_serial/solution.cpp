/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Serial solution and derivative functions
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

// Compute the exact solution
int Solution(realtype t, N_Vector u, UserData *udata)
{
  // Initialize u to zero (handles boundary conditions)
  N_VConst(ZERO, u);

  realtype x, y;
  realtype sin_sqr_x, sin_sqr_y;

  // Constants for computing solution
  realtype cos_sqr_t = cos(PI * t) * cos(PI * t);

  // Iterative over domain interior
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? udata->nx_loc : udata->nx_loc - 1;

  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? udata->ny_loc : udata->ny_loc - 1;

  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  for (sunindextype j = jstart; j < jend; j++)
  {
    for (sunindextype i = istart; i < iend; i++)
    {
      x  = (udata->is + i) * udata->dx;
      y  = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uarray[IDX(i,j,udata->nx_loc)] = sin_sqr_x * sin_sqr_y * cos_sqr_t;
    }
  }

  return 0;
}


// Compute the exact solution derivative
int SolutionDerivative(realtype t, N_Vector up, UserData *udata)
{
  // Initialize u to zero (handles boundary conditions)
  N_VConst(ZERO, up);

  realtype x, y;
  realtype sin_sqr_x, sin_sqr_y;

  // Constants for computing solution
  realtype cos_sin_t = -TWO * PI * cos(PI * t) * sin(PI * t);

  // Iterative over domain interior
  sunindextype istart = (udata->HaveNbrW) ? 0 : 1;
  sunindextype iend   = (udata->HaveNbrE) ? udata->nx_loc : udata->nx_loc - 1;

  sunindextype jstart = (udata->HaveNbrS) ? 0 : 1;
  sunindextype jend   = (udata->HaveNbrN) ? udata->ny_loc : udata->ny_loc - 1;

  realtype *uparray = N_VGetArrayPointer(up);
  if (check_flag((void *) uparray, "N_VGetArrayPointer", 0)) return -1;

  for (sunindextype j = jstart; j < jend; j++)
  {
    for (sunindextype i = istart; i < iend; i++)
    {
      x  = (udata->is + i) * udata->dx;
      y  = (udata->js + j) * udata->dy;

      sin_sqr_x = sin(PI * x) * sin(PI * x);
      sin_sqr_y = sin(PI * y) * sin(PI * y);

      uparray[IDX(i,j,udata->nx_loc)] = sin_sqr_x * sin_sqr_y * cos_sin_t;
    }
  }

  return 0;
}
