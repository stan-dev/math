/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * This is the implementation file for the MPIPlusX NVECTOR.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <nvector/nvector_mpimanyvector.h>
#include <nvector/nvector_mpiplusx.h>
#include <sundials/sundials_math.h>

#define MPIPLUSX_LOCAL_VECTOR(v) ( N_VGetSubvector_MPIManyVector(v, 0) )

N_Vector N_VMake_MPIPlusX(MPI_Comm comm, N_Vector X)
{
  N_Vector v;

  if (X == NULL) return NULL;

  v = NULL;
  v = N_VMake_MPIManyVector(comm, 1, &X);
  if (v == NULL) return NULL;

  /* override certain ops */
  v->ops->nvgetvectorid = N_VGetVectorID_MPIPlusX;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_MPIPlusX;
  v->ops->nvsetarraypointer = N_VSetArrayPointer_MPIPlusX;

  return v;
}

N_Vector_ID N_VGetVectorID_MPIPlusX(N_Vector v)
{
  return SUNDIALS_NVEC_MPIPLUSX;
}

realtype* N_VGetArrayPointer_MPIPlusX(N_Vector v)
{
  return N_VGetSubvectorArrayPointer_MPIManyVector(v, 0);
}

void N_VSetArrayPointer_MPIPlusX(realtype *vdata, N_Vector v)
{
  N_VSetSubvectorArrayPointer_MPIManyVector(vdata, v, 0);
}

N_Vector N_VGetLocalVector_MPIPlusX(N_Vector v)
{
  return MPIPLUSX_LOCAL_VECTOR(v);
}

sunindextype N_VGetLocalLength_MPIPlusX(N_Vector v)
{
  return N_VGetLength(MPIPLUSX_LOCAL_VECTOR(v));
}
