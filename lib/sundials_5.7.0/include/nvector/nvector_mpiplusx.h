/* -----------------------------------------------------------------
 * Programmer(s): Cody Balos @ LLNL
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
 * This is the header file for the MPI+X implementation of the
 * NVECTOR module. The MPIPlusX NVECTOR is really just an extension
 * of the ManyVector.
 * -----------------------------------------------------------------*/

#ifndef _NVECTOR_MPIPLUSX_H
#define _NVECTOR_MPIPLUSX_H

#include <mpi.h>
#include <nvector/nvector_mpimanyvector.h>
#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef N_VectorContent_MPIManyVector N_VectorContent_MPIPlusX;


SUNDIALS_EXPORT N_Vector N_VMake_MPIPlusX(MPI_Comm comm, N_Vector X);

SUNDIALS_EXPORT N_Vector_ID N_VGetVectorID_MPIPlusX(N_Vector v);

SUNDIALS_EXPORT realtype* N_VGetArrayPointer_MPIPlusX(N_Vector v);

SUNDIALS_EXPORT void N_VSetArrayPointer_MPIPlusX(realtype *vdata, N_Vector v);

SUNDIALS_EXPORT N_Vector N_VGetLocalVector_MPIPlusX(N_Vector v);

SUNDIALS_EXPORT sunindextype N_VGetLocalLength_MPIPlusX(N_Vector v);

SUNDIALS_STATIC_INLINE
int N_VEnableFusedOps_MPIPlusX(N_Vector v, booleantype tf)
{ return N_VEnableFusedOps_MPIManyVector(v, tf); }

#ifdef __cplusplus
}
#endif

#endif
