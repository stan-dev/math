/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * Implementation header file for the main ARKODE integrator.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_MRI_TABLES_IMPL_H
#define _ARKODE_MRI_TABLES_IMPL_H

#include "arkode_mristep_impl.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* Returns the stage type (implicit/explicit + fast/nofast) */
int mriStepCoupling_GetStageType(MRIStepCoupling MRIC, int is);

/* Returns index maps for where to store stage RHS evaluations */
int mriStepCoupling_GetStageMap(MRIStepCoupling MRIC, int* stage_map,
                                int* nstored_stages);


#ifdef __cplusplus
}
#endif

#endif
