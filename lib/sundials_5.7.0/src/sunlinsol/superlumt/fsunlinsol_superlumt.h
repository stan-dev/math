/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
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
 * This file (companion of fsunlinsol_superlumt.c) contains the
 * definitions needed for the initialization of superlumt
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_SUPERLUMT_H
#define _FSUNLINSOL_SUPERLUMT_H

#include <sunlinsol/sunlinsol_superlumt.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNSUPERLUMT_INIT            SUNDIALS_F77_FUNC(fsunsuperlumtinit,        FSUNSUPERLUMTINIT)
#define FSUNSUPERLUMT_SETORDERING     SUNDIALS_F77_FUNC(fsunsuperlumtsetordering, FSUNSUPERLUMTSETORDERING)
#define FSUNMASSSUPERLUMT_INIT        SUNDIALS_F77_FUNC(fsunmasssuperlumtinit,        FSUNMASSSUPERLUMTINIT)
#define FSUNMASSSUPERLUMT_SETORDERING SUNDIALS_F77_FUNC(fsunmasssuperlumtsetordering, FSUNMASSSUPERLUMTSETORDERING)
#else
#define FSUNSUPERLUMT_INIT            fsunsuperlumtinit_
#define FSUNSUPERLUMT_SETORDERING     fsunsuperlumtsetordering_
#define FSUNMASSSUPERLUMT_INIT        fsunmasssuperlumtinit_
#define FSUNMASSSUPERLUMT_SETORDERING fsunmasssuperlumtsetordering_
#endif


/* Declarations of global variables */

extern SUNLinearSolver F2C_CVODE_linsol;
extern SUNLinearSolver F2C_IDA_linsol;
extern SUNLinearSolver F2C_KINSOL_linsol;
extern SUNLinearSolver F2C_ARKODE_linsol;
extern SUNLinearSolver F2C_ARKODE_mass_sol;

/* 
 * Prototypes of exported functions 
 *
 * FSUNSUPERLUMT_INIT - initializes superlumt linear solver for main problem
 * FSUNSUPERLUMT_SETORDERING - sets the ordering choice used by SUPERLUMT for main problem
 * FSUNMASSSUPERLUMT_INIT - initializes superlumt linear solver for mass matrix
 * FSUNMASSSUPERLUMT_SETORDERING - sets the ordering choice used by SUPERLUMT for mass matrix
 */

void FSUNSUPERLUMT_INIT(int *code, int *num_threads, int *ier);
void FSUNSUPERLUMT_SETORDERING(int *code, int *ordering, int *ier);
void FSUNMASSSUPERLUMT_INIT(int *num_threads, int *ier);
void FSUNMASSSUPERLUMT_SETORDERING(int *ordering, int *ier);

#ifdef __cplusplus
}
#endif

#endif
