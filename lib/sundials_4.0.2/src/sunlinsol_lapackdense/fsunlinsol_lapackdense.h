/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
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
 * This file (companion of fsunlinsol_lapackdense.c) contains the
 * definitions needed for the initialization of LAPACK dense
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_LAPDENSE_H
#define _FSUNLINSOL_LAPDENSE_H

#include <sunlinsol/sunlinsol_lapackdense.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNLAPACKDENSE_INIT     SUNDIALS_F77_FUNC(fsunlapackdenseinit, FSUNLAPACKDENSEINIT)
#define FSUNMASSLAPACKDENSE_INIT SUNDIALS_F77_FUNC(fsunmasslapackdenseinit, FSUNMASSLAPACKDENSEINIT)
#else
#define FSUNLAPACKDENSE_INIT     fsunlapackdenseinit_
#define FSUNMASSLAPACKDENSE_INIT fsunmasslapackdenseinit_
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
 * FSUNLAPACKDENSE_INIT - initializes LAPACK dense linear solver for main problem
 * FSUNMASSLAPACKDENSE_INIT - initializes LAPACK dense linear solver for mass matrix solve
 */

void FSUNLAPACKDENSE_INIT(int *code, int *ier);
void FSUNMASSLAPACKDENSE_INIT(int *ier);

#ifdef __cplusplus
}
#endif

#endif
