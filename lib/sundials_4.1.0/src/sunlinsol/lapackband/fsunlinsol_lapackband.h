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
 * This file (companion of fsunlinsol_lapackband.c) contains the
 * definitions needed for the initialization of LAPACK band
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_LAPBAND_H
#define _FSUNLINSOL_LAPBAND_H

#include <sunlinsol/sunlinsol_lapackband.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNLAPACKBAND_INIT     SUNDIALS_F77_FUNC(fsunlapackbandinit, FSUNLAPACKBANDINIT)
#define FSUNMASSLAPACKBAND_INIT SUNDIALS_F77_FUNC(fsunmasslapackbandinit, FSUNMASSLAPACKBANDINIT)
#else
#define FSUNLAPACKBAND_INIT     fsunlapackbandinit_
#define FSUNMASSLAPACKBAND_INIT fsunmasslapackbandinit_
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
 * FSUNLAPACKBAND_INIT - initializes LAPACK band linear solver for main problem
 * FSUNMASSLAPACKBAND_INIT - initializes LAPACK band linear solver for mass matrix solve
 */

void FSUNLAPACKBAND_INIT(int *code, int *ier);
void FSUNMASSLAPACKBAND_INIT(int *ier);

#ifdef __cplusplus
}
#endif

#endif
