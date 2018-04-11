/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
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
