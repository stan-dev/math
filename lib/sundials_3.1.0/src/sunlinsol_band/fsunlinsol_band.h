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
 * This file (companion of fsunlinsol_band.c) contains the
 * definitions needed for the initialization of band
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_BAND_H
#define _FSUNLINSOL_BAND_H

#include <sunlinsol/sunlinsol_band.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNBANDLINSOL_INIT     SUNDIALS_F77_FUNC(fsunbandlinsolinit, FSUNBANDLINSOLINIT)
#define FSUNMASSBANDLINSOL_INIT SUNDIALS_F77_FUNC(fsunmassbandlinsolinit, FSUNMASSBANDLINSOLINIT)
#else
#define FSUNBANDLINSOL_INIT     fsunbandlinsolinit_
#define FSUNMASSBANDLINSOL_INIT fsunmassbandlinsolinit_
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
 * FSUNBANDLINSOL_INIT - initializes band linear solver for main problem
 * FSUNMASSBANDLINSOL_INIT - initializes band linear solver for mass matrix solve
 */

void FSUNBANDLINSOL_INIT(int *code, int *ier);
void FSUNMASSBANDLINSOL_INIT(int *ier);

#ifdef __cplusplus
}
#endif

#endif
