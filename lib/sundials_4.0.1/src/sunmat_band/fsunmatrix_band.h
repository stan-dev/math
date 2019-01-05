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
 * This file (companion of fsunmatrix_band.c) contains the
 * definitions needed for the initialization of band
 * matrix operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNMATRIX_BAND_H
#define _FSUNMATRIX_BAND_H

#include <sunmatrix/sunmatrix_band.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNBANDMAT_INIT     SUNDIALS_F77_FUNC(fsunbandmatinit, FSUNBANDMATINIT)
#define FSUNBANDMASSMAT_INIT SUNDIALS_F77_FUNC(fsunbandmassmatinit, FSUNBANDMASSMATINIT)
#else
#define FSUNBANDMAT_INIT     fsunbandmatinit_
#define FSUNBANDMASSMAT_INIT fsunbandmassmatinit_
#endif


/* Declarations of global variables */

extern SUNMatrix F2C_CVODE_matrix;
extern SUNMatrix F2C_IDA_matrix;
extern SUNMatrix F2C_KINSOL_matrix;
extern SUNMatrix F2C_ARKODE_matrix;
extern SUNMatrix F2C_ARKODE_mass_matrix;

/* 
 * Prototypes of exported functions 
 *
 * FSUNBANDMAT_INIT - initializes band matrix operations for main problem
 * FSUNBANDMASSMAT_INIT - initializes band matrix operations for mass matrix solve
 */

void FSUNBANDMAT_INIT(int *code, long int *N, long int *mu, long int *ml, int *ier);
void FSUNBANDMASSMAT_INIT(long int *N, long int *mu, long int *ml, int *ier);

#ifdef __cplusplus
}
#endif

#endif
