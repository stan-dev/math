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
 * This file (companion of fsunmatrix_dense.c) contains the
 * definitions needed for the initialization of dense
 * matrix operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNMATRIX_DENSE_H
#define _FSUNMATRIX_DENSE_H

#include <sunmatrix/sunmatrix_dense.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNDENSEMAT_INIT     SUNDIALS_F77_FUNC(fsundensematinit, FSUNDENSEMATINIT)
#define FSUNDENSEMASSMAT_INIT SUNDIALS_F77_FUNC(fsundensemassmatinit, FSUNDENSEMASSMATINIT)
#else
#define FSUNDENSEMAT_INIT     fsundensematinit_
#define FSUNDENSEMASSMAT_INIT fsundensemassmatinit_
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
 * FSUNDENSEMAT_INIT - initializes dense matrix operations for main problem
 * FSUNDENSEMASSMAT_INIT - initializes dense matrix operations for mass matrix solver
 */

void FSUNDENSEMAT_INIT(int *code, long int *M, long int *N, int *ier);
void FSUNDENSEMASSMAT_INIT(long int *M, long int *N, int *ier);

#ifdef __cplusplus
}
#endif

#endif
