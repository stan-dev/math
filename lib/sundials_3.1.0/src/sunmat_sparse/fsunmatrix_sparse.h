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
 * This file (companion of fsunmatrix_sparse.c) contains the
 * definitions needed for the initialization of sparse
 * matrix operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNMATRIX_SPARSE_H
#define _FSUNMATRIX_SPARSE_H

#include <sunmatrix/sunmatrix_sparse.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNSPARSEMAT_INIT     SUNDIALS_F77_FUNC(fsunsparsematinit, FSUNSPARSEMATINIT)
#define FSUNSPARSEMASSMAT_INIT SUNDIALS_F77_FUNC(fsunsparsemassmatinit, FSUNSPARSEMASSMATINIT)
#else
#define FSUNSPARSEMAT_INIT     fsunsparsematinit_
#define FSUNSPARSEMASSMAT_INIT fsunsparsemassmatinit_
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
 * FSUNSPARSEMAT_INIT - initializes sparse matrix operations for main problem
 * FSUNSPARSEMASSMAT_INIT - initializes sparse matrix operations for mass matrix solve
 */

void FSUNSPARSEMAT_INIT(int *code, long int *M, long int *N,
                        long int *NNZ, int *sparsetype, int *ier);

void FSUNSPARSEMASSMAT_INIT(long int *M, long int *N,
                            long int *NNZ, int *sparsetype, int *ier);

#ifdef __cplusplus
}
#endif

#endif
