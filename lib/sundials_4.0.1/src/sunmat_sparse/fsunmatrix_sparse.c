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
 * This file (companion of fsunmatrix_sparse.h) contains the
 * implementation needed for the Fortran initialization of sparse
 * vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fsunmatrix_sparse.h"

/* Define global matrix variables */

SUNMatrix F2C_CVODE_matrix;
SUNMatrix F2C_IDA_matrix;
SUNMatrix F2C_KINSOL_matrix;
SUNMatrix F2C_ARKODE_matrix;
SUNMatrix F2C_ARKODE_mass_matrix;

/* Fortran callable interfaces */

void FSUNSPARSEMAT_INIT(int *code, long int *M, long int *N,
                        long int *NNZ, int *sparsetype, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (F2C_CVODE_matrix)  SUNMatDestroy(F2C_CVODE_matrix);
    F2C_CVODE_matrix = NULL;
    F2C_CVODE_matrix = SUNSparseMatrix(*M, *N, *NNZ, *sparsetype);
    if (F2C_CVODE_matrix == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    if (F2C_IDA_matrix)  SUNMatDestroy(F2C_IDA_matrix);
    F2C_IDA_matrix = NULL;
    F2C_IDA_matrix = SUNSparseMatrix(*M, *N, *NNZ, *sparsetype);
    if (F2C_IDA_matrix == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    if (F2C_KINSOL_matrix)  SUNMatDestroy(F2C_KINSOL_matrix);
    F2C_KINSOL_matrix = NULL;
    F2C_KINSOL_matrix = SUNSparseMatrix(*M, *N, *NNZ, *sparsetype);
    if (F2C_KINSOL_matrix == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    if (F2C_ARKODE_matrix)  SUNMatDestroy(F2C_ARKODE_matrix);
    F2C_ARKODE_matrix = NULL;
    F2C_ARKODE_matrix = SUNSparseMatrix(*M, *N, *NNZ, *sparsetype);
    if (F2C_ARKODE_matrix == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FSUNSPARSEMASSMAT_INIT(long int *M, long int *N, long int *NNZ, 
                            int *sparsetype, int *ier)
{
  *ier = 0;
  if (F2C_ARKODE_mass_matrix)  SUNMatDestroy(F2C_ARKODE_mass_matrix);
  F2C_ARKODE_mass_matrix = NULL;
  F2C_ARKODE_mass_matrix = SUNSparseMatrix(*M, *N, *NNZ, *sparsetype);
  if (F2C_ARKODE_mass_matrix == NULL) *ier = -1;
}
