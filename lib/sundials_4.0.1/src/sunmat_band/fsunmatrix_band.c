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
 * This file (companion of fsunmatrix_band.h) contains the
 * implementation needed for the Fortran initialization of band
 * vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fsunmatrix_band.h"

/* Define global matrix variables */

SUNMatrix F2C_CVODE_matrix;
SUNMatrix F2C_IDA_matrix;
SUNMatrix F2C_KINSOL_matrix;
SUNMatrix F2C_ARKODE_matrix;
SUNMatrix F2C_ARKODE_mass_matrix;

/* Fortran callable interfaces */

void FSUNBANDMAT_INIT(int *code, long int *N, long int *mu, 
                      long int *ml, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (F2C_CVODE_matrix)  SUNMatDestroy(F2C_CVODE_matrix);
    F2C_CVODE_matrix = NULL;
    F2C_CVODE_matrix = SUNBandMatrix(*N, *mu, *ml);
    if (F2C_CVODE_matrix == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    if (F2C_IDA_matrix)  SUNMatDestroy(F2C_IDA_matrix);
    F2C_IDA_matrix = NULL;
    F2C_IDA_matrix = SUNBandMatrix(*N, *mu, *ml);
    if (F2C_IDA_matrix == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    if (F2C_KINSOL_matrix)  SUNMatDestroy(F2C_KINSOL_matrix);
    F2C_KINSOL_matrix = NULL;
    F2C_KINSOL_matrix = SUNBandMatrix(*N, *mu, *ml);
    if (F2C_KINSOL_matrix == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    if (F2C_ARKODE_matrix)  SUNMatDestroy(F2C_ARKODE_matrix);
    F2C_ARKODE_matrix = NULL;
    F2C_ARKODE_matrix = SUNBandMatrix(*N, *mu, *ml);
    if (F2C_ARKODE_matrix == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FSUNBANDMASSMAT_INIT(long int *N, long int *mu, 
                          long int *ml, int *ier)
{
  *ier = 0;
  if (F2C_ARKODE_mass_matrix)  SUNMatDestroy(F2C_ARKODE_mass_matrix);
  F2C_ARKODE_mass_matrix = NULL;
  F2C_ARKODE_mass_matrix = SUNBandMatrix(*N, *mu, *ml);
  if (F2C_ARKODE_mass_matrix == NULL) *ier = -1;
}
