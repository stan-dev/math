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
 * This file (companion of fsunlinsol_klu.h) contains the
 * implementation needed for the Fortran initialization of klu
 * linear solver operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fsunlinsol_klu.h"

/* Define global linsol variables */

SUNLinearSolver F2C_CVODE_linsol;
SUNLinearSolver F2C_IDA_linsol;
SUNLinearSolver F2C_KINSOL_linsol;
SUNLinearSolver F2C_ARKODE_linsol;
SUNLinearSolver F2C_ARKODE_mass_sol;

/* Declarations of external global variables */

extern SUNMatrix F2C_CVODE_matrix;
extern SUNMatrix F2C_IDA_matrix;
extern SUNMatrix F2C_KINSOL_matrix;
extern SUNMatrix F2C_ARKODE_matrix;
extern SUNMatrix F2C_ARKODE_mass_matrix;

extern N_Vector F2C_CVODE_vec;
extern N_Vector F2C_IDA_vec;
extern N_Vector F2C_KINSOL_vec;
extern N_Vector F2C_ARKODE_vec;

/* Fortran callable interfaces */

void FSUNKLU_INIT(int *code, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (F2C_CVODE_linsol)  SUNLinSolFree(F2C_CVODE_linsol);
    F2C_CVODE_linsol = NULL;
    F2C_CVODE_linsol = SUNLinSol_KLU(F2C_CVODE_vec, F2C_CVODE_matrix);
    if (F2C_CVODE_linsol == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    if (F2C_IDA_linsol)  SUNLinSolFree(F2C_IDA_linsol);
    F2C_IDA_linsol = NULL;
    F2C_IDA_linsol = SUNLinSol_KLU(F2C_IDA_vec, F2C_IDA_matrix);
    if (F2C_IDA_linsol == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    if (F2C_KINSOL_linsol)  SUNLinSolFree(F2C_KINSOL_linsol);
    F2C_KINSOL_linsol = NULL;
    F2C_KINSOL_linsol = SUNLinSol_KLU(F2C_KINSOL_vec, F2C_KINSOL_matrix);
    if (F2C_KINSOL_linsol == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    if (F2C_ARKODE_linsol)  SUNLinSolFree(F2C_ARKODE_linsol);
    F2C_ARKODE_linsol = NULL;
    F2C_ARKODE_linsol = SUNLinSol_KLU(F2C_ARKODE_vec, F2C_ARKODE_matrix);
    if (F2C_ARKODE_linsol == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FSUNKLU_REINIT(int *code, long int *NNZ, int *reinit_type, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    *ier = SUNLinSol_KLUReInit(F2C_CVODE_linsol, F2C_CVODE_matrix,
                               *NNZ, *reinit_type);
    break;
  case FCMIX_IDA:
    *ier = SUNLinSol_KLUReInit(F2C_IDA_linsol, F2C_IDA_matrix,
                               *NNZ, *reinit_type);
    break;
  case FCMIX_KINSOL:
    *ier = SUNLinSol_KLUReInit(F2C_KINSOL_linsol, F2C_KINSOL_matrix,
                               *NNZ, *reinit_type);
    break;
  case FCMIX_ARKODE:
    *ier = SUNLinSol_KLUReInit(F2C_ARKODE_linsol, F2C_ARKODE_matrix,
                               *NNZ, *reinit_type);
    break;
  default:
    *ier = -1;
  }
}


void FSUNKLU_SETORDERING(int *code, int *ordering_choice, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    *ier = SUNLinSol_KLUSetOrdering(F2C_CVODE_linsol, *ordering_choice);
    break;
  case FCMIX_IDA:
    *ier = SUNLinSol_KLUSetOrdering(F2C_IDA_linsol, *ordering_choice);
    break;
  case FCMIX_KINSOL:
    *ier = SUNLinSol_KLUSetOrdering(F2C_KINSOL_linsol, *ordering_choice);
    break;
  case FCMIX_ARKODE:
    *ier = SUNLinSol_KLUSetOrdering(F2C_ARKODE_linsol, *ordering_choice);
    break;
  default:
    *ier = -1;
  }
}


void FSUNMASSKLU_INIT(int *ier)
{
  *ier = 0;
  if (F2C_ARKODE_mass_sol)  SUNLinSolFree(F2C_ARKODE_mass_sol);
  F2C_ARKODE_mass_sol = NULL;
  F2C_ARKODE_mass_sol = SUNLinSol_KLU(F2C_ARKODE_vec, 
                                      F2C_ARKODE_mass_matrix);
  if (F2C_ARKODE_mass_sol == NULL) *ier = -1;
}


void FSUNMASSKLU_REINIT(long int *NNZ, int *reinit_type, int *ier)
{
  *ier = 0;
  *ier = SUNLinSol_KLUReInit(F2C_ARKODE_mass_sol, F2C_ARKODE_mass_matrix,
                             *NNZ, *reinit_type);
}


void FSUNMASSKLU_SETORDERING(int *ordering_choice, int *ier)
{
  *ier = 0;
  *ier = SUNLinSol_KLUSetOrdering(F2C_ARKODE_mass_sol, *ordering_choice);
}
