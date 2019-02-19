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
 * This file (companion of fsunlinsol_superlumt.h) contains the
 * implementation needed for the Fortran initialization of superlumt
 * linear solver operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fsunlinsol_superlumt.h"

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

void FSUNSUPERLUMT_INIT(int *code, int *num_threads, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (F2C_CVODE_linsol)  SUNLinSolFree(F2C_CVODE_linsol);
    F2C_CVODE_linsol = NULL;
    F2C_CVODE_linsol = SUNLinSol_SuperLUMT(F2C_CVODE_vec,
                                           F2C_CVODE_matrix,
                                           *num_threads);
    if (F2C_CVODE_linsol == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    if (F2C_IDA_linsol)  SUNLinSolFree(F2C_IDA_linsol);
    F2C_IDA_linsol = NULL;
    F2C_IDA_linsol = SUNLinSol_SuperLUMT(F2C_IDA_vec,
                                         F2C_IDA_matrix,
                                         *num_threads);
    if (F2C_IDA_linsol == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    if (F2C_KINSOL_linsol)  SUNLinSolFree(F2C_KINSOL_linsol);
    F2C_KINSOL_linsol = NULL;
    F2C_KINSOL_linsol = SUNLinSol_SuperLUMT(F2C_KINSOL_vec,
                                            F2C_KINSOL_matrix,
                                            *num_threads);
    if (F2C_KINSOL_linsol == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    if (F2C_ARKODE_linsol)  SUNLinSolFree(F2C_ARKODE_linsol);
    F2C_ARKODE_linsol = NULL;
    F2C_ARKODE_linsol = SUNLinSol_SuperLUMT(F2C_ARKODE_vec,
                                            F2C_ARKODE_matrix,
                                            *num_threads);
    if (F2C_ARKODE_linsol == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FSUNSUPERLUMT_SETORDERING(int *code, int *ordering_choice, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    *ier = SUNLinSol_SuperLUMTSetOrdering(F2C_CVODE_linsol, *ordering_choice);
    break;
  case FCMIX_IDA:
    *ier = SUNLinSol_SuperLUMTSetOrdering(F2C_IDA_linsol, *ordering_choice);
    break;
  case FCMIX_KINSOL:
    *ier = SUNLinSol_SuperLUMTSetOrdering(F2C_KINSOL_linsol, *ordering_choice);
    break;
  case FCMIX_ARKODE:
    *ier = SUNLinSol_SuperLUMTSetOrdering(F2C_ARKODE_linsol, *ordering_choice);
    break;
  default:
    *ier = -1;
  }
}


void FSUNMASSSUPERLUMT_INIT(int *num_threads, int *ier)
{
  *ier = 0;
  if (F2C_ARKODE_mass_sol)  SUNLinSolFree(F2C_ARKODE_mass_sol);
  F2C_ARKODE_mass_sol = NULL;
  F2C_ARKODE_mass_sol = SUNLinSol_SuperLUMT(F2C_ARKODE_vec,
                                            F2C_ARKODE_mass_matrix,
                                            *num_threads);
  if (F2C_ARKODE_mass_sol == NULL) *ier = -1;
}


void FSUNMASSSUPERLUMT_SETORDERING(int *ordering_choice, int *ier)
{
  *ier = 0;
  *ier = SUNLinSol_SuperLUMTSetOrdering(F2C_ARKODE_mass_sol, 
                                        *ordering_choice);
}
