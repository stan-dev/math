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
 * This file (companion of fsunlinsol_spfgmr.h) contains the
 * implementation needed for the Fortran initialization of SPFGMR
 * linear solver operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fsunlinsol_spfgmr.h"

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

void FSUNSPFGMR_INIT(int *code, int *pretype, int *maxl, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (F2C_CVODE_linsol)  SUNLinSolFree(F2C_CVODE_linsol);
    F2C_CVODE_linsol = NULL;
    F2C_CVODE_linsol = SUNSPFGMR(F2C_CVODE_vec, *pretype, *maxl);
    if (F2C_CVODE_linsol == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    if (F2C_IDA_linsol)  SUNLinSolFree(F2C_IDA_linsol);
    F2C_IDA_linsol = NULL;
    F2C_IDA_linsol = SUNSPFGMR(F2C_IDA_vec, *pretype, *maxl);
    if (F2C_IDA_linsol == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    if (F2C_KINSOL_linsol)  SUNLinSolFree(F2C_KINSOL_linsol);
    F2C_KINSOL_linsol = NULL;
    F2C_KINSOL_linsol = SUNSPFGMR(F2C_KINSOL_vec, *pretype, *maxl);
    if (F2C_KINSOL_linsol == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    if (F2C_ARKODE_linsol)  SUNLinSolFree(F2C_ARKODE_linsol);
    F2C_ARKODE_linsol = NULL;
    F2C_ARKODE_linsol = SUNSPFGMR(F2C_ARKODE_vec, *pretype, *maxl);
    if (F2C_ARKODE_linsol == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}


void FSUNSPFGMR_SETGSTYPE(int *code, int *gstype, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (!F2C_CVODE_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetGSType(F2C_CVODE_linsol, *gstype);
    break;
  case FCMIX_IDA:
    if (!F2C_IDA_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetGSType(F2C_IDA_linsol, *gstype);
    break;
  case FCMIX_KINSOL:
    if (!F2C_KINSOL_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetGSType(F2C_KINSOL_linsol, *gstype);
    break;
  case FCMIX_ARKODE:
    if (!F2C_ARKODE_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetGSType(F2C_ARKODE_linsol, *gstype);
    break;
  default:
    *ier = -1;
  }
}


void FSUNSPFGMR_SETPRECTYPE(int *code, int *pretype, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (!F2C_CVODE_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetPrecType(F2C_CVODE_linsol, *pretype);
    break;
  case FCMIX_IDA:
    if (!F2C_IDA_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetPrecType(F2C_IDA_linsol, *pretype);
    break;
  case FCMIX_KINSOL:
    if (!F2C_KINSOL_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetPrecType(F2C_KINSOL_linsol, *pretype);
    break;
  case FCMIX_ARKODE:
    if (!F2C_ARKODE_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetPrecType(F2C_ARKODE_linsol, *pretype);
    break;
  default:
    *ier = -1;
  }
}


void FSUNSPFGMR_SETMAXRS(int *code, int *maxrs, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    if (!F2C_CVODE_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetMaxRestarts(F2C_CVODE_linsol, *maxrs);
    break;
  case FCMIX_IDA:
    if (!F2C_IDA_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetMaxRestarts(F2C_IDA_linsol, *maxrs);
    break;
  case FCMIX_KINSOL:
    if (!F2C_KINSOL_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetMaxRestarts(F2C_KINSOL_linsol, *maxrs);
    break;
  case FCMIX_ARKODE:
    if (!F2C_ARKODE_linsol) {
      *ier = -1;
      return;
    }
    *ier = SUNSPFGMRSetMaxRestarts(F2C_ARKODE_linsol, *maxrs);
    break;
  default:
    *ier = -1;
  }
}


void FSUNMASSSPFGMR_INIT(int *pretype, int *maxl, int *ier)
{
  *ier = 0;
  if (F2C_ARKODE_mass_sol)  SUNLinSolFree(F2C_ARKODE_mass_sol);
  F2C_ARKODE_mass_sol = NULL;
  F2C_ARKODE_mass_sol = SUNSPFGMR(F2C_ARKODE_vec, *pretype, *maxl);
  if (F2C_ARKODE_mass_sol == NULL) *ier = -1;
}


void FSUNMASSSPFGMR_SETGSTYPE(int *gstype, int *ier)
{
  *ier = 0;
  if (!F2C_ARKODE_mass_sol) {
      *ier = -1;
      return;
  }
  *ier = SUNSPFGMRSetGSType(F2C_ARKODE_mass_sol, *gstype);
}


void FSUNMASSSPFGMR_SETPRECTYPE(int *pretype, int *ier)
{
  *ier = 0;
  if (!F2C_ARKODE_mass_sol) {
      *ier = -1;
      return;
  }
  *ier = SUNSPFGMRSetPrecType(F2C_ARKODE_mass_sol, *pretype);
}


void FSUNMASSSPFGMR_SETMAXRS(int *maxrs, int *ier)
{
  *ier = 0;
  if (!F2C_ARKODE_mass_sol) {
      *ier = -1;
      return;
  }
  *ier = SUNSPFGMRSetMaxRestarts(F2C_ARKODE_mass_sol, *maxrs);
}
