/*
 * -----------------------------------------------------------------
 * Programmer(s): Steven Smith @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file (companion of nvector_openmp.h) contains the
 * implementation needed for the Fortran initialization of openmp
 * vector operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fnvector_openmp.h"

/* Define global vector variables */

N_Vector F2C_CVODE_vec;
N_Vector F2C_CVODE_vecQ;
N_Vector *F2C_CVODE_vecS;
N_Vector F2C_CVODE_vecB;
N_Vector F2C_CVODE_vecQB;

N_Vector F2C_IDA_vec;
N_Vector F2C_IDA_vecQ;
N_Vector *F2C_IDA_vecS;
N_Vector F2C_IDA_vecB;
N_Vector F2C_IDA_vecQB;

N_Vector F2C_KINSOL_vec;

N_Vector F2C_ARKODE_vec;

/* Fortran callable interfaces */

void FNV_INITOMP(int *code, long int *N, int *num_threads, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vec = NULL;
    F2C_CVODE_vec = N_VNewEmpty_OpenMP((sunindextype)(*N), *num_threads);
    if (F2C_CVODE_vec == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vec = NULL;
    F2C_IDA_vec = N_VNewEmpty_OpenMP((sunindextype)(*N), *num_threads);
    if (F2C_IDA_vec == NULL) *ier = -1;
    break;
  case FCMIX_KINSOL:
    F2C_KINSOL_vec = NULL;
    F2C_KINSOL_vec = N_VNewEmpty_OpenMP((sunindextype)(*N), *num_threads);
    if (F2C_KINSOL_vec == NULL) *ier = -1;
    break;
  case FCMIX_ARKODE:
    F2C_ARKODE_vec = NULL;
    F2C_ARKODE_vec = N_VNewEmpty_OpenMP((sunindextype)(*N), *num_threads);
    if (F2C_ARKODE_vec == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITOMP_Q(int *code, long int *Nq, int *num_threads, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQ = NULL;
    F2C_CVODE_vecQ = N_VNewEmpty_OpenMP((sunindextype)(*Nq), *num_threads);
    if (F2C_CVODE_vecQ == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQ = NULL;
    F2C_IDA_vecQ = N_VNewEmpty_OpenMP((sunindextype)(*Nq), *num_threads);
    if (F2C_IDA_vecQ == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITOMP_B(int *code, long int *NB, int *num_threads, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecB = NULL;
    F2C_CVODE_vecB = N_VNewEmpty_OpenMP((sunindextype)(*NB), *num_threads);
    if (F2C_CVODE_vecB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecB = NULL;
    F2C_IDA_vecB = N_VNewEmpty_OpenMP((sunindextype)(*NB), *num_threads);
    if (F2C_IDA_vecB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITOMP_QB(int *code, long int *NqB, int *num_threads, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecQB = NULL;
    F2C_CVODE_vecQB = N_VNewEmpty_OpenMP((sunindextype)(*NqB), *num_threads);
    if (F2C_CVODE_vecQB == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecQB = NULL;
    F2C_IDA_vecQB = N_VNewEmpty_OpenMP((sunindextype)(*NqB), *num_threads);
    if (F2C_IDA_vecQB == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}

void FNV_INITOMP_S(int *code, int *Ns, int *ier)
{
  *ier = 0;

  switch(*code) {
  case FCMIX_CVODE:
    F2C_CVODE_vecS = NULL;
    F2C_CVODE_vecS = (N_Vector *) N_VCloneVectorArrayEmpty_OpenMP(*Ns, F2C_CVODE_vec);
    if (F2C_CVODE_vecS == NULL) *ier = -1;
    break;
  case FCMIX_IDA:
    F2C_IDA_vecS = NULL;
    F2C_IDA_vecS = (N_Vector *) N_VCloneVectorArrayEmpty_OpenMP(*Ns, F2C_IDA_vec);
    if (F2C_IDA_vecS == NULL) *ier = -1;
    break;
  default:
    *ier = -1;
  }
}
