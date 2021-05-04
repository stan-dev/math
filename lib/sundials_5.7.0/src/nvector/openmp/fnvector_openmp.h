/*
 * -----------------------------------------------------------------
 * Programmer(s): Steven Smith @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file (companion of nvector_openmp.h) contains the
 * definitions needed for the initialization of openmp
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_OPENMP_H
#define _FNVECTOR_OPENMP_H

#include <nvector/nvector_openmp.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FNV_INITOMP    SUNDIALS_F77_FUNC(fnvinitomp, FNVINITOMP)
#else
#define FNV_INITOMP    fnvinitomp_
#endif

#if defined(SUNDIALS_F77_FUNC_)

#define FNV_INITOMP_Q  SUNDIALS_F77_FUNC_(fnvinitomp_q, FNVINITOMP_Q)
#define FNV_INITOMP_S  SUNDIALS_F77_FUNC_(fnvinitomp_s, FNVINITOMP_S)
#define FNV_INITOMP_B  SUNDIALS_F77_FUNC_(fnvinitomp_b, FNVINITOMP_B)
#define FNV_INITOMP_QB SUNDIALS_F77_FUNC_(fnvinitomp_qb, FNVINITOMP_QB)

#else

#define FNV_INITOMP_Q  fnvinitomp_q_
#define FNV_INITOMP_S  fnvinitomp_s_
#define FNV_INITOMP_B  fnvinitomp_b_
#define FNV_INITOMP_QB fnvinitomp_qb_

#endif

/* Declarations of global variables */

extern N_Vector F2C_CVODE_vec;
extern N_Vector F2C_CVODE_vecQ;
extern N_Vector *F2C_CVODE_vecS;
extern N_Vector F2C_CVODE_vecB;
extern N_Vector F2C_CVODE_vecQB;

extern N_Vector F2C_IDA_vec;
extern N_Vector F2C_IDA_vecQ;
extern N_Vector *F2C_IDA_vecS;
extern N_Vector F2C_IDA_vecB;
extern N_Vector F2C_IDA_vecQB;

extern N_Vector F2C_KINSOL_vec;

extern N_Vector F2C_ARKODE_vec;

/*
 * Prototypes of exported functions
 *
 * FNV_INITOMP    - initializes openmp vector operations for main problem
 * FNV_INITOMP_Q  - initializes openmp vector operations for quadratures
 * FNV_INITOMP_S  - initializes openmp vector operations for sensitivities
 * FNV_INITOMP_B  - initializes openmp vector operations for adjoint problem
 * FNV_INITOMP_QB - initializes openmp vector operations for adjoint quadratures
 *
 */

void FNV_INITOMP(int *code, long int *neq, int *num_threads, int *ier);
void FNV_INITOMP_Q(int *code, long int *Nq, int *num_threads, int *ier);
void FNV_INITOMP_S(int *code, int *Ns, int *ier);
void FNV_INITOMP_B(int *code, long int *NB, int *num_threads, int *ier);
void FNV_INITOMP_QB(int *code, long int *NqB, int *num_threads, int *ier);

#ifdef __cplusplus
}
#endif

#endif
