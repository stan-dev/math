/*
 * -----------------------------------------------------------------
 * Programmer(s): Steven Smith @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * This file (companion of nvector_pthreads.h) contains the
 * definitions needed for the initialization of pthreads
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_PTHREADS_H
#define _FNVECTOR_PTHREADS_H

#include <nvector/nvector_pthreads.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FNV_INITPTS    SUNDIALS_F77_FUNC(fnvinitpts, FNVINITPTS)
#else
#define FNV_INITPTS    fnvinitpts_
#endif

#if defined(SUNDIALS_F77_FUNC_)

#define FNV_INITPTS_Q  SUNDIALS_F77_FUNC_(fnvinitpts_q, FNVINITPTS_Q)
#define FNV_INITPTS_S  SUNDIALS_F77_FUNC_(fnvinitpts_s, FNVINITPTS_S)
#define FNV_INITPTS_B  SUNDIALS_F77_FUNC_(fnvinitpts_b, FNVINITPTS_B)
#define FNV_INITPTS_QB SUNDIALS_F77_FUNC_(fnvinitpts_qb, FNVINITPTS_QB)

#else

#define FNV_INITPTS_Q  fnvinitpts_q_
#define FNV_INITPTS_S  fnvinitpts_s_
#define FNV_INITPTS_B  fnvinitpts_b_
#define FNV_INITPTS_QB fnvinitpts_qb_

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
 * FNV_INITPTS    - initializes pthreads vector operations for main problem
 * FNV_INITPTS_Q  - initializes pthreads vector operations for quadratures
 * FNV_INITPTS_S  - initializes pthreads vector operations for sensitivities
 * FNV_INITPTS_B  - initializes pthreads vector operations for adjoint problem
 * FNV_INITPTS_QB - initializes pthreads vector operations for adjoint quadratures
 *
 */

void FNV_INITPTS(int *code, long int *neq, int *num_threads, int *ier);
void FNV_INITPTS_Q(int *code, long int *Nq, int *num_threads, int *ier);
void FNV_INITPTS_S(int *code, int *Ns, int *ier);
void FNV_INITPTS_B(int *code, long int *NB, int *num_threads, int *ier);
void FNV_INITPTS_QB(int *code, long int *NqB, int *num_threads, int *ier);

#ifdef __cplusplus
}
#endif

#endif
