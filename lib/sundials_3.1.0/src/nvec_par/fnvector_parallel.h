/*
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
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
 * This file (companion of nvector_parallel.c) contains the
 * definitions needed for the initialization of parallel
 * vector operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FNVECTOR_PARALLEL_H
#define _FNVECTOR_PARALLEL_H

#include <nvector/nvector_parallel.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FNV_INITP    SUNDIALS_F77_FUNC(fnvinitp, FNVINITP)
#else
#define FNV_INITP    fnvinitp_
#endif

#if defined(SUNDIALS_F77_FUNC_)

#define FNV_INITP_Q  SUNDIALS_F77_FUNC_(fnvinitp_q, FNVINITP_Q)
#define FNV_INITP_S  SUNDIALS_F77_FUNC_(fnvinitp_s, FNVINITP_S)
#define FNV_INITP_B  SUNDIALS_F77_FUNC_(fnvinitp_b, FNVINITP_B)
#define FNV_INITP_QB SUNDIALS_F77_FUNC_(fnvinitp_qb, FNVINITP_QB)

#else

#define FNV_INITP_Q  fnvinitp_q_
#define FNV_INITP_S  fnvinitp_s_
#define FNV_INITP_B  fnvinitp_b_
#define FNV_INITP_QB fnvinitp_qb_

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
 * FNV_INITP    - initializes parallel vector operations for main problem
 * FNV_INITP_Q  - initializes parallel vector operations for quadratures
 * FNV_INITP_S  - initializes parallel vector operations for sensitivities
 * FNV_INITP_B  - initializes parallel vector operations for adjoint problem
 * FNV_INITP_QB - initializes parallel vector operations for adjoint quadratures
 *
 */

#ifndef SUNDIALS_MPI_COMM_F2C
#define MPI_Fint int
#endif

void FNV_INITP(MPI_Fint *comm, int *code, long int *L, long int *N, int *ier);
void FNV_INITP_Q(MPI_Fint *comm, int *code, long int *Lq, long int *Nq, int *ier);
void FNV_INITP_B(MPI_Fint *comm, int *code, long int *LB, long int *NB, int *ier);
void FNV_INITP_QB(MPI_Fint *comm, int *code, long int *LqB, long int *NqB, int *ier);
void FNV_INITP_S(int *code, int *Ns, int *ier);

#ifdef __cplusplus
}
#endif

#endif
