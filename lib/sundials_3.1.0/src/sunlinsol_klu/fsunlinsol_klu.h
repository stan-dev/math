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
 * This file (companion of fsunlinsol_klu.c) contains the
 * definitions needed for the initialization of klu
 * linear solver operations in Fortran.
 * -----------------------------------------------------------------
 */

#ifndef _FSUNLINSOL_KLU_H
#define _FSUNLINSOL_KLU_H

#include <sunlinsol/sunlinsol_klu.h>
#include <sundials/sundials_fnvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNKLU_INIT            SUNDIALS_F77_FUNC(fsunkluinit,        FSUNKLUINIT)
#define FSUNKLU_REINIT          SUNDIALS_F77_FUNC(fsunklureinit,      FSUNKLUREINIT)
#define FSUNKLU_SETORDERING     SUNDIALS_F77_FUNC(fsunklusetordering, FSUNKLUSETORDERING)
#define FSUNMASSKLU_INIT        SUNDIALS_F77_FUNC(fsunmasskluinit,        FSUNMASSKLUINIT)
#define FSUNMASSKLU_REINIT      SUNDIALS_F77_FUNC(fsunmassklureinit,      FSUNMASSKLUREINIT)
#define FSUNMASSKLU_SETORDERING SUNDIALS_F77_FUNC(fsunmassklusetordering, FSUNMASSKLUSETORDERING)
#else
#define FSUNKLU_INIT            fsunkluinit_
#define FSUNKLU_REINIT          fsunklureinit_
#define FSUNKLU_SETORDERING     fsunklusetordering_
#define FSUNMASSKLU_INIT        fsunmasskluinit_
#define FSUNMASSKLU_REINIT      fsunmassklureinit_
#define FSUNMASSKLU_SETORDERING fsunmassklusetordering_
#endif


/* Declarations of global variables */

extern SUNLinearSolver F2C_CVODE_linsol;
extern SUNLinearSolver F2C_IDA_linsol;
extern SUNLinearSolver F2C_KINSOL_linsol;
extern SUNLinearSolver F2C_ARKODE_linsol;
extern SUNLinearSolver F2C_ARKODE_mass_sol;

/* 
 * Prototypes of exported functions 
 *
 * FSUNKLU_INIT - initializes klu linear solver for main problem
 * FSUNKLU_REINIT - reinitializes klu linear solver for main problem
 * FSUNKLU_SETORDERING - sets the ordering choice used by KLU for main problem
 * FSUNMASSKLU_INIT - initializes klu linear solver for mass matrix solve
 * FSUNMASSKLU_REINIT - reinitializes klu linear solver for mass matrix solve
 * FSUNMASSKLU_SETORDERING - sets the ordering choice used by KLU for mass matrix solve
 */

void FSUNKLU_INIT(int *code, int *ier);
void FSUNKLU_REINIT(int *code, long int *NNZ, 
                    int *reinit_type, int *ier);
void FSUNKLU_SETORDERING(int *code, int *ordering,
                         int *ier);
void FSUNMASSKLU_INIT(int *ier);
void FSUNMASSKLU_REINIT(long int *NNZ, 
                        int *reinit_type, int *ier);
void FSUNMASSKLU_SETORDERING(int *ordering, int *ier);

#ifdef __cplusplus
}
#endif

#endif
