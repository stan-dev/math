/*-----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *-----------------------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and
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
 *-----------------------------------------------------------------------------
 * This file contains the definitions needed for initialization of the
 * SUNNonlinearSolver fixed-point moudule operations in Fortran.
 *---------------------------------------------------------------------------*/

#ifndef _FSUNNONLINSOL_FIXEDPOINT_H
#define _FSUNNONLINSOL_FIXEDPOINT_H

#include <sundials/sundials_fnvector.h>       /* FCMIX_* solver IDs */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#if defined(SUNDIALS_F77_FUNC)
#define FSUNFIXEDPOINT_INIT        SUNDIALS_F77_FUNC(fsunfixedpointinit,        FSUNFIXEDPOINTINIT)
#define FSUNFIXEDPOINT_SETMAXITERS SUNDIALS_F77_FUNC(fsunfixedpointsetmaxiters, FSUNFIXEDPOINTSETMAXITERS)
#else
#define FSUNFIXEDPOINT_INIT        fsunfixedpointinit_
#define FSUNFIXEDPOINT_SETMAXITERS fsunfixedpointsetmaxiters_
#endif

/* Declarations of global variables */

extern SUNNonlinearSolver F2C_CVODE_nonlinsol;
extern SUNNonlinearSolver F2C_IDA_nonlinsol;
extern SUNNonlinearSolver F2C_ARKODE_nonlinsol;

/*-----------------------------------------------------------------------------
  Prototypes of exported functions

  FSUNFIXEDPOINT_INIT - initializes fixed point nonlinear solver for main problem
  FSUNFIXEDPOINT_SETMAXITERS - sets the maximum number of nonlinear iterations
  ---------------------------------------------------------------------------*/

void FSUNFIXEDPOINT_INIT(int *code, int *m, int *ier);
void FSUNFIXEDPOINT_SETMAXITERS(int *code, int *maxiters, int *ier);

#ifdef __cplusplus
}
#endif

#endif
