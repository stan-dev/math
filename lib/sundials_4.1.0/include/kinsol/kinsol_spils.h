/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Scott Cohen, Alan Hindmarsh, Radu Serban,
 *                  and Aaron Collier @ LLNL
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
 * Header file for the deprecated Scaled Preconditioned Iterative
 * Linear Solver interface in KINSOL; these routines now just wrap
 * the updated KINSOL generic linear solver interface in kinsol_ls.h.
 * -----------------------------------------------------------------*/

#ifndef _KINSPILS_H
#define _KINSPILS_H

#include <kinsol/kinsol_ls.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/*===============================================================
  Function Types (typedefs for equivalent types in kinsol_ls.h)
  ===============================================================*/

typedef KINLsPrecSetupFn KINSpilsPrecSetupFn;
typedef KINLsPrecSolveFn KINSpilsPrecSolveFn;
typedef KINLsJacTimesVecFn KINSpilsJacTimesVecFn;

/*====================================================================
  Exported Functions (wrappers for equivalent routines in kinsol_ls.h)
  ====================================================================*/

int KINSpilsSetLinearSolver(void *kinmem, SUNLinearSolver LS);

int KINSpilsSetPreconditioner(void *kinmem, KINSpilsPrecSetupFn psetup,
                              KINSpilsPrecSolveFn psolve);

int KINSpilsSetJacTimesVecFn(void *kinmem, KINSpilsJacTimesVecFn jtv);

int KINSpilsGetWorkSpace(void *kinmem, long int *lenrwLS, long int *leniwLS);

int KINSpilsGetNumPrecEvals(void *kinmem, long int *npevals);

int KINSpilsGetNumPrecSolves(void *kinmem, long int *npsolves);

int KINSpilsGetNumLinIters(void *kinmem, long int *nliters);

int KINSpilsGetNumConvFails(void *kinmem, long int *nlcfails);

int KINSpilsGetNumJtimesEvals(void *kinmem, long int *njvevals);

int KINSpilsGetNumFuncEvals(void *kinmem, long int *nfevals);

int KINSpilsGetLastFlag(void *kinmem, long int *flag);

char *KINSpilsGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif
