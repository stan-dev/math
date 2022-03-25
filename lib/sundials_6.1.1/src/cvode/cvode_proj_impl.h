/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * Based on CPODES by Radu Serban @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Implementation header file for projections in CVODE.
 * ---------------------------------------------------------------------------*/

#ifndef _CVODE_PROJ_IMPL_H
#define _CVODE_PROJ_IMPL_H

#include "cvode/cvode.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* =============================================================================
 * Default Projection Constants
 *
 * PROJ_MAX_FAILS  max nunmber of projection failures in one step attempt
 * PROJ_EPS        projection solve tolerance
 * PROJ_FAIL_ETA   maximum step size decrease on projection failure
 * ===========================================================================*/

#define PROJ_MAX_FAILS 10
#define PROJ_EPS       RCONST(0.1)
#define PROJ_FAIL_ETA  RCONST(0.25)

/* =============================================================================
 * Projection Data Structure
 * ===========================================================================*/

/* -----------------------------------------------------------------------------
 * Types : struct CVodeProjMemRec, CVodeProjMem
 * -----------------------------------------------------------------------------
 * The type CVodeProjMem is type pointer to struct CVodeProjMemRec. This
 * structure contains data pertaining to the use of projection capabilities.
 * ---------------------------------------------------------------------------*/
typedef struct CVodeProjMemRec {

  booleantype internal_proj;  /* use the internal projection algorithm?      */
  booleantype err_proj;       /* is error projection enabled?                */
  booleantype first_proj;     /* is this the first time we project?          */

  long int freq;              /* projection frequency                        */
  long int nstlprj;           /* step number of last projection              */

  int max_fails;              /* maximum number of projection failures       */

  CVProjFn pfun;              /* function to perform projection              */

  realtype eps_proj;          /* projection solve tolerance                  */
  realtype eta_pfail;         /* projection failure step reduction factor    */

  long int nproj;             /* number of projections performed             */
  long int npfails;           /* number of projection failures               */

} *CVodeProjMem;

#ifdef __cplusplus
}
#endif

#endif
