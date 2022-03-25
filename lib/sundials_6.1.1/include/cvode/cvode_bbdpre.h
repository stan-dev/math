/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Michael Wittman, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the header file for the CVBBDPRE module, for a
 * band-block-diagonal preconditioner, i.e. a block-diagonal
 * matrix with banded blocks.
 * -----------------------------------------------------------------*/

#ifndef _CVBBDPRE_H
#define _CVBBDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* User-supplied function Types */

typedef int (*CVLocalFn)(sunindextype Nlocal, realtype t,
                         N_Vector y, N_Vector g, void *user_data);

typedef int (*CVCommFn)(sunindextype Nlocal, realtype t,
                        N_Vector y, void *user_data);

/* Exported Functions */

SUNDIALS_EXPORT int CVBBDPrecInit(void *cvode_mem, sunindextype Nlocal,
                                  sunindextype mudq, sunindextype mldq,
                                  sunindextype mukeep, sunindextype mlkeep,
                                  realtype dqrely, CVLocalFn gloc, CVCommFn cfn);

SUNDIALS_EXPORT int CVBBDPrecReInit(void *cvode_mem,
                                    sunindextype mudq, sunindextype mldq,
                                    realtype dqrely);


/* Optional output functions */

SUNDIALS_EXPORT int CVBBDPrecGetWorkSpace(void *cvode_mem,
                                          long int *lenrwBBDP,
                                          long int *leniwBBDP);

SUNDIALS_EXPORT int CVBBDPrecGetNumGfnEvals(void *cvode_mem,
                                            long int *ngevalsBBDP);


#ifdef __cplusplus
}
#endif

#endif
