/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
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
 * This is the header file for the ARKBANDPRE module, which provides
 * a banded difference quotient Jacobian-based preconditioner.
 * -----------------------------------------------------------------*/

#ifndef _ARKBANDPRE_H
#define _ARKBANDPRE_H

#include <sundials/sundials_nvector.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* BandPrec inititialization function */

SUNDIALS_EXPORT int ARKBandPrecInit(void *arkode_mem, sunindextype N,
                                    sunindextype mu, sunindextype ml);

/* Optional output functions */

SUNDIALS_EXPORT int ARKBandPrecGetWorkSpace(void *arkode_mem,
                                            long int *lenrwLS,
                                            long int *leniwLS);
SUNDIALS_EXPORT int ARKBandPrecGetNumRhsEvals(void *arkode_mem,
                                              long int *nfevalsBP);


#ifdef __cplusplus
}
#endif

#endif
