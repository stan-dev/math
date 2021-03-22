/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * This is the header file for the ARKODE + XBraid interface.
 * ---------------------------------------------------------------------------*/

#ifndef _ARKODE_XBRAID_H
#define _ARKODE_XBRAID_H

#include "sundials/sundials_xbraid.h"
#include "braid.h"
#include "mpi.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif


/* -------------------------------
 * Construct, initialize, and free
 * ------------------------------- */


SUNDIALS_EXPORT int ARKBraid_Create(void *arkode_mem, braid_App *app);

SUNDIALS_EXPORT int ARKBraid_BraidInit(MPI_Comm comm_w, MPI_Comm comm_t,
                                       realtype tstart, realtype tstop,
                                       sunindextype nt, braid_App app,
                                       braid_Core *core);

SUNDIALS_EXPORT int ARKBraid_Free(braid_App *app);


/* ----------------------
 * ARKBraid Set Functions
 * ---------------------- */


SUNDIALS_EXPORT int ARKBraid_SetStepFn(braid_App app, braid_PtFcnStep step);

SUNDIALS_EXPORT int ARKBraid_SetInitFn(braid_App app, braid_PtFcnInit init);

SUNDIALS_EXPORT int ARKBraid_SetSpatialNormFn(braid_App app,
                                              braid_PtFcnSpatialNorm snorm);

SUNDIALS_EXPORT int ARKBraid_SetAccessFn(braid_App app,
                                         braid_PtFcnAccess access);


/* ----------------------
 * ARKBraid Get Functions
 * ---------------------- */


SUNDIALS_EXPORT int ARKBraid_GetVecTmpl(braid_App app, N_Vector *tmpl);

SUNDIALS_EXPORT int ARKBraid_GetARKStepMem(braid_App app, void **arkode_mem);

SUNDIALS_EXPORT int ARKBraid_GetUserData(braid_App app, void **user_data);

SUNDIALS_EXPORT int ARKBraid_GetLastBraidFlag(braid_App app, int *last_flag);

SUNDIALS_EXPORT int ARKBraid_GetLastARKStepFlag(braid_App app, int *last_flag);

SUNDIALS_EXPORT int ARKBraid_GetSolution(braid_App app, realtype *tout,
                                         N_Vector yout);


/* --------------------------
 * XBraid Interface Functions
 * -------------------------- */


SUNDIALS_EXPORT int ARKBraid_Step(braid_App app, braid_Vector ustop,
                                  braid_Vector fstop, braid_Vector u,
                                  braid_StepStatus status);

SUNDIALS_EXPORT int ARKBraid_Init(braid_App app, realtype t,
                                  braid_Vector *u_ptr);

SUNDIALS_EXPORT int ARKBraid_Access(braid_App app, braid_Vector u,
                                    braid_AccessStatus astatus);


/* -----------------
 * Utility Functions
 * ----------------- */


SUNDIALS_EXPORT int ARKBraid_TakeStep(void *arkode_mem, realtype tstart,
                                      realtype tstop, N_Vector y,
                                      int *ark_flag);


#ifdef __cplusplus
}
#endif

#endif
