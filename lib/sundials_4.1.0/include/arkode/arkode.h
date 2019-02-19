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
 * This is the header file for the main ARKode infrastructure.
 * -----------------------------------------------------------------
 * ARKode is used to numerically solve the ordinary initial value
 * problems using one-step methods.  Users do not call ARKode
 * infrastructure routines directly; they instead interact with
 * one of the time stepping modules built on top of ARKode.
 * These time step modules define their supported problem types,
 * solver options, etc.
 *
 * This file serves to define constants and provide function
 * prototypes for use across ARKode-based time integration
 * modules.
 * -----------------------------------------------------------------*/

#ifndef _ARKODE_H
#define _ARKODE_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>
#include <arkode/arkode_butcher.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------
 * ARKode Constants
 * ----------------- */

/* itask */
#define ARK_NORMAL         1
#define ARK_ONE_STEP       2

/* return values */

#define ARK_SUCCESS               0
#define ARK_TSTOP_RETURN          1
#define ARK_ROOT_RETURN           2

#define ARK_WARNING              99

#define ARK_TOO_MUCH_WORK        -1
#define ARK_TOO_MUCH_ACC         -2
#define ARK_ERR_FAILURE          -3
#define ARK_CONV_FAILURE         -4

#define ARK_LINIT_FAIL           -5
#define ARK_LSETUP_FAIL          -6
#define ARK_LSOLVE_FAIL          -7
#define ARK_RHSFUNC_FAIL         -8
#define ARK_FIRST_RHSFUNC_ERR    -9
#define ARK_REPTD_RHSFUNC_ERR    -10
#define ARK_UNREC_RHSFUNC_ERR    -11
#define ARK_RTFUNC_FAIL          -12
#define ARK_LFREE_FAIL           -13
#define ARK_MASSINIT_FAIL        -14
#define ARK_MASSSETUP_FAIL       -15
#define ARK_MASSSOLVE_FAIL       -16
#define ARK_MASSFREE_FAIL        -17
#define ARK_MASSMULT_FAIL        -18

#define ARK_MEM_FAIL             -20
#define ARK_MEM_NULL             -21
#define ARK_ILL_INPUT            -22
#define ARK_NO_MALLOC            -23
#define ARK_BAD_K                -24
#define ARK_BAD_T                -25
#define ARK_BAD_DKY              -26
#define ARK_TOO_CLOSE            -27

#define ARK_POSTPROCESS_FAIL     -28
#define ARK_VECTOROP_ERR         -29

#define ARK_NLS_INIT_FAIL        -30
#define ARK_NLS_SETUP_FAIL       -31
#define ARK_NLS_SETUP_RECVR      -32
#define ARK_NLS_OP_ERR           -33

#define ARK_INNERSTEP_FAIL       -34

#define ARK_UNRECOGNIZED_ERROR   -99

/* ------------------------------
 * User-Supplied Function Types
 * ------------------------------ */

typedef int (*ARKRhsFn)(realtype t, N_Vector y,
                        N_Vector ydot, void *user_data);

typedef int (*ARKRootFn)(realtype t, N_Vector y,
                         realtype *gout, void *user_data);

typedef int (*ARKEwtFn)(N_Vector y, N_Vector ewt, void *user_data);

typedef int (*ARKRwtFn)(N_Vector y, N_Vector rwt, void *user_data);

typedef void (*ARKErrHandlerFn)(int error_code, const char *module,
                                const char *function, char *msg,
                                void *user_data);

typedef int (*ARKAdaptFn)(N_Vector y, realtype t, realtype h1,
                          realtype h2, realtype h3,
                          realtype e1, realtype e2,
                          realtype e3, int q, int p,
                          realtype *hnew, void *user_data);

typedef int (*ARKExpStabFn)(N_Vector y, realtype t,
                            realtype *hstab, void *user_data);

typedef int (*ARKVecResizeFn)(N_Vector y, N_Vector ytemplate,
                              void *user_data);

typedef int (*ARKPostProcessStepFn)(realtype t, N_Vector y,
                                    void *user_data);


#ifdef __cplusplus
}
#endif

#endif
