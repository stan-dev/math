/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * SUNDIALS context class implementation.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_CONTEXT_IMPL_H
#define _SUNDIALS_CONTEXT_IMPL_H

#include <sundials/sundials_types.h>
#include <sundials/sundials_profiler.h>
#include <sundials/sundials_context.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

struct _SUNContext {
   SUNProfiler profiler;
   booleantype own_profiler;
};

#ifdef __cplusplus
}
#endif
#endif
