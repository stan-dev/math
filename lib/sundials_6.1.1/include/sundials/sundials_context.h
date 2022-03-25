/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
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
 * SUNDIALS context class. A context object holds data that all
 * SUNDIALS objects in a simulation share. It is thread-safe provided
 * that each thread has its own context object.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_CONTEXT_H
#define _SUNDIALS_CONTEXT_H

#include "sundials/sundials_types.h"
#include "sundials/sundials_profiler.h"

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef struct _SUNContext *SUNContext;

SUNDIALS_EXPORT int SUNContext_Create(void* comm, SUNContext* ctx);
SUNDIALS_EXPORT int SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler);
SUNDIALS_EXPORT int SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler);
SUNDIALS_EXPORT int SUNContext_Free(SUNContext* ctx);

#ifdef __cplusplus
}

namespace sundials
{

class Context
{
public:
   Context(void* comm = NULL)
   {
      SUNContext_Create(comm, &sunctx_);
   }

   operator SUNContext() { return sunctx_; }

   ~Context()
   {
      SUNContext_Free(&sunctx_);
   }

private:
   SUNContext sunctx_;

};

} /* namespace sundials */

#endif
#endif
