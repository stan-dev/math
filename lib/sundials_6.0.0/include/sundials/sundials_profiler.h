/* -----------------------------------------------------------------
 * Programmer: Cody J. Balos @ LLNL
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
 * -----------------------------------------------------------------*/

#ifndef _SUNDIALS_PROFILER_H
#define _SUNDIALS_PROFILER_H

#include <stdio.h>

#include "sundials/sundials_config.h"

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && defined(SUNDIALS_CALIPER_ENABLED)
#include "caliper/cali.h"
#endif

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef struct _SUNProfiler *SUNProfiler;

SUNDIALS_EXPORT int SUNProfiler_Create(void* comm, const char* title, SUNProfiler* p);
SUNDIALS_EXPORT int SUNProfiler_Free(SUNProfiler* p);
SUNDIALS_EXPORT int SUNProfiler_Begin(SUNProfiler p, const char* name);
SUNDIALS_EXPORT int SUNProfiler_End(SUNProfiler p, const char* name);
SUNDIALS_EXPORT int SUNProfiler_Print(SUNProfiler p, FILE* fp);

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && defined(SUNDIALS_CALIPER_ENABLED)

#define SUNDIALS_MARK_FUNCTION_BEGIN(profobj) CALI_MARK_FUNCTION_BEGIN

#define SUNDIALS_MARK_FUNCTION_END(profobj) CALI_MARK_FUNCTION_END

#define SUNDIALS_WRAP_STATEMENT(profobj, name, stmt) CALI_WRAP_STATEMENT(name, stmt)

#define SUNDIALS_MARK_BEGIN(profobj, name) CALI_MARK_BEGIN(name)

#define SUNDIALS_MARK_END(profobj, name) CALI_MARK_END(name)

#ifdef __cplusplus
#define SUNDIALS_CXX_MARK_FUNCTION(projobj) CALI_CXX_MARK_FUNCTION
#endif

#elif defined(SUNDIALS_BUILD_WITH_PROFILING)

#define SUNDIALS_MARK_FUNCTION_BEGIN(profobj) SUNProfiler_Begin(profobj, __func__)

#define SUNDIALS_MARK_FUNCTION_END(profobj) SUNProfiler_End(profobj, __func__)

#define SUNDIALS_WRAP_STATEMENT(profobj, name, stmt) \
    SUNProfiler_Begin(profobj, (name)); \
    stmt; \
    SUNProfiler_End(profobj, (name));

#define SUNDIALS_MARK_BEGIN(profobj, name) SUNProfiler_Begin(profobj, (name))

#define SUNDIALS_MARK_END(profobj, name) SUNProfiler_End(profobj, (name))

#ifdef __cplusplus
#define SUNDIALS_CXX_MARK_FUNCTION(profobj) sundials::ProfilerMarkScope __ProfilerMarkScope(profobj, __func__)
#endif

#else

#define SUNDIALS_MARK_FUNCTION_BEGIN(profobj)

#define SUNDIALS_MARK_FUNCTION_END(profobj)

#define SUNDIALS_WRAP_STATEMENT(profobj, name, stmt)

#define SUNDIALS_MARK_BEGIN(profobj, name)

#define SUNDIALS_MARK_END(profobj, name)

#ifdef __cplusplus
#define SUNDIALS_CXX_MARK_FUNCTION(profobj)
#endif

#endif

#ifdef __cplusplus
}

namespace sundials
{
/* Convenience class for C++ codes.
   Allows for simpler profiler statements using C++ scoping rules. */
class ProfilerMarkScope
{
public:
  ProfilerMarkScope(SUNProfiler prof, const char* name) {
    prof_ = prof;
    name_ = name;
    SUNProfiler_Begin(prof_, name_);
  }

  ~ProfilerMarkScope() {
    SUNProfiler_End(prof_, name_);
  }
private:
  SUNProfiler prof_;
  const char* name_;
};
}

#endif
#endif /* SUNDIALS_PROFILER_H_ */
