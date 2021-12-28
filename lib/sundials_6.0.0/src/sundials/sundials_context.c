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
 * SUNDIALS context class. A context object holds data that all
 * SUNDIALS objects in a simulation share.
 * ----------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sundials/sundials_profiler.h>
#include <sundials/sundials_context.h>
#include "sundials_context_impl.h"
#include "sundials_debug.h"


int SUNContext_Create(void* comm, SUNContext* sunctx)
{
  SUNProfiler profiler = NULL;

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  if (SUNProfiler_Create(comm, "SUNContext Default", &profiler))
    return(-1);
#endif

  *sunctx = NULL;
  *sunctx = (SUNContext) malloc(sizeof(struct _SUNContext));

  if (*sunctx == NULL)
  {
#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
    SUNProfiler_Free(&profiler);
#endif
    return(-1);
  }

  (*sunctx)->profiler = profiler;
  (*sunctx)->own_profiler = SUNTRUE;

  return(0);
}

int SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler)
{
  if (sunctx == NULL) return(-1);

#ifdef SUNDIALS_BUILD_WITH_PROFILING
  /* get profiler */
  *profiler = sunctx->profiler;
#else
  *profiler = NULL;
#endif

  return(0);
}

int SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler)
{
  if (sunctx == NULL) return(-1);

#ifdef SUNDIALS_BUILD_WITH_PROFILING
  /* free any existing profiler */
  if (sunctx->profiler && sunctx->own_profiler) {
    if (SUNProfiler_Free(&(sunctx->profiler)))
      return(-1);
    sunctx->profiler = NULL;
  }

  /* set profiler */
  sunctx->profiler = profiler;
  sunctx->own_profiler = SUNFALSE;
#endif

  return(0);
}

int SUNContext_Free(SUNContext* sunctx)
{
#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  FILE *fp;
  char* sunprofiler_print_env;
#endif

  if (!sunctx) return(0);
  if (!(*sunctx)) return(0);

#if defined(SUNDIALS_BUILD_WITH_PROFILING) && !defined(SUNDIALS_CALIPER_ENABLED)
  /* Find out where we are printing to */
  sunprofiler_print_env = getenv("SUNPROFILER_PRINT");
  fp = NULL;
  if (sunprofiler_print_env)
  {
    if (!strcmp(sunprofiler_print_env, "0"))
      fp = NULL;
    else if (!strcmp(sunprofiler_print_env, "1") ||
             !strcmp(sunprofiler_print_env, "TRUE") ||
             !strcmp(sunprofiler_print_env, "stdout"))
      fp = stdout;
    else
      fp = fopen(sunprofiler_print_env, "a");
  }

  /* Enforce that the profiler is freed before finalizing,
     if it is not owned by the sunctx. */
  if ((*sunctx)->profiler)
  {
    if (fp) SUNProfiler_Print((*sunctx)->profiler, fp);
    if (fp) fclose(fp);
    if ((*sunctx)->own_profiler && SUNProfiler_Free(&(*sunctx)->profiler))
      return(-1);
  }
#endif

  free(*sunctx);
  *sunctx = NULL;

  return(0);
}
