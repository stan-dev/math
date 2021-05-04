/*
 * Copyright 2006. Rene Rivera
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#ifndef BJAM_MEM_H
#define BJAM_MEM_H

#include "config.h"

/* Standard C memory allocation. */
#include <stdlib.h>

#define bjam_malloc_x(s) malloc(s)
#define bjam_calloc_x(n,s) calloc(n,s)
#define bjam_realloc_x(p,s) realloc(p,s)
#define bjam_free_x(p) free(p)

#ifndef bjam_malloc_atomic_x
    #define bjam_malloc_atomic_x(s) bjam_malloc_x(s)
#endif
#ifndef bjam_calloc_atomic_x
    #define bjam_calloc_atomic_x(n,s) bjam_calloc_x(n,s)
#endif
#ifndef bjam_mem_init_x
    #define bjam_mem_init_x()
#endif
#ifndef bjam_mem_close_x
    #define bjam_mem_close_x()
#endif
#ifndef bjam_malloc_raw_x
    #define bjam_malloc_raw_x(s) bjam_malloc_x(s)
#endif
#ifndef bjam_calloc_raw_x
    #define bjam_calloc_raw_x(n,s) bjam_calloc_x(n,s)
#endif
#ifndef bjam_realloc_raw_x
    #define bjam_realloc_raw_x(p,s) bjam_realloc_x(p,s)
#endif
#ifndef bjam_free_raw_x
    #define bjam_free_raw_x(p) bjam_free_x(p)
#endif

#ifdef OPT_DEBUG_PROFILE
    /* Profile tracing of memory allocations. */
    #include "debug.h"

    #define BJAM_MALLOC(s) (profile_memory(s), bjam_malloc_x(s))
    #define BJAM_MALLOC_ATOMIC(s) (profile_memory(s), bjam_malloc_atomic_x(s))
    #define BJAM_CALLOC(n,s) (profile_memory(n*s), bjam_calloc_x(n,s))
    #define BJAM_CALLOC_ATOMIC(n,s) (profile_memory(n*s), bjam_calloc_atomic_x(n,s))
    #define BJAM_REALLOC(p,s) (profile_memory(s), bjam_realloc_x(p,s))

    #define BJAM_MALLOC_RAW(s) (profile_memory(s), bjam_malloc_raw_x(s))
    #define BJAM_CALLOC_RAW(n,s) (profile_memory(n*s), bjam_calloc_raw_x(n,s))
    #define BJAM_REALLOC_RAW(p,s) (profile_memory(s), bjam_realloc_raw_x(p,s))
#else
    /* No mem tracing. */
    #define BJAM_MALLOC(s) bjam_malloc_x(s)
    #define BJAM_MALLOC_ATOMIC(s) bjam_malloc_atomic_x(s)
    #define BJAM_CALLOC(n,s) bjam_calloc_x(n,s)
    #define BJAM_CALLOC_ATOMIC(n,s) bjam_calloc_atomic_x(n,s)
    #define BJAM_REALLOC(p,s) bjam_realloc_x(p,s)

    #define BJAM_MALLOC_RAW(s) bjam_malloc_raw_x(s)
    #define BJAM_CALLOC_RAW(n,s) bjam_calloc_raw_x(n,s)
    #define BJAM_REALLOC_RAW(p,s) bjam_realloc_raw_x(p,s)
#endif

#define BJAM_MEM_INIT() bjam_mem_init_x()
#define BJAM_MEM_CLOSE() bjam_mem_close_x()

#define BJAM_FREE(p) bjam_free_x(p)
#define BJAM_FREE_RAW(p) bjam_free_raw_x(p)

#endif
