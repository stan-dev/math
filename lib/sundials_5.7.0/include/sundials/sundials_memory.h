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
 * SUNDIALS memory helpers and types.
 * ----------------------------------------------------------------*/

#ifndef _SUNDIALS_MEMORY_H
#define _SUNDIALS_MEMORY_H

#include <stdlib.h>

#include <sundials/sundials_types.h>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

typedef enum
{
  SUNMEMTYPE_HOST,      /* pageable memory accessible on the host     */
  SUNMEMTYPE_PINNED,    /* page-locked memory accesible on the host   */
  SUNMEMTYPE_DEVICE,    /* memory accessible from the device          */
  SUNMEMTYPE_UVM        /* memory accessible from the host or device  */
} SUNMemoryType;


/*
 * SUNMemory is a simple abstraction of a pointer to some
 * contiguos memory, so that we can keep track of its type
 * and its ownership.
 */

typedef struct _SUNMemory *SUNMemory;

struct _SUNMemory
{
  void*         ptr;
  SUNMemoryType type;
  booleantype   own;
};

/* Creates a new SUNMemory object with a NULL ptr */
SUNDIALS_EXPORT SUNMemory SUNMemoryNewEmpty();

/*
 * SUNMemoryHelper holds ops which can allocate, deallocate,
 * and copy SUNMemory.
 */

typedef struct _SUNMemoryHelper_Ops *SUNMemoryHelper_Ops;
typedef struct _SUNMemoryHelper *SUNMemoryHelper;

struct _SUNMemoryHelper
{
  void*               content;
  SUNMemoryHelper_Ops ops;
};

struct _SUNMemoryHelper_Ops
{
  /* operations that implementations are required to provide */
  int             (*alloc)(SUNMemoryHelper, SUNMemory* memptr, size_t mem_size, SUNMemoryType mem_type);
  int             (*dealloc)(SUNMemoryHelper, SUNMemory mem);
  int             (*copy)(SUNMemoryHelper, SUNMemory dst, SUNMemory src, size_t mem_size);

  /* operations that provide default implementations */
  int             (*copyasync)(SUNMemoryHelper, SUNMemory dst, SUNMemory src,
                               size_t mem_size, void* ctx);
  SUNMemoryHelper (*clone)(SUNMemoryHelper);
  int             (*destroy)(SUNMemoryHelper);
};


/*
 * Generic SUNMemoryHelper functions that work without a SUNMemoryHelper object.
 */

/* Creates a new SUNMemory object which points to the same data as another
 * SUNMemory object.
 * The SUNMemory returned will not own the ptr, therefore, it will not free
 * the ptr in Dealloc. */
SUNDIALS_EXPORT SUNMemory SUNMemoryHelper_Alias(SUNMemory mem);

/* Creates a new SUNMemory object with ptr set to the user provided pointer
 * The SUNMemory returned will not own the ptr, therefore, it will not free
 * the ptr in Dealloc. */
SUNDIALS_EXPORT SUNMemory SUNMemoryHelper_Wrap(void* ptr, SUNMemoryType mem_type);


/*
 * Required SUNMemoryHelper operations.
 */


SUNDIALS_EXPORT int SUNMemoryHelper_Alloc(SUNMemoryHelper, SUNMemory* memptr,
                                          size_t mem_size, SUNMemoryType mem_type);

SUNDIALS_EXPORT int SUNMemoryHelper_Dealloc(SUNMemoryHelper, SUNMemory mem);

SUNDIALS_EXPORT int SUNMemoryHelper_Copy(SUNMemoryHelper, SUNMemory dst,
                                         SUNMemory src, size_t mem_size);

/*
 * Optional SUNMemoryHelper operations.
 */

SUNDIALS_EXPORT int SUNMemoryHelper_CopyAsync(SUNMemoryHelper, SUNMemory dst,
                                              SUNMemory src, size_t mem_size,
                                              void* ctx);

/* Clones the SUNMemoryHelper */
SUNDIALS_EXPORT SUNMemoryHelper SUNMemoryHelper_Clone(SUNMemoryHelper);

/* Frees the SUNMemoryHelper */
SUNDIALS_EXPORT int SUNMemoryHelper_Destroy(SUNMemoryHelper);

/*
 * Utility SUNMemoryHelper functions.
 */

/* Creates an empty SUNMemoryHelper object */
SUNDIALS_EXPORT SUNMemoryHelper SUNMemoryHelper_NewEmpty();

/* Copyies the SUNMemoryHelper ops structure from src->ops to dst->ops. */
SUNDIALS_EXPORT int SUNMemoryHelper_CopyOps(SUNMemoryHelper src,
                                            SUNMemoryHelper dst);

/* Checks that all required SUNMemoryHelper ops are provided */
SUNDIALS_EXPORT
booleantype SUNMemoryHelper_ImplementsRequiredOps(SUNMemoryHelper);


#ifdef __cplusplus
}
#endif

#endif