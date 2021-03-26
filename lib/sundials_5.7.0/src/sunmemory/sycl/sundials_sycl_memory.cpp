/* ----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------
 * SUNDIALS SYCL memory helper implementation.
 * ----------------------------------------------------------------*/

#include <cstdlib>

#include <sunmemory/sunmemory_sycl.h>
#include "sundials_debug.h"

#define SYCL_QUEUE(h) (*((sycl::queue*)(h->content)))

SUNMemoryHelper SUNMemoryHelper_Sycl(sycl::queue *Q)
{
  // Check for non-NULL queue
  if (Q == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: input queue is NULL\n");
    return NULL;
  }

  // Allocate the helper
  SUNMemoryHelper helper = SUNMemoryHelper_NewEmpty();
  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: SUNMemoryHelper_NewEmpty returned NULL\n");
    return NULL;
  }

  // Set the ops
  helper->ops->alloc     = SUNMemoryHelper_Alloc_Sycl;
  helper->ops->dealloc   = SUNMemoryHelper_Dealloc_Sycl;
  helper->ops->copy      = SUNMemoryHelper_Copy_Sycl;
  helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Sycl;
  helper->ops->clone     = SUNMemoryHelper_Clone_Sycl;
  helper->ops->destroy   = SUNMemoryHelper_Destroy_Sycl;

  // Attach the sycl queue pointer as the content
  helper->content = (void*) Q;

  return helper;
}

SUNMemoryHelper SUNMemoryHelper_Clone_Sycl(SUNMemoryHelper helper)
{
  // Check input
  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Clone_Sycl: input helper is NULL\n");
    return NULL;
  }

  // Allocate the helper
  SUNMemoryHelper new_helper = SUNMemoryHelper_NewEmpty();
  if (helper == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: SUNMemoryHelper_NewEmpty returned NULL\n");
    return NULL;
  }

  // Set the ops
  if (SUNMemoryHelper_CopyOps(helper, new_helper))
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: SUNMemoryHelper_CopyOps returned nonzero\n");
    return NULL;
  }

  // Copy the sycl queue pointer
  new_helper->content = helper->content;

  return new_helper;
}

int SUNMemoryHelper_Destroy_Sycl(SUNMemoryHelper helper)
{
  helper->content = NULL;
  free(helper->ops);
  free(helper);
  return 0;
}

int SUNMemoryHelper_Alloc_Sycl(SUNMemoryHelper helper, SUNMemory* memptr,
                               size_t mem_size, SUNMemoryType mem_type)
{
  // Allocate the memory struct
  SUNMemory mem = SUNMemoryNewEmpty();
  if (mem == NULL)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: SUNMemoryNewEmpty returned NULL\n");
    return -1;
  }

  // Initialize the memory content
  mem->ptr  = NULL;
  mem->own  = SUNTRUE;
  mem->type = mem_type;

  // Allocate the data pointer
  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr = malloc(mem_size);
    if (mem->ptr == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc returned NULL\n");
      free(mem);
      return -1;
    }
  }
  else if (mem_type == SUNMEMTYPE_PINNED)
  {
    mem->ptr = sycl::malloc_host(mem_size, SYCL_QUEUE(helper));
    if (mem->ptr == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_host returned NULL\n");
      free(mem);
      return -1;
    }
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    mem->ptr = sycl::malloc_device(mem_size, SYCL_QUEUE(helper));
    if (mem->ptr == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_device returned NULL\n");
      free(mem);
      return -1;
    }
  }
  else if (mem_type == SUNMEMTYPE_UVM)
  {
    mem->ptr = sycl::malloc_shared(mem_size, SYCL_QUEUE(helper));
    if (mem->ptr == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_shared returned NULL\n");
      free(mem);
      return -1;
    }
  }
  else
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: unknown memory type\n");
    free(mem);
    return -1;
  }

  *memptr = mem;
  return 0;
}

int SUNMemoryHelper_Dealloc_Sycl(SUNMemoryHelper helper, SUNMemory mem)
{
  if (mem == NULL) return 0;

  if (mem->ptr != NULL && mem->own)
  {
    if (mem->type == SUNMEMTYPE_HOST)
    {
      free(mem->ptr);
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_PINNED ||
             mem->type == SUNMEMTYPE_DEVICE ||
             mem->type == SUNMEMTYPE_UVM)
    {
      sycl::free(mem->ptr, SYCL_QUEUE(helper));
      mem->ptr = NULL;
    }
    else
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Sycl: unknown memory type\n");
      return -1;
    }
  }

  free(mem);
  return 0;
}


int SUNMemoryHelper_Copy_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                              SUNMemory src, size_t memory_size)
{
  if (SUNMemoryHelper_CopyAsync_Sycl(helper, dst, src, memory_size, NULL))
    return -1;
  SYCL_QUEUE(helper).wait_and_throw();
  return 0;
}


int SUNMemoryHelper_CopyAsync_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                                   SUNMemory src, size_t memory_size,
                                   void* ctx)
{
  if (src->type == SUNMEMTYPE_HOST && dst->type == SUNMEMTYPE_HOST)
  {
    memcpy(dst->ptr, src->ptr, memory_size);
  }
  else
  {
    SYCL_QUEUE(helper).memcpy(dst->ptr, src->ptr, memory_size);
  }
  return 0;
}
