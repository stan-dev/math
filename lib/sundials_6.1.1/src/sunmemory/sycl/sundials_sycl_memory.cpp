/* ----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
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

SUNMemoryHelper SUNMemoryHelper_Sycl(SUNContext sunctx)
{
  // Allocate the helper
  SUNMemoryHelper helper = SUNMemoryHelper_NewEmpty(sunctx);
  if (!helper)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: SUNMemoryHelper_NewEmpty returned NULL\n");
    return NULL;
  }

  // Set the ops
  helper->ops->alloc     = SUNMemoryHelper_Alloc_Sycl;
  helper->ops->dealloc   = SUNMemoryHelper_Dealloc_Sycl;
  helper->ops->copy      = SUNMemoryHelper_Copy_Sycl;
  helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Sycl;

  return helper;
}

int SUNMemoryHelper_Alloc_Sycl(SUNMemoryHelper helper, SUNMemory* memptr,
                               size_t mem_size, SUNMemoryType mem_type,
                               void* queue)
{
  // Check inputs
  if (!queue) {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  // Allocate the memory struct
  SUNMemory mem = SUNMemoryNewEmpty();
  if (!mem)
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Sycl: SUNMemoryNewEmpty returned NULL\n");
    return -1;
  }

  // Initialize the memory content
  mem->ptr  = nullptr;
  mem->own  = SUNTRUE;
  mem->type = mem_type;

  // Allocate the data pointer
  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr = malloc(mem_size);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc returned NULL\n");
      free(mem);
      return -1;
    }
  }
  else if (mem_type == SUNMEMTYPE_PINNED)
  {
    mem->ptr = ::sycl::malloc_host(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_host returned NULL\n");
      free(mem);
      return -1;
    }
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    mem->ptr = ::sycl::malloc_device(mem_size, *sycl_queue);
    if (!(mem->ptr))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sycl: malloc_device returned NULL\n");
      free(mem);
      return -1;
    }
  }
  else if (mem_type == SUNMEMTYPE_UVM)
  {
    mem->ptr = ::sycl::malloc_shared(mem_size, *sycl_queue);
    if (!(mem->ptr))
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

int SUNMemoryHelper_Dealloc_Sycl(SUNMemoryHelper helper, SUNMemory mem,
                                 void* queue)
{
  if (!mem) return 0;

  if (mem->ptr && mem->own)
  {
    if (!queue) {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Sycl: queue is NULL\n");
      return -1;
    }
    ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

    if (mem->type == SUNMEMTYPE_HOST)
    {
      free(mem->ptr);
      mem->ptr = nullptr;
    }
    else if (mem->type == SUNMEMTYPE_PINNED ||
             mem->type == SUNMEMTYPE_DEVICE ||
             mem->type == SUNMEMTYPE_UVM)
    {
      ::sycl::free(mem->ptr, *sycl_queue);
      mem->ptr = nullptr;
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
                              SUNMemory src, size_t memory_size, void* queue)
{
  if (!queue) {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Copy_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  if (SUNMemoryHelper_CopyAsync_Sycl(helper, dst, src, memory_size, queue))
    return -1;
  sycl_queue->wait_and_throw();
  return 0;
}


int SUNMemoryHelper_CopyAsync_Sycl(SUNMemoryHelper helper, SUNMemory dst,
                                   SUNMemory src, size_t memory_size,
                                   void* queue)
{
  if (!queue) {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_CopyAsync_Sycl: queue is NULL\n");
    return -1;
  }
  ::sycl::queue* sycl_queue = static_cast<::sycl::queue*>(queue);

  if (src->type == SUNMEMTYPE_HOST && dst->type == SUNMEMTYPE_HOST)
  {
    memcpy(dst->ptr, src->ptr, memory_size);
  }
  else
  {
    sycl_queue->memcpy(dst->ptr, src->ptr, memory_size);
  }
  return 0;
}
