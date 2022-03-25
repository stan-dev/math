/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Example of a custom SUNMemoryHelper that only supports SYCL
 * unmanaged memory only and synchronous copies.
 * -----------------------------------------------------------------*/

#include <cstdlib>
#include <CL/sycl.hpp>
#include <sundials/sundials_memory.h>

int MyMemoryHelper_Alloc(SUNMemoryHelper helper, SUNMemory* memptr,
                         size_t memsize, SUNMemoryType mem_type, void* queue)
{
  if (!queue) return -1;
  sycl::queue* sycl_queue = static_cast<sycl::queue*>(queue);

  SUNMemory mem = SUNMemoryNewEmpty();
  if (mem == NULL) return -1;

  mem->ptr = NULL;
  mem->own = SUNTRUE;

  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr  = malloc(memsize);
    if (mem->ptr == NULL) { free(mem); return -1; }
    mem->type = SUNMEMTYPE_HOST;
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    mem->ptr = sycl::malloc_device(memsize, *sycl_queue);
    if (mem->ptr == NULL) { free(mem); return -1; }
    mem->type = SUNMEMTYPE_DEVICE;
  }
  else
  {
    free(mem);
    return -1;
  }

  *memptr = mem;
  return 0;
}

int MyMemoryHelper_Dealloc(SUNMemoryHelper helper, SUNMemory mem, void* queue)
{
  if (!mem) return 0;

  if (mem->ptr && mem->own)
  {
    if (!queue) return -1;

    sycl::queue* sycl_queue = static_cast<sycl::queue*>(queue);

    if (mem->type == SUNMEMTYPE_HOST)
    {
      free(mem->ptr);
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_DEVICE)
    {
      sycl::free(mem->ptr, *sycl_queue);
      mem->ptr = NULL;
    }
    else
    {
      return -1;
    }
  }

  free(mem);

  return 0;
}

int MyMemoryHelper_Copy(SUNMemoryHelper helper, SUNMemory dst,
                        SUNMemory src, size_t memory_size, void* queue)
{
  if (!queue) return -1;
  sycl::queue* sycl_queue = static_cast<sycl::queue*>(queue);

  if (src->type == SUNMEMTYPE_HOST && dst->type == SUNMEMTYPE_HOST)
  {
    memcpy(dst->ptr, src->ptr, memory_size);
  }
  else
  {
    sycl_queue->memcpy(dst->ptr, src->ptr, memory_size);
  }
  sycl_queue->wait_and_throw();
  return 0;
}

SUNMemoryHelper MyMemoryHelper(SUNContext sunctx)
{
  SUNMemoryHelper helper = SUNMemoryHelper_NewEmpty(sunctx);
  if (helper == NULL) return NULL;

  /* Set the ops */
  helper->ops->alloc   = MyMemoryHelper_Alloc;
  helper->ops->dealloc = MyMemoryHelper_Dealloc;
  helper->ops->copy    = MyMemoryHelper_Copy;

  return helper;
}
