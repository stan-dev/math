/* -----------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
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
 * Example of a custom SUNMemoryHelper that only supports SYCL
 * unmanaged memory only and synchronous copies.
 * -----------------------------------------------------------------*/

#include <cstdlib>
#include <CL/sycl.hpp>
#include <sundials/sundials_memory.h>

#define MY_QUEUE(h) (*((sycl::queue*)(h->content)))

int MyMemoryHelper_Alloc(SUNMemoryHelper helper, SUNMemory* memptr,
                         size_t memsize, SUNMemoryType mem_type)
{
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
    mem->ptr = sycl::malloc_device(memsize, MY_QUEUE(helper));
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

int MyMemoryHelper_Dealloc(SUNMemoryHelper helper, SUNMemory mem)
{
  if (mem != NULL)
  {
    if (mem->ptr != NULL && mem->own)
    {
      if (mem->type == SUNMEMTYPE_HOST)
      {
        free(mem->ptr);
        mem->ptr = NULL;
      }
      else if (mem->type == SUNMEMTYPE_DEVICE)
      {
        sycl::free(mem->ptr, MY_QUEUE(helper));
        mem->ptr = NULL;
      }
      else
      {
        return -1;
      }
    }

    free(mem);
  }

  return 0;
}

int MyMemoryHelper_Copy(SUNMemoryHelper helper, SUNMemory dst,
                        SUNMemory src, size_t memory_size)
{
  if (src->type == SUNMEMTYPE_HOST && dst->type == SUNMEMTYPE_HOST)
  {
    memcpy(dst->ptr, src->ptr, memory_size);
  }
  else
  {
    MY_QUEUE(helper).memcpy(dst->ptr, src->ptr, memory_size);
  }
  MY_QUEUE(helper).wait_and_throw();
  return 0;
}


SUNMemoryHelper MyMemoryHelper_Clone(SUNMemoryHelper helper)
{
  if (helper == NULL) return NULL;

  SUNMemoryHelper new_helper = SUNMemoryHelper_NewEmpty();
  if (helper == NULL) return NULL;

  if (SUNMemoryHelper_CopyOps(helper, new_helper)) return NULL;

  new_helper->content = helper->content;

  return new_helper;
}

int MyMemoryHelper_Destroy(SUNMemoryHelper helper)
{
  helper->content = NULL;
  free(helper->ops);
  free(helper);
  return 0;
}

SUNMemoryHelper MyMemoryHelper(sycl::queue *Q)
{
  SUNMemoryHelper helper = SUNMemoryHelper_NewEmpty();
  if (helper == NULL) return NULL;

  /* Set the ops */
  helper->ops->alloc   = MyMemoryHelper_Alloc;
  helper->ops->dealloc = MyMemoryHelper_Dealloc;
  helper->ops->copy    = MyMemoryHelper_Copy;
  helper->ops->clone   = MyMemoryHelper_Clone;
  helper->ops->destroy = MyMemoryHelper_Destroy;

  /* Attach the sycl queue pointer as the content */
  helper->content = (void*) Q;

  return helper;
}
