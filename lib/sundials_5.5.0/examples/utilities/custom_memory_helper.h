/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example of a custom SUNMemoryHelper that only supports CUDA
 * unmanaged memory only and synchronous copies.
 * -----------------------------------------------------------------*/

#include <assert.h>
#include <cuda_runtime.h>
#include <sundials/sundials_memory.h>


#define MY_CUDACHK(ans) { cudaVerify((ans), __FILE__, __LINE__, 1); }
static void cudaVerify(cudaError_t code, const char *file, int line, int abort)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr, "CUDA ERROR: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) assert(false);
   }
}

int MyMemoryHelper_Alloc(SUNMemoryHelper helper, SUNMemory* memptr,
                         size_t memsize, SUNMemoryType mem_type)
{
  SUNMemory mem = SUNMemoryNewEmpty();

  mem->ptr = NULL;
  mem->own = SUNTRUE;

  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr  = malloc(memsize);
    if (mem->ptr == NULL) return(-1);
    mem->type = SUNMEMTYPE_HOST;
  }
  else if (mem_type == SUNMEMTYPE_UVM ||
           mem_type == SUNMEMTYPE_DEVICE)
  {
    MY_CUDACHK( cudaMalloc(&(mem->ptr), memsize) );
    mem->type = SUNMEMTYPE_DEVICE;
  }
  else
  {
    free(mem);
    return(-1);
  }

  *memptr = mem;
  return(0);
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
        MY_CUDACHK( cudaFree(mem->ptr) );
        mem->ptr = NULL;
      }
      else
      {
        return(-1);
      }
    }

    free(mem);
  }

  return(0);
}

int MyMemoryHelper_Copy(SUNMemoryHelper helper, SUNMemory dst,
                        SUNMemory src, size_t memory_size)
{
  switch(src->type)
  {
    case SUNMEMTYPE_HOST:
      if (dst->type == SUNMEMTYPE_HOST)
      {
        memcpy(dst->ptr, src->ptr, memory_size);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE)
      {
        MY_CUDACHK( cudaMemcpy(dst->ptr, src->ptr,
                               memory_size,
                               cudaMemcpyHostToDevice) );
      }
      break;
    case SUNMEMTYPE_DEVICE:
      if (dst->type == SUNMEMTYPE_HOST)
      {
        MY_CUDACHK( cudaMemcpy(dst->ptr, src->ptr,
                               memory_size,
                               cudaMemcpyDeviceToHost) );
      }
      else if (dst->type == SUNMEMTYPE_DEVICE)
      {
        MY_CUDACHK( cudaMemcpy(dst->ptr, src->ptr,
                               memory_size,
                               cudaMemcpyDeviceToDevice) );
      }
      break;
    default:
      return(-1);
  }

  return(0);
}

SUNMemoryHelper MyMemoryHelper()
{
  SUNMemoryHelper helper;

  /* Allocate helper */
  helper = SUNMemoryHelper_NewEmpty();

  /* Set the ops */
  helper->ops->alloc     = MyMemoryHelper_Alloc;
  helper->ops->dealloc   = MyMemoryHelper_Dealloc;
  helper->ops->copy      = MyMemoryHelper_Copy;
  helper->ops->copyasync = NULL;

  /* Attach user data and ops */
  helper->content = NULL;

  return helper;
}
