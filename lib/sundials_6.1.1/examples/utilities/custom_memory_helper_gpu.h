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
 * Example of a custom SUNMemoryHelper that only supports CUDA/HIP
 * unmanaged memory only and synchronous copies.
 * -----------------------------------------------------------------*/

#include <assert.h>
#include <string.h>
#if defined(__NVCC__)
#include <cuda_runtime.h>
#define EX_USES_CUDA
#elif defined(__HCC__) || (defined(__clang__) && defined(__HIP__))
#include <hip/hip_runtime.h>
#define EX_USES_HIP
#endif

#include <sundials/sundials_memory.h>

#if defined(EX_USES_CUDA)
#define MY_GPU(a) cuda ## a
#elif defined(EX_USES_HIP)
#define MY_GPU(a) hip ## a
#endif

#define MY_GPUCHK(ans) { gpuVerify((ans), __FILE__, __LINE__, 1); }

static void gpuVerify(MY_GPU(Error_t) code, const char *file, int line, int abort)
{
   if (code != MY_GPU(Success))
   {
      fprintf(stderr, "GPU ERROR: %s %s %d\n", MY_GPU(GetErrorString)(code), file, line);
      if (abort) assert(false);
   }
}

int MyMemoryHelper_Alloc(SUNMemoryHelper helper, SUNMemory* memptr,
                         size_t memsize, SUNMemoryType mem_type, void* queue)
{
  SUNMemory mem = SUNMemoryNewEmpty();

  mem->ptr = NULL;
  mem->own = SUNTRUE;

  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr  = malloc(memsize);
    if (mem->ptr == NULL) { free(mem); return(-1); }
    mem->type = SUNMEMTYPE_HOST;
  }
  else if (mem_type == SUNMEMTYPE_UVM ||
           mem_type == SUNMEMTYPE_DEVICE)
  {
    MY_GPUCHK( MY_GPU(Malloc)(&(mem->ptr), memsize) );
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

int MyMemoryHelper_Dealloc(SUNMemoryHelper helper, SUNMemory mem, void* queue)
{
  if (!mem) return 0;

  if (mem->ptr && mem->own)
  {
    if (mem->type == SUNMEMTYPE_HOST)
    {
      free(mem->ptr);
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_DEVICE)
    {
      MY_GPUCHK( MY_GPU(Free)(mem->ptr) );
      mem->ptr = NULL;
    }
    else
    {
      return(-1);
    }
  }

  free(mem);

  return(0);
}

int MyMemoryHelper_Copy(SUNMemoryHelper helper, SUNMemory dst,
                        SUNMemory src, size_t memory_size, void* queue)
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

        MY_GPUCHK( MY_GPU(Memcpy)(dst->ptr, src->ptr,
                                  memory_size,
                                  MY_GPU(MemcpyHostToDevice)) );
      }
      break;
    case SUNMEMTYPE_DEVICE:
      if (dst->type == SUNMEMTYPE_HOST)
      {
        MY_GPUCHK( MY_GPU(Memcpy)(dst->ptr, src->ptr,
                                  memory_size,
                                  MY_GPU(MemcpyDeviceToHost)) );
      }
      else if (dst->type == SUNMEMTYPE_DEVICE)
      {
        MY_GPUCHK( MY_GPU(Memcpy)(dst->ptr, src->ptr,
                                  memory_size,
                                  MY_GPU(MemcpyDeviceToDevice)) );
      }
      break;
    default:
      return(-1);
  }

  return(0);
}

SUNMemoryHelper MyMemoryHelper(SUNContext sunctx)
{
  SUNMemoryHelper helper;

  /* Allocate helper */
  helper = SUNMemoryHelper_NewEmpty(sunctx);

  /* Set the ops */
  helper->ops->alloc     = MyMemoryHelper_Alloc;
  helper->ops->dealloc   = MyMemoryHelper_Dealloc;
  helper->ops->copy      = MyMemoryHelper_Copy;
  helper->ops->copyasync = NULL;

  /* Attach user data and ops */
  helper->content = NULL;

  return helper;
}
