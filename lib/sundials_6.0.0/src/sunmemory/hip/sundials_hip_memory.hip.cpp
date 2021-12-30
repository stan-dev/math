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
 * SUNDIALS HIP memory helper implementation.
 * ----------------------------------------------------------------*/

#include <cstdlib>

#include <sunmemory/sunmemory_hip.h>
#include "sundials_debug.h"
#include "sundials_hip.h"


SUNMemoryHelper SUNMemoryHelper_Hip(SUNContext sunctx)
{
  SUNMemoryHelper helper;

  /* Allocate the helper */
  helper = SUNMemoryHelper_NewEmpty(sunctx);

  /* Set the ops */
  helper->ops->alloc     = SUNMemoryHelper_Alloc_Hip;
  helper->ops->dealloc   = SUNMemoryHelper_Dealloc_Hip;
  helper->ops->copy      = SUNMemoryHelper_Copy_Hip;
  helper->ops->copyasync = SUNMemoryHelper_CopyAsync_Hip;

  /* Attach content and ops */
  helper->content = NULL;

  return helper;
}

int SUNMemoryHelper_Alloc_Hip(SUNMemoryHelper helper, SUNMemory* memptr,
                              size_t mem_size, SUNMemoryType mem_type,
                              void* queue)
{
  SUNMemory mem = SUNMemoryNewEmpty();

  mem->ptr  = NULL;
  mem->own  = SUNTRUE;
  mem->type = mem_type;

  if (mem_type == SUNMEMTYPE_HOST)
  {
    mem->ptr = malloc(mem_size);
    if (mem->ptr == NULL)
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: malloc returned NULL\n");
      free(mem);
      return(-1);
    }
  }
  else if (mem_type == SUNMEMTYPE_PINNED)
  {
    if (!SUNDIALS_HIP_VERIFY(hipMallocHost(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: hipMallocHost failed\n");
      free(mem);
      return(-1);
    }
  }
  else if (mem_type == SUNMEMTYPE_DEVICE)
  {
    if (!SUNDIALS_HIP_VERIFY(hipMalloc(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: hipMalloc failed\n");
      free(mem);
      return(-1);
    }
  }
  else if (mem_type == SUNMEMTYPE_UVM)
  {
    if (!SUNDIALS_HIP_VERIFY(hipMallocManaged(&(mem->ptr), mem_size)))
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: hipMallocManaged failed\n");
      free(mem);
      return(-1);
    }
  }
  else
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Hip: unknown memory type\n");
    free(mem);
    return(-1);
  }

  *memptr = mem;
  return(0);
}

int SUNMemoryHelper_Dealloc_Hip(SUNMemoryHelper helper, SUNMemory mem,
                                void* queue)
{
  if (mem == NULL) return(0);

  if (mem->ptr != NULL && mem->own)
  {
    if (mem->type == SUNMEMTYPE_HOST)
    {
      free(mem->ptr);
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_PINNED)
    {
      if (!SUNDIALS_HIP_VERIFY(hipFreeHost(mem->ptr)))
      {
        SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Hip: hipFreeHost failed\n");
        return(-1);
      }
      mem->ptr = NULL;
    }
    else if (mem->type == SUNMEMTYPE_DEVICE ||
             mem->type == SUNMEMTYPE_UVM)
    {
      if (!SUNDIALS_HIP_VERIFY(hipFree(mem->ptr)))
      {
        SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Hip: hipFree failed\n");
        return(-1);
      }
      mem->ptr = NULL;
    }
    else
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Hip: unknown memory type\n");
      return(-1);
    }
  }

  free(mem);
  return(0);
}

int SUNMemoryHelper_Copy_Hip(SUNMemoryHelper helper, SUNMemory dst,
                             SUNMemory src, size_t memory_size, void* queue)
{
  int retval = 0;
  hipError_t cuerr = hipSuccess;

  switch(src->type)
  {
    case SUNMEMTYPE_HOST:
    case SUNMEMTYPE_PINNED:
      if (dst->type == SUNMEMTYPE_HOST ||
          dst->type == SUNMEMTYPE_PINNED)
      {
        memcpy(dst->ptr, src->ptr, memory_size);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE ||
               dst->type == SUNMEMTYPE_UVM)
      {
        cuerr = hipMemcpy(dst->ptr, src->ptr,
                           memory_size,
                           hipMemcpyHostToDevice);
      }
      if (!SUNDIALS_HIP_VERIFY(cuerr)) retval = -1;
      break;
    case SUNMEMTYPE_UVM:
    case SUNMEMTYPE_DEVICE:
      if (dst->type == SUNMEMTYPE_HOST ||
          dst->type == SUNMEMTYPE_PINNED)
      {
        cuerr = hipMemcpy(dst->ptr, src->ptr,
                           memory_size,
                           hipMemcpyDeviceToHost);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE ||
               dst->type == SUNMEMTYPE_UVM)
      {
        cuerr = hipMemcpy(dst->ptr, src->ptr,
                           memory_size,
                           hipMemcpyDeviceToDevice);
      }
      if (!SUNDIALS_HIP_VERIFY(cuerr)) retval = -1;
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_CopyAsync_Hip: unknown memory type\n");
      retval = -1;
  }

  return(retval);
}

int SUNMemoryHelper_CopyAsync_Hip(SUNMemoryHelper helper, SUNMemory dst,
                                  SUNMemory src, size_t memory_size,
                                  void* queue)
{
  int retval = 0;
  hipError_t cuerr = hipSuccess;
  hipStream_t stream = 0;

  if (queue != NULL)
  {
    stream = *((hipStream_t*) queue);
  }

  switch(src->type)
  {
    case SUNMEMTYPE_HOST:
    case SUNMEMTYPE_PINNED:
      if (dst->type == SUNMEMTYPE_HOST ||
          dst->type == SUNMEMTYPE_PINNED)
      {
        memcpy(dst->ptr, src->ptr, memory_size);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE ||
               dst->type == SUNMEMTYPE_UVM)
      {
        cuerr = hipMemcpyAsync(dst->ptr, src->ptr,
                                memory_size,
                                hipMemcpyHostToDevice,
                                stream);
      }
      if (!SUNDIALS_HIP_VERIFY(cuerr)) retval = -1;
      break;
    case SUNMEMTYPE_UVM:
    case SUNMEMTYPE_DEVICE:
      if (dst->type == SUNMEMTYPE_HOST ||
          dst->type == SUNMEMTYPE_PINNED)
      {
        cuerr = hipMemcpyAsync(dst->ptr, src->ptr,
                                memory_size,
                                hipMemcpyDeviceToHost,
                                stream);
      }
      else if (dst->type == SUNMEMTYPE_DEVICE ||
              dst->type == SUNMEMTYPE_UVM)
      {
        cuerr = hipMemcpyAsync(dst->ptr, src->ptr,
                                memory_size,
                                hipMemcpyDeviceToDevice,
                                stream);
      }
      if (!SUNDIALS_HIP_VERIFY(cuerr)) retval = -1;
      break;
    default:
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_CopyAsync_Hip: unknown memory type\n");
      retval = -1;
  }

  return(retval);
}
