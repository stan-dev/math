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
 * SUNDIALS memory helper implementation that uses the standard
 * system memory allocators.
 * ----------------------------------------------------------------*/

#include <string.h>
#include <stdlib.h>

#include <sunmemory/sunmemory_system.h>
#include "sundials_debug.h"

SUNMemoryHelper SUNMemoryHelper_Sys(SUNContext sunctx)
{
  SUNMemoryHelper helper;

  /* Allocate the helper */
  helper = SUNMemoryHelper_NewEmpty(sunctx);

  /* Set the ops */
  helper->ops->alloc     = SUNMemoryHelper_Alloc_Sys;
  helper->ops->dealloc   = SUNMemoryHelper_Dealloc_Sys;
  helper->ops->copy      = SUNMemoryHelper_Copy_Sys;

  /* Attach content and ops */
  helper->content = NULL;

  return helper;
}

int SUNMemoryHelper_Alloc_Sys(SUNMemoryHelper helper, SUNMemory* memptr,
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
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sys: malloc returned NULL\n");
      free(mem);
      return(-1);
    }
  }
  else
  {
    SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Alloc_Sys: unsupported memory type\n");
    free(mem);
    return(-1);
  }

  *memptr = mem;
  return(0);
}

int SUNMemoryHelper_Dealloc_Sys(SUNMemoryHelper helper, SUNMemory mem,
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
    else
    {
      SUNDIALS_DEBUG_PRINT("ERROR in SUNMemoryHelper_Dealloc_Sys: unsupported memory type\n");
      return(-1);
    }
  }

  free(mem);
  return(0);
}

int SUNMemoryHelper_Copy_Sys(SUNMemoryHelper helper, SUNMemory dst,
                             SUNMemory src, size_t memory_size, void* queue)
{
  int retval = 0;
  memcpy(dst->ptr, src->ptr, memory_size);
  return(retval);
}
