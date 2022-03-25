/* -----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------*/

#ifndef _SUNDIALS_CHECK_RETVAL_H_
#define _SUNDIALS_CHECK_RETVAL_H_

#include <stdio.h>

/* --------------------------------------------------------------
 * Function to check return values:
 *
 * opt == 0  means the function allocates memory and returns a
 *           pointer so check if a NULL pointer was returned
 * opt == 1  means the function returns an integer where a
 *           value < 0 indicates an error occured
 * --------------------------------------------------------------*/
static int check_retval(void *returnvalue, const char *funcname, int opt, int myid)
{
  int* errvalue;

  if (opt == 0 && returnvalue == NULL)
  {
    /* A NULL pointer was returned - no memory allocated */
    if (myid == 0)
      fprintf(stderr, "\nERROR: %s() failed - returned NULL pointer\n\n",
              funcname);
    return(1);
  }
  else if (opt == 1)
  {
    errvalue = (int *) returnvalue;

    /* A value < 0 was returned - function failed */
    if (*errvalue < 0)
    {
      if (myid == 0)
        fprintf(stderr, "\nERROR: %s() returned %d\n\n", funcname, *errvalue);
      return(1);
    }
  }

  /* return success */
  return(0);
}

#endif