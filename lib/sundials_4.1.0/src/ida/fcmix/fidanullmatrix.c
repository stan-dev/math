/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * File that provides a globally-defined, but NULL-valued, 
 * SUNMatrix object, to ensure that F2C_IDA_matrix is defined 
 * for cases when no matrix object is linked in with the main 
 * executable.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "fida.h"
#include "ida_impl.h"

/*=============================================================*/

/* Define global matrix variable */

SUNMatrix F2C_IDA_matrix;

/*=============================================================*/

/* C routine that is called when using matrix-free linear solvers */
void FIDANullMatrix()
{
  F2C_IDA_matrix = NULL;
}

/*===============================================================
   EOF
===============================================================*/
