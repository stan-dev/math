/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * File that provides a globally-defined, but NULL-valued, 
 * SUNMatrix object, to ensure that F2C_KINSOL_matrix is defined 
 * for cases when no matrix object is linked in with the main 
 * executable.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "fkinsol.h"
#include "kinsol_impl.h"

/*=============================================================*/

/* Define global matrix variable */

SUNMatrix F2C_KINSOL_matrix;

/*=============================================================*/

/* C routine that is called when using matrix-free linear solvers 
   or fixed-point nonlinear solver */
void FKINNullMatrix()
{
  F2C_KINSOL_matrix = NULL;
}

/*===============================================================
   EOF
===============================================================*/
