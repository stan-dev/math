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
 * File that provides globally-defined, but NULL-valued, 
 * SUNMatrix objects, to ensure that F2C_ARKODE_matrix and 
 * F2C_ARKODE_mass_matrix are defined for cases when no matrix 
 * object is linked in with the main executable.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"

/*=============================================================*/

/* Define global matrix variables */

SUNMatrix F2C_ARKODE_matrix;
SUNMatrix F2C_ARKODE_mass_matrix;

/*=============================================================*/

/* C routine that is called when solving an explicit problem, or 
   when using matrix-free linear solvers */
void FARKNullMatrix()
{
  F2C_ARKODE_matrix = NULL;
  F2C_ARKODE_mass_matrix = NULL;
}

/*===============================================================
   EOF
===============================================================*/
