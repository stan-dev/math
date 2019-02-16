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
 * SUNLinearSolver objects, to ensure that F2C_ARKODE_linsol and 
 * F2C_ARKODE_mass_sol are defined for cases when no linear 
 * solver object is linked in with the main executable.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"

/*=============================================================*/

/* Define global matrix variables */

SUNLinearSolver F2C_ARKODE_linsol;
SUNLinearSolver F2C_ARKODE_mass_sol; 

/*=============================================================*/

/* C routine that is called when solving an explicit problem */
void FARKNullLinsol()
{
  F2C_ARKODE_linsol = NULL;
  F2C_ARKODE_mass_sol = NULL;
}

/*===============================================================
   EOF
===============================================================*/
