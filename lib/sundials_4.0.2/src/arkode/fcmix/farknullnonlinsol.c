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
 * SUNNonlinearSolver object, to ensure that F2C_ARKODE_nonlinsol  
 * isdefined for cases when no Fortran-defined nonlinear solver 
 * object is linked in with the main executable.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"

/*=============================================================*/

/* Define global variable */

SUNNonlinearSolver F2C_ARKODE_nonlinsol;

/*=============================================================*/

/* C routine that is called when solving an explicit problem */
void FARKNullNonlinsol()
{
  F2C_ARKODE_nonlinsol = NULL;
}

/*===============================================================
   EOF
===============================================================*/
