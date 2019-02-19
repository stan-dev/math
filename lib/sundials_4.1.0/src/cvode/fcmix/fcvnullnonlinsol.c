/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * File that provides a globally-defined, but NULL-valued, SUNNonlinearSolver
 * object, to ensure that F2C_CVODE_nonlinsol is defined for cases when the
 * default nonlinear solver is used and thus no Fortran nonlinear solver object
 * is linked in with the main executable.
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "fcvode.h"
#include "cvode_impl.h"

/*=============================================================*/

/* Define global linear solver variable */

SUNNonlinearSolver F2C_CVODE_nonlinsol;

/*=============================================================*/

/* C routine that is called when using the default nonlinear solver */
void FCVNullNonlinSol()
{
  F2C_CVODE_nonlinsol = NULL;
}

/*===============================================================
   EOF
===============================================================*/
