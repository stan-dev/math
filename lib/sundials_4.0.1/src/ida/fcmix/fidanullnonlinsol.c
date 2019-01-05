/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------------------
 * File that provides a globally-defined, but NULL-valued, SUNNonlinearSolver
 * object, to ensure that F2C_IDA_nonlinsol is defined for cases when the
 * default nonlinear solver is used and thus no Fortran nonlinear solver object
 * is linked in with the main executable.
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "fida.h"
#include "ida_impl.h"

/*=============================================================*/

/* Define global linear solver variable */

SUNNonlinearSolver F2C_IDA_nonlinsol;

/*=============================================================*/

/* C routine that is called when using the default nonlinear solver */
void FIDANullNonlinSol()
{
  F2C_IDA_nonlinsol = NULL;
}

/*===============================================================
   EOF
===============================================================*/
