/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2018, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
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
