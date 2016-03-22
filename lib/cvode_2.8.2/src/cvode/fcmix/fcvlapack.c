/*
 * -----------------------------------------------------------------
 * $Revision: 4075 $
 * $Date: 2014-04-24 10:46:58 -0700 (Thu, 24 Apr 2014) $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVLAPACK
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_lapack.h>

/***************************************************************************/

void FCV_LAPACKDENSE(int *neq, int *ier)
{
  /* neq  is the problem size */

  *ier = CVLapackDense(CV_cvodemem, *neq);

  CV_ls = CV_LS_LAPACKDENSE;
}

/***************************************************************************/

void FCV_LAPACKBAND(int *neq, int *mupper, int *mlower, int *ier)
{
  /* 
     neq        is the problem size
     mupper     is the upper bandwidth
     mlower     is the lower bandwidth 
  */

  *ier = CVLapackBand(CV_cvodemem, *neq, *mupper, *mlower);

  CV_ls = CV_LS_LAPACKBAND;
}

/***************************************************************************/

