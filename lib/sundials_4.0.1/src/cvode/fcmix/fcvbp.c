/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *      Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2017, Southern Methodist University and 
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
 * -----------------------------------------------------------------
 * This module contains the routines necessary to interface with the
 * CVBANDPRE module and user-supplied Fortran routines.
 * The routines here call the generically named routines and provide 
 * a standard interface to the C code of the CVBANDPRE package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"                 /* actual fn. names, prototypes and global vars.*/
#include "fcvbp.h"                  /* prototypes of interfaces to CVBANDPRE        */

#include <cvode/cvode_bandpre.h>    /* prototypes of CVBANDPRE functions and macros */

/***************************************************************************/

void FCV_BPINIT(long int *N, long int *mu, long int *ml, int *ier)
{
  /* 
     Call CVBandPrecInit to initialize the CVBANDPRE module:
     N      is the vector size
     mu, ml are the half-bandwidths of the retained preconditioner blocks
  */

  *ier = CVBandPrecInit(CV_cvodemem, *N, *mu, *ml);

  return;
}

/***************************************************************************/

/* C function FCVBPOPT to access optional outputs from CVBANDPRE_Data */

void FCV_BPOPT(long int *lenrwbp, long int *leniwbp, long int *nfebp)
{
  CVBandPrecGetWorkSpace(CV_cvodemem, lenrwbp, leniwbp);
  CVBandPrecGetNumRhsEvals(CV_cvodemem, nfebp);
}
