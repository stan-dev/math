/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
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
 *---------------------------------------------------------------
 * This module contains the routines necessary to interface with
 * the ARKBANDPRE module and user-supplied Fortran routines. The
 * routines here call the generically named routines and provide
 * a standard interface to the C code of the ARKBANDPRE package.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "farkbp.h"
#include <arkode/arkode_bandpre.h>

/*=============================================================*/

/* Fortran interface to C routine ARKBandPrecInit; see farkbp.h
   for additional information */
void FARK_BPINIT(long int *N, long int *mu,
                 long int *ml, int *ier)
{
  *ier = ARKBandPrecInit(ARK_arkodemem, *N, *mu, *ml);
  return;
}

/*=============================================================*/

/* Fortran interface to C routines ARKBandPrecGetWorkSpace and
   ARKBandPrecGetNumRhsEvals; see farkbp.h for additional
   information */
void FARK_BPOPT(long int *lenrwbp, long int *leniwbp, long int *nfebp)
{
  ARKBandPrecGetWorkSpace(ARK_arkodemem, lenrwbp, leniwbp);
  ARKBandPrecGetNumRhsEvals(ARK_arkodemem, nfebp);
  return;
}

/*===============================================================
   EOF
===============================================================*/
