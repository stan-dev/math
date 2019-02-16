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
 * Fortran/C interface routines for ARKODE/ARKLS, for the case
 * of a user-supplied mass-matrix approximation routine.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_arkstep.h>
#include <sunmatrix/sunmatrix_band.h>


/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_BMASS(long int *N, long int *MU, 
                         long int *ML, long int *EBAND, 
                         realtype *T, realtype *BMASS,
                         long int *IPAR, realtype *RPAR,
                         realtype *V1, realtype *V2, realtype *V3, 
                         int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface routine to ARKStepSetMassFn; see farkode.h
   for further details */
void FARK_BANDSETMASS(int *ier)
{
  *ier = ARKStepSetMassFn(ARK_arkodemem, FARKBandMass);
}

/*=============================================================*/

/* C interface to user-supplied Fortran subroutine FARKBMASS; see
   farkode.h for further details */
int FARKBandMass(realtype t, SUNMatrix M, void *user_data,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *massdata, *v1data, *v2data, *v3data;
  long int N, mupper, mlower, smu, eband;
  FARKUserData ARK_userdata;

  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  N = SUNBandMatrix_Columns(M);
  mupper = SUNBandMatrix_UpperBandwidth(M);
  mlower = SUNBandMatrix_LowerBandwidth(M);
  smu = SUNBandMatrix_StoredUpperBandwidth(M);
  eband   = smu + mlower + 1;
  massdata = SUNBandMatrix_Column(M,0) - mupper;
  ARK_userdata = (FARKUserData) user_data;

  FARK_BMASS(&N, &mupper, &mlower, &eband, &t, massdata,
             ARK_userdata->ipar, ARK_userdata->rpar, v1data,
             v2data, v3data, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
