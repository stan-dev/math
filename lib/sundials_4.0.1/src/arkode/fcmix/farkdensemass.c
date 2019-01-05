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
 * Fortran/C interface routines for ARKODE/ARKLS, for the case
 * of a user-supplied mass-matrix approximation routine.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_arkstep.h>
#include <sunmatrix/sunmatrix_dense.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_DMASS(long int *N, realtype *T,
                         realtype *DMASS, long int *IPAR,
                         realtype *RPAR, realtype *V1,
                         realtype *V2, realtype *V3, int *ier);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface routine to ARKStepSetMassFn; see
   farkode.h for further details */
void FARK_DENSESETMASS(int *ier)
{
  *ier = ARKStepSetMassFn(ARK_arkodemem, FARKDenseMass);
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKDMASS; see
   farkode.h for additional information  */
int FARKDenseMass(realtype t, SUNMatrix M, void *user_data,
                  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *massdata, *v1data, *v2data, *v3data;
  long int N;
  FARKUserData ARK_userdata;

  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  N = SUNDenseMatrix_Columns(M);
  massdata = SUNDenseMatrix_Column(M,0);
  ARK_userdata = (FARKUserData) user_data;

  FARK_DMASS(&N, &t, massdata, ARK_userdata->ipar, ARK_userdata->rpar,
             v1data, v2data, v3data, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
