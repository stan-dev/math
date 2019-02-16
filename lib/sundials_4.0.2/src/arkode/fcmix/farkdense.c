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
 * of a user-supplied Jacobian approximation routine.
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

  extern void FARK_DJAC(long int *N, realtype *T, realtype *Y, 
                        realtype *FY, realtype *DJAC,
                        realtype *H, long int *IPAR,
                        realtype *RPAR, realtype *V1,
                        realtype *V2, realtype *V3, int *ier);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetJacFn; see
   farkode.h for additional information */
void FARK_DENSESETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKStepSetJacFn(ARK_arkodemem, NULL);
  } else {
    *ier = ARKStepSetJacFn(ARK_arkodemem, FARKDenseJac);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKDJAC; see
   farkode.h for additional information  */
int FARKDenseJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                 void *user_data, N_Vector vtemp1, N_Vector vtemp2,
                 N_Vector vtemp3)
{
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  long int N;
  FARKUserData ARK_userdata;
  int ier = 0;

  /* Initialize all pointers to NULL */
  ydata = fydata = jacdata = v1data = v2data = v3data = NULL;

  ARKStepGetLastStep(ARK_arkodemem, &h);
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  N       = SUNDenseMatrix_Columns(J);
  jacdata = SUNDenseMatrix_Column(J,0);
  ARK_userdata = (FARKUserData) user_data;

  FARK_DJAC(&N, &t, ydata, fydata, jacdata, &h,
            ARK_userdata->ipar, ARK_userdata->rpar, v1data,
            v2data, v3data, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
