/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Fortran/C interface routines for CVODE/CVLS, for the case
 * of a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "fcvode.h"     /* actual fn. names, prototypes and global vars.*/
#include "cvode_impl.h" /* definition of CVodeMem type                  */

#include <cvode/cvode_ls.h>
#include <sunmatrix/sunmatrix_dense.h>

/***************************************************************************/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_DJAC(long int *N, realtype *T, realtype *Y,
                       realtype *FY, realtype *DJAC, realtype *H,
                       long int *IPAR, realtype *RPAR, realtype *V1, 
                       realtype *V2, realtype *V3, int *ier);
#ifdef __cplusplus
}
#endif

/***************************************************************************/

void FCV_DENSESETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = CVodeSetJacFn(CV_cvodemem, NULL);
  } else {
    *ier = CVodeSetJacFn(CV_cvodemem, FCVDenseJac);
  }
}

/***************************************************************************/

/* C function CVDenseJac interfaces between CVODE and a Fortran subroutine
   FCVDJAC for solution of a linear system with dense Jacobian approximation.
   Addresses of arguments are passed to FCVDJAC, using accessor functions 
   from the SUNDenseMatrix and N_Vector modules. */

int FCVDenseJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector vtemp1, N_Vector vtemp2,
                N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *jacdata, *v1data, *v2data, *v3data;
  realtype h;
  long int N;
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);

  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  N = SUNDenseMatrix_Columns(J);
  jacdata = SUNDenseMatrix_Column(J,0);

  CV_userdata = (FCVUserData) user_data;

  FCV_DJAC(&N, &t, ydata, fydata, jacdata, &h, 
           CV_userdata->ipar, CV_userdata->rpar, v1data,
           v2data, v3data, &ier); 
  return(ier);
}

