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
 * Fortran/C interface routines for ARKODE, for the case of a 
 * user-supplied step adaptivity routine.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_ADAPT(realtype *Y, realtype *T, realtype *H1, 
                         realtype *H2, realtype *H3, realtype *E1, 
                         realtype *E2, realtype *E3, int *Q, int *P, 
                         realtype *HNEW, long int *IPAR, 
                         realtype *RPAR, int *IER);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetAdaptivityFn; see 
   farkode.h for further information */
void FARK_ADAPTSET(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = ARKStepSetAdaptivityFn(ARK_arkodemem, NULL, NULL);
  } else {
    *ier = ARKStepSetAdaptivityFn(ARK_arkodemem, FARKAdapt, 
                                  ARK_arkodemem);
  }
  return;
}

/*=============================================================*/

/* C interface to user-supplied fortran routine FARKADAPT; see 
   farkode.h for further information */
int FARKAdapt(N_Vector y, realtype t, realtype h1, realtype h2, 
              realtype h3, realtype e1, realtype e2, realtype e3, 
              int q, int p, realtype *hnew, void *user_data)
{
  int ier = 0;
  realtype *ydata;
  FARKUserData ARK_userdata;

  ydata  = N_VGetArrayPointer(y);
  ARK_userdata = (FARKUserData) user_data;

  FARK_ADAPT(ydata, &t, &h1, &h2, &h3, &e1, &e2, &e3, &q, &p, hnew, 
             ARK_userdata->ipar, ARK_userdata->rpar, &ier);
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
