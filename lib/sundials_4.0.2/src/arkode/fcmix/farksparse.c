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
 * of a user-supplied sparse Jacobian routine.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "farkode.h"
#include "arkode_impl.h"
#include <arkode/arkode_arkstep.h>
#include <sunmatrix/sunmatrix_sparse.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

 
  extern void FARK_SPJAC(realtype *T, realtype *Y, 
                         realtype *FY, long int *N, 
                         long int *NNZ, realtype *JDATA, 
                         sunindextype *JRVALS, sunindextype *JCPTRS, 
                         realtype *H,  long int *IPAR,
                         realtype *RPAR, realtype *V1,
                         realtype *V2, realtype *V3,
                         int *ier);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKStepSetJacFn; see 
   farkode.h for further information */
void FARK_SPARSESETJAC(int *ier)
{
#if defined(SUNDIALS_INT32_T)
  arkProcessError((ARKodeMem) ARK_arkodemem, ARK_ILL_INPUT, "ARKODE",
                  "FARKSPARSESETJAC", 
                  "Sparse Fortran users must configure SUNDIALS with 64-bit integers.");
  *ier = 1;
#else  
  *ier = ARKStepSetJacFn(ARK_arkodemem, FARKSparseJac);
#endif
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKSPJAC; see 
   farkode.h for additional information  */
int FARKSparseJac(realtype t, N_Vector y, N_Vector fy, 
                  SUNMatrix J, void *user_data, N_Vector vtemp1, 
                  N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *v1data, *v2data, *v3data, *Jdata;
  realtype h;
  long int NP, NNZ; 
  sunindextype *indexvals, *indexptrs; 
  FARKUserData ARK_userdata;

  ARKStepGetLastStep(ARK_arkodemem, &h);
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  ARK_userdata = (FARKUserData) user_data;
  NP = SUNSparseMatrix_NP(J);
  NNZ = SUNSparseMatrix_NNZ(J);
  Jdata = SUNSparseMatrix_Data(J);
  indexvals = SUNSparseMatrix_IndexValues(J);
  indexptrs = SUNSparseMatrix_IndexPointers(J);

  FARK_SPJAC(&t, ydata, fydata, &NP, &NNZ, Jdata, indexvals, 
             indexptrs, &h, ARK_userdata->ipar, ARK_userdata->rpar, 
             v1data, v2data, v3data, &ier); 
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
