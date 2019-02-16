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
#include <sunmatrix/sunmatrix_sparse.h>

/*=============================================================*/

/* Prototype of the Fortran routine */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FARK_SPMASS(realtype *T, long int *N, 
                          long int *NNZ, realtype *MDATA, 
                          sunindextype *MRVALS, sunindextype *MCPTRS, 
                          long int *IPAR, realtype *RPAR, 
                          realtype *V1, realtype *V2, realtype *V3, 
                          int *ier);

#ifdef __cplusplus
}
#endif

/*=============================================================*/

/* Fortran interface to C routine ARKSlsSetMassFn; see 
   farkode.h for further information */
void FARK_SPARSESETMASS(int *ier)
{
  *ier = ARKStepSetMassFn(ARK_arkodemem, FARKSparseMass);
}

/*=============================================================*/

/* C interface to user-supplied Fortran routine FARKSPMASS; see 
   farkode.h for additional information  */
int FARKSparseMass(realtype t, SUNMatrix MassMat, void *user_data, 
                   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *v1data, *v2data, *v3data, *Mdata;
  FARKUserData ARK_userdata;
  long int NP, NNZ;
  sunindextype *indexvals, *indexptrs;
  
  v1data = N_VGetArrayPointer(vtemp1);
  v2data = N_VGetArrayPointer(vtemp2);
  v3data = N_VGetArrayPointer(vtemp3);
  NP = SUNSparseMatrix_NP(MassMat);
  NNZ = SUNSparseMatrix_NNZ(MassMat);
  Mdata = SUNSparseMatrix_Data(MassMat);
  indexvals = SUNSparseMatrix_IndexValues(MassMat);
  indexptrs = SUNSparseMatrix_IndexPointers(MassMat);
  ARK_userdata = (FARKUserData) user_data;

  FARK_SPMASS(&t, &NP, &NNZ, Mdata, indexvals, indexptrs, 
              ARK_userdata->ipar, ARK_userdata->rpar, v1data, 
              v2data, v3data, &ier); 
  return(ier);
}

/*===============================================================
   EOF
===============================================================*/
