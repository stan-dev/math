/*-----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Aaron Collier @ LLNL
 *-----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *-----------------------------------------------------------------
 * Fortran/C interface routines for IDA/IDALS, for the case of
 * a user-supplied Jacobian approximation routine.
 *-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fida.h"     /* function names, prototypes, global vars.*/
#include "ida_impl.h" /* definition of IDAMem type               */

#include <ida/ida_ls.h>
#include <sunmatrix/sunmatrix_band.h>

/*************************************************/

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

  extern void FIDA_BJAC(long int* N, long int* MU, long int* ML,
                        long int* EBAND, realtype* T, realtype* Y,
                        realtype* YP, realtype* R, realtype* CJ,
                        realtype* J, realtype* EWT, realtype* H,
                        long int* IPAR, realtype* RPAR, realtype* V1, 
                        realtype* V2, realtype* V3, int* IER);
#ifdef __cplusplus
}
#endif

/*************************************************/

void FIDA_BANDSETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = IDASetJacFn(IDA_idamem, NULL);
  } else {
    if (F2C_IDA_ewtvec == NULL) {
      F2C_IDA_ewtvec = N_VClone(F2C_IDA_vec);
      if (F2C_IDA_ewtvec == NULL) {
        *ier = -1;
        return;
      }
    }
    *ier = IDASetJacFn(IDA_idamem, FIDABandJac);
  }
  return;
}

/*************************************************/

int FIDABandJac(realtype t, realtype c_j, N_Vector yy,  
		N_Vector yp, N_Vector rr, SUNMatrix J,
		void *user_data, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3)
{
  realtype *yy_data, *yp_data, *rr_data, *jacdata, *ewtdata, *v1data, *v2data, *v3data;
  realtype h;
  long int N, mupper, mlower, smu, eband;
  int ier;
  FIDAUserData IDA_userdata;

  /* Initialize all pointers to NULL */
  yy_data = yp_data = rr_data = jacdata = ewtdata = NULL;
  v1data = v2data = v3data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  IDAGetErrWeights(IDA_idamem, F2C_IDA_ewtvec);
  IDAGetLastStep(IDA_idamem, &h);

  /* Get pointers to vector data */
  yy_data = N_VGetArrayPointer(yy);
  yp_data = N_VGetArrayPointer(yp);
  rr_data = N_VGetArrayPointer(rr);
  ewtdata = N_VGetArrayPointer(F2C_IDA_ewtvec);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);

  N       = SUNBandMatrix_Columns(J);
  mupper  = SUNBandMatrix_UpperBandwidth(J);
  mlower  = SUNBandMatrix_LowerBandwidth(J);
  smu     = SUNBandMatrix_StoredUpperBandwidth(J);
  eband   = smu + mlower + 1;
  jacdata = SUNBandMatrix_Column(J,0) - mupper;

  IDA_userdata = (FIDAUserData) user_data;

  /* Call user-supplied routine */
  FIDA_BJAC(&N, &mupper, &mlower, &eband, &t, yy_data, yp_data,
            rr_data, &c_j, jacdata, ewtdata, &h, 
            IDA_userdata->ipar, IDA_userdata->rpar,
            v1data, v2data, v3data, &ier);

  return(ier);
}
