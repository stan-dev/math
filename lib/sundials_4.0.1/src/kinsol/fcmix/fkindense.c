/* -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
 *                David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 * Fortran/C interface routines for KINSOL/KINLS, for the case
 * of a user-supplied Jacobian approximation routine.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"     /* prototypes of standard interfaces and global vars.*/
#include "kinsol_impl.h" /* definition of KINMem type                         */

#include <kinsol/kinsol_ls.h>
#include <sunmatrix/sunmatrix_dense.h>

/*
 * ----------------------------------------------------------------
 * prototypes of the user-supplied fortran routines
 * ----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_DJAC(long int* N, realtype* uudata , realtype* fdata,
                    realtype* jacdata, realtype* v1, realtype* v2,
                    int* ier);

#ifdef __cplusplus
}
#endif

/*
 * ----------------------------------------------------------------
 * Function : FKIN_DENSESETJAC
 * ----------------------------------------------------------------
 */

void FKIN_DENSESETJAC(int *flag, int *ier)
{
  if (*flag == 0) {
    *ier = KINSetJacFn(KIN_kinmem, NULL);
  }
  else {
    *ier = KINSetJacFn(KIN_kinmem, FKINDenseJac);
  }
  return;
}

/*
 * ----------------------------------------------------------------
 * Function : FKINDenseJac
 * ----------------------------------------------------------------
 * C function FKINDenseJac interfaces between KINSOL and a Fortran
 * subroutine FKDJAC for solution of a linear system with dense
 * Jacobian approximation. Addresses are passed to FKDJAC, using
 * the SUNDenseMatrix_Columns function. Auxiliary data is assumed
 * to be communicated by Common.
 * ----------------------------------------------------------------
 */

int FKINDenseJac(N_Vector uu, N_Vector fval, SUNMatrix J,
                 void *user_data, N_Vector vtemp1, N_Vector vtemp2)
{
  realtype *uu_data, *fval_data, *jacdata, *v1_data, *v2_data;
  long int N;
  int ier;

  /* Initialize all pointers to NULL */
  uu_data = fval_data = jacdata = v1_data = v2_data = NULL;

  /* NOTE: The user-supplied routine should set ier to an
     appropriate value, but we preset the value to zero
     (meaning SUCCESS) so the user need only reset the
     value if an error occurred */
  ier = 0;

  /* Get pointers to vector data */
  uu_data   = N_VGetArrayPointer(uu);
  fval_data = N_VGetArrayPointer(fval);
  v1_data   = N_VGetArrayPointer(vtemp1);
  v2_data   = N_VGetArrayPointer(vtemp2);

  N       = SUNDenseMatrix_Columns(J);
  jacdata = SUNDenseMatrix_Column(J,0);

  /* Call user-supplied routine */
  FK_DJAC(&N, uu_data, fval_data, jacdata, v1_data, v2_data, &ier);

  return(ier);
}
