/* -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh and
 *                Radu Serban @ LLNL
 *                David J. Gardner @ LLNL
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
 * Routines used to interface between KINSOL and a Fortran
 * user-supplied routine FKJTIMES (Jacobian J times vector v).
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"
#include "kinsol_impl.h"

#include <kinsol/kinsol_ls.h>

/*------------------------------------------------------------------
  prototype of the user-supplied fortran routine
  ------------------------------------------------------------------*/
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

extern void FK_JTIMES(realtype* vdata, realtype* Jvdata, int* new_uu,
                      realtype* uudata, int* ier);

#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------
  Function : FKIN_LSSETJAC
  ------------------------------------------------------------------*/
void FKIN_LSSETJAC(int *flag, int *ier)
{
  if ((*flag) == 0) KINSetJacTimesVecFn(KIN_kinmem, NULL);
  else              KINSetJacTimesVecFn(KIN_kinmem, FKINJtimes);

  return;
}

/*------------------------------------------------------------------
  Function : FKIN_SPILSSETJAC -- DEPRECATED
  ------------------------------------------------------------------*/
void FKIN_SPILSSETJAC(int *flag, int *ier)
{ FKIN_LSSETJAC(flag, ier); }

/*------------------------------------------------------------------
  Function : FKINJtimes
  ------------------------------------------------------------------
  C function FKINJtimes is used to interface between
  KINSp* / KINSp*JTimes and FK_JTIMES (user-supplied Fortran
  routine).
  ------------------------------------------------------------------*/
int FKINJtimes(N_Vector v, N_Vector Jv,
               N_Vector uu, booleantype *new_uu, 
               void *user_data)
{
  int retcode;
  realtype *vdata, *Jvdata, *uudata;

  vdata = Jvdata = uudata = NULL;

  vdata  = N_VGetArrayPointer(v);
  uudata = N_VGetArrayPointer(uu);
  Jvdata = N_VGetArrayPointer(Jv);
 
  FK_JTIMES(vdata, Jvdata, (int *) new_uu, uudata, &retcode);

  return(retcode);
}
