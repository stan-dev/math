/* -----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 *                Daniel R. Reynolds @ SMU
 *                David J. Gardner @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2015, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Lawrence Livermore National Laboratory in part under 
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "fkinsol.h"
#include "kinsol_impl.h"

#include <kinsol/kinsol_ls.h>
#include <sunmatrix/sunmatrix_sparse.h>

/*=============================================================*/

/* Prototype of the Fortran routine */
 
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
extern void FKIN_SPJAC(realtype *Y, realtype *FY, long int *N,
                       long int *NNZ, realtype *JDATA,
                       sunindextype *JRVALS, sunindextype *JCPTRS,
                       realtype *V1, realtype *V2, int *ier);
 
#ifdef __cplusplus
}
#endif
 
/*=============================================================*/

/* Fortran interface to C routine KINSlsSetSparseJacFn; see
   fkinsol.h for further information */
void FKIN_SPARSESETJAC(int *ier)
{
#if defined(SUNDIALS_INT32_T)
  KINProcessError((KINMem) KIN_kinmem, KIN_ILL_INPUT, "KIN",
                  "FKINSPARSESETJAC",
                  "Sparse Fortran users must configure SUNDIALS with 64-bit integers.");
  *ier = 1;
#else
  *ier = KINSetJacFn(KIN_kinmem, FKINSparseJac);
#endif
}

/*=============================================================*/
 
/* C interface to user-supplied Fortran routine FKINSPJAC; see 
   fkinsol.h for additional information  */
int FKINSparseJac(N_Vector y, N_Vector fy, SUNMatrix J,
                  void *user_data, N_Vector vtemp1,
                  N_Vector vtemp2)
{
  int ier;
  realtype *ydata, *fydata, *v1data, *v2data, *Jdata;
  long int NP, NNZ;
  sunindextype *indexvals, *indexptrs;
 
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);

  NP = SUNSparseMatrix_NP(J);
  NNZ = SUNSparseMatrix_NNZ(J);
  Jdata = SUNSparseMatrix_Data(J);
  indexvals = SUNSparseMatrix_IndexValues(J);
  indexptrs = SUNSparseMatrix_IndexPointers(J);
 
  FKIN_SPJAC(ydata, fydata, &NP, &NNZ,
             Jdata, indexvals, indexptrs,
             v1data, v2data, &ier);
  return(ier);
}

