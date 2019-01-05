/*-----------------------------------------------------------------
 * Programmer(s): Carol Woodward @ LLNL
 *                Daniel R. Reynolds @ SMU
 *-----------------------------------------------------------------
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
 *-----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "fida.h"
#include "ida_impl.h"
#include <ida/ida_ls.h>
#include <sunmatrix/sunmatrix_sparse.h>

/*=============================================================*/

/* Prototype of the Fortran routine */
 
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
extern void FIDA_SPJAC(realtype *T, realtype *CJ, realtype *Y, 
		       realtype *YP, realtype *R, long int *N,
                       long int *NNZ, realtype *JDATA,
                       sunindextype *JRVALS, sunindextype *JCPTRS,
                       realtype *H, long int *IPAR, realtype *RPAR, 
		       realtype *V1, realtype *V2, 
		       realtype *V3, int *ier);
 
#ifdef __cplusplus
}
#endif
 
/*=============================================================*/

/* Fortran interface to C routine IDASlsSetSparseJacFn; see
   fida.h for further information */
void FIDA_SPARSESETJAC(int *ier)
{
#if defined(SUNDIALS_INT32_T)
  IDAProcessError((IDAMem) IDA_idamem, IDA_ILL_INPUT, "IDA",
                  "FIDASPARSESETJAC", 
                  "Sparse Fortran users must configure SUNDIALS with 64-bit integers.");
  *ier = 1;
#else  
  *ier = IDASetJacFn(IDA_idamem, FIDASparseJac);
#endif
}

/*=============================================================*/
 
/* C interface to user-supplied Fortran routine FIDASPJAC; see 
   fida.h for additional information  */
int FIDASparseJac(realtype t, realtype cj, N_Vector y, N_Vector yp,
		  N_Vector fval, SUNMatrix J, void *user_data, 
		  N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *ypdata, *rdata, *v1data, *v2data, *v3data, *Jdata;
  realtype h;
  long int NP, NNZ; 
  sunindextype *indexvals, *indexptrs; 
  FIDAUserData IDA_userdata;

  IDAGetLastStep(IDA_idamem, &h);
  ydata   = N_VGetArrayPointer(y);
  ypdata  = N_VGetArrayPointer(yp);
  rdata   = N_VGetArrayPointer(fval);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  IDA_userdata = (FIDAUserData) user_data;
  NP = SUNSparseMatrix_NP(J);
  NNZ = SUNSparseMatrix_NNZ(J);
  Jdata = SUNSparseMatrix_Data(J);
  indexvals = SUNSparseMatrix_IndexValues(J);
  indexptrs = SUNSparseMatrix_IndexPointers(J);

  FIDA_SPJAC(&t, &cj, ydata, ypdata, rdata, &NP, &NNZ,
	    Jdata, indexvals, indexptrs, &h, 
	    IDA_userdata->ipar, IDA_userdata->rpar, v1data, 
	    v2data, v3data, &ier); 
  return(ier);
}

