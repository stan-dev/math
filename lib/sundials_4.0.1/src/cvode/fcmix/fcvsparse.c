/*
 * -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds and Ting Yan @ SMU
 *     Carol Woodward @ LLNL
 * -----------------------------------------------------------------
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
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include "fcvode.h"
#include "cvode_impl.h"
#include <cvode/cvode_ls.h>
#include <sunmatrix/sunmatrix_sparse.h>

/* Prototype of the Fortran routine */
 
#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
 
extern void FCV_SPJAC(realtype *T, realtype *Y, 
                      realtype *FY, long int *N,
                      long int *NNZ, realtype *JDATA,
                      sunindextype *JRVALS, 
                      sunindextype *JCPTRS, realtype *H, 
		      long int *IPAR, realtype *RPAR, 
                      realtype *V1, realtype *V2, 
                      realtype *V3, int *ier);
 
#ifdef __cplusplus
}
#endif
 
/*=============================================================*/

/* Fortran interface to C routine CVSlsSetSparseJacFn; see
   fcvode.h for further information */
void FCV_SPARSESETJAC(int *ier)
{
#if defined(SUNDIALS_INT32_T)
  cvProcessError((CVodeMem) CV_cvodemem, CV_ILL_INPUT, "CVODE",
                  "FCVSPARSESETJAC", 
                  "Sparse Fortran users must configure SUNDIALS with 64-bit integers.");
  *ier = 1;
#else  
  *ier = CVodeSetJacFn(CV_cvodemem, FCVSparseJac);
#endif
}

/*=============================================================*/
 
/* C interface to user-supplied Fortran routine FCVSPJAC; see 
   fcvode.h for additional information  */
int FCVSparseJac(realtype t, N_Vector y, N_Vector fy, 
		 SUNMatrix J, void *user_data, N_Vector vtemp1, 
		 N_Vector vtemp2, N_Vector vtemp3)
{
  int ier;
  realtype *ydata, *fydata, *v1data, *v2data, *v3data, *Jdata;
  realtype h;
  long int NP, NNZ; 
  sunindextype *indexvals, *indexptrs; 
  FCVUserData CV_userdata;

  CVodeGetLastStep(CV_cvodemem, &h);
  ydata   = N_VGetArrayPointer(y);
  fydata  = N_VGetArrayPointer(fy);
  v1data  = N_VGetArrayPointer(vtemp1);
  v2data  = N_VGetArrayPointer(vtemp2);
  v3data  = N_VGetArrayPointer(vtemp3);
  CV_userdata = (FCVUserData) user_data;
  NP = SUNSparseMatrix_NP(J);
  NNZ = SUNSparseMatrix_NNZ(J);
  Jdata = SUNSparseMatrix_Data(J);
  indexvals = SUNSparseMatrix_IndexValues(J);
  indexptrs = SUNSparseMatrix_IndexPointers(J);

  FCV_SPJAC(&t, ydata, fydata, &NP, &NNZ, Jdata, indexvals, 
	    indexptrs, &h, CV_userdata->ipar, CV_userdata->rpar,
            v1data, v2data, v3data, &ier); 
  return(ier);
}

