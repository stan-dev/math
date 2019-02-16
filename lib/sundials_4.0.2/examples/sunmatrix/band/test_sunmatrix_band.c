/*
 * ----------------------------------------------------------------- 
 * Programmer(s): Daniel Reynolds @ SMU
 *                David Gardner @ LLNL
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
 * This is the testing routine to check the SUNMatrix Band module 
 * implementation. 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_band.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include "test_sunmatrix.h"


/* ----------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int          fails = 0;            /* counter for test failures  */
  sunindextype cols, uband, lband;   /* matrix columns, bandwidths */
  SUNMatrix    A, I;                 /* test matrices              */
  N_Vector     x, y;                 /* test vectors               */
  int          print_timing;
  sunindextype i, j, k, kstart, kend, jstart, jend;
  realtype     *colj, *xdata, *ydata;

  /* check input and set vector length */
  if (argc < 5){
    printf("ERROR: FOUR (4) Inputs required: matrix cols, matrux uband, matrix lband, print timing \n");
    return(-1);
  }

  cols = atol(argv[1]); 
  if (cols <= 0) {
    printf("ERROR: number of matrix columns must be a positive integer \n");
    return(-1); 
  }

  uband = atol(argv[2]); 
  if ((uband <= 0) || (uband >= cols)){
    printf("ERROR: matrix upper bandwidth must be a positive integer, less than number of columns \n");
    return(-1); 
  }

  lband = atol(argv[3]); 
  if ((lband <= 0) || (lband >= cols)){
    printf("ERROR: matrix lower bandwidth must be a positive integer, less than number of columns \n");
    return(-1); 
  }

  print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  printf("\nBand matrix test: size %ld, bandwidths %ld %ld\n\n",
         (long int) cols, (long int) uband, (long int) lband);

  /* Initialize vectors and matrices to NULL */
  x = NULL;
  y = NULL;
  A = NULL;
  I = NULL;
  
  /* Create matrices and vectors */
  A = SUNBandMatrix(cols, uband, lband);
  I = SUNBandMatrix(cols, 0, 0);
  x = N_VNew_Serial(cols);
  y = N_VNew_Serial(cols);

  /* Fill matrices */
  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);
  for (j=0; j<cols; j++) {
    
    /* identity matrix */
    colj = SUNBandMatrix_Column(I, j);
    colj[0] = RCONST(1.0);

    /* A matrix */
    colj = SUNBandMatrix_Column(A, j);
    kstart = (j<uband) ? -j : -uband;
    kend = (j>cols-1-lband) ? cols-1-j: lband;
    for (k=kstart; k<=kend; k++)
      colj[k] = j - k;    /* A(i,j) = j + (j-i) */
  }

  /* Fill vectors */
  for (i=0; i<cols; i++) {

    /* x vector */
    xdata[i] = i;

    /* y vector */
    ydata[i] = RCONST(0.0);
    jstart = SUNMAX(0, i-lband);
    jend = SUNMIN(cols-1, i+uband);
    for (j=jstart; j<=jend; j++)
      ydata[i] += (j+j-i)*(j);
  }

  /* Run Tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_BAND, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  fails += Test_SUNMatScaleAdd(A, I, 0);
  fails += Test_SUNMatScaleAddI(A, I, 0);
  fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);

  /* Print result */
  if (fails) {
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
    printf("\nA =\n");
    SUNBandMatrix_Print(A,stdout);
    printf("\nI =\n");
    SUNBandMatrix_Print(I,stdout);
    printf("\nx =\n");
    N_VPrint_Serial(x);
    printf("\ny =\n");
    N_VPrint_Serial(y);
  } else {
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");
  }

  /* Free matrices and vectors */
  SUNMatDestroy(A);
  SUNMatDestroy(I);
  N_VDestroy(x);
  N_VDestroy(y);

  return(fails);
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  sunindextype i, j, istart, iend;
  realtype *Acolj, *Bcolj;

  /* check matrix type and dimensions */
  if (SUNMatGetID(A) != SUNMatGetID(B))
    return 1;
  if (SUNBandMatrix_Columns(A) != SUNBandMatrix_Columns(B))
    return 1;
  if (SUNBandMatrix_Rows(A) != SUNBandMatrix_Rows(B))
    return 1;
  if (SUNBandMatrix_LowerBandwidth(A) != SUNBandMatrix_LowerBandwidth(B))
    return 1;
  if (SUNBandMatrix_UpperBandwidth(A) != SUNBandMatrix_UpperBandwidth(B))
    return 1;
  
  /* check matrix data */
  for (j=0; j<SUNBandMatrix_Columns(A); j++) {

    /* access matrix columns */
    Acolj = SUNBandMatrix_Column(A, j);
    Bcolj = SUNBandMatrix_Column(B, j);

    /* compare entries in this column */
    istart = (j<SUNBandMatrix_UpperBandwidth(A)) ? -j : -SUNBandMatrix_UpperBandwidth(A);
    iend = (j>SUNBandMatrix_Columns(A)-1-SUNBandMatrix_LowerBandwidth(A)) ?
      SUNBandMatrix_Columns(A)-1-j : SUNBandMatrix_LowerBandwidth(A);
    for (i=istart; i<=iend; i++)
      failure += FNEQ(Acolj[i], Bcolj[i], tol);
  }

  if (failure > ZERO) {
    printf("check_matrix failure, A = \n");
    SUNBandMatrix_Print(A, stdout);
    printf("B = \n");
    SUNBandMatrix_Print(B, stdout);
  }

    
  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  sunindextype i, j, istart, iend;
  realtype *Acolj;

  /* check matrix data */
  for (j=0; j<SUNBandMatrix_Columns(A); j++) {

    /* access matrix column */
    Acolj = SUNBandMatrix_Column(A, j);

    /* compare entries in this column */
    istart = (j<SUNBandMatrix_UpperBandwidth(A)) ? -j : -SUNBandMatrix_UpperBandwidth(A);
    iend = (j>SUNBandMatrix_Columns(A)-1-SUNBandMatrix_LowerBandwidth(A)) ?
      SUNBandMatrix_Columns(A)-1-j : SUNBandMatrix_LowerBandwidth(A);
    for (i=istart; i<=iend; i++) {
      if (FNEQ(Acolj[i], val, tol)) {
        failure++;
        printf("j = %li, Acolj[%li] = %g, val = %g\n",
               (long int) j, (long int) i, Acolj[i], val);
      }
    }
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_vector(N_Vector X, N_Vector Y, realtype tol)
{
  int failure = 0;
  sunindextype i, local_length;
  realtype *Xdata, *Ydata;
  
  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);
  local_length = N_VGetLength_Serial(X);
  
  /* check vector data */
  for(i=0; i < local_length; i++) {
    failure += FNEQ(Xdata[i], Ydata[i], tol);
  }
  
  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(SUNMatrix A)
{
  realtype *Adata = SUNBandMatrix_Data(A);
  if (Adata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  return SUNTRUE;
}

