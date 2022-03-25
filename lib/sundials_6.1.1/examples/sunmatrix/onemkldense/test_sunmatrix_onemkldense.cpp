/* ---------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * ---------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ---------------------------------------------------------------------------
 * This is the testing routine to check the oneMKL dense SUNMatrix
 * implementation.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_onemkldense.h>
#include <sunmemory/sunmemory_sycl.h>
#include <nvector/nvector_sycl.h>
#include "test_sunmatrix.h"

/* ---------------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * ---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{
  int          fails = 0;        // counter for test failures
  sunindextype i, j, k, m, n;
  SUNContext   sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  // Check inputs and set matrix dimensions
  if (argc < 5)
  {
    printf("ERROR: FOUR (4) Input required: matrix rows, matrix cols, number of matrix blocks, print timing \n");
    return -1;
  }

  // Matrix block dimensions
  sunindextype matrows = (sunindextype) atol(argv[1]);
  if (matrows <= 0)
  {
    printf("ERROR: number of rows must be a positive integer \n");
    return -1;
  }

  sunindextype matcols = (sunindextype) atol(argv[2]);
  if (matcols <= 0)
  {
    printf("ERROR: number of cols must be a positive integer \n");
    return -1;
  }

  // Number of matrix blocks
  sunindextype nblocks = (sunindextype) atol(argv[3]);
  if (nblocks <= 0)
  {
    printf("ERROR: number of blocks must be a positive integer \n");
    return -1;
  }

  // Timing flag
  int print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  int square = (matrows == matcols) ? 1 : 0;
  printf("\noneMKL dense matrix test: size %ld by %ld\n\n",
         (long int) matrows, (long int) matcols);

  // Create an in-order GPU queue
  sycl::gpu_selector selector;
  sycl::queue myQueue(selector,
                      sycl::property_list{sycl::property::queue::in_order{}});

  sycl::device dev = myQueue.get_device();
  std::cout << "Running on "
            << (dev.get_info<sycl::info::device::name>())
            << std::endl;
  std::cout << " is host? "
            << (dev.is_host() ? "Yes" : "No")
            << std::endl;
  std::cout << " is cpu? "
            << (dev.is_cpu() ? "Yes" : "No")
            << std::endl;
  std::cout << " is gpu? "
            << (dev.is_gpu() ? "Yes" : "No")
            << std::endl;
  std::cout << " is accelerator? "
            << (dev.is_accelerator() ? "Yes" : "No")
            << std::endl;
  std::cout << " is the queue in order? "
            << (myQueue.is_in_order() ? "Yes" : "No")
            << std::endl;
  std::cout << " supports usm host allocations? "
            << (dev.get_info<sycl::info::device::usm_host_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " supports usm device allocations? "
            << (dev.get_info<sycl::info::device::usm_device_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " suports usm shared allocations? "
            << (dev.get_info<sycl::info::device::usm_shared_allocations>() ?
                "Yes" : "No")
            << std::endl;
  std::cout << " max work group size: "
            << dev.get_info<sycl::info::device::max_work_group_size>()
            << std::endl;
  std::cout << " max global memory size (bytes): "
            << dev.get_info<sycl::info::device::global_mem_size>()
            << std::endl;
  std::cout << " max local memory size (bytes): "
            << dev.get_info<sycl::info::device::local_mem_size>()
            << std::endl;
  std::cout << std::endl;

  // Create Sycl memory helper
  SUNMemoryHelper memhelper = SUNMemoryHelper_Sycl(sunctx);
  if (!memhelper)
  {
    printf("Memory helper creation failed\n");
    return 1;
  }

  // Create vectors and matrices
  N_Vector x = N_VNew_Sycl(matcols * nblocks, &myQueue, sunctx);
  if (!x)
  {
    printf("Vector creation failed\n");
    return 1;
  }

  N_Vector y = N_VNew_Sycl(matrows * nblocks, &myQueue, sunctx);
  if (!y)
  {
    printf("Vector creation failed\n");
    N_VDestroy(x);
    return 1;
  }

  SUNMatrix A = SUNMatrix_OneMklDenseBlock(nblocks, matrows, matcols,
                                           SUNMEMTYPE_DEVICE, memhelper,
                                           &myQueue, sunctx);
  if (!A)
  {
    printf("Matrix creation failed\n");
    N_VDestroy(x);
    N_VDestroy(y);
  }

  myQueue.wait_and_throw();

  SUNMatrix I = NULL;
  if (square)
  {
    I = SUNMatClone(A);
    if (!I)
    {
      printf("Matrix creation failed\n");
      N_VDestroy(x);
      N_VDestroy(y);
      SUNMatDestroy(A);
    }
  }

  // Allocate host data
  realtype* Adata = (realtype*) malloc(sizeof(realtype)*SUNMatrix_OneMklDense_LData(A));
  if (!Adata)
  {
    printf("Data allocation failed\n");
    N_VDestroy(x);
    N_VDestroy(y);
    SUNMatDestroy(A);
  }

  realtype* Idata = NULL;
  if (square)
  {
    Idata = (realtype*) malloc(sizeof(realtype)*SUNMatrix_OneMklDense_LData(I));
    if (!Idata)
    {
      printf("Data allocation failed\n");
      N_VDestroy(x);
      N_VDestroy(y);
      SUNMatDestroy(A);
      free(Adata);
    }
  }

  // Fill matrices and vectors
  for(k=0; k < nblocks; k++)
    for(j=0; j < matcols; j++)
      for(i=0; i < matrows; i++)
        Adata[k * matcols * matrows + j * matrows + i] = (j + 1) * (i + j);

  myQueue.wait_and_throw();

  SUNMatrix_OneMklDense_CopyToDevice(A, Adata);

  if (square)
  {
    for(k = 0; k < nblocks; k++)
      for(j = 0; j < matcols; j++)
        for(i = 0; i < matrows; i++)
          Idata[k * matcols * matrows + j * matrows + i] =
            (j == i) ? ONE : ZERO;

    myQueue.wait_and_throw();

    SUNMatrix_OneMklDense_CopyToDevice(I, Idata);
  }

  realtype* xdata = N_VGetArrayPointer(x);
  for(k = 0; k < nblocks; k++)
    for(i = 0; i < matcols; i++)
      xdata[matcols * k + i] = ONE / (i+1);

  myQueue.wait_and_throw();

  N_VCopyToDevice_Sycl(x);

  realtype* ydata = N_VGetArrayPointer(y);
  for(k = 0; k < nblocks; k++)
  {
    for(i = 0; i < matrows; i++)
    {
      m = i;
      n = m + matcols - 1;
      ydata[matrows * k + i] = HALF*(n+1-m)*(n+m);
    }
  }

  myQueue.wait_and_throw();

  N_VCopyToDevice_Sycl(y);

  // SUNMatrix Tests
  fails += Test_SUNMatGetID(A, SUNMATRIX_ONEMKLDENSE, 0);
  fails += Test_SUNMatClone(A, 0);
  fails += Test_SUNMatCopy(A, 0);
  fails += Test_SUNMatZero(A, 0);
  if (square)
  {
    fails += Test_SUNMatScaleAdd(A, I, 0);
    fails += Test_SUNMatScaleAddI(A, I, 0);
  }
  fails += Test_SUNMatMatvecSetup(A, 0);
  fails += Test_SUNMatMatvec(A, x, y, 0);
  fails += Test_SUNMatSpace(A, 0);

  // Print result
  if (fails)
    printf("FAIL: SUNMatrix module failed %i tests \n \n", fails);
  else
    printf("SUCCESS: SUNMatrix module passed all tests \n \n");

  // Free vectors and matrices
  N_VDestroy(x);
  N_VDestroy(y);
  SUNMatDestroy(A);
  free(Adata);
  if (square)
  {
    SUNMatDestroy(I);
    free(Idata);
  }
  SUNMemoryHelper_Destroy(memhelper);
  SUNContext_Free(&sunctx);

  return fails;
}


/* ---------------------------------------------------------------------------
 * Check matrix
 * ---------------------------------------------------------------------------*/


int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  sunindextype i = 0;
  sunindextype Aldata = SUNMatrix_OneMklDense_LData(A);
  sunindextype Bldata = SUNMatrix_OneMklDense_LData(B);
  realtype *Adata = (realtype*) malloc(sizeof(realtype) * Aldata);
  realtype *Bdata = (realtype*) malloc(sizeof(realtype) * Bldata);

  // Copy data to host
  SUNMatrix_OneMklDense_CopyFromDevice(A, Adata);
  SUNMatrix_OneMklDense_CopyFromDevice(B, Bdata);

  // Check lengths
  if (Aldata != Bldata)
  {
    printf(">>> ERROR: check_matrix: Different data array lengths \n");
    return 1;
  }

  // Compare data
  for(i = 0; i < Aldata; i++)
  {
    failure += SUNRCompareTol(Adata[i], Bdata[i], tol);
  }

  free(Adata);
  free(Bdata);

  if (failure > ZERO)
    return 1;
  else
    return 0;
}


int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  sunindextype i = 0;
  sunindextype Aldata = SUNMatrix_OneMklDense_LData(A);
  realtype *Adata = (realtype*) malloc(sizeof(realtype) * Aldata);

  // copy data to host
  SUNMatrix_OneMklDense_CopyFromDevice(A, Adata);

  // compare data
  for(i = 0; i < Aldata; i++)
  {
    int check = SUNRCompareTol(Adata[i], val, tol);
    if (check)
    {
      printf("failed at %ld\n", (long int) i);
      failure += check;
    }
  }

  free(Adata);

  if (failure > ZERO)
    return 1;
  else
    return 0;
}


int check_vector(N_Vector actual, N_Vector expected, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  // copy vectors to host
  N_VCopyFromDevice_Sycl(actual);
  N_VCopyFromDevice_Sycl(expected);

  // get vector data
  xdata = N_VGetArrayPointer(actual);
  ydata = N_VGetArrayPointer(expected);

  // check data lengths
  xldata = N_VGetLength(actual);
  yldata = N_VGetLength(expected);

  if (xldata != yldata)
  {
    printf(">>> ERROR: check_vector: Different data array lengths \n");
    return 1;
  }

  // check vector data
  for(i = 0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO)
  {
    printf("Check_vector failures:\n");
    for(i = 0; i < xldata; i++)
      if (SUNRCompareTol(xdata[i], ydata[i], tol) != 0)
        printf("  actual[%ld] = %g != %e (err = %g)\n", (long int) i,
               xdata[i], ydata[i], SUNRabs(xdata[i]-ydata[i]));
  }

  if (failure > ZERO)
    return 1;
  else
    return 0;
}


booleantype has_data(SUNMatrix A)
{
  realtype *Adata = SUNMatrix_OneMklDense_Data(A);
  if (Adata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}


booleantype is_square(SUNMatrix A)
{
  if (SUNMatrix_OneMklDense_Rows(A) == SUNMatrix_OneMklDense_Columns(A))
    return SUNTRUE;
  else
    return SUNFALSE;
}


void sync_device(SUNMatrix A)
{
  ((SUNMatrixContent_OneMklDense)(A->content))->queue->wait_and_throw();
}
