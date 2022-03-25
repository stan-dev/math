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
 * This is the testing routine to check the SUNLinSol Dense module
 * implementation.
 * ---------------------------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>

#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_onemkldense.h>
#include <sunmatrix/sunmatrix_onemkldense.h>
#include <sunmemory/sunmemory_sycl.h>
#include <nvector/nvector_sycl.h>
#include "test_sunlinsol.h"

/* ---------------------------------------------------------------------------
 * SUNLinSol_OneMklDense Testing Routine
 * ---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{
  int          fails = 0; // counter for test failures
  sunindextype i, j, k;
  SUNContext   sunctx;

  if (SUNContext_Create(NULL, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  // Check inputs and set matrix dimensions
  if (argc < 4){
    printf("ERROR: THREE (3) Inputs required: matrix cols, number of blocks, print timing \n");
    return -1;
  }

  // Number of matrix columns and rows
  sunindextype cols = (sunindextype) atol(argv[1]);
  if (cols <= 0) {
    printf("ERROR: number of matrix columns must be a positive integer \n");
    return -1;
  }
  sunindextype rows = cols;

  // Number of matrix blocks
  sunindextype nblocks = (sunindextype) atol(argv[2]);
  if (nblocks <= 0) {
    printf("ERROR: number of blocks must be a positive integer \n");
    return -1;
  }

  // Timing flag
  int print_timing = atoi(argv[3]);
  SetTiming(print_timing);

  printf("\noneMKL dense linear solver test: size %ld, blocks %ld\n\n",
         (long int) cols, (long int) nblocks);

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
  N_Vector x = N_VNew_Sycl(cols * nblocks, &myQueue, sunctx);
  if (!x)
  {
    printf("Vector creation failed\n");
    return 1;
  }

  N_Vector y = N_VClone(x);
  if (!y)
  {
    printf("Vector creation failed\n");
    N_VDestroy(x);
    return 1;
  }

  N_Vector b = N_VClone(x);
  if (!b)
  {
    printf("Vector creation failed\n");
    N_VDestroy(x);
    N_VDestroy(y);
    return 1;
  }

  SUNMatrix A;
  if (nblocks > 1)
  {
    A = SUNMatrix_OneMklDenseBlock(nblocks, rows, cols, SUNMEMTYPE_DEVICE,
                                   memhelper, &myQueue, sunctx);
  }
  else
  {
    A = SUNMatrix_OneMklDense(rows, cols, SUNMEMTYPE_DEVICE, memhelper,
                              &myQueue, sunctx);
  }

  if (!A)
  {
    printf("Matrix creation failed\n");
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
  }

  SUNMatrix B = SUNMatClone(A);
  if (!B)
  {
    printf("Matrix creation failed\n");
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
    SUNMatDestroy(A);
  }

  SUNMatrix I = SUNMatClone(A);
  if (!I)
  {
    printf("Matrix creation failed\n");
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
    SUNMatDestroy(A);
    SUNMatDestroy(B);
  }

  // Allocate host data
  realtype* Adata = (realtype*) malloc(sizeof(realtype) *
                                       SUNMatrix_OneMklDense_LData(A));
  if (!Adata)
  {
    printf("Data allocation failed\n");
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
    SUNMatDestroy(A);
    SUNMatDestroy(B);
    SUNMatDestroy(I);
  }

  realtype* Idata = (realtype*) malloc(sizeof(realtype) *
                                       SUNMatrix_OneMklDense_LData(I));
  if (!Idata)
  {
    printf("Data allocation failed\n");
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
    SUNMatDestroy(A);
    SUNMatDestroy(B);
    SUNMatDestroy(I);
    free(Adata);
  }

  // Fill A matrix with uniform random data in [0,1/cols]
  for (k = 0; k < nblocks; k++)
    for (j = 0; j < cols; j++)
      for (i = 0; i < rows; i++)
        Adata[k * cols * rows + j * rows + i] =
          (realtype) rand() / (realtype) RAND_MAX / cols;

  // Create anti-identity matrix
  for (k = 0; k < nblocks; k++)
    for(j = 0; j < cols; j++)
      for (i = 0; i < rows; i++)
        Idata[k * cols * rows + j * rows + i] =
          ((rows-1-i) == j) ? RCONST(1.0) : RCONST(0.0);

  // Add anti-identity to ensure the solver needs to do row-swapping
  for (k = 0; k < nblocks; k++)
    for (i = 0; i < rows; i++)
      for(j = 0; j < cols; j++)
        Adata[k * cols * rows + j * rows + i] +=
          Idata[k * cols * rows + j * rows + i];

  SUNMatrix_OneMklDense_CopyToDevice(A, Adata);
  SUNMatrix_OneMklDense_CopyToDevice(I, Idata);

  // Fill x vector with uniform random data in [0,1]
  realtype* xdata = N_VGetArrayPointer(x);
  for (j = 0; j < cols * nblocks; j++)
    xdata[j] = (realtype) rand() / (realtype) RAND_MAX;

  N_VCopyToDevice_Sycl(x);

  // copy A and x into B and y to print in case of solver failure
  SUNMatCopy(A, B);
  N_VScale(ONE, x, y);

  // create right-hand side vector for linear solve
  fails += SUNMatMatvecSetup(A);
  fails += SUNMatMatvec(A, x, b);
  if (fails)
  {
    printf("FAIL: SUNLinSol SUNMatMatvec failure\n");

    // Free matrices and vectors
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
    SUNMatDestroy(A);
    SUNMatDestroy(B);
    SUNMatDestroy(I);
    free(Adata);
    free(Idata);
    return 1;
  }

  // Create dense linear solver
  SUNLinearSolver LS = SUNLinSol_OneMklDense(x, A, sunctx);
  if (!LS)
  {
    printf("FAIL: SUNLinSol_OneMklDense failure\n");

    // Free matrices and vectors
    N_VDestroy(x);
    N_VDestroy(y);
    N_VDestroy(b);
    SUNMatDestroy(A);
    SUNMatDestroy(B);
    SUNMatDestroy(I);
    free(Adata);
    free(Idata);

    return 1;
  }

  // Run Tests
  fails += Test_SUNLinSolInitialize(LS, 0);
  fails += Test_SUNLinSolSetup(LS, A, 0);
  fails += Test_SUNLinSolSolve(LS, A, x, b, RCONST(1e-10), SUNTRUE, 0);
  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_ONEMKLDENSE, 0);
  fails += Test_SUNLinSolLastFlag(LS, 0);
  fails += Test_SUNLinSolSpace(LS, 0);

  // Print result
  if (fails)
  {
    printf("FAIL: SUNLinSol module failed %i tests \n \n", fails);
    printf("\nx (original) =\n");
    N_VCopyFromDevice_Sycl(y);
    N_VPrint(y);
    printf("\nx (computed) =\n");
    N_VCopyFromDevice_Sycl(x);
    N_VPrint(x);
    printf("\nb =\n");
    N_VCopyFromDevice_Sycl(b);
    N_VPrint(b);
  }
  else
  {
    printf("SUCCESS: SUNLinSol module passed all tests \n \n");
  }

  // Free solver, matrix and vectors
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  SUNMatDestroy(B);
  SUNMatDestroy(I);
  N_VDestroy(x);
  N_VDestroy(y);
  N_VDestroy(b);
  free(Adata);
  free(Idata);
  SUNMemoryHelper_Destroy(memhelper);
  SUNContext_Free(&sunctx);

  return fails;
}


/* ---------------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * ---------------------------------------------------------------------------*/


int check_vector(N_Vector X, N_Vector Y, realtype tol)
{
  int failure = 0;
  sunindextype i = 0;
  sunindextype local_length = N_VGetLength(X);
  realtype* Xdata = N_VGetArrayPointer(X);
  realtype* Ydata = N_VGetArrayPointer(Y);

  // Copy data to host
  N_VCopyFromDevice_Sycl(X);
  N_VCopyFromDevice_Sycl(Y);

  // Check vector data
  for (i = 0; i < local_length; i++)
    failure += SUNRCompareTol(Xdata[i], Ydata[i], tol);

  if (failure > ZERO)
  {
    realtype maxerr = ZERO;
    for(i = 0; i < local_length; i++)
      maxerr = SUNMAX(SUNRabs(Xdata[i] - Ydata[i]), maxerr);
    printf("check err failure: maxerr = %g (tol = %g)\n", maxerr, tol);
    return 1;
  }
  else
  {
    return 0;
  }
}

void sync_device()
{
}
