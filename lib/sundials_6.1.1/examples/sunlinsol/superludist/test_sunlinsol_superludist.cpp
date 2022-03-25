/*
 * -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the testing routine to check the SUNLinSol SuperLUDIST
 * module implementation.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sunlinsol/sunlinsol_superludist.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_slunrloc.h>
#include <nvector/nvector_parallel.h>
#include <nvector/nvector_serial.h>
#include "test_sunlinsol.h"

/* 1 --> print extended output from each process to a file instead
 * of to stdout */
#define LOG_PROCESS_TO_FILE 0

/* helper functions */
int csr_from_dense(SUNMatrix Ad, realtype droptol, realtype **matdata,
                   sunindextype **colind, sunindextype **rowptrs);

/* ----------------------------------------------------------------------
 * SUNSuperLUDIST Linear Solver Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int             fails = 0;            /* local test fail counter    */
  int             globfails = 0;        /* global test fail counter   */
  sunindextype    N;                    /* matrix columns, rows       */
  sunindextype    M_loc, M_rem;         /* local matrix rows          */
  sunindextype    Mnprocs;              /* M/procs                    */
  sunindextype    NNZ_local;            /* local number of nonzeros   */
  sunindextype    fst_row;              /* row index in global matrix */
  SUNLinearSolver LS;                   /* linear solver object       */
  SUNMatrix       A, D;                 /* test matrices              */
  SuperMatrix     *Asuper;              /* SuperMatrices of A         */
  N_Vector        gx, gy, gb, x, y, b;  /* test vectors               */
  realtype        *matdata;             /* underlying data arrays     */
  realtype        *xdata, *ydata, *bdata;
  sunindextype    *rowptrs, *colind;
  gridinfo_t      grid;                 /* SuperLU-DIST process grid  */
  xLUstruct_t     lu;
  xScalePermstruct_t scaleperm;
  xSOLVEstruct_t  solve;                /* SuperLU-DIST solve struct  */
  SuperLUStat_t   stat;
  superlu_dist_options_t options;       /* SuperLU-DIST options struct*/
  sunindextype    i, j, k;

  int             nprocs;
  int             nprow, npcol;         /* process grid rows & cols   */
  int             print_timing;
  MPI_Status      mpistatus;
  MPI_Comm        comm;
  int             rank;
  FILE            *fp;
#if LOG_PROCESS_TO_FILE
  char            filename[64];
#endif
  SUNContext      sunctx;

  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;

  if (SUNContext_Create(&comm, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

#if LOG_PROCESS_TO_FILE
  /* create log files */
  sprintf(filename, "proc-%d.log", rank);
  fp = fopen(filename, "w+");
#else
  fp = stdout;
#endif

  /* check input and set matrix dimensions */
  if (argc < 5){
    if (rank == 0)
      printf("ERROR: FOUR (4) Inputs required: matrix size, nprow, npcol, print timing \n");
    return(1);
  }

  N = (sunindextype) atol(argv[1]);
  if (N <= 0) {
    if (rank == 0)
      printf("ERROR: matrix size must be a positive integer \n");
    return(1);
  }

  nprow = atoi(argv[2]);
  npcol = atoi(argv[3]);
  if (nprow <= 0 || npcol <= 0) {
    if (rank == 0)
      printf("ERROR: nprow and npcol must be positive integers \n");
    return(1);
  } else if (nprocs < (nprow*npcol)) {
    if (rank == 0)
      printf("ERROR: nprocs must be >= nprow*npcol \n");
    return(1);
  }

  if (rank == 0) {
    printf("\nSuperLUDIST linear solver test: \
           nprocs %d, nprow %d, npcol %d, size %ld\n\n",
           nprocs, nprow, npcol, (long int) N);
  }

  print_timing = atoi(argv[4]);
  SetTiming(print_timing);

  /* intiailize SuperLU-DIST process grid */
  superlu_gridinit(comm, nprow, npcol, &grid);
  /* excess processes just exit */
  if (grid.iam >= nprow*npcol) {
    superlu_gridexit(&grid);
    MPI_Finalize();
    return(0);
  }
  /* now update the number of procs available */
  nprocs = nprow*npcol;

  /* determine number of local rows and global row index of
   * the first row in the local matrix */
  M_loc   = Mnprocs = N/nprocs;
  M_rem   = N - M_loc*nprocs;
  fst_row = grid.iam*M_loc;
  /* if there are remaining rows, the last process gets them */
  if ((M_loc*nprocs) != N) {
    if (grid.iam == (nprocs-1))
      M_loc = N - M_loc*(nprocs-1);
  }

  /* initialize global matrices and vectors */
  D  = NULL;
  A  = NULL;
  x  = NULL;
  y  = NULL;
  b  = NULL;
  gx = NULL;
  gy = NULL;
  gb = NULL;

  if (grid.iam == 0) {

    /* Create matrices and vectors */
    D = SUNDenseMatrix(N, N, sunctx);

    /* Fill matrix with uniform random data in [0,1/N] */
    for (k=0; k<N; k++) {
      i = rand() % N;
      j = rand() % N;
      matdata = SUNDenseMatrix_Column(D,j);
      matdata[i] = (realtype) rand() / (realtype) RAND_MAX / N;
    }

    /* Add identity to matrix */
    fails = SUNMatScaleAddI(ONE, D);
    if (fails) {
      printf(">>> FAIL: SUNDenseMatrix SUNMatScaleAddI failure\n");
      return(1);
    }

    /* create the global vectors */
    gx = N_VNew_Serial(N, sunctx);
    gy = N_VNew_Serial(N, sunctx);
    gb = N_VNew_Serial(N, sunctx);
    xdata = N_VGetArrayPointer(gx);
    ydata = N_VGetArrayPointer(gy);
    bdata = N_VGetArrayPointer(gb);

    /* fill x with random data */
    for (i=0; i<N; i++)
      xdata[i] = (realtype) rand() / (realtype) RAND_MAX;

    /* copy x into y to print in case of solver failure */
    N_VScale(ONE, gx, gy);

    matdata = NULL; rowptrs = NULL; colind  = NULL;
    fails = csr_from_dense(D, ZERO, &matdata, &colind, &rowptrs);
    if (fails != 0) {
      printf(">>> FAIL: csr_from_dense failure \n");
      return(1);
    }

    /* distribute matrix and vectors */
    for (i=1; i<nprocs; i++) {
      sunindextype *coltemp, *rowtemp;
      realtype *datatemp;
      sunindextype fst_rowi, M_loci;

      M_loci   = (i == (nprocs-1)) ? M_loc + M_rem : Mnprocs;
      fst_rowi = Mnprocs*i;

      /* send number of local NNZ */
      NNZ_local = rowptrs[fst_rowi+M_loci] - rowptrs[fst_rowi];
      MPI_Send(&NNZ_local, 1, MPI_SUNINDEXTYPE, i, i, grid.comm);

      /* send out rowptrs */
      rowtemp = &rowptrs[fst_rowi];
      MPI_Send(rowtemp, M_loci+1, MPI_SUNINDEXTYPE, i, i, grid.comm);

      /* send corresponding column indices */
      coltemp = &colind[rowptrs[fst_rowi]];
      MPI_Send(coltemp, NNZ_local, MPI_SUNINDEXTYPE, i, i, grid.comm);

      /* send corresponding data */
      datatemp = &matdata[rowptrs[fst_rowi]];
      MPI_Send(datatemp, NNZ_local, MPI_SUNREALTYPE, i, i, grid.comm);

      /* send vector data */
      datatemp = &xdata[fst_rowi];
      MPI_Send(datatemp, M_loci, MPI_SUNREALTYPE, i, i, grid.comm);

      datatemp = &ydata[fst_rowi];
      MPI_Send(datatemp, M_loci, MPI_SUNREALTYPE, i, i, grid.comm);

      datatemp = &bdata[fst_rowi];
      MPI_Send(datatemp, M_loci, MPI_SUNREALTYPE, i, i, grid.comm);
    }

    NNZ_local = rowptrs[M_loc] - rowptrs[0];

    /* Create local SuperLU-DIST SLU_NR_loc SuperMatrix */
    Asuper = NULL;
    Asuper = (SuperMatrix*) SUPERLU_MALLOC(sizeof(SuperMatrix));

    /* Allocates Asuper->Store and sets structure members */
    dCreate_CompRowLoc_Matrix_dist(Asuper, N, N, NNZ_local, M_loc, fst_row,
                                   matdata, colind, rowptrs, SLU_NR_loc, SLU_D, SLU_GE);

    A = SUNMatrix_SLUNRloc(Asuper, &grid, sunctx);
    if (A == NULL) {
      fails++;
      printf("process %6d: FAIL: SUNMatrix_SLUNRloc returned NULL\n", grid.iam);
      Destroy_CompRowLoc_Matrix_dist(Asuper);
      SUNMatDestroy(D);
      N_VDestroy(gx); N_VDestroy(gy); N_VDestroy(gb);
      return(fails);
    }

    /* last entry in rowptrs should be the num of local nz */
    rowptrs[M_loc] = NNZ_local;

    /* make the local NVectors */
    x = N_VMake_Parallel(grid.comm, M_loc, N, xdata, sunctx);
    y = N_VMake_Parallel(grid.comm, M_loc, N, ydata, sunctx);
    b = N_VMake_Parallel(grid.comm, M_loc, N, bdata, sunctx);

  } else {

    sunindextype shift;

    /* recieve number of local nnz */
    MPI_Recv(&NNZ_local, 1, MPI_SUNINDEXTYPE, 0, grid.iam, grid.comm, &mpistatus);

    /* Allocate memory for matrix members */
    matdata = (realtype*) SUPERLU_MALLOC(NNZ_local*sizeof(realtype));
    colind  = (sunindextype*) SUPERLU_MALLOC(NNZ_local*sizeof(sunindextype));
    rowptrs = (sunindextype*) SUPERLU_MALLOC((M_loc+1)*sizeof(sunindextype));

    /* receive distributed matrix */
    MPI_Recv(rowptrs, M_loc+1, MPI_SUNINDEXTYPE, 0, grid.iam, grid.comm, &mpistatus);
    MPI_Recv(colind, NNZ_local, MPI_SUNINDEXTYPE, 0, grid.iam, grid.comm, &mpistatus);
    MPI_Recv(matdata, NNZ_local, MPI_SUNREALTYPE, 0, grid.iam, grid.comm, &mpistatus);

    /* localize rowptrs */
    shift = rowptrs[0];
    for (i=0; i<M_loc; i++)
      rowptrs[i] = rowptrs[i]-shift;
    rowptrs[M_loc] = NNZ_local;

    /* Create local SuperLU-DIST SuperMatrix */
    Asuper = NULL;
    Asuper = (SuperMatrix*) SUPERLU_MALLOC(sizeof(SuperMatrix));
    dCreate_CompRowLoc_Matrix_dist(Asuper, N, N, NNZ_local, M_loc, fst_row,
                                   matdata, colind, rowptrs, SLU_NR_loc, SLU_D, SLU_GE);

    /* Create local SuperLU-DIST SUNMatrix */
    A = SUNMatrix_SLUNRloc(Asuper, &grid, sunctx);
    if (A == NULL) {
      fails++;
      printf("process %6d: FAIL: SUNMatrix_SLUNRloc returned NULL\n", grid.iam);
      Destroy_CompRowLoc_Matrix_dist(Asuper);
      return(fails);
    }

    /* make the local NVectors */
    x = N_VNew_Parallel(grid.comm, M_loc, N, sunctx);
    y = N_VNew_Parallel(grid.comm, M_loc, N, sunctx);
    b = N_VNew_Parallel(grid.comm, M_loc, N, sunctx);
    xdata = N_VGetArrayPointer(x);
    ydata = N_VGetArrayPointer(y);
    bdata = N_VGetArrayPointer(b);

    /* recieve vectors */
    MPI_Recv(xdata, M_loc, MPI_SUNREALTYPE, 0, grid.iam, grid.comm, &mpistatus);
    MPI_Recv(ydata, M_loc, MPI_SUNREALTYPE, 0, grid.iam, grid.comm, &mpistatus);
    MPI_Recv(bdata, M_loc, MPI_SUNREALTYPE, 0, grid.iam, grid.comm, &mpistatus);
  }

  /* Initialize all of the SuperLU-DIST structures */
  set_default_options_dist(&options);
  xScalePermstructInit(N, N, &scaleperm);
  xLUstructInit(N, &lu);
  PStatInit(&stat);

  /* Dont print out stats in this test */
  options.PrintStat = NO;

  /* Create SuperLUDIST linear solver */
  LS = SUNLinSol_SuperLUDIST(y, A, &grid, &lu, &scaleperm, &solve, &stat, &options, sunctx);
  if (LS == NULL) {
    fails++;
    printf("process %6d: FAIL: SUNLinSol_SuperLUDIST returned NULL\n", grid.iam);
    Destroy_CompRowLoc_Matrix_dist(Asuper);
    SUNMatDestroy(A);
    N_VDestroy(x); N_VDestroy(y); N_VDestroy(b);
    if (rank == 0) {
      SUNMatDestroy(D);
      N_VDestroy(gx); N_VDestroy(gy); N_VDestroy(gb);
    }
    return(fails);
  }

  /* Setup A so that the Matvec can be used */
  fails = SUNMatMatvecSetup(A);
  if (fails) {
    printf("process %6d: FAIL: SUNMatrix_SLUNRloc SUNMatMatvecSetup failure\n", grid.iam);
    Destroy_CompRowLoc_Matrix_dist(Asuper);
    SUNMatDestroy(A);
    N_VDestroy(x); N_VDestroy(y); N_VDestroy(b);
    if (rank == 0) {
      SUNMatDestroy(D);
      N_VDestroy(gx); N_VDestroy(gy); N_VDestroy(gb);
    }
    return(fails);
  }
  /* create right-hand side vector for linear solve */
  fails = SUNMatMatvec(A, x, b);
  if (fails) {
    printf("process %6d: FAIL: SUNMatrix_SLUNRloc SUNMatMatvec failure\n", grid.iam);
    Destroy_CompRowLoc_Matrix_dist(Asuper);
    SUNMatDestroy(A);
    N_VDestroy(x); N_VDestroy(y); N_VDestroy(b);
    if (rank == 0) {
      SUNMatDestroy(D);
      N_VDestroy(gx); N_VDestroy(gy); N_VDestroy(gb);
    }
    return(fails);
  }

  /* Run Tests */
  fails += Test_SUNLinSolInitialize(LS, grid.iam);
  fails += Test_SUNLinSolSetup(LS, A, grid.iam);
  fails += Test_SUNLinSolSolve(LS, A, x, b, 100*UNIT_ROUNDOFF, SUNTRUE, grid.iam);

  fails += Test_SUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, grid.iam);
  fails += Test_SUNLinSolGetID(LS, SUNLINEARSOLVER_SUPERLUDIST, grid.iam);
  fails += Test_SUNLinSolLastFlag(LS, grid.iam);
  fails += Test_SUNLinSolSpace(LS, grid.iam);

  /* Print result */
  if (fails) {
    fprintf(fp, "process %6d: FAIL: SUNLinSol module failed %i tests \n \n", grid.iam, fails);
    fprintf(fp, "process %6d: \nA =\n", grid.iam);
    SUNMatrix_SLUNRloc_Print(A, fp);
    fprintf(fp, "process %6d: \nx (original) =\n", grid.iam);
    N_VPrintFile_Parallel(y, fp);
    fprintf(fp, "process %6d: \nb =\n", grid.iam);
    N_VPrintFile_Parallel(b, fp);
    fprintf(fp, "process %6d: \nx (computed)=\n", grid.iam);
    N_VPrintFile_Parallel(x, fp);
  }

  /* Free solver, matrix and vectors */
  if (grid.iam == 0) {
    SUNMatDestroy(D);
    N_VDestroy(gx); N_VDestroy(gy); N_VDestroy(gb);
  }
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(x); N_VDestroy(y); N_VDestroy(b);

  /* Free superlu-dist structures and exit */
  Destroy_CompRowLoc_Matrix_dist(Asuper);
  free(Asuper); Asuper = NULL;
  PStatFree(&stat);
  xScalePermstructFree(&scaleperm);
  xLUstructFree(&lu);

#if LOG_PROCESS_TO_FILE
  /* Close file pointer */
  fclose(fp);
#endif

  /* Check all procs for errors */
  MPI_Allreduce(&fails, &globfails, 1, MPI_INT, MPI_MAX, grid.comm);
  if (grid.iam == 0 && globfails > 0) {
    printf("FAIL: %d failures across %d processes\n", globfails, nprocs);
  } else if (grid.iam == 0) {
    printf("SUCCESS: SUNLinSol module passed all tests \n \n");
  }

  superlu_gridexit(&grid);
  SUNContext_Free(&sunctx);
  MPI_Finalize();
  return(globfails);
}

/* ----------------------------------------------------------------------
 * Implementation-specific 'check' routines
 * --------------------------------------------------------------------*/
int check_vector(N_Vector X, N_Vector Y, realtype tol)
{
  int failure = 0;
  sunindextype i, local_length, maxloc;
  realtype *Xdata, *Ydata, maxerr;

  Xdata = N_VGetArrayPointer(X);
  Ydata = N_VGetArrayPointer(Y);
  local_length = N_VGetLocalLength_Parallel(X);

  /* check vector data */
  for(i=0; i < local_length; i++)
    failure += SUNRCompareTol(Xdata[i], Ydata[i], tol);

  if (failure > ZERO) {
    maxerr = ZERO;
    maxloc = -1;
    for(i=0; i < local_length; i++) {
      if (std::abs(Xdata[i]-Ydata[i]) >  maxerr) {
        maxerr = std::abs(Xdata[i]-Ydata[i]);
        maxloc = i;
      }
    }
    printf("check err failure: maxerr = %g at loc %li (tol = %g)\n",
	   maxerr, (long int) maxloc, tol);
    return(1);
  }
  else
    return(0);
}

int csr_from_dense(SUNMatrix Ad, realtype droptol, realtype **matdata,
                   sunindextype **colind, sunindextype **rowptrs)
{
  sunindextype i, j, nnz;
  sunindextype M, N;

  if (droptol < ZERO)
    return -1;
  if (SUNMatGetID(Ad) != SUNMATRIX_DENSE)
    return -1;

  /* set size of new matrix */
  M = SUNDenseMatrix_Rows(Ad);
  N = SUNDenseMatrix_Columns(Ad);

  /* determine total number of nonzeros */
  nnz = 0;
  for (j=0; j<N; j++)
    for (i=0; i<M; i++)
      nnz += (std::abs(SM_ELEMENT_D(Ad,i,j)) > droptol);

  /* allocate */
  (*matdata) = (realtype*) malloc(nnz*sizeof(realtype));
  (*colind)  = (sunindextype*) malloc(nnz*sizeof(sunindextype));
  (*rowptrs) = (sunindextype*) malloc((M+1)*sizeof(sunindextype));

  /* copy nonzeros from Ad into As, based on CSR/CSC type */
  nnz = 0;
  for (i=0; i<M; i++) {
    (*rowptrs)[i] = nnz;
    for (j=0; j<N; j++) {
      if ( std::abs(SM_ELEMENT_D(Ad,i,j)) > droptol ) {
        (*colind)[nnz] = j;
        (*matdata)[nnz++] = SM_ELEMENT_D(Ad,i,j);
      }
    }
  }
  (*rowptrs)[M] = nnz;

  return 0;
}

void sync_device()
{
}
