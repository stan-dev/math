/*
 * ----------------------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * ----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * ----------------------------------------------------------------------------
 * This is the header file for the SUNMatrix that wraps the SuperLU-DIST
 * SLU_NR_loc type SuperMatrix.
 * ----------------------------------------------------------------------------
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <superlu_ddefs.h>

#include <nvector/nvector_serial.h>
#include <nvector/nvector_parallel.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_slunrloc.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include "test_sunmatrix.h"

/* Control whethere to print on rank 0 only or all ranks
 * when using TEST_STATUS* macros */
static int print_all_ranks;
static int print_timing;

/* helper functions */
int csr_from_dense(SUNMatrix Ad, realtype droptol, realtype **matdata,
                   sunindextype **colind, sunindextype **rowptrs);

/* ----------------------------------------------------------------------------
 * Main SUNMatrix Testing Routine
 * --------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int fails = 0;                            /* counter for test failures        */
  int globfails = 0;                        /* counter for test failures        */
  int nprocs, nprow, npcol, rank;           /* process grid size and rank       */
  MPI_Status mpistatus;                     /* MPI status                       */
  MPI_Comm comm;                            /* MPI communicator                 */
  gridinfo_t grid;                          /* SuperLU-DIST process grid        */

  sunindextype M, N;                        /* matrix size                      */
  realtype *matdata;                        /* pointer to matrix data           */
  sunindextype *rowptrs, *colind;           /* indexptrs and indexvals          */
  SUNMatrix D, A, I;                        /* global and local matrices        */
  SuperMatrix *Asuper, *Isuper;             /* SLU_NR_loc SuperMatrices         */
  NRformat_loc *Istore;                     /* SuperMatrix->Store               */
  N_Vector gx, gy, x, y;                    /* test vectors                     */
  realtype *xdata, *ydata;                  /* vector data                      */
  sunindextype M_local;                     /* num rows in local matrix         */
  sunindextype Mnprocs;                     /* M/nprocs                         */
  sunindextype Mrem;                        /* M%nprocs                         */
  sunindextype NNZ_local;                   /* num of nnz in local matrix       */
  sunindextype fst_row;                     /* global. index of 1st local row   */
  int square;                               /* is A a square matrix             */
  sunindextype i, j, k;                     /* just some iteration variables    */

  MPI_Init(&argc, &argv);

  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  if (SUNContext_Create(&comm, &sunctx)) {
    printf("ERROR: SUNContext_Create failed\n");
    return(-1);
  }

  if (argc < 5) {
    if (rank == 0)
      printf("ERROR: FIVE (5) arguments required: matrix size, nprow, npcol, print_all_ranks, print_timing\n");
    return(1);
  }

  /* should TEST_STATUS print on all ranks */
  print_all_ranks = atoi(argv[4]);
  SetPrintAllRanks(print_all_ranks);

  print_timing = atoi(argv[5]);
  SetTiming(print_timing);

  /* Validate the process grid size */
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

  /* Extract matrix size arguments */
  M = N = (sunindextype) atol(argv[1]);

  /* Setup the process grid */
  superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
  /* excess processes just exit */
  if (grid.iam >= nprow*npcol) {
    superlu_gridexit(&grid);
    MPI_Finalize();
    return(0);
  }
  /* now update the number of procs available */
  nprocs = nprow*npcol;

  /* Determine local matrix information */
  M_local = Mnprocs = M/nprocs;
  Mrem    = M - M_local*nprocs;
  fst_row = M_local*grid.iam;
  /* if there is remaining rows, the last process gets them */
  if ((M_local*nprocs) != M) {
    if (grid.iam == (nprocs-1)) /* last proc. gets all*/
      M_local = M - M_local*(nprocs-1);
  }

  /* Initialize matrices and vectors to null */
  D  = NULL;
  A  = NULL;
  I  = NULL;
  x  = NULL;
  y  = NULL;
  gx = NULL;
  gy = NULL;
  matdata = NULL;
  rowptrs = NULL;
  colind  = NULL;

  /* Setup global matrix */
  if (grid.iam == 0) {

    /* Create the matrix as dense first */
    D = SUNDenseMatrix(M, N, sunctx);

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
    xdata = N_VGetArrayPointer(gx);
    ydata = N_VGetArrayPointer(gy);

    /* fill x with random data, then find y = Dx using Matvec */
    for (i=0; i<N; i++)
      xdata[i] = (realtype) rand() / (realtype) RAND_MAX;
    if (SUNMatMatvec(D, gx, gy)) {
      TEST_STATUS(">>> FAIL: SUNMatMatvec failed for dense matrix\n", grid.iam);
      SUNMatDestroy(D); N_VDestroy(gx); N_VDestroy(gy);
      return(1);
    }

    fails = csr_from_dense(D, ZERO, &matdata, &colind, &rowptrs);
    if (fails != 0) {
      printf(">>> FAIL: csr_from_dense failure \n");
      return(1);
    }

    /* distribute matrix and vectors */
    for (i=1; i<nprocs; i++) {
      sunindextype *coltemp, *rowtemp;
      realtype *datatemp;
      sunindextype fstrow_temp, M_localtemp;

      M_localtemp = (i == (nprocs-1)) ? M_local + Mrem : Mnprocs;
      fstrow_temp = Mnprocs*i;

      /* send number of local NNZ */
      NNZ_local = rowptrs[fstrow_temp+M_localtemp] - rowptrs[fstrow_temp];
      MPI_Send(&NNZ_local, 1, MPI_SUNINDEXTYPE, i, i, grid.comm);

      /* send out rowptrs */
      rowtemp = &rowptrs[fstrow_temp];
      MPI_Send(rowtemp, M_localtemp+1, MPI_SUNINDEXTYPE, i, i, grid.comm);

      /* send corresponding column indices */
      coltemp = &colind[rowptrs[fstrow_temp]];
      MPI_Send(coltemp, NNZ_local, MPI_SUNINDEXTYPE, i, i, grid.comm);

      /* send corresponding data */
      datatemp = &matdata[rowptrs[fstrow_temp]];
      MPI_Send(datatemp, NNZ_local, MPI_SUNREALTYPE, i, i, grid.comm);

      /* send vector data */
      datatemp = &xdata[fstrow_temp];
      MPI_Send(datatemp, M_localtemp, MPI_SUNREALTYPE, i, i, grid.comm);

      datatemp = &ydata[fstrow_temp];
      MPI_Send(datatemp, M_localtemp, MPI_SUNREALTYPE, i, i, grid.comm);
    }

    NNZ_local = rowptrs[M_local] - rowptrs[0];

    /* Create local SuperLU-DIST SLU_NR_loc SuperMatrix */
    Asuper = NULL;
    Asuper = (SuperMatrix*) malloc(sizeof(SuperMatrix));

    /* Allocates Asuper->Store and sets structure members */
    dCreate_CompRowLoc_Matrix_dist(Asuper, M, N, NNZ_local, M_local, fst_row,
                                   matdata, colind, rowptrs, SLU_NR_loc, SLU_D, SLU_GE);

    A = SUNMatrix_SLUNRloc(Asuper, &grid, sunctx);
    if (A == NULL) {
      fails++;
      TEST_STATUS(">>> FAIL: Failed to create SUNMatrix_SLUNRloc\n", grid.iam);
      SUNMatDestroy(D); N_VDestroy(gx); N_VDestroy(gy);
      return(fails);
    }

    /* last entry in rowptrs should be the num of local nz */
    rowptrs[M_local] = NNZ_local;

    /* make the local NVectors */
    x = N_VMake_Parallel(grid.comm, M_local, N, xdata, sunctx);
    y = N_VMake_Parallel(grid.comm, M_local, N, ydata, sunctx);

  } else {

    sunindextype shift;

    /* recieve number of local nnz */
    MPI_Recv(&NNZ_local, 1, MPI_SUNINDEXTYPE, 0, grid.iam, grid.comm, &mpistatus);

    /* Allocate memory for matrix members */
    matdata = (realtype*) malloc(NNZ_local*sizeof(realtype));
    colind  = (sunindextype*) malloc(NNZ_local*sizeof(sunindextype));
    rowptrs = (sunindextype*) malloc((M_local+1)*sizeof(sunindextype));

    /* receive distributed matrix */
    MPI_Recv(rowptrs, M_local+1, MPI_SUNINDEXTYPE, 0, grid.iam, grid.comm, &mpistatus);
    MPI_Recv(colind, NNZ_local, MPI_SUNINDEXTYPE, 0, grid.iam, grid.comm, &mpistatus);
    MPI_Recv(matdata, NNZ_local, MPI_SUNREALTYPE, 0, grid.iam, grid.comm, &mpistatus);

    /* localize rowptrs */
    shift = rowptrs[0];
    for (i=0; i<M_local; i++)
      rowptrs[i] = rowptrs[i]-shift;
    rowptrs[M_local] = NNZ_local;

    /* Create local SuperLU-DIST SuperMatrix */
    Asuper = NULL;
    Asuper = (SuperMatrix*) malloc(sizeof(SuperMatrix));
    dCreate_CompRowLoc_Matrix_dist(Asuper, M, N, NNZ_local, M_local, fst_row,
                                   matdata, colind, rowptrs, SLU_NR_loc, SLU_D, SLU_GE);

    /* Create local SuperLU-DIST SUNMatrix */
    A = SUNMatrix_SLUNRloc(Asuper, &grid, sunctx);
    if (A == NULL) {
      fails++;
      TEST_STATUS(">>> FAIL: Failed to create SUNMatrix_SLUNRloc\n", grid.iam);
      Destroy_CompRowLoc_Matrix_dist(Asuper);
      return(fails);
    }

    /* make the local NVectors */
    x = N_VNew_Parallel(grid.comm, M_local, N, sunctx);
    y = N_VNew_Parallel(grid.comm, M_local, N, sunctx);
    xdata = N_VGetArrayPointer(x);
    ydata = N_VGetArrayPointer(y);

    /* recieve vectors */
    MPI_Recv(xdata, M_local, MPI_SUNREALTYPE, 0, grid.iam, grid.comm, &mpistatus);
    MPI_Recv(ydata, M_local, MPI_SUNREALTYPE, 0, grid.iam, grid.comm, &mpistatus);
  }

  /* Create local I matrix that is the same shape as A (which was
   * made with the diagonal) since SuperLU requires the the matrices
   * to be the same shape and have the same sparsity pattern. */
  square = is_square(A);
  if (square) {
    /* clone/copy A so we keep colind, rowptrs */
    I = SUNMatClone(A);
    if (I == NULL) {
      fails++;
      TEST_STATUS(">>> FAIL: Failed to alloc I\n", grid.iam);
    }
    fails += SUNMatCopy(A, I);

    Isuper = SUNMatrix_SLUNRloc_SuperMatrix(I);
    Istore = (NRformat_loc*) Isuper->Store;
    matdata = (realtype *) Istore->nzval;
    rowptrs = (sunindextype *) Istore->rowptr;
    colind  = (sunindextype *) Istore->colind;

    for(i=0; i<M_local; i++) {
      for (j=0; j<rowptrs[i+1]-rowptrs[i]; j++) {
        if (colind[rowptrs[i]+j] == (fst_row+i))
          matdata[rowptrs[i]+j] = ONE;
        else
          matdata[rowptrs[i]+j] = ZERO;
      }
    }
  }

  /* Commence SUNMatrix module tests */
  fails += Test_SUNMatGetID(A, SUNMATRIX_SLUNRLOC, grid.iam);
  fails += Test_SUNMatClone(A, grid.iam);
  fails += Test_SUNMatCopy(A, grid.iam);
  fails += Test_SUNMatZero(A, grid.iam);
  fails += Test_SUNMatScaleAdd(A, I, grid.iam);
  if (square)
    fails += Test_SUNMatScaleAddI(A, I, grid.iam);
  fails += Test_SUNMatMatvecSetup(A, grid.iam);
  fails += Test_SUNMatMatvec(A, x, y, grid.iam);
  fails += Test_SUNMatSpace(A, grid.iam);

  if (fails) {
    TEST_STATUS2("FAIL: SUNMatrix module SUNMatrix_SLUNRloc failed %i tests.\n\n",
                 fails, grid.iam);
  } else {
    TEST_STATUS("SUCCESS: SUNMatrix_SLUNRloc passed all tests.\n\n", grid.iam);
  }

  /* Destroy the SuperLU-DIST SuperMatrices */
  Destroy_CompRowLoc_Matrix_dist(Asuper);
  free(Asuper); Asuper = NULL;

  /* Destroy matricies and vectors */
  if (grid.iam == 0) {
    SUNMatDestroy(D);
    N_VDestroy(gx); N_VDestroy(gy);
  }
  SUNMatDestroy(A); SUNMatDestroy(I);
  N_VDestroy(x); N_VDestroy(y);

  /* Check all procs for errors */
  MPI_Allreduce(&fails, &globfails, 1, MPI_INT, MPI_MAX, grid.comm);

  superlu_gridexit(&grid);
  SUNContext_Free(&sunctx);
  MPI_Finalize();
  return(globfails);
}


/* ----------------------------------------------------------------------
 * Check matrix
 * --------------------------------------------------------------------*/
int check_matrix(SUNMatrix A, SUNMatrix B, realtype tol)
{
  int failure = 0;
  gridinfo_t *grid;
  SuperMatrix *Asuper, *Bsuper;
  NRformat_loc *Astore, *Bstore;
  realtype *Adata, *Bdata;
  sunindextype *Aindexptrs, *Bindexptrs;
  sunindextype *Aindexvals, *Bindexvals;
  sunindextype i, Arows, Brows, Annz, Bnnz;

  grid = SUNMatrix_SLUNRloc_ProcessGrid(A);

  /* get matrix pointers */
  Asuper     = SUNMatrix_SLUNRloc_SuperMatrix(A);
  Astore     = (NRformat_loc*) Asuper->Store;
  Adata      = (realtype*) Astore->nzval;
  Aindexptrs = Astore->rowptr;
  Aindexvals = Astore->colind;
  Arows      = Astore->m_loc;
  Annz       = Aindexptrs[Arows];

  Bsuper     = SUNMatrix_SLUNRloc_SuperMatrix(B);
  Bstore     = (NRformat_loc*) Bsuper->Store;
  Bdata      = (realtype*) Bstore->nzval;
  Bindexptrs = Bstore->rowptr;
  Bindexvals = Bstore->colind;
  Brows      = Bstore->m_loc;
  Bnnz       = Bindexptrs[Brows];

  /* matrices must have same sparsetype, shape and actual data lengths */
  if (SUNMatGetID(A) != SUNMatGetID(B)) {
    TEST_STATUS3(">>> ERROR: check_matrix: Different storage types (%d vs %d)\n",
                 SUNMatGetID(A), SUNMatGetID(B), grid->iam);
    return(1);
  }
  if (Asuper->nrow != Bsuper->nrow) {
    TEST_STATUS3(">>> ERROR: check_matrix: Different numbers of rows (%ld vs %ld)\n",
                 (long int) Asuper->nrow, (long int) Bsuper->nrow, grid->iam);
    return(1);
  }
  if (Asuper->ncol != Bsuper->ncol) {
    TEST_STATUS3(">>> ERROR: check_matrix: Different numbers of cols (%ld vs %ld)\n",
                 (long int) Asuper->nrow, (long int) Bsuper->nrow, grid->iam);
    return(1);
  }
  if (Arows != Brows) {
    TEST_STATUS3(">>> ERROR: check_matrix: Different numbers of local rows (%ld vs %ld)\n",
                 (long int) Arows, (long int) Brows, grid->iam);
    return(1);
  }
  if (Annz != Bnnz) {
    TEST_STATUS3(">>> ERROR: check_matrix: Different numbers of nonzeros (%ld vs %ld)\n",
                 (long int) Annz, (long int) Bnnz, grid->iam);
    return(1);
  }

  /* compare sparsity patterns */
  for (i=0; i<Arows; i++)
    failure += (Aindexptrs[i] != Bindexptrs[i]);
  if (failure > ZERO) {
    TEST_STATUS(">>> ERROR: check_matrix: Different indexptrs \n", grid->iam);
    return(1);
  }
  for (i=0; i<Annz; i++)
    failure += (Aindexvals[i] != Bindexvals[i]);
  if (failure > ZERO) {
    TEST_STATUS(">>> ERROR: check_matrix: Different indexvals \n", grid->iam);
    return(1);
  }

  /* compare global start index */
  if (Astore->fst_row != Bstore->fst_row) {
    TEST_STATUS(">>> ERRROR: check_matrix: Different globalidx \n", grid->iam);
    return(1);
  }

  /* compare matrix values */
  for(i=0; i<Annz; i++) {
    failure += SUNRCompareTol(Adata[i], Bdata[i], tol);
    if (SUNRCompareTol(Adata[i], Bdata[i], tol)) {
      TEST_STATUS3("Adata[%ld] != Bdata[%ld] \n", (long int) i, (long int) i,
                   grid->iam);
    }
  }
  if (failure > ZERO) {
    TEST_STATUS(">>> ERROR: check_matrix: Different entries \n", grid->iam);
    return(1);
  }

  return(0);
}

int check_matrix_entry(SUNMatrix A, realtype val, realtype tol)
{
  int failure = 0;
  SuperMatrix *Asuper;
  NRformat_loc *Astore;
  realtype *Adata;
  sunindextype i, nnz_loc;

  /* get matrix pointers */
  Asuper = SUNMatrix_SLUNRloc_SuperMatrix(A);
  Astore = (NRformat_loc*) Asuper->Store;
  Adata  = (realtype*) Astore->nzval;

  /* compare data */
  nnz_loc = Astore->nnz_loc;
  for(i=0; i < nnz_loc; i++) {
    if (SUNRCompareTol(Adata[i],val,tol) != 0) {
      printf("rhs=%g\n", std::abs(val)*tol);
      printf("  Adata[%ld] = %g != %g (err = %g)\n", (long int) i, Adata[i],
             val, std::abs(Adata[i]-val));
      failure++;
    }
  }

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

int check_vector(N_Vector x, N_Vector y, realtype tol)
{
  int failure = 0;
  realtype *xdata, *ydata;
  sunindextype xldata, yldata;
  sunindextype i;

  /* get vector data */
  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);

  /* check data lengths */
  xldata = N_VGetLocalLength_Parallel(x);
  yldata = N_VGetLocalLength_Parallel(y);

  if (xldata != yldata) {
    TEST_STATUS(">>> ERROR: check_vector: Different data array lengths \n", 0);
    return(1);
  }

  /* check vector data */
  for(i=0; i < xldata; i++)
    failure += SUNRCompareTol(xdata[i], ydata[i], tol);

  if (failure > ZERO)
    return(1);
  else
    return(0);
}

booleantype has_data(SUNMatrix A)
{
  SuperMatrix *Asuper;
  NRformat_loc *Astore;
  realtype *Adata;

  /* get matrix pointers */
  Asuper = SUNMatrix_SLUNRloc_SuperMatrix(A);
  Astore = (NRformat_loc*) Asuper->Store;
  Adata  = (realtype*) Astore->nzval;

  if (Adata == NULL)
    return SUNFALSE;
  else
    return SUNTRUE;
}

booleantype is_square(SUNMatrix A)
{
  SuperMatrix *Asuper = SUNMatrix_SLUNRloc_SuperMatrix(A);
  if (Asuper->nrow == Asuper->ncol)
    return SUNTRUE;
  else
    return SUNFALSE;
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

void sync_device(SUNMatrix A)
{
  /* not running on GPU, just return */
  return;
}
