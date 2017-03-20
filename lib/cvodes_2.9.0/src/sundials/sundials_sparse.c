/*
 * -----------------------------------------------------------------
 * $Revision: 4761 $
 * $Date: 2016-05-18 20:00:35 -0700 (Wed, 18 May 2016) $
 * -----------------------------------------------------------------
 * Programmers: Carol Woodward, Slaven Peles @ LLNL
 *              Daniel R. Reynolds @ SMU
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
 * This is the implementation file for operations on the SUNDIALS
 * sparse matrix structure.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_sparse.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/*
 * ==================================================================
 * Private function prototypes (functions working on SlsMat)
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Functions: SparseMatvecCSC
 * -----------------------------------------------------------------
 * This function computes the matrix-vector product, y=A*x, where A
 * is a CSC sparse matrix of dimension MxN, x is a realtype array of 
 * length N, and y is a realtype array of length M. Upon successful
 * completion, the return value is zero; otherwise 1 is returned.
 * -----------------------------------------------------------------
 */

static int SparseMatvecCSC(const SlsMat A, const realtype *x, realtype *y);

/*
 * -----------------------------------------------------------------
 * Functions: SparseMatvecCSR
 * -----------------------------------------------------------------
 * This function computes the matrix-vector product, y=A*x, where A
 * is a CSR sparse matrix of dimension MxN, x is a realtype array of 
 * length N, and y is a realtype array of length M. Upon successful
 * completion, the return value is zero; otherwise 1 is returned.
 * -----------------------------------------------------------------
 */

static int SparseMatvecCSR(const SlsMat A, const realtype *x, realtype *y);



/*
 * ==================================================================
 * Implementation of sparse matrix methods (functions on SlsMat)
 * ==================================================================
 */

/* 
 * Default Constructor
 * 
 * Creates a new (empty) sparse matrix of a desired size and nonzero density.
 * Returns NULL if a memory allocation error occurred.
 * 
 */
SlsMat SparseNewMat(int M, int N, int NNZ, int sparsetype)
{
  SlsMat A;

  if ( (M <= 0) || (N <= 0) ) return(NULL);

  A = NULL;
  A = (SlsMat) malloc(sizeof(struct _SlsMat));
  if (A==NULL) return (NULL);
  
  A->sparsetype = sparsetype;
  
  switch(A->sparsetype){
    case CSC_MAT:
      A->NP = N;
      A->rowvals = &(A->indexvals);
      A->colptrs = &(A->indexptrs);
      /* CSR indices */
      A->colvals = NULL;
      A->rowptrs = NULL;
      break;
    case CSR_MAT:
      A->NP = M;
      A->colvals = &(A->indexvals);
      A->rowptrs = &(A->indexptrs);
      /* CSC indices */
      A->rowvals = NULL;
      A->colptrs = NULL;
      break;
    default:
      free(A); 
      A = NULL;
      return(NULL);
  }

  A->data = (realtype *) malloc(NNZ * sizeof(realtype));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }
  
  A->indexvals = (int *) malloc(NNZ * sizeof(int));
  if (A->indexvals == NULL) {
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }
  A->indexptrs = (int *) malloc((A->NP + 1) * sizeof(int));
  if (A->indexptrs == NULL) {
    free(A->indexvals);
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }

  A->M = M;
  A->N = N;
  A->NNZ = NNZ;
  /* A->colptrs[N] = NNZ; */
  A->indexptrs[A->NP] = 0;

  return(A);
}

/** 
 * Constructor
 * 
 * Creates a new sparse matrix out of an existing dense or band matrix.  
 * Returns NULL if a memory allocation error occurred.
 * 
 */
SlsMat SparseFromDenseMat(const DlsMat Ad, int sparsetype)
{
  int i, j, nnz;
  int M, N;
  realtype dtmp;
  SlsMat As = NULL;
  
  switch(sparsetype) {
    case CSC_MAT:
      /* CSC is transpose of CSR */
      M = Ad->N;
      N = Ad->M;
      break;
    case CSR_MAT:
      M = Ad->M;
      N = Ad->N;
      break;
    default:
      /* Sparse matrix type not recognized */
      return NULL;
  }

  /* proceed according to A's type (dense/band) */
  if (Ad->type == SUNDIALS_DENSE) {

    /* determine total number of nonzeros */
    nnz = 0;
    for (j=0; j<Ad->N; j++)
      for (i=0; i<Ad->M; i++)
        nnz += (DENSE_ELEM(Ad,i,j) != 0.0);

    /* allocate sparse matrix */
    As = SparseNewMat(Ad->M, Ad->N, nnz, sparsetype);
    if (As == NULL)  return NULL;

    /* copy nonzeros from A into As */
    nnz = 0;
    for (i=0; i<M; i++) {
      (As->indexptrs)[i] = nnz;
      for (j=0; j<N; j++) {
        /* CSR = row major looping; CSC = column major looping */
        dtmp = (sparsetype == CSR_MAT) ? DENSE_ELEM(Ad,i,j) : DENSE_ELEM(Ad,j,i);
        if ( dtmp != ZERO ) { 
          (As->indexvals)[nnz] = j;
          As->data[nnz++] = dtmp;
        }
      }
    }
    (As->indexptrs)[M] = nnz;

  } else { /* SUNDIALS_BAND */

    /* determine total number of nonzeros */
    nnz = 0;
    for (j=0; j<Ad->N; j++)
      for (i=j-(Ad->mu); i<j+(Ad->ml); i++)
        nnz += (BAND_ELEM(Ad,i,j) != 0.0);

    /* allocate sparse matrix */
    As = SparseNewMat(Ad->M, Ad->N, nnz, sparsetype);
    if (As == NULL)  return NULL;

    /* copy nonzeros from A into As */
    nnz = 0;
    for (i=0; i<M; i++) {
      (As->indexptrs)[i] = nnz;
      for (j=i-(Ad->mu); j<i+(Ad->ml); j++) {
        /* CSR = row major looping; CSC = column major looping */
        dtmp = (sparsetype == CSR_MAT) ? BAND_ELEM(Ad,i,j) : BAND_ELEM(Ad,j,i);
        if ( dtmp != 0.0 ) { 
          (As->indexvals)[nnz] = j;
          As->data[nnz++] = dtmp;
        }
      }
    }
    (As->indexptrs)[M] = nnz;

  }

  return(As);
}


/**
 * 
 * Destructor
 * 
 * Frees memory and deletes the structure for an existing sparse matrix.
 * 
 */
int SparseDestroyMat(SlsMat A)
{
  if (A->data) {
    free(A->data);  
    A->data = NULL;
  }
  if (A->indexvals) {
    free(A->indexvals);
    A->indexvals = NULL;
    A->rowvals   = NULL;
    A->colvals   = NULL;
  }
  if (A->indexptrs) {
    free(A->indexptrs);
    A->indexptrs = NULL;
    A->colptrs   = NULL;
    A->rowptrs   = NULL;
  }
  free(A); 
  A = NULL;
  
  return 0;
}


/** 
 * Sets all sparse matrix entries to zero.
 */
int SparseSetMatToZero(SlsMat A)
{
  int i;

  for (i=0; i<A->NNZ; i++) {
    A->data[i] = ZERO;
    A->indexvals[i] = 0;
  }

  for (i=0; i<A->NP; i++) {
    A->indexptrs[i] = 0;
  }
  /* A->colptrs[A->N] = A->NNZ; */
  A->indexptrs[A->NP] = 0;
  
  return 0;
}


/** 
 * Copies the sparse matrix A into sparse matrix B.  
 * 
 * It is assumed that A and B have the same dimensions, but we account 
 * for the situation in which B has fewer nonzeros than A.
 *  
 */
int SparseCopyMat(const SlsMat A, SlsMat B)
{
  int i;
  int A_nz = A->indexptrs[A->NP];
  
  if(A->M != B->M || A->N != B->N) {
    /* fprintf(stderr, "Error: Copying sparse matrices of different size!\n"); */
    return (-1);
  }
    
  
  /* ensure B is of the same type as A */
  B->sparsetype = A->sparsetype;

  /* ensure that B is allocated with at least as 
     much memory as we have nonzeros in A */
  if (B->NNZ < A_nz) {
    B->indexvals = realloc(B->indexvals, A_nz*sizeof(int));
    B->data = realloc(B->data, A_nz*sizeof(realtype));
    B->NNZ = A_nz;
  }

  /* zero out B so that copy works correctly */
  SparseSetMatToZero(B);

  /* copy the data and row indices over */
  for (i=0; i<A_nz; i++){
    B->data[i] = A->data[i];
    B->indexvals[i] = A->indexvals[i];
  }

  /* copy the column pointers over */
  for (i=0; i<A->NP; i++) {
    B->indexptrs[i] = A->indexptrs[i];
  }
  B->indexptrs[A->NP] = A_nz;
  
  return 0;
}


/** 
 * Scales a sparse matrix A by the coefficient b.
 */
int SparseScaleMat(realtype b, SlsMat A)
{
  int i;

  for (i=0; i<A->indexptrs[A->NP]; i++){
    A->data[i] = b * (A->data[i]);
  }
  return 0;
}




/** 
 * Adds 1 to every diagonal entry of A.  
 * 
 * Works for general [rectangular] matrices and handles potentially increased 
 * size if A does not currently contain a value on the diagonal.
 * 
 * The function was developed originally for CSC matrices. To make it work for 
 * CSR, one simply need to transpose it, i.e. transpose M and N in the 
 * implementation.  
 * 
 */
int SparseAddIdentityMat(SlsMat A)
{
  int j, i, p, nz, newmat, found;
  int *w, *Ap, *Ai, *Cp, *Ci;
  realtype *x, *Ax, *Cx;
  SlsMat C;
  int M;
  int N;

  /* determine if A already contains values on the diagonal (hence 
     memory allocation necessary)*/
  newmat=0;
  for (j=0; j < SUNMIN(A->N,A->M); j++) {
    /* scan column (row if CSR) of A, searching for diagonal value */
    found = 0;
    for (i=A->indexptrs[j]; i<A->indexptrs[j+1]; i++) {
      if (A->indexvals[i] == j) {
        found = 1;
        break;
      }
    }
    /* if no diagonal found, signal new matrix */
    if (!found) {
      newmat=1;
      break;
    }
  }

  /* perform operation */

  /*   case 1: A already contains a diagonal */
  if (!newmat) {

    /* iterate through columns, adding 1.0 to diagonal */
    for (j=0; j < SUNMIN(A->N,A->M); j++)
      for (i=A->indexptrs[j]; i<A->indexptrs[j+1]; i++)
        if (A->indexvals[i] == j) 
          A->data[i] += ONE;

  /*   case 2: A does not already contain a diagonal */
  } else {
    
    if (A->sparsetype == CSC_MAT) {
      M = A->M;
      N = A->N;
    }
    else if (A->sparsetype == CSR_MAT) {
      M = A->N;
      N = A->M;
    }
    else
      return (-1);
  
    /* create work arrays for row indices and nonzero column values */
    w = (int *) malloc(A->M * sizeof(int));
    x = (realtype *) malloc(A->M * sizeof(realtype));

    /* create new matrix for sum (overestimate nnz as sum of each) */
    C = SparseNewMat(A->M, A->N, (A->indexptrs)[A->NP] + SUNMIN(A->M, A->N), A->sparsetype);

    /* access data from CSR structures (return if failure) */
    Cp = Ci = Ap = Ai = NULL;
    Cx = Ax = NULL;
    if (C->indexptrs)  Cp = C->indexptrs;
    else  return (-1);
    if (C->indexvals)  Ci = C->indexvals;
    else  return (-1);
    if (C->data)       Cx = C->data;
    else  return (-1);
    if (A->indexptrs)  Ap = A->indexptrs;
    else  return (-1);
    if (A->indexvals)  Ai = A->indexvals;
    else  return (-1);
    if (A->data)       Ax = A->data;
    else  return (-1);

    /* initialize total nonzero count */
    nz = 0;

    /* iterate through columns (rows for CSR) */
    for (j=0; j<N; j++) {

      /* set current column (row) pointer to current # nonzeros */
      Cp[j] = nz;

      /* clear out temporary arrays for this column (row) */
      for (i=0; i<M; i++) {
        w[i] = 0;
        x[i] = 0.0;
      }

      /* iterate down column (along row) of A, collecting nonzeros */
      for (p=Ap[j]; p<Ap[j+1]; p++) {
        w[Ai[p]] += 1;       /* indicate that row is filled */
        x[Ai[p]] = Ax[p];    /* collect value */
      }

      /* add identity to this column (row) */
      if (j < M) {
        w[j] += 1;     /* indicate that row is filled */
        x[j] += ONE;   /* update value */
      }

      /* fill entries of C with this column's (row's) data */
      for (i=0; i<M; i++) {
        if ( w[i] > 0 ) { 
          Ci[nz] = i;  
          Cx[nz++] = x[i];
        }
      }
    }

    /* indicate end of data */
    Cp[N] = nz;

    /* update A's structure with C's values; nullify C's pointers */
    A->NNZ = C->NNZ;

    if (A->data)
      free(A->data);  
    A->data = C->data;
    C->data = NULL;

    if (A->indexvals)
      free(A->indexvals);
    A->indexvals = C->indexvals;
    C->indexvals = NULL;

    if (A->indexptrs)
      free(A->indexptrs);
    A->indexptrs = C->indexptrs;
    C->indexptrs = NULL;

    /* clean up */
    SparseDestroyMat(C); 
    free(w);
    free(x);

    /* reallocate the new matrix to remove extra space */
    SparseReallocMat(A);
  }
  return 0;
}


/** 
 * Add two sparse matrices: A = A+B.  
 * 
 * Handles potentially increased size if matrices have different sparsity patterns.  
 * Returns 0 if successful, and 1 if unsuccessful (in which case A is left unchanged).
 * 
 * The function was developed originally for CSC matrices. To make it work for 
 * CSR, one simply need to transpose it, i.e. transpose M and N in the 
 * implementation.  
 * 
 */
int SparseAddMat(SlsMat A, const SlsMat B)
{
  int j, i, p, nz, newmat;
  int *w, *Ap, *Ai, *Bp, *Bi, *Cp, *Ci;
  realtype *x, *Ax, *Bx, *Cx;
  SlsMat C;
  int M;
  int N;

  /* ensure that matrix dimensions agree */
  if ((A->M != B->M) || (A->N != B->N)) {
    /* fprintf(stderr, "Error: Adding sparse matrices of different size!\n"); */
    return(-1);
  }
  
  /* if A is CSR matrix, transpose M and N */
  if (A->sparsetype == CSC_MAT) {
    M = A->M;
    N = A->N;
  }
  else if (A->sparsetype == CSR_MAT) {
    M = A->N;
    N = A->M;
  }
  else
    return(-1);
  
  /* create work arrays for row indices and nonzero column values */
  w = (int *) malloc(M * sizeof(int));
  x = (realtype *) malloc(M * sizeof(realtype));

  /* determine if A already contains the sparsity pattern of B */
  newmat=0;
  for (j=0; j<N; j++) {

    /* clear work array */
    for (i=0; i<M; i++)  w[i] = 0;

    /* scan column of A, incrementing w by one */
    for (i=A->indexptrs[j]; i<A->indexptrs[j+1]; i++)
      w[A->indexvals[i]] += 1;

    /* scan column of B, decrementing w by one */
    for (i=B->indexptrs[j]; i<B->indexptrs[j+1]; i++)
      w[B->indexvals[i]] -= 1;

    /* if any entry of w is negative, A doesn't contain B's sparsity */
    for (i=0; i<M; i++)
      if (w[i] < 0) {
        newmat = 1;
        break;
      }
    if (newmat) break;

  }

  /* perform operation */

  /*   case 1: A already contains sparsity pattern of B */
  if (!newmat) {

    /* iterate through columns, adding matrices */
    for (j=0; j<N; j++) {

      /* clear work array */
      for (i=0; i<M; i++)
        x[i] = ZERO;

      /* scan column of B, updating work array */
      for (i = B->indexptrs[j]; i < B->indexptrs[j+1]; i++)
        x[B->indexvals[i]] = B->data[i];

      /* scan column of A, updating entries appropriately array */
      for (i = A->indexptrs[j]; i < A->indexptrs[j+1]; i++)
        A->data[i] += x[A->indexvals[i]];

    }

  /*   case 2: A does not already contain B's sparsity */
  } else {

    /* create new matrix for sum (overestimate nnz as sum of each) */
    C = SparseNewMat(M, N, (A->indexptrs[N])+(B->indexptrs[N]), A->sparsetype);

    /* access data from CSR structures (return if failure) */
    Cp = Ci = Ap = Ai = Bp = Bi = NULL;
    Cx = Ax = Bx = NULL;
    if (C->indexptrs)  Cp = C->indexptrs;
    else  return(-1);
    if (C->indexvals)  Ci = C->indexvals;
    else  return(-1);
    if (C->data)       Cx = C->data;
    else  return(-1);
    if (A->indexptrs)  Ap = (A->indexptrs);
    else  return(-1);
    if (A->indexvals)  Ai = (A->indexvals);
    else  return(-1);
    if (A->data)       Ax = A->data;
    else  return(-1);
    if (B->indexptrs)  Bp = B->indexptrs;
    else  return(-1);
    if (B->indexvals)  Bi = B->indexvals;
    else  return(-1);
    if (B->data)       Bx = B->data;
    else  return(-1);

    /* initialize total nonzero count */
    nz = 0;

    /* iterate through columns */
    for (j=0; j<N; j++) {

      /* set current column pointer to current # nonzeros */
      Cp[j] = nz;

      /* clear out temporary arrays for this column */
      for (i=0; i<M; i++) {
        w[i] = 0;
        x[i] = 0.0;
      }

      /* iterate down column of A, collecting nonzeros */
      for (p=Ap[j]; p<Ap[j+1]; p++) {
        w[Ai[p]] += 1;       /* indicate that row is filled */
        x[Ai[p]] = Ax[p];    /* collect value */
      }

      /* iterate down column of B, collecting nonzeros */
      for (p=Bp[j]; p<Bp[j+1]; p++) {
        w[Bi[p]] += 1;       /* indicate that row is filled */
        x[Bi[p]] += Bx[p];   /* collect value */
      }

      /* fill entries of C with this column's data */
      for (i=0; i<M; i++) {
        if ( w[i] > 0 ) { 
          Ci[nz] = i;  
          Cx[nz++] = x[i];
        }
      }
    }

    /* indicate end of data */
    Cp[N] = nz;

    /* update A's structure with C's values; nullify C's pointers */
    A->NNZ = C->NNZ;

    free(A->data);  
    A->data = C->data;
    C->data = NULL;

    free(A->indexvals);
    A->indexvals = C->indexvals;
    C->indexvals = NULL;

    free(A->indexptrs);
    A->indexptrs = C->indexptrs;
    C->indexptrs = NULL;

    /* clean up */
    SparseDestroyMat(C); 

    /* reallocate the new matrix to remove extra space */
    SparseReallocMat(A);

  }

  /* clean up */
  free(w);
  free(x);

  /* return success */
  return(0);
}


/** 
 * Resizes the memory allocated for a given sparse matrix, shortening 
 * it down to the number of actual nonzero entries.
 */
int SparseReallocMat(SlsMat A)
{
  int nzmax; 

  nzmax = A->indexptrs[A->NP];
  A->indexvals = realloc(A->indexvals, nzmax*sizeof(int));
  A->data = realloc(A->data, nzmax*sizeof(realtype));
  A->NNZ = nzmax;
  
  return 0;
}


/** 
 * Computes y=A*x, where A is a sparse matrix of dimension MxN, x is a 
 * realtype array of length N, and y is a realtype array of length M. 
 * 
 * Returns 0 if successful, -1 if unsuccessful (failed memory access).
 */
int SparseMatvec(const SlsMat A, const realtype *x, realtype *y)
{
  if(A->sparsetype == CSC_MAT)
    return SparseMatvecCSC(A, x, y);
  else if (A->sparsetype == CSR_MAT)
    return SparseMatvecCSR(A, x, y);
  else
    return(-1);
}


/** 
 * Prints the nonzero entries of a sparse matrix to screen.
 */
void SparsePrintMat(const SlsMat A, FILE* outfile)
{
  int i,j, NNZ;
  char *matrixtype;
  char *indexname;

  NNZ = A->NNZ;

  switch(A->sparsetype) {
    case CSC_MAT:
      indexname = (char*) "col";
      matrixtype = (char*) "CSC";
      break;
    case CSR_MAT:
      indexname = (char*) "row";
      matrixtype = (char*) "CSR";
      break;
    default:
      /* Sparse matrix type not recognized */
      return;
  }


  fprintf(outfile, "\n");
  
  fprintf(outfile, "%d by %d %s matrix, NNZ: %d \n", A->M, A->N, matrixtype, NNZ);
  for (j=0; j < A->NP; j++) {
    fprintf(outfile, "%s %d : locations %d to %d\n", indexname, j, (A->indexptrs)[j], (A->indexptrs)[j+1]-1);
    fprintf(outfile, "  ");
    for (i = (A->indexptrs)[j]; i < (A->indexptrs)[j+1]; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(outfile, "%d: %Lg   ", A->indexvals[i], A->data[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(outfile, "%d: %g   ", A->indexvals[i], A->data[i]);
#else
      fprintf(outfile, "%d: %g   ", A->indexvals[i], A->data[i]);
#endif
    }
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "\n");
    
}



/*
 * ==================================================================
 * Private function definitions
 * ==================================================================
 */



/** 
 * Computes y=A*x, where A is a CSC matrix of dimension MxN, x is a 
 * realtype array of length N, and y is a realtype array of length M. 
 * 
 * Returns 0 if successful, -1 if unsuccessful (failed memory access).
 */
int SparseMatvecCSC(const SlsMat A, const realtype *x, realtype *y)
{
  int j, i;
  int *Ap, *Ai;
  realtype *Ax;

  /* access data from CSR structure (return if failure) */
  if (*A->colptrs)  Ap = A->indexptrs;
  else  return(-1);
  if (*A->rowvals)  Ai = A->indexvals;
  else  return(-1);
  if (A->data)      Ax = A->data;
  else  return(-1);

  /* ensure that vectors are non-NULL */
  if ((x == NULL) || (y == NULL))
    return(-1);

  /* initialize result */
  for (i=0; i<A->M; i++)
    y[i] = 0.0;

  /* iterate through matrix columns */
  for (j=0; j<A->N; j++) {

    /* iterate down column of A, performing product */
    for (i=Ap[j]; i<Ap[j+1]; i++)
      y[Ai[i]] += Ax[i]*x[j];

  }

  /* return success */
  return(0);
}


/** 
 * Computes y=A*x, where A is a CSR matrix of dimension MxN, x is a 
 * realtype array of length N, and y is a realtype array of length M. 
 * 
 * Returns 0 if successful, -1 if unsuccessful (failed memory access).
 */
int SparseMatvecCSR(const SlsMat A, const realtype *x, realtype *y)
{
  int j, i;
  int *Ap, *Aj;
  realtype *Ax;

  /* access data from CSR structure (return if failure) */
  if (*A->rowptrs)  Ap = A->indexptrs;
  else  return(-1);
  if (*A->colvals)  Aj = A->indexvals;
  else  return(-1);
  if (A->data)      Ax = A->data;
  else  return(-1);

  /* ensure that vectors are non-NULL */
  if ((x == NULL) || (y == NULL))
    return(-1);

  /* initialize result */
  for (i=0; i<A->M; i++)
    y[i] = 0.0;

  /* iterate through matrix rows */
  for (i=0; i<A->M; ++i) {

    /* iterate along row of A, performing product */
    for (j=Ap[i]; j<Ap[i+1]; ++j)
      y[i] += Ax[j]*x[Aj[j]];

  }

  /* return success */
  return(0);
}


