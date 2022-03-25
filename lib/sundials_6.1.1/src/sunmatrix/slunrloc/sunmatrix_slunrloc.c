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
 * This is the implementation file for the SuperLU SuperMatrix SLU_NR_loc
 * format compatibile SUNMatrix.
 * ----------------------------------------------------------------------------
 */

#include <mpi.h>
#include <stdarg.h>
#include <stdlib.h>
#include <superlu_ddefs.h>
#include <sunmatrix/sunmatrix_slunrloc.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_mpi_types.h>

/*
 * ----------------------------------------------------------------------------
 *  Macros for accessing the SUNMatrix_SLUNRloc content structure members
 *  and useful values stored in nested SuperLU structures.
 * ----------------------------------------------------------------------------
 */

#define SM_CONTENT_SLUNRLOC(A)       ( (SUNMatrixContent_SLUNRloc)(A->content) )

#define SM_SUPERMATRIX_SLUNRLOC(A)   ( SM_CONTENT_SLUNRLOC(A)->A_super )

#define SM_GRID_SLUNRLOC(A)          ( SM_CONTENT_SLUNRLOC(A)->grid )

#define SM_COMMPATTERN_SLUNRLOC(A)   ( SM_CONTENT_SLUNRLOC(A)->gsmv_comm )

#define SM_ROWTOPROC_SLUNRLOC(A)     ( SM_CONTENT_SLUNRLOC(A)->row_to_proc )

#define SM_OWNDATA_SLUNRLOC(A)       ( SM_CONTENT_SLUNRLOC(A)->own_data )

#define SM_COLSORTED_SLUNRLOC(A)     ( SM_CONTENT_SLUNRLOC(A)->ACS_super )

#define SM_SUPERSTORE_SLUNRLOC(A)    ( (NRformat_loc*)(SM_SUPERMATRIX_SLUNRLOC(A)->Store) )

#define SM_GLOBALROWS_SLUNRLOC(A)    ( SM_SUPERMATRIX_SLUNRLOC(A)->nrow )

#define SM_GLOBALCOLS_SLUNRLOC(A)    ( SM_SUPERMATRIX_SLUNRLOC(A)->ncol )

#define SM_LOCALROWS_SLUNRLOC(A)     ( SM_SUPERSTORE_SLUNRLOC(A)->m_loc )

#define SM_LOCALNNZ_SLUNRLOC(A)      ( SM_SUPERSTORE_SLUNRLOC(A)->nnz_loc )

#define SM_FSTROW_SLUNRLOC(A)        ( SM_SUPERSTORE_SLUNRLOC(A)->fst_row )

#define SM_DATA_SLUNRLOC(A)          ( (realtype*)SM_SUPERSTORE_SLUNRLOC(A)->nzval )

#define SM_COLIND_SLUNRLOC(A)        ( SM_SUPERSTORE_SLUNRLOC(A)->colind )

#define SM_ROWPTRS_SLUNRLOC(A)       ( SM_SUPERSTORE_SLUNRLOC(A)->rowptr )

/* constants */
#define ZERO RCONST(0.0)

/* Private function prototypes */
static booleantype SMCompatible_SLUNRloc(SUNMatrix A, SUNMatrix B);

/*
 * ----------------------------------------------------------------------------
 *  Exported functions
 * ----------------------------------------------------------------------------
 */

SUNMatrix SUNMatrix_SLUNRloc(SuperMatrix *A_super, gridinfo_t *grid)
{
  SUNMatrix A;
  SUNMatrix_Ops ops;
  SUNMatrixContent_SLUNRloc content;

  /* Check for valid intputs */
  if (A_super == NULL || grid == NULL)
    return(NULL);

  if (A_super->Stype != SLU_NR_loc ||
      A_super->Dtype != SLU_D ||
      A_super->Mtype != SLU_GE)
    return(NULL);

  /* Create the matrix */
  A = NULL;
  A = (SUNMatrix) malloc(sizeof(*A));
  if (A == NULL) return(NULL);

  /* Create the matrix operations structure */
  ops = NULL;
  ops = (SUNMatrix_Ops) malloc(sizeof(struct _generic_SUNMatrix_Ops));
  if (ops == NULL) { free(A); return(NULL); }

  /* Attach operations */
  ops->getid       = SUNMatGetID_SLUNRloc;
  ops->clone       = SUNMatClone_SLUNRloc;
  ops->destroy     = SUNMatDestroy_SLUNRloc;
  ops->zero        = SUNMatZero_SLUNRloc;
  ops->copy        = SUNMatCopy_SLUNRloc;
  ops->scaleadd    = SUNMatScaleAdd_SLUNRloc;
  ops->scaleaddi   = SUNMatScaleAddI_SLUNRloc;
  ops->matvecsetup = SUNMatMatvecSetup_SLUNRloc;
  ops->matvec      = SUNMatMatvec_SLUNRloc;
  ops->space       = SUNMatSpace_SLUNRloc;

  /* Create content */
  content = NULL;
  content = (SUNMatrixContent_SLUNRloc)
            malloc(sizeof(struct _SUNMatrixContent_SLUNRloc));
  if (content == NULL) { free(ops); free(A); return(NULL); }

  /* Attach content to SuperMatrix */
  content->A_super     = A_super;
  content->own_data    = SUNFALSE;
  content->grid        = grid;
  content->row_to_proc = NULL;
  content->gsmv_comm   = NULL;
  content->ACS_super   = NULL;

  /* Attach content and ops to SUNMatrix */
  A->content = content;
  A->ops     = ops;

  return(A);
}

void SUNMatrix_SLUNRloc_Print(SUNMatrix A, FILE *fp)
{
  STAN_SUNDIALS_FPRINTF(fp, "====== START SUNMatrix_SLUNRloc_Print %p  ======\n", (void*) A);
  STAN_SUNDIALS_FPRINTF(fp, "A->content->A_super = %p\n", (void*) SM_SUPERMATRIX_SLUNRLOC(A));

  /* Call SuperLU_DIST print routine */
  file_dPrint_CompRowLoc_Matrix_dist(fp, SM_SUPERMATRIX_SLUNRLOC(A));

  STAN_SUNDIALS_FPRINTF(fp, "======= END SUNMatrix_SLUNRloc_Print %p  =======\n", (void*) A);
}


/*
 * ----------------------------------------------------------------------------
 *  Implementation of accessor functions
 * ----------------------------------------------------------------------------
 */

SuperMatrix* SUNMatrix_SLUNRloc_SuperMatrix(SUNMatrix A)
{
  return(SM_SUPERMATRIX_SLUNRLOC(A));
}

gridinfo_t* SUNMatrix_SLUNRloc_ProcessGrid(SUNMatrix A)
{
  return(SM_GRID_SLUNRLOC(A));
}

booleantype SUNMatrix_SLUNRloc_OwnData(SUNMatrix A)
{
  return(SM_OWNDATA_SLUNRLOC(A));
}

/*
 * ----------------------------------------------------------------------------
 *  Implementation of matrix operations
 * ----------------------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_SLUNRloc(SUNMatrix A)
{
  return(SUNMATRIX_SLUNRLOC);
}

SUNMatrix SUNMatClone_SLUNRloc(SUNMatrix A)
{
  SUNMatrix B;
  SuperMatrix *B_super;

  /* allocate new SuperMatrix */
  B_super = NULL;
  B_super = malloc(sizeof(SuperMatrix));
  if (B_super == NULL) return(NULL);

  /* call the SuperLU-DIST clone function */
  dClone_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A), B_super);

  /* create the new SUNMatrix */
  B = NULL;
  B = SUNMatrix_SLUNRloc(B_super, SM_GRID_SLUNRLOC(A));
  if (B == NULL) {
    /* call SuperLU-DIST destroy function */
    Destroy_CompRowLoc_Matrix_dist(B_super);
    free(B_super);
    return(NULL);
  }

  /* Allocated the SuperMatrix ourselves, so SUNMatrix now owns the data */
  SM_OWNDATA_SLUNRLOC(B) = SUNTRUE;

  return(B);
}

void SUNMatDestroy_SLUNRloc(SUNMatrix A)
{
  if (SM_OWNDATA_SLUNRLOC(A)) {
    /* call SuperLU-DIST destroy function */
    Destroy_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A));
    free(SM_SUPERMATRIX_SLUNRLOC(A));
  }
  if (SM_COLSORTED_SLUNRLOC(A)) {
    /* if CS exists, then the Matvec has been initialized and we must finalize
       it by calling the SuperLU DIST pdgsmv_finalize routine to free up memory
       allocated */

    pdgsmv_finalize(SM_COMMPATTERN_SLUNRLOC(A));
    Destroy_CompRowLoc_Matrix_dist(SM_COLSORTED_SLUNRLOC(A));
    free(SM_COLSORTED_SLUNRLOC(A));
    SM_COLSORTED_SLUNRLOC(A) = NULL;
  }
  if (SM_ROWTOPROC_SLUNRLOC(A)) {
    free(SM_ROWTOPROC_SLUNRLOC(A));
    SM_ROWTOPROC_SLUNRLOC(A) = NULL;
  }
  if (SM_COMMPATTERN_SLUNRLOC(A)) {
    free(SM_COMMPATTERN_SLUNRLOC(A));
    SM_COMMPATTERN_SLUNRLOC(A) = NULL;
  }
  free(A->content); A->content = NULL;
  free(A->ops); A->ops = NULL;
  free(A); A = NULL;
  return;
}

int SUNMatZero_SLUNRloc(SUNMatrix A)
{
  /* call SuperLU-DIST clone function */
  dZero_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A));
  return(SUNMAT_SUCCESS);
}

int SUNMatCopy_SLUNRloc(SUNMatrix A, SUNMatrix B)
{
  /* call SuperLU-DIST copy function */
  dCopy_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A),
                               SM_SUPERMATRIX_SLUNRLOC(B));
  return(SUNMAT_SUCCESS);
}

int SUNMatScaleAdd_SLUNRloc(realtype c, SUNMatrix A, SUNMatrix B)
{
  /* check that B can be added into A */
  if (!SMCompatible_SLUNRloc(A, B)) return(SUNMAT_ILL_INPUT);

  /* call SuperLU-DIST ScaleAdd function */
  dScaleAdd_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A),
                                   SM_SUPERMATRIX_SLUNRLOC(B), c);
  return(SUNMAT_SUCCESS);
}

int SUNMatScaleAddI_SLUNRloc(realtype c, SUNMatrix A)
{
  /* call SuperLU-DIST ScaleAddI function */
  dScaleAddId_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A), c);
  return(SUNMAT_SUCCESS);
}

int SUNMatMatvec_SLUNRloc(SUNMatrix A, N_Vector x, N_Vector y)
{
  SuperMatrix *ACS;
  realtype *xdata, *ydata;

  /* Extract the column-sorted A */
  ACS = SM_COLSORTED_SLUNRLOC(A);

  /* The column-sorted matrix ACS and the communication pattern must be
     established by calling SUNMatMatvecSetup prior to calling SUNMatMatvec. */
  if (ACS == NULL || SM_ROWTOPROC_SLUNRLOC(A) == NULL ||
      SM_COMMPATTERN_SLUNRLOC(A) == NULL)
    return(SUNMAT_MATVEC_SETUP_REQUIRED);

  xdata = N_VGetArrayPointer(x);
  ydata = N_VGetArrayPointer(y);
  if (xdata == NULL || ydata == NULL)
    return(SUNMAT_MEM_FAIL);

  /* Call SuperLU-DIST Matvec routine to perform the actual Matvec. */
  pdgsmv(0, ACS, SM_GRID_SLUNRLOC(A), SM_COMMPATTERN_SLUNRLOC(A), xdata, ydata);

  return(SUNMAT_SUCCESS);
}

int SUNMatMatvecSetup_SLUNRloc(SUNMatrix A)
{
  sunindextype *temp;
  sunindextype nprocs;
  sunindextype i, j;

  SuperMatrix *ACS = SM_COLSORTED_SLUNRLOC(A);

  /* If ACS is NULL, then this is the first setup call and things must be
     allocated */
  if (ACS == NULL) {
    ACS = (SuperMatrix *) malloc(sizeof(SuperMatrix));

    /* Clone and copy A to create ACS which will be A but with column-sorted
       column indices to [internal, external]. ACS is used with the Matvec
       routine. */
    dClone_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A), ACS);
    dCopy_CompRowLoc_Matrix_dist(SM_SUPERMATRIX_SLUNRLOC(A), ACS);

    SM_ROWTOPROC_SLUNRLOC(A) = (sunindextype *)
      malloc(SM_GLOBALROWS_SLUNRLOC(A)*sizeof(sunindextype));
    if (SM_ROWTOPROC_SLUNRLOC(A) == NULL) {
      Destroy_CompRowLoc_Matrix_dist(ACS); free(ACS);
      return(SUNMAT_MEM_FAIL);
    }

    SM_COMMPATTERN_SLUNRLOC(A) = (pdgsmv_comm_t *) malloc(sizeof(pdgsmv_comm_t));
    if (SM_COMMPATTERN_SLUNRLOC(A) == NULL) {
      free(SM_ROWTOPROC_SLUNRLOC(A));
      Destroy_CompRowLoc_Matrix_dist(ACS); free(ACS);
      return(SUNMAT_MEM_FAIL);
    }

    SM_COLSORTED_SLUNRLOC(A) = ACS;
  } else {
    /* if ACS has already been created, we can reuse it to save allocations,
       but we must finalize the last matvec to avoid leaking memory.*/
    pdgsmv_finalize(SM_COMMPATTERN_SLUNRLOC(A));
  }

  /* calculate the number of processes the matrix is spread across */
  nprocs = SM_GRID_SLUNRLOC(A)->nprow*SM_GRID_SLUNRLOC(A)->npcol;

  /* establish a row number to process mapping */
  temp = (sunindextype*) malloc((nprocs+1)*sizeof(sunindextype));
  if (temp == NULL) { SUNMatDestroy(A); return(SUNMAT_MEM_FAIL); }

  MPI_Allgather(&SM_FSTROW_SLUNRLOC(A), 1, MPI_SUNINDEXTYPE,
                temp, 1, MPI_SUNINDEXTYPE, SM_GRID_SLUNRLOC(A)->comm);

  temp[nprocs] = SM_GLOBALROWS_SLUNRLOC(A);
  for (i=0; i<nprocs; i++) {
    for (j=temp[i]; j<temp[i+1]; j++)
      SM_ROWTOPROC_SLUNRLOC(A)[j] = i;
  }
  free(temp);

  /* Call SuperLU-DIST routine to establish communication pattern for the
     Matvec. WARNING: This will overwrite the matrix provided. It is modified
     with colind permuted to [internal, external]. */
  pdgsmv_init(ACS, SM_ROWTOPROC_SLUNRLOC(A), SM_GRID_SLUNRLOC(A),
              SM_COMMPATTERN_SLUNRLOC(A));

  return(SUNMAT_SUCCESS);
}

int SUNMatSpace_SLUNRloc(SUNMatrix A, long int *lenrw, long int *leniw)
{
  /* since the SuperLU_DIST structures are opaque objects, we omit those
     from these results */

  *leniw = SM_GLOBALROWS_SLUNRLOC(A); /* length(row_to_proc) */
  *lenrw = 0;

  return(SUNMAT_SUCCESS);
}

/*
 * ----------------------------------------------------------------------------
 * private functions
 * ----------------------------------------------------------------------------
 */

/* Function to check compatibility of two sparse SUNMatrix objects.
   Checks to make sure that the the matrices are both SLUNRLOC,
   have the same number of rows, cols and nonzeros.*/
static booleantype SMCompatible_SLUNRloc(SUNMatrix A, SUNMatrix B)
{
  /* both matrices must be SLUNRLOC */
  if ((SUNMatGetID(A) != SUNMATRIX_SLUNRLOC) ||
      (SUNMatGetID(B) != SUNMATRIX_SLUNRLOC))
    return SUNFALSE;

  if (SM_GLOBALCOLS_SLUNRLOC(A) != SM_GLOBALCOLS_SLUNRLOC(B))
    return SUNFALSE;

  if (SM_LOCALROWS_SLUNRLOC(A) != SM_LOCALROWS_SLUNRLOC(B))
    return SUNFALSE;

  if (SM_LOCALNNZ_SLUNRLOC(A) != SM_LOCALNNZ_SLUNRLOC(B))
    return SUNFALSE;

  return SUNTRUE;
}
