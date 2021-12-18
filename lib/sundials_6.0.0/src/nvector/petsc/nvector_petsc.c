/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on N_Vector_Parallel by Scott D. Cohen, Alan C. Hindmarsh,
 * Radu Serban, and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a PETSc implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_petsc.h>
#include <sundials/sundials_math.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Error Message */
#define BAD_N1 "N_VNewEmpty_Petsc -- Sum of local vector lengths differs from "
#define BAD_N2 "input global length. \n\n"
#define BAD_N   BAD_N1 BAD_N2

/*
 * -----------------------------------------------------------------
 * Simplifying macros NV_CONTENT_PTC, NV_OWN_DATA_PTC,
 *                    NV_LOCLENGTH_PTC, NV_GLOBLENGTH_PTC,
 *                    NV_COMM_PTC
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * sunindextype v_len, s_len, i;
 *
 * (1) NV_CONTENT_PTC
 *
 *     This routines gives access to the contents of the PETSc
 *     vector wrapper N_Vector.
 *
 *     The assignment v_cont = NV_CONTENT_PTC(v) sets v_cont to be
 *     a pointer to the N_Vector (PETSc wrapper) content structure.
 *
 * (2) NV_PVEC_PTC, NV_OWN_DATA_PTC, NV_LOCLENGTH_PTC, NV_GLOBLENGTH_PTC,
 *     and NV_COMM_PTC
 *
 *     These routines give access to the individual parts of
 *     the content structure of a PETSc N_Vector wrapper.
 *
 *     NV_PVEC_PTC(v) returns the PETSc vector (Vec) object.
 *
 *     The assignment v_llen = NV_LOCLENGTH_PTC(v) sets v_llen to
 *     be the length of the local part of the vector v. The call
 *     NV_LOCLENGTH_PTC(v) = llen_v should NOT be used! It will
 *     change the value stored in the N_Vector content structure,
 *     but it will NOT change the length of the actual PETSc vector.
 *
 *     The assignment v_glen = NV_GLOBLENGTH_PTC(v) sets v_glen to
 *     be the global length of the vector v. The call
 *     NV_GLOBLENGTH_PTC(v) = glen_v should NOT be used! It will
 *     change the value stored in the N_Vector content structure,
 *     but it will NOT change the length of the actual PETSc vector.
 *
 *     The assignment v_comm = NV_COMM_PTC(v) sets v_comm to be the
 *     MPI communicator of the vector v. The assignment
 *     NV_COMM_PTC(v) = comm_v should NOT be used! It will change
 *     the value stored in the N_Vector content structure, but it
 *     will NOT change the MPI communicator of the actual PETSc
 *     vector.
 *
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_PTC(v)    ( (N_VectorContent_Petsc)(v->content) )

#define NV_LOCLENGTH_PTC(v)  ( NV_CONTENT_PTC(v)->local_length )

#define NV_GLOBLENGTH_PTC(v) ( NV_CONTENT_PTC(v)->global_length )

#define NV_OWN_DATA_PTC(v)   ( NV_CONTENT_PTC(v)->own_data )

#define NV_PVEC_PTC(v)       ( NV_CONTENT_PTC(v)->pvec )

#define NV_COMM_PTC(v)       ( NV_CONTENT_PTC(v)->comm )



/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Petsc(N_Vector v)
{
  return SUNDIALS_NVEC_PETSC;
}


/* ----------------------------------------------------------------
 * Function to create a new N_Vector wrapper with an empty (NULL)
 * PETSc vector.
 */

N_Vector N_VNewEmpty_Petsc(MPI_Comm comm,
                           sunindextype local_length,
                           sunindextype global_length,
                           SUNContext sunctx)
{
  N_Vector v;
  N_VectorContent_Petsc content;
  sunindextype n, Nsum;
  PetscErrorCode ierr;

  /* Compute global length as sum of local lengths */
  n = local_length;
  ierr = MPI_Allreduce(&n, &Nsum, 1, MPI_SUNINDEXTYPE, MPI_SUM, comm);
  CHKERRABORT(comm,ierr);
  if (Nsum != global_length) {
    STAN_SUNDIALS_FPRINTF(stderr, BAD_N);
    return(NULL);
  }

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Petsc;
  v->ops->nvclone           = N_VClone_Petsc;
  v->ops->nvcloneempty      = N_VCloneEmpty_Petsc;
  v->ops->nvdestroy         = N_VDestroy_Petsc;
  v->ops->nvspace           = N_VSpace_Petsc;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_Petsc;
  v->ops->nvsetarraypointer = N_VSetArrayPointer_Petsc;
  v->ops->nvgetcommunicator = N_VGetCommunicator_Petsc;
  v->ops->nvgetlength       = N_VGetLength_Petsc;

  /* standard vector operations */
  v->ops->nvlinearsum    = N_VLinearSum_Petsc;
  v->ops->nvconst        = N_VConst_Petsc;
  v->ops->nvprod         = N_VProd_Petsc;
  v->ops->nvdiv          = N_VDiv_Petsc;
  v->ops->nvscale        = N_VScale_Petsc;
  v->ops->nvabs          = N_VAbs_Petsc;
  v->ops->nvinv          = N_VInv_Petsc;
  v->ops->nvaddconst     = N_VAddConst_Petsc;
  v->ops->nvdotprod      = N_VDotProd_Petsc;
  v->ops->nvmaxnorm      = N_VMaxNorm_Petsc;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Petsc;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Petsc;
  v->ops->nvmin          = N_VMin_Petsc;
  v->ops->nvwl2norm      = N_VWL2Norm_Petsc;
  v->ops->nvl1norm       = N_VL1Norm_Petsc;
  v->ops->nvcompare      = N_VCompare_Petsc;
  v->ops->nvinvtest      = N_VInvTest_Petsc;
  v->ops->nvconstrmask   = N_VConstrMask_Petsc;
  v->ops->nvminquotient  = N_VMinQuotient_Petsc;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProdLocal_Petsc;
  v->ops->nvmaxnormlocal     = N_VMaxNormLocal_Petsc;
  v->ops->nvminlocal         = N_VMinLocal_Petsc;
  v->ops->nvl1normlocal      = N_VL1NormLocal_Petsc;
  v->ops->nvinvtestlocal     = N_VInvTestLocal_Petsc;
  v->ops->nvconstrmasklocal  = N_VConstrMaskLocal_Petsc;
  v->ops->nvminquotientlocal = N_VMinQuotientLocal_Petsc;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Petsc;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Petsc;

  /* single buffer reduction operations */
  v->ops->nvdotprodmultilocal     = N_VDotProdMultiLocal_Petsc;
  v->ops->nvdotprodmultiallreduce = N_VDotProdMultiAllReduce_Petsc;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_Petsc;
  v->ops->nvbufpack   = N_VBufPack_Petsc;
  v->ops->nvbufunpack = N_VBufUnpack_Petsc;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Petsc) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_data      = SUNFALSE;
  content->pvec          = NULL;

  return(v);
}



/* ----------------------------------------------------------------
 * Function to create an N_Vector wrapper for a PETSc vector.
 */

N_Vector N_VMake_Petsc(Vec pvec, SUNContext sunctx)
{
  N_Vector v = NULL;
  MPI_Comm comm;
  PetscInt local_length;
  PetscInt global_length;

  VecGetLocalSize(pvec, &local_length);
  VecGetSize(pvec, &global_length);
  PetscObjectGetComm((PetscObject) pvec, &comm);

  v = N_VNewEmpty_Petsc(comm, local_length, global_length, sunctx);
  if (v == NULL)
     return(NULL);

  /* Attach data */
  NV_OWN_DATA_PTC(v) = SUNFALSE;
  NV_PVEC_PTC(v)     = pvec;

  return(v);
}

/* ----------------------------------------------------------------
 * Function to create an array of new PETSc vector wrappers.
 */

N_Vector *N_VCloneVectorArray_Petsc(int count, N_Vector w)
{
  return(N_VCloneVectorArray(count, w));
}

/* ----------------------------------------------------------------
 * Function to create an array of new PETSc vector wrappers with
 * empty (NULL) PETSc vectors.
 */

N_Vector *N_VCloneVectorArrayEmpty_Petsc(int count, N_Vector w)
{
  return(N_VCloneEmptyVectorArray(count, w));
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Petsc
 */

void N_VDestroyVectorArray_Petsc(N_Vector *vs, int count)
{
  N_VDestroyVectorArray(vs, count);
  return;
}

/* ----------------------------------------------------------------
 * Function to extract PETSc vector
 */

Vec N_VGetVector_Petsc(N_Vector v)
{
  return NV_PVEC_PTC(v);
}

/* ----------------------------------------------------------------
 * Function to set the PETSc vector
 */

void N_VSetVector_Petsc(N_Vector v, Vec p)
{
  NV_PVEC_PTC(v) = p;
}

/* ----------------------------------------------------------------
 * Function to print the global data in a PETSc vector to stdout
 */

void N_VPrint_Petsc(N_Vector x)
{
  Vec xv = NV_PVEC_PTC(x);
  MPI_Comm comm = NV_COMM_PTC(x);

  VecView(xv, PETSC_VIEWER_STDOUT_(comm));

  return;
}

/* ----------------------------------------------------------------
 * Function to print the global data in a PETSc vector to fname
 */

void N_VPrintFile_Petsc(N_Vector x, const char fname[])
{
  Vec xv = NV_PVEC_PTC(x);
  MPI_Comm comm = NV_COMM_PTC(x);
  PetscViewer viewer;

  PetscViewerASCIIOpen(comm, fname, &viewer);

  VecView(xv, viewer);

  PetscViewerDestroy(&viewer);

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Petsc(N_Vector w)
{
  N_Vector v;
  N_VectorContent_Petsc content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Petsc) malloc(sizeof *content);
  if (content == NULL) { N_VDestroy(v); return(NULL); }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->local_length  = NV_LOCLENGTH_PTC(w);
  content->global_length = NV_GLOBLENGTH_PTC(w);
  content->comm          = NV_COMM_PTC(w);
  content->own_data      = SUNFALSE;
  content->pvec          = NULL;

  return(v);
}

N_Vector N_VClone_Petsc(N_Vector w)
{
  N_Vector v = NULL;
  Vec pvec   = NULL;
  Vec wvec   = NV_PVEC_PTC(w);

  /* PetscErrorCode ierr; */

  v = N_VCloneEmpty_Petsc(w);
  if (v == NULL)
    return(NULL);

  /* Duplicate vector */
  VecDuplicate(wvec, &pvec);
  if (pvec == NULL) {
    N_VDestroy_Petsc(v);
    return(NULL);
  }

  /* Attach data */
  NV_OWN_DATA_PTC(v) = SUNTRUE;
  NV_PVEC_PTC(v)     = pvec;

  return(v);
}

void N_VDestroy_Petsc(N_Vector v)
{
  if (v == NULL) return;

  /* free content */
  if (v->content != NULL) {
    if (NV_OWN_DATA_PTC(v) && NV_PVEC_PTC(v) != NULL) {
      VecDestroy(&(NV_PVEC_PTC(v)));
      NV_PVEC_PTC(v) = NULL;
    }
    free(v->content);
    v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL) { free(v->ops); v->ops = NULL; }
  free(v); v = NULL;

  return;
}

void N_VSpace_Petsc(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  MPI_Comm comm;
  int npes;

  comm = NV_COMM_PTC(v);
  MPI_Comm_size(comm, &npes);

  *lrw = NV_GLOBLENGTH_PTC(v);
  *liw = 2*npes;

  return;
}

/*
 * Not implemented for PETSc wrapper.
 */
realtype *N_VGetArrayPointer_Petsc(N_Vector v)
{
  return NULL;
}

/*
 * Not implemented for PETSc wrapper.
 */
void N_VSetArrayPointer_Petsc(realtype *v_data, N_Vector v)
{
  return;
}

void *N_VGetCommunicator_Petsc(N_Vector v)
{
  return((void *) &(NV_COMM_PTC(v)));
}

sunindextype N_VGetLength_Petsc(N_Vector v)
{
  return(NV_GLOBLENGTH_PTC(v));
}

void N_VLinearSum_Petsc(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  Vec zv = NV_PVEC_PTC(z);

  if (x == y) {
    N_VScale_Petsc(a + b, x, z); /* z <~ ax+bx */
    return;
  }

  if (z == y) {
    if (b == ONE) {
      VecAXPY(yv, a, xv);   /* BLAS usage: axpy  y <- ax+y */
      return;
    }
    VecAXPBY(yv, a, b, xv); /* BLAS usage: axpby y <- ax+by */
    return;
  }

  if (z == x) {
    if (a == ONE) {
      VecAXPY(xv, b, yv);   /* BLAS usage: axpy  x <- by+x */
      return;
    }
    VecAXPBY(xv, b, a, yv); /* BLAS usage: axpby x <- by+ax */
    return;
  }


  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  VecAXPBYPCZ(zv, a, b, 0.0, xv, yv); /* PETSc, probably not optimal */

  return;
}

void N_VConst_Petsc(realtype c, N_Vector z)
{
  Vec zv = NV_PVEC_PTC(z);

  VecSet(zv, c);

  return;
}

void N_VProd_Petsc(N_Vector x, N_Vector y, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  Vec zv = NV_PVEC_PTC(z);

  VecPointwiseMult(zv, xv, yv);

  return;
}

void N_VDiv_Petsc(N_Vector x, N_Vector y, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  Vec zv = NV_PVEC_PTC(z);

  VecPointwiseDivide(zv, xv, yv); /* z = x/y */

  return;
}

void N_VScale_Petsc(realtype c, N_Vector x, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VecScale(xv, c);
    return;
  }

  VecAXPBY(zv, c, 0.0, xv);

  return;
}

void N_VAbs_Petsc(N_Vector x, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);

  if(z != x)
    VecCopy(xv, zv); /* copy x~>z */
  VecAbs(zv);

  return;
}

void N_VInv_Petsc(N_Vector x, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);

  if(z != x)
    VecCopy(xv, zv); /* copy x~>z */
  VecReciprocal(zv);

  return;
}

void N_VAddConst_Petsc(N_Vector x, realtype b, N_Vector z)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);

  if(z != x)
    VecCopy(xv, zv); /* copy x~>z */
  VecShift(zv, b);

  return;
}

realtype N_VDotProdLocal_Petsc(N_Vector x, N_Vector y)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  PetscScalar *xd;
  PetscScalar *yd;
  PetscReal sum = ZERO;

  VecGetArray(xv, &xd);
  VecGetArray(yv, &yd);
  for (i = 0; i < N; i++)
    sum += xd[i] * yd[i];
  VecRestoreArray(xv, &xd);
  VecRestoreArray(yv, &yd);
  return ((realtype) sum);
}

realtype N_VDotProd_Petsc(N_Vector x, N_Vector y)
{
  Vec xv = NV_PVEC_PTC(x);
  Vec yv = NV_PVEC_PTC(y);
  PetscScalar dotprod;

  VecDot(xv, yv, &dotprod);
  return dotprod;
}

realtype N_VMaxNormLocal_Petsc(N_Vector x)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  PetscScalar *xd;
  PetscReal max = ZERO;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++) {
    if (PetscAbsScalar(xd[i]) > max) max = PetscAbsScalar(xd[i]);
  }
  VecRestoreArray(xv, &xd);
  return ((realtype) max);
}

realtype N_VMaxNorm_Petsc(N_Vector x)
{
  Vec xv = NV_PVEC_PTC(x);
  PetscReal norm;

  VecNorm(xv, NORM_INFINITY, &norm);
  return norm;
}

realtype N_VWSqrSumLocal_Petsc(N_Vector x, N_Vector w)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  Vec wv = NV_PVEC_PTC(w);
  PetscScalar *xd;
  PetscScalar *wd;
  PetscReal sum = ZERO;

  VecGetArray(xv, &xd);
  VecGetArray(wv, &wd);
  for (i = 0; i < N; i++) {
    sum += PetscSqr(PetscAbsScalar(xd[i] * wd[i]));
  }
  VecRestoreArray(xv, &xd);
  VecRestoreArray(wv, &wd);
  return ((realtype) sum);
}

realtype N_VWrmsNorm_Petsc(N_Vector x, N_Vector w)
{
  realtype global_sum;
  sunindextype N_global = NV_GLOBLENGTH_PTC(x);
  realtype sum = N_VWSqrSumLocal_Petsc(x, w);
  (void) MPI_Allreduce(&sum, &global_sum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_COMM_PTC(x));
  return (SUNRsqrt(global_sum/N_global));
}

realtype N_VWSqrSumMaskLocal_Petsc(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  Vec wv = NV_PVEC_PTC(w);
  Vec idv = NV_PVEC_PTC(id);
  PetscScalar *xd;
  PetscScalar *wd;
  PetscScalar *idd;
  PetscReal sum = ZERO;

  VecGetArray(xv, &xd);
  VecGetArray(wv, &wd);
  VecGetArray(idv, &idd);
  for (i = 0; i < N; i++) {
    PetscReal tag = (PetscReal) idd[i];
    if (tag > ZERO) {
      sum += PetscSqr(PetscAbsScalar(xd[i] * wd[i]));
    }
  }
  VecRestoreArray(xv, &xd);
  VecRestoreArray(wv, &wd);
  VecRestoreArray(idv, &idd);
  return sum;
}

realtype N_VWrmsNormMask_Petsc(N_Vector x, N_Vector w, N_Vector id)
{
  realtype global_sum;
  sunindextype N_global = NV_GLOBLENGTH_PTC(x);
  realtype sum = N_VWSqrSumMaskLocal_Petsc(x, w, id);
  (void) MPI_Allreduce(&sum, &global_sum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_COMM_PTC(x));
  return (SUNRsqrt(global_sum/N_global));
}

realtype N_VMinLocal_Petsc(N_Vector x)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  PetscScalar *xd;
  PetscReal min = BIG_REAL;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++) {
    if (xd[i] < min) min = xd[i];
  }
  VecRestoreArray(xv, &xd);
  return ((realtype) min);
}

realtype N_VMin_Petsc(N_Vector x)
{
  Vec xv = NV_PVEC_PTC(x);
  PetscReal minval;
  PetscInt i;

  VecMin(xv, &i, &minval);
  return minval;
}

realtype N_VWL2Norm_Petsc(N_Vector x, N_Vector w)
{
  realtype global_sum;
  realtype sum = N_VWSqrSumLocal_Petsc(x, w);
  (void) MPI_Allreduce(&sum, &global_sum, 1, MPI_SUNREALTYPE, MPI_SUM, NV_COMM_PTC(x));
  return (SUNRsqrt(global_sum));
}

realtype N_VL1NormLocal_Petsc(N_Vector x)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  PetscScalar *xd;
  PetscReal sum = ZERO;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++) {
    sum += PetscAbsScalar(xd[i]);
  }
  VecRestoreArray(xv, &xd);
  return ((realtype) sum);
}

realtype N_VL1Norm_Petsc(N_Vector x)
{
  Vec xv = NV_PVEC_PTC(x);
  PetscReal norm;

  VecNorm(xv, NORM_1, &norm);
  return norm;
}

void N_VCompare_Petsc(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);
  PetscReal cpet = c; /* <~ realtype should typedef to PETScReal */
  PetscScalar *xdata;
  PetscScalar *zdata;

  VecGetArray(xv, &xdata);
  VecGetArray(zv, &zdata);
  for (i = 0; i < N; i++) {
    zdata[i] = PetscAbsScalar(xdata[i]) >= cpet ? ONE : ZERO;
  }
  VecRestoreArray(xv, &xdata);
  VecRestoreArray(zv, &zdata);

  return;
}

booleantype N_VInvTestLocal_Petsc(N_Vector x, N_Vector z)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  Vec zv = NV_PVEC_PTC(z);
  PetscScalar *xd;
  PetscScalar *zd;
  PetscReal val = ONE;

  VecGetArray(xv, &xd);
  VecGetArray(zv, &zd);
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO)
      val = ZERO;
    else
      zd[i] = ONE/xd[i];
  }
  VecRestoreArray(xv, &xd);
  VecRestoreArray(zv, &zd);

  if (val == ZERO)
    return(SUNFALSE);
  else
    return(SUNTRUE);
}

booleantype N_VInvTest_Petsc(N_Vector x, N_Vector z)
{
  realtype val2;
  realtype val = (N_VInvTestLocal_Petsc(x, z)) ? ONE : ZERO;
  (void) MPI_Allreduce(&val, &val2, 1, MPI_SUNREALTYPE, MPI_MIN, NV_COMM_PTC(x));
  if (val2 == ZERO)
    return(SUNFALSE);
  else
    return(SUNTRUE);
}

booleantype N_VConstrMaskLocal_Petsc(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  realtype temp;
  booleantype test;
  Vec xv = NV_PVEC_PTC(x);
  Vec cv = NV_PVEC_PTC(c);
  Vec mv = NV_PVEC_PTC(m);
  PetscScalar *xd;
  PetscScalar *cd;
  PetscScalar *md;

  temp = ZERO;

  VecGetArray(xv, &xd);
  VecGetArray(cv, &cd);
  VecGetArray(mv, &md);
  for (i = 0; i < N; i++) {
    PetscReal cc = (PetscReal) cd[i]; /* <~ Drop imaginary parts if any. */
    PetscReal xx = (PetscReal) xd[i]; /* <~ Constraints defined on Re{x} */
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cc == ZERO)
      continue;

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cc) > ONEPT5 && xx*cc <= ZERO) ||
           (SUNRabs(cc) > HALF   && xx*cc <  ZERO);
    if (test) {
      temp = md[i] = ONE;
    }
  }
  VecRestoreArray(xv, &xd);
  VecRestoreArray(cv, &cd);
  VecRestoreArray(mv, &md);

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

booleantype N_VConstrMask_Petsc(N_Vector c, N_Vector x, N_Vector m)
{
  realtype temp2;
  realtype temp = (N_VConstrMaskLocal_Petsc(c, x, m)) ? ZERO : ONE;
  (void) MPI_Allreduce(&temp, &temp2, 1, MPI_SUNREALTYPE, MPI_MAX, NV_COMM_PTC(x));
  return (temp2 == ONE) ? SUNFALSE : SUNTRUE;
}

realtype N_VMinQuotientLocal_Petsc(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce = SUNTRUE;
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(num);

  Vec nv = NV_PVEC_PTC(num);
  Vec dv = NV_PVEC_PTC(denom);
  PetscScalar *nd;
  PetscScalar *dd;
  PetscReal minval = BIG_REAL;

  VecGetArray(nv, &nd);
  VecGetArray(dv, &dd);
  for (i = 0; i < N; i++) {
    PetscReal nr = (PetscReal) nd[i];
    PetscReal dr = (PetscReal) dd[i];
    if (dr == ZERO)
      continue;
    else {
      if (!notEvenOnce)
        minval = SUNMIN(minval, nr/dr);
      else {
        minval = nr/dr;
        notEvenOnce = SUNFALSE;
      }
    }
  }
  VecRestoreArray(nv, &nd);
  VecRestoreArray(dv, &dd);
  return((realtype) minval);
}

realtype N_VMinQuotient_Petsc(N_Vector num, N_Vector denom)
{
  PetscReal gmin;
  realtype minval = N_VMinQuotientLocal_Petsc(num, denom);
  (void) MPI_Allreduce(&minval, &gmin, 1, MPI_SUNREALTYPE, MPI_MIN, NV_COMM_PTC(num));
  return(gmin);
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */


int N_VLinearCombination_Petsc(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int  i;
  Vec* xv;
  Vec  zv;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Petsc(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_Petsc(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  /* get petsc vectors */
  xv = (Vec*) malloc(nvec * sizeof(Vec));
  for (i=0; i<nvec; i++)
    xv[i] = NV_PVEC_PTC(X[i]);

  zv = NV_PVEC_PTC(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE)) {
    VecMAXPY(zv, nvec-1, c+1, xv+1);
    free(xv);
    return(0);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z) {
    VecScale(zv, c[0]);
    VecMAXPY(zv, nvec-1, c+1, xv+1);
    free(xv);
    return(0);
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  VecAXPBY(zv, c[0], 0.0, xv[0]);
  VecMAXPY(zv, nvec-1, c+1, xv+1);
  free(xv);

  return(0);
}


int N_VScaleAddMulti_Petsc(int nvec, realtype* a, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  PetscScalar  *xd, *yd, *zd;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Petsc(a[0], x, ONE, Y[0], Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N = NV_LOCLENGTH_PTC(x);
  VecGetArray(NV_PVEC_PTC(x), &xd);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      VecGetArray(NV_PVEC_PTC(Y[i]), &yd);
      for (j=0; j<N; j++) {
        yd[j] += a[i] * xd[j];
      }
      VecRestoreArray(NV_PVEC_PTC(Y[i]), &yd);
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    VecGetArray(NV_PVEC_PTC(Y[i]), &yd);
    VecGetArray(NV_PVEC_PTC(Z[i]), &zd);
    for (j=0; j<N; j++) {
      zd[j] = a[i] * xd[j] + yd[j];
    }
    VecRestoreArray(NV_PVEC_PTC(Y[i]), &yd);
    VecRestoreArray(NV_PVEC_PTC(Z[i]), &zd);
  }
  return(0);
}


int N_VDotProdMulti_Petsc(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  int  i;
  Vec* yv;
  Vec  xv;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VDotProd */
  if (nvec == 1) {
    dotprods[0] = N_VDotProd_Petsc(x, Y[0]);
    return(0);
  }

  /* get petsc vectors */
  yv = (Vec*) malloc(nvec * sizeof(Vec));
  for (i=0; i<nvec; i++)
    yv[i] = NV_PVEC_PTC(Y[i]);

  xv = NV_PVEC_PTC(x);

  VecMDot(xv, nvec, yv, dotprods);
  free(yv);

  return(0);
}

/*
 * -----------------------------------------------------------------------------
 * single buffer reduction operations
 * -----------------------------------------------------------------------------
 */

int N_VDotProdMultiLocal_Petsc(int nvec, N_Vector x, N_Vector* Y,
                               realtype* dotprods)
{
  int j;
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec xv = NV_PVEC_PTC(x);
  Vec yv;
  PetscScalar *xd;
  PetscScalar *yd;

  VecGetArray(xv, &xd);
  for (j = 0; j < nvec; j++) {
    yv = NV_PVEC_PTC(Y[j]);
    VecGetArray(yv, &yd);
    dotprods[j] = ZERO;
    for (i = 0; i < N; i++) {
      dotprods[j] += xd[i] * yd[i];
    }
    VecRestoreArray(yv, &yd);
  }
  VecRestoreArray(xv, &xd);

  return (0);
}

int N_VDotProdMultiAllReduce_Petsc(int nvec, N_Vector x, realtype* dotprods)
{
  int retval;
  retval = MPI_Allreduce(MPI_IN_PLACE, dotprods, nvec, MPI_SUNREALTYPE,
                         MPI_SUM, NV_COMM_PTC(x));
  return retval == MPI_SUCCESS ? 0 : -1;
}

/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Petsc(int nvec,
                                  realtype a, N_Vector* X,
                                  realtype b, N_Vector* Y,
                                  N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  PetscScalar  *xd, *yd, *zd;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_Petsc(a, X[0], b, Y[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LOCLENGTH_PTC(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i=0; i<nvec; i++) {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    VecGetArray(NV_PVEC_PTC(Y[i]), &yd);
    VecGetArray(NV_PVEC_PTC(Z[i]), &zd);
    for (j=0; j<N; j++) {
      zd[j] = a * xd[j] + b * yd[j];
    }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    VecRestoreArray(NV_PVEC_PTC(Y[i]), &yd);
    VecRestoreArray(NV_PVEC_PTC(Z[i]), &zd);
  }

  return(0);
}


int N_VScaleVectorArray_Petsc(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  PetscScalar  *xd, *zd;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_Petsc(c[0], X[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LOCLENGTH_PTC(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z) {
    for (i=0; i<nvec; i++) {
      VecGetArray(NV_PVEC_PTC(X[i]), &xd);
      for (j=0; j<N; j++) {
        xd[j] *= c[i];
      }
      VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    }
    return(0);
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i=0; i<nvec; i++) {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    VecGetArray(NV_PVEC_PTC(Z[i]), &zd);
    for (j=0; j<N; j++) {
      zd[j] = c[i] * xd[j];
    }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    VecRestoreArray(NV_PVEC_PTC(Z[i]), &zd);
  }
  return(0);
}


int N_VConstVectorArray_Petsc(int nvec, realtype c, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  PetscScalar  *zd;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VConst */
  if (nvec == 1) {
    N_VConst_Petsc(c, Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LOCLENGTH_PTC(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i=0; i<nvec; i++) {
    VecGetArray(NV_PVEC_PTC(Z[i]), &zd);
    for (j=0; j<N; j++) {
      zd[j] = c;
    }
    VecRestoreArray(NV_PVEC_PTC(Z[i]), &zd);
  }

  return(0);
}


int N_VWrmsNormVectorArray_Petsc(int nvec, N_Vector* X, N_Vector* W, realtype* nrm)
{
  int          i, retval;
  sunindextype j, Nl, Ng;
  realtype*    wd=NULL;
  realtype*    xd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNorm_Petsc(X[0], W[0]);
    return(0);
  }

  /* get vector lengths and communicator */
  Nl   = NV_LOCLENGTH_PTC(X[0]);
  Ng   = NV_GLOBLENGTH_PTC(X[0]);
  comm = NV_COMM_PTC(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    VecGetArray(NV_PVEC_PTC(W[i]), &wd);
    nrm[i] = ZERO;
    for (j=0; j<Nl; j++) {
      nrm[i] += PetscSqr(PetscAbsScalar(xd[j] * wd[j]));
    }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    VecRestoreArray(NV_PVEC_PTC(W[i]), &wd);
  }
  retval = MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE, MPI_SUM, comm);

  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/Ng);

  return retval == MPI_SUCCESS ? 0 : -1;
}


int N_VWrmsNormMaskVectorArray_Petsc(int nvec, N_Vector* X, N_Vector* W,
                                        N_Vector id, realtype* nrm)
{
  int          i, retval;
  sunindextype j, Nl, Ng;
  PetscScalar  *wd, *xd, *idd;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNormMask_Petsc(X[0], W[0], id);
    return(0);
  }

  /* get vector lengths and communicator */
  Nl   = NV_LOCLENGTH_PTC(X[0]);
  Ng   = NV_GLOBLENGTH_PTC(X[0]);
  comm = NV_COMM_PTC(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  VecGetArray(NV_PVEC_PTC(id), &idd);
  for (i=0; i<nvec; i++) {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    VecGetArray(NV_PVEC_PTC(W[i]), &wd);
    nrm[i] = ZERO;
    for (j=0; j<Nl; j++) {
      if (idd[j] > ZERO)
        nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    VecRestoreArray(NV_PVEC_PTC(W[i]), &wd);
  }
  VecRestoreArray(NV_PVEC_PTC(id), &idd);

  retval = MPI_Allreduce(MPI_IN_PLACE, nrm, nvec, MPI_SUNREALTYPE, MPI_SUM, comm);

  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/Ng);

  return retval == MPI_SUCCESS ? 0 : -1;
}


int N_VScaleAddMultiVectorArray_Petsc(int nvec, int nsum, realtype* a,
                                          N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  int          i, j;
  sunindextype k, N;
  PetscScalar  *xd, *yd, *zd;

  int          retval;
  N_Vector*    YY;
  N_Vector*    ZZ;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VLinearSum */
    if (nsum == 1) {
      N_VLinearSum_Petsc(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      return(0);
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector *) malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (j=0; j<nsum; j++) {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    retval = N_VScaleAddMulti_Petsc(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);
    return(retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1) {
    retval = N_VLinearSumVectorArray_Petsc(nvec, a[0], X, ONE, Y[0], Z[0]);
    return(retval);
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N  = NV_LOCLENGTH_PTC(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      VecGetArray(NV_PVEC_PTC(X[i]), &xd);
      for (j=0; j<nsum; j++) {
        VecGetArray(NV_PVEC_PTC(Y[j][i]), &yd);
        for (k=0; k<N; k++) {
          yd[k] += a[j] * xd[k];
        }
        VecRestoreArray(NV_PVEC_PTC(Y[j][i]), &yd);
      }
      VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    VecGetArray(NV_PVEC_PTC(X[i]), &xd);
    for (j=0; j<nsum; j++) {
      VecGetArray(NV_PVEC_PTC(Y[j][i]), &yd);
      VecGetArray(NV_PVEC_PTC(Z[j][i]), &zd);
      for (k=0; k<N; k++) {
        zd[k] = a[j] * xd[k] + yd[k];
      }
      VecRestoreArray(NV_PVEC_PTC(Y[j][i]), &yd);
      VecRestoreArray(NV_PVEC_PTC(Z[j][i]), &zd);
    }
    VecRestoreArray(NV_PVEC_PTC(X[i]), &xd);
  }

  return(0);
}


int N_VLinearCombinationVectorArray_Petsc(int nvec, int nsum,
                                             realtype* c,
                                             N_Vector** X,
                                             N_Vector* Z)
{
  int          i; /* vector arrays index in summation [0,nsum) */
  int          j; /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  PetscScalar  *zd, *xd;

  realtype*    ctmp;
  N_Vector*    Y;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);
  if (nsum < 1) return(-1);

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1) {

    /* should have called N_VScale */
    if (nsum == 1) {
      N_VScale_Petsc(c[0], X[0][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearSum */
    if (nsum == 2) {
      N_VLinearSum_Petsc(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nsum; i++) {
      Y[i] = X[i][0];
    }

    N_VLinearCombination_Petsc(nsum, c, Y, Z[0]);

    free(Y);
    return(0);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VScaleVectorArray */
  if (nsum == 1) {

    ctmp = (realtype*) malloc(nvec * sizeof(realtype));

    for (j=0; j<nvec; j++) {
      ctmp[j] = c[0];
    }

    N_VScaleVectorArray_Petsc(nvec, ctmp, X[0], Z);

    free(ctmp);
    return(0);
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2) {
    N_VLinearSumVectorArray_Petsc(nvec, c[0], X[0], c[1], X[1], Z);
    return(0);
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LOCLENGTH_PTC(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE)) {
    for (j=0; j<nvec; j++) {
      VecGetArray(NV_PVEC_PTC(Z[j]), &zd);
      for (i=1; i<nsum; i++) {
        VecGetArray(NV_PVEC_PTC(X[i][j]), &xd);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
        VecRestoreArray(NV_PVEC_PTC(X[i][j]), &xd);
      }
      VecRestoreArray(NV_PVEC_PTC(Z[j]), &zd);
    }
    return(0);
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z) {
    for (j=0; j<nvec; j++) {
      VecGetArray(NV_PVEC_PTC(Z[j]), &zd);
      for (k=0; k<N; k++) {
        zd[k] *= c[0];
      }
      for (i=1; i<nsum; i++) {
        VecGetArray(NV_PVEC_PTC(X[i][j]), &xd);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
        VecRestoreArray(NV_PVEC_PTC(X[i][j]), &xd);
      }
      VecRestoreArray(NV_PVEC_PTC(Z[j]), &zd);
    }
    return(0);
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j=0; j<nvec; j++) {
    VecGetArray(NV_PVEC_PTC(X[0][j]), &xd);
    VecGetArray(NV_PVEC_PTC(Z[j]), &zd);
    for (k=0; k<N; k++) {
      zd[k] = c[0] * xd[k];
    }
    VecRestoreArray(NV_PVEC_PTC(X[0][j]), &xd);
    for (i=1; i<nsum; i++) {
      VecGetArray(NV_PVEC_PTC(X[i][j]), &xd);
      for (k=0; k<N; k++) {
        zd[k] += c[i] * xd[k];
      }
      VecRestoreArray(NV_PVEC_PTC(X[i][j]), &xd);
    }
    VecRestoreArray(NV_PVEC_PTC(Z[j]), &zd);
  }
  return(0);
}


/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */


int N_VBufSize_Petsc(N_Vector x, sunindextype *size)
{
  if (x == NULL) return(-1);
  *size = NV_LOCLENGTH_PTC(x) * ((sunindextype)sizeof(PetscScalar));
  return(0);
}


int N_VBufPack_Petsc(N_Vector x, void *buf)
{
  Vec          xv;
  sunindextype i, N;
  PetscScalar  *xd = NULL;
  PetscScalar  *bd = NULL;

  if (x == NULL || buf == NULL) return(-1);

  xv = NV_PVEC_PTC(x);
  N  = NV_LOCLENGTH_PTC(x);
  bd = (PetscScalar*) buf;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++)
    bd[i] = xd[i];
  VecRestoreArray(xv, &xd);

  return(0);
}


int N_VBufUnpack_Petsc(N_Vector x, void *buf)
{
  Vec          xv;
  sunindextype i, N;
  PetscScalar  *xd = NULL;
  PetscScalar  *bd = NULL;

  if (x == NULL || buf == NULL) return(-1);

  xv = NV_PVEC_PTC(x);
  N  = NV_LOCLENGTH_PTC(x);
  bd = (PetscScalar*) buf;

  VecGetArray(xv, &xd);
  for (i = 0; i < N; i++)
    xd[i] = bd[i];
  VecRestoreArray(xv, &xd);

  return(0);
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Petsc;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Petsc;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Petsc;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Petsc;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Petsc;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Petsc;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Petsc;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Petsc;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Petsc;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Petsc;
    /* enable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = N_VDotProdMultiLocal_Petsc;
  } else {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
    /* disable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Petsc;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Petsc;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_Petsc;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Petsc;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Petsc;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Petsc;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Petsc;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Petsc;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Petsc;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Petsc;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMultiLocal_Petsc(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmultilocal = N_VDotProdMultiLocal_Petsc;
  else
    v->ops->nvdotprodmultilocal = NULL;

  /* return success */
  return(0);
}
