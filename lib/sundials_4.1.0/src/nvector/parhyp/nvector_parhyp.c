/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL and Jean M. Sexton @ SMU
 * -----------------------------------------------------------------
 * Based on work by Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                  and Aaron Collier @ LLNL
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
 * This is the implementation file for a HYPRE ParVector wrapper
 * for the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_mpi.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

/* Error Message */
#define BAD_N1 "N_VNew_ParHyp -- Sum of local vector lengths differs from "
#define BAD_N2 "input global length. \n\n"
#define BAD_N   BAD_N1 BAD_N2

/*
 * -----------------------------------------------------------------
 * Simplifying macros NV_CONTENT_PH, NV_DATA_PH, NV_LOCLENGTH_PH,
 *                    NV_GLOBLENGTH_PH, and NV_COMM_PH
 * -----------------------------------------------------------------
 * In the descriptions below, the following user declarations
 * are assumed:
 *
 * N_Vector v;
 * sunindextype v_len, s_len, i;
 *
 * (1) NV_CONTENT_PH
 *
 *     This routines gives access to the contents of the HYPRE
 *     vector wrapper (the N_Vector).
 *
 *     The assignment v_cont = NV_CONTENT_PH(v) sets v_cont to be
 *     a pointer to the N_Vector content structure.
 *
 * (2) NV_DATA_PH, NV_LOCLENGTH_PH, NV_GLOBLENGTH_PH, and NV_COMM_PH
 *
 *     These routines give access to the individual parts of
 *     the content structure of a parhyp N_Vector.
 *
 *     The assignment v_llen = NV_LOCLENGTH_PH(v) sets v_llen to
 *     be the length of the local part of the vector v. The call
 *     NV_LOCLENGTH_PH(v) = llen_v generally should NOT be used! It
 *     will change locally stored value with the HYPRE local vector
 *     length, but it will NOT change the length of the actual HYPRE
 *     local vector.
 *
 *     The assignment v_glen = NV_GLOBLENGTH_PH(v) sets v_glen to
 *     be the global length of the vector v. The call
 *     NV_GLOBLENGTH_PH(v) = glen_v generally should NOT be used! It
 *     will change locally stored value with the HYPRE parallel vector
 *     length, but it will NOT change the length of the actual HYPRE
 *     parallel vector.
 *
 *     The assignment v_comm = NV_COMM_PH(v) sets v_comm to be the
 *     MPI communicator of the vector v. The assignment
 *     NV_COMM_C(v) = comm_v sets the MPI communicator of v to be
 *     NV_COMM_PH(v) = comm_v generally should NOT be used! It
 *     will change locally stored value with the HYPRE parallel vector
 *     communicator, but it will NOT change the communicator of the
 *     actual HYPRE parallel vector.
 *
 * (3) NV_DATA_PH, NV_HYPRE_PARVEC_PH
 *
 *     The assignment v_data = NV_DATA_PH(v) sets v_data to be
 *     a pointer to the first component of the data inside the
 *     local vector of the HYPRE_parhyp vector for the vector v.
 *     The assignment NV_DATA_PH(v) = data_v should NOT be used.
 *     Instead, use NV_HYPRE_PARVEC_PH to obtain pointer to HYPRE
 *     vector and then use HYPRE functions to manipulate vector data.
 *
 *     The assignment v_parhyp = NV_HYPRE_PARVEC_PH(v) sets v_parhyp
 *     to be a pointer to HYPRE_ParVector of vector v. The assignment
 *     NV_HYPRE_PARVEC_PH(v) = parhyp_v sets pointer to
 *     HYPRE_ParVector of vector v to be parhyp_v.
 *
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_PH(v)    ( (N_VectorContent_ParHyp)(v->content) )

#define NV_LOCLENGTH_PH(v)  ( NV_CONTENT_PH(v)->local_length )

#define NV_GLOBLENGTH_PH(v) ( NV_CONTENT_PH(v)->global_length )

#define NV_OWN_PARVEC_PH(v) ( NV_CONTENT_PH(v)->own_parvector )

#define NV_HYPRE_PARVEC_PH(v) ( NV_CONTENT_PH(v)->x )

#define NV_DATA_PH(v)       ( NV_HYPRE_PARVEC_PH(v) == NULL ? NULL : hypre_VectorData(hypre_ParVectorLocalVector(NV_HYPRE_PARVEC_PH(v))) )

#define NV_COMM_PH(v)       ( NV_CONTENT_PH(v)->comm )


/* Private function prototypes */

/* z=x+y */
static void VSum_ParHyp(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_ParHyp(N_Vector x, N_Vector y, N_Vector z);
/* z=c(x+y) */
static void VScaleSum_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=ax+y */
static void VLin1_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_ParHyp(N_Vector v)
{
  return SUNDIALS_NVEC_PARHYP;
}


/* ----------------------------------------------------------------
 * Function to create a new parhyp vector without underlying
 * HYPRE vector.
 */
N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm,
                            sunindextype local_length,
                            sunindextype global_length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_ParHyp content;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = N_VGetVectorID_ParHyp;
  ops->nvclone           = N_VClone_ParHyp;
  ops->nvcloneempty      = N_VCloneEmpty_ParHyp;
  ops->nvdestroy         = N_VDestroy_ParHyp;
  ops->nvspace           = N_VSpace_ParHyp;
  ops->nvgetarraypointer = N_VGetArrayPointer_ParHyp;
  ops->nvsetarraypointer = N_VSetArrayPointer_ParHyp;

  /* standard vector operations */
  ops->nvlinearsum    = N_VLinearSum_ParHyp;
  ops->nvconst        = N_VConst_ParHyp;
  ops->nvprod         = N_VProd_ParHyp;
  ops->nvdiv          = N_VDiv_ParHyp;
  ops->nvscale        = N_VScale_ParHyp;
  ops->nvabs          = N_VAbs_ParHyp;
  ops->nvinv          = N_VInv_ParHyp;
  ops->nvaddconst     = N_VAddConst_ParHyp;
  ops->nvdotprod      = N_VDotProd_ParHyp;
  ops->nvmaxnorm      = N_VMaxNorm_ParHyp;
  ops->nvwrmsnormmask = N_VWrmsNormMask_ParHyp;
  ops->nvwrmsnorm     = N_VWrmsNorm_ParHyp;
  ops->nvmin          = N_VMin_ParHyp;
  ops->nvwl2norm      = N_VWL2Norm_ParHyp;
  ops->nvl1norm       = N_VL1Norm_ParHyp;
  ops->nvcompare      = N_VCompare_ParHyp;
  ops->nvinvtest      = N_VInvTest_ParHyp;
  ops->nvconstrmask   = N_VConstrMask_ParHyp;
  ops->nvminquotient  = N_VMinQuotient_ParHyp;

  /* fused vector operations (optional, NULL means disabled by default) */
  ops->nvlinearcombination = NULL;
  ops->nvscaleaddmulti     = NULL;
  ops->nvdotprodmulti      = NULL;

  /* vector array operations (optional, NULL means disabled by default) */
  ops->nvlinearsumvectorarray         = NULL;
  ops->nvscalevectorarray             = NULL;
  ops->nvconstvectorarray             = NULL;
  ops->nvwrmsnormvectorarray          = NULL;
  ops->nvwrmsnormmaskvectorarray      = NULL;
  ops->nvscaleaddmultivectorarray     = NULL;
  ops->nvlinearcombinationvectorarray = NULL;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_ParHyp) malloc(sizeof(struct _N_VectorContent_ParHyp));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_parvector = SUNFALSE;
  content->x             = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}


/* ----------------------------------------------------------------
 * Function to create a parhyp N_Vector wrapper around user
 * supplie HYPRE vector.
 */

N_Vector N_VMake_ParHyp(HYPRE_ParVector x)
{
  N_Vector v;
  MPI_Comm comm = hypre_ParVectorComm(x);
  HYPRE_Int global_length = hypre_ParVectorGlobalSize(x);
  HYPRE_Int local_begin = hypre_ParVectorFirstIndex(x);
  HYPRE_Int local_end = hypre_ParVectorLastIndex(x);
  HYPRE_Int local_length = local_end - local_begin + 1;

  v = NULL;
  v = N_VNewEmpty_ParHyp(comm, local_length, global_length);
  if (v == NULL)
    return(NULL);

  NV_OWN_PARVEC_PH(v)   = SUNFALSE;
  NV_HYPRE_PARVEC_PH(v) = x;

  return(v);
}


/* ----------------------------------------------------------------
 * Function to create an array of new parhyp vectors.
 */

N_Vector *N_VCloneVectorArray_ParHyp(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_ParHyp(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_ParHyp(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to create an array of new parhyp vector wrappers
 * without uderlying HYPRE vectors.
 */

N_Vector *N_VCloneVectorArrayEmpty_ParHyp(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_ParHyp(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_ParHyp(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_ParHyp
 */

void N_VDestroyVectorArray_ParHyp(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++)
    N_VDestroy_ParHyp(vs[j]);

  free(vs);
  vs = NULL;

  return;
}


/* ----------------------------------------------------------------
 * Extract HYPRE vector
 */

HYPRE_ParVector N_VGetVector_ParHyp(N_Vector v)
{
  return NV_HYPRE_PARVEC_PH(v);
}

/* ----------------------------------------------------------------
 * Function to print a parhyp vector.
 * TODO: Consider using a HYPRE function for this.
 */

void N_VPrint_ParHyp(N_Vector x)
{
  N_VPrintFile_ParHyp(x, stdout);
}

/* ----------------------------------------------------------------
 * Function to print a parhyp vector.
 * TODO: Consider using a HYPRE function for this.
 */

void N_VPrintFile_ParHyp(N_Vector x, FILE *outfile)
{
  sunindextype i, N;
  realtype *xd;

  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);

  for (i = 0; i < N; i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    STAN_SUNDIALS_FPRINTF(outfile, "%Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    STAN_SUNDIALS_FPRINTF(outfile, "%g\n", xd[i]);
#else
    STAN_SUNDIALS_FPRINTF(outfile, "%g\n", xd[i]);
#endif
  }
  STAN_SUNDIALS_FPRINTF(outfile, "\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_ParHyp(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_ParHyp content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Added variables for hypre_parhyp intialization */
  int nprocs, myid;
  MPI_Comm_size(NV_COMM_PH(w), &nprocs);
  MPI_Comm_rank(NV_COMM_PH(w), &myid);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = w->ops->nvgetvectorid;
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;

  /* standard vector operations */
  ops->nvlinearsum    = w->ops->nvlinearsum;
  ops->nvconst        = w->ops->nvconst;
  ops->nvprod         = w->ops->nvprod;
  ops->nvdiv          = w->ops->nvdiv;
  ops->nvscale        = w->ops->nvscale;
  ops->nvabs          = w->ops->nvabs;
  ops->nvinv          = w->ops->nvinv;
  ops->nvaddconst     = w->ops->nvaddconst;
  ops->nvdotprod      = w->ops->nvdotprod;
  ops->nvmaxnorm      = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm     = w->ops->nvwrmsnorm;
  ops->nvmin          = w->ops->nvmin;
  ops->nvwl2norm      = w->ops->nvwl2norm;
  ops->nvl1norm       = w->ops->nvl1norm;
  ops->nvcompare      = w->ops->nvcompare;
  ops->nvinvtest      = w->ops->nvinvtest;
  ops->nvconstrmask   = w->ops->nvconstrmask;
  ops->nvminquotient  = w->ops->nvminquotient;

  /* fused vector operations */
  ops->nvlinearcombination = w->ops->nvlinearcombination;
  ops->nvscaleaddmulti     = w->ops->nvscaleaddmulti;
  ops->nvdotprodmulti      = w->ops->nvdotprodmulti;

  /* vector array operations */
  ops->nvlinearsumvectorarray         = w->ops->nvlinearsumvectorarray;
  ops->nvscalevectorarray             = w->ops->nvscalevectorarray;
  ops->nvconstvectorarray             = w->ops->nvconstvectorarray;
  ops->nvwrmsnormvectorarray          = w->ops->nvwrmsnormvectorarray;
  ops->nvwrmsnormmaskvectorarray      = w->ops->nvwrmsnormmaskvectorarray;
  ops->nvscaleaddmultivectorarray     = w->ops->nvscaleaddmultivectorarray;
  ops->nvlinearcombinationvectorarray = w->ops->nvlinearcombinationvectorarray;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_ParHyp) malloc(sizeof(struct _N_VectorContent_ParHyp));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  /* Attach lengths and communicator */
  content->local_length  = NV_LOCLENGTH_PH(w);
  content->global_length = NV_GLOBLENGTH_PH(w);
  content->comm          = NV_COMM_PH(w);
  content->own_parvector = SUNFALSE;
  content->x             = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/*
 * Clone HYPRE vector wrapper.
 *
 */
N_Vector N_VClone_ParHyp(N_Vector w)
{
  N_Vector v;
  HYPRE_ParVector vx;
  const HYPRE_ParVector wx = NV_HYPRE_PARVEC_PH(w);

  v = NULL;
  v = N_VCloneEmpty_ParHyp(w);
  if (v==NULL)
    return(NULL);

  vx = hypre_ParVectorCreate(wx->comm, wx->global_size, wx->partitioning);
  hypre_ParVectorInitialize(vx);

  hypre_ParVectorSetPartitioningOwner(vx, 0);
  hypre_ParVectorSetDataOwner(vx, 1);
  hypre_SeqVectorSetDataOwner(hypre_ParVectorLocalVector(vx), 1);

  NV_HYPRE_PARVEC_PH(v) = vx;
  NV_OWN_PARVEC_PH(v) = SUNTRUE;

  return(v);
}

void N_VDestroy_ParHyp(N_Vector v)
{
  if ((NV_OWN_PARVEC_PH(v) == SUNTRUE)) {
    hypre_ParVectorDestroy(NV_HYPRE_PARVEC_PH(v));
  }

  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}


void N_VSpace_ParHyp(N_Vector v, sunindextype *lrw, sunindextype *liw)
{
  MPI_Comm comm;
  int npes;

  comm = NV_COMM_PH(v);
  MPI_Comm_size(comm, &npes);

  *lrw = NV_GLOBLENGTH_PH(v);
  *liw = 2*npes;

  return;
}


/*
 * This function is disabled in ParHyp implementation and returns NULL.
 * The user should extract HYPRE vector using N_VGetVector_ParHyp and
 * then use HYPRE functions to get pointer to raw data of the local HYPRE
 * vector.
 */
realtype *N_VGetArrayPointer_ParHyp(N_Vector v)
{
  return NULL; /* ((realtype *) NV_DATA_PH(v)); */
}


/*
 * This method is not implemented for HYPRE vector wrapper.
 * TODO: Put error handler in the function body.
 */
void N_VSetArrayPointer_ParHyp(realtype *v_data, N_Vector v)
{
  /* Not implemented for Hypre vector */
}

/*
 * Computes z[i] = a*x[i] + b*y[i]
 *
 */
void N_VLinearSum_ParHyp(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  booleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    HYPRE_Complex   alpha=a;
    HYPRE_ParVectorAxpy( alpha, (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(x),
                                (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(y));
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    HYPRE_Complex   beta=b;
    HYPRE_ParVectorAxpy( beta, (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(y),
                               (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(x));
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum_ParHyp(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_ParHyp(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_ParHyp(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_ParHyp(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum_ParHyp(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff_ParHyp(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+(b*yd[i]);

  return;
}

void N_VConst_ParHyp(realtype c, N_Vector z)
{
  HYPRE_Complex value = c;
  HYPRE_ParVectorSetConstantValues( (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(z), value);
  return;
}

/* ----------------------------------------------------------------------------
 * Compute componentwise product z[i] = x[i]*y[i]
 */

void N_VProd_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]*yd[i];

  return;
}


/* ----------------------------------------------------------------------------
 * Compute componentwise division z[i] = x[i]/y[i]
 */

void N_VDiv_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]/yd[i];

  return;
}


void N_VScale_ParHyp(realtype c, N_Vector x, N_Vector z)
{
  HYPRE_Complex value = c;

  if (x != z) {
     HYPRE_ParVectorCopy((HYPRE_ParVector) NV_HYPRE_PARVEC_PH(x), (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(z));
  }
  HYPRE_ParVectorScale(value, (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(z));

  return;
}


void N_VAbs_ParHyp(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = SUNRabs(xd[i]);

  return;
}

void N_VInv_ParHyp(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = ONE/xd[i];

  return;
}

void N_VAddConst_ParHyp(N_Vector x, realtype b, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
     zd[i] = xd[i] + b;

  return;
}

realtype N_VDotProd_ParHyp(N_Vector x, N_Vector y)
{

  HYPRE_Real gsum;
  HYPRE_ParVectorInnerProd( (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(x),
                            (HYPRE_ParVector) NV_HYPRE_PARVEC_PH(y), &gsum);

  return(gsum);
}

realtype N_VMaxNorm_ParHyp(N_Vector x)
{
  sunindextype i, N;
  realtype max, *xd, gmax;
  MPI_Comm comm;

  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  comm = NV_COMM_PH(x);

  max = ZERO;

  for (i = 0; i < N; i++) {
    if (SUNRabs(xd[i]) > max) max = SUNRabs(xd[i]);
  }

  gmax = SUNMPI_Allreduce_scalar(max, 2, comm);

  return(gmax);
}

realtype N_VWrmsNorm_ParHyp(N_Vector x, N_Vector w)
{
  sunindextype i, N, N_global;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  N        = NV_LOCLENGTH_PH(x);
  N_global = NV_GLOBLENGTH_PH(x);
  xd       = NV_DATA_PH(x);
  wd       = NV_DATA_PH(w);
  comm     = NV_COMM_PH(x);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(SUNRsqrt(gsum/N_global));
}

realtype N_VWrmsNormMask_ParHyp(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N, N_global;
  realtype sum, prodi, *xd, *wd, *idd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = idd = NULL;

  N        = NV_LOCLENGTH_PH(x);
  N_global = NV_GLOBLENGTH_PH(x);
  xd       = NV_DATA_PH(x);
  wd       = NV_DATA_PH(w);
  idd      = NV_DATA_PH(id);
  comm = NV_COMM_PH(x);

  for (i = 0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i]*wd[i];
      sum += SUNSQR(prodi);
    }
  }

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(SUNRsqrt(gsum/N_global));
}

realtype N_VMin_ParHyp(N_Vector x)
{
  sunindextype i, N;
  realtype min, *xd, gmin;
  MPI_Comm comm;

  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  comm = NV_COMM_PH(x);

  min = BIG_REAL;

  if (N > 0) {

    xd = NV_DATA_PH(x);

    min = xd[0];

    for (i = 1; i < N; i++) {
      if (xd[i] < min)
        min = xd[i];
    }

  }

  gmin = SUNMPI_Allreduce_scalar(min, 3, comm);

  return(gmin);
}

realtype N_VWL2Norm_ParHyp(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  realtype sum, prodi, *xd, *wd, gsum;
  MPI_Comm comm;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  wd = NV_DATA_PH(w);
  comm = NV_COMM_PH(x);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SUNSQR(prodi);
  }

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(SUNRsqrt(gsum));
}

realtype N_VL1Norm_ParHyp(N_Vector x)
{
  sunindextype i, N;
  realtype sum, gsum, *xd;
  MPI_Comm comm;

  sum = ZERO;
  xd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  comm = NV_COMM_PH(x);

  for (i = 0; i<N; i++)
    sum += SUNRabs(xd[i]);

  gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return(gsum);
}

void N_VCompare_ParHyp(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++) {
    zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO;
  }

  return;
}

booleantype N_VInvTest_ParHyp(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *zd, val, gval;
  MPI_Comm comm;

  xd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  zd = NV_DATA_PH(z);
  comm = NV_COMM_PH(x);

  val = ONE;
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO)
      val = ZERO;
    else
      zd[i] = ONE/xd[i];
  }

  gval = SUNMPI_Allreduce_scalar(val, 3, comm);

  if (gval == ZERO)
    return(SUNFALSE);
  else
    return(SUNTRUE);
}

booleantype N_VConstrMask_ParHyp(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  realtype temp;
  realtype *cd, *xd, *md;
  booleantype test;
  MPI_Comm comm;

  cd = xd = md = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  cd = NV_DATA_PH(c);
  md = NV_DATA_PH(m);
  comm = NV_COMM_PH(x);

  temp = ZERO;

  for (i = 0; i < N; i++) {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO)
      continue;

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cd[i]) > ONEPT5 && xd[i]*cd[i] <= ZERO) ||
           (SUNRabs(cd[i]) > HALF   && xd[i]*cd[i] <  ZERO);
    if (test) {
      temp = md[i] = ONE;
    }
  }

  /* Find max temp across all MPI ranks */
  temp = SUNMPI_Allreduce_scalar(temp, 2, comm);

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

realtype N_VMinQuotient_ParHyp(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce;
  sunindextype i, N;
  realtype *nd, *dd, min;
  MPI_Comm comm;

  nd = dd = NULL;

  N  = NV_LOCLENGTH_PH(num);
  nd = NV_DATA_PH(num);
  dd = NV_DATA_PH(denom);
  comm = NV_COMM_PH(num);

  notEvenOnce = SUNTRUE;
  min = BIG_REAL;

  for (i = 0; i < N; i++) {
    if (dd[i] == ZERO) continue;
    else {
      if (!notEvenOnce) min = SUNMIN(min, nd[i]/dd[i]);
      else {
        min = nd[i]/dd[i];
        notEvenOnce = SUNFALSE;
      }
    }
  }

  return(SUNMPI_Allreduce_scalar(min, 3, comm));
}


/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */


int N_VLinearCombination_ParHyp(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_ParHyp(c[0], X[0], z);
    return(0);
  }

  /* should have called N_VLinearSum */
  if (nvec == 2) {
    N_VLinearSum_ParHyp(c[0], X[0], c[1], X[1], z);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_LOCLENGTH_PH(z);
  zd = NV_DATA_PH(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE)) {
    for (i=1; i<nvec; i++) {
      xd = NV_DATA_PH(X[i]);
      for (j=0; j<N; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z) {
    for (j=0; j<N; j++) {
      zd[j] *= c[0];
    }
    for (i=1; i<nvec; i++) {
      xd = NV_DATA_PH(X[i]);
      for (j=0; j<N; j++) {
        zd[j] += c[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = NV_DATA_PH(X[0]);
  for (j=0; j<N; j++) {
    zd[j] = c[0] * xd[j];
  }
  for (i=1; i<nvec; i++) {
    xd = NV_DATA_PH(X[i]);
    for (j=0; j<N; j++) {
      zd[j] += c[i] * xd[j];
    }
  }
  return(0);
}


int N_VScaleAddMulti_ParHyp(int nvec, realtype* a, N_Vector x, N_Vector* Y,
                             N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_ParHyp(a[0], x, ONE, Y[0], Z[0]);
    return(0);
  }

  /* get vector length and data array */
  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      yd = NV_DATA_PH(Y[i]);
      for (j=0; j<N; j++) {
        yd[j] += a[i] * xd[j];
      }
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    yd = NV_DATA_PH(Y[i]);
    zd = NV_DATA_PH(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = a[i] * xd[j] + yd[j];
    }
  }
  return(0);
}


int N_VDotProdMulti_ParHyp(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VDotProd */
  if (nvec == 1) {
    dotprods[0] = N_VDotProd_ParHyp(x, Y[0]);
    return(0);
  }

  /* get vector length, data array, and communicator */
  N    = NV_LOCLENGTH_PH(x);
  xd   = NV_DATA_PH(x);
  comm = NV_COMM_PH(x);

  /* compute multiple dot products */
  for (i=0; i<nvec; i++) {
    yd = NV_DATA_PH(Y[i]);
    dotprods[i] = ZERO;
    for (j=0; j<N; j++) {
      dotprods[i] += xd[j] * yd[j];
    }
  }
  SUNMPI_Allreduce(dotprods, nvec, 1, comm);

  return(0);
}


/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */


int N_VLinearSumVectorArray_ParHyp(int nvec,
                                   realtype a, N_Vector* X,
                                   realtype b, N_Vector* Y,
                                   N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VLinearSum */
  if (nvec == 1) {
    N_VLinearSum_ParHyp(a, X[0], b, Y[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LOCLENGTH_PH(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_PH(X[i]);
    yd = NV_DATA_PH(Y[i]);
    zd = NV_DATA_PH(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = a * xd[j] + b * yd[j];
    }
  }

  return(0);
}


int N_VScaleVectorArray_ParHyp(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    xd=NULL;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VScale */
  if (nvec == 1) {
    N_VScale_ParHyp(c[0], X[0], Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LOCLENGTH_PH(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_PH(X[i]);
      for (j=0; j<N; j++) {
        xd[j] *= c[i];
      }
    }
    return(0);
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_PH(X[i]);
    zd = NV_DATA_PH(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = c[i] * xd[j];
    }
  }

  return(0);
}


int N_VConstVectorArray_ParHyp(int nvec, realtype c, N_Vector* Z)
{
  int          i;
  sunindextype j, N;
  realtype*    zd=NULL;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VConst */
  if (nvec == 1) {
    N_VConst_ParHyp(c, Z[0]);
    return(0);
  }

  /* get vector length */
  N = NV_LOCLENGTH_PH(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i=0; i<nvec; i++) {
    zd = NV_DATA_PH(Z[i]);
    for (j=0; j<N; j++) {
      zd[j] = c;
    }
  }

  return(0);
}


int N_VWrmsNormVectorArray_ParHyp(int nvec, N_Vector* X, N_Vector* W, realtype* nrm)
{
  int          i;
  sunindextype j, Nl, Ng;
  realtype*    wd=NULL;
  realtype*    xd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNorm_ParHyp(X[0], W[0]);
    return(0);
  }

  /* get vector lengths and communicator */
  Nl   = NV_LOCLENGTH_PH(X[0]);
  Ng   = NV_GLOBLENGTH_PH(X[0]);
  comm = NV_COMM_PH(X[0]);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_PH(X[i]);
    wd = NV_DATA_PH(W[i]);
    nrm[i] = ZERO;
    for (j=0; j<Nl; j++) {
      nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
  }
  SUNMPI_Allreduce(nrm, nvec, 1, comm);

  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/Ng);

  return(0);
}


int N_VWrmsNormMaskVectorArray_ParHyp(int nvec, N_Vector* X, N_Vector* W,
                                        N_Vector id, realtype* nrm)
{
  int          i;
  sunindextype j, Nl, Ng;
  realtype*    wd=NULL;
  realtype*    xd=NULL;
  realtype*    idd=NULL;
  MPI_Comm     comm;

  /* invalid number of vectors */
  if (nvec < 1) return(-1);

  /* should have called N_VWrmsNorm */
  if (nvec == 1) {
    nrm[0] = N_VWrmsNormMask_ParHyp(X[0], W[0], id);
    return(0);
  }

  /* get vector lengths, communicator, and mask data */
  Nl   = NV_LOCLENGTH_PH(X[0]);
  Ng   = NV_GLOBLENGTH_PH(X[0]);
  comm = NV_COMM_PH(X[0]);
  idd  = NV_DATA_PH(id);

  /* compute the WRMS norm for each vector in the vector array */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_PH(X[i]);
    wd = NV_DATA_PH(W[i]);
    nrm[i] = ZERO;
    for (j=0; j<Nl; j++) {
      if (idd[j] > ZERO)
        nrm[i] += SUNSQR(xd[j] * wd[j]);
    }
  }
  SUNMPI_Allreduce(nrm, nvec, 1, comm);

  for (i=0; i<nvec; i++)
    nrm[i] = SUNRsqrt(nrm[i]/Ng);

  return(0);
}


int N_VScaleAddMultiVectorArray_ParHyp(int nvec, int nsum, realtype* a,
                                          N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  int          i, j;
  sunindextype k, N;
  realtype*    xd=NULL;
  realtype*    yd=NULL;
  realtype*    zd=NULL;

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
      N_VLinearSum_ParHyp(a[0], X[0], ONE, Y[0][0], Z[0][0]);
      return(0);
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector *) malloc(nsum * sizeof(N_Vector));
    ZZ = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (j=0; j<nsum; j++) {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    retval = N_VScaleAddMulti_ParHyp(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);
    return(retval);
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1) {
    retval = N_VLinearSumVectorArray_ParHyp(nvec, a[0], X, ONE, Y[0], Z[0]);
    return(retval);
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N  = NV_LOCLENGTH_PH(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z) {
    for (i=0; i<nvec; i++) {
      xd = NV_DATA_PH(X[i]);
      for (j=0; j<nsum; j++){
        yd = NV_DATA_PH(Y[j][i]);
        for (k=0; k<N; k++) {
          yd[k] += a[j] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i=0; i<nvec; i++) {
    xd = NV_DATA_PH(X[i]);
    for (j=0; j<nsum; j++) {
      yd = NV_DATA_PH(Y[j][i]);
      zd = NV_DATA_PH(Z[j][i]);
      for (k=0; k<N; k++) {
        zd[k] = a[j] * xd[k] + yd[k];
      }
    }
  }
  return(0);
}


int N_VLinearCombinationVectorArray_ParHyp(int nvec, int nsum,
                                             realtype* c,
                                             N_Vector** X,
                                             N_Vector* Z)
{
  int          i; /* vector arrays index in summation [0,nsum) */
  int          j; /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  realtype*    zd=NULL;
  realtype*    xd=NULL;

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
      N_VScale_ParHyp(c[0], X[0][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearSum */
    if (nsum == 2) {
      N_VLinearSum_ParHyp(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return(0);
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector *) malloc(nsum * sizeof(N_Vector));

    for (i=0; i<nsum; i++) {
      Y[i] = X[i][0];
    }

    N_VLinearCombination_ParHyp(nsum, c, Y, Z[0]);

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

    N_VScaleVectorArray_ParHyp(nvec, ctmp, X[0], Z);

    free(ctmp);
    return(0);
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2) {
    N_VLinearSumVectorArray_ParHyp(nvec, c[0], X[0], c[1], X[1], Z);
    return(0);
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = NV_LOCLENGTH_PH(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && (c[0] == ONE)) {
    for (j=0; j<nvec; j++) {
      zd = NV_DATA_PH(Z[j]);
      for (i=1; i<nsum; i++) {
        xd = NV_DATA_PH(X[i][j]);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z) {
    for (j=0; j<nvec; j++) {
      zd = NV_DATA_PH(Z[j]);
      for (k=0; k<N; k++) {
        zd[k] *= c[0];
      }
      for (i=1; i<nsum; i++) {
        xd = NV_DATA_PH(X[i][j]);
        for (k=0; k<N; k++) {
          zd[k] += c[i] * xd[k];
        }
      }
    }
    return(0);
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j=0; j<nvec; j++) {
    xd = NV_DATA_PH(X[0][j]);
    zd = NV_DATA_PH(Z[j]);
    for (k=0; k<N; k++) {
      zd[k] = c[0] * xd[k];
    }
    for (i=1; i<nsum; i++) {
      xd = NV_DATA_PH(X[i][j]);
      for (k=0; k<N; k++) {
        zd[k] += c[i] * xd[k];
      }
    }
  }
  return(0);
}


/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static void VSum_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]+yd[i];

  return;
}

static void VDiff_ParHyp(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]-yd[i];

  return;
}


static void VScaleSum_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]+yd[i]);

  return;
}

static void VScaleDiff_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = c*(xd[i]-yd[i]);

  return;
}

static void VLin1_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])+yd[i];

  return;
}

static void VLin2_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  realtype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  yd = NV_DATA_PH(y);
  zd = NV_DATA_PH(z);

  for (i = 0; i < N; i++)
    zd[i] = (a*xd[i])-yd[i];

  return;
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_ParHyp;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_ParHyp;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_ParHyp;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_ParHyp;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_ParHyp;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_ParHyp;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_ParHyp;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_ParHyp;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_ParHyp;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_ParHyp;
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
  }

  /* return success */
  return(0);
}


int N_VEnableLinearCombination_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_ParHyp;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_ParHyp;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_ParHyp;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_ParHyp;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_ParHyp;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_ParHyp;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_ParHyp;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_ParHyp;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_ParHyp;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_ParHyp(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_ParHyp;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}
