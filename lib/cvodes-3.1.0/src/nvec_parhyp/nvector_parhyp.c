/* -----------------------------------------------------------------
 * Programmer(s): Jean M. Sexton @ SMU
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * Based on work by Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                  and Aaron Collier @ LLNL
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
 * This is the implementation file for a parhyp MPI implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_parhyp.h>
#include <sundials/sundials_math.h>

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
 *     to be a pointer to hypre_ParVector of vector v. The assignment
 *     NV_HYPRE_PARVEC_PH(v) = parhyp_v sets pointer to 
 *     hypre_ParVector of vector v to be parhyp_v.
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

/* Reduction operations add/max/min over the processor group */
static realtype VAllReduce_ParHyp(realtype d, int op, MPI_Comm comm);
/* z=x */
/* static void VCopy_ParHyp(N_Vector x, N_Vector z); */
/* z=x+y */
static void VSum_ParHyp(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_ParHyp(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
/* static void VNeg_ParHyp(N_Vector x, N_Vector z); */
/* z=c(x+y) */
static void VScaleSum_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_ParHyp(realtype c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_ParHyp(realtype a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
/* static void Vaxpy_ParHyp(realtype a, N_Vector x, N_Vector y); */
/* x <- ax */
/* static void VScaleBy_ParHyp(realtype a, N_Vector x); */

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
  ops->nvlinearsum       = N_VLinearSum_ParHyp;
  ops->nvconst           = N_VConst_ParHyp;
  ops->nvprod            = N_VProd_ParHyp;
  ops->nvdiv             = N_VDiv_ParHyp;
  ops->nvscale           = N_VScale_ParHyp;
  ops->nvabs             = N_VAbs_ParHyp;
  ops->nvinv             = N_VInv_ParHyp;
  ops->nvaddconst        = N_VAddConst_ParHyp;
  ops->nvdotprod         = N_VDotProd_ParHyp;
  ops->nvmaxnorm         = N_VMaxNorm_ParHyp;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_ParHyp;
  ops->nvwrmsnorm        = N_VWrmsNorm_ParHyp;
  ops->nvmin             = N_VMin_ParHyp;
  ops->nvwl2norm         = N_VWL2Norm_ParHyp;
  ops->nvl1norm          = N_VL1Norm_ParHyp;
  ops->nvcompare         = N_VCompare_ParHyp;
  ops->nvinvtest         = N_VInvTest_ParHyp;
  ops->nvconstrmask      = N_VConstrMask_ParHyp;
  ops->nvminquotient     = N_VMinQuotient_ParHyp;

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

N_Vector N_VMake_ParHyp(hypre_ParVector *x)
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

hypre_ParVector* N_VGetVector_ParHyp(N_Vector v)
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
    fprintf(outfile, "%Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%g\n", xd[i]);
#else
    fprintf(outfile, "%g\n", xd[i]);
#endif
  }
  fprintf(outfile, "\n");

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
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;  
  ops->nvprod            = w->ops->nvprod;   
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale; 
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;    
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

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
  hypre_ParVector *vx;
  const hypre_ParVector *wx = NV_HYPRE_PARVEC_PH(w);
  
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
   
  gmax = VAllReduce_ParHyp(max, 2, comm);

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

  gsum = VAllReduce_ParHyp(sum, 1, comm);

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

  gsum = VAllReduce_ParHyp(sum, 1, comm);

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

  gmin = VAllReduce_ParHyp(min, 3, comm);

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

  gsum = VAllReduce_ParHyp(sum, 1, comm);

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

  gsum = VAllReduce_ParHyp(sum, 1, comm);

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

  gval = VAllReduce_ParHyp(val, 3, comm);

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
  MPI_Comm comm;

  cd = xd = md = NULL;

  N  = NV_LOCLENGTH_PH(x);
  xd = NV_DATA_PH(x);
  cd = NV_DATA_PH(c);
  md = NV_DATA_PH(m);
  comm = NV_COMM_PH(x);

  temp = ONE;

  for (i = 0; i < N; i++) {
    md[i] = ZERO;
    if (cd[i] == ZERO) continue;
    if (cd[i] > ONEPT5 || cd[i] < -ONEPT5) {
      if (xd[i]*cd[i] <= ZERO) { 
        temp = ZERO; 
        md[i] = ONE;
      }
      continue;
    }
    if (cd[i] > HALF || cd[i] < -HALF) {
      if (xd[i]*cd[i] < ZERO ) { 
        temp = ZERO; 
        md[i] = ONE;
      }
    }
  }

  temp = VAllReduce_ParHyp(temp, 3, comm);

  if (temp == ONE) 
    return(SUNTRUE);
  else 
    return(SUNFALSE);
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

  return(VAllReduce_ParHyp(min, 3, comm));
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static realtype VAllReduce_ParHyp(realtype d, int op, MPI_Comm comm)
{
  /* 
   * This function does a global reduction.  The operation is
   *   sum if op = 1,
   *   max if op = 2,
   *   min if op = 3.
   * The operation is over all processors in the communicator 
   */

  realtype out;

  switch (op) {
   case 1: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }

  return(out);
}


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

