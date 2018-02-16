/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * 
 * Based on N_Vector_Parallel by Scott D. Cohen, Alan C. Hindmarsh, 
 * Radu Serban, and Aaron Collier @ LLNL
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
 *     the content structure of a parallel N_Vector.
 *
 *     NV_PVEC_PTC(v) returns pointer to the PETSc vector. 
 *
 *     The assignment v_llen = NV_LOCLENGTH_PTC(v) sets v_llen to
 *     be the length of the local part of the vector v. The call
 *     NV_LOCLENGTH_PTC(v) = llen_v sets the local length
 *     of v to be llen_v.
 *
 *     The assignment v_glen = NV_GLOBLENGTH_PTC(v) sets v_glen to
 *     be the global length of the vector v. The call
 *     NV_GLOBLENGTH_PTC(v) = glen_v sets the global length of v to
 *     be glen_v.
 *
 *     The assignment v_comm = NV_COMM_PTC(v) sets v_comm to be the
 *     MPI communicator of the vector v. The assignment
 *     NV_COMM_PTC(v) = comm_v sets the MPI communicator of v to be
 *     comm_v.
 *
 * -----------------------------------------------------------------
 */

#define NV_CONTENT_PTC(v)    ( (N_VectorContent_Petsc)(v->content) )

#define NV_LOCLENGTH_PTC(v)  ( NV_CONTENT_PTC(v)->local_length )

#define NV_GLOBLENGTH_PTC(v) ( NV_CONTENT_PTC(v)->global_length )

#define NV_OWN_DATA_PTC(v)   ( NV_CONTENT_PTC(v)->own_data )

#define NV_PVEC_PTC(v)       ( NV_CONTENT_PTC(v)->pvec )

#define NV_COMM_PTC(v)       ( NV_CONTENT_PTC(v)->comm )


/* Private function prototypes */

/* Reduction operations add/max/min over the processor group */
static realtype VAllReduce_Petsc(realtype d, int op, MPI_Comm comm);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation 
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Petsc(N_Vector v)
{
  return SUNDIALS_NVEC_PETSC;
}


/* ----------------------------------------------------------------
 * Function to create a new parallel vector with empty data array
 */

N_Vector N_VNewEmpty_Petsc(MPI_Comm comm, 
                           sunindextype local_length,
                           sunindextype global_length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Petsc content;
  sunindextype n, Nsum;
  PetscErrorCode ierr;

  /* Compute global length as sum of local lengths */
  n = local_length;
  ierr = MPI_Allreduce(&n, &Nsum, 1, PVEC_INTEGER_MPI_TYPE, MPI_SUM, comm);
  CHKERRABORT(comm,ierr);
  if (Nsum != global_length) {
    fprintf(stderr, BAD_N);
    return(NULL);
  } 

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = N_VGetVectorID_Petsc;
  ops->nvclone           = N_VClone_Petsc;
  ops->nvcloneempty      = N_VCloneEmpty_Petsc;
  ops->nvdestroy         = N_VDestroy_Petsc;
  ops->nvspace           = N_VSpace_Petsc;
  ops->nvgetarraypointer = N_VGetArrayPointer_Petsc;
  ops->nvsetarraypointer = N_VSetArrayPointer_Petsc;
  ops->nvlinearsum       = N_VLinearSum_Petsc;
  ops->nvconst           = N_VConst_Petsc;
  ops->nvprod            = N_VProd_Petsc;
  ops->nvdiv             = N_VDiv_Petsc;
  ops->nvscale           = N_VScale_Petsc;
  ops->nvabs             = N_VAbs_Petsc;
  ops->nvinv             = N_VInv_Petsc;
  ops->nvaddconst        = N_VAddConst_Petsc;
  ops->nvdotprod         = N_VDotProd_Petsc;
  ops->nvmaxnorm         = N_VMaxNorm_Petsc;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Petsc;
  ops->nvwrmsnorm        = N_VWrmsNorm_Petsc;
  ops->nvmin             = N_VMin_Petsc;
  ops->nvwl2norm         = N_VWL2Norm_Petsc;
  ops->nvl1norm          = N_VL1Norm_Petsc;
  ops->nvcompare         = N_VCompare_Petsc;
  ops->nvinvtest         = N_VInvTest_Petsc;
  ops->nvconstrmask      = N_VConstrMask_Petsc;
  ops->nvminquotient     = N_VMinQuotient_Petsc;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Petsc) malloc(sizeof(struct _N_VectorContent_Petsc));
  if (content == NULL) { 
    free(ops); 
    free(v); 
    return(NULL); 
  }

  /* Attach lengths and communicator */
  content->local_length  = local_length;
  content->global_length = global_length;
  content->comm          = comm;
  content->own_data      = SUNFALSE;
  content->pvec          = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}



/* ---------------------------------------------------------------- 
 * Function to create a parallel N_Vector with user data component
 * This function is NOT implemented for PETSc wrapper!
 */

N_Vector N_VMake_Petsc(Vec *pvec)
{
  N_Vector v = NULL;
  MPI_Comm comm;
  PetscInt local_length;
  PetscInt global_length;

  VecGetLocalSize(*pvec, &local_length);
  VecGetSize(*pvec, &global_length);
  PetscObjectGetComm((PetscObject) (*pvec), &comm);
  
  v = N_VNewEmpty_Petsc(comm, local_length, global_length);
  if (v == NULL) 
     return(NULL);

  /* Attach data */
  NV_OWN_DATA_PTC(v) = SUNFALSE;
  NV_PVEC_PTC(v)     = pvec;

  return(v);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new parallel vectors. 
 */

N_Vector *N_VCloneVectorArray_Petsc(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_Petsc(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Petsc(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ---------------------------------------------------------------- 
 * Function to create an array of new parallel vectors with empty
 * (NULL) data array.
 */

N_Vector *N_VCloneVectorArrayEmpty_Petsc(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_Petsc(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Petsc(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Petsc
 */

void N_VDestroyVectorArray_Petsc(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Petsc(vs[j]);

  free(vs); 
  vs = NULL;

  return;
}

/* ---------------------------------------------------------------- 
 * Function to extract PETSc vector 
 */

Vec *N_VGetVector_Petsc(N_Vector v)
{
  return NV_PVEC_PTC(v);
}

/* ---------------------------------------------------------------- 
 * Function to print the global data in a PETSc parallel vector to
 * stdout
 */

void N_VPrint_Petsc(N_Vector x)
{
  Vec *xv = NV_PVEC_PTC(x);
  MPI_Comm comm = NV_COMM_PTC(x);
  
  VecView(*xv, PETSC_VIEWER_STDOUT_(comm));

  return;
}

/* ---------------------------------------------------------------- 
 * Function to print the global data in a PETSc parallel vector to
 * fname
 */

void N_VPrintFile_Petsc(N_Vector x, const char fname[])
{
  Vec *xv = NV_PVEC_PTC(x);
  MPI_Comm comm = NV_COMM_PTC(x);
  PetscViewer viewer;

  PetscViewerASCIIOpen(comm, fname, &viewer);

  VecView(*xv, viewer);

  PetscViewerDestroy(&viewer);

  return;
}

/* ---------------------------------------------------------------- 
 * Function to print the local data in a PETSc parallel vector to
 * outfile
 */

/*
void N_VPrintFileLocal_Petsc(N_Vector x, FILE *outfile)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec *xv = NV_PVEC_PTC(x);
  PetscScalar *xd;

  VecGetArray(*xv, &xd);

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

  VecRestoreArray(*xv, &xd);

  return;
}
*/

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Petsc(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Petsc content;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { 
    free(v); 
    return(NULL); 
  }
  
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
  content = (N_VectorContent_Petsc) malloc(sizeof(struct _N_VectorContent_Petsc));
  if (content == NULL) { 
    free(ops); 
    free(v); 
    return(NULL); 
  }

  /* Attach lengths and communicator */
  content->local_length  = NV_LOCLENGTH_PTC(w);
  content->global_length = NV_GLOBLENGTH_PTC(w);
  content->comm          = NV_COMM_PTC(w);
  content->own_data      = SUNFALSE;
  content->pvec          = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

N_Vector N_VClone_Petsc(N_Vector w)
{
  N_Vector v     = NULL;
  Vec *pvec      = NULL;
  Vec *wvec      = NV_PVEC_PTC(w);
  
  /* PetscErrorCode ierr; */
  
  v = N_VCloneEmpty_Petsc(w);
  if (v == NULL) 
    return(NULL);

  /* Create data */

  /* Allocate empty PETSc vector */
  pvec = (Vec*) malloc(sizeof(Vec));
  if(pvec == NULL) {
    N_VDestroy_Petsc(v); 
    return(NULL);
  }
    
  /* ierr = */ 
  VecDuplicate(*wvec, pvec);
  if(pvec == NULL) {
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
  if (NV_OWN_DATA_PTC(v) == SUNTRUE) {
    VecDestroy((NV_PVEC_PTC(v)));
    NV_PVEC_PTC(v) = NULL;
  }
  
  free(v->content); 
  v->content = NULL;
  free(v->ops); 
  v->ops = NULL;
  free(v); 
  v = NULL;

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

void N_VLinearSum_Petsc(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  Vec *xv = NV_PVEC_PTC(x);
  Vec *yv = NV_PVEC_PTC(y);
  Vec *zv = NV_PVEC_PTC(z);
  
  if (x == y) {
    N_VScale_Petsc(a + b, x, z); /* z <~ ax+bx */
    return;
  }

  if (z == y) {
    if (b == ONE) { 
      VecAXPY(*yv, a, *xv);   /* BLAS usage: axpy  y <- ax+y */
      return;
    }
    VecAXPBY(*yv, a, b, *xv); /* BLAS usage: axpby y <- ax+by */
    return;
  }

  if (z == x) {
    if (a == ONE) { 
      VecAXPY(*xv, b, *yv);   /* BLAS usage: axpy  x <- by+x */
      return;
    }
    VecAXPBY(*xv, b, a, *yv); /* BLAS usage: axpby x <- by+ax */
    return;
  }


  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
  
  VecAXPBYPCZ(*zv, a, b, 0.0, *xv, *yv); // PETSc, probably not optimal 

  return;
}

void N_VConst_Petsc(realtype c, N_Vector z)
{
  Vec *zv      = NV_PVEC_PTC(z);

  VecSet(*zv, c);
  
  return;
}

void N_VProd_Petsc(N_Vector x, N_Vector y, N_Vector z)
{
  Vec *xv = NV_PVEC_PTC(x);
  Vec *yv = NV_PVEC_PTC(y);
  Vec *zv = NV_PVEC_PTC(z);
  
  VecPointwiseMult(*zv, *xv, *yv);
  
  return;
}

void N_VDiv_Petsc(N_Vector x, N_Vector y, N_Vector z)
{
  Vec *xv = NV_PVEC_PTC(x);
  Vec *yv = NV_PVEC_PTC(y);
  Vec *zv = NV_PVEC_PTC(z);

  VecPointwiseDivide(*zv, *xv, *yv); /* z = x/y */

  return;
}

void N_VScale_Petsc(realtype c, N_Vector x, N_Vector z)
{
  Vec *xv = NV_PVEC_PTC(x);
  Vec *zv = NV_PVEC_PTC(z);

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VecScale(*xv, c);
    return;
  }
  
  VecAXPBY(*zv, c, 0.0, *xv); 

  return;
}

void N_VAbs_Petsc(N_Vector x, N_Vector z)
{
  Vec *xv = NV_PVEC_PTC(x);
  Vec *zv = NV_PVEC_PTC(z);

  if(z != x)
    VecCopy(*xv, *zv); /* copy x~>z */
  VecAbs(*zv); 
  
  return;
}

void N_VInv_Petsc(N_Vector x, N_Vector z)
{
  Vec *xv = NV_PVEC_PTC(x);
  Vec *zv = NV_PVEC_PTC(z);

  if(z != x)
    VecCopy(*xv, *zv); /* copy x~>z */
  VecReciprocal(*zv);

  return;
}

void N_VAddConst_Petsc(N_Vector x, realtype b, N_Vector z)
{
  Vec *xv = NV_PVEC_PTC(x);
  Vec *zv = NV_PVEC_PTC(z);

  if(z != x)
    VecCopy(*xv, *zv); /* copy x~>z */
  VecShift(*zv, b);

  return;
}

realtype N_VDotProd_Petsc(N_Vector x, N_Vector y)
{
  Vec *xv = NV_PVEC_PTC(x);
  Vec *yv = NV_PVEC_PTC(y);
  PetscScalar dotprod;
  
  VecDot(*xv, *yv, &dotprod);
  
  return dotprod;
}

realtype N_VMaxNorm_Petsc(N_Vector x)
{
  Vec *xv = NV_PVEC_PTC(x);
  PetscReal norm;
  
  VecNorm(*xv, NORM_INFINITY, &norm);
  
  return norm;
}

realtype N_VWrmsNorm_Petsc(N_Vector x, N_Vector w)
{
  sunindextype i;
  sunindextype N        = NV_LOCLENGTH_PTC(x);
  sunindextype N_global = NV_GLOBLENGTH_PTC(x);
  MPI_Comm comm     = NV_COMM_PTC(x);
  Vec *xv = NV_PVEC_PTC(x);
  Vec *wv = NV_PVEC_PTC(w);
  PetscScalar *xd;
  PetscScalar *wd;
  PetscReal sum = ZERO;
  realtype global_sum;
  
  VecGetArray(*xv, &xd);
  VecGetArray(*wv, &wd);
  for (i = 0; i < N; i++) {
    sum += PetscSqr(PetscAbsScalar(xd[i] * wd[i]));
  }
  VecRestoreArray(*xv, &xd);
  VecRestoreArray(*wv, &wd);
  
  global_sum = VAllReduce_Petsc(sum, 1, comm);
  return (SUNRsqrt(global_sum/N_global)); 
}

realtype N_VWrmsNormMask_Petsc(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i;
  sunindextype N        = NV_LOCLENGTH_PTC(x);
  sunindextype N_global = NV_GLOBLENGTH_PTC(x);
  MPI_Comm comm     = NV_COMM_PTC(x);
  
  Vec *xv = NV_PVEC_PTC(x);
  Vec *wv = NV_PVEC_PTC(w);
  Vec *idv = NV_PVEC_PTC(id);
  PetscScalar *xd;
  PetscScalar *wd;
  PetscScalar *idd;
  PetscReal sum = ZERO;
  realtype global_sum;
  
  VecGetArray(*xv, &xd);
  VecGetArray(*wv, &wd);
  VecGetArray(*idv, &idd);
  for (i = 0; i < N; i++) {
    PetscReal tag = (PetscReal) idd[i];
    if (tag > ZERO) {
      sum += PetscSqr(PetscAbsScalar(xd[i] * wd[i]));
    }
  }
  VecRestoreArray(*xv, &xd);
  VecRestoreArray(*wv, &wd);
  VecRestoreArray(*idv, &idd);

  global_sum = VAllReduce_Petsc(sum, 1, comm);
  return (SUNRsqrt(global_sum/N_global)); 
}

realtype N_VMin_Petsc(N_Vector x)
{
  Vec *xv = NV_PVEC_PTC(x);
  PetscReal minval;
  PetscInt i;
  
  VecMin(*xv, &i, &minval);
  
  return minval;
}

realtype N_VWL2Norm_Petsc(N_Vector x, N_Vector w)
{
  sunindextype i;
  sunindextype N        = NV_LOCLENGTH_PTC(x);
  MPI_Comm comm     = NV_COMM_PTC(x);

  Vec *xv = NV_PVEC_PTC(x);
  Vec *wv = NV_PVEC_PTC(w);
  PetscScalar *xd;
  PetscScalar *wd;
  PetscReal sum = ZERO;
  realtype global_sum;
  
  VecGetArray(*xv, &xd);
  VecGetArray(*wv, &wd);
  for (i = 0; i < N; i++) {
    sum += PetscSqr(PetscAbsScalar(xd[i] * wd[i]));
  }
  VecRestoreArray(*xv, &xd);
  VecRestoreArray(*wv, &wd);

  global_sum = VAllReduce_Petsc(sum, 1, comm);
  return (SUNRsqrt(global_sum)); 
}

realtype N_VL1Norm_Petsc(N_Vector x)
{
  Vec *xv = NV_PVEC_PTC(x);
  PetscReal norm;
  
  VecNorm(*xv, NORM_1, &norm);
  
  return norm;
}

void N_VCompare_Petsc(realtype c, N_Vector x, N_Vector z)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  Vec *xv = NV_PVEC_PTC(x);
  Vec *zv = NV_PVEC_PTC(z);
  PetscReal cpet = c; // <~ realtype should typedef to PETScReal
  PetscScalar *xdata;
  PetscScalar *zdata;

  VecGetArray(*xv, &xdata);
  VecGetArray(*zv, &zdata);
  for (i = 0; i < N; i++) {
    zdata[i] = PetscAbsScalar(xdata[i]) >= cpet ? ONE : ZERO;
  }
  VecRestoreArray(*xv, &xdata);
  VecRestoreArray(*zv, &zdata);

  return;
}

booleantype N_VInvTest_Petsc(N_Vector x, N_Vector z)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  MPI_Comm comm = NV_COMM_PTC(x);
  Vec *xv = NV_PVEC_PTC(x);
  Vec *zv = NV_PVEC_PTC(z);
  PetscScalar *xd;
  PetscScalar *zd;
  PetscReal val = ONE;
  
  VecGetArray(*xv, &xd);
  VecGetArray(*zv, &zd);
  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO) 
      val = ZERO;
    else
      zd[i] = ONE/xd[i];
  }
  VecRestoreArray(*xv, &xd);
  VecRestoreArray(*zv, &zd);

  val = VAllReduce_Petsc(val, 3, comm);

  if (val == ZERO)
    return(SUNFALSE);
  else
    return(SUNTRUE);
}

booleantype N_VConstrMask_Petsc(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i;
  sunindextype N = NV_LOCLENGTH_PTC(x);
  MPI_Comm comm = NV_COMM_PTC(x);
  realtype minval = ONE;
  Vec *xv = NV_PVEC_PTC(x);
  Vec *cv = NV_PVEC_PTC(c);
  Vec *mv = NV_PVEC_PTC(m);
  PetscScalar *xd;
  PetscScalar *cd;
  PetscScalar *md;

  VecGetArray(*xv, &xd);
  VecGetArray(*cv, &cd);
  VecGetArray(*mv, &md);
  for (i = 0; i < N; i++) {
    PetscReal cc = (PetscReal) cd[i]; /* <~ Drop imaginary parts if any. */
    PetscReal xx = (PetscReal) xd[i]; /* <~ Constraints defined on Re{x} */
    md[i] = ZERO;
    if (cc == ZERO) continue;
    if (cc > ONEPT5 || cc < -ONEPT5) {
      if (xx*cc <= ZERO) { minval = ZERO; md[i] = ONE; }
      continue;
    }
    if (cc > HALF || cc < -HALF) {
      if (xx*cc < ZERO ) { minval = ZERO; md[i] = ONE; }
    }
  }
  VecRestoreArray(*xv, &xd);
  VecRestoreArray(*cv, &cd);
  VecRestoreArray(*mv, &md);

  minval = VAllReduce_Petsc(minval, 3, comm);

  if (minval == ONE) 
    return(SUNTRUE);
  else
    return(SUNFALSE);
}

realtype N_VMinQuotient_Petsc(N_Vector num, N_Vector denom)
{
  booleantype notEvenOnce = SUNTRUE;
  sunindextype i; 
  sunindextype N    = NV_LOCLENGTH_PTC(num);
  MPI_Comm comm = NV_COMM_PTC(num);

  Vec *nv = NV_PVEC_PTC(num);
  Vec *dv = NV_PVEC_PTC(denom);
  PetscScalar *nd;
  PetscScalar *dd;
  PetscReal minval = BIG_REAL;

  VecGetArray(*nv, &nd);
  VecGetArray(*dv, &dd);
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
  VecRestoreArray(*nv, &nd);
  VecRestoreArray(*dv, &dd);

  return(VAllReduce_Petsc(minval, 3, comm));
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static realtype VAllReduce_Petsc(realtype d, int op, MPI_Comm comm)
{
  /* 
   * This function does a global reduction.  The operation is
   *   sum if op = 1,
   *   max if op = 2,
   *   min if op = 3.
   * The operation is over all processors in the communicator 
   */

  PetscErrorCode ierr;
  realtype out;

  switch (op) {
   case 1: ierr = MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_SUM, comm);
           break;

   case 2: ierr = MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MAX, comm);
           break;

   case 3: ierr = MPI_Allreduce(&d, &out, 1, PVEC_REAL_MPI_TYPE, MPI_MIN, comm);
           break;

   default: break;
  }
  CHKERRABORT(comm, ierr);

  return(out);
}

