/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the implementation file for a MPI+RAJA implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/raja/Vector.hpp>
#include <sundials/sundials_mpi.h>
#include <RAJA/RAJA.hpp>


#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

// RAJA defines
#define CUDA_BLOCK_SIZE 256
#define RAJA_NODE_TYPE RAJA::cuda_exec< CUDA_BLOCK_SIZE >
#define RAJA_REDUCE_TYPE RAJA::cuda_reduce< CUDA_BLOCK_SIZE >
#define RAJA_LAMBDA [=] __device__

extern "C" {

using namespace sunrajavec;

// Type defines
typedef sunrajavec::Vector<realtype, sunindextype> vector_type;

// Static constants
static constexpr sunindextype zeroIdx = 0;

/*
 * ----------------------------------------------------------------
 * private accessor/helper functions
 * ----------------------------------------------------------------
 */

static inline sunindextype getLocalLength(N_Vector v)
{
  vector_type* vp = static_cast<vector_type*>(v->content);
  return vp->size();
}

static inline SUNMPI_Comm getMPIComm(N_Vector v)
{
  vector_type* vp = static_cast<vector_type*>(v->content);
  return vp->comm();
}

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Raja(N_Vector v)
{
  return SUNDIALS_NVEC_RAJA;
}

N_Vector N_VNewEmpty_Raja()
{
  N_Vector v;
  N_Vector_Ops ops;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = N_VGetVectorID_Raja;
  ops->nvclone           = N_VClone_Raja;
  ops->nvcloneempty      = N_VCloneEmpty_Raja;
  ops->nvdestroy         = N_VDestroy_Raja;
  ops->nvspace           = N_VSpace_Raja;
  ops->nvgetarraypointer = NULL; //N_VGetArrayPointer_Raja;
  ops->nvsetarraypointer = NULL; //N_VSetArrayPointer_Raja;

  /* standard vector operations */
  ops->nvlinearsum    = N_VLinearSum_Raja;
  ops->nvconst        = N_VConst_Raja;
  ops->nvprod         = N_VProd_Raja;
  ops->nvdiv          = N_VDiv_Raja;
  ops->nvscale        = N_VScale_Raja;
  ops->nvabs          = N_VAbs_Raja;
  ops->nvinv          = N_VInv_Raja;
  ops->nvaddconst     = N_VAddConst_Raja;
  ops->nvdotprod      = N_VDotProd_Raja;
  ops->nvmaxnorm      = N_VMaxNorm_Raja;
  ops->nvwrmsnormmask = N_VWrmsNormMask_Raja;
  ops->nvwrmsnorm     = N_VWrmsNorm_Raja;
  ops->nvmin          = N_VMin_Raja;
  ops->nvwl2norm      = N_VWL2Norm_Raja;
  ops->nvl1norm       = N_VL1Norm_Raja;
  ops->nvcompare      = N_VCompare_Raja;
  ops->nvinvtest      = N_VInvTest_Raja;
  ops->nvconstrmask   = N_VConstrMask_Raja;
  ops->nvminquotient  = N_VMinQuotient_Raja;

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

  /* Attach ops and set content to NULL */
  v->content = NULL;
  v->ops     = ops;

  return(v);
}


#if SUNDIALS_MPI_ENABLED
N_Vector N_VNew_Raja(MPI_Comm comm,
                     sunindextype local_length,
                     sunindextype global_length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  v->content = new vector_type(comm, local_length, global_length);

  return(v);
}
#else
N_Vector N_VNew_Raja(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  v->content = new vector_type(length);

  return(v);
}
#endif


N_Vector N_VMake_Raja(N_VectorContent_Raja c)
{
  N_Vector v;
  vector_type* x = static_cast<vector_type*>(c);
  sunindextype length = x->size();

  v = NULL;
  v = N_VNewEmpty_Raja();
  if (v == NULL) return(NULL);

  v->content = c;

  return(v);
}


/* -----------------------------------------------------------------
 * Function to return the global length of the vector.
 */
sunindextype N_VGetLength_Raja(N_Vector v)
{
  vector_type* xd = static_cast<vector_type*>(v->content);
  return xd->sizeGlobal();
}

#if SUNDIALS_MPI_ENABLED
/* -----------------------------------------------------------------
 * Function to return the local length of the vector.
 */
sunindextype N_VGetLocalLength_Raja(N_Vector v)
{
  vector_type* xd = static_cast<vector_type*>(v->content);
  return xd->size();
}

/* -----------------------------------------------------------------
 * Function to return the MPI communicator for the vector.
 */
MPI_Comm N_VGetMPIComm_Raja(N_Vector v)
{
  vector_type* xd = static_cast<vector_type*>(v->content);
  return (xd->comm());
}
#endif

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

realtype *N_VGetHostArrayPointer_Raja(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->host());
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Raja(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->device());
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Raja(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  xv->copyToDev();
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Raja(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  xv->copyFromDev();
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

void N_VPrint_Raja(N_Vector X)
{
  N_VPrintFile_Raja(X, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void N_VPrintFile_Raja(N_Vector X, FILE *outfile)
{
  const realtype *xd = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);
  sunindextype i;

  for (i = 0; i < N; ++i) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", xd[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", xd[i]);
#else
    fprintf(outfile, "%11.8g\n", xd[i]);
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

N_Vector N_VCloneEmpty_Raja(N_Vector w)
{
  N_Vector v;
  N_Vector_Ops ops;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

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
  v->content = NULL;
  v->ops  = ops;

  return(v);
}

N_Vector N_VClone_Raja(N_Vector w)
{
  N_Vector v;
  vector_type* wdat = static_cast<vector_type*>(w->content);
  vector_type* vdat = new vector_type(*wdat);
  v = NULL;
  v = N_VCloneEmpty_Raja(w);
  if (v == NULL) return(NULL);

  v->content = vdat;

  return(v);
}


void N_VDestroy_Raja(N_Vector v)
{
  vector_type* x = static_cast<vector_type*>(v->content);
  if (x != NULL) {
    delete x;
    v->content = NULL;
  }

  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VSpace_Raja(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  SUNMPI_Comm comm = getMPIComm(X);
  int npes;

  SUNMPI_Comm_size(comm, &npes);

  *lrw = N_VGetLength_Raja(X);
  *liw = 2*npes;
}

void N_VConst_Raja(realtype c, N_Vector Z)
{
  const sunindextype N = getLocalLength(Z);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N), RAJA_LAMBDA(sunindextype i) {
     zdata[i] = c;
  });
}

void N_VLinearSum_Raja(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const realtype *ydata = N_VGetDeviceArrayPointer_Raja(Y);
  const sunindextype N = getLocalLength(X);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = a*xdata[i] + b*ydata[i];
    }
  );
}

void N_VProd_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const realtype *ydata = N_VGetDeviceArrayPointer_Raja(Y);
  const sunindextype N = getLocalLength(X);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] * ydata[i];
    }
  );
}

void N_VDiv_Raja(N_Vector X, N_Vector Y, N_Vector Z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const realtype *ydata = N_VGetDeviceArrayPointer_Raja(Y);
  const sunindextype N = getLocalLength(X);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] / ydata[i];
    }
  );
}

void N_VScale_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = c * xdata[i];
    }
  );
}

void N_VAbs_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = abs(xdata[i]);
    }
  );
}

void N_VInv_Raja(N_Vector X, N_Vector Z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall<RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = ONE / xdata[i];
    }
  );
}

void N_VAddConst_Raja(N_Vector X, realtype b, N_Vector Z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = xdata[i] + b;
    }
  );
}

realtype N_VDotProd_Raja(N_Vector X, N_Vector Y)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const realtype *ydata = N_VGetDeviceArrayPointer_Raja(Y);
  const sunindextype N = getLocalLength(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += xdata[i] * ydata[i] ;
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(X);
  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return gsum;
}

realtype N_VMaxNorm_Raja(N_Vector X)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);

  RAJA::ReduceMax< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result.max(abs(xdata[i]));
    }
  );

  /* Reduce across MPI processes */
  realtype maximum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(X);
  return SUNMPI_Allreduce_scalar(maximum, 2, comm);
}

realtype N_VWrmsNorm_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const realtype *wdata = N_VGetDeviceArrayPointer_Raja(W);
  const sunindextype N = getLocalLength(X);
  const sunindextype Nglobal = N_VGetLength_Raja(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(X);
  return std::sqrt(SUNMPI_Allreduce_scalar(sum, 1, comm)/Nglobal);
}

realtype N_VWrmsNormMask_Raja(N_Vector X, N_Vector W, N_Vector ID)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const realtype *wdata = N_VGetDeviceArrayPointer_Raja(W);
  const realtype *iddata = N_VGetDeviceArrayPointer_Raja(ID);
  const sunindextype N = getLocalLength(X);
  const sunindextype Nglobal = N_VGetLength_Raja(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (iddata[i] > ZERO)
        gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(X);
  return std::sqrt(SUNMPI_Allreduce_scalar(sum, 1, comm)/Nglobal);
}

realtype N_VMin_Raja(N_Vector X)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);

  RAJA::ReduceMin< RAJA_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result.min(xdata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype minumum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(X);
  return SUNMPI_Allreduce_scalar(minumum, 3, comm);
}

realtype N_VWL2Norm_Raja(N_Vector X, N_Vector W)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const realtype *wdata = N_VGetDeviceArrayPointer_Raja(W);
  const sunindextype N = getLocalLength(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (xdata[i] * wdata[i] * xdata[i] * wdata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(X);
  return std::sqrt(SUNMPI_Allreduce_scalar(sum, 1, comm));
}

realtype N_VL1Norm_Raja(N_Vector X)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(0.0);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      gpu_result += (abs(xdata[i]));
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(X);
  return SUNMPI_Allreduce_scalar(sum, 1, comm);
}

void N_VCompare_Raja(realtype c, N_Vector X, N_Vector Z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(X);
  const sunindextype N = getLocalLength(X);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(Z);

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      zdata[i] = abs(xdata[i]) >= c ? ONE : ZERO;
    }
  );
}

booleantype N_VInvTest_Raja(N_Vector x, N_Vector z)
{
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(x);
  const sunindextype N = getLocalLength(x);
  realtype *zdata = N_VGetDeviceArrayPointer_Raja(z);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(ZERO);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (xdata[i] == ZERO) {
        gpu_result += ONE;
      } else {
        zdata[i] = ONE/xdata[i];
      }
    }
  );

  /* Reduce across MPI processes */
  realtype minimum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(x);
  realtype global_minimum = SUNMPI_Allreduce_scalar(minimum, 3, comm);

  return (global_minimum < HALF);
}

booleantype N_VConstrMask_Raja(N_Vector c, N_Vector x, N_Vector m)
{
  const realtype *cdata = N_VGetDeviceArrayPointer_Raja(c);
  const realtype *xdata = N_VGetDeviceArrayPointer_Raja(x);
  const sunindextype N = getLocalLength(x);
  realtype *mdata = N_VGetDeviceArrayPointer_Raja(m);

  RAJA::ReduceSum< RAJA_REDUCE_TYPE, realtype> gpu_result(ZERO);
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      bool test = (abs(cdata[i]) > ONEPT5 && cdata[i]*xdata[i] <= ZERO) ||
                  (abs(cdata[i]) > HALF   && cdata[i]*xdata[i] <  ZERO);
      mdata[i] = test ? ONE : ZERO;
      gpu_result += mdata[i];
    }
  );

  /* Reduce across MPI processes */
  realtype sum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(x);
  realtype global_sum = SUNMPI_Allreduce_scalar(sum, 1, comm);

  return (global_sum < HALF);
}

realtype N_VMinQuotient_Raja(N_Vector num, N_Vector denom)
{
  const realtype *ndata = N_VGetDeviceArrayPointer_Raja(num);
  const realtype *ddata = N_VGetDeviceArrayPointer_Raja(denom);
  const sunindextype N = getLocalLength(num);

  RAJA::ReduceMin< RAJA_REDUCE_TYPE, realtype> gpu_result(std::numeric_limits<realtype>::max());
  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      if (ddata[i] != ZERO)
        gpu_result.min(ndata[i]/ddata[i]);
    }
  );

  /* Reduce across MPI processes */
  realtype minimum = static_cast<realtype>(gpu_result);
  SUNMPI_Comm comm = getMPIComm(num);
  return SUNMPI_Allreduce_scalar(minimum, 3, comm);
}


/*
 * -----------------------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearCombination_Raja(int nvec, realtype* c, N_Vector* X, N_Vector z)
{
  cudaError_t  err;

  sunindextype N = getLocalLength(z);
  realtype* d_zd = N_VGetDeviceArrayPointer_Raja(z);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = N_VGetDeviceArrayPointer_Raja(X[j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      d_zd[i] = d_c[0] * d_Xd[0][i];
      for (int j=1; j<nvec; j++)
        d_zd[i] += d_c[j] * d_Xd[j][i];
    }
  );

  // Free host array
  delete[] h_Xd;

  // Free device arrays
  err = cudaFree(d_c);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VScaleAddMulti_Raja(int nvec, realtype* c, N_Vector x, N_Vector* Y, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getLocalLength(x);
  realtype* d_xd = N_VGetDeviceArrayPointer_Raja(x);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Yd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Yd[j] = N_VGetDeviceArrayPointer_Raja(Y[j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = N_VGetDeviceArrayPointer_Raja(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = d_c[j] * d_xd[i] + d_Yd[j][i];
    }
  );

  // Free host array
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_c);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Raja(int nvec,
                                 realtype a, N_Vector* X,
                                 realtype b, N_Vector* Y,
                                 N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getLocalLength(Z[0]);

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = N_VGetDeviceArrayPointer_Raja(X[j]);

  realtype** h_Yd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Yd[j] = N_VGetDeviceArrayPointer_Raja(Y[j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = N_VGetDeviceArrayPointer_Raja(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = a * d_Xd[j][i] + b * d_Yd[j][i];
    }
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VScaleVectorArray_Raja(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getLocalLength(Z[0]);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nvec*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nvec*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = N_VGetDeviceArrayPointer_Raja(X[j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = N_VGetDeviceArrayPointer_Raja(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = d_c[j] * d_Xd[j][i];
    }
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VConstVectorArray_Raja(int nvec, realtype c, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getLocalLength(Z[0]);

  // Create array of device pointers on host
  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = N_VGetDeviceArrayPointer_Raja(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        d_Zd[j][i] = c;
    }
  );

  // Free host array
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VScaleAddMultiVectorArray_Raja(int nvec, int nsum, realtype* c,
                                     N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  cudaError_t err;

  sunindextype N = getLocalLength(X[0]);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Xd[j] = N_VGetDeviceArrayPointer_Raja(X[j]);

  realtype** h_Yd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Yd[j*nsum+k] = N_VGetDeviceArrayPointer_Raja(Y[k][j]);

  realtype** h_Zd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Zd[j*nsum+k] = N_VGetDeviceArrayPointer_Raja(Z[k][j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Yd;
  err = cudaMalloc((void**) &d_Yd, nsum*nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Yd, h_Yd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nsum*nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++)
        for (int k=0; k<nsum; k++)
          d_Zd[j*nsum+k][i] = d_c[k] * d_Xd[j][i] + d_Yd[j*nsum+k][i];
    }
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Yd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Yd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


int N_VLinearCombinationVectorArray_Raja(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  cudaError_t err;

  sunindextype N = getLocalLength(Z[0]);

  // Copy c array to device
  realtype* d_c;
  err = cudaMalloc((void**) &d_c, nsum*sizeof(realtype));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_c, c, nsum*sizeof(realtype), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  // Create array of device pointers on host
  realtype** h_Xd = new realtype*[nsum*nvec];
  for (int j=0; j<nvec; j++)
    for (int k=0; k<nsum; k++)
      h_Xd[j*nsum+k] = N_VGetDeviceArrayPointer_Raja(X[k][j]);

  realtype** h_Zd = new realtype*[nvec];
  for (int j=0; j<nvec; j++)
    h_Zd[j] = N_VGetDeviceArrayPointer_Raja(Z[j]);

  // Copy array of device pointers to device from host
  realtype** d_Xd;
  err = cudaMalloc((void**) &d_Xd, nsum*nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Xd, h_Xd, nsum*nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  realtype** d_Zd;
  err = cudaMalloc((void**) &d_Zd, nvec*sizeof(realtype*));
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaMemcpy(d_Zd, h_Zd, nvec*sizeof(realtype*), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) return cudaGetLastError();

  RAJA::forall< RAJA_NODE_TYPE >(RAJA::RangeSegment(zeroIdx, N),
    RAJA_LAMBDA(sunindextype i) {
      for (int j=0; j<nvec; j++) {
        d_Zd[j][i] = d_c[0] * d_Xd[j*nsum][i];
        for (int k=1; k<nsum; k++) {
          d_Zd[j][i] += d_c[k] * d_Xd[j*nsum+k][i];
        }
      }
    }
  );

  // Free host array
  delete[] h_Xd;
  delete[] h_Zd;

  // Free device arrays
  err = cudaFree(d_Xd);
  if (err != cudaSuccess) return cudaGetLastError();
  err = cudaFree(d_Zd);
  if (err != cudaSuccess) return cudaGetLastError();

  return(0);
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Raja;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Raja;
    v->ops->nvdotprodmulti      = NULL;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Raja;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Raja;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Raja;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Raja;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Raja;
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


int N_VEnableLinearCombination_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Raja;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Raja;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Raja;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Raja;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Raja;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Raja;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Raja(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Raja;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}

} // extern "C"
