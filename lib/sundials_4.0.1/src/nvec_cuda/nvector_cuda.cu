/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles, and Cody J. Balos @ LLNL
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
 * This is the implementation file for a MPI+CUDA implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#include <nvector/cuda/Vector.hpp>
#include <nvector/cuda/VectorKernels.cuh>
#include <nvector/cuda/VectorArrayKernels.cuh>
#include <sundials/sundials_mpi.h>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)

extern "C" {

using namespace suncudavec;

/*
 * Type definitions
 */

typedef suncudavec::Vector<realtype, sunindextype> vector_type;

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Cuda(N_Vector v)
{
  return SUNDIALS_NVEC_CUDA;
}

N_Vector N_VNewEmpty_Cuda()
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

  ops->nvgetvectorid     = N_VGetVectorID_Cuda;
  ops->nvclone           = N_VClone_Cuda;
  ops->nvcloneempty      = N_VCloneEmpty_Cuda;
  ops->nvdestroy         = N_VDestroy_Cuda;
  ops->nvspace           = N_VSpace_Cuda;
  ops->nvgetarraypointer = NULL;
  ops->nvsetarraypointer = NULL;

  /* standard vector operations */
  ops->nvlinearsum    = N_VLinearSum_Cuda;
  ops->nvconst        = N_VConst_Cuda;
  ops->nvprod         = N_VProd_Cuda;
  ops->nvdiv          = N_VDiv_Cuda;
  ops->nvscale        = N_VScale_Cuda;
  ops->nvabs          = N_VAbs_Cuda;
  ops->nvinv          = N_VInv_Cuda;
  ops->nvaddconst     = N_VAddConst_Cuda;
  ops->nvdotprod      = N_VDotProd_Cuda;
  ops->nvmaxnorm      = N_VMaxNorm_Cuda;
  ops->nvwrmsnormmask = N_VWrmsNormMask_Cuda;
  ops->nvwrmsnorm     = N_VWrmsNorm_Cuda;
  ops->nvmin          = N_VMin_Cuda;
  ops->nvwl2norm      = N_VWL2Norm_Cuda;
  ops->nvl1norm       = N_VL1Norm_Cuda;
  ops->nvcompare      = N_VCompare_Cuda;
  ops->nvinvtest      = N_VInvTest_Cuda;
  ops->nvconstrmask   = N_VConstrMask_Cuda;
  ops->nvminquotient  = N_VMinQuotient_Cuda;

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
N_Vector N_VNew_Cuda(MPI_Comm comm,
                     sunindextype local_length,
                     sunindextype global_length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  v->content = new vector_type(comm, local_length, global_length);

  return(v);
}

N_Vector N_VNewManaged_Cuda(MPI_Comm comm,
                            sunindextype local_length,
                            sunindextype global_length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* create suncudavec::Vector with managed memory */
  v->content = new vector_type(comm, local_length, global_length, true);

  return(v);
}

N_Vector N_VMake_Cuda(MPI_Comm comm,
                      sunindextype local_length, sunindextype global_length,
                      realtype *h_vdata, realtype *d_vdata)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* create suncudavec::Vector using the user-provided data arrays */
  v->content = new vector_type(comm, local_length, global_length, false, false, h_vdata, d_vdata);

  return(v);
}

N_Vector N_VMakeManaged_Cuda(MPI_Comm comm,
                             sunindextype local_length, sunindextype global_length,
                             realtype *vdata)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);
  
  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* create suncudavec::Vector with managed memory using the user-provided data arrays */
  v->content = new vector_type(comm, local_length, global_length, true, false, vdata, vdata);

  return(v);
}
#else
N_Vector N_VNew_Cuda(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  v->content = new vector_type(length);

  return(v);
}

N_Vector N_VNewManaged_Cuda(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* create suncudavec::Vector with managed memory */
  v->content = new vector_type(length, true);

  return(v);
}

N_Vector N_VMake_Cuda(sunindextype length, realtype *h_vdata, realtype *d_vdata)
{
  N_Vector v;

  if (h_vdata == NULL || d_vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* create suncudavec::Vector using the user-provided data arrays */
  v->content = new vector_type(length, false, false, h_vdata, d_vdata);

  return(v);
}

N_Vector N_VMakeManaged_Cuda(sunindextype length, realtype *vdata)
{
  N_Vector v;

  if (vdata == NULL) return(NULL);

  v = NULL;
  v = N_VNewEmpty_Cuda();
  if (v == NULL) return(NULL);

  /* create suncudavec::Vector with managed memory using the user-provided data arrays */
  v->content = new vector_type(length, true, false, vdata, vdata);

  return(v);
}
#endif

/* -----------------------------------------------------------------
 * Function to return the global length of the vector.
 */
sunindextype N_VGetLength_Cuda(N_Vector v)
{
  vector_type* xd = static_cast<vector_type*>(v->content);
  return (xd->sizeGlobal());
}


#if SUNDIALS_MPI_ENABLED
/* -----------------------------------------------------------------
 * Function to return the local length of the vector.
 */
sunindextype N_VGetLocalLength_Cuda(N_Vector v)
{
  vector_type* xd = static_cast<vector_type*>(v->content);
  return (xd->size());
}

/* -----------------------------------------------------------------
 * Function to return the MPI communicator for the vector.
 */
MPI_Comm N_VGetMPIComm_Cuda(N_Vector v)
{
  vector_type* xd = static_cast<vector_type*>(v->content);
  return (xd->comm());
}
#endif

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

realtype *N_VGetHostArrayPointer_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->host());
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->device());
}

/* ----------------------------------------------------------------------------
 * Return a flag indicating if the memory for the vector data is managed
 */
booleantype N_VIsManagedMemory_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  return (xv->isManaged());
}

/*
 * ----------------------------------------------------------------------------
 * Sets the cudaStream_t to use for execution of the CUDA kernels.
 */
void N_VSetCudaStream_Cuda(N_Vector x, cudaStream_t *stream)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  xv->partStream().setStream(*stream);
  xv->partReduce().setStream(*stream);
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  xv->copyToDev();
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Cuda(N_Vector x)
{
  vector_type* xv = static_cast<vector_type*>(x->content);
  xv->copyFromDev();
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to stdout
 */

void N_VPrint_Cuda(N_Vector x)
{
  N_VPrintFile_Cuda(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a CUDA-based vector to outfile
 */

void N_VPrintFile_Cuda(N_Vector x, FILE *outfile)
{
  sunindextype i;
  vector_type* xd = static_cast<vector_type*>(x->content);

  for (i = 0; i < xd->size(); i++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Lg\n", xd->host()[i]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16g\n", xd->host()[i]);
#else
    fprintf(outfile, "%11.8g\n", xd->host()[i]);
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

N_Vector N_VCloneEmpty_Cuda(N_Vector w)
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

N_Vector N_VClone_Cuda(N_Vector w)
{
  N_Vector v;
  vector_type* wdat = static_cast<vector_type*>(w->content);
  vector_type* vdat = new vector_type(*wdat);
  v = NULL;
  v = N_VCloneEmpty_Cuda(w);
  if (v == NULL) return(NULL);

  v->content = vdat;

  return(v);
}


void N_VDestroy_Cuda(N_Vector v)
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

void N_VSpace_Cuda(N_Vector X, sunindextype *lrw, sunindextype *liw)
{
  int npes;
  vector_type* x = static_cast<vector_type*>(X->content);

  SUNMPI_Comm_size(x->comm(), &npes);

  *lrw = x->sizeGlobal();
  *liw = 2*npes;
}

void N_VConst_Cuda(realtype a, N_Vector X)
{
  vector_type *xvec = static_cast<vector_type*>(X->content);
  setConst(a, *xvec);
}

void N_VLinearSum_Cuda(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  const vector_type *yvec = static_cast<vector_type*>(Y->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  linearSum(a, *xvec, b, *yvec, *zvec);
}

void N_VProd_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  const vector_type *yvec = static_cast<vector_type*>(Y->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  prod(*xvec, *yvec, *zvec);
}

void N_VDiv_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  const vector_type *yvec = static_cast<vector_type*>(Y->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  div(*xvec, *yvec, *zvec);
}

void N_VScale_Cuda(realtype a, N_Vector X, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  scale(a, *xvec, *zvec);
}

void N_VAbs_Cuda(N_Vector X, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  absVal(*xvec, *zvec);
}

void N_VInv_Cuda(N_Vector X, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  inv(*xvec, *zvec);
}

void N_VAddConst_Cuda(N_Vector X, realtype b, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  addConst(b, *xvec, *zvec);
}

realtype N_VDotProd_Cuda(N_Vector X, N_Vector Y)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  const vector_type *yvec = static_cast<vector_type*>(Y->content);
  SUNMPI_Comm comm = xvec->comm();

  realtype sum = dotProd(*xvec, *yvec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return gsum;
}

realtype N_VMaxNorm_Cuda(N_Vector X)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  SUNMPI_Comm comm = xvec->comm();

  realtype locmax = maxNorm(*xvec);

  realtype globmax = SUNMPI_Allreduce_scalar(locmax, 2, comm);
  return globmax;
}

realtype N_VWrmsNorm_Cuda(N_Vector X, N_Vector W)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  const vector_type *wvec = static_cast<vector_type*>(W->content);
  const sunindextype Nglob = xvec->sizeGlobal();
  SUNMPI_Comm comm = xvec->comm();

  realtype sum = wL2NormSquare(*xvec, *wvec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return std::sqrt(gsum/Nglob);
}

realtype N_VWrmsNormMask_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  const vector_type *wvec = static_cast<vector_type*>(W->content);
  const vector_type *ivec = static_cast<vector_type*>(Id->content);
  const sunindextype Nglob = xvec->sizeGlobal();
  SUNMPI_Comm comm = xvec->comm();

  realtype sum = wL2NormSquareMask(*xvec, *wvec, *ivec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return std::sqrt(gsum/Nglob);
}

realtype N_VMin_Cuda(N_Vector X)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  SUNMPI_Comm comm = xvec->comm();

  realtype locmin = findMin(*xvec);

  realtype globmin = SUNMPI_Allreduce_scalar(locmin, 3, comm);
  return globmin;
}

realtype N_VWL2Norm_Cuda(N_Vector X, N_Vector W)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  const vector_type *wvec = static_cast<vector_type*>(W->content);
  SUNMPI_Comm comm = xvec->comm();

  realtype sum = wL2NormSquare(*xvec, *wvec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return std::sqrt(gsum);
}

realtype N_VL1Norm_Cuda(N_Vector X)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  SUNMPI_Comm comm = xvec->comm();

  realtype sum = L1Norm(*xvec);

  realtype gsum = SUNMPI_Allreduce_scalar(sum, 1, comm);
  return gsum;
}

void N_VCompare_Cuda(realtype c, N_Vector X, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  compare(c, *xvec, *zvec);
}

booleantype N_VInvTest_Cuda(N_Vector X, N_Vector Z)
{
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  vector_type *zvec = static_cast<vector_type*>(Z->content);
  SUNMPI_Comm comm = xvec->comm();
  
  realtype locmin = invTest(*xvec, *zvec);

  realtype globmin = SUNMPI_Allreduce_scalar(locmin, 3, comm);
  return (globmin < HALF);
}

/*
 * Creates mask for variables violating constraints
 */
booleantype N_VConstrMask_Cuda(N_Vector C, N_Vector X, N_Vector M)
{
  const vector_type *cvec = static_cast<vector_type*>(C->content);
  const vector_type *xvec = static_cast<vector_type*>(X->content);
  vector_type *mvec = static_cast<vector_type*>(M->content);
  SUNMPI_Comm comm = xvec->comm();

  realtype locsum = constrMask(*cvec, *xvec, *mvec);

  realtype globsum = SUNMPI_Allreduce_scalar(locsum, 1, comm);
  return (globsum < HALF);
}

realtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom)
{
  const vector_type *numvec = static_cast<vector_type*>(num->content);
  const vector_type *denvec = static_cast<vector_type*>(denom->content);
  SUNMPI_Comm comm = numvec->comm();

  realtype locmin = minQuotient(*numvec, *denvec);

  realtype globmin = SUNMPI_Allreduce_scalar(locmin, 3, comm);
  return globmin;
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

int N_VLinearCombination_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type*  Zv;

  Zv = static_cast<vector_type*>(Z->content);

  Xv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Xv[i] = static_cast<vector_type*>(X[i]->content);

  err = linearCombination(nvec, c, Xv, Zv);

  delete[] Xv;

  return err == cudaSuccess ? 0 : -1;
}

int N_VScaleAddMulti_Cuda(int nvec, realtype* c, N_Vector X, N_Vector* Y,
                           N_Vector* Z)
{
  cudaError_t err;
  vector_type*  Xv;
  vector_type** Yv;
  vector_type** Zv;

  Xv = static_cast<vector_type*>(X->content);

  Yv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Yv[i] = static_cast<vector_type*>(Y[i]->content);

  Zv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Zv[i] = static_cast<vector_type*>(Z[i]->content);

  err = scaleAddMulti(nvec, c, Xv, Yv, Zv);

  delete[] Yv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VDotProdMulti_Cuda(int nvec, N_Vector x, N_Vector* Y, realtype* dotprods)
{
  cudaError_t err;
  SUNMPI_Comm comm;
  vector_type*  Xv;
  vector_type** Yv;

  Xv = static_cast<vector_type*>(x->content);
  comm = Xv->comm();

  Yv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Yv[i] = static_cast<vector_type*>(Y[i]->content);

  err = dotProdMulti(nvec, Xv, Yv, dotprods);

  delete[] Yv;

  SUNMPI_Allreduce(dotprods, nvec, 1, comm);

  return err == cudaSuccess ? 0 : -1;
}



/*
 * -----------------------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------------------
 */

int N_VLinearSumVectorArray_Cuda(int nvec, realtype a, N_Vector* X, realtype b,
                                 N_Vector* Y, N_Vector* Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type** Yv;
  vector_type** Zv;

  Xv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Xv[i] = static_cast<vector_type*>(X[i]->content);

  Yv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Yv[i] = static_cast<vector_type*>(Y[i]->content);

  Zv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Zv[i] = static_cast<vector_type*>(Z[i]->content);

  err = linearSumVectorArray(nvec, a, Xv, b, Yv, Zv);

  delete[] Xv;
  delete[] Yv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VScaleVectorArray_Cuda(int nvec, realtype* c, N_Vector* X, N_Vector* Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type** Zv;

  Xv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Xv[i] = static_cast<vector_type*>(X[i]->content);

  Zv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Zv[i] = static_cast<vector_type*>(Z[i]->content);

  err = scaleVectorArray(nvec, c, Xv, Zv);

  delete[] Xv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VConstVectorArray_Cuda(int nvec, realtype c, N_Vector* Z)
{
  cudaError_t err;
  vector_type** Zv;

  Zv = new vector_type*[nvec];
  for (int i=0; i<nvec; i++)
    Zv[i] = static_cast<vector_type*>(Z[i]->content);

  err = constVectorArray(nvec, c, Zv);

  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VWrmsNormVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                realtype* norms)
{
  cudaError_t err;
  const vector_type* xvec = static_cast<vector_type*>(X[0]->content);
  vector_type** Xv;
  vector_type** Wv;
  
  SUNMPI_Comm comm = xvec->comm();
  sunindextype N = xvec->sizeGlobal();

  Xv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Xv[k] = static_cast<vector_type*>(X[k]->content);

  Wv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Wv[k] = static_cast<vector_type*>(W[k]->content);

  err = wL2NormSquareVectorArray(nvec, Xv, Wv, norms);

  delete[] Xv;
  delete[] Wv;

  SUNMPI_Allreduce(norms, nvec, 1, comm);

  for (int k=0; k<nvec; ++k) {
    norms[k] = std::sqrt(norms[k]/N);
  }

  return err == cudaSuccess ? 0 : -1;
}


int N_VWrmsNormMaskVectorArray_Cuda(int nvec, N_Vector* X, N_Vector* W,
                                    N_Vector id, realtype* norms)
{
  cudaError_t err;
  const vector_type* xvec = static_cast<vector_type*>(X[0]->content);
  vector_type** Xv;
  vector_type** Wv;
  vector_type*  IDv;
  
  SUNMPI_Comm comm = xvec->comm();
  sunindextype N = xvec->sizeGlobal();

  Xv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Xv[k] = static_cast<vector_type*>(X[k]->content);

  Wv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Wv[k] = static_cast<vector_type*>(W[k]->content);

  IDv = static_cast<vector_type*>(id->content);

  err = wL2NormSquareMaskVectorArray(nvec, Xv, Wv, IDv, norms);

  delete[] Xv;
  delete[] Wv;

  SUNMPI_Allreduce(norms, nvec, 1, comm);

  for (int k=0; k<nvec; ++k) {
    norms[k] = std::sqrt(norms[k]/N);
  }

  return err == cudaSuccess ? 0 : -1;
}


int N_VScaleAddMultiVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                      N_Vector* X, N_Vector** Y, N_Vector** Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type** Yv;
  vector_type** Zv;

  Xv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Xv[k] = static_cast<vector_type*>(X[k]->content);

  Yv = new vector_type*[nsum*nvec];
  for (int k=0; k<nvec; k++)
    for (int j=0; j<nsum; j++)
      Yv[k*nsum+j] = static_cast<vector_type*>(Y[j][k]->content);

  Zv = new vector_type*[nsum*nvec];
  for (int k=0; k<nvec; k++)
    for (int j=0; j<nsum; j++)
      Zv[k*nsum+j] = static_cast<vector_type*>(Z[j][k]->content);

  err = scaleAddMultiVectorArray(nvec, nsum, c, Xv, Yv, Zv);

  delete[] Xv;
  delete[] Yv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


int N_VLinearCombinationVectorArray_Cuda(int nvec, int nsum, realtype* c,
                                         N_Vector** X, N_Vector* Z)
{
  cudaError_t err;
  vector_type** Xv;
  vector_type** Zv;

  Xv = new vector_type*[nsum*nvec];
  for (int k=0; k<nvec; k++)
    for (int j=0; j<nsum; j++)
      Xv[k*nsum+j] = static_cast<vector_type*>(X[j][k]->content);

  Zv = new vector_type*[nvec];
  for (int k=0; k<nvec; k++)
    Zv[k] = static_cast<vector_type*>(Z[k]->content);

  err = linearCombinationVectorArray(nvec, nsum, c, Xv, Zv);

  delete[] Xv;
  delete[] Zv;

  return err == cudaSuccess ? 0 : -1;
}


/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

int N_VEnableFusedOps_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  if (tf) {
    /* enable all fused vector operations */
    v->ops->nvlinearcombination = N_VLinearCombination_Cuda;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Cuda;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Cuda;
    /* enable all vector array operations */
    v->ops->nvlinearsumvectorarray         = N_VLinearSumVectorArray_Cuda;
    v->ops->nvscalevectorarray             = N_VScaleVectorArray_Cuda;
    v->ops->nvconstvectorarray             = N_VConstVectorArray_Cuda;
    v->ops->nvwrmsnormvectorarray          = N_VWrmsNormVectorArray_Cuda;
    v->ops->nvwrmsnormmaskvectorarray      = N_VWrmsNormMaskVectorArray_Cuda;
    v->ops->nvscaleaddmultivectorarray     = N_VScaleAddMultiVectorArray_Cuda;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Cuda;
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


int N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombination = N_VLinearCombination_Cuda;
  else
    v->ops->nvlinearcombination = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmulti = N_VScaleAddMulti_Cuda;
  else
    v->ops->nvscaleaddmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvdotprodmulti = N_VDotProdMulti_Cuda;
  else
    v->ops->nvdotprodmulti = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearsumvectorarray = N_VLinearSumVectorArray_Cuda;
  else
    v->ops->nvlinearsumvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscalevectorarray = N_VScaleVectorArray_Cuda;
  else
    v->ops->nvscalevectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvconstvectorarray = N_VConstVectorArray_Cuda;
  else
    v->ops->nvconstvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormvectorarray = N_VWrmsNormVectorArray_Cuda;
  else
    v->ops->nvwrmsnormvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvwrmsnormmaskvectorarray = N_VWrmsNormMaskVectorArray_Cuda;
  else
    v->ops->nvwrmsnormmaskvectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Cuda;
  else
    v->ops->nvscaleaddmultivectorarray = NULL;

  /* return success */
  return(0);
}

int N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf)
{
  /* check that vector is non-NULL */
  if (v == NULL) return(-1);

  /* check that ops structure is non-NULL */
  if (v->ops == NULL) return(-1);

  /* enable/disable operation */
  if (tf)
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Cuda;
  else
    v->ops->nvlinearcombinationvectorarray = NULL;

  /* return success */
  return(0);
}

} // extern "C"
