/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
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
 * This is the implementation file for a serial implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/cuda/Vector.hpp>
#include <nvector/cuda/VectorKernels.cuh>

extern "C" {

using namespace suncudavec;

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Cuda(N_Vector v)
{
  return SUNDIALS_NVEC_CUDA;
}

N_Vector N_VNewEmpty_Cuda(sunindextype length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Cuda content;

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
  ops->nvlinearsum       = N_VLinearSum_Cuda;
  ops->nvconst           = N_VConst_Cuda;
  ops->nvprod            = N_VProd_Cuda;
  ops->nvdiv             = N_VDiv_Cuda;
  ops->nvscale           = N_VScale_Cuda;
  ops->nvabs             = N_VAbs_Cuda;
  ops->nvinv             = N_VInv_Cuda;
  ops->nvaddconst        = N_VAddConst_Cuda;
  ops->nvdotprod         = N_VDotProd_Cuda;
  ops->nvmaxnorm         = N_VMaxNorm_Cuda;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Cuda;
  ops->nvwrmsnorm        = N_VWrmsNorm_Cuda;
  ops->nvmin             = N_VMin_Cuda;
  ops->nvwl2norm         = N_VWL2Norm_Cuda;
  ops->nvl1norm          = N_VL1Norm_Cuda;
  ops->nvcompare         = N_VCompare_Cuda;
  ops->nvinvtest         = N_VInvTest_Cuda;
  ops->nvconstrmask      = N_VConstrMask_Cuda;
  ops->nvminquotient     = N_VMinQuotient_Cuda;

  /* Create content */
  content = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}


N_Vector N_VNew_Cuda(sunindextype length)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Cuda(length);
  if (v == NULL)
    return(NULL);

  v->content = new Vector<realtype, sunindextype>(length);

  return(v);
}


N_Vector N_VMake_Cuda(N_VectorContent_Cuda c)
{
  N_Vector v;
  Vector<realtype, sunindextype>* x = static_cast<Vector<realtype, sunindextype>*>(c);
  sunindextype length = x->size();

  v = NULL;
  v = N_VNewEmpty_Cuda(length);
  if (v == NULL) return(NULL);

  v->content = c;

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new CUDA-based vectors.
 */

N_Vector *N_VCloneVectorArray_Cuda(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_Cuda(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Cuda(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new serial vectors with NULL data array.
 */

N_Vector *N_VCloneVectorArrayEmpty_Cuda(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc(count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_Cuda(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Cuda(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* -----------------------------------------------------------------
 * Function to return the length of the vector.
 */
sunindextype N_VGetLength_Cuda(N_Vector v)
{
  Vector<realtype, sunindextype>* xd = static_cast<Vector<realtype, sunindextype>*>(v->content);
  return xd->size();
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Cuda
 */

void N_VDestroyVectorArray_Cuda(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Cuda(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw host data
 */

realtype *N_VGetHostArrayPointer_Cuda(N_Vector x)
{
  Vector<realtype, sunindextype>* xv = static_cast<Vector<realtype, sunindextype>*>(x->content);
  return (xv->host());
}

/* ----------------------------------------------------------------------------
 * Return pointer to the raw device data
 */

realtype *N_VGetDeviceArrayPointer_Cuda(N_Vector x)
{
  Vector<realtype, sunindextype>* xv = static_cast<Vector<realtype, sunindextype>*>(x->content);
  return (xv->device());
}

/* ----------------------------------------------------------------------------
 * Copy vector data to the device
 */

void N_VCopyToDevice_Cuda(N_Vector x)
{
  Vector<realtype, sunindextype>* xv = static_cast<Vector<realtype, sunindextype>*>(x->content);
  xv->copyToDev();
}

/* ----------------------------------------------------------------------------
 * Copy vector data from the device to the host
 */

void N_VCopyFromDevice_Cuda(N_Vector x)
{
  Vector<realtype, sunindextype>* xv = static_cast<Vector<realtype, sunindextype>*>(x->content);
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
  Vector<realtype, sunindextype>* xd = static_cast<Vector<realtype, sunindextype>*>(x->content);

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
  v->content = NULL;
  v->ops  = ops;

  return(v);
}

N_Vector N_VClone_Cuda(N_Vector w)
{
  N_Vector v;
  Vector<realtype, sunindextype>* wdat = static_cast<Vector<realtype, sunindextype>*>(w->content);
  Vector<realtype, sunindextype>* vdat = new Vector<realtype, sunindextype>(*wdat);
  v = NULL;
  v = N_VCloneEmpty_Cuda(w);
  if (v == NULL) return(NULL);

  v->content = vdat;

  return(v);
}


void N_VDestroy_Cuda(N_Vector v)
{
  Vector<realtype, sunindextype>* x = static_cast<Vector<realtype, sunindextype>*>(v->content);
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
  *lrw = (extract<realtype, sunindextype>(X))->size();
  *liw = 1;
}

void N_VConst_Cuda(realtype a, N_Vector X)
{
  setConst(a, *extract<realtype, sunindextype>(X));
}

void N_VLinearSum_Cuda(realtype a, N_Vector X, realtype b, N_Vector Y, N_Vector Z)
{
  linearSum(a, *extract<realtype, sunindextype>(X), b, *extract<realtype, sunindextype>(Y), *extract<realtype, sunindextype>(Z));
}

void N_VProd_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  prod(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Y), *extract<realtype, sunindextype>(Z));
}

void N_VDiv_Cuda(N_Vector X, N_Vector Y, N_Vector Z)
{
  div(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Y), *extract<realtype, sunindextype>(Z));
}

void N_VScale_Cuda(realtype a, N_Vector X, N_Vector Z)
{
  scale(a, *extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Z));
}

void N_VAbs_Cuda(N_Vector X, N_Vector Z)
{
  absVal(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Z));
}

void N_VInv_Cuda(N_Vector X, N_Vector Z)
{
  inv(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Z));
}

void N_VAddConst_Cuda(N_Vector X, realtype b, N_Vector Z)
{
  addConst(b, *extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Z));
}

realtype N_VDotProd_Cuda(N_Vector X, N_Vector Y)
{
  return (dotProd(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Y)));
}

realtype N_VMaxNorm_Cuda(N_Vector X)
{
  return (maxNorm(*extract<realtype, sunindextype>(X)));
}

realtype N_VWrmsNorm_Cuda(N_Vector X, N_Vector W)
{
  return (wrmsNorm(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(W)));
}

realtype N_VWrmsNormMask_Cuda(N_Vector X, N_Vector W, N_Vector Id)
{
  return (wrmsNormMask(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(W), *extract<realtype, sunindextype>(Id)));
}

realtype N_VMin_Cuda(N_Vector X)
{
  return (findMin(*extract<realtype, sunindextype>(X)));
}

realtype N_VWL2Norm_Cuda(N_Vector X, N_Vector W)
{
  return (wL2Norm(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(W)));
}

realtype N_VL1Norm_Cuda(N_Vector X)
{
  return (L1Norm(*extract<realtype, sunindextype>(X)));
}

void N_VCompare_Cuda(realtype c, N_Vector X, N_Vector Z)
{
  compare(c, *extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Z));
}

booleantype N_VInvTest_Cuda(N_Vector X, N_Vector Z)
{
  return (booleantype) (invTest(*extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(Z)));
}

booleantype N_VConstrMask_Cuda(N_Vector C, N_Vector X, N_Vector M)
{
  return (booleantype) (constrMask(*extract<realtype, sunindextype>(C), *extract<realtype, sunindextype>(X), *extract<realtype, sunindextype>(M)));
}

realtype N_VMinQuotient_Cuda(N_Vector num, N_Vector denom)
{
  return (minQuotient(*extract<realtype, sunindextype>(num), *extract<realtype, sunindextype>(denom)));
}

} // extern "C"
