/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 *
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
 * This is the implementation file for a Trilinos implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <nvector/nvector_trilinos.h>
#include <nvector/trilinos/SundialsTpetraVectorInterface.hpp>
#include <nvector/trilinos/SundialsTpetraVectorKernels.hpp>

#define ZERO   RCONST(0.0)
#define HALF   RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)



/*
 * -----------------------------------------------------------------
 * using statements
 * -----------------------------------------------------------------
 */

using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::outArg;
using Teuchos::REDUCE_SUM;
using Teuchos::reduceAll;

/*
 * -----------------------------------------------------------------
 * type definitions
 * -----------------------------------------------------------------
 */

typedef Sundials::TpetraVectorInterface::vector_type vector_type;

/* ----------------------------------------------------------------
 * Returns vector type ID. Used to identify vector implementation
 * from abstract N_Vector interface.
 */
N_Vector_ID N_VGetVectorID_Trilinos(N_Vector v)
{
  return SUNDIALS_NVEC_TRILINOS;
}


/* ----------------------------------------------------------------
 * Function to create a new Trilinos vector with empty data array
 */

N_Vector N_VNewEmpty_Trilinos()
{
  N_Vector v;

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty();
  if (v == NULL) return(NULL);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Trilinos;
  v->ops->nvclone           = N_VClone_Trilinos;
  v->ops->nvcloneempty      = N_VCloneEmpty_Trilinos;
  v->ops->nvdestroy         = N_VDestroy_Trilinos;
  v->ops->nvspace           = N_VSpace_Trilinos;
  v->ops->nvgetcommunicator = N_VGetCommunicator_Trilinos;
  v->ops->nvgetlength       = N_VGetLength_Trilinos;

  /* standard vector operations */
  v->ops->nvlinearsum       = N_VLinearSum_Trilinos;
  v->ops->nvconst           = N_VConst_Trilinos;
  v->ops->nvprod            = N_VProd_Trilinos;
  v->ops->nvdiv             = N_VDiv_Trilinos;
  v->ops->nvscale           = N_VScale_Trilinos;
  v->ops->nvabs             = N_VAbs_Trilinos;
  v->ops->nvinv             = N_VInv_Trilinos;
  v->ops->nvaddconst        = N_VAddConst_Trilinos;
  v->ops->nvdotprod         = N_VDotProd_Trilinos;
  v->ops->nvmaxnorm         = N_VMaxNorm_Trilinos;
  v->ops->nvwrmsnorm        = N_VWrmsNorm_Trilinos;
  v->ops->nvwrmsnormmask    = N_VWrmsNormMask_Trilinos;
  v->ops->nvmin             = N_VMin_Trilinos;
  v->ops->nvwl2norm         = N_VWL2Norm_Trilinos;
  v->ops->nvl1norm          = N_VL1Norm_Trilinos;
  v->ops->nvcompare         = N_VCompare_Trilinos;
  v->ops->nvinvtest         = N_VInvTest_Trilinos;
  v->ops->nvconstrmask      = N_VConstrMask_Trilinos;
  v->ops->nvminquotient     = N_VMinQuotient_Trilinos;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  v->ops->nvdotprodlocal     = N_VDotProdLocal_Trilinos;
  v->ops->nvmaxnormlocal     = N_VMaxNormLocal_Trilinos;
  v->ops->nvminlocal         = N_VMinLocal_Trilinos;
  v->ops->nvl1normlocal      = N_VL1NormLocal_Trilinos;
  v->ops->nvinvtestlocal     = N_VInvTestLocal_Trilinos;
  v->ops->nvconstrmasklocal  = N_VConstrMaskLocal_Trilinos;
  v->ops->nvminquotientlocal = N_VMinQuotientLocal_Trilinos;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Trilinos;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Trilinos;

  return(v);
}



/* ----------------------------------------------------------------
 * Function to create an N_Vector attachment to Tpetra vector.
 * void* argument is to allow for calling this method from C code.
 *
 */

N_Vector N_VMake_Trilinos(Teuchos::RCP<vector_type> vec)
{
  N_Vector v = NULL;

  // Create an N_Vector with operators attached and empty content
  v = N_VNewEmpty_Trilinos();
  if (v == NULL) return(NULL);

  // Create vector content using a pointer to Tpetra vector
  v->content = new Sundials::TpetraVectorInterface(vec);
  if (v->content == NULL) { N_VDestroy(v); return NULL; }

  return(v);
}


/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Trilinos(N_Vector w)
{
  N_Vector v;

  if (w == NULL) return(NULL);

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty();
  if (v == NULL) return(NULL);

  /* Attach operations */
  if (N_VCopyOps(w, v)) { N_VDestroy(v); return(NULL); }

  return(v);
}

N_Vector N_VClone_Trilinos(N_Vector w)
{
  N_Vector v = N_VCloneEmpty_Trilinos(w);
  if (v == NULL) return(NULL);

  // Get raw pointer to Tpetra vector
  Teuchos::RCP<vector_type> wvec = N_VGetVector_Trilinos(w);

  // Clone wvec and get raw pointer to the clone
  Teuchos::RCP<vector_type> tvec =
    Teuchos::rcp(new vector_type(*wvec, Teuchos::Copy));

  // Create vector content using the raw pointer to the cloned Tpetra vector
  v->content = new Sundials::TpetraVectorInterface(tvec);
  if (v->content == NULL) { N_VDestroy(v); return NULL; }

  return(v);
}

void N_VDestroy_Trilinos(N_Vector v)
{
  if (v == NULL) return;

  if(v->content != NULL) {
    Sundials::TpetraVectorInterface* iface =
      reinterpret_cast<Sundials::TpetraVectorInterface*>(v->content);

      // iface was created with 'new', so use 'delete' to destroy it.
      delete iface;
      v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL) { free(v->ops); v->ops = NULL; }
  free(v); v = NULL;

  return;
}

void N_VSpace_Trilinos(N_Vector x, sunindextype *lrw, sunindextype *liw)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm = xv->getMap()->getComm();
  int npes = comm->getSize();

  *lrw = (sunindextype)(xv->getGlobalLength());
  *liw = 2*npes;
}

/*
 * MPI communicator accessor
 */
void *N_VGetCommunicator_Trilinos(N_Vector x)
{
  using namespace Sundials;

#ifdef SUNDIALS_TRILINOS_HAVE_MPI
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  /* Access Teuchos::Comm* (which is actually a Teuchos::MpiComm*) */
  auto comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int>>(xv->getMap()->getComm());

  return((void*) comm->getRawMpiComm().get());   /* extract raw pointer to MPI_Comm */
#else
  return(NULL);
#endif
}

/*
 * Global vector length accessor
 */
sunindextype N_VGetLength_Trilinos(N_Vector x)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return ((sunindextype) xv->getGlobalLength());
}

/*
 * Linear combination of two vectors: z = a*x + b*y
 */
void N_VLinearSum_Trilinos(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  if (x == z) {
    zv->update(b, *yv, a);
  } else if (y == z) {
    zv->update(a, *xv, b);
  } else {
    zv->update(a, *xv, b, *yv, ZERO);
  }

}

/*
 * Set all vector elements to a constant: z[i] = c
 */
void N_VConst_Trilinos(realtype c, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<vector_type> zv = N_VGetVector_Trilinos(z);

  zv->putScalar(c);
}

/*
 * Elementwise multiply vectors: z[i] = x[i]*y[i]
 */
void N_VProd_Trilinos(N_Vector x, N_Vector y, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  zv->elementWiseMultiply(ONE, *xv, *yv, ZERO);
}

/*
 * Elementwise divide vectors: z[i] = x[i]/y[i]
 */
void N_VDiv_Trilinos(N_Vector x, N_Vector y, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  TpetraVector::elementWiseDivide(*xv, *yv, *zv);
}

/*
 * Scale vector: z = c*x
 */
void N_VScale_Trilinos(realtype c, N_Vector x, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  zv->scale(c, *xv);
}

/*
 * Elementwise absolute value: z[i] = |x[i]|
 */
void N_VAbs_Trilinos(N_Vector x, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  zv->abs(*xv);
}

/*
 * Elementwise inverse: z[i] = 1/x[i]
 */
void N_VInv_Trilinos(N_Vector x, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  zv->reciprocal(*xv);
}

/*
 * Add constant: z = x + b
 */
void N_VAddConst_Trilinos(N_Vector x, realtype b, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  TpetraVector::addConst(*xv, b, *zv);
}

/*
 * Scalar product of vectors x and y
 */
realtype N_VDotProd_Trilinos(N_Vector x, N_Vector y)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);

  return xv->dot(*yv);
}

/*
 * Max norm (L infinity) of vector x
 */
realtype N_VMaxNorm_Trilinos(N_Vector x)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return xv->normInf();
}

/*
 * Weighted RMS norm
 */
realtype N_VWrmsNorm_Trilinos(N_Vector x, N_Vector w)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv = N_VGetVector_Trilinos(w);

  return TpetraVector::normWrms(*xv, *wv);
}

/*
 * Masked weighted RMS norm
 */
realtype N_VWrmsNormMask_Trilinos(N_Vector x, N_Vector w, N_Vector id)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv  = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv  = N_VGetVector_Trilinos(w);
  Teuchos::RCP<const vector_type> idv = N_VGetVector_Trilinos(id);

  return TpetraVector::normWrmsMask(*xv, *wv, *idv);
}

/*
 * Returns minimum vector element
 */
realtype N_VMin_Trilinos(N_Vector x)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv  = N_VGetVector_Trilinos(x);

  return TpetraVector::minElement(*xv);
}

/*
 * Weighted L2 norm
 */
realtype N_VWL2Norm_Trilinos(N_Vector x, N_Vector w)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv = N_VGetVector_Trilinos(w);

  return TpetraVector::normWL2(*xv, *wv);
}

/*
 * L1 norm
 */
realtype N_VL1Norm_Trilinos(N_Vector x)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return xv->norm1();
}

/*
 * Elementwise z[i] = |x[i]| >= c ? 1 : 0
 */
void N_VCompare_Trilinos(realtype c, N_Vector x, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  TpetraVector::compare(c, *xv, *zv);
}

/*
 * Elementwise inverse with zero checking: z[i] = 1/x[i], x[i] != 0
 */
booleantype N_VInvTest_Trilinos(N_Vector x, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv       = N_VGetVector_Trilinos(z);

  return TpetraVector::invTest(*xv, *zv) ? SUNTRUE : SUNFALSE;
}

/*
 * Checks constraint violations for vector x. Constraints are defined in
 * vector c, and constraint violation flags are stored in vector m.
 */
booleantype N_VConstrMask_Trilinos(N_Vector c, N_Vector x, N_Vector m)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> cv = N_VGetVector_Trilinos(c);
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> mv       = N_VGetVector_Trilinos(m);

  return TpetraVector::constraintMask(*cv, *xv, *mv) ? SUNTRUE : SUNFALSE;
}

/*
 * Find minimum quotient: minq  = min ( num[i]/denom[i]), denom[i] != 0.
 */
realtype N_VMinQuotient_Trilinos(N_Vector num, N_Vector denom)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> numv = N_VGetVector_Trilinos(num);
  Teuchos::RCP<const vector_type> denv = N_VGetVector_Trilinos(denom);

  return TpetraVector::minQuotient(*numv, *denv);
}

/*
 * MPI task-local dot product
 */
realtype N_VDotProdLocal_Trilinos(N_Vector x, N_Vector y)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> yv = N_VGetVector_Trilinos(y);

  return TpetraVector::dotProdLocal(*xv, *yv);
}

/*
 * MPI task-local maximum norm
 */
realtype N_VMaxNormLocal_Trilinos(N_Vector x)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return TpetraVector::maxNormLocal(*xv);
}

/*
 * MPI task-local minimum element
 */
realtype N_VMinLocal_Trilinos(N_Vector x)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return TpetraVector::minLocal(*xv);
}

/*
 * MPI task-local L1 norm
 */
realtype N_VL1NormLocal_Trilinos(N_Vector x)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);

  return TpetraVector::L1NormLocal(*xv);
}

/*
 * MPI task-local weighted squared sum
 */
realtype N_VWSqrSumLocal_Trilinos(N_Vector x, N_Vector w)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv = N_VGetVector_Trilinos(w);

  return TpetraVector::WSqrSumLocal(*xv, *wv);
}

/*
 * MPI task-local weighted masked squared sum
 */
realtype N_VWSqrSumMaskLocal_Trilinos(N_Vector x, N_Vector w, N_Vector id)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv  = N_VGetVector_Trilinos(x);
  Teuchos::RCP<const vector_type> wv  = N_VGetVector_Trilinos(w);
  Teuchos::RCP<const vector_type> idv = N_VGetVector_Trilinos(id);

  return TpetraVector::WSqrSumMaskLocal(*xv, *wv, *idv);
}

/*
 * MPI task-local elementwise inverse with zero checking: z[i] = 1/x[i], x[i] != 0
 */
booleantype N_VInvTestLocal_Trilinos(N_Vector x, N_Vector z)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> zv = N_VGetVector_Trilinos(z);

  return TpetraVector::invTestLocal(*xv, *zv) ? SUNTRUE : SUNFALSE;
}

/*
 * MPI task-local constraint checking for vector x. Constraints are defined in
 * vector c, and constraint violation flags are stored in vector m.
 */
booleantype N_VConstrMaskLocal_Trilinos(N_Vector c, N_Vector x, N_Vector m)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> cv = N_VGetVector_Trilinos(c);
  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  Teuchos::RCP<vector_type> mv       = N_VGetVector_Trilinos(m);

  return TpetraVector::constraintMaskLocal(*cv, *xv, *mv) ? SUNTRUE : SUNFALSE;
}

/*
 * MPI task-local minimum quotient: minq  = min ( num[i]/denom[i]), denom[i] != 0.
 */
realtype N_VMinQuotientLocal_Trilinos(N_Vector num, N_Vector denom)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> numv = N_VGetVector_Trilinos(num);
  Teuchos::RCP<const vector_type> denv = N_VGetVector_Trilinos(denom);

  return TpetraVector::minQuotientLocal(*numv, *denv);
}
