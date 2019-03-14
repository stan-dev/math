/* -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 *
 * Based on N_Vector_Parallel by Scott D. Cohen, Alan C. Hindmarsh,
 * Radu Serban, and Aaron Collier @ LLNL
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
  N_Vector_Ops ops;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvgetvectorid     = N_VGetVectorID_Trilinos;
  ops->nvclone           = N_VClone_Trilinos;
  ops->nvcloneempty      = N_VCloneEmpty_Trilinos;
  ops->nvdestroy         = N_VDestroy_Trilinos;
  ops->nvspace           = N_VSpace_Trilinos;
  ops->nvgetarraypointer = NULL;
  ops->nvsetarraypointer = NULL;
  ops->nvlinearsum       = N_VLinearSum_Trilinos;
  ops->nvconst           = N_VConst_Trilinos;
  ops->nvprod            = N_VProd_Trilinos;
  ops->nvdiv             = N_VDiv_Trilinos;
  ops->nvscale           = N_VScale_Trilinos;
  ops->nvabs             = N_VAbs_Trilinos;
  ops->nvinv             = N_VInv_Trilinos;
  ops->nvaddconst        = N_VAddConst_Trilinos;
  ops->nvdotprod         = N_VDotProd_Trilinos;
  ops->nvmaxnorm         = N_VMaxNorm_Trilinos;
  ops->nvwrmsnorm        = N_VWrmsNorm_Trilinos;
  ops->nvwrmsnormmask    = N_VWrmsNormMask_Trilinos;
  ops->nvmin             = N_VMin_Trilinos;
  ops->nvwl2norm         = N_VWL2Norm_Trilinos;
  ops->nvl1norm          = N_VL1Norm_Trilinos;
  ops->nvcompare         = N_VCompare_Trilinos;
  ops->nvinvtest         = N_VInvTest_Trilinos;
  ops->nvconstrmask      = N_VConstrMask_Trilinos;
  ops->nvminquotient     = N_VMinQuotient_Trilinos;

  /* fused vector operations */
  ops->nvlinearcombination = NULL;
  ops->nvscaleaddmulti     = NULL;
  ops->nvdotprodmulti      = NULL;

  /* vector array operations */
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
  if (v == NULL)
     return(NULL);

  // Create vector content using a pointer to Tpetra vector
  v->content = new Sundials::TpetraVectorInterface(vec);
  if (v->content == NULL) {
    free(v->ops);
    free(v);
    return NULL;
  }

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
  N_Vector_Ops ops;

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
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvwrmsnormmask    = w->ops->nvwrmsnormmask;
  ops->nvmin             = w->ops->nvmin;
  ops->nvwl2norm         = w->ops->nvwl2norm;
  ops->nvl1norm          = w->ops->nvl1norm;
  ops->nvcompare         = w->ops->nvcompare;
  ops->nvinvtest         = w->ops->nvinvtest;
  ops->nvconstrmask      = w->ops->nvconstrmask;
  ops->nvminquotient     = w->ops->nvminquotient;

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

  /* Attach ops and set content to NULL */
  v->content = NULL;
  v->ops     = ops;

  return(v);
}

N_Vector N_VClone_Trilinos(N_Vector w)
{
  N_Vector v = N_VCloneEmpty_Trilinos(w);
  if (v == NULL)
    return(NULL);

  // Get raw pointer to Tpetra vector
  Teuchos::RCP<vector_type> wvec = N_VGetVector_Trilinos(w);

  // Clone wvec and get raw pointer to the clone
  Teuchos::RCP<vector_type> tvec =
    Teuchos::rcp(new vector_type(*wvec, Teuchos::Copy));

  // Create vector content using the raw pointer to the cloned Tpetra vector
  v->content = new Sundials::TpetraVectorInterface(tvec);
  if (v->content == NULL) {
    free(v->ops);
    free(v);
    return NULL;
  }

  return(v);
}

void N_VDestroy_Trilinos(N_Vector v)
{
  if(v->content != NULL) {
    Sundials::TpetraVectorInterface* iface =
      reinterpret_cast<Sundials::TpetraVectorInterface*>(v->content);

      // iface was created with 'new', so use 'delete' to destroy it.
      delete iface;
      v->content = NULL;
  }

  free(v->ops);
  v->ops = NULL;

  free(v);
  v = NULL;
}

void N_VSpace_Trilinos(N_Vector x, sunindextype *lrw, sunindextype *liw)
{
  using namespace Sundials;

  Teuchos::RCP<const vector_type> xv = N_VGetVector_Trilinos(x);
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm = xv->getMap()->getComm();
  int npes = comm->getSize();

  *lrw = xv->getGlobalLength();
  *liw = 2*npes;
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
 * vector c, and constrain violation flags are stored in vector m.
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
