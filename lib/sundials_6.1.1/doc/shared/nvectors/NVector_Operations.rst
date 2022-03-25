..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors.Ops:

Description of the NVECTOR operations
=====================================

.. _NVectors.Ops.Standard:

Standard vector operations
--------------------------

The standard vector operations defined by the generic ``N_Vector``
module are defined as follows.  For each of these operations, we give
the name, usage of the function, and a description of its mathematical
operations below.


.. c:function:: N_Vector_ID N_VGetVectorID(N_Vector w)

   Returns the vector type identifier for the vector ``w``.  It is
   used to determine the vector implementation type (e.g. serial,
   parallel, ...) from the abstract ``N_Vector`` interface.  Returned
   values are given in :numref:`NVectors.Description.vectorIDs`.

   Usage:

   .. code-block:: c

      id = N_VGetVectorID(w);


.. c:function:: N_Vector N_VClone(N_Vector w)

   Creates a new ``N_Vector`` of the same type as an existing vector
   *w* and sets the *ops* field.  It does not copy the vector, but
   rather allocates storage for the new vector.

   Usage:

   .. code-block:: c

      v = N_VClone(w);


.. c:function:: N_Vector N_VCloneEmpty(N_Vector w)

   Creates a new ``N_Vector`` of the same type as an existing vector
   *w* and sets the *ops* field.  It does not allocate storage for the
   new vector's data.

   Usage:

   .. code-block:: c

      v = N VCloneEmpty(w);


.. c:function:: void N_VDestroy(N_Vector v)

   Destroys the ``N_Vector`` *v* and frees memory allocated for its
   internal data.

   Usage:

   .. code-block:: c

      N_VDestroy(v);


.. c:function:: void N_VSpace(N_Vector v, sunindextype* lrw, sunindextype* liw)

   Returns storage requirements for the ``N_Vector`` *v*:

   * *lrw* contains the number of ``realtype`` words
   * *liw* contains the number of integer words.

   This function is advisory only, for use in
   determining a user's total space requirements; it could be a dummy
   function in a user-supplied NVECTOR module if that information is
   not of interest.

   Usage:

   .. code-block:: c

      N_VSpace(nvSpec, &lrw, &liw);


.. c:function:: realtype* N_VGetArrayPointer(N_Vector v)

   Returns a pointer to a ``realtype`` array from the ``N_Vector``
   *v*.  Note that this assumes that the internal data in the
   ``N_Vector`` is a contiguous array of ``realtype`` and is
   accesible from the CPU.

   This routine is
   only used in the solver-specific interfaces to the dense and banded
   (serial) linear solvers, and in the interfaces to the banded
   (serial) and band-block-diagonal (parallel) preconditioner modules
   provided with SUNDIALS.

   Usage:

   .. code-block:: c

      vdata = N_VGetArrayPointer(v);


.. c:function:: realtype* N_VGetDeviceArrayPointer(N_Vector v)

   Returns a device pointer to a ``realtype`` array from the ``N_Vector``
   ``v``. Note that this assumes that the internal data in ``N_Vector`` is a
   contiguous array of ``realtype`` and is accessible from the device (e.g.,
   GPU).

   This operation is *optional* except when using the GPU-enabled direct
   linear solvers.

   Usage:

   .. code-block:: c

      vdata = N_VGetArrayPointer(v);


.. c:function:: void N_VSetArrayPointer(realtype* vdata, N_Vector v)

   Replaces the data array pointer in an ``N_Vector`` with a given
   array of ``realtype``.  Note that this assumes that the internal
   data in the ``N_Vector`` is a contiguous array of
   ``realtype``. This routine is only used in the interfaces to the
   dense (serial) linear solver, hence need not exist in a
   user-supplied NVECTOR module.

   Usage:

   .. code-block:: c

      N_VSetArrayPointer(vdata,v);


.. c:function:: void* N_VGetCommunicator(N_Vector v)

   Returns a pointer to the ``MPI_Comm`` object associated with the
   vector (if applicable).  For MPI-unaware vector implementations, this
   should return ``NULL``.

   Usage:

   .. code-block:: c

      commptr = N_VGetCommunicator(v);


.. c:function:: sunindextype N_VGetLength(N_Vector v)

   Returns the global length (number of "active" entries) in the
   NVECTOR *v*.  This value should be cumulative across all processes
   if the vector is used in a parallel environment.  If *v*
   contains additional storage, e.g., for parallel communication, those
   entries should not be included.

   Usage:

   .. code-block:: c

      global_length = N_VGetLength(v);


.. c:function:: void N_VLinearSum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z)

   Performs the operation *z = ax + by*, where *a* and *b* are
   ``realtype`` scalars and *x* and *y* are of type ``N_Vector``:

   .. math::
      z_i = a x_i + b y_i, \quad i=0,\ldots,n-1.

   The output vector *z* can be the same as either of the input vectors (*x* or *y*).


   Usage:

   .. code-block:: c

      N_VLinearSum(a, x, b, y, z);


.. c:function:: void N_VConst(realtype c, N_Vector z)

   Sets all components of the ``N_Vector`` *z* to ``realtype`` *c*:

   .. math::
      z_i = c, \quad i=0,\ldots,n-1.

   Usage:

   .. code-block:: c

      N_VConst(c, z);


.. c:function:: void N_VProd(N_Vector x, N_Vector y, N_Vector z)

   Sets the ``N_Vector`` *z* to be the component-wise product of the
   ``N_Vector`` inputs *x* and *y*:

   .. math::
      z_i = x_i y_i, \quad i=0,\ldots,n-1.

   Usage:

   .. code-block:: c

      N_VProd(x, y, z);


.. c:function:: void N_VDiv(N_Vector x, N_Vector y, N_Vector z)

   Sets the ``N_Vector`` *z* to be the component-wise ratio of the
   ``N_Vector`` inputs *x* and *y*:

   .. math::
      z_i = \frac{x_i}{y_i}, \quad i=0,\ldots,n-1.

   The :math:`y_i` may not be tested for 0 values. It should only be
   called with a *y* that is guaranteed to have all nonzero components.

   Usage:

   .. code-block:: c

      N_VDiv(x, y, z);


.. c:function:: void N_VScale(realtype c, N_Vector x, N_Vector z)

   Scales the ``N_Vector`` *x* by the ``realtype`` scalar *c* and
   returns the result in *z*:

   .. math::
      z_i = c x_i, \quad i=0,\ldots,n-1.

   Usage:

   .. code-block:: c

      N_VScale(c, x, z);


.. c:function:: void N_VAbs(N_Vector x, N_Vector z)

   Sets the components of the ``N_Vector`` *z* to be the absolute
   values of the components of the ``N_Vector`` *x*:

   .. math::
      z_i = |x_i|, \quad i=0,\ldots,n-1.

   Usage:

   .. code-block:: c

      N_VAbs(x, z);


.. c:function:: void N_VInv(N_Vector x, N_Vector z)

   Sets the components of the ``N_Vector`` *z* to be the inverses of
   the components of the ``N_Vector`` *x*:

   .. math::
      z_i = \frac{1}{x_i}, \quad i=0,\ldots,n-1.

   This routine may not check for division by 0.  It should be called
   only with an *x* which is guaranteed to have all nonzero components.

   Usage:

   .. code-block:: c

      N_VInv(x, z);


.. c:function:: void N_VAddConst(N_Vector x, realtype b, N_Vector z)

   Adds the ``realtype`` scalar *b* to all components of *x* and
   returns the result in the ``N_Vector`` *z*:

   .. math::
      z_i = x_i+b, \quad i=0,\ldots,n-1.

   Usage:

   .. code-block:: c

      N_VAddConst(x, b, z);


.. c:function:: realtype N_VDotProd(N_Vector x, N_Vector z)

   Returns the value of the dot-product of the ``N_Vectors`` *x* and *y*:

   .. math::
      d = \sum_{i=0}^{n-1} x_i y_i.

   Usage:

   .. code-block:: c

      d = N_VDotProd(x, y);


.. c:function:: realtype N_VMaxNorm(N_Vector x)

   Returns the value of the :math:`l_{\infty}` norm of the
   ``N_Vector`` *x*:

   .. math::
      m = \max_{0\le i< n} |x_i|.

   Usage:

   .. code-block:: c

      m = N_VMaxNorm(x);


.. c:function:: realtype N_VWrmsNorm(N_Vector x, N_Vector w)

   Returns the weighted root-mean-square norm of the ``N_Vector`` *x*
   with (positive) ``realtype`` weight vector *w*:

   .. math::
      m = \sqrt{\left( \sum_{i=0}^{n-1} (x_i w_i)^2 \right) / n}

   Usage:

   .. code-block:: c

      m = N_VWrmsNorm(x, w);


.. c:function:: realtype N_VWrmsNormMask(N_Vector x, N_Vector w, N_Vector id)

   Returns the weighted root mean square norm of the ``N_Vector`` *x*
   with ``realtype`` weight vector *w* built using only the
   elements of *x* corresponding to positive elements of the
   ``N_Vector`` *id*:

   .. math::
      m = \sqrt{\left( \sum_{i=0}^{n-1} (x_i w_i H(id_i))^2 \right) / n},

   where :math:`H(\alpha)=\begin{cases} 1 & \alpha>0\\ 0 & \alpha \leq 0\end{cases}`.

   Usage:

   .. code-block:: c

      m = N_VWrmsNormMask(x, w, id);

.. c:function:: realtype N_VMin(N_Vector x)

   Returns the smallest element of the ``N_Vector`` *x*:

   .. math::
      m = \min_{0\le i< n} x_i.

   Usage:

   .. code-block:: c

      m = N_VMin(x);

.. c:function:: realtype N_VWl2Norm(N_Vector x, N_Vector w)

   Returns the weighted Euclidean :math:`l_2` norm of the ``N_Vector``
   *x* with ``realtype`` weight vector *w*:

   .. math::
      m = \sqrt{\sum_{i=0}^{n-1}\left(x_i w_i\right)^2}.

   Usage:

   .. code-block:: c

      m = N_VWL2Norm(x, w);

.. c:function:: realtype N_VL1Norm(N_Vector x)

   Returns the :math:`l_1` norm of the ``N_Vector`` *x*:

   .. math::
      m = \sum_{i=0}^{n-1} |x_i|.

   Usage:

   .. code-block:: c

      m = N_VL1Norm(x);


.. c:function:: void N_VCompare(realtype c, N_Vector x, N_Vector z)

   Compares the components of the ``N_Vector`` *x* to the ``realtype``
   scalar *c* and returns an ``N_Vector`` *z* such that for all
   :math:`0\le i< n`,

   .. math::
      z_i = \begin{cases} 1.0 &\quad\text{if}\; |x_i| \ge c,\\
                          0.0 &\quad\text{otherwise}\end{cases}.

   Usage:

   .. code-block:: c

      N_VCompare(c, x, z);

.. c:function:: booleantype N_VInvTest(N_Vector x, N_Vector z)

   Sets the components of the ``N_Vector`` *z* to be the inverses of
   the components of the ``N_Vector`` *x*, with prior testing for
   zero values:

   .. math::
      z_i = \frac{1}{x_i}, \quad i=0,\ldots,n-1.

   This routine returns a boolean assigned to ``SUNTRUE`` if all
   components of *x* are nonzero (successful inversion) and returns
   ``SUNFALSE`` otherwise.

   Usage:

   .. code-block:: c

      t = N_VInvTest(x, z);

.. c:function:: booleantype N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)

   Performs the following constraint tests based on the values in
   :math:`c_i`:

   .. math::
      \begin{array}{rllll}
      x_i &>& 0 \;&\text{if}\; &c_i = 2, \\
      x_i &\ge& 0 \;&\text{if}\; &c_i = 1, \\
      x_i &<& 0 \;&\text{if}\; &c_i = -2, \\
      x_i &\le& 0 \;&\text{if}\; &c_i = -1.
      \end{array}

   There is no constraint on :math:`x_i` if :math:`c_i = 0`. This
   routine returns a boolean assigned to ``SUNFALSE`` if any element
   failed the constraint test and assigned to ``SUNTRUE`` if all
   passed. It also sets a mask vector *m*, with elements equal to 1.0
   where the constraint test failed, and 0.0 where the test
   passed. This routine is used only for constraint checking.

   Usage:

   .. code-block:: c

      t = N_VConstrMask(c, x, m);

.. c:function:: realtype N_VMinQuotient(N_Vector num, N_Vector denom)

   This routine returns the minimum of the quotients obtained by
   termwise dividing the elements of *n* by the elements in *d*:

   .. math::
      \min_{0\le i< n} \frac{\text{num}_i}{\text{denom}_i}.

   A zero element in *denom* will be skipped.  If no such quotients
   are found, then the large value ``BIG_REAL`` (defined in the header
   file ``sundials_types.h``) is returned.

   Usage:

   .. code-block:: c

      minq = N_VMinQuotient(num, denom);



.. _NVectors.Ops.Fused:

Fused operations
----------------

The following fused vector operations are *optional*. These
operations are intended to increase data reuse, reduce parallel
communication on distributed memory systems, and lower the number of
kernel launches on systems with accelerators. If a particular NVECTOR
implementation defines one of the fused vector operations as
``NULL``, the NVECTOR interface will call one of the above standard
vector operations as necessary.  As above, for each operation, we give
the name, usage of the function, and a description of its mathematical
operations below.


.. c:function:: int N_VLinearCombination(int nv, realtype* c, N_Vector* X, N_Vector z)

   This routine computes the linear combination of *nv* vectors with :math:`n` elements:

   .. math::
      z_i = \sum_{j=0}^{nv-1} c_j x_{j,i}, \quad i=0,\ldots,n-1,

   where :math:`c` is an array of :math:`nv` scalars, :math:`x_j` is a
   vector in the vector array *X*, and *z* is the output
   vector. If the output vector *z* is one of the vectors in *X*, then
   it *must* be the first vector in the vector array. The operation
   returns 0 for success and a non-zero value otherwise.

   Usage:

   .. code-block:: c

      retval = N_VLinearCombination(nv, c, X, z);


.. c:function:: int N_VScaleAddMulti(int nv, realtype* c, N_Vector x, N_Vector* Y, N_Vector* Z)

   This routine scales and adds one vector to *nv* vectors with :math:`n` elements:

   .. math::
      z_{j,i} = c_j x_i + y_{j,i}, \quad j=0,\ldots,nv-1 \quad i=0,\ldots,n-1,

   where *c* is an array of scalars, *x* is a vector, :math:`y_j` is a
   vector in the vector array *Y*, and :math:`z_j` is an output vector
   in the vector array *Z*. The operation returns 0 for success and a
   non-zero value otherwise.

   Usage:

   .. code-block:: c

      retval = N_VScaleAddMulti(nv, c, x, Y, Z);


.. c:function:: int N_VDotProdMulti(int nv, N_Vector x, N_Vector* Y, realtype* d)

   This routine computes the dot product of a vector with *nv* vectors
   having :math:`n` elements:

   .. math::
      d_j = \sum_{i=0}^{n-1} x_i y_{j,i}, \quad j=0,\ldots,nv-1,

   where *d* is an array of scalars containing the computed dot
   products, *x* is a vector, and :math:`y_j` is a vector the vector
   array *Y*. The operation returns 0 for success and a non-zero value
   otherwise.

   Usage:

   .. code-block:: c

      retval = N_VDotProdMulti(nv, x, Y, d);


.. _NVectors.Ops.Array:

Vector array operations
-----------------------

The following vector array operations are also *optional*. As with the
fused vector operations, these are intended to increase data reuse,
reduce parallel communication on distributed memory systems, and lower
the number of kernel launches on systems with accelerators. If a
particular NVECTOR implementation defines one of the fused or vector
array operations as ``NULL``, the NVECTOR interface will call one of
the above standard vector operations as necessary.  As above, for each
operation, we give the name, usage of the function, and a description
of its mathematical operations below.


.. c:function:: int N_VLinearSumVectorArray(int nv, realtype a, N_Vector X, realtype b, N_Vector* Y, N_Vector* Z)

   This routine computes the linear sum of two vector arrays of *nv* vectors with :math:`n` elements:

   .. math::
      z_{j,i} = a x_{j,i} + b y_{j,i}, \quad i=0,\ldots,n-1 \quad j=0,\ldots,nv-1,

   where *a* and *b* are scalars, :math:`x_j` and :math:`y_j` are
   vectors in the vector arrays *X* and *Y* respectively, and
   :math:`z_j` is a vector in the output vector array *Z*. The
   operation returns 0 for success and a non-zero value otherwise.

   Usage:

   .. code-block:: c

      retval = N_VLinearSumVectorArray(nv, a, X, b, Y, Z);


.. c:function:: int N_VScaleVectorArray(int nv, realtype* c, N_Vector* X, N_Vector* Z)

   This routine scales each element in a vector of :math:`n` elements
   in a vector array of *nv* vectors by a potentially different constant:

   .. math::
      z_{j,i} = c_j x_{j,i}, \quad i=0,\ldots,n-1 \quad j=0,\ldots,nv-1,

   where *c* is an array of scalars, :math:`x_j` is a vector in the
   vector array *X*, and :math:`z_j` is a vector in the output vector
   array *Z*. The operation returns 0 for success and a non-zero value otherwise.

   Usage:

   .. code-block:: c

      retval = N_VScaleVectorArray(nv, c, X, Z);


.. c:function:: int N_VConstVectorArray(int nv, realtype c, N_Vector* Z)

   This routine sets each element in a vector of :math:`n` elements in
   a vector array of *nv* vectors to the same value:

   .. math::
      z_{j,i} = c, \quad i=0,\ldots,n-1 \quad j=0,\ldots,nv-1,

   where *c* is a scalar and :math:`z_j` is a vector in the vector
   array *Z*. The operation returns 0 for success and a non-zero value otherwise.

   Usage:

   .. code-block:: c

      retval = N_VConstVectorArray(nv, c, Z);


.. c:function:: int N_VWrmsNormVectorArray(int nv, N_Vector* X, N_Vector* W, realtype* m)

   This routine computes the weighted root mean square norm of each
   vector in a vector array:

   .. math::
      m_j = \left( \frac1n \sum_{i=0}^{n-1} \left(x_{j,i} w_{j,i}\right)^2\right)^{1/2}, \quad j=0,\ldots,nv-1,

   where :math:`x_j` is a vector in the vector array *X*, :math:`w_j`
   is a weight vector in the vector array *W*, and *m* is the output
   array of scalars containing the computed norms. The operation
   returns 0 for success and a non-zero value otherwise.

   Usage:

   .. code-block:: c

      retval = N_VWrmsNormVectorArray(nv, X, W, m);


.. c:function:: int N_VWrmsNormMaskVectorArray(int nv, N_Vector* X, N_Vector* W, N_Vector id, realtype* m)

   This routine computes the masked weighted root mean square norm of
   each vector in a vector array:

   .. math::
      m_j = \left( \frac1n \sum_{i=0}^{n-1} \left(x_{j,i} w_{j,i} H(id_i)\right)^2 \right)^{1/2}, \quad j=0,\ldots,nv-1,

   where :math:`H(id_i)=1` if :math:`id_i > 0` and is zero otherwise,
   :math:`x_j` is a vector in the vector array *X*, :math:`w_j` is a
   weight vector in the vector array *W*, *id* is the mask vector, and
   *m* is the output array of scalars containing the computed
   norms. The operation returns 0 for success and a non-zero value
   otherwise.

   Usage:

   .. code-block:: c

      retval = N_VWrmsNormMaskVectorArray(nv, X, W, id, m);


.. c:function:: int N_VScaleAddMultiVectorArray(int nv, int nsum, realtype* c, N_Vector* X, N_Vector** YY, N_Vector** ZZ)

   This routine scales and adds a vector array of *nv* vectors to
   *nsum* other vector arrays:

   .. math::
      z_{k,j,i} = c_k x_{j,i} + y_{k,j,i}, \quad i=0,\ldots,n-1 \quad j=0,\ldots,nv-1, \quad k=0,\ldots,nsum-1

   where *c* is an array of scalars, :math:`x_j` is a vector in the
   vector array *X*, :math:`y_{k,j}` is a vector in the array of
   vector arrays *YY*, and :math:`z_{k,j}` is an output vector in the
   array of vector arrays *ZZ*. The operation returns 0 for success
   and a non-zero value otherwise.

   Usage:

   .. code-block:: c

      retval = N_VScaleAddMultiVectorArray(nv, nsum, c, x, YY, ZZ);


.. c:function:: int N_VLinearCombinationVectorArray(int nv, int nsum, realtype* c, N_Vector** XX, N_Vector* Z)

   This routine computes the linear combination of *nsum* vector
   arrays containing *nv* vectors:

   .. math::
      z_{j,i} = \sum_{k=0}^{nsum-1} c_k x_{k,j,i}, \quad i=0,\ldots,n-1 \quad j=0,\ldots,nv-1,

   where *c* is an array of scalars, :math:`x_{k,j}` is a vector in
   array of vector arrays *XX*, and :math:`z_{j,i}` is an output
   vector in the vector array *Z*. If the output vector array is one
   of the vector arrays in *XX*, it *must* be the first vector array
   in *XX*. The operation returns 0 for success and a non-zero value
   otherwise.

   Usage:

   .. code-block:: c

      retval = N_VLinearCombinationVectorArray(nv, nsum, c, XX, Z);


.. _NVectors.Ops.Local:

Local reduction operations
--------------------------

The following local reduction operations are also *optional*. As with
the fused and vector array operations, these are intended to reduce
parallel communication on distributed memory systems. If a particular
NVECTOR implementation defines one of the local reduction operations
as ``NULL``, the NVECTOR interface will call one of the above standard
vector operations as necessary.  As above, for each operation, we give
the name, usage of the function, and a description of its mathematical
operations below.


.. c:function:: realtype N_VDotProdLocal(N_Vector x, N_Vector y)

   This routine computes the MPI task-local portion of the ordinary
   dot product of *x* and *y*:

   .. math::
      d=\sum_{i=0}^{n_{local}-1} x_i y_i,

   where :math:`n_{local}` corresponds to the number of components in
   the vector on this MPI task (or :math:`n_{local}=n` for MPI-unaware
   applications).

   Usage:

   .. code-block:: c

      d = N_VDotProdLocal(x, y);


.. c:function:: realtype N_VMaxNormLocal(N_Vector x)

   This routine computes the MPI task-local portion of the maximum
   norm of the NVECTOR *x*:

   .. math::
      m = \max_{0\le i< n_{local}} | x_i |,

   where :math:`n_{local}` corresponds to the number of components in
   the vector on this MPI task (or :math:`n_{local}=n` for MPI-unaware
   applications).

   Usage:

   .. code-block:: c

      m = N_VMaxNormLocal(x);


.. c:function:: realtype N_VMinLocal(N_Vector x)

   This routine computes the smallest element of the MPI task-local
   portion of the NVECTOR *x*:

   .. math::
      m = \min_{0\le i< n_{local}} x_i,

   where :math:`n_{local}` corresponds to the number of components in
   the vector on this MPI task (or :math:`n_{local}=n` for MPI-unaware
   applications).

   Usage:

   .. code-block:: c

      m = N_VMinLocal(x);


.. c:function:: realtype N_VL1NormLocal(N_Vector x)

   This routine computes the MPI task-local portion of the :math:`l_1`
   norm of the ``N_Vector`` *x*:

   .. math::
      n = \sum_{i=0}^{n_{local}-1} | x_i |,

   where :math:`n_{local}` corresponds to the number of components in
   the vector on this MPI task (or :math:`n_{local}=n` for MPI-unaware
   applications).

   Usage:

   .. code-block:: c

      n = N_VL1NormLocal(x);


.. c:function:: realtype N_VWSqrSumLocal(N_Vector x, N_Vector w)

   This routine computes the MPI task-local portion of the weighted
   squared sum of the NVECTOR *x* with weight vector *w*:

   .. math::
      s = \sum_{i=0}^{n_{local}-1} (x_i w_i)^2,

   where :math:`n_{local}` corresponds to the number of components in
   the vector on this MPI task (or :math:`n_{local}=n` for MPI-unaware
   applications).

   Usage:

   .. code-block:: c

      s = N_VWSqrSumLocal(x, w);


.. c:function:: realtype N_VWSqrSumMaskLocal(N_Vector x, N_Vector w, N_Vector id)

   This routine computes the MPI task-local portion of the weighted
   squared sum of the NVECTOR *x* with weight vector *w* built using
   only the elements of *x* corresponding to positive elements of the NVECTOR *id*:

   .. math::
      m = \sum_{i=0}^{n_{local}-1} (x_i w_i H(id_i))^2,

   where

   .. math::
      H(\alpha) = \begin{cases} 1 & \alpha > 0 \\ 0 & \alpha \leq 0 \end{cases}

   and :math:`n_{local}` corresponds to the number of components in
   the vector on this MPI task (or :math:`n_{local}=n` for MPI-unaware
   applications).

   Usage:

   .. code-block:: c

      s = N_VWSqrSumMaskLocal(x, w, id);


.. c:function:: booleantype N_VInvTestLocal(N_Vector x)

   This routine sets the MPI task-local components of the
   NVECTOR *z* to be the inverses of the components of the NVECTOR
   *x*, with prior testing for zero values:

   .. math::
      z_i = \frac{1}{x_i}, \: i=0,\ldots,n_{local}-1

   where :math:`n_{local}` corresponds to the number of components in
   the vector on this MPI task (or :math:`n_{local}=n` for MPI-unaware
   applications).  This routine returns a boolean assigned to
   ``SUNTRUE`` if all task-local components of *x* are nonzero
   (successful inversion) and returns ``SUNFALSE`` otherwise.

   Usage:

   .. code-block:: c

      t = N_VInvTestLocal(x);


.. c:function:: booleantype N_VConstrMaskLocal(N_Vector c, N_Vector x, N_Vector m)

   Performs the following constraint tests based on the values in
   :math:`c_i`:

   .. math::
      \begin{array}{rllll}
      x_i &>& 0 \;&\text{if}\; &c_i = 2, \\
      x_i &\ge& 0 \;&\text{if}\; &c_i = 1, \\
      x_i &<& 0 \;&\text{if}\; &c_i = -2, \\
      x_i &\le& 0 \;&\text{if}\; &c_i = -1.
      \end{array}

   for all MPI task-local components of the vectors.
   This routine returns a boolean assigned to ``SUNFALSE`` if any
   task-local element failed the constraint test and assigned to
   ``SUNTRUE`` if all passed.  It also sets a mask vector *m*, with
   elements equal to 1.0 where the constraint test failed, and 0.0
   where the test passed.  This routine is used only for constraint
   checking.

   Usage:

   .. code-block:: c

      t = N_VConstrMaskLocal(c, x, m);


.. c:function:: realtype N_VMinQuotientLocal(N_Vector num, N_Vector denom)

   This routine returns the minimum of the quotients obtained by
   term-wise dividing :math:`num_i` by :math:`denom_i`, for all MPI
   task-local components of the vectors.  A zero element in *denom*
   will be skipped. If no such quotients are found, then the large value
   ``BIG_REAL`` (defined in the header file ``sundials_types.h``)
   is returned.

   Usage:

   .. code-block:: c

      minq = N_VMinQuotientLocal(num, denom);


.. _NVectors.Ops.SingleBufferReduction:

Single Buffer Reduction Operations
----------------------------------

The following *optional* operations are used to combine separate reductions into
a single MPI call by splitting the local computation and communication into
separate functions. These operations are used in low-synchronization
orthogonalization methods to reduce the number of MPI ``Allreduce`` calls. If a
particular NVECTOR implementation does not define these operations additional
communication will be required.

.. c:function:: int N_VDotProdMultiLocal(int nv, N_Vector x, N_Vector* Y, realtype* d)

   This routine computes the MPI task-local portion of the dot product of a
   vector :math:`x` with *nv* vectors :math:`y_j`:

   .. math::
      d_j = \sum_{i=0}^{n_{local}-1} x_i y_{j,i}, \quad j=0,\ldots,nv-1,

   where :math:`d` is an array of scalars containing the computed dot products,
   :math:`x` is a vector, :math:`y_j` is a vector in the vector array *Y*, and
   :math:`n_{local}` corresponds to the number of components in the vector on
   this MPI task. The operation returns 0 for success and a non-zero value
   otherwise.

   Usage:

   .. code-block:: c

      retval = N_VDotProdMultiLocal(nv, x, Y, d);


.. c:function:: int N_VDotProdMultiAllReduce(int nv, N_Vector x, realtype* d)

   This routine combines the MPI task-local portions of the dot product of a
   vector :math:`x` with *nv* vectors:

   .. code-block:: c

      retval = MPI_Allreduce(MPI_IN_PLACE, d, nv, MPI_SUNREALTYPE, MPI_SUM, comm)

   where *d* is an array of *nv* scalars containing the local contributions to
   the dot product and *comm* is the MPI communicator associated with the vector
   *x*. The operation returns 0 for success and a non-zero value otherwise.

   Usage:

   .. code-block:: c

      retval = N_VDotProdMultiAllReduce(nv, x, d);


.. _NVectors.Ops.Exchange:

Exchange operations
-------------------

The following vector exchange operations are also *optional* and are
intended only for use when interfacing with the XBraid library for
parallel-in-time integration. In that setting these operations are
required but are otherwise unused by SUNDIALS packages and may be set
to ``NULL``. For each operation, we give the function signature, a
description of the expected behavior, and an example of the function
usage.



.. c:function:: int N_VBufSize(N_Vector x, sunindextype *size)

   This routine returns the buffer size need to exchange in the data in the
   vector *x* between computational nodes.

   Usage:

   .. code-block:: c

      flag = N_VBufSize(x, &buf_size)



.. c:function:: int N_VBufPack(N_Vector x, void *buf)

   This routine fills the exchange buffer *buf* with the vector data in *x*.

   Usage:

   .. code-block:: c

      flag = N_VBufPack(x, &buf)


.. c:function:: int N_VBufUnpack(N_Vector x, void *buf)

   This routine unpacks the data in the exchange buffer *buf* into the vector
   *x*.

   Usage:

   .. code-block:: c

      flag = N_VBufUnpack(x, buf)
