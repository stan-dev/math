..
   Programmer(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------


.. _SUNMatrix.OneMklDense:

The SUNMATRIX_ONEMKLDENSE Module
================================

The SUNMATRIX_ONEMKLDENSE module is intended for interfacing with direct linear
solvers from the `Intel oneAPI Math Kernel Library (oneMKL) <https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html>`_
using the SYCL (DPC++) programming model. The implementation currently supports
a standard LAPACK column-major storage format as well as a low-storage format
for block-diagonal matrices,

.. math::

   \mathbf{A} =
   \begin{bmatrix}
      \mathbf{A_0} & 0            & \cdots & 0 \\
      0            & \mathbf{A_2} & \cdots & 0 \\
      \vdots       & \vdots       & \ddots & \vdots \\
      0            & 0            & \cdots & \mathbf{A_{n-1}}
   \end{bmatrix}

This matrix implementation is best paired with the
:ref:`SUNLinearSolver_OneMklDense <SUNLinSol.OneMklDense>` linear solver.

The header file to include when using this class is
``sunmatrix/sunmatrix_onemkldense.h``. The installed library to link to is
``libsundials_sunmatrixonemkldense.lib`` where ``lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

.. warning::

   The SUNMATRIX_ONEMKLDENSE class is experimental and subject to change.


SUNMATRIX_ONEMKLDENSE Functions
-------------------------------

The SUNMATRIX_ONEMKLDENSE class defines implementations of the following matrix
operations listed in :numref:`SUNMatrix.Ops`.

* ``SUNMatGetID_OneMklDense`` -- returns ``SUNMATRIX_ONEMKLDENSE``
* ``SUNMatClone_OneMklDense``
* ``SUNMatDestroy_OneMklDense``
* ``SUNMatZero_OneMklDense``
* ``SUNMatCopy_OneMklDense``
* ``SUNMatScaleAdd_OneMklDense``
* ``SUNMatScaleAddI_OneMklDense``
* ``SUNMatMatvec_OneMklDense``
* ``SUNMatSpace_OneMklDense``

In addition, the SUNMATRIX_ONEMKLDENSE class defines the following
implementation specific functions.


Constructors
^^^^^^^^^^^^


.. cpp:function:: SUNMatrix SUNMatrix_OneMklDense(sunindextype M, sunindextype N, SUNMemoryType memtype, SUNMemoryHelper memhelper, sycl::queue* queue, SUNContext sunctx)

   This constructor function creates and allocates memory for an
   :math:`M \times N` SUNMATRIX_ONEMKLDENSE ``SUNMatrix``.

   **Arguments:**
      * *M* -- the number of matrix rows.
      * *N* -- the number of matrix columns.
      * *memtype* -- the type of memory to use for the matrix data; can be
        ``SUNMEMTYPE_UVM`` or ``SUNMEMTYPE_DEVICE``.
      * *memhelper* -- the memory helper used for allocating data.
      * *queue* -- the SYCL queue to which operations will be submitted.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNMatrix`` object otherwise ``NULL``.


.. cpp:function:: SUNMatrix SUNMatrix_OneMklDenseBlock(sunindextype nblocks, sunindextype M_block, sunindextype N_block, SUNMemoryType memtype, SUNMemoryHelper memhelper, sycl::queue* queue, SUNContext sunctx)

   This constructor function creates and allocates memory for a block diagonal
   SUNMATRIX_ONEMKLDENSE ``SUNMatrix`` with *nblocks* of size
   :math:`M_{block} \times N_{block}`.

   **Arguments:**
      * *nblocks* -- the number of matrix rows.
      * *M_block* -- the number of matrix rows in each block.
      * *N_block* -- the number of matrix columns in each block.
      * *memtype* -- the type of memory to use for the matrix data; can be
        ``SUNMEMTYPE_UVM`` or ``SUNMEMTYPE_DEVICE``.
      * *memhelper* -- the memory helper used for allocating data.
      * *queue* -- the SYCL queue to which operations will be submitted.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNMatrix`` object otherwise ``NULL``.


Access Matrix Dimensions
^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: sunindextype SUNMatrix_OneMklDense_Rows(SUNMatrix A)

   This function returns the number of rows in the ``SUNMatrix`` object. For
   block diagonal matrices, the number of rows is computed as
   :math:`M_{\text{block}} \times \text{nblocks}`.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, the number of rows in the ``SUNMatrix`` object otherwise
      ``SUNMATRIX_ILL_INPUT``.


.. c:function:: sunindextype SUNMatrix_OneMklDense_Columns(SUNMatrix A)

   This function returns the number of columns in the ``SUNMatrix`` object. For
   block diagonal matrices, the number of columns is computed as
   :math:`N_{\text{block}} \times \text{nblocks}`.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, the number of columns in the ``SUNMatrix`` object otherwise
      ``SUNMATRIX_ILL_INPUT``.


Access Matrix Block Dimensions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: sunindextype SUNMatrix_OneMklDense_NumBlocks(SUNMatrix A)

   This function returns the number of blocks in the ``SUNMatrix`` object.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, the number of blocks in the ``SUNMatrix`` object otherwise
      ``SUNMATRIX_ILL_INPUT``.


.. c:function:: sunindextype SUNMatrix_OneMklDense_BlockRows(SUNMatrix A)

   This function returns the number of rows in a block of the ``SUNMatrix``
   object.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, the number of rows in a block of the ``SUNMatrix`` object
      otherwise ``SUNMATRIX_ILL_INPUT``.


.. c:function:: sunindextype SUNMatrix_OneMklDense_BlockColumns(SUNMatrix A)

   This function returns the number of columns in a block of the ``SUNMatrix``
   object.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, the number of columns in a block of the ``SUNMatrix``
      object otherwise ``SUNMATRIX_ILL_INPUT``.


Access Matrix Data
^^^^^^^^^^^^^^^^^^


.. c:function:: sunindextype SUNMatrix_OneMklDense_LData(SUNMatrix A)

   This function returns the length of the ``SUNMatrix`` data array.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, the length of the ``SUNMatrix`` data array otherwise
      ``SUNMATRIX_ILL_INPUT``.


.. c:function:: realtype* SUNMatrix_OneMklDense_Data(SUNMatrix A)

   This function returns the ``SUNMatrix`` data array.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, the ``SUNMatrix`` data array otherwise ``NULL``.


.. c:function:: realtype* SUNMatrix_OneMklDense_Column(SUNMatrix A, sunindextype j)

   This function returns a pointer to the data array for column *j* in the
   ``SUNMatrix``.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.
      * *j* -- the column index.

   **Return value:**
      If successful, a pointer to the data array for the ``SUNMatrix`` column
      otherwise ``NULL``.

   .. note::

      No bounds-checking is performed by this function, *j* should be strictly
      less than :math:`nblocks * N_{\text{block}}`.


Access Matrix Block Data
^^^^^^^^^^^^^^^^^^^^^^^^


.. c:function:: sunindextype SUNMatrix_OneMklDense_BlockLData(SUNMatrix A)

   This function returns the length of the ``SUNMatrix`` data array for each
   block of the ``SUNMatrix`` object.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, the length of the ``SUNMatrix`` data array for each block
      otherwise ``SUNMATRIX_ILL_INPUT``.


.. c:function:: realtype** SUNMatrix_OneMklDense_BlockData(SUNMatrix A)

   This function returns an array of pointers that point to the start of the
   data array for each block in the ``SUNMatrix``.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.

   **Return value:**
      If successful, an array of data pointers to each of the ``SUNMatrix``
      blocks otherwise ``NULL``.


.. c:function:: realtype* SUNMatrix_OneMklDense_Block(SUNMatrix A, sunindextype k)

   This function returns a pointer to the data array for block *k* in the
   ``SUNMatrix``.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.
      * *k* -- the block index.

   **Return value:**
      If successful, a pointer to the data array for the ``SUNMatrix`` block
      otherwise ``NULL``.

   .. note::

      No bounds-checking is performed by this function, *j* should be strictly
      less than *nblocks*.


.. c:function:: realtype* SUNMatrix_OneMklDense_BlockColumn(SUNMatrix A, sunindextype k, sunindextype j)

   This function returns a pointer to the data array for column *j* of block *k*
   in the ``SUNMatrix``.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object.
      * *k* -- the block index.
      * *j* -- the column index.

   **Return value:**
      If successful, a pointer to the data array for the ``SUNMatrix`` column
      otherwise ``NULL``.

   .. note::

      No bounds-checking is performed by this function, *k* should be strictly
      less than *nblocks* and *j* should be strictly less than
      :math:`N_{\text{block}}`.


Copy Data
^^^^^^^^^


.. c:function:: int SUNMatrix_OneMklDense_CopyToDevice(SUNMatrix A, realtype* h_data)

   This function copies the matrix data to the GPU device from the provided host
   array.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object
      * *h_data* -- a host array pointer to copy data from.

   **Return value:**
      * ``SUNMAT_SUCCESS`` -- if the copy is successful.
      * ``SUNMAT_ILL_INPUT`` -- if either the ``SUNMatrix`` is not a
        ``SUNMATRIX_ONEMKLDENSE`` matrix.
      * ``SUNMAT_MEM_FAIL`` -- if the copy fails.


.. c:function:: int SUNMatrix_OneMklDense_CopyFromDevice(SUNMatrix A, realtype* h_data)

   This function copies the matrix data from the GPU device to the provided host
   array.

   **Arguments:**
      * *A* -- a ``SUNMatrix`` object
      * *h_data* -- a host array pointer to copy data to.

   **Return value:**
      * ``SUNMAT_SUCCESS`` -- if the copy is successful.
      * ``SUNMAT_ILL_INPUT`` -- if either the ``SUNMatrix`` is not a
        ``SUNMATRIX_ONEMKLDENSE`` matrix.
      * ``SUNMAT_MEM_FAIL`` -- if the copy fails.


SUNMATRIX_ONEMKLDENSE Usage Notes
---------------------------------

.. warning::

   The SUNMATRIX_ONEMKLDENSE class only supports 64-bit indexing, thus SUNDIALS
   must be built for 64-bit indexing to use this class.

   When using the SUNMATRIX_ONEMKLDENSE class with a SUNDIALS package (e.g.
   CVODE), the queue given to matrix should be the same stream used for the
   NVECTOR object that is provided to the package, and the NVECTOR object given
   to the :c:func:`SUNMatMatvec` operation. If different streams are utilized,
   synchronization issues may occur.
