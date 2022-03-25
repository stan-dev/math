..
   Programmer(s): Cody J. Balos @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMatrix.cuSparse:

The SUNMATRIX_CUSPARSE Module
======================================

The SUNMATRIX_CUSPARSE module is an interface to the NVIDIA cuSPARSE matrix for
use on NVIDIA GPUs :cite:p:`cuSPARSE_site`.  All data stored by this matrix
implementation resides on the GPU at all times.

The header file to be included when using this module is ``sunmatrix/sunmatrix_cusparse.h``.
The installed library to link to is ``libsundials_sunmatrixcusparse.lib`` where ``.lib`` is
typically ``.so`` for shared libraries and ``.a`` for static libraries.

.. _SUNMatrix.cuSparse.Description:

SUNMATRIX_CUSPARSE Description
------------------------------

The implementation currently supports the cuSPARSE CSR matrix format described
in the cuSPARSE documentation, as well as a unique low-storage format for
block-diagonal matrices of the form

.. math::

   \mathbf{A} =
   \begin{bmatrix}
      \mathbf{A_0} & 0 & \cdots & 0\\
      0 & \mathbf{A_2} & \cdots & 0\\
      \vdots & \vdots & \ddots & \vdots\\
      0 & 0 & \cdots & \mathbf{A_{n-1}}\\
   \end{bmatrix},

where all the block matrices :math:`\mathbf{A_j}` share the same sparsity pattern.
We will refer to this format as BCSR (not to be confused with the canonical BSR format where
each block is stored as dense). In this format, the CSR column indices and row pointers
are only stored for the first block and are computed only as necessary for other blocks.
This can drastically reduce the amount of storage required compared to the regular CSR
format when the number of blocks is large. This format is well-suited for, and
intended to be used with, the ``SUNLinearSolver_cuSolverSp_batchQR`` linear solver
(see :numref:`SUNLinSol.cuSolverSp`).

**The SUNMATRIX_CUSPARSE module is experimental and subject to change.**

.. _SUNMatrix.cuSparse.Functions:

SUNMATRIX_CUSPARSE Functions
----------------------------------

The SUNMATRIX_CUSPARSE module defines GPU-enabled sparse implementations of all matrix
operations listed in :numref:`SUNMatrix.Ops` except for the :c:func:`SUNMatSpace`
and :c:func:`SUNMatMatvecSetup` operations:

* ``SUNMatGetID_cuSparse`` -- returns ``SUNMATRIX_CUSPARSE``

* ``SUNMatClone_cuSparse``

* ``SUNMatDestroy_cuSparse``

* ``SUNMatZero_cuSparse``

* ``SUNMatCopy_cuSparse``

* ``SUNMatScaleAdd_cuSparse`` -- performs :math:`A = cA + B`, where :math:`A` and :math:`B`
  must have the same sparsity pattern

* ``SUNMatScaleAddI_cuSparse`` -- performs :math:`A = cA + I`, where the diagonal of :math:`A`
  must be present

* ``SUNMatMatvec_cuSparse``


In addition, the SUNMATRIX_CUSPARSE module defines the following implementation specific
functions:

.. c:function:: SUNMatrix SUNMatrix_cuSparse_NewCSR(int M, int N, int NNZ, cusparseHandle_t cusp, SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNMATRIX_CUSPARSE
   ``SUNMatrix`` that uses the CSR storage format. Its arguments are the
   number of rows and columns of the matrix, ``M`` and ``N``, the number of
   nonzeros to be stored in the matrix, ``NNZ``, and a valid ``cusparseHandle_t``.


.. c:function:: SUNMatrix SUNMatrix_cuSparse_NewBlockCSR(int nblocks, int blockrows, int blockcols, int blocknnz, cusparseHandle_t cusp, SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNMATRIX_CUSPARSE
   ``SUNMatrix`` object that leverages the ``SUNMAT_CUSPARSE_BCSR`` storage
   format to store a block diagonal matrix where each block shares the same
   sparsity pattern. The blocks must be square. The function arguments are the
   number of blocks, ``nblocks``, the number of rows, ``blockrows``, the number of
   columns, ``blockcols``, the number of nonzeros in each each block, ``blocknnz``,
   and a valid ``cusparseHandle_t``.

   .. warning::

      The ``SUNMAT_CUSPARSE_BCSR`` format currently only supports square matrices, i.e.,
      ``blockrows == blockcols``.


.. c:function:: SUNMatrix SUNMatrix_cuSparse_MakeCSR(cusparseMatDescr_t mat_descr, int M, int N, int NNZ, int *rowptrs , int *colind , realtype *data, cusparseHandle_t cusp, SUNContext sunctx)

   This constructor function creates a SUNMATRIX_CUSPARSE ``SUNMatrix``
   object from user provided pointers. Its arguments are a ``cusparseMatDescr_t``
   that must have index base ``CUSPARSE_INDEX_BASE_ZERO``, the number of rows
   and columns of the matrix, ``M`` and ``N``, the number of nonzeros to be stored
   in the matrix, ``NNZ``, and a valid ``cusparseHandle_t``.


.. c:function:: int SUNMatrix_cuSparse_Rows(SUNMatrix A)

   This function returns the number of rows in the sparse ``SUNMatrix``.


.. c:function:: int SUNMatrix_cuSparse_Columns(SUNMatrix A)

   This function returns the number of columns in the sparse ``SUNMatrix``.


.. c:function:: int SUNMatrix_cuSparse_NNZ(SUNMatrix A)

   This function returns the number of entries allocated for nonzero
   storage for the sparse ``SUNMatrix``.


.. c:function:: int SUNMatrix_cuSparse_SparseType(SUNMatrix A)

   This function returns the storage type (``SUNMAT_CUSPARSE_CSR``
   or ``SUNMAT_CUSPARSE_BCSR``) for the sparse ``SUNMatrix``.


.. c:function:: realtype* SUNMatrix_cuSparse_Data(SUNMatrix A)

   This function returns a pointer to the data array for the
   sparse ``SUNMatrix``.


.. c:function:: int* SUNMatrix_cuSparse_IndexValues(SUNMatrix A)

   This function returns a pointer to the index value array for the sparse
   ``SUNMatrix`` -- for the CSR format this is an array of column indices for
   each nonzero entry. For the BCSR format this is an array of the column indices
   for each nonzero entry in the first block only.


.. c:function:: int* SUNMatrix_cuSparse_IndexPointers(SUNMatrix A)

   This function returns a pointer to the index pointer array for the
   sparse ``SUNMatrix`` -- for the CSR format this is an array of the locations
   of the first entry of each row in the ``data`` and ``indexvalues`` arrays,
   for the BCSR format this is an array of the locations of each row in the
   ``data`` and ``indexvalues`` arrays in the first block only.


.. c:function:: int SUNMatrix_cuSparse_NumBlocks(SUNMatrix A)

   This function returns the number of matrix blocks.


.. c:function:: int SUNMatrix_cuSparse_BlockRows(SUNMatrix A)

   This function returns the number of rows in a matrix block.


.. c:function:: int SUNMatrix_cuSparse_BlockColumns(SUNMatrix A)

   This function returns the number of columns in a matrix block.


.. c:function:: int SUNMatrix_cuSparse_BlockNNZ(SUNMatrix A)

   This function returns the number of nonzeros in each
   matrix block.


.. c:function:: realtype* SUNMatrix_cuSparse_BlockData(SUNMatrix A, int blockidx)

   This function returns a pointer to the location in the ``data`` array
   where the data for the block, ``blockidx``, begins. Thus, ``blockidx``
   must be less than ``SUNMatrix_cuSparse_NumBlocks(A)``. The first block
   in the SUNMatrix is index 0, the second block is index 1, and so on.


.. c:function:: cusparseMatDescr_t SUNMatrix_cuSparse_MatDescr(SUNMatrix A)

   This function returns the ``cusparseMatDescr_t`` object associated with
   the matrix.


.. c:function:: int SUNMatrix_cuSparse_CopyToDevice(SUNMatrix A, realtype* h_data, int* h_idxptrs, int* h_idxvals)

   This functions copies the matrix information to the GPU device from the provided
   host arrays. A user may provide ``NULL`` for any of ``h_data``, ``h_idxptrs``, or
   ``h_idxvals`` to avoid copying that information.

   The function returns ``SUNMAT_SUCCESS`` if the copy operation(s) were successful,
   or a nonzero error code otherwise.

.. c:function:: int SUNMatrix_cuSparse_CopyFromDevice(SUNMatrix A, realtype* h_data, int* h_idxptrs, int* h_idxvals)

   This functions copies the matrix information from the GPU device to the provided
   host arrays. A user may provide ``NULL`` for any of ``h_data``, ``h_idxptrs``, or
   ``h_idxvals`` to avoid copying that information. Otherwise:

   * The ``h_data`` array must be at least ``SUNMatrix_cuSparse_NNZ(A)*sizeof(realtype)``
     bytes.

   * The ``h_idxptrs`` array must be at least
     ``(SUNMatrix_cuSparse_BlockDim(A)+1)*sizeof(int)`` bytes.

   * The ``h_idxvals`` array must be at least
     ``(SUNMatrix_cuSparse_BlockNNZ(A))*sizeof(int)`` bytes.

   The function returns ``SUNMAT_SUCCESS`` if the copy operation(s) were successful,
   or a nonzero error code otherwise.


.. c:function:: int SUNMatrix_cuSparse_SetFixedPattern(SUNMatrix A, booleantype yesno)

   This function changes the behavior of the the ``SUNMatZero`` operation on the object
   ``A``.  By default the matrix sparsity pattern is not considered to be fixed, thus,
   the ``SUNMatZero`` operation zeros out all ``data`` array as well as the ``indexvalues``
   and ``indexpointers`` arrays. Providing a value of ``1`` or ``SUNTRUE`` for the
   ``yesno`` argument changes the behavior of ``SUNMatZero`` on ``A`` so that only the
   data is zeroed out, but not the ``indexvalues`` or ``indexpointers`` arrays.
   Providing a value of ``0`` or ``SUNFALSE`` for the ``yesno`` argument is equivalent
   to the default behavior.

.. c:function:: int SUNMatrix_cuSparse_SetKernelExecPolicy(SUNMatrix A, SUNCudaExecPolicy* exec_policy)

   This function sets the execution policies which control the kernel parameters
   utilized when launching the CUDA kernels. By default the matrix is setup to use
   a policy which tries to leverage the structure of the matrix. See
   :numref:`NVectors.CUDA.SUNCudaExecPolicy` for more information about the
   :cpp:type:`SUNCudaExecPolicy` class.


.. _SUNMatrix.cuSparse.Notes:

SUNMATRIX_CUSPARSE Usage Notes
----------------------------------

The SUNMATRIX_CUSPARSE module only supports 32-bit indexing, thus SUNDIALS must be built
for 32-bit indexing to use this module.

The SUNMATRIX_CUSPARSE module can be used with CUDA streams by calling the cuSPARSE
function ``cusparseSetStream`` on the ``cusparseHandle_t`` that is provided to the
SUNMATRIX_CUSPARSE constructor.

.. warning::

   When using the SUNMATRIX_CUSPARSE module with a SUNDIALS package (e.g. ARKODE), the
   stream given to cuSPARSE should be the same stream used for the NVECTOR object that
   is provided to the package, and the NVECTOR object given to the ``SUNMatvec`` operation.
   If different streams are utilized, synchronization issues may occur.
