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

.. _SUNMatrix.Sparse:

The SUNMATRIX_SPARSE Module
======================================

The sparse implementation of the ``SUNMatrix`` module, SUNMATRIX_SPARSE,
is designed to work with either *compressed-sparse-column* (CSC) or
*compressed-sparse-row* (CSR) sparse matrix formats.  To this end, it
defines the *content* field of ``SUNMatrix`` to be the following
structure:

.. code-block:: c

   struct _SUNMatrixContent_Sparse {
     sunindextype M;
     sunindextype N;
     sunindextype NNZ;
     sunindextype NP;
     realtype *data;
     int sparsetype;
     sunindextype *indexvals;
     sunindextype *indexptrs;
     /* CSC indices */
     sunindextype **rowvals;
     sunindextype **colptrs;
     /* CSR indices */
     sunindextype **colvals;
     sunindextype **rowptrs;
   };

A diagram of the underlying data representation in a sparse matrix is
shown in :numref:`SUNSparseMatrix`.  A more
complete description of the parts of this *content* field is given below:

* ``M`` - number of rows

* ``N`` - number of columns

* ``NNZ`` - maximum number of nonzero entries in the matrix
  (allocated length of ``data`` and ``indexvals`` arrays)

* ``NP`` - number of index pointers (e.g. number of column pointers
  for CSC matrix). For CSC matrices ``NP=N``, and for CSR matrices
  ``NP=M``. This value is set automatically at construction based the
  input choice for ``sparsetype``.

* ``data`` - pointer to a contiguous block of ``realtype``
  variables (of length ``NNZ``), containing the values of the
  nonzero entries in the matrix

* ``sparsetype`` - type of the sparse matrix (``CSC_MAT`` or ``CSR_MAT``)

* ``indexvals`` - pointer to a contiguous block of ``int`` variables
  (of length ``NNZ``), containing the row indices (if CSC) or column
  indices (if CSR) of each nonzero matrix entry held in ``data``

* ``indexptrs`` - pointer to a contiguous block of ``int``
  variables (of length ``NP+1``). For CSC matrices each entry provides
  the index of the first column entry into the ``data`` and
  ``indexvals`` arrays, e.g. if ``indexptr[3]=7``, then the first
  nonzero entry in the fourth column of the matrix is located in
  ``data[7]``, and is located in row ``indexvals[7]`` of the matrix.
  The last entry contains the total number of nonzero values in the
  matrix and hence points one past the end of the active data in the
  ``data`` and ``indexvals`` arrays. For CSR matrices, each entry
  provides the index of the first row entry into the ``data`` and
  ``indexvals`` arrays.

The following pointers are added to the SUNMATRIX_SPARSE content
structure for user convenience, to provide a more intuitive interface
to the CSC and CSR sparse matrix data structures. They are set
automatically when creating a sparse ``SUNMatrix``, based on the
sparse matrix storage type.

* ``rowvals`` - pointer to ``indexvals`` when ``sparsetype`` is
  ``CSC_MAT``, otherwise set to ``NULL``.

* ``colptrs`` - pointer to ``indexptrs`` when ``sparsetype`` is
  ``CSC_MAT``, otherwise set to ``NULL``.

* ``colvals`` - pointer to ``indexvals`` when ``sparsetype`` is
  ``CSR_MAT``, otherwise set to ``NULL``.

* ``rowptrs`` - pointer to ``indexptrs`` when ``sparsetype`` is
  ``CSR_MAT``, otherwise set to ``NULL``.

For example, the :math:`5\times 4` matrix

.. math::
   \left[\begin{array}{cccc}
     0 & 3 & 1 & 0\\
     3 & 0 & 0 & 2\\
     0 & 7 & 0 & 0\\
     1 & 0 & 0 & 9\\
     0 & 0 & 0 & 5
   \end{array}\right]

could be stored as a CSC matrix in this structure as either

.. code-block:: c

   M = 5;
   N = 4;
   NNZ = 8;
   NP = N;
   data = {3.0, 1.0, 3.0, 7.0, 1.0, 2.0, 9.0, 5.0};
   sparsetype = CSC_MAT;
   indexvals = {1, 3, 0, 2, 0, 1, 3, 4};
   indexptrs = {0, 2, 4, 5, 8};

or

.. code-block:: c

   M = 5;
   N = 4;
   NNZ = 10;
   NP = N;
   data = {3.0, 1.0, 3.0, 7.0, 1.0, 2.0, 9.0, 5.0, *, *};
   sparsetype = CSC_MAT;
   indexvals = {1, 3, 0, 2, 0, 1, 3, 4, *, *};
   indexptrs = {0, 2, 4, 5, 8};

where the first has no unused space, and the second has additional
storage (the entries marked with ``*`` may contain any values).
Note in both cases that the final value in ``indexptrs`` is 8,
indicating the total number of nonzero entries in the matrix.

Similarly, in CSR format, the same matrix could be stored as

.. code-block:: c

   M = 5;
   N = 4;
   NNZ = 8;
   NP = M;
   data = {3.0, 1.0, 3.0, 2.0, 7.0, 1.0, 9.0, 5.0};
   sparsetype = CSR_MAT;
   indexvals = {1, 2, 0, 3, 1, 0, 3, 3};
   indexptrs = {0, 2, 4, 5, 7, 8};


.. _SUNSparseMatrix:
.. figure:: /figs/cscmat.png

   Diagram of the storage for a compressed-sparse-column matrix of
   type SUNMATRIX_SPARSE: Here ``A`` is an :math:`M \times N` sparse
   CSC matrix with storage for up to ``NNZ`` nonzero entries (the
   allocated length of both ``data`` and ``indexvals``).  The entries
   in ``indexvals`` may assume values from ``0`` to ``M-1``,
   corresponding to the row index (zero-based) of
   each nonzero value.  The entries in ``data`` contain the values of
   the nonzero entries, with the row ``i``, column ``j`` entry of
   ``A`` (again, zero-based) denoted as ``A(i,j)``.  The ``indexptrs``
   array contains ``N+1`` entries; the first ``N`` denote the starting
   index of each column within the ``indexvals`` and ``data`` arrays,
   while the final entry points one past the final nonzero entry.
   Here, although ``NNZ`` values are allocated, only ``nz`` are
   actually filled in; the greyed-out portions of ``data`` and
   ``indexvals`` indicate extra allocated space.


The header file to be included when using this module is
``sunmatrix/sunmatrix_sparse.h``.

The following macros are provided to access the content of a
SUNMATRIX_SPARSE matrix. The prefix ``SM_`` in the names
denotes that these macros are for *SUNMatrix* implementations,
and the suffix ``_S`` denotes that these are specific to
the *sparse* version.


.. c:macro:: SM_CONTENT_S(A)

   This macro gives access to the contents of the sparse ``SUNMatrix`` *A*.

   The assignment ``A_cont = SM_CONTENT_S(A)`` sets
   ``A_cont`` to be a pointer to the sparse ``SUNMatrix`` content
   structure.

   Implementation:

   .. code-block:: c

      #define SM_CONTENT_S(A)   ( (SUNMatrixContent_Sparse)(A->content) )


.. c:macro:: SM_ROWS_S(A)

   Access the number of rows in the sparse ``SUNMatrix`` *A*.

   This may be used either to retrieve or to set the value.  For
   example, the assignment ``A_rows = SM_ROWS_S(A)`` sets ``A_rows``
   to be the number of rows in the matrix *A*.  Similarly, the
   assignment ``SM_ROWS_S(A) = A_rows`` sets the number of
   columns in *A* to equal ``A_rows``.

   Implementation:

   .. code-block:: c

      #define SM_ROWS_S(A)   ( SM_CONTENT_S(A)->M )


.. c:macro:: SM_COLUMNS_S(A)

   Access the number of columns in the sparse ``SUNMatrix`` *A*.  As
   with ``SM_ROWS_S``, this may be used either to retrieve or to set
   the value.

   Implementation:

   .. code-block:: c

      #define SM_COLUMNS_S(A)   ( SM_CONTENT_S(A)->N )


.. c:macro:: SM_NNZ_S(A)

   Access the allocated number of nonzeros in the sparse ``SUNMatrix``
   *A*.  As with ``SM_ROWS_S``, this may be used either to retrieve or
   to set the value.

   Implementation:

   .. code-block:: c

      #define SM_NNZ_S(A)   ( SM_CONTENT_S(A)->NNZ )


.. c:macro:: SM_NP_S(A)

   Access the number of index pointers ``NP`` in the sparse
   ``SUNMatrix`` *A*.  As with ``SM_ROWS_S``, this may be used either
   to retrieve or to set the value.

   Implementation:

   .. code-block:: c

      #define SM_NP_S(A)   ( SM_CONTENT_S(A)->NP )


.. c:macro:: SM_SPARSETYPE_S(A)

   Access the sparsity type parameter in the sparse ``SUNMatrix`` *A*.
   As with ``SM_ROWS_S``, this may be used either to retrieve or to
   set the value.

   Implementation:

   .. code-block:: c

      #define SM_SPARSETYPE_S(A)   ( SM_CONTENT_S(A)->sparsetype )


.. c:macro:: SM_DATA_S(A)

   This macro gives access to the ``data`` pointer for the matrix
   entries.

   The assignment ``A_data = SM_DATA_S(A)`` sets ``A_data`` to be
   a pointer to the first component of the data array for the sparse
   ``SUNMatrix A``.  The assignment ``SM_DATA_S(A) = A_data``
   sets the data array of *A* to be ``A_data`` by storing the
   pointer ``A_data``.

   Implementation:

   .. code-block:: c

      #define SM_DATA_S(A)   ( SM_CONTENT_S(A)->data )


.. c:macro:: SM_INDEXVALS_S(A)

   This macro gives access to the ``indexvals`` pointer for the matrix
   entries.

   The assignment ``A_indexvals = SM_INDEXVALS_S(A)``
   sets ``A_indexvals`` to be a pointer to the array of index values
   (i.e. row indices for a CSC matrix, or column indices for a CSR
   matrix) for the sparse ``SUNMatrix`` *A*.

   Implementation:

   .. code-block:: c

      #define SM_INDEXVALS_S(A)   ( SM_CONTENT_S(A)->indexvals )


.. c:macro:: SM_INDEXPTRS_S(A)

   This macro gives access to the ``indexptrs`` pointer for the matrix entries.

   The assignment ``A_indexptrs = SM_INDEXPTRS_S(A)``
   sets ``A_indexptrs`` to be a pointer to the array of index
   pointers (i.e. the starting indices in the data/indexvals arrays for
   each row or column in CSR or CSC formats, respectively).

   Implementation:

   .. code-block:: c

      #define SM_INDEXPTRS_S(A)   ( SM_CONTENT_S(A)->indexptrs )


The SUNMATRIX_SPARSE module defines sparse implementations of all matrix
operations listed in :numref:`SUNMatrix.Ops`. Their names are
obtained from those in that section by appending the suffix ``_Sparse``
(e.g. ``SUNMatCopy_Sparse``).  The module SUNMATRIX_SPARSE provides the
following additional user-callable routines:


.. c:function:: SUNMatrix SUNSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ, int sparsetype, SUNContext sunctx)

   This constructor function creates and allocates memory for a sparse
   ``SUNMatrix``.  Its arguments are the number of rows and columns of
   the matrix, *M* and *N*, the maximum number of nonzeros to be
   stored in the matrix, *NNZ*, and a flag *sparsetype* indicating
   whether to use CSR or CSC format (valid choices are ``CSR_MAT`` or
   ``CSC_MAT``).



.. c:function:: SUNMatrix SUNSparseFromDenseMatrix(SUNMatrix A, realtype droptol, int sparsetype)

   This constructor function creates a new sparse matrix from an
   existing SUNMATRIX_DENSE object by copying all values with
   magnitude larger than *droptol* into the sparse matrix structure.

   Requirements:

   * *A* must have type ``SUNMATRIX_DENSE``

   * *droptol* must be non-negative

   * *sparsetype* must be either ``CSC_MAT`` or ``CSR_MAT``

   The function returns ``NULL`` if any requirements are violated, or if
   the matrix storage request cannot be satisfied.



.. c:function:: SUNMatrix SUNSparseFromBandMatrix(SUNMatrix A, realtype droptol, int sparsetype)

   This constructor function creates a new sparse matrix from an
   existing SUNMATRIX_BAND object by copying all values with
   magnitude larger than *droptol* into the sparse matrix structure.

   Requirements:

   * *A* must have type ``SUNMATRIX_BAND``

   * *droptol* must be non-negative

   * *sparsetype* must be either ``CSC_MAT`` or ``CSR_MAT``.

   The function returns ``NULL`` if any requirements are violated, or if
   the matrix storage request cannot be satisfied.



.. c:function:: int SUNSparseMatrix_Realloc(SUNMatrix A)

   This function reallocates internal storage arrays in a sparse matrix
   so that the resulting sparse matrix has no wasted space (i.e. the
   space allocated for nonzero entries equals the actual number of
   nonzeros, ``indexptrs[NP]``). Returns 0 on success and
   1 on failure (e.g. if the input matrix is not sparse).



.. c:function:: void SUNSparseMatrix_Print(SUNMatrix A, FILE* outfile)

   This function prints the content of a sparse ``SUNMatrix`` to the
   output stream specified by ``outfile``.  Note: ``stdout``
   or ``stderr`` may be used as arguments for ``outfile`` to print
   directly to standard output or standard error, respectively.


.. c:function:: sunindextype SUNSparseMatrix_Rows(SUNMatrix A)

   This function returns the number of rows in the sparse ``SUNMatrix``.


.. c:function:: sunindextype SUNSparseMatrix_Columns(SUNMatrix A)

   This function returns the number of columns in the sparse ``SUNMatrix``.


.. c:function:: sunindextype SUNSparseMatrix_NNZ(SUNMatrix A)

   This function returns the number of entries allocated for nonzero
   storage for the sparse ``SUNMatrix``.


.. c:function:: sunindextype SUNSparseMatrix_NP(SUNMatrix A)

   This function returns the number of index pointers for the
   sparse ``SUNMatrix`` (the ``indexptrs`` array has ``NP+1``
   entries).


.. c:function:: int SUNSparseMatrix_SparseType(SUNMatrix A)

   This function returns the storage type (``CSR_MAT``
   or ``CSC_MAT``) for the sparse  ``SUNMatrix``.


.. c:function:: realtype* SUNSparseMatrix_Data(SUNMatrix A)

   This function returns a pointer to the data array for the
   sparse ``SUNMatrix``.


.. c:function:: sunindextype* SUNSparseMatrix_IndexValues(SUNMatrix A)

   This function returns a pointer to index value array for the sparse
   ``SUNMatrix`` -- for CSR format this is the column index for each nonzero
   entry, for CSC format this is the row index for each nonzero entry.


.. c:function:: sunindextype* SUNSparseMatrix_IndexPointers(SUNMatrix A)

   This function returns a pointer to the index pointer array for the
   sparse ``SUNMatrix`` -- for CSR format this is the location of the first
   entry of each row in the ``data`` and ``indexvalues`` arrays, for
   CSC format this is the location of the first entry of each column.


.. note:: Within the ``SUNMatMatvec_Sparse`` routine, internal
          consistency checks are performed to ensure that the matrix
          is called with consistent ``N_Vector`` implementations.
          These are currently limited to: NVECTOR_SERIAL,
          NVECTOR_OPENMP, NVECTOR_PTHREADS, and NVECTOR_CUDA when using
          managed memory. As additional compatible vector implementations
          are added to SUNDIALS, these will be included within this
          compatibility check.
