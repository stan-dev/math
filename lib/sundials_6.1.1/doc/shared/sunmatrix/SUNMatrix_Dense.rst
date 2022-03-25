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

.. _SUNMatrix.Dense:

The SUNMATRIX_DENSE Module
======================================

The dense implementation of the ``SUNMatrix`` module, SUNMATRIX_DENSE,
defines the *content* field of ``SUNMatrix`` to be the following structure:

.. code-block:: c

   struct _SUNMatrixContent_Dense {
     sunindextype M;
     sunindextype N;
     realtype *data;
     sunindextype ldata;
     realtype **cols;
   };

These entries of the *content* field contain the following information:

* ``M`` - number of rows

* ``N`` - number of columns

* ``data`` - pointer to a contiguous block of ``realtype`` variables.
  The elements of the dense matrix are stored columnwise, i.e. the
  :math:`(i,j)` element of a dense ``SUNMatrix`` object
  (with :math:`0 \le i < M` and :math:`0 \le j < N`) may be accessed
  via ``data[j*M+i]``.

* ``ldata`` - length of the data array (:math:`= M\, N`).

* ``cols`` - array of pointers. ``cols[j]`` points to the first
  element of the j-th column of the matrix in the array ``data``.
  The :math:`(i,j)` element of a dense ``SUNMatrix``
  (with :math:`0 \le i < M` and :math:`0 \le j < N`) may be accessed
  may be accessed via ``cols[j][i]``.


The header file to be included when using this module is
``sunmatrix/sunmatrix_dense.h``.

The following macros are provided to access the content of a
SUNMATRIX_DENSE matrix. The prefix ``SM_`` in the names denotes that
these macros are for *SUNMatrix* implementations, and the suffix
``_D`` denotes that these are specific to the *dense* version.


.. c:macro:: SM_CONTENT_D(A)

   This macro gives access to the contents of the dense ``SUNMatrix`` *A*.

   The assignment ``A_cont = SM_CONTENT_D(A)`` sets
   ``A_cont`` to be a pointer to the dense ``SUNMatrix`` content
   structure.

   Implementation:

   .. code-block:: c

      #define SM_CONTENT_D(A)   ( (SUNMatrixContent_Dense)(A->content) )


.. c:macro:: SM_ROWS_D(A)

   Access the number of rows in the dense ``SUNMatrix`` *A*.

   This may be used either to retrieve or to set the value.  For
   example, the assignment ``A_rows = SM_ROWS_D(A)`` sets ``A_rows`` to be
   the number of rows in the matrix ``A``.  Similarly, the
   assignment ``SM_ROWS_D(A) = A_rows`` sets the number of
   columns in ``A`` to equal ``A_rows``.

   Implementation:

   .. code-block:: c

      #define SM_ROWS_D(A)   ( SM_CONTENT_D(A)->M )


.. c:macro:: SM_COLUMNS_D(A)

   Access the number of columns in the dense ``SUNMatrix`` *A*.

   This may be used either to retrieve or to set the value.  For
   example, the assignment ``A_columns = SM_COLUMNS_D(A)`` sets
   ``A_columns`` to be the number of columns in the matrix ``A``.
   Similarly, the assignment ``SM_COLUMNS_D(A) = A_columns`` sets the
   number of columns in ``A`` to equal ``A_columns``

   Implementation:

   .. code-block:: c

      #define SM_COLUMNS_D(A)   ( SM_CONTENT_D(A)->N )


.. c:macro:: SM_LDATA_D(A)

   Access the total data length in the dense ``SUNMatrix`` *A*.

   This may be used either to retrieve or to set the value.  For
   example, the assignment ``A_ldata = SM_LDATA_D(A)`` sets
   ``A_ldata`` to be the length of the data array in the matrix ``A``.
   Similarly, the assignment ``SM_LDATA_D(A) = A_ldata`` sets the
   parameter for the length of the data array in ``A`` to equal
   ``A_ldata``.

   Implementation:

   .. code-block:: c

      #define SM_LDATA_D(A)   ( SM_CONTENT_D(A)->ldata )


.. c:macro:: SM_DATA_D(A)

   This macro gives access to the ``data`` pointer for the matrix entries.

   The assignment ``A_data = SM_DATA_D(A)`` sets ``A_data`` to be
   a pointer to the first component of the data array for the dense
   ``SUNMatrix A``.  The assignment ``SM_DATA_D(A) = A_data``
   sets the data array of ``A`` to be ``A_data`` by storing the
   pointer ``A_data``.

   Implementation:

   .. code-block:: c

      #define SM_DATA_D(A)   ( SM_CONTENT_D(A)->data )


.. c:macro:: SM_COLS_D(A)

   This macro gives access to the ``cols`` pointer for the matrix entries.

   The assignment ``A_cols = SM_COLS_D(A)`` sets ``A_cols`` to be
   a pointer to the array of column pointers for the dense ``SUNMatrix A``.
   The assignment ``SM_COLS_D(A) = A_cols`` sets the column pointer
   array of ``A`` to be ``A_cols`` by storing the pointer
   ``A_cols``.

   Implementation:

   .. code-block:: c

      #define SM_COLS_D(A)   ( SM_CONTENT_D(A)->cols )


.. c:macro:: SM_COLUMN_D(A)

   This macros gives access to the individual columns of the data
   array of a dense ``SUNMatrix``.

   The assignment ``col_j = SM_COLUMN_D(A,j)`` sets ``col_j`` to be
   a pointer to the first entry of the ``j``-th column of the :math:`M \times N`
   dense matrix ``A`` (with :math:`0 \le j < N`).  The type of the
   expression ``SM_COLUMN_D(A,j)`` is ``realtype *``.  The pointer
   returned by the call ``SM_COLUMN_D(A,j)`` can be treated as
   an array which is indexed from 0 to ``M-1``.

   Implementation:

   .. code-block:: c

      #define SM_COLUMN_D(A,j)    ( (SM_CONTENT_D(A)->cols)[j] )


.. c:macro:: SM_ELEMENT_D(A)

   This macro gives access to the individual entries of the data array
   of a dense ``SUNMatrix``.

   The assignments ``SM_ELEMENT_D(A,i,j) = a_ij`` and ``a_ij =
   SM_ELEMENT_D(A,i,j)`` reference the :math:`A_{i,j}` element of the
   :math:`M \times N` dense matrix ``A`` (with :math:`0 \le i < M` and
   :math:`0 \le j < N`).

   Implementation:

   .. code-block:: c

      #define SM_ELEMENT_D(A,i,j) ( (SM_CONTENT_D(A)->cols)[j][i] )



The SUNMATRIX_DENSE module defines dense implementations of all matrix
operations listed in :numref:`SUNMatrix.Ops`. Their names are obtained
from those in that section by appending the suffix ``_Dense``
(e.g. ``SUNMatCopy_Dense``).  The module SUNMATRIX_DENSE provides the
following additional user-callable routines:


.. c:function:: SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N, SUNContext sunctx)

   This constructor function creates and allocates memory for a dense
   ``SUNMatrix``.  Its arguments are the number of rows, ``M``, and
   columns, ``N``, for the dense matrix.


.. c:function:: void SUNDenseMatrix_Print(SUNMatrix A, FILE* outfile)

   This function prints the content of a dense ``SUNMatrix`` to the
   output stream specified by ``outfile``.  Note: ``stdout``
   or ``stderr`` may be used as arguments for ``outfile`` to print
   directly to standard output or standard error, respectively.


.. c:function:: sunindextype SUNDenseMatrix_Rows(SUNMatrix A)

   This function returns the number of rows in the dense ``SUNMatrix``.


.. c:function:: sunindextype SUNDenseMatrix_Columns(SUNMatrix A)

   This function returns the number of columns in the dense ``SUNMatrix``.


.. c:function:: sunindextype SUNDenseMatrix_LData(SUNMatrix A)

   This function returns the length of the data array for the dense ``SUNMatrix``.


.. c:function:: realtype* SUNDenseMatrix_Data(SUNMatrix A)

   This function returns a pointer to the data array for the dense ``SUNMatrix``.


.. c:function:: realtype** SUNDenseMatrix_Cols(SUNMatrix A)

   This function returns a pointer to the cols array for the dense ``SUNMatrix``.


.. c:function:: realtype* SUNDenseMatrix_Column(SUNMatrix A, sunindextype j)

   This function returns a pointer to the first entry of the jth
   column of the dense ``SUNMatrix``.  The resulting pointer should
   be indexed over the range ``0`` to ``M-1``.



**Notes**

* When looping over the components of a dense ``SUNMatrix A``,
  the most efficient approaches are to:

  * First obtain the component array via ``A_data = SUNDenseMatrix_Data(A)``,
    or equivalently ``A_data = SM_DATA_D(A)``, and then access ``A_data[i]``
    within the loop.

  * First obtain the array of column pointers via
    ``A_cols = SUNDenseMatrix_Cols(A)``, or equivalently
    ``A_cols = SM_COLS_D(A)``, and then access ``A_cols[j][i]`` within the loop.

  * Within a loop over the columns, access the column pointer via
    ``A_colj = SUNDenseMatrix_Column(A,j)`` and then to access the
    entries within that column using ``A_colj[i]`` within the loop.

  All three of these are more efficient than
  using ``SM_ELEMENT_D(A,i,j)`` within a double loop.

* Within the ``SUNMatMatvec_Dense`` routine, internal consistency
  checks are performed to ensure that the matrix is called with
  consistent ``N_Vector`` implementations.  These are currently
  limited to: NVECTOR_SERIAL, NVECTOR_OPENMP, and NVECTOR_PTHREADS.
  As additional compatible vector implementations are added to
  SUNDIALS, these will be included within this compatibility check.
