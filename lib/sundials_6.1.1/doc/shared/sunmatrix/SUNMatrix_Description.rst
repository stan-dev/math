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

.. _SUNMatrix.Description:

Description of the SUNMATRIX Modules
====================================

For problems that involve direct methods for solving linear systems,
the SUNDIALS packages not only operate on generic vectors, but also
on generic matrices (of type ``SUNMatrix``), through a set of
operations defined by the particular SUNMATRIX implementation.
Users can provide their own specific implementation of the
SUNMATRIX module, particularly in cases where they provide their
own ``N_Vector`` and/or linear solver modules, and require matrices
that are compatible with those implementations.  The generic
``SUNMatrix`` operations are described below, and descriptions of
the SUNMATRIX implementations provided with SUNDIALS follow.

The generic ``SUNMatrix`` type has been modeled after the
object-oriented style of the generic :c:type:`N_Vector` type.
Specifically, a generic ``SUNMatrix`` is a pointer to a structure
that has an implementation-dependent *content* field containing
the description and actual data of the matrix, and an *ops* field
pointing to a structure with generic matrix operations.
The type ``SUNMatrix`` is defined as:

.. c:type:: struct _generic_SUNMatrix *SUNMatrix

and the generic structure is defined as

.. code-block:: c

   struct _generic_SUNMatrix {
       void *content;
       struct _generic_SUNMatrix_Ops *ops;
   };

Here, the ``_generic_SUNMatrix_Ops`` structure is essentially a list of
function pointers to the various actual matrix operations, and is
defined as

.. code-block:: c

   struct _generic_SUNMatrix_Ops {
     SUNMatrix_ID (*getid)(SUNMatrix);
     SUNMatrix    (*clone)(SUNMatrix);
     void         (*destroy)(SUNMatrix);
     int          (*zero)(SUNMatrix);
     int          (*copy)(SUNMatrix, SUNMatrix);
     int          (*scaleadd)(realtype, SUNMatrix, SUNMatrix);
     int          (*scaleaddi)(realtype, SUNMatrix);
     int          (*matvecsetup)(SUNMatrix);
     int          (*matvec)(SUNMatrix, N_Vector, N_Vector);
     int          (*space)(SUNMatrix, long int*, long int*);
   };


The generic SUNMATRIX module defines and implements the matrix
operations acting on a ``SUNMatrix``. These routines are nothing but
wrappers for the matrix operations defined by a particular SUNMATRIX
implementation, which are accessed through the *ops* field of the
``SUNMatrix`` structure. To illustrate this point we show below the
implementation of a typical matrix operation from the generic
SUNMATRIX module, namely ``SUNMatZero``, which sets all values of a
matrix ``A`` to zero, returning a flag denoting a successful/failed
operation:

.. code-block:: c

   int SUNMatZero(SUNMatrix A)
   {
     return((int) A->ops->zero(A));
   }

:numref:`SUNMatrix.Ops` contains a complete list of all
matrix operations defined by the generic SUNMATRIX module.  A
particular implementation of the SUNMATRIX module must:

* Specify the *content* field of the ``SUNMatrix`` object.

* Define and implement a minimal subset of the matrix operations.
  See the documentation for each SUNDIALS package and/or linear solver
  to determine which SUNMATRIX operations they require.

  Note that the names of these routines should be unique to that
  implementation in order to permit using more than one SUNMATRIX
  module (each with different ``SUNMatrix`` internal data
  representations) in the same code.

* Define and implement user-callable constructor and destructor
  routines to create and free a ``SUNMatrix`` with the new *content*
  field and with *ops* pointing to the new matrix operations.

* Optionally, define and implement additional user-callable routines
  acting on the newly defined ``SUNMatrix`` (e.g., a routine to print the
  *content* for debugging purposes).

* Optionally, provide accessor macros as needed for that particular
  implementation to be used to access different parts in the content
  field of the newly defined ``SUNMatrix``.

To aid in the creation of custom SUNMATRIX modules the generic SUNMATRIX module
provides three utility functions :c:func:`SUNMatNewEmpty`,  :c:func:`SUNMatCopyOps()`,
and :c:func:`SUNMatFreeEmpty`. When used in custom SUNMATRIX constructors and clone
routines these functions will ease the introduction of any new optional matrix
operations to the SUNMATRIX API by ensuring only required operations need to be
set and all operations are copied when cloning a matrix.

.. c:function:: SUNMatrix SUNMatNewEmpty()

  This function allocates a new generic ``SUNMatrix`` object and initializes its
  content pointer and the function pointers in the operations structure to ``NULL``.

  **Return value:**
     If successful, this function returns a ``SUNMatrix`` object. If an error
     occurs when allocating the object, then this routine will return ``NULL``.

.. c:function:: int SUNMatCopyOps(SUNMatrix A, SUNMatrix B)

  This function copies the function pointers in the ``ops`` structure of ``A``
  into the ``ops`` structure of ``B``.

   **Arguments:**
      * *A* -- the matrix to copy operations from.
      * *B* -- the matrix to copy operations to.

   **Return value:**
      If successful, this function returns ``0``. If either of the inputs
      are ``NULL`` or the ``ops`` structure of either input is ``NULL``,
      then is function returns a non-zero value.

.. c:function:: void SUNMatFreeEmpty(SUNMatrix A)

  This routine frees the generic ``SUNMatrix`` object, under the assumption that any
  implementation-specific data that was allocated within the underlying content structure
  has already been freed. It will additionally test whether the ops pointer is ``NULL``,
  and, if it is not, it will free it as well.

   **Arguments:**
      * *A* -- the SUNMatrix object to free


Each SUNMATRIX implementation included in SUNDIALS has a unique
identifier specified in enumeration and shown in
:numref:`SUNMatrix.Description.matrixIDs`. It is recommended that a
user-supplied SUNMATRIX implementation use the ``SUNMATRIX_CUSTOM``
identifier.


.. _SUNMatrix.Description.matrixIDs:
.. table:: Identifiers associated with matrix kernels supplied with SUNDIALS
   :align: center

   ======================  =================================================  ========
   Matrix ID               Matrix type                                        ID Value
   ======================  =================================================  ========
   SUNMATRIX_DENSE         Dense :math:`M \times N` matrix                    0
   SUNMATRIX_MAGMADENSE    Magma dense :math:`M \times N` matrix              1
   SUNMATRIX_BAND          Band :math:`M \times M` matrix                     2
   SUNMATRIX_SPARSE        Sparse (CSR or CSC) :math:`M\times N` matrix       3
   SUNMATRIX_SLUNRLOC      SUNMatrix wrapper for SuperLU_DIST SuperMatrix     4
   SUNMATRIX_CUSPARSE      CUDA sparse CSR matrix                             5
   SUNMATRIX_CUSTOM        User-provided custom matrix                        6
   ======================  =================================================  ========
