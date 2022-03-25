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

.. _NVectors.NVTrilinos:

The NVECTOR_TRILINOS Module
===========================

The NVECTOR_TRILINOS module is an NVECTOR wrapper around the
`Trilinos <https://github.com/trilinos/Trilinos>`_ Tpetra vector.
The interface to Tpetra is implemented in the
``sundials::trilinos::nvector_tpetra::TpetraVectorInterface`` class. This class simply stores
a reference counting pointer to a Tpetra vector and inherits from
an empty structure

.. code-block:: c

   struct _N_VectorContent_Trilinos {};

to interface the C++ class with the NVECTOR C code.
A pointer to an instance of this class is kept in the *content* field
of the ``N_Vector`` object, to ensure that the Tpetra vector
is not deleted for as long as the ``N_Vector`` object exists.

The Tpetra vector type in the ``sundials::trilinos::nvector_tpetra::TpetraVectorInterface``
class is defined as:

.. cpp:type:: Tpetra::Vector<realtype, int, sunindextype> vector_type;

The Tpetra vector will use the SUNDIALS-specified ``realtype`` as its scalar
type, ``int`` as the local ordinal type, and ``sunindextype`` as the global ordinal type.
This type definition will use Tpetra's default node type. Available Kokkos node
types as of the Trilinos 12.14 release are serial (single thread), OpenMP, Pthread,
and CUDA. The default node type is selected when building the Kokkos package.
For example, the Tpetra vector will use a CUDA node if Tpetra was built with
CUDA support and the CUDA node was selected as the default when Tpetra was
built.

The header file to include when using this module is ``nvector_trilinos.h``.
The installed module library to link to is ``libsundials_nvectrilinos.lib``
where ``.lib`` is typically ``.so`` for shared libraries and ``.a``
for static libraries.


NVECTOR_TRILINOS functions
-----------------------------------

The NVECTOR_TRILINOS module defines implementations of all vector
operations listed in :numref:`NVectors.Ops`,
:numref:`NVectors.Ops.Fused`, :numref:`NVectors.Ops.Array`, and
:numref:`NVectors.Ops.Local`, except for
:c:func:`N_VGetArrayPointer` and :c:func:`N_VSetArrayPointer`.  As
such, this vector cannot be used with the SUNDIALS direct solvers
and preconditioners.  When access to raw
vector data is needed, it is recommended to extract the Trilinos
Tpetra vector first, and then use Tpetra vector methods to access the
data.  Usage examples of NVECTOR_TRILINOS are provided in example
programs for IDA.

The names of vector operations are obtained from those in
:numref:`NVectors.Ops` by appending the suffice ``_Trilinos``
(e.g. ``N_VDestroy_Trilinos``).  Vector operations call existing
``Tpetra::Vector`` methods when available. Vector operations specific
to SUNDIALS are implemented as standalone functions in the namespace
``sundials::trilinos::nvector_tpetra::TpetraVector``, located in the file ``SundialsTpetraVectorKernels.hpp``.
The module NVECTOR_TRILINOS provides the following additional user-callable routines:


.. cpp:function:: Teuchos::RCP<vector_type> N_VGetVector_Trilinos(N_Vector v)

   This C++ function takes an ``N_Vector`` as the argument and returns
   a reference counting pointer to the underlying Tpetra vector. This
   is a standalone function defined in the global namespace.


.. cpp:function:: N_Vector N_VMake_Trilinos(Teuchos::RCP<vector_type> v)

   This C++ function creates and allocates memory for an
   NVECTOR_TRILINOS wrapper around a user-provided Tpetra vector. This
   is a standalone function defined in the global namespace.



**Notes**

* The template parameter ``vector_type`` should be set as:

  .. code-block:: cpp

     typedef sundials::trilinos::nvector_tpetra::TpetraVectorInterface::vector_type vector_type

  This will ensure that data types used in Tpetra vector match those
  in SUNDIALS.

* When there is a need to access components of an ``N_Vector_Trilinos v``,
  it is recommeded to extract the Trilinos vector object via ``x_vec =
  N_VGetVector_Trilinos(v)`` and then access components using the
  appropriate Trilinos functions.

* The functions ``N_VDestroy_Trilinos`` and
  ``N_VDestroyVectorArray_Trilinos`` only delete the ``N_Vector``
  wrapper. The underlying Tpetra vector object will exist for as long
  as there is at least one reference to it.
