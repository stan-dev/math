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

.. _NVectors.NVPETSc:

The NVECTOR_PETSC Module
========================

The NVECTOR_PETSC module is an NVECTOR wrapper around the PETSc vector. It
defines the *content* field of a ``N_Vector`` to be a structure
containing the global and local lengths of the vector, a pointer to
the PETSc vector, an MPI communicator, and a boolean flag  *own_data*
indicating ownership of the wrapped PETSc vector.

.. code-block:: c

   struct _N_VectorContent_Petsc {
      sunindextype local_length;
      sunindextype global_length;
      booleantype own_data;
      Vec *pvec;
      MPI_Comm comm;
   };

The header file to be included when using this module is
``nvector_petsc.h``.  The installed module library to link to is
``libsundials_nvecpetsc.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

Unlike native SUNDIALS vector types, NVECTOR_PETSC does not provide
macros to access its member variables.  Note that NVECTOR_PETSC
requires SUNDIALS to be built with MPI support.


NVECTOR_PETSC functions
-----------------------------------

The NVECTOR_PETSC module defines implementations of all vector operations listed
in :numref:`NVectors.Ops` except for :c:func:`N_VGetArrayPointer` and
:c:func:`N_VSetArrayPointer`.  As such, this vector cannot be used with SUNDIALS
Fortran interfaces.  When access to raw vector data is needed, it is recommended
to extract the PETSc vector first, and then use PETSc methods to access the
data.  Usage examples of NVECTOR_PETSC is provided in example programs for IDA.

The names of vector operations are obtained from those in
:numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`, :numref:`NVectors.Ops.Array`, and
:numref:`NVectors.Ops.Local` by appending the suffice ``_Petsc``
(e.g. ``N_VDestroy_Petsc``).  The module NVECTOR_PETSC provides the
following additional user-callable routines:


.. c:function:: N_Vector N_VNewEmpty_Petsc(MPI_Comm comm, sunindextype local_length, sunindextype global_length, SUNContext sunctx)

   This function creates a new PETSC ``N_Vector`` with the pointer to
   the wrapped PETSc vector set to ``NULL``. It is used by the
   ``N_VMake_Petsc`` and ``N_VClone_Petsc`` implementations.  It
   should be used only with great caution.


.. c:function:: N_Vector N_VMake_Petsc(Vec* pvec, SUNContext sunctx)

   This function creates and allocates memory for an NVECTOR_PETSC
   wrapper with a user-provided PETSc vector.  It does *not* allocate
   memory for the vector ``pvec`` itself.


.. c:function:: Vec *N_VGetVector_Petsc(N_Vector v)

   This function returns a pointer to the underlying PETSc vector.


.. c:function:: void N_VPrint_Petsc(N_Vector v)

   This function prints the global content of a wrapped PETSc vector to ``stdout``.


.. c:function:: void N_VPrintFile_Petsc(N_Vector v, const char fname[])

   This function prints the global content of a wrapped PETSc vector to ``fname``.


By default all fused and vector array operations are disabled in the NVECTOR_PETSC
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VMake_Petsc`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from while vectors created with
:c:func:`N_VMake_Petsc` will have the default settings for the NVECTOR_PETSC module.

.. c:function:: int N_VEnableFusedOps_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the PETSc vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the PETSc vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the PETSc vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the PETSc vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the PETSc vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the PETSc vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the PETSc vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the PETSc vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the PETSc vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the PETSc vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_Petsc(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the PETSc vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an ``N_Vector_Petsc v``, it
  is recommeded to extract the PETSc vector via ``x_vec = N_VGetVector_Petsc(v);``
  and then access components using appropriate PETSc functions.

* The functions :c:func:`N_VNewEmpty_Petsc`, :c:func:`N_VMake_Petsc`,
  and :c:func:`N_VCloneVectorArrayEmpty_Petsc()` set the field
  *own_data* to ``SUNFALSE``. The routines :c:func:`N_VDestroy_Petsc()` and
  :c:func:`N_VDestroyVectorArray_Petsc()` will not attempt to free the
  pointer ``pvec`` for any ``N_Vector`` with *own_data* set to
  ``SUNFALSE``. In such a case, it is the user's responsibility to
  deallocate the ``pvec`` pointer.

* To maximize efficiency, vector operations in the NVECTOR_PETSC
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representations of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations.
