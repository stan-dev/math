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


.. _NVectors.ParHyp:

The NVECTOR_PARHYP Module
=========================

The NVECTOR_PARHYP implementation of the NVECTOR  module provided with
SUNDIALS is a wrapper around HYPRE's ParVector class.
Most of the vector kernels simply call HYPRE vector operations.
The implementation defines the *content* field of ``N_Vector`` to
be a structure containing the global and local lengths of the vector, a
pointer to an object of type ``hypre_ParVector``, an MPI communicator,
and a boolean flag *own_parvector* indicating ownership of the
HYPRE parallel vector object *x*.


.. code-block:: c

   struct _N_VectorContent_ParHyp {
     sunindextype local_length;
     sunindextype global_length;
     booleantype own_data;
     booleantype own_parvector;
     realtype *data;
     MPI_Comm comm;
     hypre_ParVector *x;
   };

The header file to be included when using this module is ``nvector_parhyp.h``.
The installed module library to link to is
``libsundials_nvecparhyp.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

Unlike native SUNDIALS vector types, NVECTOR_PARHYP does not provide macros
to access its member variables.
Note that NVECTOR_PARHYP requires SUNDIALS to be built with MPI support.



NVECTOR_PARHYP functions
-----------------------------------

The NVECTOR_PARHYP module defines implementations of all vector operations
listed in :numref:`NVectors.Ops` except for :c:func:`N_VSetArrayPointer` and
:c:func:`N_VGetArrayPointer` because accessing raw vector data is handled by
low-level HYPRE functions.  As such, this vector is not available for use with
SUNDIALS Fortran interfaces.  When access to raw vector data is needed, one
should extract the HYPRE vector first, and then use HYPRE methods to access the
data.  Usage examples of NVECTOR_PARHYP are provided in the
``cvAdvDiff_non_ph.c`` example programs for CVODE and the
``ark_diurnal_kry_ph.c`` example program for ARKODE.

The names of parhyp methods are obtained from those in
:numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`, :numref:`NVectors.Ops.Array`, and
:numref:`NVectors.Ops.Local` by appending the suffix ``_ParHyp``
(e.g. ``N_VDestroy_ParHyp``).  The module NVECTOR_PARHYP provides the
following additional user-callable routines:


.. c:function:: N_Vector N_VNewEmpty_ParHyp(MPI_Comm comm, sunindextype local_length, sunindextype global_length, SUNContext sunctx)

   This function creates a new parhyp ``N_Vector`` with the pointer to the
   HYPRE vector set to ``NULL``.


.. c:function:: N_Vector N_VMake_ParHyp(hypre_ParVector *x, SUNContext sunctx)

   This function creates an ``N_Vector`` wrapper around an existing
   HYPRE parallel vector.  It does *not* allocate memory for ``x`` itself.


.. c:function:: hypre_ParVector *N_VGetVector_ParHyp(N_Vector v)

   This function returns a pointer to the underlying HYPRE vector.


.. c:function:: void N_VPrint_ParHyp(N_Vector v)

   This function prints the local content of a parhyp vector to ``stdout``.


.. c:function:: void N_VPrintFile_ParHyp(N_Vector v, FILE *outfile)

   This function prints the local content of a parhyp vector to ``outfile``.


By default all fused and vector array operations are disabled in the NVECTOR_PARHYP
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VMake_ParHyp`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from while vectors created with
:c:func:`N_VMake_ParHyp` will have the default settings for the NVECTOR_PARHYP module.

.. c:function:: int N_VEnableFusedOps_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the parhyp vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the parhyp vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the parhyp vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the parhyp vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the parhyp vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the parhyp vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the parhyp vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the parhyp vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the parhyp vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the parhyp vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_ParHyp(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the parhyp vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an ``N_Vector_ParHyp v``,
  it is recommended to extract the HYPRE vector via
  ``x_vec = N_VGetVector_ParHyp(v)`` and then access components using
  appropriate HYPRE functions.

* :c:func:`N_VNewEmpty_ParHyp`, :c:func:`N_VMake_ParHyp`, and
  :c:func:`N_VCloneVectorArrayEmpty_ParHyp()` set the field *own_parvector*
  to ``SUNFALSE``.  The functions :c:func:`N_VDestroy_ParHyp()` and
  :c:func:`N_VDestroyVectorArray_ParHyp()` will not attempt to delete an
  underlying HYPRE vector for any ``N_Vector`` with *own_parvector*
  set to ``SUNFALSE``.  In such a case, it is the user's responsibility
  to delete the underlying vector.

* To maximize efficiency, vector operations in the NVECTOR_PARHYP
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representations of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations.
