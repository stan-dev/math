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

.. _NVectors.NVParallel:

The NVECTOR_PARALLEL Module
===========================

The NVECTOR_PARALLEL implementation of the NVECTOR module provided with
SUNDIALS is based on MPI.  It defines the *content* field of an
``N_Vector`` to be a structure containing the global and local lengths
of the vector, a pointer to the beginning of a contiguous local data
array, an MPI communicator, an a boolean flag *own_data* indicating
ownership of the data array *data*.

.. code-block:: c

   struct _N_VectorContent_Parallel {
      sunindextype local_length;
      sunindextype global_length;
      booleantype own_data;
      realtype *data;
      MPI_Comm comm;
   };

The header file to be included when using this module is
``nvector_parallel.h``. The installed module library to link to is
``libsundials_nvecparallel.lib`` where ``.lib`` is typically ``.so``
for shared libraries and ``.a`` for static libraries.


NVECTOR_PARALLEL accessor macros
-----------------------------------

The following seven macros are provided to access the content of a
NVECTOR_PARALLEL vector. The suffix ``_P`` in the names denotes the
distributed memory parallel version.

.. c:macro:: NV_CONTENT_P(v)

   This macro gives access to the contents of the parallel
   ``N_Vector`` *v*.

   The assignment ``v_cont = NV_CONTENT_P(v)`` sets ``v_cont`` to be a
   pointer to the ``N_Vector`` *content* structure of type ``struct
   N_VectorContent_Parallel``.

   Implementation:

   .. code-block:: c

      #define NV_CONTENT_P(v) ( (N_VectorContent_Parallel)(v->content) )


.. c:macro:: NV_OWN_DATA_P(v)

   Access the *own_data* component of the parallel ``N_Vector`` *v*.

   Implementation:

   .. code-block:: c

      #define NV_OWN_DATA_P(v)   ( NV_CONTENT_P(v)->own_data )

.. c:macro:: NV_DATA_P(v)

   The assignment ``v_data = NV_DATA_P(v)`` sets ``v_data`` to be a
   pointer to the first component of the *local_data* for the
   ``N_Vector v``.

   The assignment ``NV_DATA_P(v) = v_data`` sets the component array of
   ``v`` to be ``v_data`` by storing the pointer ``v_data`` into
   *data*.

   Implementation:

   .. code-block:: c

      #define NV_DATA_P(v)       ( NV_CONTENT_P(v)->data )

.. c:macro:: NV_LOCLENGTH_P(v)

   The assignment ``v_llen = NV_LOCLENGTH_P(v)`` sets ``v_llen`` to be
   the length of the local part of ``v``.

   The call ``NV_LOCLENGTH_P(v) = llen_v`` sets the *local_length* of
   ``v`` to be ``llen_v``.

   Implementation:

   .. code-block:: c

      #define NV_LOCLENGTH_P(v)  ( NV_CONTENT_P(v)->local_length )

.. c:macro:: NV_GLOBLENGTH_P(v)

   The assignment ``v_glen = NV_GLOBLENGTH_P(v)`` sets ``v_glen`` to be
   the *global_length* of the vector ``v``.

   The call ``NV_GLOBLENGTH_P(v) = glen_v`` sets the *global_length*
   of ``v`` to be ``glen_v``.

   Implementation:

   .. code-block:: c

      #define NV_GLOBLENGTH_P(v) ( NV_CONTENT_P(v)->global_length )

.. c:macro:: NV_COMM_P(v)

   This macro provides access to the MPI communicator used by the
   parallel ``N_Vector`` *v*.

   Implementation:

   .. code-block:: c

      #define NV_COMM_P(v) ( NV_CONTENT_P(v)->comm )

.. c:macro:: NV_Ith_P(v,i)

   This macro gives access to the individual components of the
   *local_data* array of an ``N_Vector``.

   The assignment ``r = NV_Ith_P(v,i)`` sets ``r`` to be the value of
   the ``i``-th component of the local part of ``v``.

   The assignment ``NV_Ith_P(v,i) = r`` sets the value of the ``i``-th
   component of the local part of ``v`` to be ``r``.

   Here ``i`` ranges from 0 to :math:`n-1`, where :math:`n` is the
   *local_length*.

   Implementation:

   .. code-block:: c

      #define NV_Ith_P(v,i) ( NV_DATA_P(v)[i] )



NVECTOR_PARALLEL functions
-----------------------------------

The NVECTOR_PARALLEL module defines parallel implementations of all
vector operations listed in :numref:`NVectors.Ops`.  Their names are
obtained from the generic names by appending the suffix ``_Parallel``
(e.g. ``N_VDestroy_Parallel``).   The module NVECTOR_PARALLEL provides
the following additional user-callable routines:


.. c:function:: N_Vector N_VNew_Parallel(MPI_Comm comm, sunindextype local_length, sunindextype global_length, SUNContext sunctx)

   This function creates and allocates memory for a parallel vector
   having global length *global_length*, having processor-local length
   *local_length*, and using the MPI communicator *comm*.


.. c:function:: N_Vector N_VNewEmpty_Parallel(MPI_Comm comm, sunindextype local_length, sunindextype global_length, SUNContext sunctx)

   This function creates a new parallel ``N_Vector`` with an empty
   (``NULL``) data array.


.. c:function:: N_Vector N_VMake_Parallel(MPI_Comm comm, sunindextype local_length, sunindextype global_length, realtype* v_data, SUNContext sunctx)

   This function creates and allocates memory for a parallel vector
   with user-provided data array.

   (This function does *not* allocate memory for ``v_data`` itself.)


.. c:function:: sunindextype N_VGetLocalLength_Parallel(N_Vector v)

   This function returns the local vector length.


.. c:function:: void N_VPrint_Parallel(N_Vector v)

   This function prints the local content of a parallel vector to ``stdout``.


.. c:function:: void N_VPrintFile_Parallel(N_Vector v, FILE *outfile)

   This function prints the local content of a parallel vector to ``outfile``.


By default all fused and vector array operations are disabled in the NVECTOR_PARALLEL
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VNew_Parallel`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees that the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from, while vectors created with
:c:func:`N_VNew_Parallel` will have the default settings for the NVECTOR_PARALLEL module.

.. c:function:: int N_VEnableFusedOps_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the parallel vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the parallel vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the parallel vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the parallel vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the parallel vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the parallel vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the parallel vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the parallel vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the parallel vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the parallel vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_Parallel(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the parallel vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When looping over the components of an ``N_Vector v``, it is
  more efficient to first obtain the local component array via
  ``v_data = N_VGetArrayPointer(v)``, or equivalently
  ``v_data = NV_DATA_P(v)``, and then access ``v_data[i]`` within the loop
  than it is to use ``NV_Ith_P(v,i)`` within the loop.

* :c:func:`N_VNewEmpty_Parallel`, :c:func:`N_VMake_Parallel`, and
  :c:func:`N_VCloneVectorArrayEmpty_Parallel()` set the field *own_data* to
  ``SUNFALSE``. The routines :c:func:`N_VDestroy_Parallel()` and
  :c:func:`N_VDestroyVectorArray_Parallel()` will not attempt to free the
  pointer data for any ``N_Vector`` with *own_data* set to
  ``SUNFALSE``. In such a case, it is the user's responsibility to
  deallocate the data pointer.

* To maximize efficiency, vector operations in the NVECTOR_PARALLEL
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representation of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations.



NVECTOR_PARALLEL Fortran Interface
------------------------------------

The NVECTOR_PARALLEL module provides a Fortran 2003 module for use from Fortran applications.

The ``fnvector_parallel_mod`` Fortran module defines interfaces to all
NVECTOR_PARALLEL C functions using the intrinsic ``iso_c_binding``
module which provides a standardized mechanism for interoperating with C. As
noted in the C function descriptions above, the interface functions are
named after the corresponding C function, but with a leading ``F``. For
example, the function ``N_VNew_Parallel`` is interfaced as
``FN_VNew_Parallel``.

The Fortran 2003 NVECTOR_PARALLEL interface module can be accessed with the ``use``
statement, i.e. ``use fnvector_parallel_mod``, and linking to the library
``libsundials_fnvectorparallel_mod.lib`` in addition to the C library.
For details on where the library and module file
``fnvector_parallel_mod.mod`` are installed see :numref:`Installation`.
We note that the module is accessible from the Fortran 2003 SUNDIALS integrators
*without* separately linking to the ``libsundials_fnvectorparallel_mod`` library.
