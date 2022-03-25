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

.. _NVectors.OpenMP:

The NVECTOR_OPENMP Module
=========================

In situations where a user has a multi-core processing unit capable of
running multiple parallel threads with shared memory, SUNDIALS provides
an implementation of NVECTOR using OpenMP, called NVECTOR_OPENMP, and
an implementation using Pthreads, called NVECTOR_PTHREADS. Testing has
shown that vectors should be of length at least :math:`100,000` before
the overhead associated with creating and using the threads is made up
by the parallelism in the vector calculations.

The OpenMP NVECTOR implementation provided with SUNDIALS,
NVECTOR_OPENMP, defines the *content* field of ``N_Vector`` to be a structure
containing the length of the vector, a pointer to the beginning of a contiguous
data array, a boolean flag *own_data* which specifies the ownership of
*data*, and the number of threads.  Operations on the vector are
threaded using OpenMP, the number of threads used is based on the
supplied argument in the vector constructor.

.. code-block:: c

   struct _N_VectorContent_OpenMP {
     sunindextype length;
     booleantype own_data;
     realtype *data;
     int num_threads;
   };

The header file to be included when using this module is ``nvector_openmp.h``.
The installed module library to link to is
``libsundials_nvecopenmp.lib`` where ``.lib`` is typically ``.so`` for shared libraries and ``.a``
for static libraries.
The Fortran module file to use when using the Fortran 2003 interface to
this module is ``fnvector_openmp_mod.mod``.


NVECTOR_OPENMP accessor macros
-----------------------------------

The following six macros are provided to access the content of an NVECTOR_OPENMP
vector. The suffix ``_OMP`` in the names denotes the OpenMP version.


.. c:macro:: NV_CONTENT_OMP(v)

   This macro gives access to the contents of the OpenMP vector
   ``N_Vector`` *v*.

   The assignment ``v_cont = NV_CONTENT_OMP(v)`` sets ``v_cont`` to be
   a pointer to the OpenMP ``N_Vector`` content structure.

   Implementation:

   .. code-block:: c

      #define NV_CONTENT_OMP(v) ( (N_VectorContent_OpenMP)(v->content) )


.. c:macro:: NV_OWN_DATA_OMP(v)

   Access the *own_data* component of the OpenMP ``N_Vector`` *v*.

   Implementation:

   .. code-block:: c

      #define NV_OWN_DATA_OMP(v) ( NV_CONTENT_OMP(v)->own_data )


.. c:macro:: NV_DATA_OMP(v)

   The assignment ``v_data = NV_DATA_OMP(v)`` sets ``v_data`` to be a
   pointer to the first component of the *data* for the ``N_Vector``
   ``v``.

   Similarly, the assignment ``NV_DATA_OMP(v) = v_data`` sets the component
   array of ``v`` to be ``v_data`` by storing the pointer ``v_data``.

   Implementation:

   .. code-block:: c

      #define NV_DATA_OMP(v) ( NV_CONTENT_OMP(v)->data )


.. c:macro:: NV_LENGTH_OMP(v)

   Access the *length* component of the OpenMP ``N_Vector`` *v*.

   The assignment ``v_len = NV_LENGTH_OMP(v)`` sets ``v_len`` to be the
   *length* of ``v``. On the other hand, the call ``NV_LENGTH_OMP(v) =
   len_v`` sets the *length* of ``v`` to be ``len_v``.

   Implementation:

   .. code-block:: c

      #define NV_LENGTH_OMP(v) ( NV_CONTENT_OMP(v)->length )


.. c:macro:: NV_NUM_THREADS_OMP(v)

   Access the *num_threads* component of the OpenMP ``N_Vector`` *v*.

   The assignment ``v_threads = NV_NUM_THREADS_OMP(v)`` sets
   ``v_threads`` to be the *num_threads* of ``v``. On the other hand,
   the call ``NV_NUM_THREADS_OMP(v) = num_threads_v`` sets the
   *num_threads* of ``v`` to be ``num_threads_v``.

   Implementation:

   .. code-block:: c

      #define NV_NUM_THREADS_OMP(v) ( NV_CONTENT_OMP(v)->num_threads )


.. c:macro:: NV_Ith_OMP(v,i)

   This macro gives access to the individual components of the *data*
   array of an ``N_Vector``, using standard 0-based C indexing.

   The assignment ``r = NV_Ith_OMP(v,i)`` sets ``r`` to be the value of
   the ``i``-th component of ``v``.

   The assignment ``NV_Ith_OMP(v,i) = r`` sets the value of the ``i``-th
   component of ``v`` to be ``r``.

   Here ``i`` ranges from 0 to :math:`n-1` for a vector of length
   :math:`n`.

   Implementation:

   .. code-block:: c

      #define NV_Ith_OMP(v,i) ( NV_DATA_OMP(v)[i] )



NVECTOR_OPENMP functions
-----------------------------------

The NVECTOR_OPENMP module defines OpenMP implementations of all vector
operations listed in :numref:`NVectors.Ops`,
:numref:`NVectors.Ops.Fused`, :numref:`NVectors.Ops.Array`, and
:numref:`NVectors.Ops.Local`.  Their names are obtained from those in
those sections by appending the suffix ``_OpenMP``
(e.g. ``N_VDestroy_OpenMP``).  All the standard vector operations
listed in :numref:`NVectors.Ops` with the suffix ``_OpenMP``
appended are callable via the Fortran 2003 interface by prepending an
`F' (e.g. ``FN_VDestroy_OpenMP``).

The module NVECTOR_OPENMP provides the following additional user-callable routines:


.. c:function:: N_Vector N_VNew_OpenMP(sunindextype vec_length, int num_threads, SUNContext sunctx)

   This function creates and allocates memory for a OpenMP
   ``N_Vector``. Arguments are the vector length and number of threads.


.. c:function:: N_Vector N_VNewEmpty_OpenMP(sunindextype vec_length, int num_threads, SUNContext sunctx)

   This function creates a new OpenMP ``N_Vector`` with an empty
   (``NULL``) data array.


.. c:function:: N_Vector N_VMake_OpenMP(sunindextype vec_length, realtype* v_data, int num_threads, SUNContext sunctx)

   This function creates and allocates memory for a OpenMP vector with
   user-provided data array, *v_data*.

   (This function does *not* allocate memory for ``v_data`` itself.)


.. c:function:: void N_VPrint_OpenMP(N_Vector v)

   This function prints the content of an OpenMP vector to ``stdout``.


.. c:function:: void N_VPrintFile_OpenMP(N_Vector v, FILE *outfile)

   This function prints the content of an OpenMP vector to ``outfile``.


By default all fused and vector array operations are disabled in the NVECTOR_OPENMP
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VNew_OpenMP`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from while vectors created with
:c:func:`N_VNew_OpenMP` will have the default settings for the NVECTOR_OPENMP module.

.. c:function:: int N_VEnableFusedOps_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the OpenMP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the OpenMP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the OpenMP vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the OpenMP vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the OpenMP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the OpenMP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the OpenMP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the OpenMP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the OpenMP vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the OpenMP vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_OpenMP(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the OpenMP vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When looping over the components of an ``N_Vector v``, it is more
  efficient to first obtain the component array via
  ``v_data = N_VGetArrayPointer(v)``, or equivalently
  ``v_data = NV_DATA_OMP(v)`` and then access ``v_data[i]`` within the
  loop than it is to use ``NV_Ith_OMP(v,i)`` within the loop.

* :c:func:`N_VNewEmpty_OpenMP`, :c:func:`N_VMake_OpenMP`, and
  :c:func:`N_VCloneVectorArrayEmpty_OpenMP()` set the field *own_data*
  to ``SUNFALSE``.  The functions :c:func:`N_VDestroy_OpenMP()` and
  :c:func:`N_VDestroyVectorArray_OpenMP()` will not attempt to free the
  pointer data for any ``N_Vector`` with *own_data* set to ``SUNFALSE``.
  In such a case, it is the user's responsibility to deallocate the
  data pointer.

* To maximize efficiency, vector operations in the NVECTOR_OPENMP
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representation of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  internal representations.


NVECTOR_OPENMP Fortran Interface
------------------------------------

The NVECTOR_OPENMP module provides a Fortran 2003 module for use from Fortran applications.

The ``fnvector_openmp_mod`` Fortran module defines interfaces to all
NVECTOR_OPENMP C functions using the intrinsic ``iso_c_binding``
module which provides a standardized mechanism for interoperating with C. As
noted in the C function descriptions above, the interface functions are
named after the corresponding C function, but with a leading ``F``. For
example, the function ``N_VNew_OpenMP`` is interfaced as
``FN_VNew_OpenMP``.

The Fortran 2003 NVECTOR_OPENMP interface module can be accessed with the ``use``
statement, i.e. ``use fnvector_openmp_mod``, and linking to the library
``libsundials_fnvectoropenmp_mod.lib`` in addition to the C library.
For details on where the library and module file
``fnvector_openmp_mod.mod`` are installed see :numref:`Installation`.
