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

.. _NVectors.OpenMPDEV:

The NVECTOR_OPENMPDEV Module
============================

In situations where a user has access to a device such as a GPU for
offloading computation, SUNDIALS provides an NVECTOR implementation using
OpenMP device offloading, called NVECTOR_OPENMPDEV.

The NVECTOR_OPENMPDEV implementation defines the *content* field
of the ``N_Vector`` to be a structure  containing the length of the vector, a pointer
to the beginning of a contiguousdata array on the host, a pointer to the beginning of
a contiguous data array on the device, and a boolean flag ``own_data`` which specifies
the ownership of host and device data arrays.

.. code-block:: c

  struct _N_VectorContent_OpenMPDEV
  {
    sunindextype length;
    booleantype  own_data;
    realtype     *host_data;
    realtype     *dev_data;
  };

The header file to include when using this module is ``nvector_openmpdev.h``.
The installed module library to link to is ``libsundials_nvecopenmpdev.lib``
where ``.lib`` is typically ``.so`` for shared libraries and ``.a``
for static libraries.


NVECTOR_OPENMPDEV accessor macros
-----------------------------------

The following macros are provided to access the content of an NVECTOR_OPENMPDEV
vector.

.. c:macro:: NV_CONTENT_OMPDEV(v)

   This macro gives access to the contents of the NVECTOR_OPENMPDEV
   ``N_Vector v``.

   The assignment ``v_cont = NV_CONTENT_S(v)`` sets ``v_cont`` to be a
   pointer to the NVECTOR_OPENMPDEV  content structure.

   Implementation:

   .. code-block:: c

      #define NV_CONTENT_OMPDEV(v) ( (N_VectorContent_OpenMPDEV)(v->content) )


.. c:macro:: NV_OWN_DATA_OMPDEV(v)

   Access the *own_data* component of the OpenMPDEV ``N_Vector v``.

   The assignment ``v_data = NV_DATA_HOST_OMPDEV(v)`` sets ``v_data`` to be
   a pointer to the first component of the data on the host for the ``N_Vector v``.

   Implementation:

   .. code-block:: c

      #define NV_OWN_DATA_OMPDEV(v) ( NV_CONTENT_OMPDEV(v)->own_data )


.. c:macro:: NV_DATA_HOST_OMPDEV(v)

   The assignment ``NV_DATA_HOST_OMPDEV(v) = v_data`` sets the host component array
   of ``v`` to  be ``v_data`` by storing the pointer ``v_data``.

   Implementation:

   .. code-block:: c

      #define NV_DATA_HOST_OMPDEV(v) ( NV_CONTENT_OMPDEV(v)->host_data )


.. c:macro:: NV_DATA_DEV_OMPDEV(v)

   The assignment ``v_dev_data = NV_DATA_DEV_OMPDEV(v)`` sets ``v_dev_data`` to be
   a pointer to the first component of the data on the device for the ``N_Vector v``.
   The assignment ``NV_DATA_DEV_OMPDEV(v) = v_dev_data`` sets the device component
   array of ``v`` to be ``v_dev_data`` by storing the pointer ``v_dev_data``.

   Implementation:

   .. code-block:: c

      #define NV_DATA_DEV_OMPDEV(v) ( NV_CONTENT_OMPDEV(v)->dev_data )


.. c:macro:: NV_LENGTH_OMPDEV(V)

   Access the *length* component of the OpenMPDEV ``N_Vector v``.

   The assignment ``v_len = NV_LENGTH_OMPDEV(v)`` sets ``v_len`` to be
   the length of ``v``. On the other hand, the call ``NV_LENGTH_OMPDEV(v) = len_v``
   sets the length of ``v`` to be ``len_v``.

   .. code-block:: c

      #define NV_LENGTH_OMPDEV(v) ( NV_CONTENT_OMPDEV(v)->length )


NVECTOR_OPENMPDEV functions
-----------------------------------

The NVECTOR_OPENMPDEV module defines OpenMP device offloading implementations of all vector
operations listed in :numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`,
:numref:`NVectors.Ops.Array`, and :numref:`NVectors.Ops.Local`, except for
:c:func:`N_VSetArrayPointer`.
As such, this vector cannot be used with the SUNDIALS direct solvers and preconditioners.
It also provides methods for copying from the host to the device and vice versa.

The names of the vector operations are obtained from those in
:numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`, :numref:`NVectors.Ops.Array`, and
:numref:`NVectors.Ops.Local` by appending the suffix ``_OpenMPDEV`` (e.g.
``N_VDestroy_OpenMPDEV``).  The module NVECTOR_OPENMPDEV provides the following additional
user-callable routines:

.. c:function:: N_Vector N_VNew_OpenMPDEV(sunindextype vec_length, SUNContext sunctx)

   This function creates and allocates memory for an NVECTOR_OPENMPDEV ``N_Vector``.


.. c:function:: N_Vector N_VNewEmpty_OpenMPDEV(sunindextype vec_length, SUNContext sunctx)

   This function creates a new NVECTOR_OPENMPDEV ``N_Vector`` with an empty
   (``NULL``) data array.


.. c:function:: N_Vector N_VMake_OpenMPDEV(sunindextype vec_length, realtype *h_vdata, realtype *d_vdata, SUNContext sunctx)

   This function creates an NVECTOR_OPENMPDEV vector with user-supplied vector data
   arrays ``h_vdata`` and ``d_vdata``. This function does not allocate memory for
   data itself.


.. c:function:: realtype *N_VGetHostArrayPointer_OpenMPDEV(N_Vector v)

   This function returns a pointer to the host data array.


.. c:function:: realtype *N_VGetDeviceArrayPointer_OpenMPDEV(N_Vector v)

   This function returns a pointer to the device data array.


.. c:function:: void N_VPrint_OpenMPDEV(N_Vector v)

   This function prints the content of an NVECTOR_OPENMPDEV vector to ``stdout``.


.. c:function:: void N_VPrintFile_OpenMPDEV(N_Vector v, FILE *outfile)

   This function prints the content of an NVECTOR_OPENMPDEV vector to ``outfile``.


.. c:function:: void N_VCopyToDevice_OpenMPDEV(N_Vector v)

   This function copies the content of an NVECTOR_OPENMPDEV vector's host data array
   to the device data array.


.. c:function:: void N_VCopyFromDevice_OpenMPDEV(N_Vector v)

   This function copies the content of an NVECTOR_OPENMPDEV vector's device data array
   to the host data array.

By default all fused and vector array operations are disabled in the NVECTOR_OPENMPDEV
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with ``N_VNew_OpenMPDEV``, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using ``N_VClone``. This guarantees the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from while vectors created with
``N_VNew_OpenMPDEV`` will have the default settings for the NVECTOR_OPENMPDEV module.

.. c:function::  int N_VEnableFusedOps_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
  vector array operations in the NVECTOR_OPENMPDEV vector. The return value is ``0`` for
  success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.


.. c:function:: int N_VEnableLinearCombination_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
  combination fused operation in the NVECTOR_OPENMPDEV vector. The return value is ``0`` for
  success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.


.. c:function:: int N_VEnableScaleAddMulti_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
  add a vector to multiple vectors fused operation in the NVECTOR_OPENMPDEV vector. The
  return value is ``0`` for success and ``-1`` if the input vector or its
  ``ops`` structure are ``NULL``.


.. c:function:: int N_VEnableDotProdMulti_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
  dot products fused operation in the NVECTOR_OPENMPDEV vector. The return value is ``0``
  for success and ``-1`` if the input vector or its ``ops`` structure are
  ``NULL``.


.. c:function:: int N_VEnableLinearSumVectorArray_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
  operation for vector arrays in the NVECTOR_OPENMPDEV vector. The return value is ``0`` for
  success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.


.. c:function:: int N_VEnableScaleVectorArray_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
  operation for vector arrays in the NVECTOR_OPENMPDEV vector. The return value is ``0`` for
  success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.


.. c:function:: int N_VEnableConstVectorArray_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
  operation for vector arrays in the NVECTOR_OPENMPDEV vector. The return value is ``0`` for
  success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.


.. c:function:: int N_VEnableWrmsNormVectorArray_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
  operation for vector arrays in the NVECTOR_OPENMPDEV vector. The return value is ``0`` for
  success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.


.. c:function:: int N_VEnableWrmsNormMaskVectorArray_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
  norm operation for vector arrays in the NVECTOR_OPENMPDEV vector. The return value is
  ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
  ``NULL``.


.. c:function:: int N_VEnableScaleAddMultiVectorArray_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
  add a vector array to multiple vector arrays operation in the NVECTOR_OPENMPDEV vector. The
  return value is ``0`` for success and ``-1`` if the input vector or its
  ``ops`` structure are ``NULL``.


.. c:function:: int N_VEnableLinearCombinationVectorArray_OpenMPDEV(N_Vector v, booleantype tf)

  This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
  combination operation for vector arrays in the NVECTOR_OPENMPDEV vector. The return value
  is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
  are ``NULL``.


**Notes**

* When looping over the components of an ``N_Vector v``, it is
  most efficient to first obtain the component array via
  ``h_data = N_VGetArrayPointer(v)`` for the host array or
  ``v_data = N_VGetDeviceArrayPointer(v)`` for the device array,
  or equivalently to use the macros
  ``h_data = NV_DATA_HOST_OMPDEV(v)`` for the host array or
  ``v_data = NV_DATA_DEV_OMPDEV(v)`` for the device array, and then
  access ``h_data[i]`` or ``v_data[i]`` within the loop.

* When accessing individual components of an ``N_Vector v`` on
  the host remember to first copy the array
  back from the device with ``N_VCopyFromDevice_OpenMPDEV(v)``
  to ensure the array is up to date.

* :c:func:`N_VNewEmpty_OpenMPDEV`, :c:func:`N_VMake_OpenMPDEV`, and
  :c:func:`N_VCloneVectorArrayEmpty_OpenMPDEV()` set the field *own_data*
  to ``SUNFALSE``.  The functions :c:func:`N_VDestroy_OpenMPDEV()` and
  :c:func:`N_VDestroyVectorArray_OpenMPDEV()` will not attempt to free the
  pointer data for any ``N_Vector`` with *own_data* set to ``SUNFALSE``.
  In such a case, it is the user's responsibility to deallocate the
  data pointers.

* To maximize efficiency, vector operations in the NVECTOR_OPENMPDEV
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representation of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  length.
