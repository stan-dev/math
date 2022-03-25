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

.. _NVectors.NVSerial:

The NVECTOR_SERIAL Module
=========================

The serial implementation of the NVECTOR module provided with
SUNDIALS, NVECTOR_SERIAL, defines the *content* field of an
``N_Vector`` to be a structure containing the length of the vector, a
pointer to the beginning of a contiguous data array, and a boolean
flag *own_data* which specifies the ownership of data.

.. code-block:: c

   struct _N_VectorContent_Serial {
      sunindextype length;
      booleantype own_data;
      realtype *data;
   };

The header file to be included when using this module is ``nvector_serial.h``.
The installed module library to link to is
``libsundials_nvecserial.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

.. _NVectors.NVSerial.Macros:

NVECTOR_SERIAL accessor macros
------------------------------

The following five macros are provided to access the content of an
NVECTOR_SERIAL vector. The suffix ``_S`` in the names denotes the serial
version.


.. c:macro:: NV_CONTENT_S(v)

   This macro gives access to the contents of the serial vector
   ``N_Vector`` *v*.

   The assignment ``v_cont = NV_CONTENT_S(v)`` sets ``v_cont`` to be a
   pointer to the serial ``N_Vector`` `content` structure.

   Implementation:

   .. code-block:: c

      #define NV_CONTENT_S(v) ( (N_VectorContent_Serial)(v->content) )


.. c:macro:: NV_OWN_DATA_S(v)

   Access the *own_data* component of the serial ``N_Vector`` *v*.

   Implementation:

   .. code-block:: c

      #define NV_OWN_DATA_S(v) ( NV_CONTENT_S(v)->own_data )


.. c:macro:: NV_DATA_S(v)

   The assignment ``v_data = NV_DATA_S(v)`` sets ``v_data`` to be a
   pointer to the first component of the *data* for the ``N_Vector``
   ``v``.

   Similarly, the assignment ``NV_DATA_S(v) = v_data`` sets the component
   array of ``v`` to be ``v_data`` by storing the pointer ``v_data``.

   Implementation:

   .. code-block:: c

      #define NV_DATA_S(v) ( NV_CONTENT_S(v)->data )


.. c:macro:: NV_LENGTH_S(v)

   Access the *length* component of the serial ``N_Vector`` *v*.

   The assignment ``v_len = NV_LENGTH_S(v)`` sets ``v_len`` to be the
   *length* of ``v``. On the other hand, the call ``NV_LENGTH_S(v) =
   len_v`` sets the *length* of ``v`` to be ``len_v``.

   Implementation:

   .. code-block:: c

      #define NV_LENGTH_S(v) ( NV_CONTENT_S(v)->length )


.. c:macro:: NV_Ith_S(v,i)

   This macro gives access to the individual components of the *data*
   array of an ``N_Vector``, using standard 0-based C indexing.

   The assignment ``r = NV_Ith_S(v,i)`` sets ``r`` to be the value of
   the ``i``-th component of ``v``.

   The assignment ``NV_Ith_S(v,i) = r`` sets the value of the ``i``-th
   component of ``v`` to be ``r``.

   Here ``i`` ranges from 0 to :math:`n-1` for a vector of length
   :math:`n`.

   Implementation:

   .. code-block:: c

      #define NV_Ith_S(v,i) ( NV_DATA_S(v)[i] )


.. _NVectors.NVSerial.Functions:

NVECTOR_SERIAL functions
------------------------

The NVECTOR_SERIAL module defines serial implementations of all vector
operations listed in :numref:`NVectors.Ops.Standard`,
:numref:`NVectors.Ops.Fused`, :numref:`NVectors.Ops.Array`, and
:numref:`NVectors.Ops.Local`.  Their names are obtained from those in
those sections by appending the suffix ``_Serial``
(e.g. ``N_VDestroy_Serial``).  All the standard vector operations
listed in :numref:`NVectors.Ops.Standard` with the suffix ``_Serial``
appended are callable via the Fortran 2003 interface by prepending an
``F`` (e.g. ``FN_VDestroy_Serial``).

The module NVECTOR_SERIAL provides the following additional
user-callable routines:

.. c:function:: N_Vector N_VNew_Serial(sunindextype vec_length, SUNContext sunctx)

   This function creates and allocates memory for a serial
   ``N_Vector``. Its only argument is the vector length.


.. c:function:: N_Vector N_VNewEmpty_Serial(sunindextype vec_length, SUNContext sunctx)

   This function creates a new serial ``N_Vector`` with an empty
   (``NULL``) data array.


.. c:function:: N_Vector N_VMake_Serial(sunindextype vec_length, realtype* v_data, SUNContext sunctx)

   This function creates and allocates memory for a serial vector with
   user-provided data array, *v_data*.

   (This function does *not* allocate memory for ``v_data`` itself.)


.. c:function:: void N_VPrint_Serial(N_Vector v)

   This function prints the content of a serial vector to ``stdout``.


.. c:function:: void N_VPrintFile_Serial(N_Vector v, FILE *outfile)

   This function prints the content of a serial vector to ``outfile``.


By default all fused and vector array operations are disabled in the NVECTOR_SERIAL
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VNew_Serial`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees that the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned, from while vectors created with
:c:func:`N_VNew_Serial` will have the default settings for the NVECTOR_SERIAL module.

.. c:function:: int N_VEnableFusedOps_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the serial vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the serial vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the serial vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the serial vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the serial vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the serial vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the serial vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the serial vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the serial vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the serial vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_Serial(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the serial vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When looping over the components of an ``N_Vector v``, it is more
  efficient to first obtain the component array via ``v_data =
  NV_DATA_S(v)``, or equivalently ``v_data = N_VGetArrayPointer(v)``,
  and then access ``v_data[i]`` within the loop than it
  is to use ``NV_Ith_S(v,i)`` within the loop.

* :c:func:`N_VNewEmpty_Serial`, :c:func:`N_VMake_Serial`, and
  :c:func:`N_VCloneVectorArrayEmpty_Serial()` set the field *own_data*
  to ``SUNFALSE``.  The functions :c:func:`N_VDestroy_Serial()` and
  :c:func:`N_VDestroyVectorArray_Serial()` will not attempt to free the
  pointer data for any ``N_Vector`` with *own_data* set to ``SUNFALSE``.
  In such a case, it is the user's responsibility to deallocate the
  data pointer.

* To maximize efficiency, vector operations in the NVECTOR_SERIAL
  implementation that have more than one ``N_Vector`` argument do not
  check for consistent internal representation of these vectors. It is
  the user's responsibility to ensure that such routines are called
  with ``N_Vector`` arguments that were all created with the same
  length.


.. _NVectors.NVSerial.Fortran:

NVECTOR_SERIAL Fortran Interface
------------------------------------

The NVECTOR_SERIAL module provides a Fortran 2003 module for use from Fortran applications.

The ``fnvector_serial_mod`` Fortran module defines interfaces to all
NVECTOR_SERIAL C functions using the intrinsic ``iso_c_binding``
module which provides a standardized mechanism for interoperating with C. As
noted in the C function descriptions above, the interface functions are
named after the corresponding C function, but with a leading ``F``. For
example, the function ``N_VNew_Serial`` is interfaced as
``FN_VNew_Serial``.

The Fortran 2003 NVECTOR_SERIAL interface module can be accessed with the ``use``
statement, i.e. ``use fnvector_serial_mod``, and linking to the library
``libsundials_fnvectorserial_mod.lib`` in addition to the C library.
For details on where the library and module file
``fnvector_serial_mod.mod`` are installed see :numref:`Installation`.
We note that the module is accessible from the Fortran 2003 SUNDIALS integrators
*without* separately linking to the ``libsundials_fnvectorserial_mod`` library.
