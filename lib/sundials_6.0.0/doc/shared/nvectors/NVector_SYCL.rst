..
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors.SYCL:

The NVECTOR_SYCL Module
=======================

The NVECTOR_SYCL module is an experimental NVECTOR implementation using the
`SYCL <https://www.khronos.org/sycl/>`_  abstraction layer. At present the only
supported SYCL compiler is the DPC++ (Intel oneAPI) compiler. This module allows
for SUNDIALS vector kernels to run on Intel GPU devices. The module is intended
for users who are already familiar with SYCL and GPU programming.

The vector content layout is as follows:

.. code-block:: c++

   struct _N_VectorContent_Sycl
   {
      sunindextype       length;
      booleantype        own_exec;
      booleantype        own_helper;
      SUNMemory          host_data;
      SUNMemory          device_data;
      SUNSyclExecPolicy* stream_exec_policy;
      SUNSyclExecPolicy* reduce_exec_policy;
      SUNMemoryHelper    mem_helper;
      sycl::queue*       queue;
      void*              priv; /* 'private' data */
   };

   typedef struct _N_VectorContent_Sycl *N_VectorContent_Sycl;


The content members are the vector length (size), boolean flags that indicate
if the vector owns the execution policies and memory helper objects (i.e., it is
in charge of freeing the objects), :c:type:`SUNMemory` objects for the vector data on
the host and device, pointers to execution policies that control how streaming
and reduction kernels are launched, a :c:type:`SUNMemoryHelper` for performing memory
operations, the SYCL queue, and a private data structure which holds additional
members that should not be accessed directly.

When instantiated with :cpp:func:`N_VNew_Sycl`, the underlying data will be
allocated on both the host and the device. Alternatively, a user can provide
host and device data arrays by using the :cpp:func:`N_VMake_Sycl` constructor.
To use managed (shared) memory, the constructors :cpp:func:`N_VNewManaged_Sycl`
and :cpp:func:`N_VMakeManaged_Sycl` are provided. Additionally, a user-defined
``SUNMemoryHelper`` for allocating/freeing data can be provided with the
constructor :cpp:func:`N_VNewWithMemHelp_Sycl`. Details on each of these
constructors are provided below.

The header file to include when using this is ``nvector_sycl.h``. The installed
module library to link to is ``libsundials_nvecsycl.lib``. The extension
``.lib`` is typically ``.so`` for shared libraries ``.a`` for static libraries.


NVECTOR_SYCL functions
-----------------------------------

The NVECTOR_SYCL module implementations of all vector operations listed in
:numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`,
:numref:`NVectors.Ops.Array`, and :numref:`NVectors.Ops.Local`, except for
:c:func:`N_VDotProdMulti()`, :c:func:`N_VWrmsNormVectorArray()`,
:c:func:`N_VWrmsNormMaskVectorArray()` as support for arrays of reduction
vectors is not yet supported.  These functions will be added to the NVECTOR_SYCL
implementation in the future. The names of vector operations are obtained from
those in the aforementioned sections by appending the suffix ``_Sycl`` (e.g.,
``N_VDestroy_Sycl``).

Additionally, the NVECTOR_SYCL module provides the following user-callable
constructors for creating a new NVECTOR_SYCL:


.. cpp:function:: N_Vector N_VNew_Sycl(sunindextype vec_length, sycl::queue* Q, SUNContext sunctx)

   This function creates and allocates memory for an NVECTOR_SYCL. Vector data
   arrays are allocated on both the host and the device associated with the
   input queue. All operation are launched in the provided queue.


.. cpp:function:: N_Vector N_VNewManaged_Sycl(sunindextype vec_length, sycl::queue* Q, SUNContext sunctx)

   This function creates and allocates memory for a NVECTOR_SYCL. The vector
   data array is allocated in managed (shared) memory using the input queue. All
   operation are launched in the provided queue.


.. cpp:function:: N_Vector N_VMake_Sycl(sunindextype length, realtype *h_vdata, realtype *d_vdata, sycl::queue* Q, SUNContext sunctx)

   This function creates an NVECTOR_SYCL with user-supplied host and device
   data arrays. This function does not allocate memory for data itself. All
   operation are launched in the provided queue.


.. cpp:function:: N_Vector N_VMakeManaged_Sycl(sunindextype length, realtype *vdata, sycl::queue *Q, SUNContext sunctx)

   This function creates an NVECTOR_SYCL with a user-supplied managed (shared)
   data array. This function does not allocate memory for data itself. All
   operation are launched in the provided queue.


.. cpp:function:: N_Vector N_VNewWithMemHelp_Sycl(sunindextype length, booleantype use_managed_mem, SUNMemoryHelper helper, sycl::queue *Q, SUNContext sunctx)

   This function creates an NVECTOR_SYCL with a user-supplied SUNMemoryHelper
   for allocating/freeing memory. All operation are launched in the provided
   queue.


.. cpp:function:: N_Vector N_VNewEmpty_Sycl()

   This function creates a new ``N_Vector`` where the members of the content
   structure have not been allocated.  This utility function is used by the
   other constructors to create a new vector.


The following user-callable functions are provided for accessing the vector data
arrays on the host and device and copying data between the two memory spaces.
Note the generic NVECTOR operations :c:func:`N_VGetArrayPointer()` and
:c:func:`N_VSetArrayPointer()` are mapped to the corresponding ``HostArray``
functions given below. To ensure memory coherency, a user will need to call the
``CopyTo`` or ``CopyFrom`` functions as necessary to transfer data between the
host and device, unless managed (shared) memory is used.


.. cpp:function:: realtype* N_VGetHostArrayPointer_Sycl(N_Vector v)

   This function returns a pointer to the vector host data array.


.. cpp:function:: realtype* N_VGetDeviceArrayPointer_Sycl(N_Vector v)

   This function returns a pointer to the vector device data array.


.. cpp:function:: void N_VSetHostArrayPointer_Sycl(realtype* h_vdata, N_Vector v)

   This function sets the host array pointer in the vector ``v``.


.. cpp:function:: void N_VSetDeviceArrayPointer_Sycl(realtype* d_vdata, N_Vector v)

   This function sets the device array pointer in the vector ``v``.


.. cpp:function:: void N_VCopyToDevice_Sycl(N_Vector v)

   This function copies host vector data to the device.


.. cpp:function:: void N_VCopyFromDevice_Sycl(N_Vector v)

   This function copies vector data from the device to the host.


.. cpp:function:: booleantype N_VIsManagedMemory_Sycl(N_Vector v)

   This function returns ``SUNTRUE`` if the vector data is allocated as managed
   (shared) memory otherwise it returns ``SUNFALSE``.


The following user-callable function is provided to set the execution policies
for how SYCL kernels are launched on a device.


.. cpp:function:: int N_VSetKernelExecPolicy_Sycl(N_Vector v, SUNSyclExecPolicy *stream_exec_policy, SUNSyclExecPolicy *reduce_exec_policy)

   This function sets the execution policies which control the kernel parameters
   utilized when launching the streaming and reduction kernels. By default the
   vector is setup to use the :cpp:func:`SUNSyclThreadDirectExecPolicy` and
   :cpp:func:`SUNSyclBlockReduceExecPolicy`. See
   :numref:`NVectors.SYCL.SUNSyclExecPolicy` below for more information about the
   :cpp:type:`SUNSyclExecPolicy` class.

   .. note::

      All vectors used in a single instance of a SUNDIALS package must use the
      same execution policy. It is **strongly recommended** that this function
      is called immediately after constructing the vector, and any subsequent
      vector be created by cloning to ensure consistent execution policies
      across vectors.


The following user-callable functions are provided to print the host vector data
array. Unless managed memory is used, a user may need to call
:cpp:func:`N_VCopyFromDevice_Sycl()` to ensure consistency between the host and
device array.


.. cpp:function:: void N_VPrint_Sycl(N_Vector v)

   This function prints the host data array to ``stdout``.


.. cpp:function:: void N_VPrintFile_Sycl(N_Vector v, FILE *outfile)

   This function prints the host data array to ``outfile``.


By default all fused and vector array operations are disabled in the
NVECTOR_SYCL module. The following additional user-callable routines are
provided to enable or disable fused and vector array operations for a specific
vector. To ensure consistency across vectors it is recommended to first create a
vector with one of the above constructors, enable/disable the desired operations
on that vector with the functions below, and then use this vector in conjunction
with :c:func:`N_VClone()` to create any additional vectors. This guarantees the
new vectors will have the same operations enabled/disabled as cloned vectors
inherit the same enable/disable options as the vector they are cloned from while
vectors created by any of the constructors above will have the default settings
for the NVECTOR_SYCL module.


.. cpp:function:: int N_VEnableFusedOps_Sycl(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the SYCL vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. cpp:function:: int N_VEnableLinearCombination_Sycl(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the SYCL vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. cpp:function:: int N_VEnableScaleAddMulti_Sycl(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the SYCL vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

..
   .. cpp:function:: int N_VEnableDotProdMulti_Sycl(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
      dot products fused operation in the SYCL vector. The return value is ``0``
      for success and ``-1`` if the input vector or its ``ops`` structure are
      ``NULL``.

.. cpp:function:: int N_VEnableLinearSumVectorArray_Sycl(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the SYCL vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. cpp:function:: int N_VEnableScaleVectorArray_Sycl(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the SYCL vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. cpp:function:: int N_VEnableConstVectorArray_Sycl(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the SYCL vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

..
   .. cpp:function:: int N_VEnableWrmsNormVectorArray_Sycl(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
      operation for vector arrays in the SYCL vector. The return value is ``0`` for
      success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

   .. cpp:function:: int N_VEnableWrmsNormMaskVectorArray_Sycl(N_Vector v, booleantype tf)

      This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
      norm operation for vector arrays in the SYCL vector. The return value is
      ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
      ``NULL``.

.. cpp:function:: int N_VEnableScaleAddMultiVectorArray_Sycl(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the SYCL vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. cpp:function:: int N_VEnableLinearCombinationVectorArray_Sycl(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the SYCL vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an NVECTOR_SYCL, ``v``, it is
  recommended to use :c:func:`N_VGetDeviceArrayPointer()` to access the device
  array or :c:func:`N_VGetArrayPointer()` for the host array. When using managed
  (shared) memory, either function may be used. To ensure memory coherency, a
  user may need to call the ``CopyTo`` or ``CopyFrom`` functions as necessary to
  transfer data between the host and device, unless managed (shared) memory is
  used.

* To maximize efficiency, vector operations in the NVECTOR_SYCL implementation
  that have more than one ``N_Vector`` argument do not check for consistent
  internal representations of these vectors. It is the user's responsibility to
  ensure that such routines are called with ``N_Vector`` arguments that were all
  created with the same internal representations.


.. _NVectors.SYCL.SUNSyclExecPolicy:

The ``SUNSyclExecPolicy`` Class
--------------------------------


In order to provide maximum flexibility to users, the SYCL kernel execution
parameters used by kernels within SUNDIALS are defined by objects of the
``sundials::sycl::ExecPolicy`` abstract class type (this class can be accessed in
the global namespace as ``SUNSyclExecPolicy``). Thus, users may provide custom
execution policies that fit the needs of their problem. The ``SUNSyclExecPolicy``
class is defined as

.. cpp:type:: sundials::sycl::ExecPolicy SUNSyclExecPolicy

where the ``sundials::sycl::ExecPolicy`` class is defined in the header file
``sundials_sycl_policies.hpp``, as follows:

.. code-block:: c++

   class ExecPolicy
   {
   public:
      virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const = 0;
      virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const = 0;
      virtual ExecPolicy* clone() const = 0;
      virtual ~ExecPolicy() {}
   };

For consistency the function names and behavior mirror the execution policies
for the CUDA and HIP vectors. In the SYCL case the ``blockSize`` is the local
work-group range in a one-dimensional ``nd_range`` (threads per group). The
``gridSize`` is the number of local work groups so the global work-group range
in a one-dimensional ``nd_range`` is ``blockSize * gridSize`` (total number of
threads). All vector kernels are written with a many-to-one mapping where work
units (vector elements) are mapped in a round-robin manner across the global
range. As such, the ``blockSize`` and ``gridSize`` can be set to any positive
value.

To define a custom execution policy, a user simply needs to create a class that
inherits from the abstract class and implements the methods. The SUNDIALS
provided ``sundials::sycl::ThreadDirectExecPolicy`` (aka in the global namespace
as ``SUNSyclThreadDirectExecPolicy``) class is a good example of a what a custom
execution policy may look like:

.. code-block:: c++

   class ThreadDirectExecPolicy : public ExecPolicy
   {
   public:
      ThreadDirectExecPolicy(const size_t blockDim)
         : blockDim_(blockDim)
      {}

      ThreadDirectExecPolicy(const ThreadDirectExecPolicy& ex)
         : blockDim_(ex.blockDim_)
      {}

      virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const
      {
         return (numWorkUnits + blockSize() - 1) / blockSize();
      }

      virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const
      {
         return blockDim_;
      }

      virtual ExecPolicy* clone() const
      {
         return static_cast<ExecPolicy*>(new ThreadDirectExecPolicy(*this));
      }

   private:
      const size_t blockDim_;
   };


SUNDIALS provides the following execution policies:

   .. cpp:function:: SUNSyclThreadDirectExecPolicy(const size_t blockDim)

      Is for kernels performing streaming operations and maps each work unit
      (vector element) to a work-item (thread). Based on the local work-group range
      (number of threads per group, ``blockSize``) the number of local work-groups
      (``gridSize``) is computed so there are enough work-items in the global
      work-group range ( total number of threads, ``blockSize * gridSize``) for one
      work unit per work-item (thread).

   .. cpp:function:: SUNSyclGridStrideExecPolicy(const size_t blockDim, const size_t gridDim)

      Is for kernels performing streaming operations and maps each work unit
      (vector element) to a work-item (thread) in a round-robin manner so the local
      work-group range (number of threads per group, ``blockSize``) and the number
      of local work-groups (``gridSize``) can be set to any positive value. In this
      case the global work-group range (total number of threads,
      ``blockSize * gridSize``) may be less than the number of work units (vector
      elements).

   .. cpp:function:: SUNSyclBlockReduceExecPolicy(const size_t blockDim)

      Is for kernels performing a reduction, the local work-group range (number
      of threads per group, ``blockSize``) and the number of local work-groups
      (``gridSize``) can be set to any positive value or the ``gridSize`` may be
      set to ``0`` in which case the global range is chosen so that there are
      enough threads for at most two work units per work-item.

By default the NVECTOR_SYCL module uses the ``SUNSyclThreadDirectExecPolicy``
and ``SUNSyclBlockReduceExecPolicy`` where the default ``blockDim`` is
determined by querying the device for the ``max_work_group_size``. User may
specify different policies by constructing a new ``SyclExecPolicy`` and
attaching it with :cpp:func:`N_VSetKernelExecPolicy_Sycl()`. For example, a policy
that uses 128 work-items (threads) per group can be created and attached like
so:

.. code-block:: c++

   N_Vector v = N_VNew_Sycl(length, SUNContext sunctx);
   SUNSyclThreadDirectExecPolicy thread_direct(128);
   SUNSyclBlockReduceExecPolicy  block_reduce(128);
   flag = N_VSetKernelExecPolicy_Sycl(v, &thread_direct, &block_reduce);


These default policy objects can be reused for multiple SUNDIALS data structures
(e.g. a :c:type:`SUNMatrix` and an :c:type:`N_Vector`) since they do not hold any modifiable
state information.
