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

.. _NVectors.HIP:

The NVECTOR_HIP Module
======================

The NVECTOR_HIP module is an NVECTOR implementation using the AMD ROCm HIP
library :cite:p:`rocm_site`. The module allows for SUNDIALS vector kernels
to run on AMD or NVIDIA GPU devices. It is intended for users who are already
familiar with HIP and GPU programming. Building this vector module requires
the HIP-clang compiler. The vector content layout is as follows:

.. code-block:: c++

   struct _N_VectorContent_Hip
   {
      sunindextype       length;
      booleantype        own_helper;
      SUNMemory          host_data;
      SUNMemory          device_data;
      SUNHipExecPolicy*  stream_exec_policy;
      SUNHipExecPolicy*  reduce_exec_policy;
      SUNMemoryHelper    mem_helper;
      void*              priv; /* 'private' data */
   };

   typedef struct _N_VectorContent_Hip *N_VectorContent_Hip;


The content members are the vector length (size), a boolean flag that signals if
the vector owns the data (i.e. it is in charge of freeing the data), pointers to
vector data on the host and the device, pointers to :cpp:type:`SUNHipExecPolicy`
implementations that control how the HIP kernels are launched for streaming and
reduction vector kernels, and a private data structure which holds additonal members
that should not be accessed directly.

When instantiated with :c:func:`N_VNew_Hip`, the underlying data will be
allocated on both the host and the device. Alternatively, a user can provide
host and device data arrays by using the :c:func:`N_VMake_Hip` constructor.
To use managed memory, the constructors :c:func:`N_VNewManaged_Hip` and
:c:func:`N_VMakeManaged_Hip` are provided. Additionally, a user-defined
``SUNMemoryHelper`` for allocating/freeing data can be provided with the
constructor :c:func:`N_VNewWithMemHelp_Hip`. Details on each of these
constructors are provided below.

To use the NVECTOR_HIP module, include ``nvector_hip.h`` and link to
the library ``libsundials_nvechip.lib``. The extension, ``.lib``, is
typically ``.so`` for shared libraries and ``.a`` for static libraries.


NVECTOR_HIP functions
-----------------------------------

Unlike other native SUNDIALS vector types, the NVECTOR_HIP module does not
provide macros to access its member variables. Instead, user should use the
accessor functions:


.. c:function:: realtype* N_VGetHostArrayPointer_Hip(N_Vector v)

   This function returns pointer to the vector data on the host.


.. c:function:: realtype* N_VGetDeviceArrayPointer_Hip(N_Vector v)

   This function returns pointer to the vector data on the device.


.. c:function:: booleantype N_VIsManagedMemory_Hip(N_Vector v)

   This function returns a boolean flag indiciating if the vector
   data array is in managed memory or not.


The NVECTOR_HIP module defines implementations of all standard vector
operations defined in :numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`,
:numref:`NVectors.Ops.Array`, and :numref:`NVectors.Ops.Local`, except for
:c:func:`N_VSetArrayPointer`.
The names of vector operations are obtained from those in
:numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`, :numref:`NVectors.Ops.Array`, and
:numref:`NVectors.Ops.Local` by appending the suffix ``_Hip``
(e.g. :c:func:`N_VDestroy_Hip`).  The module NVECTOR_HIP provides the
following additional user-callable routines:



.. c:function:: N_Vector N_VNew_Hip(sunindextype length, SUNContext sunctx)

   This function creates and allocates memory for a HIP ``N_Vector``.
   The vector data array is allocated on both the host and device.


.. c:function:: N_Vector N_VNewManaged_Hip(sunindextype vec_length, SUNContext sunctx)

   This function creates and allocates memory for a HIP
   ``N_Vector``. The vector data array is allocated in managed memory.


.. c:function:: N_Vector N_VNewWithMemHelp_Hip(sunindextype length, booleantype use_managed_mem, SUNMemoryHelper helper, SUNContext sunctx)

   This function creates a new HIP ``N_Vector`` with a user-supplied
   SUNMemoryHelper for allocating/freeing memory.


.. c:function:: N_Vector N_VNewEmpty_Hip(sunindextype vec_length, SUNContext sunctx)

   This function creates a new HIP ``N_Vector`` where the members of the content
   structure have not been allocated. This utility function is used by the
   other constructors to create a new vector.


.. c:function:: N_Vector N_VMake_Hip(sunindextype vec_length, realtype *h_vdata, realtype *d_vdata, SUNContext sunctx)


   This function creates a HIP ``N_Vector`` with user-supplied vector data arrays
   for the host and the device.


.. c:function:: N_Vector N_VMakeManaged_Hip(sunindextype vec_length, realtype *vdata, SUNContext sunctx)

   This function creates a HIP ``N_Vector`` with a user-supplied
   managed memory data array.



The module NVECTOR_HIP also provides the following user-callable routines:

.. c:function:: void N_VSetKernelExecPolicy_Hip(N_Vector v, SUNHipExecPolicy* stream_exec_policy, SUNHipExecPolicy* reduce_exec_policy)

   This function sets the execution policies which control the kernel parameters
   utilized when launching the streaming and reduction HIP kernels. By default
   the vector is setup to use the :cpp:func:`SUNHipThreadDirectExecPolicy` and
   :cpp:func:`SUNHipBlockReduceExecPolicy`. Any custom execution policy for
   reductions must ensure that the grid dimensions (number of thread blocks) is
   a multiple of the HIP warp size (32 for NVIDIA GPUs, 64 for AMD GPUs). See
   :numref:`NVectors.HIP.SUNHipExecPolicy` below for more information about the
   :cpp:type:`SUNHipExecPolicy` class. Providing ``NULL`` for an argument will
   result in the default policy being restored.

   .. note::

      Note: All vectors used in a single instance of a SUNDIALS package must use
      the same execution policy. It is **strongly recommended** that this
      function is called immediately after constructing the vector, and any
      subsequent vector be created by cloning to ensure consistent execution
      policies across vectors*


.. c:function:: realtype* N_VCopyToDevice_Hip(N_Vector v)

   This function copies host vector data to the device.


.. c:function:: realtype* N_VCopyFromDevice_Hip(N_Vector v)

   This function copies vector data from the device to the host.


.. c:function:: void N_VPrint_Hip(N_Vector v)

   This function prints the content of a HIP vector to ``stdout``.


.. c:function:: void N_VPrintFile_Hip(N_Vector v, FILE *outfile)

   This function prints the content of a HIP vector to ``outfile``.


By default all fused and vector array operations are disabled in the NVECTOR_HIP
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VNew_Hip`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from while vectors created with
:c:func:`N_VNew_Hip` will have the default settings for the NVECTOR_HIP module.

.. c:function:: int N_VEnableFusedOps_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the HIP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the HIP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the HIP vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the HIP vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the HIP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the HIP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the HIP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the HIP vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the HIP vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the HIP vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_Hip(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the HIP vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an ``N_Vector_Hip``, ``v``,
  it is recommeded to use functions :c:func:`N_VGetDeviceArrayPointer_Hip()` or
  :c:func:`N_VGetHostArrayPointer_Hip()`. However, when using managed memory,
  the function :c:func:`N_VGetArrayPointer` may also be used.

* To maximize efficiency, vector operations in the NVECTOR_HIP implementation
  that have more than one ``N_Vector`` argument do not check for
  consistent internal representations of these vectors. It is the user's
  responsibility to ensure that such routines are called with ``N_Vector``
  arguments that were all created with the same internal representations.


.. _NVectors.HIP.SUNHipExecPolicy:

The ``SUNHipExecPolicy`` Class
--------------------------------


In order to provide maximum flexibility to users, the HIP kernel execution parameters used
by kernels within SUNDIALS are defined by objects of the ``sundials::hip::ExecPolicy``
abstract class type (this class can be accessed in the global namespace as ``SUNHipExecPolicy``).
Thus, users may provide custom execution policies that fit the needs of their problem. The
``SUNHipExecPolicy`` class is defined as

.. cpp:type:: sundials::hip::ExecPolicy SUNHipExecPolicy

where the ``sundials::hip::ExecPolicy`` class is defined in the header file
``sundials_hip_policies.hpp``, as follows:

.. code-block:: c++

   class ExecPolicy
   {
   public:
      ExecPolicy(hipStream_t stream = 0) : stream_(stream) { }
      virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const = 0;
      virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const = 0;
      virtual const hipStream_t* stream() const { return (&stream_); }
      virtual ExecPolicy* clone() const = 0;
      ExecPolicy* clone_new_stream(hipStream_t stream) const {
         ExecPolicy* ex = clone();
         ex->stream_ = stream;
         return ex;
      }
      virtual bool atomic() const { return false; }
      virtual ~ExecPolicy() {}
   protected:
      hipStream_t stream_;
   };



To define a custom execution policy, a user simply needs to create a class that inherits from
the abstract class and implements the methods. The SUNDIALS provided
``sundials::hip::ThreadDirectExecPolicy`` (aka in the global namespace as
``SUNHipThreadDirectExecPolicy``) class is a good example of a what a custom execution policy
may look like:

.. code-block:: c++

   class ThreadDirectExecPolicy : public ExecPolicy
   {
   public:
      ThreadDirectExecPolicy(const size_t blockDim, hipStream_t stream = 0)
         : blockDim_(blockDim), ExecPolicy(stream)
      {}

      ThreadDirectExecPolicy(const ThreadDirectExecPolicy& ex)
         : blockDim_(ex.blockDim_), ExecPolicy(ex.stream_)
      {}

      virtual size_t gridSize(size_t numWorkUnits = 0, size_t /*blockDim*/ = 0) const
      {
         /* ceil(n/m) = floor((n + m - 1) / m) */
         return (numWorkUnits + blockSize() - 1) / blockSize();
      }

      virtual size_t blockSize(size_t /*numWorkUnits*/ = 0, size_t /*gridDim*/ = 0) const
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


In total, SUNDIALS provides 4 execution policies:


   .. cpp:function:: SUNHipThreadDirectExecPolicy(const size_t blockDim, const hipStream_t stream = 0)

      Maps each HIP thread to a work unit. The number of threads per block
      (blockDim) can be set to anything. The grid size will be calculated so
      that there are enough threads for one thread per element. If a HIP stream
      is provided, it will be used to execute the kernel.

   .. cpp:function:: SUNHipGridStrideExecPolicy(const size_t blockDim, const size_t gridDim, const hipStream_t stream = 0)

      Is for kernels that use grid stride loops. The number of threads per block (blockDim)
      can be set to anything. The number of blocks (gridDim) can be set to
      anything. If a HIP stream is provided, it will be used to execute the
      kernel.

   .. cpp:function:: SUNHipBlockReduceExecPolicy(const size_t blockDim, const hipStream_t stream = 0)

      Is for kernels performing a reduction across indvidual thread blocks. The
      number of threads per block (blockDim) can be set to any valid multiple of
      the HIP warp size. The grid size (gridDim) can be set to any value greater
      than 0. If it is set to 0, then the grid size will be chosen so that there
      is enough threads for one thread per work unit. If a HIP stream is
      provided, it will be used to execute the kernel.

   .. cpp:function:: SUNHipBlockReduceAtomicExecPolicy(const size_t blockDim, const hipStream_t stream = 0)

      Is for kernels performing a reduction across indvidual thread blocks using
      atomic operations. The number of threads per block (blockDim) can be set
      to any valid multiple of the HIP warp size. The grid size (gridDim) can be
      set to any value greater than 0. If it is set to 0, then the grid size
      will be chosen so that there is enough threads for one thread per work
      unit. If a HIP stream is provided, it will be used to execute the kernel.


For example, a policy that uses 128 threads per block and a user provided stream can be
created like so:

.. code-block:: c++

   hipStream_t stream;
   hipStreamCreate(&stream);
   SUNHipThreadDirectExecPolicy thread_direct(128, stream);


These default policy objects can be reused for multiple SUNDIALS data structures
(e.g. a :c:type:`SUNMatrix` and an :c:type:`N_Vector`) since they do not hold any
modifiable state information.
