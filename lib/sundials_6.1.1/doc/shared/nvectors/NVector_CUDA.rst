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

.. _NVectors.CUDA:

The NVECTOR_CUDA Module
=======================

The NVECTOR_CUDA module is an NVECTOR implementation in the CUDA language.
The module allows for SUNDIALS vector kernels to run on NVIDIA GPU devices. It is intended
for users who are already familiar with CUDA and GPU programming. Building this vector
module requires a CUDA compiler and, by extension, a C++ compiler. The vector content layout
is as follows:

.. code-block:: c++

   struct _N_VectorContent_Cuda
   {
      sunindextype       length;
      booleantype        own_helper;
      SUNMemory          host_data;
      SUNMemory          device_data;
      SUNCudaExecPolicy* stream_exec_policy;
      SUNCudaExecPolicy* reduce_exec_policy;
      SUNMemoryHelper    mem_helper;
      void*              priv; /* 'private' data */
   };

   typedef struct _N_VectorContent_Cuda *N_VectorContent_Cuda;


The content members are the vector length (size), boolean flags that indicate
if the vector owns the execution policies and memory helper objects (i.e., it is
in change of freeing the objects), :c:type:`SUNMemory` objects for the vector data on
the host and device, pointers to execution policies that control how streaming
and reduction kernels are launched, a :c:type:`SUNMemoryHelper` for performing memory
operations, and a private data structure which holds additonal members that
should not be accessed directly.

When instantiated with :c:func:`N_VNew_Cuda`, the underlying data will be
allocated on both the host and the device. Alternatively, a user can provide
host and device data arrays by using the :c:func:`N_VMake_Cuda` constructor.
To use CUDA managed memory, the constructors :c:func:`N_VNewManaged_Cuda` and
:c:func:`N_VMakeManaged_Cuda` are provided. Additionally, a user-defined
``SUNMemoryHelper`` for allocating/freeing data can be provided with the
constructor :c:func:`N_VNewWithMemHelp_Cuda`. Details on each of these
constructors are provided below.

To use the NVECTOR_CUDA module, include ``nvector_cuda.h`` and link to
the library ``libsundials_nveccuda.lib``. The extension, ``.lib``, is
typically ``.so`` for shared libraries and ``.a`` for static libraries.


.. _NVectors.CUDA.Functions:

NVECTOR_CUDA functions
-----------------------------------

Unlike other native SUNDIALS vector types, the NVECTOR_CUDA module does not
provide macros to access its member variables. Instead, user should use the
accessor functions:


.. c:function:: realtype* N_VGetHostArrayPointer_Cuda(N_Vector v)

   This function returns pointer to the vector data on the host.


.. c:function:: realtype* N_VGetDeviceArrayPointer_Cuda(N_Vector v)

   This function returns pointer to the vector data on the device.


.. c:function:: booleantype N_VIsManagedMemory_Cuda(N_Vector v)

   This function returns a boolean flag indiciating if the vector
   data array is in managed memory or not.


The NVECTOR_CUDA module defines implementations of all standard vector
operations defined in :numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`,
:numref:`NVectors.Ops.Array`, and :numref:`NVectors.Ops.Local`, except for
:c:func:`N_VSetArrayPointer`, and, if using unmanaged memory,
:c:func:`N_VGetArrayPointer`.  As such, this vector can only be used with
SUNDIALS direct solvers and preconditioners when using managed memory.
The NVECTOR_CUDA module provides separate functions to access data on the host
and on the device for the unmanaged memory use case. It also provides methods for
copying from the host to the device and vice versa. Usage examples of NVECTOR_CUDA
are provided in example programs for CVODE :cite:p:`cvode_ex`.

The names of vector operations are obtained from those in
:numref:`NVectors.Ops`, :numref:`NVectors.Ops.Fused`, :numref:`NVectors.Ops.Array`, and
:numref:`NVectors.Ops.Local` by appending the suffix ``_Cuda``
(e.g. ``N_VDestroy_Cuda``).  The module NVECTOR_CUDA provides the
following additional user-callable routines:



.. c:function:: N_Vector N_VNew_Cuda(sunindextype length, SUNContext sunctx)

   This function creates and allocates memory for a CUDA ``N_Vector``.
   The vector data array is allocated on both the host and device.


.. c:function:: N_Vector N_VNewManaged_Cuda(sunindextype vec_length, SUNContext sunctx)

   This function creates and allocates memory for a CUDA
   ``N_Vector``. The vector data array is allocated in managed memory.


.. c:function:: N_Vector N_VNewWithMemHelp_Cuda(sunindextype length, booleantype use_managed_mem, SUNMemoryHelper helper, SUNContext sunctx)

   This function creates a new CUDA ``N_Vector`` with a user-supplied
   SUNMemoryHelper for allocating/freeing memory.


.. c:function:: N_Vector N_VNewEmpty_Cuda(sunindextype vec_length, SUNContext sunctx)

   This function creates a new CUDA ``N_Vector`` where the members of the
   content structure have not been allocated. This utility function is used by
   the other constructors to create a new vector.


.. c:function:: N_Vector N_VMake_Cuda(sunindextype vec_length, realtype *h_vdata, realtype *d_vdata, SUNContext sunctx)


   This function creates a CUDA ``N_Vector`` with user-supplied vector data arrays
   for the host and the device.


.. c:function:: N_Vector N_VMakeManaged_Cuda(sunindextype vec_length, realtype *vdata, SUNContext sunctx)

   This function creates a CUDA ``N_Vector`` with a user-supplied
   managed memory data array.


.. c:function:: N_Vector N_VMakeWithManagedAllocator_Cuda(sunindextype length, void* (*allocfn)(size_t size), void (*freefn)(void* ptr))

   This function creates a CUDA ``N_Vector`` with a user-supplied memory allocator.
   It requires the user to provide a corresponding free function as well.
   The memory allocated by the allocator function must behave like CUDA managed memory.



The module NVECTOR_CUDA also provides the following user-callable routines:

.. c:function:: void N_VSetKernelExecPolicy_Cuda(N_Vector v, SUNCudaExecPolicy* stream_exec_policy, SUNCudaExecPolicy* reduce_exec_policy)

   This function sets the execution policies which control the kernel parameters
   utilized when launching the streaming and reduction CUDA kernels. By default
   the vector is setup to use the :cpp:func:`SUNCudaThreadDirectExecPolicy` and
   :cpp:func:`SUNCudaBlockReduceAtomicExecPolicy`. Any custom execution policy
   for reductions must ensure that the grid dimensions (number of thread blocks)
   is a multiple of the CUDA warp size (32). See
   :numref:`NVectors.CUDA.SUNCudaExecPolicy` below for more information about
   the :cpp:type:`SUNCudaExecPolicy` class. Providing ``NULL`` for an argument
   will result in the default policy being restored.

   .. note::

      Note: All vectors used in a single instance of a SUNDIALS package must use
      the same execution policy. It is **strongly recommended** that this
      function is called immediately after constructing the vector, and any
      subsequent vector be created by cloning to ensure consistent execution
      policies across vectors


.. c:function:: realtype* N_VCopyToDevice_Cuda(N_Vector v)

   This function copies host vector data to the device.


.. c:function:: realtype* N_VCopyFromDevice_Cuda(N_Vector v)

   This function copies vector data from the device to the host.


.. c:function:: void N_VPrint_Cuda(N_Vector v)

   This function prints the content of a CUDA vector to ``stdout``.


.. c:function:: void N_VPrintFile_Cuda(N_Vector v, FILE *outfile)

   This function prints the content of a CUDA vector to ``outfile``.


By default all fused and vector array operations are disabled in the NVECTOR_CUDA
module. The following additional user-callable routines are provided to
enable or disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a vector
with :c:func:`N_VNew_Cuda`, enable/disable the desired operations for that vector
with the functions below, and create any additional vectors from that vector
using :c:func:`N_VClone`. This guarantees the new vectors will have the same
operations enabled/disabled as cloned vectors inherit the same enable/disable
options as the vector they are cloned from while vectors created with
:c:func:`N_VNew_Cuda` will have the default settings for the NVECTOR_CUDA module.

.. c:function:: int N_VEnableFusedOps_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the CUDA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the CUDA vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the CUDA vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the CUDA vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableScaleAddMultiVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector array to multiple vector arrays operation in the CUDA vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombinationVectorArray_Cuda(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination operation for vector arrays in the CUDA vector. The return value
   is ``0`` for success and ``-1`` if the input vector or its ``ops`` structure
   are ``NULL``.


**Notes**

* When there is a need to access components of an ``N_Vector_Cuda``, ``v``,
  it is recommeded to use functions :c:func:`N_VGetDeviceArrayPointer_Cuda()` or
  :c:func:`N_VGetHostArrayPointer_Cuda()`. However, when using managed memory,
  the function :c:func:`N_VGetArrayPointer` may also be used.

* To maximize efficiency, vector operations in the NVECTOR_CUDA implementation
  that have more than one ``N_Vector`` argument do not check for
  consistent internal representations of these vectors. It is the user's
  responsibility to ensure that such routines are called with ``N_Vector``
  arguments that were all created with the same internal representations.


.. _NVectors.CUDA.SUNCudaExecPolicy:

The ``SUNCudaExecPolicy`` Class
--------------------------------


In order to provide maximum flexibility to users, the CUDA kernel execution parameters used
by kernels within SUNDIALS are defined by objects of the ``sundials::cuda::ExecPolicy``
abstract class type (this class can be accessed in the global namespace as ``SUNCudaExecPolicy``).
Thus, users may provide custom execution policies that fit the needs of their problem. The
``SUNCudaExecPolicy`` class is defined as

.. cpp:type:: sundials::cuda::ExecPolicy SUNCudaExecPolicy

where the ``sundials::cuda::ExecPolicy`` class is defined in the header file
``sundials_cuda_policies.hpp``, as follows:

.. code-block:: c++

   class ExecPolicy
   {
   public:
      ExecPolicy(cudaStream_t stream = 0) : stream_(stream) { }
      virtual size_t gridSize(size_t numWorkUnits = 0, size_t blockDim = 0) const = 0;
      virtual size_t blockSize(size_t numWorkUnits = 0, size_t gridDim = 0) const = 0;
      virtual const cudaStream_t* stream() const { return (&stream_); }
      virtual ExecPolicy* clone() const = 0;
      ExecPolicy* clone_new_stream(cudaStream_t stream) const {
         ExecPolicy* ex = clone();
         ex->stream_ = stream;
         return ex;
      }
      virtual bool atomic() const { return false; }
      virtual ~ExecPolicy() {}
   protected:
      cudaStream_t stream_;
   };


To define a custom execution policy, a user simply needs to create a class that
inherits from the abstract class and implements the methods. The SUNDIALS
provided ``sundials::cuda::ThreadDirectExecPolicy`` (aka in the global namespace
as ``SUNCudaThreadDirectExecPolicy``) class is a good example of a what a custom
execution policy may look like:

.. code-block:: c++

   class ThreadDirectExecPolicy : public ExecPolicy
   {
   public:
      ThreadDirectExecPolicy(const size_t blockDim, cudaStream_t stream = 0)
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


In total, SUNDIALS provides 3 execution policies:

   .. cpp:function:: SUNCudaThreadDirectExecPolicy(const size_t blockDim, const cudaStream_t stream = 0)

      Maps each CUDA thread to a work unit. The number of threads per block
      (blockDim) can be set to anything. The grid size will be calculated so
      that there are enough threads for one thread per element. If a CUDA stream
      is provided, it will be used to execute the kernel.

   .. cpp:function:: SUNCudaGridStrideExecPolicy(const size_t blockDim, const size_t gridDim, const cudaStream_t stream = 0)

      Is for kernels that use grid stride loops. The number of threads per block
      (blockDim) can be set to anything. The number of blocks (gridDim) can be
      set to anything. If a CUDA stream is provided, it will be used to execute
      the kernel.

   .. cpp:function:: SUNCudaBlockReduceExecPolicy(const size_t blockDim, const cudaStream_t stream = 0)

      Is for kernels performing a reduction across indvidual thread blocks. The
      number of threads per block (blockDim) can be set to any valid multiple of
      the CUDA warp size. The grid size (gridDim) can be set to any value
      greater than 0. If it is set to 0, then the grid size will be chosen so
      that there is enough threads for one thread per work unit. If a CUDA
      stream is provided, it will be used to execute the kernel.

   .. cpp:function:: SUNCudaBlockReduceAtomicExecPolicy(const size_t blockDim, const cudaStream_t stream = 0)

      Is for kernels performing a reduction across indvidual thread blocks using
      atomic operations. The number of threads per block (blockDim) can be set
      to any valid multiple of the CUDA warp size. The grid size (gridDim) can be
      set to any value greater than 0. If it is set to 0, then the grid size
      will be chosen so that there is enough threads for one thread per work
      unit. If a CUDA stream is provided, it will be used to execute the kernel.


For example, a policy that uses 128 threads per block and a user provided stream can be
created like so:

.. code-block:: c++

   cudaStream_t stream;
   cudaStreamCreate(&stream);
   SUNCudaThreadDirectExecPolicy thread_direct(128, stream);


These default policy objects can be reused for multiple SUNDIALS data structures
(e.g. a :c:type:`SUNMatrix` and an :c:type:`N_Vector`) since they do not hold any
modifiable state information.
