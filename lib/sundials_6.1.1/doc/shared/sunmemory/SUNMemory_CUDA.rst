..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMemory.CUDA:

The SUNMemoryHelper_Cuda Implementation
=======================================

The SUNMemoryHelper_Cuda module is an implementation of the ``SUNMemoryHelper``
API that interfaces to the NVIDIA :cite:p:`cuda_site` library.  The
implementation defines the constructor

.. c:function:: SUNMemoryHelper SUNMemoryHelper_Cuda(SUNContext sunctx)

   Allocates and returns a ``SUNMemoryHelper`` object for handling CUDA memory
   if successful. Otherwise it returns ``NULL``.


.. _SUNMemory.CUDA.Operations:

SUNMemoryHelper_Cuda API Functions
----------------------------------

The implementation provides the following operations defined by the
``SUNMemoryHelper`` API:

.. c:function:: SUNMemory SUNMemoryHelper_Alloc_Cuda(SUNMemoryHelper helper, \
                                                     SUNMemory memptr, \
                                                     size_t mem_size, \
                                                     SUNMemoryType mem_type, \
                                                     void* queue)

   Allocates a ``SUNMemory`` object whose ``ptr`` field is allocated for
   ``mem_size`` bytes and is of type ``mem_type``. The new object will have
   ownership of ``ptr`` and will be deallocated when
   :c:func:`SUNMemoryHelper_Dealloc` is called.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object.
   * ``memptr`` -- pointer to the allocated ``SUNMemory``.
   * ``mem_size`` -- the size in bytes of the ``ptr``.
   * ``mem_type`` -- the ``SUNMemoryType`` of the ``ptr``. Supported values are:

     * ``SUNMEMTYPE_HOST`` -- memory is allocated with a call to ``malloc``.

     * ``SUNMEMTYPE_PINNED`` -- memory is allocated with a call to
       ``cudaMallocHost``.

     * ``SUNMEMTYPE_DEVICE`` -- memory is allocated with a call to
       ``cudaMalloc``.

     * ``SUNMEMTYPE_UVM`` -- memory is allocated with a call to
       ``cudaMallocManaged``.

   * ``queue`` -- currently unused.


   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_Dealloc_Cuda(SUNMemoryHelper helper, \
                                                 SUNMemory mem, void* queue)

   Deallocates the ``mem->ptr`` field if it is owned by ``mem``, and then
   deallocates the ``mem`` object.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object.
   * ``mem`` -- the ``SUNMemory`` object.
   * ``queue`` -- currently unused.

   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_Copy_Cuda(SUNMemoryHelper helper, \
                                              SUNMemory dst, SUNMemory src, \
                                              size_t mem_size, void* queue)

   Synchronously copies ``mem_size`` bytes from the the source memory to the
   destination memory.  The copy can be across memory spaces, e.g. host to
   device, or within a memory space, e.g. host to host.  The ``helper`` object
   will use the memory types of ``dst`` and ``src`` to determine the appropriate
   transfer type necessary.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object.
   * ``dst`` -- the destination memory to copy to.
   * ``src`` -- the source memory to copy from.
   * ``mem_size`` -- the number of bytes to copy.
   * ``queue`` -- currently unused.

   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_CopyAsync_Cuda(SUNMemoryHelper helper, \
                                                   SUNMemory dst, \
                                                   SUNMemory src, \
                                                   size_t mem_size, void* queue)

   Asynchronously copies ``mem_size`` bytes from the the source memory to the
   destination memory.  The copy can be across memory spaces, e.g. host to
   device, or within a memory space, e.g. host to host.  The ``helper`` object
   will use the memory types of ``dst`` and ``src`` to determine the appropriate
   transfer type necessary.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object.
   * ``dst`` -- the destination memory to copy to.
   * ``src`` -- the source memory to copy from.
   * ``mem_size`` -- the number of bytes to copy.
   * ``queue`` -- the ``cudaStream_t`` handle for the stream that the copy will
     be performed on.

   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).
