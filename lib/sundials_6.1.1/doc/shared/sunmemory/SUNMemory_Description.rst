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

.. _SUNMemory.Description:

The SUNMemoryHelper API
=======================

This API consists of three new SUNDIALS types: :c:type:`SUNMemoryType`,
:c:type:`SUNMemory`, and :c:type:`SUNMemoryHelper`:


.. c:type:: struct _SUNMemory *SUNMemory

   The ``SUNMemory`` type is a pointer a structure containing a pointer to
   actual data (``ptr``), the data memory type, and a flag indicating ownership
   of that data pointer. This structure is defined as

   .. code-block:: c

      struct _SUNMemory
      {
        void*         ptr;
        SUNMemoryType type;
        booleantype   own;
      };


.. c:enum:: SUNMemoryType

   The ``SUNMemoryType`` type is an enumeration that defines the supported
   memory types:

   .. code-block:: c

      typedef enum
      {
        SUNMEMTYPE_HOST,      /* pageable memory accessible on the host     */
        SUNMEMTYPE_PINNED,    /* page-locked memory accesible on the host   */
        SUNMEMTYPE_DEVICE,    /* memory accessible from the device          */
        SUNMEMTYPE_UVM        /* memory accessible from the host or device  */
      } SUNMemoryType;


.. c:type:: struct _SUNMemoryHelper *SUNMemoryHelper

   The ``SUNMemoryHelper`` type is a pointer to a structure containing a pointer
   to the implementation-specific member data (``content``) and a virtual method
   table of member functions (``ops``). This strucutre is defined as

   .. code-block:: c

      struct _SUNMemoryHelper
      {
        void*               content;
        SUNMemoryHelper_Ops ops;
      };


.. c:type:: struct _SUNMemoryHelper_Ops *SUNMemoryHelper_Ops

   The ``SUNMemoryHelper_Ops`` type is defined as a pointer to the structure
   containing the function pointers to the member function implementations. This
   structure is define as

   .. code-block:: c

      struct _SUNMemoryHelper_Ops
      {
        /* operations that implementations are required to provide */
        int (*alloc)(SUNMemoryHelper, SUNMemory* memptr size_t mem_size,
                     SUNMemoryType mem_type, void* queue);
        int (*dealloc)(SUNMemoryHelper, SUNMemory mem, void* queue);
        int (*copy)(SUNMemoryHelper, SUNMemory dst, SUNMemory src,
                    size_t mem_size, void* queue);

        /* operations that provide default implementations */
        int             (*copyasync)(SUNMemoryHelper, SUNMemory dst,
                                     SUNMemory src, size_t mem_size,
                                     void* queue);
        SUNMemoryHelper (*clone)(SUNMemoryHelper);
        int             (*destroy)(SUNMemoryHelper);
      };


.. _SUNMemory.Description.Required:

Implementation defined operations
---------------------------------

The SUNMemory API defines the following operations that an implementation to
must define:

.. c:function:: SUNMemory SUNMemoryHelper_Alloc(SUNMemoryHelper helper, \
                                                SUNMemory* memptr, \
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
   * ``mem_type`` -- the ``SUNMemoryType`` of the ``ptr``.
   * ``queue`` -- typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.

   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_Dealloc(SUNMemoryHelper helper, \
                                            SUNMemory mem, void* queue)

   Deallocates the ``mem->ptr`` field if it is owned by ``mem``, and then
   deallocates the ``mem`` object.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object.
   * ``mem`` -- the ``SUNMemory`` object.
   * ``queue`` -- typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.

   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).


.. c:function:: int SUNMemoryHelper_Copy(SUNMemoryHelper helper, \
                                         SUNMemory dst, SUNMemory src, \
                                         size_t mem_size, void* queue)

   Synchronously copies ``mem_size`` bytes from the the source memory to the
   destination memory.  The copy can be across memory spaces, e.g. host to
   device, or within a memory space, e.g. host to host.  The ``helper``
   object should use the memory types of ``dst`` and ``src`` to determine
   the appropriate transfer type necessary.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object.
   * ``dst`` -- the destination memory to copy to.
   * ``src`` -- the source memory to copy from.
   * ``mem_size`` -- the number of bytes to copy.
   * ``queue`` -- typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.

   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).



.. _SUNMemory.Description.Utilities:

Utility Functions
-----------------

The SUNMemoryHelper API defines the following functions which do not
require a SUNMemoryHelper instance:

.. c:function:: SUNMemory SUNMemoryHelper_Alias(SUNMemory mem1)

   Returns a ``SUNMemory`` object whose ``ptr`` field points to the same address
   as ``mem1``. The new object *will not* have ownership of ``ptr``, therefore,
   it will not free ``ptr`` when :c:func:`SUNMemoryHelper_Dealloc` is called.

   **Arguments:**

   * ``mem1`` -- a ``SUNMemory`` object.

   **Returns:**

   * A ``SUNMemory`` object or ``NULL`` if an error occurs.


.. c:function:: SUNMemory SUNMemoryHelper_Wrap(void* ptr, \
                                               SUNMemoryType mem_type)

   Returns a ``SUNMemory`` object whose ``ptr`` field points to the ``ptr``
   argument passed to the function. The new object *will not* have ownership of
   ``ptr``, therefore, it will not free ``ptr`` when
   :c:func:`SUNMemoryHelper_Dealloc` is called.

   **Arguments:**

   * ``ptr`` -- the data pointer to wrap in a ``SUNMemory`` object.
   * ``mem_type`` -- the ``SUNMemoryType`` of the ``ptr``.

   **Returns:**

   * A ``SUNMemory`` object or ``NULL`` if an error occurs.


.. c:function:: SUNMemoryHelper SUNMemoryHelper_NewEmpty()

   Returns an empty ``SUNMemoryHelper``. This is useful for building custom
   ``SUNMemoryHelper`` implementations.

   **Returns:**

   * A ``SUNMemoryHelper`` object or ``NULL`` if an error occurs.


.. c:function:: int SUNMemoryHelper_CopyOps(SUNMemoryHelper src, \
                                            SUNMemoryHelper dst)

   Copies the ``ops`` field of ``src`` to the ``ops`` field of ``dst``.
   This is useful for building custom ``SUNMemoryHelper`` implementations.

   **Arguments:**

   * ``src`` -- the object to copy from.
   * ``dst`` -- the object to copy to.

   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).



.. _SUNMemory.Description.Overridable:

Implementation overridable operations with defaults
---------------------------------------------------

In addition, the SUNMemoryHelper API defines the following *optionally
overridable* operations which an implementation may define:


.. c:function:: int SUNMemoryHelper_CopyAsync(SUNMemoryHelper helper, \
                                              SUNMemory dst, SUNMemory src, \
                                              size_t mem_size, void* queue)

   Asynchronously copies ``mem_size`` bytes from the the source memory to the
   destination memory.  The copy can be across memory spaces, e.g. host to
   device, or within a memory space, e.g. host to host.  The ``helper`` object
   should use the memory types of ``dst`` and ``src`` to determine the
   appropriate transfer type necessary.  The ``ctx`` argument is used when a
   different execution stream needs to be provided to perform the copy in,
   e.g. with ``CUDA`` this would be a ``cudaStream_t``.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object.
   * ``dst`` -- the destination memory to copy to.
   * ``src`` -- the source memory to copy from.
   * ``mem_size`` -- the number of bytes to copy.
   * ``queue`` -- typically a handle for an object representing an alternate
     execution stream (e.g., a CUDA/HIP stream or SYCL queue), but it can
     also be any implementation specific data.

   **Returns:**

   An ``int`` flag indicating success (zero) or failure (non-zero).

   .. note::

      If this operation is not defined by the implementation, then
      :c:func:`SUNMemoryHelper_Copy` will be used.


.. c:function:: SUNMemoryHelper SUNMemoryHelper_Clone(SUNMemoryHelper helper)

   Clones the ``SUNMemoryHelper`` object itself.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object to clone.

   **Returns:**

   * A ``SUNMemoryHelper`` object.

   .. note::

      If this operation is not defined by the implementation, then the default
      clone will only copy the ``SUNMemoryHelper_Ops`` structure stored in
      ``helper->ops``, and not the ``helper->content`` field.


.. c:function:: int SUNMemoryHelper_Destroy(SUNMemoryHelper helper)

   Destroys (frees) the ``SUNMemoryHelper`` object itself.

   **Arguments:**

   * ``helper`` -- the ``SUNMemoryHelper`` object to destroy.

   **Returns:**

   * An ``int`` flag indicating success (zero) or failure (non-zero).

   .. note::

      If this operation is not defined by the implementation, then the default
      destroy will only free the ``helper->ops`` field and the ``helper``
      itself. The ``helper->content`` field will not be freed.


.. _SUNMemory.Description.Custom:

Implementing a custom SUNMemoryHelper
-------------------------------------

A particular implementation of the SUNMemoryHelper API must:

*  Define and implement the required operations. Note that the names of
   these routines should be unique to that implementation in order to
   permit using more than one SUNMemoryHelper module in the same code.

*  Optionally, specify the *content* field of SUNMemoryHelper.

*  Optionally, define and implement additional user-callable routines
   acting on the newly defined SUNMemoryHelper.

An example of a custom SUNMemoryHelper is given in
``examples/utilities/custom_memory_helper.h``.
