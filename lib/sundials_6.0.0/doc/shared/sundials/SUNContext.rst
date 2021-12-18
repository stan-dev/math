.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDIALS.SUNContext:

The SUNContext Type
=====================

In SUNDIALS v6.0.0, the concept of a SUNDIALS context object was introduced.
Every SUNDIALS object holds a reference to a common ``SUNContext`` object, which
is in turn unique to a thread of execution.

The ``SUNContext`` type is defined in the header file ``sundials/sundials_context.h`` as

.. c:type:: struct _SUNContext *SUNContext

Users should create a ``SUNContext`` object prior to any other calls to SUNDIALS library
functions by calling:

.. c:function:: int SUNContext_Create(void* comm, SUNContext* ctx)

   Creates a ``SUNContext`` object associated with the thread of execution.
   The data of the ``SUNContext`` class is private.

   **Arguments**:
      * ``comm`` -- a pointer to the MPI communicator or NULL if not using MPI.
      * ``ctx`` --  [in,out] upon successful exit, a pointer to the newly
        created ``SUNContext`` object.

   **Returns**:
      * Will return < 0 if an error occurs, and zero otherwise.

The created ``SUNContext`` object should be provided to the constructor routines
for different SUNDIALS classes/modules. E.g.,

.. code-block:: C

   SUNContext sunctx;
   void* package_mem;
   N_Vector x;

   SUNContext_Create(NULL, &sunctx);

   package_mem = CVodeCreate(..., sunctx);
   package_mem = IDACreate(..., sunctx);
   package_mem = KINCreate(..., sunctx);
   package_mem = ARKStepCreate(..., sunctx);

   x = N_VNew_<SomeVector>(..., sunctx);

After all other SUNDIALS code, the ``SUNContext`` object should be freed with a call to:

.. c:function:: int SUNContext_Free(SUNContext* ctx)

   Frees the ``SUNContext`` object.

   **Arguments**:
      * ``ctx`` -- pointer to a valid ``SUNContext`` object, ``NULL`` upon successful return.

   **Returns**:
      * Will return < 0 if an error occurs, and zero otherwise.


.. warning::

   When MPI is being used, the :c:func:`SUNContext_Free` must be called prior to ``MPI_Finalize``.



The ``SUNContext`` API further consists of the following functions:

.. c:function:: int SUNContext_GetProfiler(SUNContext ctx, SUNProfiler* profiler)

   Gets the ``SUNProfiler`` object associated with the ``SUNContext`` object.

   **Arguments**:
      * ``ctx`` -- a valid ``SUNContext`` object.
      * ``profiler`` -- [in,out] a pointer to the ``SUNProfiler`` object associated with this context; will be ``NULL`` if profiling is not enabled.

   **Returns**:
      * Will return < 0 if an error occurs, and zero otherwise.


.. c:function:: int SUNContext_SetProfiler(SUNContext ctx, SUNProfiler profiler)

   Sets the ``SUNProfiler`` object associated with the ``SUNContext`` object.

   **Arguments**:
      * ``ctx`` -- a valid ``SUNContext`` object.
      * ``profiler`` -- a ``SUNProfiler`` object to associate with this context; this is ignored if profiling is not enabled.

   **Returns**:
      * Will return < 0 if an error occurs, and zero otherwise.


.. _SUNDIALS.SUNContext.Threads:

Implications for Multi-Threading
--------------------------------

In multi-threading applications where multiple SUNDIALS simulations are conducted concurrently,
e.g. by having one instance of an integrator per thread, a ``SUNContext`` object needs to be created
for each thread.

.. _SUNDIALS.SUNContext.CPP:

Convenience class for C++ Users
-------------------------------

For C++ users, a class, ``sundials::Context``, that follows RAII is provided:

.. code-block:: cpp

   namespace sundials
   {

   class Context
   {
   public:
      Context(void* comm = NULL)
      {
         SUNContext_Create(comm, &sunctx_);
      }

      operator SUNContext() { return sunctx_; }

      ~Context()
      {
         SUNContext_Free(&sunctx_);
      }

   private:
      SUNContext sunctx_;

   };

   } // namespace sundials
