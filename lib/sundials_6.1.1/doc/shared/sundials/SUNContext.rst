.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDIALS.SUNContext:

The SUNContext Type
=====================

In SUNDIALS v6.0.0, the concept of a SUNDIALS simulation context was introduced,
in particular the ``SUNContext`` class. All of the SUNDIALS objects (vectors,
linear and nonlinear solvers, matrices, etc) that collectively form a SUNDIALS
simulation, hold a reference to a common ``SUNContext`` object.

The ``SUNContext`` class/type is defined in the header file
``sundials/sundials_context.h`` as

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

Implications for task-based programming and multi-threading
------------------------------------------------------------

Applications that need to have *concurrently initialized* SUNDIALS simulations
need to take care to understand the following:

#. A ``SUNContext`` object must only be associated with *one* SUNDIALS simulation
(a solver object and its associated vectors etc.) at a time.

   - Concurrently initialized is not the same as concurrently executing. Even if
     two SUNDIALS simulations execute sequentially, if both are initialized
     at the same time with the same ``SUNContext``, behavior is undefined.

   - It is OK to reuse a ``SUNContext`` object with another SUNDIALS simulation
     after the first simulation has completed and all of the simulation's
     associated objects (vectors, matrices, algebraic solvers, etc.) have been
     destroyed.

#. The creation and destruction of a ``SUNContext`` object is cheap, especially
in comparison to the cost of creating/destroying a SUNDIALS solver object.

The following (incomplete) code examples demonstrate these points using CVODE as
the example SUNDIALS package.

.. code-block:: c

   SUNContext sunctxs[num_threads];
   int cvode_initialized[num_threads];
   void* cvode_mem[num_threads];

   // Create
   for (int i = 0; i < num_threads; i++) {
      sunctxs[i] = SUNContext_Create(...);
      cvode_mem[i] = CVodeCreate(..., sunctxs[i]);
      cvode_initialized[i] = 0; // not yet initialized
      // set optional cvode inputs...
   }

   // Solve
   #pragma omp parallel for
   for (int i = 0; i < num_problems; i++) {
      int retval = 0;
      int tid = omp_get_thread_num();
      if (!cvode_initialized[tid]) {
         retval = CVodeInit(cvode_mem[tid], ...);
         cvode_initialized[tid] = 1;
      } else {
         retval = CVodeReInit(cvode_mem[tid], ...);
      }
      CVode(cvode_mem[i], ...);
   }

   // Destroy
   for (int i = 0; i < num_threads; i++) {
      // get optional cvode outputs...
      CVodeFree(&cvode_mem[i]);
      SUNContext_Free(&sunctxs[i]);
   }

Since each thread has its own unique CVODE and SUNContext object pair, there
should be no thread-safety issues. Users should be sure that you apply the same
idea to the other SUNDIALS objects needed as well (e.g. an ``N_Vector``).

The variation of the above code example demonstrates another possible approach:

.. code-block:: c

   // Create, Solve, Destroy
   #pragma omp parallel for
   for (int i = 0; i < num_problems; i++) {
      int retval = 0;
      void* cvode_mem;
      SUNContext sunctx;

      sunctx = SUNContext_Create(...);
      cvode_mem = CVodeCreate(..., sunctx);
      retval = CVodeInit(cvode_mem, ...);

      // set optional cvode inputs...

      CVode(cvode_mem, ...);

      // get optional cvode outputs...

      CVodeFree(&cvode_mem);
      SUNContext_Free(&sunctx);
   }

So long as the overhead of creating/destroying the CVODE object is small
compared to the cost of solving the ODE, this approach is a fine alternative to
the first approach since :c:func:`SUNContext_Create` and
:c:func:`SUNContext_Free` are much cheaper than the CVODE create/free routines.


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
