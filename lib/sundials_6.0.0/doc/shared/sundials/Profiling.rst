.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDIALS.Profiling:

Performance Profiling
=====================

SUNDIALS includes a lightweight performance profiling layer that can be enabled
at compile-time. Optionally, this profiling layer can leverage Caliper
:cite:p:`Caliper:2016` for more advanced instrumentation and profiling. By
default, only SUNDIALS library code is profiled. However, a public profiling API
can be utilized to leverage the SUNDIALS profiler to time user code regions as
well (see :numref:`SUNDIALS.Profiling.API`).

.. _SUNDIALS.Profiling.Enabling:

Enabling Profiling
------------------

To enable profiling, SUNDIALS must be built with the CMake option
:cmakeop:`SUNDIALS_BUILD_WITH_PROFILING` set to ``ON``. To utilize Caliper
support, the CMake option :cmakeop:`ENABLE_CALIPER` must also be set to ``ON``.
More details in regards to configuring SUNDIALS with CMake can be found in
:numref:`Installation`.

When SUNDIALS is built with profiling enabled and **without Caliper**, then the
environment variable ``SUNPROFILER_PRINT`` can be utilized to enable/disable the
printing of profiler information. Setting ``SUNPROFILER_PRINT=1`` will cause the
profiling information to be printed to stdout when the SUNDIALS simulation context is
freed. Setting ``SUNPROFILER_PRINT=0`` will result in no profiling information
being printed unless the :c:func:`SUNProfiler_Print` function is called
explicitly. By default, ``SUNPROFILER_PRINT`` is assumed to be ``0``.
``SUNPROFILER_PRINT`` can also be set to a file path where the output should be printed.

If Caliper is enabled, then users should refer to the `Caliper documentation <https://software.llnl.gov/Caliper/>`_
for information on getting profiler output. In most cases, this involves
setting the ``CALI_CONFIG`` environment variable.

.. warning::

   While the SUNDIALS profiling scheme is relatively lightweight, enabling
   profiling can still negatively impact performance. As such, it is recommended
   that profiling is enabled judiciously.


.. _SUNDIALS.Profiling.API:

Profiler API
------------

The primary way of interacting with the SUNDIALS profiler is through the following
macros:

.. code-block:: C

   SUNDIALS_MARK_FUNCTION_BEGIN(profobj)
   SUNDIALS_MARK_FUNCTION_END(profobj)
   SUNDIALS_WRAP_STATEMENT(profobj, name, stmt)
   SUNDIALS_MARK_BEGIN(profobj, name)
   SUNDIALS_MARK_END(profobj, name)

Additionally, in C++ applications, the follow macro is available:

.. code-block:: C++

   SUNDIALS_CXX_MARK_FUNCTION(profobj)

These macros can be used to time specific functions or code regions. When using
the ``*_BEGIN`` macros, it is important that a matching ``*_END`` macro is
placed at all exit points for the scope/function. The
``SUNDIALS_CXX_MARK_FUNCTION`` macro only needs to be placed at the beginning of
a function, and leverages RAII to implicitly end the region.

The ``profobj`` argument to the macro should be a ``SUNProfiler`` object, i.e.
an instance of the struct

.. c:type:: struct _SUNProfiler *SUNProfiler

When SUNDIALS is built with profiling, a default profiling object is stored in the
``SUNContext`` object and can be accessed with a call to
:c:func:`SUNContext_GetProfiler`.

The ``name`` argument should be a unique string indicating the name of the
region/function. It is important that the name given to the ``*_BEGIN`` macros
matches the name given to the ``*_END`` macros.


In addition to the macros, the following methods of the ``SUNProfiler`` class
are available.


.. c:function:: int SUNProfiler_Create(void* comm, const char* title, SUNProfiler* p)

   Creates a new ``SUNProfiler`` object.

   **Arguments:**
      * ``comm`` -- a pointer to the MPI communicator if MPI is enabled, otherwise can be ``NULL``
      * ``title`` -- a title or description of the profiler
      * ``p`` -- [in,out] On input this is a pointer to a ``SUNProfiler``, on output it will point to a new ``SUNProfiler`` instance

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred


.. c:function:: int SUNProfiler_Free(SUNProfiler* p)

   Frees a ``SUNProfiler`` object.

   **Arguments:**
      * ``p`` -- [in,out] On input this is a pointer to a ``SUNProfiler``, on output it will be ``NULL``

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred


.. c:function:: int SUNProfiler_Begin(SUNProfiler p, const char* name)

   Starts timing the region indicated by the ``name``.

   **Arguments:**
      * ``p`` -- a ``SUNProfiler`` object
      * ``name`` -- a name for the profiling region

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred


.. c:function:: int SUNProfiler_End(SUNProfiler p, const char* name)

   Ends the timing of a region indicated by the ``name``.

   **Arguments:**
      * ``p`` -- a ``SUNProfiler`` object
      * ``name`` -- a name for the profiling region

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred


.. c:function:: int SUNProfiler_Print(SUNProfiler p)

   Prints out a profiling summary. When constructed with an MPI comm the summary
   will include the average and maximum time per rank (in seconds) spent in each
   marked up region.

   **Arguments:**
      * ``p`` -- a ``SUNProfiler`` object

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred


.. _SUNDIALS.Profiling.Example:

Example Usage
-------------

The following is an excerpt from the CVODE example code ``examples/cvode/serial/cvAdvDiff_bnd.c``.
It is applicable to any of the SUNDIALS solver packages.

.. code-block:: c

   SUNContext ctx;
   SUNProfiler profobj;

   /* Create the SUNDIALS context */
   retval = SUNContext_Create(NULL, &ctx);

   /* Get a reference to the profiler */
   retval = SUNContext_GetProfiler(ctx, &profobj);

   /* ... */

   SUNDIALS_MARK_BEGIN(profobj, "Integration loop");
   umax = N_VMaxNorm(u);
   PrintHeader(reltol, abstol, umax);
   for(iout=1, tout=T1; iout <= NOUT; iout++, tout += DTOUT) {
      retval = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
      umax = N_VMaxNorm(u);
      retval = CVodeGetNumSteps(cvode_mem, &nst);
      PrintOutput(t, umax, nst);
   }
   SUNDIALS_MARK_END(profobj, "Integration loop");
   PrintFinalStats(cvode_mem);  /* Print some final statistics   */


.. _SUNDIALS.Profiling.Other:

Other Considerations
--------------------

If many regions are being timed, it may be necessary to increase the maximum
number of profiler entries (the default is ``2560``). This can be done
by setting the environment variable ``SUNPROFILER_MAX_ENTRIES``.
