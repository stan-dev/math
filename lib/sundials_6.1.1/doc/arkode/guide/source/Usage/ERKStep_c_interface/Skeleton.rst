.. ----------------------------------------------------------------
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

.. _ARKODE.Usage.ERKStep.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling
program) for the integration of an ODE IVP using the ERKStep module.
Most of the steps are independent of the NVECTOR implementation used.
For the steps that are not, refer to :numref:`NVectors` for
the specific name of the function to be called or macro to be
referenced.

.. index:: User main program

#. Initialize parallel or multi-threaded environment, if appropriate.

   For example, call ``MPI_Init`` to initialize MPI if used, or set
   ``num_threads``, the number of threads to use within the threaded
   vector functions, if used.

#. Create the SUNDIALS simulation context object.

   Call :c:func:`SUNContext_Create` to allocate the ``SUNContext`` object.

#. Set problem dimensions, etc.

   This generally includes the problem size, ``N``, and may include
   the local vector length ``Nlocal``.

   .. note::

      The variables ``N`` and ``Nlocal`` should be of type
      ``sunindextype``.

#. Set vector of initial values

   To set the vector ``y0`` of initial values, use the appropriate
   functions defined by the particular NVECTOR implementation.

   For native SUNDIALS vector implementations (except the CUDA and
   RAJA based ones), use a call of the form

   .. code-block:: c

      y0 = N_VMake_***(..., ydata);

   if the ``realtype`` array ``ydata`` containing the initial values of
   :math:`y` already exists.  Otherwise, create a new vector by making
   a call of the form

   .. code-block:: c

      y0 = N_VNew_***(...);

   and then set its elements by accessing the underlying data where it
   is located with a call of the form

   .. code-block:: c

      ydata = N_VGetArrayPointer_***(y0);

   For details on each of SUNDIALS' provided vector implementations, see
   the corresponding sections in :numref:`NVectors` for details.

#. Create ERKStep object

   Call ``arkode_mem = ERKStepCreate(...)`` to create the ERKStep memory
   block. :c:func:`ERKStepCreate` returns a ``void*`` pointer to
   this memory structure. See :numref:`ARKODE.Usage.ERKStep.Initialization` for
   details.

#. Specify integration tolerances

   Call :c:func:`ERKStepSStolerances()` or
   :c:func:`ERKStepSVtolerances()` to specify either a scalar relative
   tolerance and scalar absolute tolerance, or a scalar relative
   tolerance and a vector of absolute tolerances,
   respectively.  Alternatively, call :c:func:`ERKStepWFtolerances()`
   to specify a function which sets directly the weights used in
   evaluating WRMS vector norms. See :numref:`ARKODE.Usage.ERKStep.Tolerances`
   for details.

#. Set optional inputs

   Call ``ERKStepSet*`` functions to change any optional inputs that
   control the behavior of ERKStep from their default values. See
   :numref:`ARKODE.Usage.ERKStep.OptionalInputs` for details.

#. Specify rootfinding problem

   Optionally, call :c:func:`ERKStepRootInit()` to initialize a rootfinding
   problem to be solved during the integration of the ODE system. See
   :numref:`ARKODE.Usage.ERKStep.RootFinding` for general details, and
   :numref:`ARKODE.Usage.ERKStep.OptionalInputs` for relevant optional
   input calls.

#. Advance solution in time

   For each point at which output is desired, call

   .. code-block:: c

      ier = ERKStepEvolve(arkode_mem, tout, yout, &tret, itask);

   Here, ``itask`` specifies the return mode. The vector ``yout``
   (which can be the same as the vector ``y0`` above) will contain
   :math:`y(t_\text{out})`. See :numref:`ARKODE.Usage.ERKStep.Integration`
   for details.

#. Get optional outputs

   Call ``ERKStepGet*`` functions to obtain optional output. See
   :numref:`ARKODE.Usage.ERKStep.OptionalOutputs` for details.

#. Deallocate memory for solution vector

    Upon completion of the integration, deallocate memory for the
    vector ``y`` (or ``yout``) by calling the NVECTOR destructor
    function:

    .. code-block:: c

       N_VDestroy(y);

#. Free solver memory

    Call :c:func:`ERKStepFree()` to free the memory allocated for
    the ERKStep module.

#. Finalize MPI, if used

    Call ``MPI_Finalize`` to terminate MPI.
