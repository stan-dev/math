.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   Based on ERKStep by Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.MRIStep.Skeleton:

A skeleton of the user's main program
============================================

The following is a skeleton of the user's main program (or calling program) for
the integration of an ODE IVP using the MRIStep module.  Most of the steps are
independent of the NVECTOR, SUNMATRIX, SUNLINSOL and SUNNONLINSOL
implementations used.  For the steps that are not, refer to :numref:`NVectors`,
:numref:`SUNMatrix`, :numref:`SUNLinSol`, and :numref:`SUNNonlinSol` for the
specific name of the function to be called or macro to be referenced.

.. index:: User main program

#. Initialize parallel or multi-threaded environment, if appropriate.

   For example, call ``MPI_Init`` to initialize MPI if used, or set
   ``num_threads``, the number of threads to use within the threaded
   vector functions, if used.

#. Create the SUNDIALS context object

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

#. Create an inner stepper object to solve the fast (inner) IVP

   * If using ARKStep as the fast (inner) integrator, create the ARKStep object
     with :c:func:`ARKStepCreate` and configure the integrator as desired for
     evolving the fast time scale. See sections :numref:`ARKODE.Usage.ARKStep.Skeleton`
     and :numref:`ARKODE.Usage.ARKStep.OptionalInputs` for details on configuring
     ARKStep.

     Once the ARKStep object is setup, create an ``MRIStepInnerStepper`` object
     with :c:func:`ARKStepCreateMRIStepInnerStepper`.

   * If supplying a user-defined fast (inner) integrator, create the
     ``MRIStepInnerStepper`` object as described in section
     :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper`.

   .. note::

      When using ARKStep as a fast (inner) integrator it is the user's
      responsibility to create, configure, and attach the integrator to the
      MRIStep module. User-specified options regarding how this fast integration
      should be performed (e.g., adaptive vs. fixed time step,
      explicit/implicit/ImEx partitioning, algebraic solvers, etc.) will be
      respected during evolution of the fast time scale during MRIStep
      integration.

      Due to the algorithms supported in MRIStep, the ARKStep module used for
      the fast time scale must be configured with an identity mass matrix.

      If a *user_data* pointer needs to be passed to user functions called by
      the fast (inner) integrator then it should be attached here by calling
      :c:func:`ARKStepSetUserData()`. This *user_data* pointer will only be
      passed to user-supplied functions that are attached to the fast (inner)
      integrator. To supply a *user_data* pointer to user-supplied functions
      called by the slow (outer) integrator the desired pointer should be
      attached by calling :c:func:`MRIStepSetUserData()` after creating the
      MRIStep memory below. The *user_data* pointers attached to the inner and
      outer integrators may be the same or different depending on what is
      required by the user code.

      Specifying a rootfinding problem for the fast integration is not
      supported. Rootfinding problems should be created and initialized with
      the slow integrator. See the steps below and :c:func:`MRIStepRootInit()`
      for more details.

#. Create an MRIStep object for the slow (outer) integration

   Create the MRIStep object by calling  :c:func:`MRIStepCreate`. One of the
   inputs to :c:func:`MRIStepCreate` is the ``MRIStepInnerStepper`` object for
   solving the fast (inner) IVP created in the previous step.

#. Set the slow step size

   Call :c:func:`MRIStepSetFixedStep()` to specify the slow time step
   size.

#. Create and configure implicit solvers (*as appropriate*)

   Specifically, if MRIStep is configured with an implicit slow right-hand side
   function in the prior step, then the following steps are recommended:

   #. Specify integration tolerances

      Call :c:func:`MRIStepSStolerances()` or :c:func:`MRIStepSVtolerances()` to
      specify either a scalar relative tolerance and scalar absolute tolerance,
      or a scalar relative tolerance and a vector of absolute tolerances,
      respectively.  Alternatively, call :c:func:`MRIStepWFtolerances()`
      to specify a function which sets directly the weights used in
      evaluating WRMS vector norms. See :numref:`ARKODE.Usage.MRIStep.Tolerances` for
      details.

   #. Create nonlinear solver object

      If a non-default nonlinear solver object is desired for implicit
      MRI stage solves (see :numref:`ARKODE.Usage.MRIStep.NonlinearSolvers`),
      then that nonlinear solver object must be created by using
      the appropriate functions defined by the particular SUNNONLINSOL
      implementation (e.g., ``NLS = SUNNonlinSol_***(...);`` where
      ``***`` is the name of the nonlinear solver (see
      :numref:`SUNNonlinSol` for details).

      For the SUNDIALS-supplied SUNNONLINSOL implementations, the
      nonlinear solver object may be created using a call of the form

      .. code-block:: c

         SUNNonlinearSolver NLS = SUNNonlinSol_*(...);

      where ``*`` can be replaced with "Newton", "FixedPoint", or other
      options, as discussed in the sections
      :numref:`ARKODE.Usage.ARKStep.NonlinearSolvers` and :numref:`SUNNonlinSol`.

      Note: by default, MRIStep will use the Newton nonlinear solver
      (see section :numref:`SUNNonlinSol.Newton`), so a custom nonlinear solver
      object is only needed when using a *different* solver, or for the user
      to exercise additional controls over the Newton solver.

   #. Attach nonlinear solver module

      If a nonlinear solver object was created above, then it must be
      attached to MRIStep using the call (for details see
      :numref:`ARKODE.Usage.MRIStep.NonlinearSolvers`):

      .. code-block:: c

         ier = MRIStepSetNonlinearSolver(...);

   #. Set nonlinear solver optional inputs

      Call the appropriate set functions for the selected nonlinear
      solver module to change optional inputs specific to that nonlinear
      solver.  These *must* be called after attaching the nonlinear
      solver to MRIStep, otherwise the optional inputs will be
      overridden by MRIStep defaults.  See :numref:`SUNNonlinSol` for more
      information on optional inputs.

   #. Create matrix object

      If a nonlinear solver requiring a linear solver will be used (e.g.,
      a Newton iteration) and if that linear solver will be matrix-based,
      then a template Jacobian matrix must be created by using the
      appropriate functions defined by the particular SUNMATRIX
      implementation.

      For the SUNDIALS-supplied SUNMATRIX implementations, the
      matrix object may be created using a call of the form

      .. code-block:: c

         SUNMatrix A = SUNBandMatrix(...);

      or similar for other matrix modules (see :numref:`SUNMatrix` for
      further information).

   #. Create linear solver object

      If a nonlinear solver requiring a linear solver will be used (e.g.,
      a Newton iteration), then the desired linear solver object(s) must be
      created by using the appropriate functions defined by the particular
      SUNLINSOL implementation.

      For any of the SUNDIALS-supplied SUNLINSOL implementations, the
      linear solver object may be created using a call of the form

      .. code-block:: c

         SUNLinearSolver LS = SUNLinSol_*(...);

      where ``*`` can be replaced with "Dense", "SPGMR", or other
      options, as discussed in :numref:`SUNLinSol`.

   #. Set linear solver optional inputs

      Call ``*Set*`` functions from the selected linear solver module
      to change optional inputs specific to that linear solver.  See the
      documentation for each SUNLINSOL module in :numref:`SUNLinSol` for details.

   #. Attach linear solver module

      If a linear solver was created above for implicit MRI stage solves,
      initialize the ARKLS linear solver interface by attaching the
      linear solver object (and Jacobian matrix object, if applicable)
      with the call (for details see :numref:`ARKODE.Usage.MRIStep.LinearSolvers`):

      .. code-block:: c

         ier = MRIStepSetLinearSolver(...);

#. Set optional inputs

   Call ``MRIStepSet*`` functions to change any optional inputs that
   control the behavior of MRIStep from their default values. See
   :numref:`ARKODE.Usage.MRIStep.OptionalInputs` for details.

#. Specify rootfinding problem

   Optionally, call :c:func:`MRIStepRootInit()` to initialize a rootfinding
   problem to be solved during the integration of the ODE system. See
   :numref:`ARKODE.Usage.MRIStep.RootFinding` for general details, and
   :numref:`ARKODE.Usage.MRIStep.OptionalInputs` for relevant optional input calls.

#. Advance solution in time

   For each point at which output is desired, call

   .. code-block:: c

      ier = MRIStepEvolve(arkode_mem, tout, yout, &tret, itask);

   Here, ``itask`` specifies the return mode. The vector ``yout``
   (which can be the same as the vector ``y0`` above) will contain
   :math:`y(t_\text{out})`. See :numref:`ARKODE.Usage.MRIStep.Integration` for details.

#. Get optional outputs

   Call ``MRIStepGet*`` and/or ``ARKStepGet*`` functions to obtain optional
   output from the slow or fast integrators respectively. See
   :numref:`ARKODE.Usage.MRIStep.OptionalOutputs` and
   :numref:`ARKODE.Usage.ARKStep.OptionalOutputs` for details.

#. Deallocate memory for solution vector

   Upon completion of the integration, deallocate memory for the
   vector ``y`` (or ``yout``) by calling the NVECTOR destructor
   function:

   .. code-block:: c

      N_VDestroy(y);

#. Free solver memory

   * If ARKStep was used as the fast (inner) IVP integrator, call
     :c:func:`MRIStepInnerStepper_Free` and :c:func:`ARKStepFree` to free the
     memory allocated for the fast (inner) integrator.

   * If a user-defined fast (inner) integrator was supplied, free the integrator
     content and call :c:func:`MRIStepInnerStepper_Free` to free the
     ``MRIStepInnerStepper`` object.

   * Call :c:func:`MRIStepFree` to free the memory allocated for the slow
     integration object.

#. Free linear solver and matrix memory (*as appropriate*)

    Call :c:func:`SUNLinSolFree()` and (possibly)
    :c:func:`SUNMatDestroy()` to free any memory allocated for any
    linear solver and/or matrix objects created above for either the fast or
    slow integrators.

#. Free nonlinear solver memory (*as appropriate*)

   If a user-supplied ``SUNNonlinearSolver`` was provided to MRIStep,
   then call :c:func:`SUNNonlinSolFree()` to free any memory allocated
   for the nonlinear solver object created above.

#. **Free the SUNContext object**
   Call :c:func:`SUNContext_Free` to free the memory allocated for the ``SUNContext`` object.

 #. Finalize MPI, if used

    Call ``MPI_Finalize`` to terminate MPI.
