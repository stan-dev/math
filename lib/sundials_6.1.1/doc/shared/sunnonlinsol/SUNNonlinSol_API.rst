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

.. _SUNNonlinSol.API:

===============================
The SUNNonlinearSolver API
===============================

The SUNNonlinSol API defines several nonlinear solver operations that enable
SUNDIALS integrators to utilize any SUNNonlinSol implementation that
provides the required functions. These functions can be divided into three
categories. The first are the core nonlinear solver functions. The second
consists of "set" routines to supply the nonlinear solver with
functions provided by the SUNDIALS time integrators and to modify solver
parameters. The final group consists of "get" routines for retrieving nonlinear
solver statistics. All of these functions are defined in the header file
``sundials/sundials_nonlinearsolver.h``.



.. _SUNNonlinSol.API.CoreFn:

SUNNonlinearSolver core functions
-----------------------------------------------------

The core nonlinear solver functions consist of two required functions to get the
nonlinear solver type (:c:func:`SUNNonlinsSolGetType`) and solve the nonlinear system
(:c:func:`SUNNonlinSolSolve`). The remaining three functions for nonlinear solver
initialization (:c:func:`SUNNonlinSolInitialization`), setup
(:c:func:`SUNNonlinSolSetup`), and destruction (:c:func:`SUNNonlinSolFree`) are optional.


.. c:function:: SUNNonlinearSolver_Type SUNNonlinSolGetType(SUNNonlinearSolver NLS)

   This *required* function returns the nonlinear solver type.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.

   **Return value:**
      The SUNNonlinSol type identifier (of type ``int``) will be one
      of the following:

      * ``SUNNONLINEARSOLVER_ROOTFIND`` -- ``0``, the SUNNonlinSol module
        solves :math:`F(y) = 0`.

      * ``SUNNONLINEARSOLVER_FIXEDPOINT`` -- ``1``, the SUNNonlinSol
        module solves :math:`G(y) = y`.



.. c:function:: int SUNNonlinSolInitialize(SUNNonlinearSolver NLS)

   This *optional* function handles nonlinear solver initialization
   and may perform any necessary memory allocations.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.

   **Return value:**
      The return value is zero for a successful call and a
      negative value for a failure.

   **Notes:**
      It is assumed all solver-specific options have been set
      prior to calling :c:func:`SUNNonlinSolInitialize`. SUNNonlinSol
      implementations that do not require initialization may set this
      operation to ``NULL``.


.. c:function:: int SUNNonlinSolSetup(SUNNonlinearSolver NLS, N_Vector y, void* mem)

   This *optional* function performs any solver setup needed for a nonlinear solve.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *y* -- the initial guess passed to the nonlinear solver.
      * *mem* -- the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successful call and a
      negative value for a failure.

   **Notes:**
      SUNDIALS integrators call :c:func:`SUNonlinSolSetup` before
      each step attempt. SUNNonlinSol implementations that do not
      require setup may set this operation to ``NULL``.


.. c:function:: int SUNNonlinSolSolve(SUNNonlinearSolver NLS, N_Vector y0, N_Vector ycor, N_Vector w, realtype tol, booleantype callLSetup, void *mem)

   This *required* function solves the nonlinear system
   :math:`F(y)=0` or :math:`G(y)=y`.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *y0* -- the predicted value for the new solution state. This
        *must* remain unchanged throughout the solution process.
      * *ycor* -- on input the initial guess for the correction to the predicted
        state (zero) and on output the final correction to the predicted
        state.
      * *w* -- the solution error weight vector used for computing weighted error norms.
      * *tol* -- the requested solution tolerance in the weighted root-mean-squared norm.
      * *callLSetup* -- a flag indicating that the integrator
        recommends for the linear solver setup function to be called.
      * *mem* -- the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successul solve, a positive value
      for a recoverable error (i.e., the solve failed and the integrator
      should reduce the step size and reattempt the step), and a negative
      value for an unrecoverable error (i.e., the solve failed the and
      the integrator should halt and return an error to the user).


.. c:function:: int SUNNonlinSolFree(SUNNonlinearSolver NLS)

   This *optional* function frees any memory allocated by the
   nonlinear solver.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure. SUNNonlinSol implementations that
      do not allocate data may set this operation to ``NULL``.




.. _SUNNonlinSol.API.SetFn:

SUNNonlinearSolver "set" functions
-------------------------------------

The following functions are used to supply nonlinear solver modules with
functions defined by the SUNDIALS integrators and to modify solver
parameters. Only the routine for setting the nonlinear system defining function
(:c:func:`SUNNonlinSolSetSysFn`) is required. All other set functions are optional.


.. c:function:: int SUNNonlinSolSetSysFn(SUNNonlinearSolver NLS, SUNNonlinSolSysFn SysFn)

   This *required* function is used to provide the nonlinear solver
   with the function defining the nonlinear system. This is the function
   :math:`F(y)` in :math:`F(y)=0` for ``SUNNONLINEARSOLVER_ROOTFIND`` modules or
   :math:`G(y)` in :math:`G(y)=y` for ``SUNNONLINEARSOLVER_FIXEDPOINT`` modules.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *SysFn* -- the function defining the nonlinear system. See
        :numref:`SUNNonlinSol.API.SUNSuppliedFn` for the definition of
        :c:type:`SUNNonlinSolSysFn`.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure.


.. c:function:: int SUNNonlinSolSetLSetupFn(SUNNonlinearSolver NLS, SUNNonlinSolLSetupFn SetupFn)

   This *optional* function is called by SUNDIALS integrators to provide
   the nonlinear solver with access to its linear solver setup function.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *SetupFn* -- a wrapper function to the SUNDIALS integrator's linear solver setup
        function. See :numref:`SUNNonlinSol.API.SUNSuppliedFn`  for the
        definition of :c:type:`SUNNonlinSolLSetupFn`.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure.

   **Notes:**
      The :c:type:`SUNNonlinSolLSetupFn` function sets up the
      linear system :math:`Ax=b` where :math:`A = \frac{\partial
      F}{\partial y}` is the linearization of the nonlinear residual
      function :math:`F(y) = 0` (when using SUNLinSol direct linear
      solvers) or calls the user-defined preconditioner setup function
      (when using SUNLinSol iterative linear solvers). SUNNonlinSol
      implementations that do not require solving this system, do not
      utilize SUNLinSol linear solvers, or use SUNLinSol linear solvers
      that do not require setup may set this operation to ``NULL``.



.. c:function:: int SUNNonlinSolSetLSolveFn(SUNNonlinearSolver NLS, SUNNonlinSolLSolveFn SolveFn)

   This *optional* function is called by SUNDIALS integrators to provide
   the nonlinear solver with access to its linear solver solve function.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *SolveFn* -- a wrapper function to the SUNDIALS integrator's
        linear solver solve function. See
        :numref:`SUNNonlinSol.API.SUNSuppliedFn` for the definition of
        :c:type:`SUNNonlinSolLSolveFn`.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure.

   **Notes:**
      The :c:type:`SUNNonlinSolLSolveFn` function solves the
      linear system :math:`Ax=b` where :math:`A = \frac{\partial
      F}{\partial y}` is the linearization of the nonlinear residual
      function :math:`F(y) = 0`.  SUNNonlinSol implementations that do
      not require solving this system or do not use SUNLinSol linear
      solvers may set this operation to ``NULL``.



.. c:function:: int SUNNonlinSolSetConvTestFn(SUNNonlinearSolver NLS, SUNNonlinSolConvTestFn CTestFn, void* ctest_data)

   This *optional* function is used to provide the nonlinear solver
   with a function for determining if the nonlinear solver iteration
   has converged. This is typically called by SUNDIALS integrators to
   define their nonlinear convergence criteria, but may be replaced by
   the user.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *CTestFn* -- a SUNDIALS integrator's nonlinear solver
        convergence test function. See
        :numref:`SUNNonlinSol.API.SUNSuppliedFn` for the definition of
        :c:type:`SUNNonlinSolConvTestFn`.
      * *ctest_data* -- is a data pointer passed to *CTestFn* every time it is
        called.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure.

   **Notes:**
      SUNNonlinSol implementations utilizing their own convergence test
      criteria may set this function to ``NULL``.



.. c:function:: int SUNNonlinSolSetMaxIters(SUNNonlinearSolver NLS, int maxiters)

   This *optional* function sets the maximum number of nonlinear solver
   iterations. This is typically called by SUNDIALS integrators to
   define their default iteration limit, but may be adjusted by the user.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *maxiters* -- the maximum number of nonlinear iterations.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure (e.g., :math:`maxiters < 1`).




.. _SUNNonlinSol.API.GetFn:

SUNNonlinearSolver "get" functions
----------------------------------

The following functions allow SUNDIALS integrators to retrieve nonlinear
solver statistics. The routines to get the number of iterations in the most
recent solve (:c:func:`SUNNonlinSolGetNumIters`) and number of convergence failures
are optional. The routine to get the current nonlinear solver iteration
(:c:func:`SUNNonlinSolGetCurIter`) is required when using the convergence test
provided by the SUNDIALS integrator or when using an iterative SUNLinSol
linear solver module; otherwise :c:func:`SUNNonlinSolGetCurIter` is optional.


.. c:function:: int SUNNonlinSolGetNumIters(SUNNonlinearSolver NLS, long int *niters)

   This *optional* function returns the number of nonlinear solver iterations
   in the most recent solve. This is typically called by the SUNDIALS
   integrator to store the nonlinear solver statistics, but may also be
   called by the user.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *niters* -- the total number of nonlinear solver iterations.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure.


.. c:function:: int SUNNonlinSolGetCurIter(SUNNonlinearSolver NLS, int *iter)

   This function returns the iteration index of the current nonlinear
   solve. This function is *required* when using SUNDIALS
   integrator-provided convergence tests or when using an iterative
   SUNLinSol linear solver module; otherwise it is *optional*.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *iter* -- the nonlinear solver iteration in the current solve
        starting from zero.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure.


.. c:function:: int SUNNonlinSolGetNumConvFails(SUNNonlinearSolver NLS, long int *nconvfails)

   This *optional* function returns the number of nonlinear solver convergence
   failures in the most recent solve. This is typically called by the SUNDIALS
   integrator to store the nonlinear solver statistics, but may also be called
   by the user.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *nconvfails* -- the total number of nonlinear solver convergence failures.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure.


.. _SUNNonlinSol.API.SUNSuppliedFn:

Functions provided by SUNDIALS integrators
--------------------------------------------

To interface with SUNNonlinSol modules, the SUNDIALS integrators
supply a variety of routines for evaluating the nonlinear system,
calling the SUNLinSol setup and solve functions, and testing the
nonlinear iteration for convergence.  These integrator-provided routines
translate between the user-supplied ODE or DAE systems and the generic
interfaces to the nonlinear or linear systems of equations that result
in their solution. The functions provided to a SUNNonlinSol
module have types defined in the header file
``sundials/sundials_nonlinearsolver.h``; these are also described below.


.. c:type:: int (*SUNNonlinSolSysFn)(N_Vector ycor, N_Vector F, void* mem)

   These functions evaluate the nonlinear system :math:`F(y)`
   for ``SUNNONLINEARSOLVER_ROOTFIND`` type modules or :math:`G(y)`
   for ``SUNNONLINEARSOLVER_FIXEDPOINT`` type modules. Memory
   for *F* must by be allocated prior to calling this function. The
   vector *ycor* will be left unchanged.

   **Arguments:**
      * *ycor* -- is the current correction to the predicted state at which the
        nonlinear system should be evaluated.
      * *F* -- is the output vector containing :math:`F(y)` or
        :math:`G(y)`, depending on the solver type.
      * *mem* -- is the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successul solve, a positive value for
      a recoverable error, and a negative value for an unrecoverable error.

   **Notes:**
      SUNDIALS integrators formulate nonlinear systems as a function of the
      correction to the predicted solution. On each call to the nonlinear system
      function the integrator will compute and store the current solution based on
      the input correction. Additionally, the residual will store the value of the
      ODE right-hand side function or DAE residual used in computing the nonlinear
      system. These stored values are then directly used in the integrator-supplied
      linear solver setup and solve functions as applicable.


.. c:type:: int (*SUNNonlinSolLSetupFn)(booleantype jbad, booleantype* jcur, void* mem)

   These functions are wrappers to the SUNDIALS integrator's function
   for setting up linear solves with SUNLinSol modules.

   **Arguments:**
      * *jbad* -- is an input indicating whether the nonlinear solver
        believes that :math:`A` has gone stale (``SUNTRUE``) or not (``SUNFALSE``).
      * *jcur* -- is an output indicating whether the routine has updated the
        Jacobian :math:`A` (``SUNTRUE``) or not (``SUNFALSE``).
      * *mem* -- is the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successul solve, a positive value for
      a recoverable error, and a negative value for an unrecoverable error.

   **Notes:**
      The :c:type:`SUNNonlinSolLSetupFn` function sets up the linear
      system :math:`Ax=b` where :math:`A = \frac{\partial F}{\partial y}`
      is the linearization of the nonlinear residual function
      :math:`F(y) = 0` (when using SUNLinSol direct linear solvers) or
      calls the user-defined preconditioner setup function (when using
      SUNLinSol iterative linear solvers). SUNNonlinSol implementations
      that do not require solving this system, do not utilize SUNLinSol
      linear solvers, or use SUNLinSol linear solvers that do not
      require setup may ignore these functions.

      As discussed in the description of :c:type:`SUNNonlinSolSysFn`, the linear
      solver setup function assumes that the nonlinear system function has been
      called prior to the linear solver setup function as the setup will utilize
      saved values from the nonlinear system evaluation (e.g., the updated
      solution).


.. c:type:: int (*SUNNonlinSolLSolveFn)(N_Vector b, void* mem)

   These functions are wrappers to the SUNDIALS integrator's function
   for solving linear systems with SUNLinSol modules.

   **Arguments:**
      * *b* -- contains the right-hand side vector for the linear
        solve on input and the solution to the linear system on output.
      * *mem* -- is the SUNDIALS integrator memory structure.

   **Return value:**
      The return value is zero for a successul solve, a positive value for
      a recoverable error, and a negative value for an unrecoverable error.

   **Notes:**
      The :c:type:`SUNNonlinSolLSolveFn` function solves the linear
      system :math:`Ax=b` where :math:`A = \frac{\partial F}{\partial y}`
      is the linearization of the nonlinear residual function
      :math:`F(y) = 0`. SUNNonlinSol implementations that do not
      require solving this system or do not use SUNLinSol linear solvers
      may ignore these functions.

      As discussed in the description of :c:type:`SUNNonlinSolSysFn`, the linear
      solver solve function assumes that the nonlinear system function has been
      called prior to the linear solver solve function as the setup may utilize
      saved values from the nonlinear system evaluation (e.g., the updated
      solution).


.. c:type:: int (*SUNNonlinSolConvTestFn)(SUNNonlinearSolver NLS, N_Vector ycor, N_Vector del, realtype tol, N_Vector ewt, void* ctest_data)

   These functions are SUNDIALS integrator-specific convergence tests for
   nonlinear solvers and are typically supplied by each SUNDIALS integrator,
   but users may supply custom problem-specific versions as desired.

   **Arguments:**
      * *NLS* -- is the SUNNonlinSol object.
      * *ycor* -- is the current correction (nonlinear iterate).
      * *del* -- is the difference between the current and prior nonlinear iterates.
      * *tol* -- is the nonlinear solver tolerance.
      * *ewt* -- is the weight vector used in computing weighted norms.
      * *ctest_data* -- is the data pointer provided to
        :c:func:`SUNNonlinSolSetConvTestFn()`.

   **Return value:**
      The return value of this routine will be a negative value if an
      unrecoverable error occurred or one of the following:

      * ``SUN_NLS_SUCCESS`` -- the iteration is converged.

      * ``SUN_NLS_CONTINUE`` -- the iteration has not converged, keep
        iterating.

      * ``SUN_NLS_CONV_RECVR`` -- the iteration appears to be
        diverging, try to recover.

   **Notes:**
      The tolerance passed to this routine by SUNDIALS integrators is
      the tolerance in a weighted root-mean-squared norm with error
      weight vector ``ewt``.  SUNNonlinSol modules utilizing their
      own convergence criteria may ignore these functions.



.. _SUNNonlinSol.API.ReturnCodes:

SUNNonlinearSolver return codes
---------------------------------

The functions provided to SUNNonlinSol modules by each SUNDIALS
integrator, and functions within the SUNDIALS-provided SUNNonlinSol
implementations, utilize a common set of return codes shown in
:numref:`SUNNonlinSol.API.CodeTable`.  Here, negative values correspond to non-recoverable
failures, positive values to recoverable failures, and zero to a
successful call.

.. _SUNNonlinSol.API.CodeTable:
.. table:: Description of the ``SUNNonlinearSolver`` return codes.
   :align: center

   +-----------------------+---------+---------------------------------------------------------------+
   | Name                  | Value   | Description                                                   |
   +=======================+=========+===============================================================+
   | SUN_NLS_SUCCESS       |    0    | successful call or converged solve                            |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_CONTINUE      |  901    | the nonlinear solver is not converged, keep iterating         |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_CONV_RECVR    |  902    | the nonlinear solver appears to be diverging, try to recover  |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_MEM_NULL      | -901    | a memory argument is ``NULL``                                 |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_MEM_FAIL      | -902    | a memory access or allocation failed                          |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_ILL_INPUT     | -903    | an illegal input option was provided                          |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_VECTOROP_ERR  | -904    | a NVECTOR operation failed                                    |
   +-----------------------+---------+---------------------------------------------------------------+
   | SUN_NLS_EXT_FAIL      | -905    | an external library call returned an error                    |
   +-----------------------+---------+---------------------------------------------------------------+



.. _SUNNonlinSol.API.Generic:

The generic SUNNonlinearSolver module
-----------------------------------------

SUNDIALS integrators interact with specific SUNNonlinSol
implementations through the generic SUNNonlinSol module on which all
other SUNNonlinSol implementations are built. The
``SUNNonlinearSolver`` type is a pointer to a structure containing an
implementation-dependent *content* field and an *ops*
field. The type ``SUNNonlinearSolver`` is defined as follows:

.. c:type:: struct _generic_SUNNonlinearSolver *SUNNonlinearSolver

and the generic structure is defined as

.. code-block:: c

   struct _generic_SUNNonlinearSolver {
     void *content;
     struct _generic_SUNNonlinearSolver_Ops *ops;
   };

where the ``_generic_SUNNonlinearSolver_Ops`` structure is a list of
pointers to the various actual nonlinear solver operations provided by a
specific implementation. The ``_generic_SUNNonlinearSolver_Ops``
structure is defined as

.. code-block:: c

   struct _generic_SUNNonlinearSolver_Ops {
     SUNNonlinearSolver_Type (*gettype)(SUNNonlinearSolver);
     int                     (*initialize)(SUNNonlinearSolver);
     int                     (*setup)(SUNNonlinearSolver, N_Vector, void*);
     int                     (*solve)(SUNNonlinearSolver, N_Vector, N_Vector,
                                      N_Vector, realtype, booleantype, void*);
     int                     (*free)(SUNNonlinearSolver);
     int                     (*setsysfn)(SUNNonlinearSolver, SUNNonlinSolSysFn);
     int                     (*setlsetupfn)(SUNNonlinearSolver, SUNNonlinSolLSetupFn);
     int                     (*setlsolvefn)(SUNNonlinearSolver, SUNNonlinSolLSolveFn);
     int                     (*setctestfn)(SUNNonlinearSolver, SUNNonlinSolConvTestFn,
                                           void*);
     int                     (*setmaxiters)(SUNNonlinearSolver, int);
     int                     (*getnumiters)(SUNNonlinearSolver, long int*);
     int                     (*getcuriter)(SUNNonlinearSolver, int*);
     int                     (*getnumconvfails)(SUNNonlinearSolver, long int*);
   };

The generic SUNNonlinSol module defines and implements the nonlinear
solver operations defined in
:numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`.
These routines are in fact only wrappers to the nonlinear solver
operations provided by a particular SUNNonlinSol implementation,
which are accessed through the ops field of the ``SUNNonlinearSolver``
structure. To illustrate this point we show below the implementation
of a typical nonlinear solver operation from the generic SUNNonlinSol
module, namely :c:func:`SUNNonlinSolSolve`, which solves the nonlinear
system and returns a flag denoting a successful or failed solve:

.. code-block:: c

   int SUNNonlinSolSolve(SUNNonlinearSolver NLS,
                         N_Vector y0, N_Vector y,
                         N_Vector w, realtype tol,
                         booleantype callLSetup, void* mem)
   {
     return((int) NLS->ops->solve(NLS, y0, y, w, tol, callLSetup, mem));
   }



.. _SUNNonlinSol.API.Custom:

Implementing a Custom SUNNonlinearSolver Module
--------------------------------------------------

A SUNNonlinSol implementation *must* do the following:

* Specify the content of the SUNNonlinSol module.

* Define and implement the required nonlinear solver operations defined
  in :numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`.
  Note that the names of the module routines should be unique to that
  implementation in order to permit using more than one SUNNonlinSol
  module (each with different ``SUNNonlinearSolver`` internal data
  representations) in the same code.

* Define and implement a user-callable constructor to create a
  ``SUNNonlinearSolver`` object.

To aid in the creation of custom ``SUNNonlinearSolver`` modules, the generic
``SUNNonlinearSolver`` module provides the utility functions
:c:func:`SUNNonlinSolNewEmpty` and :c:func:`SUNNonlinsolFreeEmpty`. When used
in custom ``SUNNonlinearSolver`` constructors these functions will ease the
introduction of any new optional nonlinear solver operations to the
``SUNNonlinearSolver`` API by ensuring that only required operations need to
be set.

.. c:function:: SUNNonlinearSolver SUNNonlinSolNewEmpty()

  This function allocates a new generic ``SUNNonlinearSolver`` object and
  initializes its content pointer and the function pointers in the operations
  structure to ``NULL``.

  **Return value:**
     If successful, this function returns a ``SUNNonlinearSolver`` object.
     If an error occurs when allocating the object, then this routine will
     return ``NULL``.

.. c:function:: void SUNNonlinSolFreeEmpty(SUNNonlinearSolver NLS)

  This routine frees the generic ``SUNNonlinearSolver`` object, under the assumption that any
  implementation-specific data that was allocated within the underlying content structure
  has already been freed. It will additionally test whether the ops pointer is ``NULL``,
  and, if it is not, it will free it as well.

   **Arguments:**
      * *NLS* -- a SUNNonlinearSolver object


Additionally, a ``SUNNonlinearSolver`` implementation *may* do
the following:

* Define and implement additional user-callable "set" routines
  acting on the ``SUNNonlinearSolver`` object, e.g., for setting
  various configuration options to tune the performance of the
  nonlinear solve algorithm.

* Provide additional user-callable "get" routines acting on the
  ``SUNNonlinearSolver`` object, e.g., for returning various solve
  statistics.
