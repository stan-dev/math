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

.. _SUNNonlinSol.Newton:

==============================================
The SUNNonlinSol_Newton implementation
==============================================

This section describes the SUNNonlinSol implementation of Newton's method. To
access the SUNNonlinSol_Newton module, include the header file
``sunnonlinsol/sunnonlinsol_newton.h``. We note that the SUNNonlinSol_Newton
module is accessible from SUNDIALS integrators *without* separately
linking to the ``libsundials_sunnonlinsolnewton`` module library.


.. _SUNNonlinSol.Newton.Math:

SUNNonlinSol_Newton description
----------------------------------------

To find the solution to

.. math::
   F(y) = 0
   :label: e:newton_sys

given an initial guess :math:`y^{(0)}`, Newton's method computes a series of
approximate solutions

.. math::
   y^{(m+1)} = y^{(m)} + \delta^{(m+1)}

where :math:`m` is the Newton iteration index, and the Newton update :math:`\delta^{(m+1)}`
is the solution of the linear system

.. math::
   A(y^{(m)}) \delta^{(m+1)} = -F(y^{(m)}) \, ,
   :label: e:newton_linsys


in which :math:`A` is the Jacobian matrix

.. math::
   A \equiv \partial F / \partial y \, .
   :label: e:newton_mat

Depending on the linear solver used, the SUNNonlinSol_Newton module
will employ either a Modified Newton method or an Inexact Newton
method :cite:p:`Bro:87,BrSa:90,DES:82,DeSc:96,Kel:95`. When used
with a direct linear solver, the Jacobian matrix :math:`A` is held
constant during the Newton iteration, resulting in a Modified Newton
method. With a matrix-free iterative linear solver, the iteration is
an Inexact Newton method.

In both cases, calls to the integrator-supplied :c:type:`SUNNonlinSolLSetupFn`
function are made infrequently to amortize the increased cost of
matrix operations (updating :math:`A` and its factorization within direct
linear solvers, or updating the preconditioner within iterative linear
solvers).  Specifically, SUNNonlinSol_Newton will call the
:c:type:`SUNNonlinSolLSetupFn` function in two instances:

(a) when requested by the integrator (the input ``callLSetSetup`` is
    ``SUNTRUE``) before attempting the Newton iteration, or

(b) when reattempting the nonlinear solve after a recoverable failure
    occurs in the Newton iteration with stale Jacobian information
    (``jcur`` is ``SUNFALSE``).  In this case, SUNNonlinSol_Newton
    will set ``jbad`` to ``SUNTRUE`` before calling the
    :c:type:`SUNNonlinSolLSetupFn()` function.

Whether the Jacobian matrix :math:`A` is fully or partially updated depends
on logic unique to each integrator-supplied :c:type:`SUNNonlinSolSetupFn`
routine. We refer to the discussion of nonlinear solver strategies
provided in the package-specific Mathematics section of the documentation for details.

The default maximum number of iterations and the stopping criteria for
the Newton iteration are supplied by the SUNDIALS integrator when
SUNNonlinSol_Newton is attached to it.  Both the maximum number of
iterations and the convergence test function may be modified by the
user by calling the :c:func:`SUNNonlinSolSetMaxIters` and/or
:c:func:`SUNNonlinSolSetConvTestFn` functions after attaching the
SUNNonlinSol_Newton object to the integrator.


.. _SUNNonlinSol.Newton.Functions:

SUNNonlinSol_Newton functions
---------------------------------------

The SUNNonlinSol_Newton module provides the following constructor
for creating the ``SUNNonlinearSolver`` object.


.. c:function:: SUNNonlinearSolver SUNNonlinSol_Newton(N_Vector y, SUNContext sunctx)

   This creates a ``SUNNonlinearSolver`` object for use with SUNDIALS
   integrators to solve nonlinear systems of the form :math:`F(y) = 0`
   using Newton's method.

   **Arguments:**
      * *y* -- a template for cloning vectors needed within the solver.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      A SUNNonlinSol object if the constructor exits successfully,
      otherwise it will be ``NULL``.


The SUNNonlinSol_Newton module implements all of the functions
defined in :numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`
except for :c:func:`SUNNonlinSolSetup`. The SUNNonlinSol_Newton functions
have the same names as those defined by the generic SUNNonlinSol API with
``_Newton`` appended to the function name. Unless using the SUNNonlinSol_Newton
module as a standalone nonlinear solver the generic functions defined
in :numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn`
should be called in favor of the SUNNonlinSol_Newton-specific implementations.

The SUNNonlinSol_Newton module also defines the following
user-callable function.


.. c:function:: int SUNNonlinSolGetSysFn_Newton(SUNNonlinearSolver NLS, SUNNonlinSolSysFn *SysFn)

   This returns the residual function that defines the nonlinear system.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *SysFn* -- the function defining the nonlinear system.

   **Return value:**
      The return value should be zero for a successful call, and a
      negative value for a failure.

   **Notes:**
      This function is intended for users that wish to evaluate the
      nonlinear residual in a custom convergence test function for the
      SUNNonlinSol_Newton module.  We note that SUNNonlinSol_Newton
      will not leverage the results from any user calls to *SysFn*.


.. c:function:: int SUNNonlinSolSetInfoFile_Newton(SUNNonlinearSolver NLS, FILE* info_file)

   This sets the output file where all informative (non-error) messages should be directed.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *info_file* -- pointer to output file (``stdout`` by default);
         a ``NULL`` input will disable output.

   **Return value:**
      * ``SUN_NLS_SUCCESS`` if successful.
      * ``SUN_NLS_MEM_NULL`` if the SUNNonlinSol memory was ``NULL``.
      * ``SUN_NLS_ILL_INPUT`` if SUNDIALS was not built with monitoring enabled.

   **Notes:**
      This function is intended for users that wish to monitor the nonlinear
      solver progress. By default, the file pointer is set to ``stdout``.

   .. warning::

      SUNDIALS must be built with the CMake option
      ``SUNDIALS_BUILD_WITH_MONITORING`` to utilize this function.
      See :numref:`Installation.CMake.Options` for more information.


.. c:function:: int SUNNonlinSolSetPrintLevel_Newton(SUNNonlinearSolver NLS, int print_level)

   This specifies the level of verbosity of the output.

   **Arguments:**
      * *NLS* -- a SUNNonlinSol object.
      * *print_level* -- flag indicating level of verbosity;
        must be one of:

         * 0, no information is printed (default).
         * 1, for each nonlinear iteration the residual norm is printed.

   **Return value:**
      * ``SUN_NLS_SUCCESS`` if successful.
      * ``SUN_NLS_MEM_NULL`` if the SUNNonlinearSolver memory was ``NULL``.
      * ``SUN_NLS_ILL_INPUT`` if SUNDIALS was not built with monitoring enabled,
        or the print level value was invalid.

   **Notes:**
      This function is intended for users that wish to monitor the nonlinear
      solver progress. By default, the print level is 0.

   .. warning::

      SUNDIALS must be built with the CMake option
      ``SUNDIALS_BUILD_WITH_MONITORING`` to utilize this function.
      See :numref:`Installation.CMake.Options` for more information.


.. _SUNNonlinSol.Newton.Content:

SUNNonlinSol_Newton content
------------------------------------------------

The *content* field of the SUNNonlinSol_Newton module is the
following structure.

.. code-block:: c

   struct _SUNNonlinearSolverContent_Newton {

     SUNNonlinSolSysFn      Sys;
     SUNNonlinSolLSetupFn   LSetup;
     SUNNonlinSolLSolveFn   LSolve;
     SUNNonlinSolConvTestFn CTest;

     N_Vector    delta;
     booleantype jcur;
     int         curiter;
     int         maxiters;
     long int    niters;
     long int    nconvfails;
     void*       ctest_data;

     int         print_level;
     FILE*       info_file;
   };

These entries of the *content* field contain the following
information:

* ``Sys`` -- the function for evaluating the nonlinear system,

* ``LSetup`` -- the package-supplied function for setting up the
  linear solver,

* ``LSolve`` -- the package-supplied function for performing a linear
  solve,

* ``CTest`` -- the function for checking convergence of the Newton iteration,

* ``delta`` -- the Newton iteration update vector,

* ``jcur`` -- the Jacobian status (``SUNTRUE`` = current, ``SUNFALSE`` = stale),

* ``curiter``  -- the current number of iterations in the solve attempt,

* ``maxiters`` -- the maximum number of Newton iterations allowed in a solve,

* ``niters`` -- the total number of nonlinear iterations across all solves,

* ``nconvfails`` -- the total number of nonlinear convergence failures across
  all solves,

* ``ctest_data`` -- the data pointer passed to the convergence test function,

* ``print_level`` - controls the amount of information to be printed to the info file,

* ``info_file``   - the file where all informative (non-error) messages will be directed.
