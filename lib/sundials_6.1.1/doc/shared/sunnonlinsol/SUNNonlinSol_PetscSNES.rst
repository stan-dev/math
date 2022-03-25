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

.. _SUNNonlinSol.PetscSNES:

================================================
The SUNNonlinSol_PetscSNES implementation
================================================

This section describes the SUNNonlinSol interface to the
`PETSc SNES nonlinear solver(s) <https://petsc.org/release/docs/manual/snes/>`_.
To enable the SUNonlinSol_PetscSNES module, SUNDIALS must be
configured to use PETSc. Instructions on how to do this are given in
:numref:`Installation.CMake.ExternalLibraries.PETSc`. To access the
SUNNonlinSol_PetscSNES module, include the header file
``sunnonlinsol/sunnonlinsol_petscsnes.h``. The library to link to is
``libsundials_sunnonlinsolpetsc.lib`` where ``.lib`` is typically ``.so`` for
shared libaries and ``.a`` for static libraries. Users of the
SUNNonlinSol_PetscSNES module should also see :numref:`NVectors.NVPETSc`
which discusses the NVECTOR interface to the PETSc ``Vec`` API.

.. _SUNNonlinSol.PetscSNES.Description:

SUNNonlinSol_PetscSNES description
----------------------------------------

The SUNNonlinSol_PetscSNES implementation allows users to utilize a
PETSc SNES nonlinear solver to solve the nonlinear systems that arise in the
SUNDIALS integrators. Since SNES uses the KSP linear solver interface underneath
it, the SUNNonlinSol_PetscSNES implementation does not interface with
SUNDIALS linear solvers. Instead, users should set nonlinear solver options,
linear solver options, and preconditioner options through the PETSc SNES, KSP,
and PC APIs.

*Important usage notes for the SUNNonlinSol_PetscSNES implementation:*

* The SUNNonlinSol_PetscSNES implementation handles calling
  ``SNESSetFunction`` at construction. The actual residual function :math:`F(y)`
  is set by the SUNDIALS integrator when the SUNNonlinSol_PetscSNES
  object is attached to it. Therefore, a user should not call ``SNESSetFunction``
  on a ``SNES`` object that is being used with SUNNonlinSol_PetscSNES.
  For these reasons it is recommended, although not always necessary, that the
  user calls :c:func:`SUNNonlinSol_PetscSNES` with the new ``SNES`` object immediately
  after calling ``SNESCreate``.

* The number of nonlinear iterations is tracked by SUNDIALS separately from the
  count kept by SNES. As such, the function :c:func:`SUNNonlinSolGetNumIters` reports
  the cumulative number of iterations across the lifetime of the
  :c:type:`SUNNonlinearSolver` object.

* Some "converged" and "diverged" convergence reasons returned by SNES are
  treated as recoverable convergence failures by SUNDIALS. Therefore, the count of
  convergence failures returned by :c:func:`SUNNonlinSolGetNumConvFails` will reflect
  the number of recoverable convergence failures as determined by SUNDIALS, and
  may differ from the count returned by ``SNESGetNonlinearStepFailures``.

* The SUNNonlinSol_PetscSNES module is not currently compatible with
  the CVODES or IDAS staggered or simultaneous sensitivity strategies.


.. _SUNNonlinSolPetscSNES.functions:

SUNNonlinearSolver_PetscSNES functions
--------------------------------------

The SUNNonlinSol_PetscSNES module provides the following constructor
for creating a :c:type:`SUNNonlinearSolver` object.

.. c:function:: SUNNonlinearSolver SUNNonlinSol_PetscSNES(N_Vector y, SNES snes, SUNContext sunctx)

  This creates a SUNNonlinSol object that wraps a PETSc ``SNES`` object for
  use with SUNDIALS. This will call ``SNESSetFunction`` on the provided
  ``SNES`` object.

  **Arguments:**
    * *snes* -- a PETSc ``SNES`` object.
    * *y* -- a ``N_Vector`` object of type NVECTOR_PETSC that is used as a template
      for the residual vector.
    * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

  **Return value:**
     A SUNNonlinSol object if the constructor exits successfully,
     otherwise it will be ``NULL``.

  .. warning::

     This function calls ``SNESSetFunction`` and will overwrite whatever
     function was previously set. Users should not call ``SNESSetFunction``
     on the ``SNES`` object provided to the constructor.


The SUNNonlinSol_PetscSNES module implements all of the functions defined in
:numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn` except for
:c:func:`SUNNonlinSolSetup`, :c:func:`SUNNonlinSolSetLSetupFn`,
:c:func:`SUNNonlinSolSetLSolveFn`, :c:func:`SUNNonlinSolSetConvTestFn`, and
:c:func:`SUNNonlinSolSetMaxIters`.

The SUNNonlinSol_PetscSNES functions have the same names as those defined by
the generic SUNNonlinSol API with ``_PetscSNES`` appended to the
function name. Unless using the SUNNonlinSol_PetscSNES module as a
standalone nonlinear solver the generic functions defined in
:numref:`SUNNonlinSol.API.CoreFn`--:numref:`SUNNonlinSol.API.GetFn` should
be called in favor of the SUNNonlinSol_PetscSNES specific implementations.

The SUNNonlinSol_PetscSNES module also defines the following
user-callable functions.

.. c:function:: int SUNNonlinSolGetSNES_PetscSNES(SUNNonlinearSolver NLS, SNES* snes)

  This gets the ``SNES`` object that was wrapped.

  **Arguments:**
    * *NLS* -- a SUNNonlinSol object.
    * *snes* -- a pointer to a PETSc ``SNES`` object that will be set upon return.

  **Return value:**
     The return value (of type ``int``) should be zero for a successful call,
     and a negative value for a failure.


.. c:function:: int SUNNonlinSolGetPetscError_PetscSNES(SUNNonlinearSolver NLS, PestcErrorCode* error)

  This gets the last error code returned by the last internal call to a PETSc API function.

  **Arguments:**
    * *NLS* -- a SUNNonlinSol object.
    * *error* -- a pointer to a PETSc error integer that will be set upon return.

  **Return value:**
     The return value (of type ``int``) should be zero for a successful call,
     and a negative value for a failure.


.. c:function:: int SUNNonlinSolGetSysFn_PetscSNES(SUNNonlinearSolver NLS, SUNNonlinSolSysFn* SysFn)

  This returns the residual function that defines the nonlinear system.

  **Arguments:**
    * *NLS* -- a SUNNonlinSol object.
    * *SysFn* -- the function defining the nonlinear system.

  **Return value:**
     The return value (of type ``int``) should be zero for a successful call,
     and a negative value for a failure.


.. _SUNNonlinSolPetscSNES.Content:

SUNNonlinearSolver_PetscSNES content
------------------------------------

The *content* field of the SUNNonlinSol_PetscSNES module is the following
structure.

.. code-block:: c

  struct _SUNNonlinearSolverContent_PetscSNES {
    int sysfn_last_err;
    PetscErrorCode petsc_last_err;
    long int nconvfails;
    long int nni;
    void *imem;
    SNES snes;
    Vec r;
    N_Vector y, f;
    SUNNonlinSolSysFn Sys;
  };

These entries of the *content* field contain the following information:

* ``sysfn_last_err``  -- last error returned by the system defining function,
* ``petsc_last_err``  -- last error returned by PETSc,
* ``nconvfails``      -- number of nonlinear converge failures (recoverable or not),
* ``nni``             -- number of nonlinear iterations,
* ``imem``            -- SUNDIALS integrator memory,
* ``snes``            -- PETSc ``SNES`` object,
* ``r``               -- the nonlinear residual,
* ``y``               -- wrapper for PETSc vectors used in the system function,
* ``f``               -- wrapper for PETSc vectors used in the system function,
* ``Sys``             -- nonlinear system definining function.
