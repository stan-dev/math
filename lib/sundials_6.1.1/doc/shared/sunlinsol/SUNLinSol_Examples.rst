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

.. _SUNLinSol.Examples:

SUNLinearSolver Examples
======================================

There are ``SUNLinearSolver`` examples that may be installed for each
implementation; these make use of the functions in ``test_sunlinsol.c``.
These example functions show simple usage of the ``SUNLinearSolver`` family
of modules.  The inputs to the examples depend on the linear solver type,
and are output to ``stdout`` if the example is run without the
appropriate number of command-line arguments.

The following is a list of the example functions in ``test_sunlinsol.c``:

* ``Test_SUNLinSolGetType``: Verifies the returned solver type against
  the value that should be returned.

* ``Test_SUNLinSolGetID``: Verifies the returned solver identifier against
  the value that should be returned.

* ``Test_SUNLinSolInitialize``: Verifies that ``SUNLinSolInitialize``
  can be called and returns successfully.

* ``Test_SUNLinSolSetup``: Verifies that ``SUNLinSolSetup`` can
  be called and returns successfully.

* ``Test_SUNLinSolSolve``: Given a ``SUNMatrix`` object :math:`A`,
  ``N_Vector`` objects :math:`x` and :math:`b` (where :math:`Ax=b`)
  and a desired solution tolerance ``tol``, this routine clones
  :math:`x` into a new vector :math:`y`, calls ``SUNLinSolSolve`` to
  fill :math:`y` as the solution to :math:`Ay=b` (to the input
  tolerance), verifies that each entry in :math:`x` and :math:`y`
  match to within ``10*tol``, and overwrites :math:`x` with :math:`y`
  prior to returning (in case the calling routine would like to
  investigate further).

* ``Test_SUNLinSolSetATimes`` (iterative solvers only): Verifies that
  ``SUNLinSolSetATimes`` can be called and returns successfully.

* ``Test_SUNLinSolSetPreconditioner`` (iterative solvers only):
  Verifies that ``SUNLinSolSetPreconditioner`` can be called and
  returns successfully.

* ``Test_SUNLinSolSetScalingVectors`` (iterative solvers only):
  Verifies that ``SUNLinSolSetScalingVectors`` can be called and
  returns successfully.

* ``Test_SUNLinSolSetZeroGuess`` (iterative solvers only): Verifies that
  ``SUNLinSolSetZeroGuess`` can be called and returns successfully.

* ``Test_SUNLinSolLastFlag``: Verifies that ``SUNLinSolLastFlag`` can
  be called, and outputs the result to ``stdout``.

* ``Test_SUNLinSolNumIters`` (iterative solvers only): Verifies that
  ``SUNLinSolNumIters`` can be called, and outputs the result to
  ``stdout``.

* ``Test_SUNLinSolResNorm`` (iterative solvers only): Verifies that
  ``SUNLinSolResNorm`` can be called, and that the result is
  non-negative.

* ``Test_SUNLinSolResid`` (iterative solvers only): Verifies that
  ``SUNLinSolResid`` can be called.

* ``Test_SUNLinSolSpace`` verifies that ``SUNLinSolSpace`` can be
  called, and outputs the results to ``stdout``.

We'll note that these tests should be performed in a particular
order.  For either direct or iterative linear
solvers, ``Test_SUNLinSolInitialize`` must be called
before ``Test_SUNLinSolSetup``, which must be called
before ``Test_SUNLinSolSolve``.  Additionally, for iterative linear
solvers ``Test_SUNLinSolSetATimes``, ``Test_SUNLinSolSetPreconditioner``
and ``Test_SUNLinSolSetScalingVectors`` should be called
before ``Test_SUNLinSolInitialize``;
similarly ``Test_SUNLinSolNumIters``, ``Test_SUNLinSolResNorm``
and ``Test_SUNLinSolResid`` should be called
after ``Test_SUNLinSolSolve``.  These are called in the appropriate
order in all of the example problems.
