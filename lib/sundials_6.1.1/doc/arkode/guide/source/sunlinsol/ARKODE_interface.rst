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

.. _SUNLinSol.ARKODE:

ARKODE SUNLinearSolver interface
==============================================

In :numref:`SUNLinSol.ARKODE.Usage`, we list the SUNLinSol module functions used
within the ARKLS interface.  As with the SUNMATRIX module, we
emphasize that the ARKODE user does not need to know detailed usage of linear
solver functions by the ARKODE code modules in order to use ARKODE. The
information is presented as an implementation detail for the interested reader.

.. _SUNLinSol.ARKODE.Usage:
.. table:: List of SUNLinSol functions called by the ARKODE linear solver
           interface, depending on the self-identified "type" reported from
           :c:func:`SUNLinSolGetType`.  Functions marked with "X" are required;
           functions marked with "O" are only called if they are non-``NULL``
           in the ``SUNLinearSolver`` implementation that is being used.
   :align: center

   +---------------------------------------+--------+-----------+-----------+----------+
   | Routine                               | DIRECT | ITERATIVE | MATRIX    | MATRIX   |
   |                                       |        |           | ITERATIVE | EMBEDDED |
   +=======================================+========+===========+===========+==========+
   | :c:func:`SUNLinSolGetType`            | X      | X         | X         | X        |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolSetATimes`          | O      | X         | O         |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolSetPreconditioner`  | O      | O         | O         |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolSetScalingVectors`  | O      | O         | O         |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolInitialize`         | X      | X         | X         |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolSetup`              | X      | X         | X         |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolSolve`              | X      | X         | X         | X        |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolNumIters`\ :sup:`1` |        | O         | O         |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolResNorm`\ :sup:`2`  |        | O         | O         |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolLastFlag`\ :sup:`3` |        |           |           |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolFree`\ :sup:`4`     |        |           |           |          |
   +---------------------------------------+--------+-----------+-----------+----------+
   | :c:func:`SUNLinSolSpace`              | O      | O         | O         | O        |
   +---------------------------------------+--------+-----------+-----------+----------+


Notes:

1. :c:func:`SUNLinSolNumIters()` is only used to accumulate overall
   iterative linear solver statistics.  If it is not implemented by
   the ``SUNLinearSolver`` module, then ARKLS will consider all
   solves as requiring zero iterations.

2. Although :c:func:`SUNLinSolResNorm()` is optional, if it is not
   implemented by the ``SUNLinearSolver`` then ARKLS will consider all
   solves a being *exact*.

3. Although ARKLS does not call :c:func:`SUNLinSolLastFlag()`
   directly, this routine is available for users to query linear
   solver failure modes.

4. Although ARKLS does not call :c:func:`SUNLinSolFree()`
   directly, this routine should be available for users to call when
   cleaning up from a simulation.


Since there are a wide range of potential SUNLinSol use cases, the following
subsections describe some details of the ARKLS interface, in the case that
interested users wish to develop custom SUNLinSol modules.


.. _SUNLinSol.Lagged_matrix:

Lagged matrix information
---------------------------------------------------

If the SUNLinSol module identifies as having type
``SUNLINEARSOLVER_DIRECT`` or ``SUNLINEARSOLVER_MATRIX_ITERATIVE``,
then it solves a linear system *defined* by a SUNMATRIX object. ARKLS
will update the matrix information infrequently according to the strategies
outlined in :numref:`ARKODE.Mathematics.Linear.Setup`.  To this end, we
differentiate between the *desired* linear system
:math:`\mathcal A x = b` with :math:`\mathcal A = (M-\gamma J)`
and the *actual* linear system

.. math::
   \tilde{\mathcal A} \tilde{x} = b \quad\Leftrightarrow\quad (M-\tilde{\gamma} J)\tilde{x} = b.

Since ARKLS updates the SUNMATRIX object infrequently, it is likely
that :math:`\gamma\ne\tilde{\gamma}`, and in turn :math:`\mathcal
A\ne\tilde{\mathcal A}`.  Therefore, after calling the
SUNLinSol-provided :c:func:`SUNLinSolSolve()` routine, we test whether
:math:`\gamma / \tilde{\gamma} \ne 1`, and if this is the case we
scale the solution :math:`\tilde{x}` to obtain the desired linear
system solution :math:`x` via

.. math::
   x = \frac{2}{1 + \gamma / \tilde{\gamma}} \tilde{x}.
   :label: ARKODE_rescaling

The motivation for this selection of the scaling factor
:math:`c = 2/(1 + \gamma/\tilde{\gamma})` follows the derivation in
:cite:p:`BBH:89,Hin:00`.  In short, if we consider a stationary
iteration for the linear system as consisting of a solve with
:math:`\tilde{\mathcal A}` followed with a scaling by :math:`c`,
then for a linear constant-coefficient problem, the error in the
solution vector will be reduced at each iteration by the error matrix
:math:`E = I - c \tilde{\mathcal A}^{-1} \mathcal A`, with a
convergence rate given by the spectral radius of :math:`E`.  Assuming
that stiff systems have a spectrum spread widely over the left
half-plane, :math:`c` is chosen to minimize the magnitude of the
eigenvalues of :math:`E`.


.. _SUNLinSol.Iterative.Tolerance:

Iterative linear solver tolerance
---------------------------------------------------

If the SUNLinSol object self-identifies as having type
``SUNLINEARSOLVER_ITERATIVE`` or ``SUNLINEARSOLVER_MATRIX_ITERATIVE``,
then ARKLS will set the input tolerance ``delta`` as described in
:numref:`ARKODE.Mathematics.Error.Linear`.  However, if the iterative linear
solver does not support scaling matrices (i.e., the
:c:func:`SUNLinSolSetScalingVectors()` routine is ``NULL``), then
ARKLS will attempt to adjust the linear solver tolerance to account
for this lack of functionality.  To this end, the following
assumptions are made:

* All solution components have similar magnitude; hence the residual
  weight vector :math:`w` used in the WRMS norm (see
  :numref:`ARKODE.Mathematics.Error.Norm`), corresponding to the left scaling
  matrix :math:`S_1`, should satisfy the assumption

  .. math::
     w_i \approx w_{mean},\quad \text{for}\quad i=0,\ldots,n-1.

* The SUNLinSol object uses a standard 2-norm to measure convergence.

Under these assumptions, ARKLS adjusts the linear solver
convergence requirement as follows
(using the notation from :eq:`eq:transformed_linear_system_components`):

.. math::
   &\| \tilde{b} - \tilde{A} \tilde{x} \|_2  <  \text{tol}\\
   \Leftrightarrow \quad & \| S_1 P_1^{-1} b - S_1 P_1^{-1} A x \|_2  <  \text{tol}\\
   \Leftrightarrow \quad & \sum_{i=0}^{n-1} \left[w_i \left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \text{tol}^2\\
   \Leftrightarrow \quad & w_{mean}^2 \sum_{i=0}^{n-1} \left[\left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \text{tol}^2\\
   \Leftrightarrow \quad & \sum_{i=0}^{n-1} \left[\left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \left(\frac{\text{tol}}{w_{mean}}\right)^2\\
   \Leftrightarrow \quad & \| P_1^{-1} (b - A x)\|_2  <  \frac{\text{tol}}{w_{mean}}

Therefore we compute the tolerance scaling factor

.. math::
   w_{mean} = \|w\|_2 / \sqrt{n}

and supply the scaled tolerance ``delta`` :math:`= \text{tol} / w_{mean}` to the SUNLinSol object.



.. _SUNLinSol.Custom:

Providing a custom SUNLinearSolver
-------------------------------------

In certain instances, users may wish to provide a custom SUNLinSol
implementation to ARKODE in order to leverage the structure of a problem.  While
the "standard" API for these routines is typically sufficient for most users,
others may need additional ARKODE-specific information on top of what is
provided.  For these purposes, we note the following advanced ouptut functions
available in ARKStep and MRIStep:


**ARKStep advanced outputs**: when solving the Newton nonlinear system of
equations in predictor-corrector form,

.. math::
   \begin{array}{ll}
   G(z_{cor}) \equiv z_{cor} - \gamma f^I\left(t^I_{n,i}, z_{i} \right) - \tilde{a}_i = 0 &\qquad  \text{[$M=I$]},\\
   G(z_{cor}) \equiv M z_{cor} - \gamma f^I\left(t^I_{n,i}, z_{i} \right) - \tilde{a}_i = 0 &\qquad  \text{[$M$ static]},\\
   G(z_{cor}) \equiv M(t^I_{n,i}) (z_{cor} - \tilde{a}_i) - \gamma f^I\left(t^I_{n,i}, z_{i}\right) = 0 &\qquad \text{[$M$ time-dependent]}.
   \end{array}

* :c:func:`ARKStepGetCurrentTime()` -- when called within the computation of a
  step (i.e., within a solve) this returns :math:`t^I_{n,i}`. Otherwise the
  current internal solution time is returned.
* :c:func:`ARKStepGetCurrentState()` -- when called within the computation of a
  step (i.e., within a solve) this returns the current stage vector
  :math:`z_{i} = z_{cor} + z_{pred}`. Otherwise the current internal solution
  is returned.
* :c:func:`ARKStepGetCurrentGamma()` -- returns :math:`\gamma`.
* :c:func:`ARKStepGetCurrentMassMatrix()` -- returns :math:`M(t)`.
* :c:func:`ARKStepGetNonlinearSystemData()` -- returns
  :math:`z_{i}`, :math:`z_{pred}`, :math:`f^I(t^I_{n,i}, y_{cur})`,
  :math:`\tilde{a}_i`, and :math:`\gamma`.


**MRIStep advanced outputs**: when solving the Newton nonlinear system of
equations in predictor-corrector form,

.. math::
   G(z_{cor}) \equiv z_{cor} - \gamma f^I\left(t^S_{n,i}, z_{i}\right) - \tilde{a}_i = 0

* :c:func:`MRIStepGetCurrentTime()` -- when called within the computation of a
  step (i.e., within a solve) this returns :math:`t^S_{n,i}`. Otherwise the
  current internal solution time is returned.
* :c:func:`MRIStepGetCurrentState()` -- when called within the computation of a
  step (i.e., within a solve) this returns the current stage vector
  :math:`z_{i} = z_{cor} + z_{pred}`. Otherwise the current internal solution
  is returned.
* :c:func:`MRIStepGetCurrentGamma()` -- returns :math:`\gamma`.
* :c:func:`MRIStepGetNonlinearSystemData()` -- returns
  :math:`z_{i}`, :math:`z_{pred}`, :math:`f^I(t^I_{n,i}, y_{cur})`,
  :math:`\tilde{a}_i`, and :math:`\gamma`.
