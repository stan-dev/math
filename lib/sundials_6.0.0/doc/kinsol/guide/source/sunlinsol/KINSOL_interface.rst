.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNLinSol.KINSOL:

KINSOL SUNLinearSolver interface
================================

:numref:`SUNLinSol.KINSOL.Table` below lists the ``SUNLinearSolver`` module linear
solver functions used within the KINLS interface. As with the ``SUNMatrix``
module, we emphasize that the KINSOL user does not need to know detailed usage
of linear solver functions by the KINSOL code modules in order to use KINSOL.
The information is presented as an implementation detail for the interested
reader.

The linear solver functions listed below are marked with "x" to indicate that they
are required, or with ":math:`\dagger`" to indicate that they are only called if
they are non-``NULL`` in the ``SUNLinearSolver`` implementation that is being
used. Note:

#. :c:func:`SUNLinSolNumIters` is only used to accumulate overall iterative linear
   solver statistics. If it is not implemented by the ``SUNLinearSolver``
   module, then KINLS will consider all solves as requiring zero iterations.

#. Although :c:func:`SUNLinSolResNorm` is optional, if it is not implemented by the
   ``SUNLinearSolver`` then KINLS will consider all solves a being *exact*.

#. Although KINLS does not call :c:func:`SUNLinSolLastFlag` directly, this routine is
   available for users to query linear solver issues directly.

#. Although KINLS does not call :c:func:`SUNLinSolFree` directly, this routine should
   be available for users to call when cleaning up from a simulation.

.. _SUNLinSol.KINSOL.Table:
.. table:: List of linear solver function usage in the KINLS interface

      ================================= =============== =============== ================
                                        DIRECT          ITERATIVE       MATRIX_ITERATIVE
      ================================= =============== =============== ================
      ``SUNLinSolGetType``              x               x               x
      ``SUNLinSolSetATimes``            :math:`\dagger` x               :math:`\dagger`
      ``SUNLinSolSetPreconditioner``    :math:`\dagger` :math:`\dagger` :math:`\dagger`
      ``SUNLinSolSetScalingVectors``    :math:`\dagger` :math:`\dagger` :math:`\dagger`
      ``SUNLinSolInitialize``           x               x               x
      ``SUNLinSolSetup``                x               x               x
      ``SUNLinSolSolve``                x               x               x
      :math:`^1`\ ``SUNLinSolNumIters``                 :math:`\dagger` :math:`\dagger`
      :math:`^2`\ ``SUNLinSolResNorm``                  :math:`\dagger` :math:`\dagger`
      :math:`^3`\ ``SUNLinSolLastFlag``
      :math:`^4`\ ``SUNLinSolFree``
      ``SUNLinSolSpace``                :math:`\dagger` :math:`\dagger` :math:`\dagger`
      ================================= =============== =============== ================


Since there are a wide range of potential ``SUNLinearSolver`` use cases, the following subsections describe some details of
the KINLS interface, in the case that interested users wish to develop custom ``SUNLinearSolver`` modules.

.. _SUNLinSol.KINSOL.lagged:

Lagged matrix information
-------------------------

If the ``SUNLinearSolver`` object self-identifies as having type
``SUNLINEARSOLVER_DIRECT`` or ``SUNLINEARSOLVER_MATRIX_ITERATIVE``, then the
``SUNLinearSolver`` object solves a linear system *defined* by a
``SUNMatrix`` object. As a result, KINSOL can perform its optional residual
monitoring scheme, described in :numref:`KINSOL.Mathematics.ModifiedNewtonResidualMon`.

.. _SUNLinSol.KINSOL.Iterative.tolerance:

Iterative linear solver tolerance
---------------------------------

If the ``SUNLinearSolver`` object self-identifies as having type
``SUNLINEARSOLVER_ITERATIVE`` or ``SUNLINEARSOLVER_MATRIX_ITERATIVE`` then KINLS
will adjust the linear solver tolerance ``delta`` as described in
:numref:`KINSOL.Mathematics.InexactNewtonStopCrit` during the course of
the nonlinear solve process. However, if the iterative linear solver does not
support scaling matrices (i.e., the ``SUNLinSolSetScalingVectors`` routine is
``NULL``), then KINLS will be unable to fully handle ill-conditioning in the
nonlinear solve process through the solution and residual scaling operators
described in :numref:`KINSOL.Mathematics.Scaling`. In this case, KINLS will attempt
to adjust the linear solver tolerance to account for this lack of functionality.
To this end, the following assumptions are made:

#. All residual components have similar magnitude; hence the scaling matrix
   :math:`D_F` used in computing the linear residual norm (see
   :numref:`KINSOL.Mathematics.Scaling`) should satisfy the assumption

   .. math:: (D_F)_{i,i} \approx D_{F,mean},\quad \text{for}\quad i=0,\ldots,n-1.

#. The ``SUNLinearSolver`` object uses a standard 2-norm to measure convergence.

Since KINSOL uses :math:`D_F` as the left-scaling matrix, :math:`S_1 = D_F`, then the linear solver convergence
requirement is converted as follows (using the notation from equations
:eq:`eq:transformed_linear_system` -- :eq:`eq:transformed_linear_system_components`:

.. math::

   &\| \tilde{b} - \tilde{A} \tilde{x} \|_2  <  \text{tol}\\
   \Leftrightarrow \quad & \| D_F P_1^{-1} b - D_F P_1^{-1} A x \|_2  <  \text{tol}\\
   \Leftrightarrow \quad & \sum_{i=0}^{n-1} \left[(D_F)_{i,i} \left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \text{tol}^2\\
   \Leftrightarrow \quad & D_{F,mean}^2 \sum_{i=0}^{n-1} \left[\left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \text{tol}^2\\
   \Leftrightarrow \quad & \sum_{i=0}^{n-1} \left[\left(P_1^{-1} (b - A x)\right)_i\right]^2  <  \left(\frac{\text{tol}}{D_{F,mean}}\right)^2\\
   \Leftrightarrow \quad & \| P_1^{-1} (b - A x)\|_2  <  \frac{\text{tol}}{D_{F,mean}}

Therefore the tolerance scaling factor

.. math:: D_{F,mean} = \frac{1}{\sqrt{n}}\left(\sum_{i=0}^{n-1} (D_F)_{i,i}^2\right)^{1/2}

is computed and the scaled tolerance ``delta``\ :math:`= \text{tol} / D_{F,mean}` is supplied to the ``SUNLinearSolver``
object.

.. _SUNLinSol.KINSOL.matrix_embedded:

Matrix-embedded solver incompatibility
--------------------------------------

At present, KINLS is incompatible with ``SUNLinearSolver`` objects that
self-identify as having type ``SUNLINEARSOLVER_MATRIX_EMBEDDED``. Support for
such user-supplied linear solvers may be added in a future release. Users
interested in such support are recommended to contact the SUNDIALS development
team.
