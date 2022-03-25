.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _KINSOL.Organization:

*****************
Code Organization
*****************

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/SundialsOrganization.rst


.. _KINSOL.Organization.KINSOL:

KINSOL organization
===================

The KINSOL package is written in ANSI C. The following summarizes the basic
structure of the package, although knowledge of this structure is not necessary
for its use.

.. _KINSOL.Organization.KINSOL.Figure:
.. figure:: /figs/kinsol/kinorg.png
   :align: center

   Overall structure diagram of the KINSOL package. Components specific to KINSOL
   begin with “KINSOL” (KINLS and KINSOLBBDPRE), all other items correspond
   to generic SUNDIALS vector, matrix, and solver interfaces.

The overall organization of the KINSOL package is shown in
:numref:`KINSOL.Organization.KINSOL.Figure`. KINSOL utilizes generic linear solvers
defined by the ``SUNLinearSolver`` (see :numref:`SUNLinSol`). As such, KINSOL
has no knowledge of the method being used to solve the linear and nonlinear
systems that arise. For any given user problem, there exists a single nonlinear
solver interface and, if necessary, one of the linear system solver interfaces
is specified, and invoked as needed during the integration.

KINSOL has a single unified linear solver interface, KINSOLLS, supporting both direct
and iterative linear solvers built using the generic ``SUNLinearSolver``
interface (see :numref:`SUNLinSol`). These solvers may utilize a ``SUNMatrix``
object (see :numref:`SUNMatrix`) for storing Jacobian information, or they may
be matrix-free. Since KINSOL can operate on any valid ``SUNLinearSolver``, the set
of linear solver modules available to KINSOL will expand as new ``SUNLinearSolver``
implementations are developed.

For users employing :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or
:ref:`SUNMATRIX_BAND <SUNMatrix.Band>` Jacobian matrices, KINSOL includes algorithms
for their approximation through difference quotients, although the user also has
the option of supplying a routine to compute the Jacobian (or an approximation
to it) directly. This user-supplied routine is required when using sparse or
user-supplied Jacobian matrices.

For users employing matrix-free iterative linear solvers, KINSOL includes an
algorithm for the approximation by difference quotients of the product
:math:`Jv`. Again, the user has the option of providing routines for this
operation, in two phases: setup (preprocessing of Jacobian data) and
multiplication.

For preconditioned iterative methods, the preconditioning must be supplied by
the user, again in two phases: setup and solve. While there is no default choice
of preconditioner analogous to the difference-quotient approximation in the
direct case, the references :cite:p:`BrHi:89,Byr:92`, together with the example
and demonstration programs included with KINSOL, offer considerable assistance in
building preconditioners.

KINSOL’s linear solver interface consists of four primary phases, devoted to (1)
memory allocation and initialization, (2) setup of the matrix data involved, (3)
solution of the system, and (4) freeing of memory. The setup and solution phases
are separate because the evaluation of Jacobians and preconditioners is done
only periodically during the integration, and only as required to achieve
convergence.  The call list within the central KINSOL module to each of the four
associated functions is fixed, thus allowing the central module to be completely
independent of the linear system method.

KINSOL also provides a preconditioner module, for use with any of the Krylov
iterative linear solvers. It works in conjunction with the
:ref:`NVECTOR_PARALLEL <NVectors.NVParallel>` and generates a preconditioner
that is a block-diagonal matrix with each block being a banded matrix.

All state information used by KINSOL to solve a given problem is stored in
``N_Vector`` instances. There is no global data in the KINSOL package, and so, in
this respect, it is reentrant. State information specific to the linear and
nonlinear solver are saved in the ``SUNLinearSolver`` and ``SUNNonlinearSolver``
instances respectively. The reentrancy of KINSOL enables the setting where two or
more problems are solved by intermixed or parallel calls to different instances
of the package from within a single user program.
