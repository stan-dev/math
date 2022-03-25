.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _CVODE.Organization:

*****************
Code Organization
*****************

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/SundialsOrganization.rst


.. _CVODE.Organization.CVODE:

CVODE organization
==================

The CVODE package is written in ANSI C. The following summarizes
the basic structure of the package, although knowledge of this structure
is not necessary for its use.

The overall organization of the CVODE package is shown in Figure
:numref:`CVODE.Organization.CVODE.Figure`.

.. _CVODE.Organization.CVODE.Figure:
.. figure:: /figs/cvode/cvorg.png
   :align: center

   Overall structure diagram of the CVODE package. Modules
   specific to CVODE begin with “CV” (CVLS, CVNLS, CVDIAG,
   CVBBDPRE, and CVBANDPRE), all other items correspond to generic
   SUNDIALS vector, matrix, and solver modules.


The central integration module, implemented in the files ``cvode.h``,
``cvode_impl.h``, and ``cvode.c``, deals with the evaluation of integration
coefficients, estimation of local error, selection of stepsize and order, and
interpolation to user output points, among other issues.

CVODE utilizes generic linear and nonlinear solver modules defined
by the ``SUNLinearSolver`` API (see :numref:`SUNLinSol`)
and ``SUNNonlinearSolver`` API (see :numref:`SUNNonlinSol`) respectively. As such, CVODE
has no knowledge of the method being used to solve the linear and
nonlinear systems that arise. For any given user problem, there exists a
single nonlinear solver interface and, if necessary, one of the linear
system solver interfaces is specified, and invoked as needed during the
integration.

At present, the package includes two linear solver interfaces. The
primary linear solver interface, CVLS, supports both direct and
iterative linear solvers built using the generic ``SUNLinearSolver`` API (see
:numref:`SUNLinSol`). These solvers may utilize a
``SUNMatrix`` object (see :numref:`SUNMatrix`) for
storing Jacobian information, or they may be matrix-free. Since
CVODE can operate on any valid ``SUNLinearSolver`` implementation, the set
of linear solver modules available to CVODE will expand as new
``SUNLinearSolver`` modules are developed.

Additionally, CVODE includes the *diagonal* linear solver interface,
CVDIAG, that creates an internally generated diagonal approximation
to the Jacobian.

For users employing :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or
:ref:`SUNMATRIX_BAND <SUNMatrix.Band>` Jacobian matrices, CVODE
includes algorithms for their approximation through difference
quotients, although the user also has the option of supplying a routine
to compute the Jacobian (or an approximation to it) directly. This
user-supplied routine is required when using sparse or user-supplied
Jacobian matrices.

For users employing matrix-free iterative linear solvers, CVODE
includes an algorithm for the approximation by difference quotients of
the product :math:`Mv`. Again, the user has the option of providing
routines for this operation, in two phases: setup (preprocessing of
Jacobian data) and multiplication.

For preconditioned iterative methods, the preconditioning must be
supplied by the user, again in two phases: setup and solve. While there
is no default choice of preconditioner analogous to the
difference-quotient approximation in the direct case, the references
:cite:p:`BrHi:89,Byr:92`, together with the example and
demonstration programs included with CVODE, offer considerable
assistance in building preconditioners.

CVODE’s linear solver interface consists of four primary phases,
devoted to (1) memory allocation and initialization, (2) setup of the
matrix data involved, (3) solution of the system, and (4) freeing of
memory. The setup and solution phases are separate because the
evaluation of Jacobians and preconditioners is done only periodically
during the integration, and only as required to achieve convergence.

CVODE also provides two preconditioner modules, for use with any of
the Krylov iterative linear solvers. The first one, CVBANDPRE, is
intended to be used with ``NVECTOR_SERIAL``, ``NVECTOR_OPENMP`` or ``NVECTOR_PTHREADS``
and provides a banded difference-quotient Jacobian-based preconditioner,
with corresponding setup and solve routines. The second preconditioner
module, CVBBDPRE, works in conjunction with ``NVECTOR_PARALLEL`` and generates a
preconditioner that is a block-diagonal matrix with each block being a
banded matrix.

All state information used by CVODE to solve a given problem is
saved in a structure, and a pointer to that structure is returned to the
user. There is no global data in the CVODE package, and so, in this
respect, it is reentrant. State information specific to the linear
solver is saved in a separate structure, a pointer to which resides in
the CVODE memory structure. The reentrancy of CVODE was
motivated by the anticipated multicomputer extension, but is also
essential in a uniprocessor setting where two or more problems are
solved by intermixed calls to the package from within a single user
program.
