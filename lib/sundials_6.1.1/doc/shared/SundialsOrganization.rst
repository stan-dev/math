.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

SUNDIALS consists of the solvers CVODE and ARKODE for ordinary differential
equation (ODE) systems, IDA for differential-algebraic (DAE) systems, and KINSOL
for nonlinear algebraic systems. In addition, SUNDIALS also includes variants of
CVODE and IDA with sensitivity analysis capabilities (using either forward or
adjoint methods), called CVODES and IDAS, respectively. The following is a list
summarizes the basic functionality of each SUNDIALS package:

* CVODE, a solver for stiff and nonstiff ODE systems :math:`\dot{y} = f(t,y)`
  based on Adams and BDF methods;

* CVODES, a solver for stiff and nonstiff ODE systems with sensitivity analysis
  capabilities;

* ARKODE, a solver for stiff, nonstiff, mixed stiff-nonstiff, and multirate ODE
  systems :math:`M(t)\, \dot{y} = f_1(t,y) + f_2(t,y)` based on Runge-Kutta
  methods;

* IDA, a solver for differential-algebraic systems :math:`F(t,y,\dot{y}) = 0`
  based on BDF methods;

* IDAS, a solver for differential-algebraic systems with sensitivity analysis
  capabilities;

* KINSOL, a solver for nonlinear algebraic systems :math:`F(u) = 0`.

The various packages in the suite share many common components and are organized
as a family. :numref:`Organization.Sundials.HighLevelDiagram` gives a high-level
overview of solver packages, the shared vector, matrix, linear solver, and
nonlinear solver interfaces (abstract base classes), and the corresponding class
implementations provided with SUNDIALS. For classes that provide interfaces to
third-party libraries (i.e., LAPACK, KLU, SuperLU_MT, SuperLU_DIST, *hypre*,
PETSc, Trilinos, and Raja) users will need to download and compile those
packages independently of SUNDIALS. The directory structure is shown in
:numref:`Organization.Sundials.DirectoryStructure`.

.. _Organization.Sundials.HighLevelDiagram:
.. figure:: /figs/sunorg1.png
   :align: center

   High-level diagram of the SUNDIALS suite.

.. _Organization.Sundials.DirectoryStructure:
.. figure:: /figs/sunorg2.png
   :align: center

   Directory structure of the SUNDIALS source tree.
