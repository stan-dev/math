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

.. _ARKODE.Usage.ARKStep:

======================================
Using the ARKStep time-stepping module
======================================

This chapter is concerned with the use of the ARKStep time-stepping
module for the solution of initial value problems (IVPs) in a C or C++
language setting.  The following sections discuss the header files and
the layout of the user's main program, and provide descriptions of the
ARKStep user-callable functions and user-supplied functions.

The example programs located in the source code ``examples/arkode``
folder, including those described in the companion document :cite:p:`arkode_ex`,
may be helpful as templates for new codes.

Users with applications written in Fortran should see the chapter
:numref:`SUNDIALS.Fortran`, which describes the Fortran/C interface
module for ARKStep, and may look to the Fortran example programs also
provided in the ARKODE ``examples`` directory.

The user should be aware that not all SUNLINSOL, SUNMATRIX, and
preconditioning modules are compatible with all NVECTOR
implementations.  Details on compatibility are given in the
documentation for each SUNMATRIX (see :numref:`SUNMatrix`) and each
SUNLINSOL module (see :numref:`SUNLinSol`). For example, NVECTOR_PARALLEL
is not compatible with the dense, banded, or sparse SUNMATRIX types,
or with the corresponding dense, banded, or sparse SUNLINSOL modules.
Please check :numref:`SUNMatrix` and :numref:`SUNLinSol` to
verify compatibility between these modules.  In addition to that
documentation, we note that the ARKBANDPRE preconditioning module is
only compatible with the NVECTOR_SERIAL, NVECTOR_OPENMP or
NVECTOR_PTHREADS vector implementations, and the preconditioner module
ARKBBDPRE can only be used with NVECTOR_PARALLEL.

ARKStep uses various input and output constants from the shared ARKODE
infrastructure. These are defined as needed in this chapter, but for
convenience the full list is provided separately in :numref:`ARKODE.Constants`.

The relevant information on using ARKStep's C and C++ interfaces is
detailed in the following subsections.

.. toctree::
   :maxdepth: 1

   Skeleton
   User_callable
   Preconditioners
   XBraid
