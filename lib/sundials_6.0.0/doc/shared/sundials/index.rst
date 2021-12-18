.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNDIALS:

==============
Using SUNDIALS
==============

As discussed in :numref:`Organization`, the six solver packages
(CVODE(S), IDA(S), ARKODE, KINSOL) that make up SUNDIALS are built upon common
classes/modules for vectors, matrices, and algebraic solvers. In addition, the
six packages all leverage some other common infrastructure, which we discuss
in this section.

.. toctree::
   SUNContext
   Profiling
   version_information
   GPU
