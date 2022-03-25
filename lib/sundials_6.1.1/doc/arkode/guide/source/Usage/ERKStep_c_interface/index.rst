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

.. _ARKODE.Usage.ERKStep:

======================================
Using the ERKStep time-stepping module
======================================

This chapter is concerned with the use of the ERKStep time-stepping
module for the solution of nonstiff initial value problems (IVPs) in a
C or C++ language setting.  The following sections discuss the header
files and the layout of the user's main program, and provide
descriptions of the ERKStep user-callable functions and user-supplied
functions.

The example programs described in the companion document :cite:p:`arkode_ex` may
be helpful. Those codes may be used as templates for new codes and are
included in the ARKODE package ``examples`` subdirectory.

ERKStep uses the input and output constants from the shared ARKODE
infrastructure. These are defined as needed in this chapter, but for
convenience the full list is provided separately in
:numref:`ARKODE.Constants`.

The relevant information on using ERKStep's C and C++ interfaces is
detailed in the following sub-sections.

.. toctree::
   :maxdepth: 1

   Skeleton
   User_callable
