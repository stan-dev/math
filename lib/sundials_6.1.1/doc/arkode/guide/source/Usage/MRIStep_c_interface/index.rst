.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.MRIStep:

==========================================
Using the MRIStep time-stepping module
==========================================

This chapter is concerned with the use of the MRIStep time-stepping module for
the solution of multirate initial value problems (IVPs) of the form
:eq:`ARKODE_IVP_two_rate` in a C or C++ language setting. The following sections
discuss the header files and the layout of the user's main program, and provide
descriptions of the MRIStep user-callable functions and user-supplied functions.

The example programs located in the source code ``examples/arkode``
folder, including those described in the companion document :cite:p:`arkode_ex`,
may be helpful as templates for new codes.

MRIStep uses the input and output constants from the shared ARKODE
infrastructure. These are defined as needed in this chapter, but for
convenience the full list is provided separately in
:numref:`ARKODE.Constants`.

The relevant information on using MRIStep's C and C++ interfaces is
detailed in the following subsections.

.. toctree::
   :maxdepth: 1

   Skeleton
   User_callable
   MRIStepCoupling
   Custom_Inner_Stepper/index
