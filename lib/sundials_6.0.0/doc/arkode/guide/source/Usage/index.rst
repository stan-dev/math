.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage:

************
Using ARKODE
************

This chapter discusses usage for ARKODE from C, C++ and Fortran applications.
The chapter builds upon :numref:`SUNDIALS`. We first discuss commonalities to
each of ARKODE's time-stepping modules, including locations and naming
conventions for the library and header files, and discussion of data types in
SUNDIALS.  We then separately discuss the C and C++ interfaces to each of
ARKODE's time stepping modules: :ref:`ARKStep <ARKODE.Usage.ARKStep>`,
:ref:`ERKStep <ARKODE.Usage.ERKStep>`, and :ref:`MRIStep
<ARKODE.Usage.MRIStep>`. Following these, we describe set of :ref:`user-supplied
routines <ARKODE.Usage.UserSupplied>` (both required and optional) that can be
supplied to ARKODE.

.. toctree::
   :maxdepth: 1

   General.rst
   ARKStep_c_interface/index.rst
   ERKStep_c_interface/index.rst
   MRIStep_c_interface/index.rst
   User_supplied.rst
