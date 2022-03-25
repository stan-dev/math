.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMemory:

###########################
Tools for Memory Management
###########################

To support applications which leverage memory pools, or utilize a memory
abstraction layer, SUNDIALS provides a set of utilities that we
collectively refer to as the SUNMemoryHelper API. The goal of this API
is to allow users to leverage operations defined by native SUNDIALS
data structures while allowing the user to have finer-grained control of
the memory management.

.. toctree::
   :maxdepth: 1

   SUNMemory_links.rst
