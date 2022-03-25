.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors:

######################
Vector Data Structures
######################

The SUNDIALS library comes packaged with a variety of NVECTOR
implementations, designed for simulations in serial, shared-memory
parallel, and distributed-memory parallel environments, as well as
interfaces to vector data structures used within external linear
solver libraries.  All native implementations assume that the
process-local data is stored contiguously, and they in turn provide a
variety of standard vector algebra operations that may be performed on
the data.

In addition, SUNDIALS provides a simple interface for generic vectors
(akin to a C++ *abstract base class*).  All of the SUNDIALS packages
(CVODE(s), IDA(s), KINSOL, ARKODE) in turn are constructed to
only depend on these generic vector operations, making them immediately
extensible to new user-defined vector objects.  The only exceptions to
this rule relate to the direct linear solver modules (and associated
matrices), since they rely on particular data storage and access
patterns in the NVECTORS used.



.. toctree::
   :maxdepth: 1

   NVector_API_link.rst
   CVODE_requirements.rst
   NVector_links.rst
