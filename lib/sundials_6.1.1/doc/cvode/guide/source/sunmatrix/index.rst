.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMatrix:

**********************
Matrix Data Structures
**********************

The SUNDIALS library comes packaged with a variety of  :c:type:`SUNMatrix`
implementations, designed for simulations requiring direct linear
solvers for problems in serial or shared-memory parallel
environments.  SUNDIALS additionally provides a simple interface for
generic matrices (akin to a C++ *abstract base class*).  All of the
major SUNDIALS packages (CVODE(s), IDA(s), KINSOL, ARKODE), are
constructed to only depend on these generic matrix operations, making
them immediately extensible to new user-defined matrix objects.  For
each of the SUNDIALS-provided matrix types, SUNDIALS also provides
:c:type:`SUNLinearSolver` implementations that factor these
matrix objects and use them in the solution of linear systems.



.. toctree::
   :maxdepth: 1

   SUNMatrix_links.rst
   CVODE_requirements
