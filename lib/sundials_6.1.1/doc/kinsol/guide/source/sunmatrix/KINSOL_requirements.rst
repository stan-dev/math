.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMatrix.KINSOL:

SUNMatrix functions used by KINSOL
==================================

In :numref:`SUNMatrix.KINSOL.Table`, we list the matrix functions in the
``SUNMatrix`` module used within the KINSOL package. The table also shows, for
each function, which of the code modules uses the function. The main KINSOL
integrator does not call any ``SUNMatrix`` functions directly, so the table
columns are specific to the KINLS and KINBBDPRE preconditioner modules. We
further note that the KINLS interface only utilizes these routines when supplied
with a *matrix-based* linear solver, i.e., the ``SUNMatrix`` object passed to
:c:func:`KINSetLinearSolver` was not ``NULL``.

At this point, we should emphasize that the KINSOL user does not need to know
anything about the usage of matrix functions by the KINSOL code modules in order
to use KINSOL. The information is presented as an implementation detail for the
interested reader.

.. _SUNMatrix.KINSOL.Table:
.. table:: List of matrix functions usage by KINSOL code modules

   +---------------------------+-----------------+-----------------+
   |                           |      KINLS      |    KINBBDPRE    |
   +===========================+=================+=================+
   | :c:func:`SUNMatGetID`     | x               |                 |
   +---------------------------+-----------------+-----------------+
   | :c:func:`SUNMatDestroy`   |                 | x               |
   +---------------------------+-----------------+-----------------+
   | :c:func:`SUNMatZero`      | x               | x               |
   +---------------------------+-----------------+-----------------+
   | :c:func:`SUNMatSpace`     |                 | :math:`\dagger` |
   +---------------------------+-----------------+-----------------+

The matrix functions listed with a :math:`\dagger` symbol are optionally used,
in that these are only called if they are implemented in the ``SUNMatrix``
module that is being used (i.e. their function pointers are non-``NULL``). The
matrix functions listed in :numref:`SUNMatrix.Description` that are *not* used
by KINSOL are: :c:func:`SUNMatCopy`, :c:func:`SUNMatClone`,
:c:func:`SUNMatScaleAdd`, :c:func:`SUNMatScaleAddI` and :c:func:`SUNMatMatvec`.
Therefore a user-supplied ``SUNMatrix`` module for KINSOL could omit these
functions.

We note that the KINBBDPRE preconditioner module is hard-coded to use the
SUNDIALS-supplied band ``SUNMatrix`` type, so the most useful information above
for user-supplied ``SUNMatrix`` implementations is the column relating the KINLS
requirements.
