.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNMatrix.CVODES:

SUNMatrix functions used by CVODES
==================================

In :numref:`SUNMatrix.CVODES.Table`, we list the matrix functions in the
``SUNMatrix`` module used within the CVODES package.
The table also shows, for each function, which of the code modules uses
the function. The main CVODES integrator does not call any
``SUNMatrix`` functions directly, so the table columns are specific to
the CVLS interface and the CVBANDPRE and
CVBBDPRE preconditioner modules. We further note that the CVLS
interface only utilizes these routines when supplied with a
*matrix-based* linear solver, i.e., the ``SUNMatrix`` object
passed to :c:func:`CVodeSetLinearSolver` was not ``NULL``.

At this point, we should emphasize that the CVODES user does not need to know
anything about the usage of matrix functions by the CVODES code modules in order
to use CVODES. The information is presented as an implementation detail for the
interested reader.

.. _SUNMatrix.CVODES.Table:
.. table:: List of matrix functions usage by CVODES code modules
   :align: center

   +---------------------------+-----------------+-----------------+-----------------+
   |                           |      CVLS       |    CVBANDPRE    |    CVBBDPRE     |
   +===========================+=================+=================+=================+
   | :c:func:`SUNMatClone`     | x               |                 |                 |
   +---------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatDestroy`   | x               | x               | x               |
   +---------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatZero`      | x               | x               | x               |
   +---------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatGetID`     | x               |                 |                 |
   +---------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatCopy`      | x               | x               | x               |
   +---------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatScaleAddI` | x               | x               | x               |
   +---------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatSpace`     | :math:`\dagger` | :math:`\dagger` | :math:`\dagger` |
   +---------------------------+-----------------+-----------------+-----------------+

The matrix functions listed with a :math:`\dagger` symbol are optionally used,
in that these are only called if they are implemented in the ``SUNMatrix``
module that is being used (i.e. their function pointers are non-``NULL``). The
matrix functions listed in :numref:`SUNMatrix.Description` that are *not* used by CVODES are:
:c:func:`SUNMatScaleAdd` and :c:func:`SUNMatMatvec`. Therefore a user-supplied ``SUNMatrix``
module for CVODES could omit these functions.

We note that the CVBANDPRE and CVBBDPRE preconditioner modules
are hard-coded to use the SUNDIALS-supplied band ``SUNMatrix`` type,
so the most useful information above for user-supplied ``SUNMatrix``
implementations is the column relating the CVLS requirements.
