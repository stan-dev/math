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


.. _SUNMatrix.ARKODE:

SUNMATRIX functions used by ARKODE
==========================================

In Table :numref:`SUNMatrix.ARKODE.Table`, we list the matrix functions in
the ``SUNMatrix`` module used within the ARKODE package.  The table
also shows, for each function, which of the code modules uses the
function.  The main ARKODE time step modules, ARKStep, ERKStep, and MRIStep,
do not call any ``SUNMatrix`` functions directly, so the table columns
are specific to the ARKLS interface and the ARKBANDPRE and ARKBBDPRE
preconditioner modules.   We further note that the ARKLS interface
only utilizes these routines when supplied with a *matrix-based*
linear solver, i.e. the ``SUNMatrix`` object (*J* or *M*) passed to
:c:func:`ARKStepSetLinearSolver()` or
:c:func:`ARKStepSetMassLinearSolver()` was not ``NULL``.

At this point, we should emphasize that the ARKODE user does not need
to know anything about the usage of matrix functions by the ARKODE
code modules in order to use ARKODE.  The information is presented as
an implementation detail for the interested reader.


.. _SUNMatrix.ARKODE.Table:
.. table:: List of matrix functions usage by ARKODE code modules
   :align: center

   +-----------------------------+-----------------+-----------------+-----------------+
   |                             |      ARKLS      |    ARKBANDPRE   |    ARKBBDPRE    |
   +=============================+=================+=================+=================+
   | :c:func:`SUNMatGetID`       |      X          |                 |                 |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatClone`       |      X          |                 |                 |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatDestroy`     |      X          |        X        |        X        |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatZero`        |      X          |        X        |        X        |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatCopy`        |      X          |        X        |        X        |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatScaleAddI`   |      X          |        X        |        X        |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatScaleAdd`    |      1          |                 |                 |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatMatvec`      |      1          |                 |                 |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatMatvecSetup` |      1,2        |                 |                 |
   +-----------------------------+-----------------+-----------------+-----------------+
   | :c:func:`SUNMatSpace`       |      2          |        2        |        2        |
   +-----------------------------+-----------------+-----------------+-----------------+

1. These matrix functions are only used for problems involving a
   non-identity mass matrix.

2. These matrix functions are optionally used, in that these are only
   called if they are implemented in the ``SUNMatrix`` module that is
   being used (i.e. their function pointers are non-``NULL``).  If not
   supplied, these modules will assume that the matrix requires no
   storage.


We note that both the ARKBANDPRE and ARKBBDPRE preconditioner modules
are hard-coded to use the SUNDIALS-supplied band ``SUNMatrix`` type,
so the most useful information above for user-supplied ``SUNMatrix``
implementations is the column relating to ARKLS requirements.
