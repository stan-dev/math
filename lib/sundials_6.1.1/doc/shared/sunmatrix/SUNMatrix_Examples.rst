..
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

.. _SUNMatrix.Examples:

SUNMATRIX Examples
==================

There are ``SUNMatrix`` examples that may be installed for each
implementation, that make use of the functions in ``test_sunmatrix.c``.
These example functions show simple usage of the ``SUNMatrix`` family
of functions.  The inputs to the examples depend on the matrix type,
and are output to ``stdout`` if the example is run without the
appropriate number of command-line arguments.

The following is a list of the example functions in ``test_sunmatrix.c``:

* ``Test_SUNMatGetID``: Verifies the returned matrix ID against
  the value that should be returned.

* ``Test_SUNMatClone``: Creates clone of an existing matrix,
  copies the data, and checks that their values match.

* ``Test_SUNMatZero``: Zeros out an existing matrix and checks
  that each entry equals 0.0.

* ``Test_SUNMatCopy``: Clones an input matrix, copies its data
  to a clone, and verifies that all values match.

* ``Test_SUNMatScaleAdd``: Given an input matrix :math:`A` and an
  input identity matrix :math:`I`, this test clones and copies
  :math:`A` to a new matrix :math:`B`, computes :math:`B = -B+B`, and
  verifies that the resulting matrix entries equal 0.  Additionally,
  if the matrix is square, this test clones and copies :math:`A` to a
  new matrix :math:`D`, clones and copies :math:`I` to a new matrix
  :math:`C`, computes :math:`D = D+I` and :math:`C = C+A` using
  :c:func:`SUNMatScaleAdd`, and then verifies that :math:`C=D`.

* ``Test_SUNMatScaleAddI``: Given an input matrix :math:`A` and an
  input identity matrix :math:`I`, this clones and copies :math:`I` to
  a new matrix :math:`B`, computes :math:`B = -B+I` using
  :c:func:`SUNMatScaleAddI`, and verifies that the resulting matrix entries
  equal 0.

* ``Test_SUNMatMatvecSetup``: verifies that :c:func:`SUNMatMatvecSetup`
  can be called.

* ``Test_SUNMatMatvec`` Given an input matrix :math:`A` and input
  vectors :math:`x` and :math:`y` such that :math:`y=Ax`, this test
  has different behavior depending on whether :math:`A` is square.  If
  it is square, it clones and copies :math:`A` to a new matrix
  :math:`B`, computes :math:`B = 3B+I` using :c:func:`SUNMatScaleAddI`,
  clones :math:`y` to new vectors :math:`w` and :math:`z`, computes
  :math:`z = Bx` using :c:func:`SUNMatMatvec`, computes :math:`w = 3y+x`
  using ``N_VLinearSum``, and verifies that :math:`w==z`.  If
  :math:`A` is not square, it just clones :math:`y` to a new vector
  :math:`z`, `computes :math:`z=Ax` using :c:func:`SUNMatMatvec`, and
  verifies that :math:`y=z`.

* ``Test_SUNMatSpace``: verifies that :c:func:`SUNMatSpace` can be
  called, and outputs the results to ``stdout``.
