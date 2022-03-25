..
   Programmer(s): David J. Gardner @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _SUNLinSol.OneMklDense:

The SUNLinSol_OneMklDense Module
================================

The SUNLinearSolver_OneMklDense implementation of the ``SUNLinearSolver`` class
interfaces to the direct linear solvers from the
`Intel oneAPI Math Kernel Library (oneMKL) <https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html>`_
for solving dense systems or block-diagonal systems with dense blocks. This
linear solver is best paired with the SUNMatrix_OneMklDense matrix.

The header file to include when using this class is
``sunlinsol/sunlinsol_onemkldense.h``. The installed library to link to is
``libsundials_sunlinsolonemkldense.lib`` where ``lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

.. warning::

   The SUNLinearSolver_OneMklDense class is experimental and subject to change.


SUNLinearSolver_OneMklDense Functions
-------------------------------------

The SUNLinearSolver_OneMklDense class defines implementations of all "direct"
linear solver operations listed in :numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_OneMklDense`` -- returns ``SUNLINEARSOLVER_ONEMKLDENSE``
* ``SUNLinSolInitialize_OneMklDense``
* ``SUNLinSolSetup_OneMklDense``
* ``SUNLinSolSolve_OneMklDense``
* ``SUNLinSolLastFlag_OneMklDense``
* ``SUNLinSolFree_OneMklDense``

In addition, the class provides the following user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_OneMklDense(N_Vector y, SUNMatrix A, SUNContext sunctx)

   This constructor function creates and allocates memory for a
   ``SUNLinearSolver`` object.

   **Arguments:**
      * *y* -- a vector for checking compatibility with the solver.
      * *A* -- a SUNMatrix_OneMklDense matrix for checking compatibility with
        the solver.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNLinearSolver`` object. If either *A* or *y* are
      incompatible then this routine will return ``NULL``. This routine analyzes
      the input matrix and vector to determine the linear system size and to
      assess compatibility with the solver.


SUNLinearSolver_OneMklDense Usage Notes
---------------------------------------

.. warning::

   The SUNLinearSolver_OneMklDense class only supports 64-bit indexing, thus
   SUNDIALS must be built for 64-bit indexing to use this class.

   When using the SUNLinearSolver_OneMklDense class with a SUNDIALS package
   (e.g. CVODE), the queue given to the matrix is also used for the linear solver.
