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

.. _SUNLinSol_Dense:

The SUNLinSol_Dense Module
======================================

The SUNLinSol_Dense implementation of the ``SUNLinearSolver`` class
is designed to be used with the corresponding SUNMATRIX_DENSE matrix type,
and one of the serial or shared-memory ``N_Vector`` implementations
(NVECTOR_SERIAL, NVECTOR_OPENMP or NVECTOR_PTHREADS).

.. _SUNLinSol_Dense.Usage:

SUNLinSol_Dense Usage
------------------------

The header file to be included when using this module is
``sunlinsol/sunlinsol_dense.h``.  The SUNLinSol_Dense module is
accessible from all SUNDIALS solvers *without*
linking to the ``libsundials_sunlinsoldense`` module library.


The module SUNLinSol_Dense provides the following user-callable constructor routine:


.. c:function:: SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A, SUNContext sunctx)

   This function creates and allocates memory for a dense ``SUNLinearSolver``.

   **Arguments:**
      * *y* -- vector used to determine the linear system size.
      * *A* -- matrix used to assess compatibility.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      New SUNLinSol_Dense object, or ``NULL`` if either ``A`` or ``y`` are incompatible.

   **Notes:**
      This routine will perform consistency checks to ensure that it is
      called with consistent ``N_Vector`` and ``SUNMatrix`` implementations.
      These are currently limited to the SUNMATRIX_DENSE matrix type and
      the NVECTOR_SERIAL, NVECTOR_OPENMP, and NVECTOR_PTHREADS vector types.  As
      additional compatible matrix and vector implementations are added to
      SUNDIALS, these will be included within this compatibility check.


For backwards compatibility, we also provide the following wrapper function:

.. c:function:: SUNLinearSolver SUNDenseLinearSolver(N_Vector y, SUNMatrix A)

   Wrapper function for :c:func:`SUNLinSol_Dense`, with identical input and
   output arguments



.. _SUNLinSol_Dense.Description:

SUNLinSol_Dense Description
-----------------------------


The SUNLinSol_Dense module defines the *content*
field of a ``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_Dense {
     sunindextype N;
     sunindextype *pivots;
     sunindextype last_flag;
   };

These entries of the *content* field contain the following
information:

* ``N`` - size of the linear system,

* ``pivots`` - index array for partial pivoting in LU factorization,

* ``last_flag`` - last error return flag from internal function evaluations.


This solver is constructed to perform the following operations:

* The "setup" call performs an :math:`LU` factorization with
  partial (row) pivoting (:math:`\mathcal O(N^3)` cost),
  :math:`PA=LU`, where :math:`P` is a permutation matrix, :math:`L` is
  a lower triangular matrix with 1's on the diagonal, and :math:`U` is
  an upper triangular matrix.  This factorization is stored in-place
  on the input SUNMATRIX_DENSE object :math:`A`, with pivoting
  information encoding :math:`P` stored in the ``pivots`` array.

* The "solve" call performs pivoting and forward and
  backward substitution using the stored ``pivots`` array and the
  :math:`LU` factors held in the SUNMATRIX_DENSE object
  (:math:`\mathcal O(N^2)` cost).


The SUNLinSol_Dense module defines dense implementations of all
"direct" linear solver operations listed in
:numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_Dense``

* ``SUNLinSolInitialize_Dense`` -- this does nothing, since all
  consistency checks are performed at solver creation.

* ``SUNLinSolSetup_Dense`` -- this performs the :math:`LU` factorization.

* ``SUNLinSolSolve_Dense`` -- this uses the :math:`LU` factors
  and ``pivots`` array to perform the solve.

* ``SUNLinSolLastFlag_Dense``

* ``SUNLinSolSpace_Dense`` -- this only returns information for
  the storage *within* the solver object, i.e. storage
  for ``N``, ``last_flag``, and ``pivots``.

* ``SUNLinSolFree_Dense``
