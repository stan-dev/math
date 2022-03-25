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

.. _SUNLinSol_LapackBand:

The SUNLinSol_LapackBand Module
======================================

The SUNLinSol_LapackBand implementation of the ``SUNLinearSolver`` class
is designed to be used with the corresponding SUNMATRIX_BAND matrix type,
and one of the serial or shared-memory ``N_Vector`` implementations
(NVECTOR_SERIAL, NVECTOR_OPENMP, or NVECTOR_PTHREADS).  The

.. _SUNLinSol_LapackBand.Usage:

SUNLinSol_LapackBand Usage
-------------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_lapackband.h``.  The installed module
library to link to is
``libsundials_sunlinsollapackband`` *.lib*
where *.lib* is typically ``.so`` for shared libraries and
``.a`` for static libraries.

The module SUNLinSol_LapackBand provides the following
user-callable routine:


.. c:function:: SUNLinearSolver SUNLinSol_LapackBand(N_Vector y, SUNMatrix A, SUNContext sunctx)

   This function creates and allocates memory for a LAPACK band
   ``SUNLinearSolver``.

   **Arguments:**
      * *y* -- vector used to determine the linear system size.
      * *A* -- matrix used to assess compatibility.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      New SUNLinSol_LapackBand object, or ``NULL`` if either
      ``A`` or ``y`` are incompatible.

   **Notes:**
      This routine will perform consistency checks to ensure that it is
      called with consistent ``N_Vector`` and ``SUNMatrix`` implementations.
      These are currently limited to the SUNMATRIX_BAND matrix type and
      the NVECTOR_SERIAL, NVECTOR_OPENMP, and NVECTOR_PTHREADS vector
      types.  As additional compatible matrix and vector implementations
      are added to SUNDIALS, these will be included within this
      compatibility check.

      Additionally, this routine will verify that the input matrix ``A``
      is allocated with appropriate upper bandwidth storage for the
      :math:`LU` factorization.


For backwards compatibility, we also provide the following wrapper function:

.. c:function:: SUNLinearSolver SUNLapackBand(N_Vector y, SUNMatrix A)

   Wrapper function for :c:func:`SUNLinSol_LapackBand`, with
   identical input and output arguments.




.. _SUNLinSol_LapackBand.Description:

SUNLinSol_LapackBand Description
-----------------------------------


SUNLinSol_LapackBand module defines the
*content* field of a ``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_Band {
     sunindextype N;
     sunindextype *pivots;
     sunindextype last_flag;
   };

These entries of the *content* field contain the following
information:

* ``N`` - size of the linear system,

* ``pivots`` - index array for partial pivoting in LU
  factorization,

* ``last_flag`` - last error return flag from internal function
  evaluations.


The SUNLinSol_LapackBand module is a ``SUNLinearSolver`` wrapper for
the LAPACK band matrix factorization and solve routines, ``*GBTRF``
and ``*GBTRS``, where ``*`` is either ``D`` or ``S``, depending on
whether SUNDIALS was configured to have :c:type:`realtype` set to
``double`` or ``single``, respectively (see
:numref:`Usage.CC.DataTypes` for details).
In order to use the SUNLinSol_LapackBand module it is assumed
that LAPACK has been installed on the system prior to installation of
SUNDIALS, and that SUNDIALS has been configured appropriately to
link with LAPACK (see
:numref:`Installation.CMake.ExternalLibraries` for details).  We note
that since there do not exist 128-bit floating-point factorization and
solve routines in LAPACK, this interface cannot be compiled when
using ``extended`` precision for :c:type:`realtype`.  Similarly, since
there do not exist 64-bit integer LAPACK routines, the
SUNLinSol_LapackBand module also cannot be compiled when using
``int64_t`` for the :c:type:`sunindextype`.

This solver is constructed to perform the following operations:

* The "setup" call performs an :math:`LU` factorization with
  partial (row) pivoting, :math:`PA=LU`, where :math:`P` is a
  permutation matrix, :math:`L` is a lower triangular matrix with 1's
  on the diagonal, and :math:`U` is an upper triangular matrix.  This
  factorization is stored in-place on the input SUNMATRIX_BAND object
  :math:`A`, with pivoting information encoding :math:`P` stored in
  the ``pivots`` array.

* The "solve" call performs pivoting and forward and
  backward substitution using the stored ``pivots`` array and the
  :math:`LU` factors held in the SUNMATRIX_BAND object.

* :math:`A` must be allocated to accommodate the increase in upper
  bandwidth that occurs during factorization.  More precisely, if
  :math:`A` is a band matrix with upper bandwidth ``mu`` and lower
  bandwidth ``ml``, then the upper triangular factor :math:`U` can
  have upper bandwidth as big as ``smu = MIN(N-1,mu+ml)``. The lower
  triangular factor :math:`L` has lower bandwidth ``ml``.

The SUNLinSol_LapackBand module defines band implementations of all
"direct" linear solver operations listed in
:numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_LapackBand``

* ``SUNLinSolInitialize_LapackBand`` -- this does nothing, since all
  consistency checks are performed at solver creation.

* ``SUNLinSolSetup_LapackBand`` -- this calls either
  ``DGBTRF`` or ``SGBTRF`` to perform the :math:`LU` factorization.

* ``SUNLinSolSolve_LapackBand`` -- this calls either
  ``DGBTRS`` or ``SGBTRS`` to use the :math:`LU` factors and
  ``pivots`` array to perform the solve.

* ``SUNLinSolLastFlag_LapackBand``

* ``SUNLinSolSpace_LapackBand`` -- this only returns information for
  the storage *within* the solver object, i.e. storage for ``N``,
  ``last_flag``, and ``pivots``.

* ``SUNLinSolFree_LapackBand``
