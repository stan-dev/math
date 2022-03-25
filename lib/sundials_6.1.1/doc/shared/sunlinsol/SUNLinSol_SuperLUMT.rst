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

.. _SUNLinSol.SuperLUMT:

The SUNLinSol_SuperLUMT Module
======================================

The SUNLinSol_SuperLUMT implementation of the ``SUNLinearSolver`` class
interfaces with the SuperLU_MT library.  This is designed to be used
with the corresponding SUNMATRIX_SPARSE matrix type, and one of the
serial or shared-memory ``N_Vector`` implementations (NVECTOR_SERIAL,
NVECTOR_OPENMP, or NVECTOR_PTHREADS).  While these are compatible, it
is not recommended to use a threaded vector module with
SUNLinSol_SuperLUMT unless it is the NVECTOR_OPENMP module and the
SuperLU_MT library has also been compiled with OpenMP.


.. _SUNLinSol.SuperLUMT.Usage:

SUNLinSol_SuperLUMT Usage
-----------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol.SuperLUMT.h``.  The installed module
library to link to is ``libsundials_sunlinsolsuperlumt`` *.lib*
where *.lib* is typically ``.so`` for shared libraries and
``.a`` for static libraries.

The module SUNLinSol_SuperLUMT provides the following user-callable routines:

.. c:function:: SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y, SUNMatrix A, int num_threads, SUNContext sunctx)

   This constructor function creates and allocates memory for a
   SUNLinSol_SuperLUMT object.

   **Arguments:**
      * *y* -- a template vector.
      * *A* -- a template matrix
      * *num_threads* -- desired number of threads (OpenMP or Pthreads,
        depending on how SuperLU_MT was installed) to use during the
        factorization steps.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNLinearSolver`` object; otherwise this routine will return ``NULL``.

   **Notes:**
      This routine analyzes the input matrix and vector to determine the linear
      system size and to assess compatibility with the SuperLU_MT library.

      This routine will perform consistency checks to ensure that it is
      called with consistent ``N_Vector`` and ``SUNMatrix``
      implementations.  These are currently limited to the
      SUNMATRIX_SPARSE matrix type (using either CSR or CSC storage
      formats) and the NVECTOR_SERIAL, NVECTOR_OPENMP, and
      NVECTOR_PTHREADS vector types.  As additional compatible matrix and
      vector implementations are added to SUNDIALS, these will be
      included within this compatibility check.

      The ``num_threads`` argument is not checked
      and is passed directly to SuperLU_MT routines.


.. c:function:: int SUNLinSol_SuperLUMTSetOrdering(SUNLinearSolver S, int ordering_choice)

   This function sets the ordering used by SuperLU_MT for reducing fill in
   the linear solve.

   **Arguments:**
      * *S* -- the SUNLinSol_SuperLUMT object to update.
      * *ordering_choice*:

        0. natural ordering

        1. minimal degree ordering on :math:`A^TA`

        2. minimal degree ordering on :math:`A^T+A`

        3. COLAMD ordering for unsymmetric matrices

      The default is 3 for COLAMD.

   **Return value:**
      * ``SUNLS_SUCCESS`` -- option successfully set
      * ``SUNLS_MEM_NULL`` -- ``S`` is ``NULL``
      * ``SUNLS_ILL_INPUT`` -- invalid ``ordering_choice``


For backwards compatibility, we also provide the following wrapper functions,
each with identical input and output arguments to the routines that
they wrap:

.. c:function:: SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A, int num_threads)

   Wrapper for :c:func:`SUNLinSol_SuperLUMT`.

and

.. c:function:: int SUNSuperLUMTSetOrdering(SUNLinearSolver S, int ordering_choice)

   Wrapper for :c:func:`SUNLinSol_SuperLUMTSetOrdering()`.




.. _SUNLinSol.SuperLUMT.Description:

SUNLinSol_SuperLUMT Description
----------------------------------

The SUNLinSol_SuperLUMT module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_SuperLUMT {
     int          last_flag;
     int          first_factorize;
     SuperMatrix  *A, *AC, *L, *U, *B;
     Gstat_t      *Gstat;
     sunindextype *perm_r, *perm_c;
     sunindextype N;
     int          num_threads;
     realtype     diag_pivot_thresh;
     int          ordering;
     superlumt_options_t *options;
   };

These entries of the *content* field contain the following
information:

* ``last_flag`` - last error return flag from internal function
  evaluations,

* ``first_factorize`` - flag indicating whether the factorization
  has ever been performed,

* ``A, AC, L, U, B`` - ``SuperMatrix`` pointers used in solve,

* ``Gstat`` - ``GStat_t`` object used in solve,

* ``perm_r, perm_c`` - permutation arrays used in solve,

* ``N`` - size of the linear system,

* ``num_threads`` - number of OpenMP/Pthreads threads to use,

* ``diag_pivot_thresh`` - threshold on diagonal pivoting,

* ``ordering`` - flag for which reordering algorithm to use,

* ``options`` - pointer to SuperLU_MT options structure.

The SUNLinSol_SuperLUMT module is a ``SUNLinearSolver`` wrapper for
the SuperLU_MT sparse matrix factorization and solver library
written by X. Sherry Li and collaborators
:cite:p:`SuperLUMT_site,Li:05,DGL:99`.  The
package performs matrix factorization using threads to enhance
efficiency in shared memory parallel environments.  It should be noted
that threads are only used in the factorization step.  In
order to use the SUNLinSol_SuperLUMT interface to SuperLU_MT, it is
assumed that SuperLU_MT has been installed on the system prior to
installation of SUNDIALS, and that SUNDIALS has been configured
appropriately to link with SuperLU_MT (see
:numref:`Installation.CMake.ExternalLibraries` for details).
Additionally, this wrapper only supports single- and
double-precision calculations, and therefore cannot be compiled if
SUNDIALS is configured to have :c:type:`realtype` set to ``extended``
(see :numref:`Usage.CC.DataTypes` for details).  Moreover,
since the SuperLU_MT library may be installed to support either 32-bit
or 64-bit integers, it is assumed that the SuperLU_MT library is
installed using the same integer precision as the SUNDIALS
:c:type:`sunindextype` option.

The SuperLU_MT library has a symbolic factorization routine that
computes the permutation of the linear system matrix to reduce fill-in
on subsequent :math:`LU` factorizations (using COLAMD, minimal degree
ordering on :math:`A^T*A`, minimal degree ordering on :math:`A^T+A`,
or natural ordering).  Of these ordering choices, the default value in
the SUNLinSol_SuperLUMT module is the COLAMD ordering.

Since the linear systems that arise within the context of SUNDIALS
calculations will typically have identical sparsity patterns, the
SUNLinSol_SuperLUMT module is constructed to perform the
following operations:

* The first time that the "setup" routine is called, it
  performs the symbolic factorization, followed by an initial
  numerical factorization.

* On subsequent calls to the "setup" routine, it skips the
  symbolic factorization, and only refactors the input matrix.

* The "solve" call performs pivoting and forward and
  backward substitution using the stored SuperLU_MT data
  structures.  We note that in this solve SuperLU_MT operates on the
  native data arrays for the right-hand side and solution vectors,
  without requiring costly data copies.


The SUNLinSol_SuperLUMT module defines implementations of all
"direct" linear solver operations listed in
:numref:`SUNLinSol.API`:


* ``SUNLinSolGetType_SuperLUMT``

* ``SUNLinSolInitialize_SuperLUMT`` -- this sets the
  ``first_factorize`` flag to 1 and resets the internal SuperLU_MT
  statistics variables.

* ``SUNLinSolSetup_SuperLUMT`` -- this performs either a :math:`LU`
  factorization or refactorization of the input matrix.

* ``SUNLinSolSolve_SuperLUMT`` -- this calls the appropriate
  SuperLU_MT solve routine to utilize the :math:`LU` factors to solve the
  linear system.

* ``SUNLinSolLastFlag_SuperLUMT``

* ``SUNLinSolSpace_SuperLUMT`` -- this only returns information for
  the storage within the solver *interface*, i.e. storage for the
  integers ``last_flag`` and ``first_factorize``.  For additional
  space requirements, see the SuperLU_MT documentation.

* ``SUNLinSolFree_SuperLUMT``
