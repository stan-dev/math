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

.. _SUNLinSol.KLU:

The SUNLinSol_KLU Module
======================================

The SUNLinSol_KLU implementation of the ``SUNLinearSolver`` class
is designed to be used with the corresponding SUNMATRIX_SPARSE matrix type,
and one of the serial or shared-memory ``N_Vector`` implementations
(NVECTOR_SERIAL, NVECTOR_OPENMP, or NVECTOR_PTHREADS).

.. _SUNLinSol.KLU.Usage:

SUNLinSol_KLU Usage
------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_klu.h``.  The installed module
library to link to is ``libsundials_sunlinsolklu`` *.lib*
where *.lib* is typically ``.so`` for shared libraries and
``.a`` for static libraries.

The module SUNLinSol_KLU provides the following additional
user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_KLU(N_Vector y, SUNMatrix A, SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNLinSol_KLU
   object.

   **Arguments:**
      * *y* -- vector used to determine the linear system size.
      * *A* -- matrix used to assess compatibility.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      New SUNLinSol_KLU object, or ``NULL`` if either ``A`` or ``y`` are incompatible.

   **Notes:**
      This routine will perform consistency checks to ensure that it is
      called with consistent ``N_Vector`` and ``SUNMatrix`` implementations.
      These are currently limited to the SUNMATRIX_SPARSE matrix type
      (using either CSR or CSC storage formats) and the NVECTOR_SERIAL,
      NVECTOR_OPENMP, and NVECTOR_PTHREADS vector types.  As additional
      compatible matrix and vector implementations are added to
      SUNDIALS, these will be included within this compatibility
      check.


.. c:function:: int SUNLinSol_KLUReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz, int reinit_type)

   This function reinitializes memory and flags for a new factorization
   (symbolic and numeric) to be conducted at the next solver setup
   call.  This routine is useful in the cases where the number of
   nonzeroes has changed or if the structure of the linear system has
   changed which would require a new symbolic (and numeric
   factorization).

   **Arguments:**
      * *S* -- existing SUNLinSol_KLU object to reinitialize.
      * *A* -- sparse ``SUNMatrix`` matrix (with updated structure)
        to use for reinitialization.
      * *nnz* -- maximum number of nonzeros expected for Jacobian matrix.
      * *reinit_type* -- governs the level of reinitialization.  The allowed values are:

         1. The Jacobian matrix will be destroyed and a new one will be
            allocated based on the ``nnz`` value passed to this call.  New
            symbolic and numeric factorizations will be completed at the next
            solver setup.

         2. Only symbolic and numeric factorizations will be completed.
            It is assumed that the Jacobian size has not exceeded the size of
            ``nnz`` given in the sparse matrix provided to the original
            constructor routine (or the previous ``SUNKLUReInit`` call).

   **Return value:**
      * ``SUNLS_SUCCESS`` -- reinitialization successful.
      * ``SUNLS_MEM_NULL`` -- either ``S`` or ``A`` are ``NULL``.
      * ``SUNLS_ILL_INPUT`` -- ``A`` does not have type ``SUNMATRIX_SPARSE`` or
         ``reinit_type`` is invalid.
      * ``SUNLS_MEM_FAIL`` reallocation of the sparse matrix failed.

   **Notes:**
      This routine assumes no other changes to solver use are necessary.


.. c:function:: int SUNLinSol_KLUSetOrdering(SUNLinearSolver S, int ordering_choice)

   This function sets the ordering used by KLU for reducing fill in
   the linear solve.

   **Arguments:**
      * *S* -- existing SUNLinSol_KLU object to update.
      * *ordering_choice* -- type of ordering to use, options are:

         0. AMD,

         1. COLAMD, and

         2. the natural ordering.

         The default is 1 for COLAMD.

   **Return value:**
      * ``SUNLS_SUCCESS`` -- ordering choice successfully updated.
      * ``SUNLS_MEM_NULL`` -- ``S`` is ``NULL``.
      * ``SUNLS_ILL_INPUT`` -- ``ordering_choice``.


.. c:function:: sun_klu_symbolic* SUNLinSol_KLUGetSymbolic(SUNLinearSolver S)

   This function returns a pointer to the KLU symbolic factorization
   stored in the SUNLinSol_KLU ``content`` structure.

   When SUNDIALS is compiled with 32-bit indices (``SUNDIALS_INDEX_SIZE=32``),
   ``sun_klu_symbolic`` is mapped to the KLU type ``klu_symbolic``; when
   SUNDIALS compiled with 64-bit indices (``SUNDIALS_INDEX_SIZE=64``) this is
   mapped to the KLU type ``klu_l_symbolic``.


.. c:function:: sun_klu_numeric* SUNLinSol_KLUGetNumeric(SUNLinearSolver S)

   This function returns a pointer to the KLU numeric factorization
   stored in the SUNLinSol_KLU ``content`` structure.

   When SUNDIALS is compiled with 32-bit indices (``SUNDIALS_INDEX_SIZE=32``),
   ``sun_klu_numeric`` is mapped to the KLU type ``klu_numeric``; when
   SUNDIALS is compiled with 64-bit indices (``SUNDIALS_INDEX_SIZE=64``) this is
   mapped to the KLU type ``klu_l_numeric``.


.. c:function:: sun_klu_common* SUNLinSol_KLUGetCommon(SUNLinearSolver S)

   This function returns a pointer to the KLU common structure
   stored in the SUNLinSol_KLU ``content`` structure.

   When SUNDIALS is compiled with 32-bit indices (``SUNDIALS_INDEX_SIZE=32``),
   ``sun_klu_common`` is mapped to the KLU type ``klu_common``; when
   SUNDIALS is compiled with 64-bit indices  (``SUNDIALS_INDEX_SIZE=64``) this is
   mapped to the KLU type ``klu_l_common``.


For backwards compatibility, we also provide the following wrapper functions,
each with identical input and output arguments to the routines that
they wrap:

.. c:function:: SUNLinearSolver SUNKLU(N_Vector y, SUNMatrix A)

   Wrapper function for :c:func:`SUNLinSol_KLU`

.. c:function:: int SUNKLUReInit(SUNLinearSolver S, SUNMatrix A, sunindextype nnz, int reinit_type)

   Wrapper function for :c:func:`SUNLinSol_KLUReInit()`

.. c:function:: int SUNKLUSetOrdering(SUNLinearSolver S, int ordering_choice)

   Wrapper function for :c:func:`SUNLinSol_KLUSetOrdering()`




.. _SUNLinSol.KLU.Description:

SUNLinSol_KLU Description
--------------------------


The SUNLinSol_KLU module defines the *content*
field of a ``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_KLU {
     int              last_flag;
     int              first_factorize;
     sun_klu_symbolic *symbolic;
     sun_klu_numeric  *numeric;
     sun_klu_common   common;
     sunindextype     (*klu_solver)(sun_klu_symbolic*, sun_klu_numeric*,
                                    sunindextype, sunindextype,
                                    double*, sun_klu_common*);
   };

These entries of the *content* field contain the following
information:

* ``last_flag`` - last error return flag from internal function
  evaluations,

* ``first_factorize`` - flag indicating whether the factorization
  has ever been performed,

* ``symbolic`` - KLU storage structure for symbolic
  factorization components, with underlying type ``klu_symbolic``
  or ``klu_l_symbolic``, depending on whether SUNDIALS was
  installed with 32-bit versus 64-bit indices, respectively,

* ``numeric`` - KLU storage structure for numeric factorization
  components, with underlying type ``klu_numeric``
  or ``klu_l_numeric``, depending on whether SUNDIALS was
  installed with 32-bit versus 64-bit indices, respectively,

* ``common`` - storage structure for common KLU solver
  components, with underlying type ``klu_common``
  or ``klu_l_common``, depending on whether SUNDIALS was
  installed with 32-bit versus 64-bit indices, respectively,

* ``klu_solver`` -- pointer to the appropriate KLU solver function
  (depending on whether it is using a CSR or CSC sparse matrix, and
  on whether SUNDIALS was installed with 32-bit or 64-bit indices).


The SUNLinSol_KLU module is a ``SUNLinearSolver`` wrapper for
the KLU sparse matrix factorization and solver library written by Tim
Davis and collaborators (:cite:p:`KLU_site,DaPa:10`).  In order to use the
SUNLinSol_KLU interface to KLU, it is assumed that KLU has
been installed on the system prior to installation of SUNDIALS, and
that SUNDIALS has been configured appropriately to link with KLU
(see :numref:`Installation.CMake.ExternalLibraries` for details).
Additionally, this wrapper only supports double-precision
calculations, and therefore cannot be compiled if SUNDIALS is
configured to have :c:type:`realtype` set to either ``extended`` or
``single`` (see :ref:`Usage.CC.DataTypes` for
details). Since the KLU library supports both 32-bit and 64-bit
integers, this interface will be compiled for either of the available
:c:type:`sunindextype` options.

The KLU library has a symbolic factorization routine that computes
the permutation of the linear system matrix to block triangular form
and the permutations that will pre-order the diagonal blocks (the only
ones that need to be factored) to reduce fill-in (using AMD, COLAMD,
CHOLAMD, natural, or an ordering given by the user).  Of these
ordering choices, the default value in the SUNLinSol_KLU
module is the COLAMD ordering.

KLU breaks the factorization into two separate parts.  The first is
a symbolic factorization and the second is a numeric factorization
that returns the factored matrix along with final pivot information.
KLU also has a refactor routine that can be called instead of the numeric
factorization.  This routine will reuse the pivot information.  This routine
also returns diagnostic information that a user can examine to determine if
numerical stability is being lost and a full numerical factorization should
be done instead of the refactor.

Since the linear systems that arise within the context of SUNDIALS
calculations will typically have identical sparsity patterns, the
SUNLinSol_KLU module is constructed to perform the
following operations:

* The first time that the "setup" routine is called, it
  performs the symbolic factorization, followed by an initial
  numerical factorization.

* On subsequent calls to the "setup" routine, it calls the
  appropriate KLU "refactor" routine, followed by estimates of
  the numerical conditioning using the relevant "rcond", and if
  necessary "condest", routine(s).  If these estimates of the
  condition number are larger than :math:`\varepsilon^{-2/3}` (where
  :math:`\varepsilon` is the double-precision unit roundoff), then a new
  factorization is performed.

* The module includes the routine ``SUNKLUReInit``, that
  can be called by the user to force a full refactorization at the
  next "setup" call.

* The "solve" call performs pivoting and forward and
  backward substitution using the stored KLU data structures.  We
  note that in this solve KLU operates on the native data arrays
  for the right-hand side and solution vectors, without requiring
  costly data copies.


The SUNLinSol_KLU module defines implementations of all
"direct" linear solver operations listed in
:numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_KLU``

* ``SUNLinSolInitialize_KLU`` -- this sets the
  ``first_factorize`` flag to 1, forcing both symbolic and numerical
  factorizations on the subsequent "setup" call.

* ``SUNLinSolSetup_KLU`` -- this performs either a :math:`LU`
  factorization or refactorization of the input matrix.

* ``SUNLinSolSolve_KLU`` -- this calls the appropriate KLU
  solve routine to utilize the :math:`LU` factors to solve the linear
  system.

* ``SUNLinSolLastFlag_KLU``

* ``SUNLinSolSpace_KLU`` -- this only returns information for
  the storage within the solver *interface*, i.e. storage for the
  integers ``last_flag`` and ``first_factorize``.  For additional
  space requirements, see the KLU documentation.

* ``SUNLinSolFree_KLU``
