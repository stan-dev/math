..
   Programmer(s): Cody J. Balos @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNLinSol.SuperLUDIST:

The SUNLinSol_SuperLUDIST Module
======================================

The SUNLinsol_SuperLUDIST implementation of the ``SUNLinearSolver`` class interfaces
with the SuperLU_DIST library.  This is designed to be used with the
SUNMatrix_SLUNRloc :c:type:`SUNMatrix`, and one of the serial, threaded or parallel
N_Vector implementations (NVECTOR_SERIAL, NVECTOR_OPENMP, NVECTOR_PTHREADS,
NVECTOR_PARALLEL, NVECTOR_PARHYP).


.. _SUNLinSol.SuperLUDIST.Usage:

SUNLinSol_SuperLUDIST Usage
-----------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_superludist.h``.  The installed module
library to link to is ``libsundials_sunlinsolsuperludist`` *.lib*
where *.lib* is typically ``.so`` for shared libraries and
``.a`` for static libraries.

The module SUNLinSol_SuperLUDIST provides the following user-callable routines:

.. warning::

  Starting with SuperLU_DIST version 6.3.0, some structures were
  renamed to have a prefix for the floating point type. The double precision API
  functions have the prefix 'd'. To maintain backwards compatibility with the
  unprefixed types, SUNDIALS provides macros to these SuperLU_DIST types with an
  'x' prefix that expand to the correct prefix. E.g., the SUNDIALS macro
  ``xLUstruct_t`` expands to ``dLUstruct_t`` or ``LUstruct_t`` based on the
  SuperLU_DIST version.


.. c:function:: SUNLinearSolver SUNLinSol_SuperLUDIST(N_Vector y, SuperMatrix *A, gridinfo_t *grid, xLUstruct_t *lu, xScalePermstruct_t *scaleperm, xSOLVEstruct_t *solve, SuperLUStat_t *stat, superlu_dist_options_t *options, SUNContext sunctx)

   This constructor function creates and allocates memory for a SUNLinSol_SuperLUDIST
   object.

   **Arguments:**
      * *y* -- a template vector.
      * *A* -- a template matrix
      * *grid*, *lu*, *scaleperm*, *solve*, *stat*, *options* -- SuperLU_DIST object pointers.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNLinearSolver`` object; otherwise this routine will return ``NULL``.

   **Notes:**
      This routine analyzes the input matrix and vector to determine the linear
      system size and to assess the compatibility with the SuperLU_DIST library.

      This routine will perform consistency checks to ensure that it is called with
      consistent N_Vector and SUNMatrix implementations. These are currently limited
      to the SUNMatrix_SLUNRloc matrix type and the NVECTOR_SERIAL, NVECTOR_OPENMP,
      NVECTOR_PTHREADS, NVECTOR_PARALLEL, and NVECTOR_PARHYP vector types. As
      additional compatible matrix and vector implementations are added to SUNDIALS,
      these will be included within this compatibility check.

      The ``grid``, ``lu``, ``scaleperm``, ``solve``, and ``options`` arguments are
      not checked and are passed directly to SuperLU_DIST routines.

      Some struct members of the ``options`` argument are modified internally by
      the SUNLinSol_SuperLUDIST solver. Specifically, the member ``Fact``
      is modified in the setup and solve routines.


.. c:function:: realtype SUNLinSol_SuperLUDIST_GetBerr(SUNLinearSolver LS)

   This function returns the componentwise relative backward error of the
   computed solution.   It takes one argument, the ``SUNLinearSolver`` object.
   The return type is ``realtype``.


.. c:function:: gridinfo_t* SUNLinSol_SuperLUDIST_GetGridinfo(SUNLinearSolver LS)

   This function returns a pointer to the SuperLU_DIST structure that
   contains the 2D process grid. It takes one argument, the ``SUNLinearSolver``
   object.


.. c:function:: xLUstruct_t* SUNLinSol_SuperLUDIST_GetLUstruct(SUNLinearSolver LS)

   This function returns a pointer to the SuperLU_DIST structure that contains
   the distributed ``L`` and ``U`` structures. It takes one argument, the
   ``SUNLinearSolver`` object.


.. c:function:: superlu_dist_options_t* SUNLinSol_SuperLUDIST_GetSuperLUOptions(SUNLinearSolver LS)

   This function returns a pointer to the SuperLU_DIST structure that contains the
   options which control how the linear system is factorized and solved. It takes
   one argument, the ``SUNLinearSolver`` object.


.. c:function:: xScalePermstruct_t* SUNLinSol_SuperLUDIST_GetScalePermstruct(SUNLinearSolver LS)

   This function returns a pointer to the SuperLU_DIST structure that contains
   the vectors that describe the transformations done to the matrix ``A``. It
   takes one argument, the ``SUNLinearSolver`` object.


.. c:function:: xSOLVEstruct_t* SUNLinSol_SuperLUDIST_GetSOLVEstruct(SUNLinearSolver LS)

   This function returns a pointer to the SuperLU_DIST structure that contains
   information for communication during the solution phase. It takes one argument
   the ``SUNLinearSolver`` object.

.. c:function:: SuperLUStat_t* SUNLinSol_SuperLUDIST_GetSuperLUStat(SUNLinearSolver LS)

   This function returns a pointer to the SuperLU_DIST structure that stores
   information about runtime and flop count. It takes one argument, the
   ``SUNLinearSolver`` object.



.. _SUNLinSol.SuperLUDIST.Description:

SUNLinSol_SuperLUDIST Description
----------------------------------

The SUNLinSol_SuperLUDIST module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_SuperLUDIST {
     booleantype             first_factorize;
     int                     last_flag;
     realtype                berr;
     gridinfo_t              *grid;
     xLUstruct_t             *lu;
     superlu_dist_options_t  *options;
     xScalePermstruct_t      *scaleperm;
     xSOLVEstruct_t          *solve;
     SuperLUStat_t           *stat;
     sunindextype            N;
   };

These entries of the *content* field contain the following
information:

* ``first_factorize`` -- flag indicating whether the factorization
  has ever been performed,

* ``last_flag`` -- last error return flag from internal function
  evaluations,

* ``berr`` -- the componentwise relative backward error of the computed solution,

* ``grid`` -- pointer to the SuperLU_DIST structure that strores the 2D process grid

* ``lu`` -- pointer to the SuperLU_DIST structure that stores the distributed ``L``
  and ``U`` factors,

* ``scaleperm`` -- pointer to the SuperLU_DIST structure that stores vectors describing
  the transformations done to the matrix ``A``,

* ``options`` -- pointer to the SuperLU_DIST stucture which contains options that control
  how the linear system is factorized and solved,

* ``solve`` -- pointer to the SuperLU_DIST solve structure,

* ``stat`` -- pointer to the SuperLU_DIST structure that stores information about runtime
  and flop count,

* ``N`` -- the number of equations in the system.


The SUNLinSol_SuperLUDIST module is a SUNLinearSolver adapter for the
SuperLU_DIST sparse matrix factorization and solver library written by
X. Sherry Li and collaborators :cite:p:`SuperLUDIST_site,GDL:07,LD:03,SLUUG:99`.
The package uses a SPMD parallel programming model and multithreading
to enhance efficiency in distributed-memory parallel environments with
multicore nodes and possibly GPU accelerators. It uses MPI for communication,
OpenMP for threading, and CUDA for GPU support. In order to use the
SUNLinSol_SuperLUDIST interface to SuperLU_DIST, it is assumed that SuperLU_DIST
has been installed on the system prior to installation of SUNDIALS, and
that SUNDIALS has been configured appropriately to link with SuperLU_DIST
(see :numref:`Installation.CMake.ExternalLibraries` for details).
Additionally, the wrapper only
supports double-precision calculations, and therefore cannot be compiled if SUNDIALS
is configured to use single or extended precision. Moreover, since the SuperLU_DIST
library may be installed to support either 32-bit or 64-bit integers,
it is assumed that the SuperLU_DIST library is installed using the same
integer size as SUNDIALS.

The SuperLU_DIST library provides many options to control how a linear
system will be factorized and solved. These options may be set by a user
on an instance of the ``superlu_dist_options_t`` struct, and then it may be provided
as an argument to the SUNLinSol_SuperLUDIST constructor. The SUNLinSol_SuperLUDIST
module will respect all options set except for ``Fact`` -- this option is
necessarily modified by the SUNLinSol_SuperLUDIST module in the setup and solve routines.

Since the linear systems that arise within the context of SUNDIALS calculations will
typically have identical sparsity patterns, the SUNLinSol_SuperLUDIST module is
constructed to perform the following operations:

* The first time that the "setup" routine is called, it
  sets the SuperLU_DIST option ``Fact`` to ``DOFACT`` so that a subsequent
  call to the "solve" routine will perform a symbolic factorization,
  followed by an initial numerical factorization before continuing
  to solve the system.

* On subsequent calls to the "setup" routine, it sets the
  SuperLU_DIST option ``Fact`` to ``SamePattern`` so that
  a subsequent call to "solve" will perform factorization assuming
  the same sparsity pattern as prior, i.e. it will reuse the column
  permutation vector.

* If "setup" is called prior to the "solve" routine, then the "solve" routine
  will perform a symbolic factorization, followed by an initial
  numerical factorization before continuing to the sparse triangular
  solves, and, potentially, iterative refinement. If "setup" is not
  called prior, "solve" will skip to the triangular solve step. We
  note that in this solve SuperLU_DIST operates on the native data arrays
  for the right-hand side and solution vectors, without requiring costly data copies.


The SUNLinSol_SuperLUDIST module defines implementations of all
"direct" linear solver operations listed in
:numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_SuperLUDIST``

* ``SUNLinSolInitialize_SuperLUDIST`` -- this sets the
  ``first_factorize`` flag to 1 and resets the internal SuperLU_DIST
  statistics variables.

* ``SUNLinSolSetup_SuperLUDIST`` -- this sets the appropriate
  SuperLU_DIST options so that a subsequent solve will perform a
  symbolic and numerical factorization before proceeding with the
  triangular solves

* ``SUNLinSolSolve_SuperLUDIST`` -- this calls the SuperLU_DIST
  solve routine to perform factorization (if the setup routine
  was called prior) and then use the $LU$ factors to solve the
  linear system.

* ``SUNLinSolLastFlag_SuperLUDIST``

* ``SUNLinSolSpace_SuperLUDIST`` -- this only returns information for
  the storage within the solver *interface*, i.e. storage for the
  integers ``last_flag`` and ``first_factorize``.  For additional
  space requirements, see the SuperLU_DIST documentation.

* ``SUNLinSolFree_SuperLUDIST``
