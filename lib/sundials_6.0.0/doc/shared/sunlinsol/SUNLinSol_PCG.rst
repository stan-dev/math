..
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNLinSol.PCG:

The SUNLinSol_PCG Module
======================================

The SUNLinSol_PCG implementation of the ``SUNLinearSolver`` class performs
the PCG (Preconditioned Conjugate Gradient :cite:p:`HeSt:52`) method; this
is an iterative linear solver that is designed to be compatible with
any ``N_Vector`` implementation that supports a minimal subset of operations
(:c:func:`N_VClone()`, :c:func:`N_VDotProd()`, :c:func:`N_VScale()`,
:c:func:`N_VLinearSum()`, :c:func:`N_VProd()`, and
:c:func:`N_VDestroy()`).  Unlike the SPGMR and SPFGMR algorithms, PCG
requires a fixed amount of memory that does not increase with the
number of allowed iterations.

Unlike all of the other iterative linear solvers supplied with
SUNDIALS, PCG should only be used on *symmetric* linear
systems (e.g. mass matrix linear systems encountered in
ARKODE). As a result, the explanation of the role of scaling and
preconditioning matrices given in general must be modified in this
scenario.  The PCG algorithm solves a linear system :math:`Ax = b` where
:math:`A` is a symmetric (:math:`A^T=A`), real-valued matrix.  Preconditioning is
allowed, and is applied in a symmetric fashion on both the right and
left.  Scaling is also allowed and is applied symmetrically.  We
denote the preconditioner and scaling matrices as follows:

* :math:`P` is the preconditioner (assumed symmetric),

* :math:`S` is a diagonal matrix of scale factors.

The matrices :math:`A` and :math:`P` are not required explicitly; only routines
that provide :math:`A` and :math:`P^{-1}` as operators are required.  The diagonal
of the matrix :math:`S` is held in a single ``N_Vector``, supplied by the user.

In this notation, PCG applies the underlying CG algorithm to the
equivalent transformed system

.. math::
   :label: eq:transformed_linear_systemPCG

   \tilde{A} \tilde{x} = \tilde{b}

where

.. math::
   :label: eq:transformed_linear_system_componentsPCG

   \tilde{A} &= S P^{-1} A P^{-1} S,\\
   \tilde{b} &= S P^{-1} b,\\
   \tilde{x} &= S^{-1} P x.

The scaling matrix must be chosen so that the vectors :math:`SP^{-1}b` and
:math:`S^{-1}Px` have dimensionless components.

The stopping test for the PCG iterations is on the L2 norm of the
scaled preconditioned residual:

.. math::

   &\| \tilde{b} - \tilde{A} \tilde{x} \|_2  <  \delta\\
   \Leftrightarrow\quad &\\
   &\| S P^{-1} b - S P^{-1} A x \|_2  <  \delta\\
   \Leftrightarrow\quad &\\
   &\| P^{-1} b - P^{-1} A x \|_S  <  \delta

where :math:`\| v \|_S = \sqrt{v^T S^T S v}`, with an input tolerance
:math:`\delta`.

.. _SUNLinSol.PCG.Usage:

SUNLinSol_PCG Usage
---------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_pcg.h``.  The SUNLinSol_PCG module
is accessible from all SUNDIALS solvers *without*
linking to the ``libsundials_sunlinsolpcg`` module library.

The module SUNLinSol_PCG provides the following user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_PCG(N_Vector y, int pretype, int maxl, SUNContext sunctx)

   This constructor function creates and allocates memory for a PCG
   ``SUNLinearSolver``.

   **Arguments:**
      * *y* -- a template vector.
      * *pretype* -- a flag indicating the type of preconditioning to use:

        * ``SUN_PREC_NONE``
        * ``SUN_PREC_LEFT``
        * ``SUN_PREC_RIGHT``
        * ``SUN_PREC_BOTH``

      * *maxl* -- the maximum number of linear iterations to allow.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNLinearSolver`` object.  If either *y* is
      incompatible then this routine will return ``NULL``.

   **Notes:**
      This routine will perform consistency checks to ensure that it is
      called with a consistent ``N_Vector`` implementation (i.e. that it
      supplies the requisite vector operations).

      A ``maxl`` argument that is :math:`\le0` will result in the default
      value (5).

      Since the PCG algorithm is designed to only support symmetric
      preconditioning, then any of the ``pretype`` inputs ``SUN_PREC_LEFT``,
      ``SUN_PREC_RIGHT``, or ``SUN_PREC_BOTH`` will result in use
      of the symmetric preconditioner;  any other integer input will
      result in the default (no preconditioning).  Although some SUNDIALS
      solvers are designed to only work with left preconditioning (IDA
      and IDAS) and others with only right preconditioning (KINSOL), PCG
      should *only* be used with these packages when the linear systems
      are known to be *symmetric*.  Since the scaling of matrix rows and
      columns must be identical in a symmetric matrix, symmetric
      preconditioning should work appropriately even for packages
      designed with one-sided preconditioning in mind.

.. c:function:: int SUNLinSol_PCGSetPrecType(SUNLinearSolver S, int pretype)

   This function updates the flag indicating use of preconditioning.

   **Arguments:**
      * *S* -- SUNLinSol_PCG object to update.
      * *pretype* -- a flag indicating the type of preconditioning to use:

        * ``SUN_PREC_NONE``
        * ``SUN_PREC_LEFT``
        * ``SUN_PREC_RIGHT``
        * ``SUN_PREC_BOTH``

   **Return value:**
      * ``SUNLS_SUCCESS`` -- successful update.
      * ``SUNLS_ILL_INPUT`` -- illegal ``pretype``
      * ``SUNLS_MEM_NULL`` -- ``S`` is ``NULL``

   **Notes:**
      As above, any one of the input values, ``SUN_PREC_LEFT``,
      ``SUN_PREC_RIGHT``, or ``SUN_PREC_BOTH`` will enable
      preconditioning; ``SUN_PREC_NONE`` disables preconditioning.



.. c:function:: int SUNLinSol_PCGSetMaxl(SUNLinearSolver S, int maxl)

   This function updates the number of linear solver iterations to allow.

   **Arguments:**
      * *S* -- SUNLinSol_PCG object to update.
      * *maxl* -- maximum number of linear iterations to allow.  Any
        non-positive input will result in the default value (5).

   **Return value:**
      * ``SUNLS_SUCCESS`` -- successful update.
      * ``SUNLS_MEM_NULL`` -- ``S`` is ``NULL``



.. c:function:: int SUNLinSolSetInfoFile_PCG(SUNLinearSolver LS, FILE* info_file)

   The function :c:func:`SUNLinSolSetInfoFile_PCG()` sets the
   output file where all informative (non-error) messages should be directed.

   **Arguments:**
      * *LS* -- a SUNLinSol object
      * *info_file* -- pointer to output file (``stdout`` by default);
         a ``NULL`` input will disable output

   **Return value:**
      * *SUNLS_SUCCESS* if successful
      * *SUNLS_MEM_NULL* if the SUNLinearSolver memory was ``NULL``
      * *SUNLS_ILL_INPUT* if SUNDIALS was not built with monitoring enabled

   **Notes:**
      This function is intended for users that wish to monitor the linear
      solver progress. By default, the file pointer is set to ``stdout``.

      **SUNDIALS must be built with the CMake option**
      ``SUNDIALS_BUILD_WITH_MONITORING`` **to utilize this function.**
      See :numref:`Installation.CMake.Options` for more information.


.. c:function:: int SUNLinSolSetPrintLevel_PCG(SUNLinearSolver LS, int print_level)

   The function :c:func:`SUNLinSolSetPrintLevel_PCG()` specifies the
   level of verbosity of the output.

   **Arguments:**
      * *LS* -- a SUNLinSol object
      * *print_level* -- flag indicating level of verbosity;
        must be one of:

         * 0, no information is printed (default)
         * 1, for each linear iteration the residual norm is printed

   **Return value:**
      * *SUNLS_SUCCESS* if successful
      * *SUNLS_MEM_NULL* if the SUNLinearSolver memory was ``NULL``
      * *SUNLS_ILL_INPUT* if SUNDIALS was not built with monitoring enabled, or
        if the print level value was invalid

   **Notes:**
      This function is intended for users that wish to monitor the linear
      solver progress. By default, the print level is 0.

      **SUNDIALS must be built with the CMake option**
      ``SUNDIALS_BUILD_WITH_MONITORING`` **to utilize this function.**
      See :numref:`Installation.CMake.Options` for more information.

For backwards compatibility, we also provide the following wrapper functions,
each with identical input and output arguments to the routines that
they wrap:

.. c:function:: SUNLinearSolver SUNPCG(N_Vector y, int pretype, int maxl)

   Wrapper function for :c:func:`SUNLinSol_PCG`

.. c:function:: int SUNPCGSetPrecType(SUNLinearSolver S, int pretype)

   Wrapper function for :c:func:`SUNLinSol_PCGSetPrecType()`

.. c:function:: int SUNPCGSetMaxl(SUNLinearSolver S, int maxl)

   Wrapper function for :c:func:`SUNLinSol_PCGSetMaxl()`



.. _SUNLinSol.PCG.Description:

SUNLinSol_PCG Description
---------------------------


The SUNLinSol_PCG module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_PCG {
     int maxl;
     int pretype;
     booleantype zeroguess;
     int numiters;
     realtype resnorm;
     int last_flag;
     SUNATimesFn ATimes;
     void* ATData;
     SUNPSetupFn Psetup;
     SUNPSolveFn Psolve;
     void* PData;
     N_Vector s;
     N_Vector r;
     N_Vector p;
     N_Vector z;
     N_Vector Ap;
     int      print_level;
     FILE*    info_file;
   };

These entries of the *content* field contain the following
information:

* ``maxl`` - number of PCG iterations to allow (default is 5),

* ``pretype`` - flag for use of preconditioning (default is none),

* ``numiters`` - number of iterations from the most-recent solve,

* ``resnorm`` - final linear residual norm from the most-recent
  solve,

* ``last_flag`` - last error return flag from an internal
  function,

* ``ATimes`` - function pointer to perform :math:`Av` product,

* ``ATData`` - pointer to structure for ``ATimes``,

* ``Psetup`` - function pointer to preconditioner setup routine,

* ``Psolve`` - function pointer to preconditioner solve routine,

* ``PData`` - pointer to structure for ``Psetup`` and ``Psolve``,

* ``s`` - vector pointer for supplied scaling matrix
  (default is ``NULL``),

* ``r`` - a ``N_Vector`` which holds the preconditioned linear system
  residual,

* ``p, z, Ap`` - ``N_Vector`` used for workspace by the
  PCG algorithm.

* ``print_level`` - controls the amount of information to be printed to the info file

* ``info_file``   - the file where all informative (non-error) messages will be directed


This solver is constructed to perform the following operations:

* During construction all ``N_Vector`` solver data is allocated, with
  vectors cloned from a template ``N_Vector`` that is input, and
  default solver parameters are set.

* User-facing "set" routines may be called to modify default
  solver parameters.

* Additional "set" routines are called by the SUNDIALS solver
  that interfaces with SUNLinSol_PCG to supply the
  ``ATimes``, ``PSetup``, and ``Psolve`` function pointers and
  ``s`` scaling vector.

* In the "initialize" call, the solver parameters are checked
  for validity.

* In the "setup" call, any non-``NULL`` ``PSetup`` function is
  called.  Typically, this is provided by the SUNDIALS solver
  itself, that translates between the generic ``PSetup`` function and
  the solver-specific routine (solver-supplied or user-supplied).

* In the "solve" call the PCG iteration is performed.  This
  will include scaling and preconditioning if those options have been
  supplied.

The SUNLinSol_PCG module defines implementations of all
"iterative" linear solver operations listed in
:numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_PCG``

* ``SUNLinSolInitialize_PCG``

* ``SUNLinSolSetATimes_PCG``

* ``SUNLinSolSetPreconditioner_PCG``

* ``SUNLinSolSetScalingVectors_PCG`` -- since PCG only supports
  symmetric scaling, the second ``N_Vector`` argument to this function
  is ignored.

* ``SUNLinSolSetZeroGuess_PCG`` -- note the solver assumes a non-zero guess by
  default and the zero guess flag is reset to ``SUNFALSE`` after each call to
  :c:func:`SUNLinSolSolve_PCG`.

* ``SUNLinSolSetup_PCG``

* ``SUNLinSolSolve_PCG``

* ``SUNLinSolNumIters_PCG``

* ``SUNLinSolResNorm_PCG``

* ``SUNLinSolResid_PCG``

* ``SUNLinSolLastFlag_PCG``

* ``SUNLinSolSpace_PCG``

* ``SUNLinSolFree_PCG``
