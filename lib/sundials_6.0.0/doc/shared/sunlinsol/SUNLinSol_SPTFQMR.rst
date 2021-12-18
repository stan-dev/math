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

.. _SUNLinSol.SPTFQMR:

The SUNLinSol_SPTFQMR Module
======================================

The SUNLinSol_SPTFQMR implementation of the ``SUNLinearSolver`` class performs
a Scaled, Preconditioned, Transpose-Free Quasi-Minimum Residual :cite:p:`Fre:93`
method; this is an iterative linear solver that is designed to be compatible
with any ``N_Vector`` implementation that supports a minimal subset of operations
(:c:func:`N_VClone()`, :c:func:`N_VDotProd()`, :c:func:`N_VScale()`,
:c:func:`N_VLinearSum()`, :c:func:`N_VProd()`, :c:func:`N_VConst()`,
:c:func:`N_VDiv()`, and :c:func:`N_VDestroy()`).  Unlike the SPGMR and
SPFGMR algorithms, SPTFQMR requires a fixed amount of memory that does
not increase with the number of allowed iterations.



.. _SUNLinSol.SPTFQMR.Usage:

SUNLinSol_SPTFQMR Usage
------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_sptfqmr.h``.  The SUNLinSol_SPTFQMR module
is accessible from all SUNDIALS solvers *without*
linking to the ``libsundials_sunlinsolsptfqmr`` module library.

The module SUNLinSol_SPTFQMR provides the following user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_SPTFQMR(N_Vector y, int pretype, int maxl, SUNContext sunctx)

   This constructor function creates and allocates memory for a SPTFQMR
   ``SUNLinearSolver``.

   **Arguments:**
      * *y* -- a template vector.
      * *pretype* -- a flag indicating the type of preconditioning to use:

        * ``SUN_PREC_NONE``
        * ``SUN_PREC_LEFT``
        * ``SUN_PREC_RIGHT``
        * ``SUN_PREC_BOTH``

      * *maxl* -- the number of Krylov basis vectors to use.
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

      Some SUNDIALS solvers are designed to only work with left
      preconditioning (IDA and IDAS) and others with only right
      preconditioning (KINSOL). While it is possible to configure a
      SUNLinSol_SPTFQMR object to use any of the preconditioning options
      with these solvers, this use mode is not supported and may result
      in inferior performance.

   .. note::

      With ``SUN_PREC_RIGHT`` or ``SUN_PREC_BOTH`` the initial guess must be zero (use
      :c:func:`SUNLinSolSetZeroGuess` to indicate the initial guess is zero).



.. c:function:: int SUNLinSol_SPTFQMRSetPrecType(SUNLinearSolver S, int pretype)

   This function updates the flag indicating use of preconditioning.

   **Arguments:**
      * *S* -- SUNLinSol_SPGMR object to update.
      * *pretype* -- a flag indicating the type of preconditioning to use:

        * ``SUN_PREC_NONE``
        * ``SUN_PREC_LEFT``
        * ``SUN_PREC_RIGHT``
        * ``SUN_PREC_BOTH``

   **Return value:**
      * ``SUNLS_SUCCESS`` -- successful update.
      * ``SUNLS_ILL_INPUT`` -- illegal ``pretype``
      * ``SUNLS_MEM_NULL`` -- ``S`` is ``NULL``


.. c:function:: int SUNLinSol_SPTFQMRSetMaxl(SUNLinearSolver S, int maxl)

   This function updates the number of linear solver iterations to allow.

   **Arguments:**
      * *S* -- SUNLinSol_SPTFQMR object to update.
      * *maxl* -- maximum number of linear iterations to allow.  Any
        non-positive input will result in the default value (5).

   **Return value:**
      * ``SUNLS_SUCCESS`` -- successful update.
      * ``SUNLS_MEM_NULL`` -- ``S`` is ``NULL``


.. c:function:: int SUNLinSolSetInfoFile_SPTFQMR(SUNLinearSolver LS, FILE* info_file)

   The function :c:func:`SUNLinSolSetInfoFile_SPTFQMR()` sets the
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


.. c:function:: int SUNLinSolSetPrintLevel_SPTFQMR(SUNLinearSolver LS, int print_level)

   The function :c:func:`SUNLinSolSetPrintLevel_SPTFQMR()` specifies the
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

.. c:function:: SUNLinearSolver SUNSPTFQMR(N_Vector y, int pretype, int maxl)

   Wrapper function for :c:func:`SUNLinSol_SPTFQMR`

.. c:function:: int SUNSPTFQMRSetPrecType(SUNLinearSolver S, int pretype)

   Wrapper function for :c:func:`SUNLinSol_SPTFQMRSetPrecType()`

.. c:function:: int SUNSPTFQMRSetMaxl(SUNLinearSolver S, int maxl)

   Wrapper function for :c:func:`SUNLinSol_SPTFQMRSetMaxl()`





.. _SUNLinSol.SPTFQMR.Description:

SUNLinSol_SPTFQMR Description
---------------------------------


The SUNLinSol_SPTFQMR module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_SPTFQMR {
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
     N_Vector s1;
     N_Vector s2;
     N_Vector r_star;
     N_Vector q;
     N_Vector d;
     N_Vector v;
     N_Vector p;
     N_Vector *r;
     N_Vector u;
     N_Vector vtemp1;
     N_Vector vtemp2;
     N_Vector vtemp3;
     int      print_level;
     FILE*    info_file;
   };

These entries of the *content* field contain the following
information:

* ``maxl`` - number of TFQMR iterations to allow (default is 5),

* ``pretype`` - flag for type of preconditioning to employ
  (default is none),

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

* ``s1, s2`` - vector pointers for supplied scaling matrices
  (default is ``NULL``),

* ``r_star`` - a ``N_Vector`` which holds the initial scaled,
  preconditioned linear system residual,

* ``q, d, v, p, u`` - ``N_Vector`` used for workspace by the SPTFQMR
  algorithm,

* ``r`` - array of two ``N_Vector`` used for workspace within the
  SPTFQMR algorithm,

* ``vtemp1, vtemp2, vtemp3`` - temporary vector storage.

* ``print_level`` - controls the amount of information to be printed to the info file

* ``info_file``   - the file where all informative (non-error) messages will be directed


This solver is constructed to perform the following operations:

* During construction all ``N_Vector`` solver data is allocated,
  with vectors cloned from a template ``N_Vector`` that is input, and
  default solver parameters are set.

* User-facing "set" routines may be called to modify default
  solver parameters.

* Additional "set" routines are called by the SUNDIALS solver
  that interfaces with SUNLinSol_SPTFQMR to supply the
  ``ATimes``, ``PSetup``, and ``Psolve`` function pointers and
  ``s1`` and ``s2`` scaling vectors.

* In the "initialize" call, the solver parameters are checked
  for validity.

* In the "setup" call, any non-``NULL`` ``PSetup`` function is
  called.  Typically, this is provided by the SUNDIALS solver itself,
  that translates between the generic ``PSetup`` function and the
  solver-specific routine (solver-supplied or user-supplied).

* In the "solve" call the TFQMR iteration is performed.  This
  will include scaling and preconditioning if those options have been
  supplied.


The SUNLinSol_SPTFQMR module defines implementations of all
"iterative" linear solver operations listed in
:numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_SPTFQMR``

* ``SUNLinSolInitialize_SPTFQMR``

* ``SUNLinSolSetATimes_SPTFQMR``

* ``SUNLinSolSetPreconditioner_SPTFQMR``

* ``SUNLinSolSetScalingVectors_SPTFQMR``

* ``SUNLinSolSetZeroGuess_SPTFQMR`` -- note the solver assumes a non-zero guess
  by default and the zero guess flag is reset to ``SUNFALSE`` after each call to
  :c:func:`SUNLinSolSolve_SPTFQMR`.

* ``SUNLinSolSetup_SPTFQMR``

* ``SUNLinSolSolve_SPTFQMR``

* ``SUNLinSolNumIters_SPTFQMR``

* ``SUNLinSolResNorm_SPTFQMR``

* ``SUNLinSolResid_SPTFQMR``

* ``SUNLinSolLastFlag_SPTFQMR``

* ``SUNLinSolSpace_SPTFQMR``

* ``SUNLinSolFree_SPTFQMR``
