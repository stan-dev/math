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

.. _SUNLinSol.SPGMR:

The SUNLinSol_SPGMR Module
======================================

The SUNLinSol_SPGMR implementation of the ``SUNLinearSolver`` class performs
a Scaled, Preconditioned, Generalized Minimum Residual :cite:p:`SaSc:86` method;
this is an iterative linear solver that is designed to be compatible with any
``N_Vector`` implementation that supports a minimal subset of operations
(:c:func:`N_VClone()`, :c:func:`N_VDotProd()`, :c:func:`N_VScale()`,
:c:func:`N_VLinearSum()`, :c:func:`N_VProd()`, :c:func:`N_VConst()`,
:c:func:`N_VDiv()`, and :c:func:`N_VDestroy()`).



.. _SUNLinSol.SPGMR.Usage:

SUNLinSol_SPGMR Usage
--------------------------

The header file to be included when using this module
is ``sunlinsol/sunlinsol_spgmr.h``.  The SUNinSol_SPGMR module
is accessible from all SUNDIALS solvers *without*
linking to the ``libsundials_sunlinsolspgmr`` module library.


The module SUNLinSol_SPGMR provides the following
user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_SPGMR(N_Vector y, int pretype, int maxl, SUNContext sunctx)

   This constructor function creates and allocates memory for a SPGMR
   ``SUNLinearSolver``.

   **Arguments:**
      * *y* -- a template vector.
      * *pretype* -- a flag indicating the type of preconditioning to use:

        * ``SUN_PREC_NONE``
        * ``SUN_PREC_LEFT``
        * ``SUN_PREC_RIGHT``
        * ``SUN_PREC_BOTH``

      * *maxl* -- the number of Krylov basis vectors to use.

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
      SUNLinSol_SPGMR object to use any of the preconditioning options
      with these solvers, this use mode is not supported and may result
      in inferior performance.


.. c:function:: int SUNLinSol_SPGMRSetPrecType(SUNLinearSolver S, int pretype)

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


.. c:function:: int SUNLinSol_SPGMRSetGSType(SUNLinearSolver S, int gstype)

   This function sets the type of Gram-Schmidt orthogonalization to use.

   **Arguments:**
      * *S* -- SUNLinSol_SPGMR object to update.
      * *gstype* -- a flag indicating the type of orthogonalization to use:

        * ``SUN_MODIFIED_GS``
        * ``SUN_CLASSICAL_GS``

   **Return value:**
      * ``SUNLS_SUCCESS`` -- successful update.
      * ``SUNLS_ILL_INPUT`` -- illegal ``gstype``
      * ``SUNLS_MEM_NULL`` -- ``S`` is ``NULL``


.. c:function:: int SUNLinSol_SPGMRSetMaxRestarts(SUNLinearSolver S, int maxrs)

   This function sets the number of GMRES restarts to allow.

   **Arguments:**
      * *S* -- SUNLinSol_SPGMR object to update.
      * *maxrs* -- maximum number of restarts to allow.  A negative input will
        result in the default of 0.

   **Return value:**
      * ``SUNLS_SUCCESS`` -- successful update.
      * ``SUNLS_MEM_NULL`` -- ``S`` is ``NULL``


.. c:function:: int SUNLinSolSetInfoFile_SPGMR(SUNLinearSolver LS, FILE* info_file)

   The function :c:func:`SUNLinSolSetInfoFile_SPGMR()` sets the
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


.. c:function:: int SUNLinSolSetPrintLevel_SPGMR(SUNLinearSolver LS, int print_level)

   The function :c:func:`SUNLinSolSetPrintLevel_SPGMR()` specifies the
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


For backwards compatibility, we also provide the wrapper functions,
each with identical input and output arguments to the routines that
they wrap:

.. c:function:: SUNLinearSolver SUNSPGMR(N_Vector y, int pretype, int maxl)

   Wrapper function for :c:func:`SUNLinSol_SPGMR`

.. c:function:: int SUNSPGMRSetPrecType(SUNLinearSolver S, int pretype)

   Wrapper function for :c:func:`SUNLinSol_SPGMRSetPrecType()`

.. c:function:: int SUNSPGMRSetGSType(SUNLinearSolver S, int gstype)

   Wrapper function for :c:func:`SUNLinSol_SPGMRSetGSType()`

.. c:function:: int SUNSPGMRSetMaxRestarts(SUNLinearSolver S, int maxrs)

   Wrapper function for :c:func:`SUNLinSol_SPGMRSetMaxRestarts()`





.. _SUNLinSol.SPGMR.Description:

SUNLinSol_SPGMR Description
-----------------------------


The SUNLinSol_SPGMR module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_SPGMR {
     int maxl;
     int pretype;
     int gstype;
     int max_restarts;
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
     N_Vector *V;
     realtype **Hes;
     realtype *givens;
     N_Vector xcor;
     realtype *yg;
     N_Vector vtemp;
     int      print_level;
     FILE*    info_file;
   };

These entries of the *content* field contain the following
information:

* ``maxl`` - number of GMRES basis vectors to use (default is 5),

* ``pretype`` - flag for type of preconditioning to employ
  (default is none),

* ``gstype`` - flag for type of Gram-Schmidt orthogonalization
  (default is modified Gram-Schmidt),

* ``max_restarts`` - number of GMRES restarts to allow
  (default is 0),

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

* ``V`` - the array of Krylov basis vectors
  :math:`v_1, \ldots, v_{\text{maxl}+1}`, stored in
  ``V[0], ... V[maxl]``. Each :math:`v_i` is a vector of type
  ``N_Vector``,

* ``Hes`` - the :math:`(\text{maxl}+1)\times\text{maxl}`
  Hessenberg matrix. It is stored row-wise so that the (i,j)th
  element is given by ``Hes[i][j]``,

* ``givens`` - a length :math:`2\,\text{maxl}` array which represents
  the Givens rotation matrices that arise in the GMRES
  algorithm. These matrices are :math:`F_0, F_1, \ldots, F_j`, where

  .. math::

     F_i = \begin{bmatrix}
        1 &        &   &     &      &   &        &   \\
          & \ddots &   &     &      &   &        &   \\
          &        & 1 &     &      &   &        &   \\
          &        &   & c_i & -s_i &   &        &   \\
          &        &   & s_i &  c_i &   &        &   \\
          &        &   &     &      & 1 &        &   \\
          &        &   &     &      &   & \ddots &   \\
          &        &   &     &      &   &        & 1\end{bmatrix},

  are represented in the ``givens`` vector as
  ``givens[0]`` :math:`= c_0`,
  ``givens[1]`` :math:`= s_0`,
  ``givens[2]`` :math:`= c_1`,
  ``givens[3]`` :math:`= s_1`, :math:`\ldots`,
  ``givens[2j]`` :math:`= c_j`,
  ``givens[2j+1]`` :math:`= s_j`,

* ``xcor`` - a vector which holds the scaled, preconditioned
  correction to the initial guess,

* ``yg`` - a length :math:`(\text{maxl}+1)` array of ``realtype``
  values used to hold "short" vectors (e.g. :math:`y` and :math:`g`),

* ``vtemp`` - temporary vector storage.

* ``print_level`` - controls the amount of information to be printed to the info file

* ``info_file``   - the file where all informative (non-error) messages will be directed


This solver is constructed to perform the following operations:

* During construction, the ``xcor`` and ``vtemp`` arrays are
  cloned from a template ``N_Vector`` that is input, and default
  solver parameters are set.

* User-facing "set" routines may be called to modify default
  solver parameters.

* Additional "set" routines are called by the SUNDIALS solver
  that interfaces with SUNLinSol_SPGMR to supply the
  ``ATimes``, ``PSetup``, and ``Psolve`` function pointers and
  ``s1`` and ``s2`` scaling vectors.

* In the "initialize" call, the remaining solver data is
  allocated (``V``, ``Hes``, ``givens``, and ``yg`` )

* In the "setup" call, any non-``NULL``
  ``PSetup`` function is called.  Typically, this is provided by
  the SUNDIALS solver itself, that translates between the generic
  ``PSetup`` function and the solver-specific routine (solver-supplied
  or user-supplied).

* In the "solve" call, the GMRES iteration is performed.  This
  will include scaling, preconditioning, and restarts if those options
  have been supplied.

The SUNLinSol_SPGMR module defines implementations of all
"iterative" linear solver operations listed in
:numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_SPGMR``

* ``SUNLinSolInitialize_SPGMR``

* ``SUNLinSolSetATimes_SPGMR``

* ``SUNLinSolSetPreconditioner_SPGMR``

* ``SUNLinSolSetScalingVectors_SPGMR``

* ``SUNLinSolSetZeroGuess_SPGMR`` -- note the solver assumes a non-zero guess by
  default and the zero guess flag is reset to ``SUNFALSE`` after each call to
  :c:func:`SUNLinSolSolve_SPGMR`.

* ``SUNLinSolSetup_SPGMR``

* ``SUNLinSolSolve_SPGMR``

* ``SUNLinSolNumIters_SPGMR``

* ``SUNLinSolResNorm_SPGMR``

* ``SUNLinSolResid_SPGMR``

* ``SUNLinSolLastFlag_SPGMR``

* ``SUNLinSolSpace_SPGMR``

* ``SUNLinSolFree_SPGMR``
