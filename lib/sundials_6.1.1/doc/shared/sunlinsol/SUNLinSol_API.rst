..
   Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNLinSol.API:

The SUNLinearSolver API
=============================

The SUNLinSol API defines several linear solver operations that enable
SUNDIALS packages to utilize this API. These functions can be divided into
three categories. The first are the core linear solver functions. The
second consist of "set" routines to supply the linear solver with functions
provided by the SUNDIALS packages and to modify solver parameters. The
final group consists of "get" routines for retrieving linear solver
statistics. All of these functions are defined in the header file
``sundials/sundials_linearsolver.h``.



.. _SUNLinSol.CoreFn:

SUNLinearSolver core functions
-----------------------------------------------------

The core linear solver functions consist of two **required**
functions: :c:func:`SUNLinSolGetType` returns the linear solver
type, and :c:func:`SUNLinSolSolve` solves the linear system :math:`Ax=b`.

The remaining **optional** functions return the solver ID
(:c:func:`SUNLinSolGetID`), initialize the linear solver object once
all solver-specific options have been set (:c:func:`SUNLinSolInitialize`),
set up the linear solver object to utilize an updated matrix :math:`A`
(:c:func:`SUNLinSolSetup`), and destroy a linear solver object
(:c:func:`SUNLinSolFree`).


.. c:function:: SUNLinearSolver_Type SUNLinSolGetType(SUNLinearSolver LS)

   Returns the type identifier for the linear solver *LS*.

   **Return value:**

      * ``SUNLINEARSOLVER_DIRECT (0)`` -- the SUNLinSol module
        requires a matrix, and computes an "exact" solution to the linear
        system defined by that matrix.

      * ``SUNLINEARSOLVER_ITERATIVE (1)`` -- the SUNLinSol module does
        not require a matrix (though one may be provided), and computes
        an inexact solution to the linear system using a matrix-free
        iterative algorithm. That is it solves the linear system defined
        by the package-supplied ``ATimes`` routine (see
        :c:func:`SUNLinSolSetATimes()` below), even if that linear system
        differs from the one encoded in the matrix object (if one is
        provided). As the solver computes the solution only inexactly (or
        may diverge), the linear solver should check for solution
        convergence/accuracy as appropriate.

      * ``SUNLINEARSOLVER_MATRIX_ITERATIVE (2)`` -- the SUNLinSol
        module requires a matrix, and computes an inexact solution to the
        linear system defined by that matrix using an iterative
        algorithm. That is it solves the linear system defined by the
        matrix object even if that linear system differs from that
        encoded by the package-supplied ``ATimes`` routine. As the solver
        computes the solution only inexactly (or may diverge), the linear
        solver should check for solution convergence/accuracy as
        appropriate.

      * ``SUNLINEARSOLVER_MATRIX_EMBEDDED (3)`` -- the SUNLinSol module sets up
        and solves the specified linear system at each linear solve call.  Any
        matrix-related data structures are held internally to the linear solver itself,
        and are not provided by the SUNDIALS package.

   **Usage:**

      .. code-block:: c

         type = SUNLinSolGetType(LS);

   .. note::

      See :numref:`SUNLinSol.Intended` for more information on
      intended use cases corresponding to the linear solver type.


.. c:function:: SUNLinearSolver_ID SUNLinSolGetID(SUNLinearSolver LS)

   Returns a non-negative linear solver identifier (of type ``int``)
   for the linear solver *LS*.

   **Return value:**

      Non-negative linear solver identifier (of type ``int``), defined by the
      enumeration ``SUNLinearSolver_ID``, with values shown in
      :numref:`SUNLinSol.API.IDs` and defined in the ``sundials_linearsolver.h``
      header file.

   **Usage:**

      .. code-block:: c

         id = SUNLinSolGetID(LS);

   .. note::

      It is recommended that a user-supplied ``SUNLinearSolver`` return the
      ``SUNLINEARSOLVER_CUSTOM`` identifier.


.. c:function:: int SUNLinSolInitialize(SUNLinearSolver LS)

   Performs linear solver initialization (assuming that all
   solver-specific options have been set).

   **Return value:**

      Zero for a successful call, and a negative value for a failure.
      Ideally, this should return one of the generic error codes listed in
      :numref:`SUNLinSol.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolInitialize(LS);


.. c:function:: int SUNLinSolSetup(SUNLinearSolver LS, SUNMatrix A)

   Performs any linear solver setup needed, based on an updated system
   ``SUNMatrix`` *A*.  This may be called frequently (e.g., with a full
   Newton method) or infrequently (for a modified Newton method), based
   on the type of integrator and/or nonlinear solver requesting the
   solves.

   **Return value:**

      Zero for a successful call, a positive value for a recoverable failure,
      and a negative value for an unrecoverable failure.  Ideally this should
      return one of the generic error codes listed in
      :numref:`SUNLinSol.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolSetup(LS, A);


.. c:function:: int SUNLinSolSolve(SUNLinearSolver LS, SUNMatrix A, N_Vector x, N_Vector b, realtype tol)

   This *required* function solves a linear system :math:`Ax = b`.

   **Arguments:**

      * *LS* -- a SUNLinSol object.
      * *A* -- a ``SUNMatrix`` object.
      * *x* -- an ``N_Vector`` object containing the initial guess for
        the solution of the linear system on input, and the solution to the
        linear system upon return.
      * *b* -- an ``N_Vector`` object containing the linear system
        right-hand side.
      * *tol* -- the desired linear solver tolerance.

   **Return value:**

      Zero for a successful call, a positive value for a recoverable failure,
      and a negative value for an unrecoverable failure.  Ideally this should
      return one of the generic error codes listed in
      :numref:`SUNLinSol.ErrorCodes`.

   **Notes:**

      **Direct solvers:** can ignore the *tol* argument.

      **Matrix-free solvers:** (those that identify as
      ``SUNLINEARSOLVER_ITERATIVE``) can ignore the ``SUNMatrix`` input
      *A*, and should rely on the matrix-vector product function supplied
      through the routine :c:func:`SUNLinSolSetATimes()`.

      **Iterative solvers:** (those that identify as
      ``SUNLINEARSOLVER_ITERATIVE`` or
      ``SUNLINEARSOLVER_MATRIX_ITERATIVE``) should attempt to solve to
      the specified tolerance *tol* in a weighted 2-norm. If the solver
      does not support scaling then it should just use a 2-norm.

      **Matrix-embedded solvers:** should ignore the ``SUNMatrix`` input *A*
      as this will be ``NULL``.  It is assumed that within this function, the
      solver will call interface routines from the relevant SUNDIALS package to
      directly form the linear system matrix :math:`A`, and then solve
      :math:`Ax=b` before returning with the solution :math:`x`.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolSolve(LS, A, x, b, tol);


.. c:function:: int SUNLinSolFree(SUNLinearSolver LS)

   Frees memory allocated by the linear solver.

   **Return value:**

      Zero for a successful call, and a negative value for a failure.
      Ideally, this should return one of the generic error codes listed in
      :numref:`SUNLinSol.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolFree(LS);




.. _SUNLinSol.SetFn:

SUNLinearSolver "set" functions
-------------------------------------

The following functions supply linear solver modules with functions defined
by the SUNDIALS packages and modify solver parameters.  Only the routine
for setting the matrix-vector product routine is required, and even then is
only required for matrix-free linear solver modules.  Otherwise, all other
set functions are optional.  SUNLinSol implementations that do not provide
the functionality for any optional routine should leave the corresponding
function pointer ``NULL`` instead of supplying a dummy routine.


.. c:function:: int SUNLinSolSetATimes(SUNLinearSolver LS, void* A_data, SUNATimesFn ATimes)

   *Required for matrix-free linear solvers* (otherwise optional).

   Provides a :c:type:`SUNATimesFn` function pointer, as well as a ``void*``
   pointer to a data structure used by this routine, to the linear
   solver object *LS*.  SUNDIALS packages call this function to set the
   matrix-vector product function to either a solver-provided
   difference-quotient via vector operations or a user-supplied
   solver-specific routine.

   **Return value:**

      Zero for a successful call, and a negative value for a failure.
      Ideally, this should return one of the generic error codes listed in
      :numref:`SUNLinSol.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolSetATimes(LS, A_data, ATimes);


.. c:function:: int SUNLinSolSetPreconditioner(SUNLinearSolver LS, void* P_data, SUNPSetupFn Pset, SUNPSolveFn Psol)

   This *optional* routine provides :c:type:`SUNPSetupFn` and
   :c:type:`SUNPSolveFn` function pointers that implement the
   preconditioner solves :math:`P_1^{-1}` and :math:`P_2^{-1}` from
   :eq:`eq:transformed_linear_system_components`. This
   routine is called by a SUNDIALS package, which provides
   translation between the generic *Pset* and *Psol* calls and the
   package- or user-supplied routines.

   **Return value:**

      Zero for a successful call, and a negative value for a failure.
      Ideally, this should return one of the generic error codes listed in
      :numref:`SUNLinSol.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolSetPreconditioner(LS, Pdata, Pset, Psol);


.. c:function:: int SUNLinSolSetScalingVectors(SUNLinearSolver LS, N_Vector s1, N_Vector s2)

   This *optional* routine provides left/right scaling vectors for the
   linear system solve.  Here, *s1* and *s2* are ``N_Vectors`` of positive
   scale factors containing the diagonal of the matrices :math:`S_1`
   and :math:`S_2` from :eq:`eq:transformed_linear_system_components`, respectively.
   Neither vector needs to be tested for positivity, and a ``NULL`` argument for either
   indicates that the corresponding scaling matrix is the
   identity.

   **Return value:**

      Zero for a successful call, and a negative value for a failure.
      Ideally, this should return one of the generic error codes listed in
      :numref:`SUNLinSol.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolSetScalingVectors(LS, s1, s2);


.. c:function:: int SUNLinSolSetZeroGuess(SUNLinearSolver LS, booleantype onoff)

   This *optional* routine indicates if the upcoming :c:func:`SUNlinSolSolve` call
   will be made with a zero initial guess (``SUNTRUE``) or a non-zero initial
   guess (``SUNFALSE``).

   **Return value:**

      Zero for a successful call, and a negative value for a failure.
      Ideally, this should return one of the generic error codes listed in
      :numref:`SUNLinSol.ErrorCodes`.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolSetZeroGuess(LS, onoff);

   **Notes:**

      It is assumed that the initial guess status is not retained across
      calls to :c:func:`SUNLinSolSolve`. As such, the linear solver interfaces in
      each of the SUNDIALS packages call :c:func:`SUNLinSolSetZeroGuess` prior to
      each call to :c:func:`SUNLinSolSolve`.


.. _SUNLinSol.GetFn:

SUNLinearSolver "get" functions
----------------------------------

The following functions allow SUNDIALS packages to retrieve results from a
linear solve.  *All routines are optional.*


.. c:function:: int SUNLinSolNumIters(SUNLinearSolver LS)

   This *optional* routine should return the number of linear
   iterations performed in the most-recent "solve" call.

   **Usage:**

      .. code-block:: c

         its = SUNLinSolNumIters(LS);


.. c:function:: realtype SUNLinSolResNorm(SUNLinearSolver LS)

   This *optional* routine should return the final residual norm from
   the most-recent "solve" call.

   **Usage:**

      .. code-block:: c

         rnorm = SUNLinSolResNorm(LS);


.. c:function:: N_Vector SUNLinSolResid(SUNLinearSolver LS)

   If an iterative method computes the preconditioned initial residual
   and returns with a successful solve without performing any
   iterations (i.e., either the initial guess or the preconditioner is
   sufficiently accurate), then this *optional* routine may be called
   by the SUNDIALS package.  This routine should return the ``N_Vector``
   containing the preconditioned initial residual vector.

   **Usage:**

      .. code-block:: c

         rvec = SUNLinSolResid(LS);

   **Notes:**

      Since ``N_Vector`` is actually a pointer, and the results are
      not modified, this routine should *not* require additional memory
      allocation.  If the SUNLinSol object does not retain a vector for
      this purpose, then this function pointer should be set to ``NULL``
      in the implementation.


.. c:function:: sunindextype SUNLinSolLastFlag(SUNLinearSolver LS)

   This *optional* routine should return the last error flag
   encountered within the linear solver.  Although not called by the
   SUNDIALS packages directly, this may be called by the user to
   investigate linear solver issues after a failed solve.

   **Usage:**

      .. code-block:: c

         lflag = SUNLinLastFlag(LS);


.. c:function:: int SUNLinSolSpace(SUNLinearSolver LS, long int *lenrwLS, long int *leniwLS)

   This *optional* routine should return the storage requirements for
   the linear solver *LS*:

   * *lrw* is a ``long int`` containing the number of realtype words
   * *liw* is a ``long int`` containing the number of integer words.

   The return value is an integer flag denoting success/failure of the operation.

   This function is advisory only, for use by users to help determine
   their total space requirements.

   **Usage:**

      .. code-block:: c

         retval = SUNLinSolSpace(LS, &lrw, &liw);





.. _SUNLinSol.SUNSuppliedFn:

Functions provided by SUNDIALS packages
---------------------------------------------

To interface with SUNLinSol modules, the SUNDIALS packages supply a
variety of routines for evaluating the matrix-vector product, and
setting up and applying the preconditioner.  These package-provided
routines translate between the user-supplied ODE, DAE, or nonlinear
systems and the generic linear solver API. The function types for
these routines are defined in the header file
``sundials/sundials_iterative.h``, and are described below.


.. c:type:: int (*SUNATimesFn)(void *A_data, N_Vector v, N_Vector z)

   Computes the action of a matrix on a vector, performing the
   operation :math:`z \gets Av`.  Memory for *z* will already be
   allocated prior to calling this function.  The parameter
   *A_data* is a pointer to any information about :math:`A` which
   the function needs in order to do its job. The vector :math:`v`
   should be left unchanged.

   **Return value:**

      Zero for a successful call, and non-zero upon failure.


.. c:type:: int (*SUNPSetupFn)(void *P_data)

   Sets up any requisite problem data in preparation for calls
   to the corresponding :c:type:`SUNPSolveFn`.


   **Return value:**

      Zero for a successful call, and non-zero upon failure.


.. c:type:: int (*SUNPSolveFn)(void *P_data, N_Vector r, N_Vector z, realtype tol, int lr)

   Solves the preconditioner equation :math:`Pz = r` for the vector :math:`z`.
   Memory for *z* will already be allocated prior to calling this function.
   The parameter *P_data* is a pointer to any information about :math:`P`
   which the function needs in order to do its job (set up by the corresponding
   :c:type:`SUNPSetupFn`). The parameter *lr* is input, and indicates
   whether :math:`P` is to be taken as the left or right
   preconditioner: *lr* = 1 for left and *lr* = 2 for right.  If
   preconditioning is on one side only, *lr* can be ignored.  If the
   preconditioner is iterative, then it should strive to solve the
   preconditioner equation so that

   .. math::

      \| Pz - r \|_{\text{wrms}} < tol

   where the error weight vector for the WRMS norm may be accessed
   from the main package memory structure.  The vector *r* should not
   be modified by the *SUNPSolveFn*.

   **Return value:**

      Zero for a successful call, a negative value for an
      unrecoverable failure condition, or a positive value for a
      recoverable failure condition (thus the calling routine may
      reattempt the solution after updating preconditioner data).


.. _SUNLinSol.ReturnCodes:

SUNLinearSolver return codes
------------------------------------

The functions provided to SUNLinSol modules by each SUNDIALS package,
and functions within the SUNDIALS-provided SUNLinSol implementations,
utilize a common set of return codes, listed in
:numref:`SUNLinSol.ErrorCodes`.  These adhere to a common pattern:

* 0 indicates success
* a positive value corresponds to a recoverable failure, and
* a negative value indicates a non-recoverable failure.

Aside from this pattern, the actual values of each error code
provide additional information to the user in case of a linear solver failure.


.. _SUNLinSol.ErrorCodes:
.. table:: SUNLinSol error codes
   :align: center

   +------------------------------+-------+---------------------------------------------------+
   | Error code                   | Value | Meaning                                           |
   +==============================+=======+===================================================+
   | ``SUNLS_SUCCESS``            | 0     | successful call or converged solve                |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_MEM_NULL``           | -801  | the memory argument to the function is ``NULL``   |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_ILL_INPUT``          | -802  | an illegal input has been provided to the         |
   |                              |       | function                                          |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_MEM_FAIL``           | -803  | failed memory access or allocation                |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_ATIMES_NULL``        | -804  | the ``Atimes`` function is ``NULL``               |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_ATIMES_FAIL_UNREC``  | -805  | an unrecoverable failure occurred in the          |
   |                              |       | ``ATimes`` routine                                |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_PSET_FAIL_UNREC``    | -806  | an unrecoverable failure occurred in the ``Pset`` |
   |                              |       | routine                                           |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_PSOLVE_NULL``        | -807  | the preconditioner solve function is ``NULL``     |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_PSOLVE_FAIL_UNREC``  | -808  | an unrecoverable failure occurred in the          |
   |                              |       | ``Psolve`` routine                                |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_PACKAGE_FAIL_UNREC`` | -809  | an unrecoverable failure occurred in an external  |
   |                              |       | linear solver package                             |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_GS_FAIL``            | -810  | a failure occurred during Gram-Schmidt            |
   |                              |       | orthogonalization (SPGMR/SPFGMR)                  |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_QRSOL_FAIL``         | -811  | a singular $R$ matrix was encountered in a QR     |
   |                              |       | factorization (SPGMR/SPFGMR)                      |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_VECTOROP_ERR``       | -812  | a vector operation error occurred                 |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_RES_REDUCED``        | 801   | an iterative solver reduced the residual, but did |
   |                              |       | not converge to the desired tolerance             |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_CONV_FAIL``          | 802   | an iterative solver did not converge (and the     |
   |                              |       | residual was not reduced)                         |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_ATIMES_FAIL_REC``    | 803   | a recoverable failure occurred in the ``ATimes``  |
   |                              |       | routine                                           |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_PSET_FAIL_REC``      | 804   | a recoverable failure occurred in the ``Pset``    |
   |                              |       | routine                                           |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_PSOLVE_FAIL_REC``    | 805   | a recoverable failure occurred in the ``Psolve``  |
   |                              |       | routine                                           |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_PACKAGE_FAIL_REC``   | 806   | a recoverable failure occurred in an external     |
   |                              |       | linear solver package                             |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_QRFACT_FAIL``        | 807   | a singular matrix was encountered during a QR     |
   |                              |       | factorization (SPGMR/SPFGMR)                      |
   +------------------------------+-------+---------------------------------------------------+
   | ``SUNLS_LUFACT_FAIL``        | 808   | a singular matrix was encountered during a LU     |
   |                              |       | factorization                                     |
   +------------------------------+-------+---------------------------------------------------+



.. _SUNLininSol.Generic:

The generic SUNLinearSolver module
-----------------------------------------

SUNDIALS packages interact with specific SUNLinSol implementations
through the generic SUNLinearSolver abstract base class.  The
``SUNLinearSolver`` type is a pointer to a structure containing an
implementation-dependent *content* field, and an *ops* field, and is
defined as

.. c:type:: struct _generic_SUNLinearSolver *SUNLinearSolver

and the generic structure is defined as

.. code-block:: c

   struct _generic_SUNLinearSolver {
     void *content;
     struct _generic_SUNLinearSolver_Ops *ops;
   };

where the ``_generic_SUNLinearSolver_Ops`` structure is a list of
pointers to the various actual linear solver operations provided by a
specific implementation.  The ``_generic_SUNLinearSolver_Ops``
structure is defined as

.. code-block:: c

   struct _generic_SUNLinearSolver_Ops {
     SUNLinearSolver_Type (*gettype)(SUNLinearSolver);
     SUNLinearSolver_ID   (*getid)(SUNLinearSolver);
     int                  (*setatimes)(SUNLinearSolver, void*, SUNATimesFn);
     int                  (*setpreconditioner)(SUNLinearSolver, void*,
                                               SUNPSetupFn, SUNPSolveFn);
     int                  (*setscalingvectors)(SUNLinearSolver,
                                               N_Vector, N_Vector);
     int                  (*setzeroguess)(SUNLinearSolver, booleantype);
     int                  (*initialize)(SUNLinearSolver);
     int                  (*setup)(SUNLinearSolver, SUNMatrix);
     int                  (*solve)(SUNLinearSolver, SUNMatrix, N_Vector,
                                   N_Vector, realtype);
     int                  (*numiters)(SUNLinearSolver);
     realtype             (*resnorm)(SUNLinearSolver);
     sunindextype         (*lastflag)(SUNLinearSolver);
     int                  (*space)(SUNLinearSolver, long int*, long int*);
     N_Vector             (*resid)(SUNLinearSolver);
     int                  (*free)(SUNLinearSolver);
   };


The generic SUNLinSol class defines and implements the linear solver
operations defined in :numref:`SUNLinSol.CoreFn` -- :numref:`SUNLinSol.GetFn`.
These routines are in fact only wrappers to the linear solver operations
defined by a particular SUNLinSol implementation, which are accessed through
the *ops* field of the ``SUNLinearSolver`` structure.  To illustrate this
point we show below the implementation of a typical linear solver operation
from the ``SUNLinearSolver`` base class, namely :c:func:`SUNLinSolInitialize`,
that initializes a ``SUNLinearSolver`` object for use after it has been
created and configured, and returns a flag denoting a successful or failed
operation:

.. code-block:: c

   int SUNLinSolInitialize(SUNLinearSolver S)
   {
     return ((int) S->ops->initialize(S));
   }



.. _SUNLinSol.API.Compatibility:

Compatibility of SUNLinearSolver modules
---------------------------------------------

Not all ``SUNLinearSolver`` implementations are compatible with all
``SUNMatrix`` and ``N_Vector`` implementations provided in SUNDIALS.
More specifically, all of the SUNDIALS iterative linear solvers
(:ref:`SPGMR <SUNLinSol.SPGMR>`, :ref:`SPFGMR <SUNLinSol.SPFGMR>`,
:ref:`SPBCGS <SUNLinSol.SPBCGS>`, :ref:`SPTFQMR <SUNLinSol.SPTFQMR>`, and
:ref:`PCG <SUNLinSol.PCG>`) are compatible with all of the SUNDIALS
``N_Vector`` modules, but the matrix-based direct SUNLinSol modules
are specifically designed to work with distinct ``SUNMatrix`` and
``N_Vector`` modules.  In the list below, we summarize the
compatibility of each matrix-based ``SUNLinearSolver``
module with the various ``SUNMatrix`` and ``N_Vector`` modules.  For
a more thorough discussion of these compatibilities, we defer to the
documentation for each individual SUNLinSol module in the sections
that follow.

* :ref:`Dense <SUNLinSol_Dense>`

  * ``SUNMatrix``: :ref:`Dense <SUNMatrix.Dense>` or user-supplied

  * ``N_Vector``: :ref:`Serial <NVectors.NVSerial>`,
    :ref:`OpenMP <NVectors.OpenMP>`, :ref:`Pthreads <NVectors.Pthreads>`,
    or user-supplied

* :ref:`LapackDense <SUNLinSol_LapackDense>`

  * ``SUNMatrix``: :ref:`Dense <SUNMatrix.Dense>` or user-supplied

  * ``N_Vector``: :ref:`Serial <NVectors.NVSerial>`,
    :ref:`OpenMP <NVectors.OpenMP>`, :ref:`Pthreads <NVectors.Pthreads>`,
    or user-supplied

* :ref:`Band <SUNLinSol_Band>`

  * ``SUNMatrix``: :ref:`Band <SUNMatrix.Band>` or user-supplied

  * ``N_Vector``: :ref:`Serial <NVectors.NVSerial>`,
    :ref:`OpenMP <NVectors.OpenMP>`, :ref:`Pthreads <NVectors.Pthreads>`,
    or user-supplied

* :ref:`LapackBand <SUNLinSol_LapackBand>`

  * ``SUNMatrix``: :ref:`Band <SUNMatrix.Band>` or user-supplied

  * ``N_Vector``: :ref:`Serial <NVectors.NVSerial>`,
    :ref:`OpenMP <NVectors.OpenMP>`, :ref:`Pthreads <NVectors.Pthreads>`,
    or user-supplied

* :ref:`KLU <SUNLinSol.KLU>`

  * ``SUNMatrix``: :ref:`Sparse <SUNMatrix.Sparse>` or user-supplied

  * ``N_Vector``: :ref:`Serial <NVectors.NVSerial>`,
    :ref:`OpenMP <NVectors.OpenMP>`, :ref:`Pthreads <NVectors.Pthreads>`,
    or user-supplied

* :ref:`SuperLU_MT <SUNLinSol.SuperLUMT>`

  * ``SUNMatrix``: :ref:`Sparse <SUNMatrix.Sparse>` or user-supplied

  * ``N_Vector``: :ref:`Serial <NVectors.NVSerial>`,
    :ref:`OpenMP <NVectors.OpenMP>`, :ref:`Pthreads <NVectors.Pthreads>`,
    or user-supplied

* :ref:`SuperLU_Dist <SUNLinSol.SuperLUDIST>`

  * ``SUNMatrix``: :ref:`SLUNRLOC <SUNMatrix.SLUNRloc>` or user-supplied

  * ``N_Vector``: :ref:`Serial <NVectors.NVSerial>`,
    :ref:`OpenMP <NVectors.OpenMP>`, :ref:`Pthreads <NVectors.Pthreads>`,
    :ref:`Parallel <NVectors.NVParallel>`, :ref:`*hypre* <NVectors.ParHyp>`,
    :ref:`PETSc <NVectors.NVPETSc>`, or user-supplied

* :ref:`Magma Dense <SUNLinSol.MagmaDense>`

  * ``SUNMatrix``: :ref:`Magma Dense <SUNMatrix.MagmaDense>` or user-supplied

  * ``N_Vector``: :ref:`HIP <NVectors.HIP>`, :ref:`RAJA <NVectors.RAJA>`, or user-supplied

* :ref:`OneMKL Dense <SUNLinSol.OneMklDense>`

  * ``SUNMatrix``: :ref:`One MKL Dense <SUNMatrix.OneMklDense>` or user-supplied

  * ``N_Vector``: :ref:`SYCL <NVectors.SYCL>`, :ref:`RAJA <NVectors.RAJA>`, or user-supplied

* :ref:`cuSolverSp batchQR <SUNLinSol.cuSolverSp>`

  * ``SUNMatrix``: :ref:`cuSparse <SUNMatrix.cuSparse>` or user-supplied

  * ``N_Vector``: :ref:`CUDA <NVectors.CUDA>`, :ref:`RAJA <NVectors.RAJA>`, or user-supplied



.. _SUNLinSol.API.Custom:

Implementing a custom SUNLinearSolver module
--------------------------------------------------

A particular implementation of the ``SUNLinearSolver`` module must:

* Specify the *content* field of the SUNLinSol module.

* Define and implement the required linear solver operations.

  .. note::

     The names of these routines should be unique to that
     implementation in order to permit using more than one
     SUNLinSol module (each with different ``SUNLinearSolver``
     internal data representations) in the same code.

* Define and implement user-callable constructor and destructor
  routines to create and free a ``SUNLinearSolver`` with
  the new *content* field and with *ops* pointing to the
  new linear solver operations.

We note that the function pointers for all unsupported optional
routines should be set to ``NULL`` in the *ops* structure.  This
allows the SUNDIALS package that is using the SUNLinSol object
to know whether the associated functionality is supported.

To aid in the creation of custom ``SUNLinearSolver`` modules the generic
``SUNLinearSolver`` module provides the utility function
:c:func:`SUNLinSolNewEmpty`. When used in custom ``SUNLinearSolver``
constructors this function will ease the introduction of any new optional linear
solver operations to the ``SUNLinearSolver`` API by ensuring that only required
operations need to be set.

.. c:function:: SUNLinearSolver SUNLinSolNewEmpty()

   This function allocates a new generic ``SUNLinearSolver`` object and
   initializes its content pointer and the function pointers in the operations
   structure to ``NULL``.

   **Return value:**

      If successful, this function returns a ``SUNLinearSolver`` object.
      If an error occurs when allocating the object, then this routine will
      return ``NULL``.

.. c:function:: void SUNLinSolFreeEmpty(SUNLinearSolver LS)

   This routine frees the generic ``SUNLinearSolver`` object, under the
   assumption that any implementation-specific data that was allocated
   within the underlying content structure has already been freed.
   It will additionally test whether the ops pointer is ``NULL``,
   and, if it is not, it will free it as well.

   **Arguments:**

      * *LS* -- a SUNLinearSolver object


Additionally, a ``SUNLinearSolver`` implementation *may* do the following:

* Define and implement additional user-callable "set" routines
  acting on the ``SUNLinearSolver``, e.g., for setting various
  configuration options to tune the linear solver for a particular
  problem.

* Provide additional user-callable "get" routines acting on the
  ``SUNLinearSolver`` object, e.g., for returning various solve
  statistics.



Each SUNLinSol implementation included in SUNDIALS has a unique
identifier specified in enumeration and shown in
:numref:`SUNLinSol.API.IDs`. It is recommended that a
user-supplied SUNLinSol implementation use the
``SUNLINEARSOLVER_CUSTOM`` identifier.

.. _SUNLinSol.API.IDs:
.. table:: Identifiers associated with :c:type:`SUNLinearSolver`
           modules supplied with SUNDIALS
   :align: center

   ==================================  ===================================================  ========
   SUNLinSol ID                        Linear solver type                                   ID Value
   ==================================  ===================================================  ========
   SUNLINEARSOLVER_BAND                Banded direct linear solver (internal)               0
   SUNLINEARSOLVER_DENSE               Dense direct linear solver (internal)                1
   SUNLINEARSOLVER_KLU                 Sparse direct linear solver (KLU)                    2
   SUNLINEARSOLVER_LAPACKBAND          Banded direct linear solver (LAPACK)                 3
   SUNLINEARSOLVER_LAPACKDENSE         Dense direct linear solver (LAPACK)                  4
   SUNLINEARSOLVER_PCG                 Preconditioned conjugate gradient iterative solver   5
   SUNLINEARSOLVER_SPBCGS              Scaled-preconditioned BiCGStab iterative solver      6
   SUNLINEARSOLVER_SPFGMR              Scaled-preconditioned FGMRES iterative solver        7
   SUNLINEARSOLVER_SPGMR               Scaled-preconditioned GMRES iterative solver         8
   SUNLINEARSOLVER_SPTFQMR             Scaled-preconditioned TFQMR iterative solver         9
   SUNLINEARSOLVER_SUPERLUDIST         Parallel sparse direct linear solver (SuperLU_Dist)  10
   SUNLINEARSOLVER_SUPERLUMT           Threaded sparse direct linear solver (SuperLU_MT)    11
   SUNLINEARSOLVER_CUSOLVERSP_BATCHQR  Sparse direct linear solver (CUDA)                   12
   SUNLINEARSOLVER_MAGMADENSE          Dense or block-dense direct linear solver (MAGMA)    13
   SUNLINEARSOLVER_ONEMKLDENSE         Dense or block-dense direct linear solver (OneMKL)   14
   SUNLINEARSOLVER_CUSTOM              User-provided custom linear solver                   15
   ==================================  ===================================================  ========


.. _SUNLinSol.Intended:


Intended use cases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The SUNLinSol and SUNMATRIX APIs are designed to require a minimal set
of routines to ease interfacing with custom or third-party linear solver
libraries. Many external solvers provide routines with similar functionality
and thus may require minimal effort to wrap within custom SUNMATRIX and
SUNLinSol implementations. As SUNDIALS packages utilize generic
SUNLinSol modules they may naturally leverage user-supplied
``SUNLinearSolver`` implementations, thus there exist a wide range of
possible linear solver combinations. Some intended use cases for both the
SUNDIALS-provided and user-supplied SUNLinSol modules are discussd in the
sections below.


Direct linear solvers
""""""""""""""""""""""""""""""""

Direct linear solver modules require a matrix and compute an "exact" solution to
the linear system *defined by the matrix*.  SUNDIALS packages strive to
amortize the high cost of matrix construction by reusing matrix information for
multiple nonlinear iterations or time steps. As a result, each package's linear
solver interface recomputes matrix information as infrequently as possible.

Alternative matrix storage formats and compatible linear solvers that are not
currently provided by, or interfaced with, SUNDIALS can leverage this
infrastructure with minimal effort. To do so, a user must implement custom
SUNMATRIX and SUNLinSol wrappers for the desired matrix format and/or linear
solver following the APIs described in :numref:`SUNMatrix`
and :numref:`SUNLinSol`.  *This user-supplied SUNLinSol module must then
self-identify as having* ``SUNLINEARSOLVER_DIRECT`` *type*.


Matrix-free iterative linear solvers
""""""""""""""""""""""""""""""""""""""

Matrix-free iterative linear solver modules do not require a matrix, and instead
compute an inexact solution to the linear system *defined by the
package-supplied* ``ATimes`` *routine*. SUNDIALS supplies multiple scaled,
preconditioned iterative SUNLinSol modules that support scaling, allowing
packages to handle non-dimensionalization, and users to define variables and
equations as natural in their applications. However, for linear solvers that do
not support left/right scaling, SUNDIALS packages must instead adjust the
tolerance supplied to the linear solver to compensate (see the iterative linear
tolerance section that follows for more details) -- this strategy may be
non-optimal since it cannot handle situations where the magnitudes of different
solution components or equations vary dramatically within a single application.

To utilize alternative linear solvers that are not currently provided by, or
interfaced with, SUNDIALS a user must implement a custom SUNLinSol wrapper
for the linear solver following the API described in
:numref:`SUNLinSol`.  *This user-supplied SUNLinSol module must then
self-identify as having* ``SUNLINEARSOLVER_ITERATIVE`` *type*.


Matrix-based iterative linear solvers (reusing :math:`A`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Matrix-based iterative linear solver modules require a matrix and compute an
inexact solution to the linear system *defined by the matrix*.  This
matrix will be updated infrequently and resued across multiple solves
to amortize the cost of matrix construction. As in the direct linear
solver case, only thin SUNMATRIX and SUNLinSol wrappers for the underlying
matrix and linear solver structures need to be created to utilize
such a linear solver. *This user-supplied SUNLinSol module must then
self-identify as having* ``SUNLINEARSOLVER_MATRIX_ITERATIVE`` *type*.

At present, SUNDIALS has one example problem that uses this approach for
wrapping a structured-grid matrix, linear solver, and preconditioner from the
*hypre* library; this may be used as a template for other customized
implementations (see ``examples/arkode/CXX_parhyp/ark_heat2D_hypre.cpp``).


Matrix-based iterative linear solvers (current :math:`A`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For users who wish to utilize a matrix-based iterative linear solver
where the matrix is *purely for preconditioning* and the linear system is
*defined by the package-supplied* ``ATimes`` *routine*, we envision two
current possibilities.

The preferred approach is for users to employ one of the SUNDIALS
scaled, preconditioned iterative linear solver implementations
(:c:func:`SUNLinSol_SPGMR`, :c:func:`SUNLinSol_SPFGMR`,
:c:func:`SUNLinSol_SPBCGS`, :c:func:`SUNLinSol_SPTFQMR`, or
:c:func:`SUNLinSol_PCG`) as the outer solver. The creation and storage of the
preconditioner matrix, and interfacing with the corresponding matrix-based
linear solver, can be handled through a package's preconditioner "setup" and
"solve" functionality without creating SUNMATRIX and SUNLinSol implementations.
This usage mode is recommended primarily because the SUNDIALS-provided modules
support variable and equation scaling as described above.

A second approach supported by the linear solver APIs is as follows. If the
SUNLinSol implementation is matrix-based, *self-identifies
as having* ``SUNLINEARSOLVER_ITERATIVE`` *type*, and *also provides a non-NULL*
:c:func:`SUNLinSolSetATimes` *routine*, then each SUNDIALS package
will call that routine to attach its package-specific matrix-vector
product routine to the SUNLinSol object. The SUNDIALS package will
then call the SUNLinSol-provided :c:func:`SUNLinSolSetup()` routine
(infrequently) to update matrix information, but will provide current
matrix-vector products to the SUNLinSol implementation through the
package-supplied ``SUNATimesFn`` routine.


Application-specific linear solvers with embedded matrix structure
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Many applications can exploit additional linear system structure arising
from to the implicit couplings in their model equations.  In certain
circumstances, the linear solve :math:`Ax=b` may be performed without
the need for a global system matrix :math:`A`, as the unformed :math:`A`
may be block diagonal or block triangular, and thus the overall linear
solve may be performed through a sequence of smaller linear solves.
In other circumstances, a linear system solve may be accomplished via
specialized fast solvers, such as the fast Fourier transform, fast
multipole method, or treecode, in which case no matrix structure
may be explicitly necessary.  In many of the above situations,
construction and preprocessing of the linear system matrix :math:`A` may be
inexpensive, and thus increased performance may be possible if the current
linear system information is used within every solve (instead of being lagged,
as occurs with matrix-based solvers that reuse :math:`A`).

To support such application-specific situations, SUNDIALS supports user-provided
linear solvers with the ``SUNLINEARSOLVER_MATRIX_EMBEDDED`` type.  For an
application to leverage this support, it should define a custom SUNLinSol
implementation having this type, that only needs to implement the required
:c:func:`SUNLinSolGetType` and :c:func:`SUNLinSolSolve` operations.
Within :c:func:`SUNLinSolSolve`, the linear solver implementation
should call package-specific interface routines (e.g.,
``ARKStepGetNonlinearSystemData``, ``CVodeGetNonlinearSystemData``,
``IDAGetNonlinearSystemData``, ``ARKStepGetCurrentGamma``,
``CVodeGetCurrentGamma``, ``IDAGetCurrentCj``, or
``MRIStepGetCurrentGamma``) to construct the relevant system matrix
:math:`A` (or portions thereof), solve the linear system :math:`Ax=b`, and
return the solution vector :math:`x`.

We note that when attaching this custom SUNLinearSolver object with the relevant
SUNDIALS package ``SetLinearSolver`` routine, the input :c:type:`SUNMatrix`
``A`` should be set to ``NULL``.

For templates of such user-provided "matrix-embedded" SUNLinSol implementations,
see the SUNDIALS examples ``ark_analytic_mels.c``, ``cvAnalytic_mels.c``,
``cvsAnalytic_mels.c``, ``idaAnalytic_mels.c``, and ``idasAnalytic_mels.c``.
