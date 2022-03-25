.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _KINSOL.Usage.CC:

**************************************************
Using KINSOL for the Solution of Nonlinear Systems
**************************************************

This section is concerned with the use of KINSOL for the solution of nonlinear
systems.

The following sections treat the header files and the layout of the user’s main
program, and provide descriptions of the KINSOL user-callable functions and
user-supplied functions. The sample programs described in the companion document
:cite:p:`kinsol_ex` may also be helpful. Those codes may be used as templates (with
the removal of some lines used in testing) and are included in the KINSOL package.

KINSOL uses various constants for both input and output. These are defined as
needed in this chapter, but for convenience are also listed separately in
:numref:`KINSOL.Constants`.

The user should be aware that not all ``SUNLinearSolver`` and ``SUNMatrix``
objects are compatible with all ``N_Vector`` implementations. Details on
compatibility are given in the documentation for each ``SUNMatrix`` (Chapter
:numref:`SUNMatrix`) and ``SUNLinearSolver`` (Chapter :numref:`SUNLinSol`)
implementation. For example, ``NVECTOR_PARALLEL`` is not compatible with the
dense, banded, or sparse ``SUNMatrix`` types, or with the corresponding dense,
banded, or sparse ``SUNLinearSolver`` objects. Please check Chapters
:numref:`SUNMatrix` and :numref:`SUNLinSol` to verify compatibility between
these objects. In addition to that documentation, we note that the KINBBDPRE
preconditioner can only be used with ``NVECTOR_PARALLEL``. It is not recommended
to use a threaded vector object with SuperLU_MT unless it is the
``NVECTOR_OPENMP`` module, and SuperLU_MT is also compiled with OpenMP.

.. _KINSOL.Usage.CC.file_access:

Access to library and header files
----------------------------------

At this point, it is assumed that the installation of KINSOL, following the
procedure described in :numref:`Installation`, has been completed successfully.

Regardless of where the user’s application program resides, its associated
compilation and load commands must make reference to the appropriate locations
for the library and header files required by KINSOL. The relevant library files are

.. code-block::

  <libdir>/libsundials_kinsol.<so|a>
  <libdir>/libsundials_nvec*.<so|a>
  <libdir>/libsundials_sunmat*.<so|a>
  <libdir>/libsundials_sunlinsol*.<so|a>
  <libdir>/libsundials_sunnonlinsol*.<so|a>

where the file extension ``.so`` is typically for shared libraries and ``.a``
for static libraries. The relevant header files are located in the
subdirectories

.. code-block::

  <incdir>/kinsol
  <incdir>/sundials
  <incdir>/nvector
  <incdir>/sunmatrix
  <incdir>/sunlinsol
  <incdir>/sunnonlinsol

The directories ``libdir`` and ``incdir`` are the install library and include
directories, respectively. For a default installation, these are
``<instdir>/lib`` or ``<instdir>/lib64`` and ``<instdir>/include``,
respectively, where ``instdir`` is the directory where SUNDIALS was installed
(see :numref:`Installation`).


.. include:: ../../../../shared/Types.rst


.. _KINSOL.Usage.CC.header_sim:

Header files
------------

The calling program must include several header files so that various macros and
data types can be used. The header file that is always required is:

* ``kinsol/kinsol.h`` the main header file for kinsol, which defines the types and
  various constants, and includes function prototypes. This includes the
  header file for KINLS, ``kinsol/kinsol_ls.h``.

Note that ``kinsol.h`` includes ``sundials_types.h``, which defines the types,
``realtype``, ``sunindextype``, and ``booleantype`` and the constants
``SUNFALSE`` and ``SUNTRUE``.

The calling program must also include an ``N_Vector`` implementation
header file, of the form ``nvector/nvector_*.h`` (see :numref:`NVectors`
for more information). This file in turn includes the header file
``sundials_nvector.h`` which defines the abstract vector data type.

If using a Newton or Picard nonlinear solver that requires the solution of a
linear system, then a linear solver module header file will be required. If the
linear solver is matrix-based, the linear solver header will also include a
header file of the from ``sunmatrix/sunmatrix_*.h`` where ``*`` is the name of
the matrix implementation compatible with the linear solver. The matrix header
file provides access to the relevant matrix functions/macros and in turn
includes the header file ``sundials_matrix.h`` which defines the abstract matrix
data type.

Other headers may be needed, according to the choice of preconditioner, etc. For
example, in the example ``kinFoodWeb_kry_p`` (see :cite:p:`kinsol_ex`),
preconditioning is done with a block-diagonal matrix. For this, even though the
``SUNLINSOL_SPGMR`` linear solver is used, the header
``sundials/sundials_dense.h`` is included for access to the underlying generic
dense matrix arithmetic routines.

.. _KINSOL.Usage.CC.skeleton_sim:

A skeleton of the user’s main program
-------------------------------------

The following is a skeleton of the user’s main program (or calling program) for
the solution of a nonlinear system problem.. Most of the steps are independent
of the ``N_Vector``, ``SUNMatrix``, and ``SUNLinearSolver`` implementations
used. For the steps that are not, refer to :numref:`NVectors`,
:numref:`SUNMatrix`, and :numref:`SUNLinSol` for the specific name of the
function to be called or macro to be referenced.

#. **Initialize parallel or multi-threaded environment** (*if appropriate*)

   For example, call ``MPI_Init`` to initialize MPI if used.

#. **Create the SUNDIALS context object**

   Call :c:func:`SUNContext_Create` to allocate the ``SUNContext`` object.

#. **Set the problem dimensions etc.**

   This generally includes the problem size ``N``, and may include the local
   vector length `Nlocal`.

#. **Create the vector with the initial guess**

   Construct an ``N_Vector`` of initial guess values using the appropriate functions
   defined by the particular ``N_Vector`` implementation (see
   :numref:`NVectors` for details).

   For native SUNDIALS vector implementations, use a call of the form
   ``y0 = N_VMake_***(..., ydata)`` if the array containing the initial values
   of :math:`y` already exists. Otherwise, create a new vector by making a call
   of the form ``N_VNew_***(...)``, and then set its elements by accessing the
   underlying data with a call of the form ``ydata = N_VGetArrayPointer(y0)``.
   Here, ``***`` is the name of the vector implementation.

   For *hypre*, PETSc, and Trilinos vector wrappers, first create and initialize
   the underlying vector, and then create an ``N_Vector`` wrapper with a call
   of the form ``y0 = N_VMake_***(yvec)``, where ``yvec`` is a *hypre*, PETSc,
   or Trilinos vector.  Note that calls like ``N_VNew_***(...)`` and
   ``N_VGetArrayPointer(...)`` are not available for these vector wrappers.

#. **Create matrix object** (*if appropriate*)

   If a linear solver is required (e.g., when using the default Newton solver)
   and the linear solver will be a matrix-based linear solver, then a template
   Jacobian matrix must be created by calling the appropriate constructor
   defined by the particular ``SUNMatrix`` implementation.

   For the native SUNDIALS ``SUNMatrix`` implementations, the matrix object may
   be created using a call of the form ``SUN***Matrix(...)`` where ``***`` is
   the name of the matrix (see :numref:`SUNMatrix` for details).

#. **Create linear solver object** (*if appropriate*)

   If a linear solver is required (e.g., when using the default Newton solver),
   then the desired linear solver object must be created by calling the
   appropriate constructor defined by the particular ``SUNLinearSolver``
   implementation.

   For any of the native SUNDIALS ``SUNLinearSolver`` implementations, the
   linear solver object may be created using a call of the form
   ``SUNLinearSolver LS = SUNLinSol_***(...);`` where ``***`` is the name of
   the linear solver (see :numref:`SUNLinSol` for details).

#. **Create KINSOL object**

   Call :c:func:`KINCreate` to create the KINSOL solver object.

#. **Initialize KINSOL solver**

   Call :c:func:`KINInit` to allocate internal memory.

#. **Attach the linear solver** (*if appropriate*)

   If a linear solver was created above, initialize the KINLS linear solver
   interface by attaching the linear solver object (and matrix object,
   if applicable) with :c:func:`KINSetLinearSolver`.

#. **Set linear solver optional inputs** (*if appropriate*)

   See :numref:`KINSOL.Usage.CC.optional_input.Table` for KINLS optional inputs
   and Chapter :numref:`SUNLinSol` for linear solver specific optional inputs.

#. **Set optional inputs**

   Call ``KINSet***`` functions to change any optional inputs that control the
   behavior of KINSOL from their default values. See :numref:`KINSOL.Usage.CC.optional_input` for details.

#. **Solve problem**

   Call ``ier = KINSol(...)`` to solve the nonlinear problem for a given
   initial guess.

   See :c:func:`KINSol` for details.

#. **Get optional outputs**

   Call ``KINGet***`` functions to obtain optional output. See
   :numref:`KINSOL.Usage.CC.optional_output` for details.

#. **Deallocate memory**

   Upon completion of the integration call the following, as necessary, to free
   any objects or memory allocated above:

   * Call :c:func:`N_VDestroy` to free vector objects.
   * Call :c:func:`SUNMatDestroy` to free matrix objects.
   * Call :c:func:`SUNLinSolFree` to free linear solvers objects.
   * Call :c:func:`SUNNonlinSolFree` to free nonlinear solvers objects.
   * Call :c:func:`KINFree` to free the memory allocated by KINSOL.
   * Call :c:func:`SUNContext_Free` to free the ``SUNContext`` object

#. **Finalize MPI, if used**

   Call ``MPI_Finalize`` to terminate MPI.


.. _KINSOL.Usage.CC.callable_fct_sim:

User-callable functions
-----------------------

This section describes the KINSOL functions that are called by the user to setup
and then solve an IVP. Some of these are required.  However, starting with
:numref:`KINSOL.Usage.CC.optional_input`, the functions listed involve optional
inputs/outputs or restarting, and those paragraphs may be skipped for a casual
use of KINSOL. In any case, refer to :numref:`KINSOL.Usage.CC.skeleton_sim` for the
correct order of these calls.

On an error, each user-callable function returns a negative value and sends an
error message to the error handler routine, which prints the message on
``stderr`` by default. However, the user can set a file as error output or can
provide his own error handler function (see :numref:`KINSOL.Usage.CC.optional_input`).

.. _KINSOL.Usage.CC.callable_fct_sim.kinmalloc:

KINSOL initialization and deallocation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. c:function:: void KINCreate(SUNContext sunctx)

   The function :c:func:`KINCreate` instantiates a KINSOL solver object.

   **Arguments:**
     - ``sunctx`` -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
     * ``void``


.. c:function:: int KINInit(void * kin_mem, KINSysFn func, N_Vector tmpl)

   The function :c:func:`KINInit` specifies the problem-defining  function,
   allocates internal memory, and initializes KINSOL.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block returned by :c:func:`KINCreate`.
     * ``func`` -- is the CC function which computes the system function :math:`F(u)` (or :math:`G(u)` for fixed-point iteration) in the nonlinear problem. This function has the form ``func(u, fval, user_data)``. (For full details see :numref:`KINSOL.Usage.CC.user_fct_sim.resFn`).
     * ``tmpl`` -- is any ``N_Vector`` (e.g. the initial guess vector ``u``) which is used as a template to create (by cloning) necessary vectors in ``kin_mem``.

   **Return value:**
     * ``KIN_SUCCESS`` -- The call to :c:func:`KINInit` was successful.
     * ``KIN_MEM_NULL`` -- The KINSOL memory block was not initialized through a previous call to :c:func:`KINCreate`.
     * ``KIN_MEM_FAIL`` -- A memory allocation request has failed.
     * ``KIN_ILL_INPUT`` -- An input argument to :c:func:`KINInit` has an illegal value.

   **Notes:**
      If an error occurred, :c:func:`KINInit` sends an error message to the
      error handler function.


.. c:function:: void KINFree(void** kin_mem)

   The function :c:func:`KINFree` frees the pointer allocated by a previous call
   to :c:func:`KINCreate`.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.

   **Return value:**
      * ``void``


.. _KINSOL.Usage.CC.callable_fct_sim.lin_solv_init:

Linear solver specification functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As previously explained, Newton and Picard iterations require the solution of
linear systems of the form :math:`J\delta = -F`. Solution of these linear
systems is handled using the KINLS linear solver interface. This interface
supports all valid ``SUNLinearSolver`` modules. Here, matrix-based
``SUNLinearSolver`` modules utilize ``SUNMatrix`` objects to store the Jacobian
matrix :math:`J = F'(u)` and factorizations used throughout
the solution process. Conversely, matrix-free ``SUNLinearSolver`` modules
instead use iterative methods to solve the linear systems of equations, and only
require the *action* of the Jacobian on a vector, :math:`Jv`.

With most iterative linear solvers, preconditioning can be done on the left
only, on the right only, on both the left and the right, or not at all. However,
only right preconditioning is supported within KINLS. If preconditioning
is done, user-supplied functions define the linear operator corresponding to a
right preconditioner matrix :math:`P`, which should approximate the system
Jacobian matrix :math:`J`. For the specification of a preconditioner, see the
iterative linear solver sections in :numref:`KINSOL.Usage.CC.optional_input`
and :numref:`KINSOL.Usage.CC.user_fct_sim`. A
preconditioner matrix :math:`P` must approximate the Jacobian :math:`J`, at
least crudely.

To specify a generic linear solver to KINSOL, after the call to :c:func:`KINCreate`
but before any calls to :c:func:`KINSol`, the user’s program must create the
appropriate ``SUNLinearSolver`` object and call the function
:c:func:`KINSetLinearSolver`, as documented below. To create the ``SUNLinearSolver``
object, the user may call one of the SUNDIALS-packaged ``SUNLinearSolver``
module constructor routines via a call of the form

.. code-block:: c

         SUNLinearSolver LS = SUNLinSol_*(...);

For a current list of such constructor routines see :numref:`SUNLinSol`.

Alternately, a user-supplied ``SUNLinearSolver`` module may be created and used
instead. The use of each of the generic linear solvers involves certain
constants, functions and possibly some macros, that are likely to be needed in
the user code. These are available in the corresponding header file associated
with the specific ``SUNMatrix`` or ``SUNLinearSolver`` module in question, as
described in Chapters :numref:`SUNMatrix` and :numref:`SUNLinSol`.

Once this solver object has been constructed, the user should attach it to
KINSOL via a call to :c:func:`KINSetLinearSolver`. The first argument passed to this
function is the KINSOL memory pointer returned by :c:func:`KINCreate`; the second
argument is the desired ``SUNLinearSolver`` object to use for solving Newton or
Picard systems. The third argument is an optional ``SUNMatrix`` object to
accompany matrix-based ``SUNLinearSolver`` inputs (for matrix-free linear
solvers, the third argument should be ``NULL``). A call to this function
initializes the KINLS linear solver interface, linking it to the main
KINSOL solver, and allows the user to specify additional parameters and routines
pertinent to their choice of linear solver.


.. c:function:: int KINSetLinearSolver(void * kin_mem, SUNLinearSolver LS, SUNMatrix J)

   The function :c:func:`KINSetLinearSolver` attaches a generic ``SUNLinSol``
   object ``LS`` and corresponding template Jacobian ``SUNMatrix``  object ``J``
   (if applicable) to KINSOL, initializing the  KINLS linear solver interface.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``LS`` -- SUNLINSOL object to use for solving Newton linear systems.
     * ``J`` -- SUNMATRIX object for used as a template for the Jacobian (or ``NULL`` if not applicable).

   **Return value:**
     * ``KINLS_SUCCESS`` -- The KINLS initialization was successful.
     * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KINLS_ILL_INPUT`` -- The KINLS interface is not compatible with the ``LS`` or ``J`` input objects or is incompatible with the current NVECTOR module.
     * ``KINLS_SUNLS_FAIL`` -- A call to the ``LS`` object failed.
     * ``KINLS_MEM_FAIL`` -- A memory allocation request failed.

   **Notes:**
      If ``LS`` is a matrix-based linear solver, then the template  Jacobian
      matrix ``J`` will be used in the solve process, so if  additional storage
      is required within the ``SUNMatrix`` object  (e.g. for factorization of a
      banded matrix), ensure that the input  object is allocated with sufficient
      size (see the documentation of  the particular ``SUNMatrix`` type in
      Chapter :numref:`SUNMatrix` for  further information).

      The previous routines :c:func:`KINDlsSetLinearSolver` and
      :c:func:`KINSpilsSetLinearSolver` are
      now wrappers for this routine, and may  still be used for
      backward-compatibility.  However, these will be  deprecated in future
      releases, so we recommend that users transition  to the new routine name
      soon.

.. _KINSOL.Usage.CC.kin:

KINSOL solver function
~~~~~~~~~~~~~~~~~~~~~~

This is the central step in the solution process, the call to solve the
nonlinear algebraic system.

.. c:function:: int KINSol(void * kin_mem, N_Vector u, int strategy, N_Vector u_scale, N_Vector f_scale)

   The function :c:func:`KINSol` computes an approximate solution to the
   nonlinear  system.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``u`` -- vector set to initial guess by user before calling :c:func:`KINSol` , but which upon return contains an approximate solution of the nonlinear system :math:`F(u) = 0`.
     * ``strategy`` -- strategy used to solve the nonlinear system. It must be of the following:

       - ``KIN_NONE`` basic Newton iteration
       - ``KIN_LINESEARCH`` Newton with globalization
       - ``KIN_FP`` fixed-point iteration with Anderson Acceleration (no linear solver needed)
       - ``KIN_PICARD`` Picard iteration with Anderson Acceleration (uses a linear solver)

     * ``u_scale`` -- vector containing diagonal elements of scaling matrix :math:`D_u` for vector ``u`` chosen so that the components of :math:`D_u\ u` (as a matrix multiplication) all have roughly the same magnitude when ``u`` is close to a root of :math:`F(u)`.
     * ``f_scale`` -- vector containing diagonal elements of scaling matrix :math:`D_F` for :math:`F(u)` chosen so that the components of :math:`D_F\ F(u)` (as a matrix multiplication) all have roughly the same magnitude when ``u`` is not too near a root of :math:`F(u)`. In the case of a fixed-point iteration, consider :math:`F(u) = G(u) - u`.

   **Return value:**
     * ``KIN_SUCCESS`` -- :c:func:`KINSol` succeeded; the scaled norm of :math:`F(u)` is less than ``fnormtol``.
     * ``KIN_INITIAL_GUESS_OK`` -- The guess ``u`` :math:`=u_0` satisfied the system :math:`F(u)=0` within the tolerances specified (the scaled norm of :math:`F(u_0)` is less than ``0.01*fnormtol``).
     * ``KIN_STEP_LT_STPTOL`` -- KINSOL stopped based on scaled step length. This means that the current iterate may be an approximate solution of the given nonlinear system, but it is also quite possible that the algorithm is "stalled" (making insufficient progress) near an invalid solution, or that the scalar ``scsteptol`` is too large (see :c:func:`KINSetScaledStepTol` in :numref:`KINSOL.Usage.CC.optional_input` to change ``scsteptol`` from its default value).
     * ``KIN_MEM_NULL`` -- The KINSOL memory block pointer was ``NULL``.
     * ``KIN_ILL_INPUT`` -- An input parameter was invalid.
     * ``KIN_NO_MALLOC`` -- The KINSOL memory was not allocated by a call to :c:func:`KINCreate`.
     * ``KIN_MEM_FAIL`` -- A memory allocation failed.
     * ``KIN_LINESEARCH_NONCONV`` -- The line search algorithm was unable to find an iterate sufficiently distinct from the current iterate, or could not find an iterate satisfying the sufficient decrease condition. Failure to satisfy the sufficient decrease condition could mean the current iterate is "close" to an approximate solution of the given nonlinear system, the difference approximation of the matrix-vector product :math:`J(u)\ v` is inaccurate, or the real scalar ``scsteptol`` is too large.
     * ``KIN_MAXITER_REACHED`` --  The maximum number of nonlinear iterations has been reached.
     * ``KIN_MXNEWT_5X_EXCEEDED`` -- Five consecutive steps have been taken that satisfy the inequality
       :math:`\|D_u p\|_{L2} > 0.99\ \texttt{mxnewtstep}` , where :math:`p` denotes the current step and ``mxnewtstep`` is a scalar upper bound on the scaled step length. Such a failure may mean that :math:`\|D_F F(u)\|_{L2}` asymptotes from above to a positive value, or the real scalar ``mxnewtstep`` is too small.
     * ``KIN_LINESEARCH_BCFAIL`` -- The line search algorithm was unable to satisfy the "beta-condition" for ``MXNBCF+1`` nonlinear iterations (not necessarily consecutive), which may indicate the algorithm is making poor progress.
     * ``KIN_LINSOLV_NO_RECOVERY`` -- The user-supplied routine ``psolve`` encountered a recoverable error, but the preconditioner is already current.
     * ``KIN_LINIT_FAIL`` -- The KINLS initialization routine (``linit``) encountered an error.
     * ``KIN_LSETUP_FAIL`` -- The KINLS setup routine (``lsetup``) encountered an error; e.g., the user-supplied routine ``pset`` (used to set up the preconditioner data) encountered an unrecoverable error.
     * ``KIN_LSOLVE_FAIL`` -- The KINLS solve routine (``lsolve``) encountered an error; e.g., the user-supplied routine ``psolve`` (used to to solve the preconditioned linear system) encountered an unrecoverable error.
     * ``KIN_SYSFUNC_FAIL`` -- The system function failed in an unrecoverable manner.
     * ``KIN_FIRST_SYSFUNC_ERR`` -- The system function failed recoverably at the first call.
     * ``KIN_REPTD_SYSFUNC_ERR`` -- The system function had repeated recoverable errors. No recovery is possible.

   **Notes:**
      The components of vectors ``u_scale`` and ``f_scale`` should be strictly
      positive.  ``KIN_SUCCESS=0``, ``KIN_INITIAL_GUESS_OK=1``, and ``KIN_STEP_LT_STPTOL=2``.
      All remaining return values are negative and therefore a test ``flag`` :math:`< 0`  will trap
      all :c:func:`KINSol` failures.


.. _KINSOL.Usage.CC.optional_input:

Optional input functions
~~~~~~~~~~~~~~~~~~~~~~~~

There are numerous optional input parameters that control the behavior of the
KINSOL solver. KINSOL provides functions that can be used to change these from
their default values.  :numref:`KINSOL.Usage.CC.optional_input.Table` lists all
optional input functions in KINSOL which are then described in detail in the
remainder of this section, beginning with those for the main KINSOL solver and
continuing with those for the KINLS linear solver interface.

We note that, on error return, all of these functions also send an error message
to the error handler function. We also note that all error return values are
negative, so a test ``retval`` :math:`<0` will catch any error.

.. _KINSOL.Usage.CC.optional_input.Table:
.. table:: Optional inputs for KINSOL and KINLS

  +--------------------------------------------------------+----------------------------------+------------------------------+
  |                   **Optional input**                   |        **Function name**         |         **Default**          |
  +========================================================+==================================+==============================+
  | **KINSOL main solver**                                 |                                  |                              |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Error handler function                                 | :c:func:`KINSetErrHandlerFn`     | internal fn.                 |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Pointer to an error file                               | :c:func:`KINSetErrFile`          | ``stderr``                   |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Info handler function                                  | :c:func:`KINSetInfoHandlerFn`    | internal fn.                 |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Pointer to an info file                                | :c:func:`KINSetInfoFile`         | ``stdout``                   |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Data for problem-defining function                     | :c:func:`KINSetUserData`         | ``NULL``                     |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Verbosity level of output                              | :c:func:`KINSetPrintLevel`       | 0                            |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Max. number of nonlinear iterations                    | :c:func:`KINSetNumMaxIters`      | 200                          |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | No initial matrix setup                                | :c:func:`KINSetNoInitSetup`      | ``SUNFALSE``                 |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | No residual monitoring\ :math:`{}^{*}`                 | :c:func:`KINSetNoResMon`         | ``SUNFALSE``                 |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Max. iterations without matrix setup                   | :c:func:`KINSetMaxSetupCalls`    | 10                           |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Max. iterations without residual check\ :math:`{}^{*}` | :c:func:`KINSetMaxSubSetupCalls` | 5                            |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Form of :math:`\eta` coefficient                       | :c:func:`KINSetEtaForm`          | ``KIN_ETACHOICE1``           |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Constant value of :math:`\eta`                         | :c:func:`KINSetEtaConstValue`    | 0.1                          |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Values of :math:`\gamma` and :math:`\alpha`            | :c:func:`KINSetEtaParams`        | 0.9 and 2.0                  |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Values of :math:`\omega_{min}` and                     | :c:func:`KINSetResMonParams`     | 0.00001 and 0.9              |
  | :math:`\omega_{max}`\ :math:`{}^{*}`                   |                                  |                              |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Constant value of :math:`\omega`\ :math:`{}^{*}`       | :c:func:`KINSetResMonConstValue` | 0.9                          |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Lower bound on :math:`\epsilon`                        | :c:func:`KINSetNoMinEps`         | ``SUNFALSE``                 |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Max. scaled length of Newton step                      | :c:func:`KINSetMaxNewtonStep`    | :math:`1000|D_u u_0|_2`      |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Max. number of :math:`\beta`-condition failures        | :c:func:`KINSetMaxBetaFails`     | 10                           |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Rel. error for D.Q. :math:`Jv`                         | :c:func:`KINSetRelErrFunc`       | :math:`\sqrt{\text{uround}}` |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Function-norm stopping tolerance                       | :c:func:`KINSetFuncNormTol`      | uround\ :math:`^{1/3}`       |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Scaled-step stopping tolerance                         | :c:func:`KINSetScaledStepTol`    | :math:`\text{uround}^{2/3}`  |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Inequality constraints on solution                     | :c:func:`KINSetConstraints`      | ``NULL``                     |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Nonlinear system function                              | :c:func:`KINSetSysFunc`          | none                         |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Return the newest fixed point iteration                | :c:func:`KINSetReturnNewest`     | ``SUNFALSE``                 |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Fixed point/Picard damping parameter                   | :c:func:`KINSetDamping`          | 1.0                          |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Anderson Acceleration subspace size                    | :c:func:`KINSetMAA`              | 0                            |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Anderson Acceleration damping parameter                | :c:func:`KINSetDampingAA`        | 1.0                          |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Anderson Acceleration delay                            | :c:func:`KINSetDelayAA`          | 0                            |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Anderson Acceleration orthogonalization routine        | :c:func:`KINSetOrthAA`           | ``KIN_ORTH_MGS``             |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | **KINLS linear solver interface**                      |                                  |                              |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Jacobian function                                      | :c:func:`KINSetJacFn`            | DQ                           |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Preconditioner functions and data                      | :c:func:`KINSetPreconditioner`   | ``NULL``, ``NULL``, ``NULL`` |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Jacobian-times-vector function and data                | :c:func:`KINSetJacTimesVecFn`    | internal DQ, ``NULL``        |
  +--------------------------------------------------------+----------------------------------+------------------------------+
  | Jacobian-times-vector system function                  | :c:func:`KINSetJacTimesVecSysFn` | ``NULL``                     |
  +--------------------------------------------------------+----------------------------------+------------------------------+


.. c:function:: int KINSetErrFile(void * kin_mem, FILE * errfp)

   The function :c:func:`KINSetErrFile` specifies the pointer to the file  where
   all KINSOL messages should be directed when the default  KINSOL error handler
   function is used.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``errfp`` -- pointer to output file.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      The default value for ``errfp`` is ``stderr``.

      Passing a value of
      ``NULL`` disables all future error message output  (except for the case in
      which the KINSOL memory pointer is ``NULL``).  This use of
      :c:func:`KINSetErrFile` is strongly discouraged.

   .. warning::
      If :c:func:`KINSetErrFile` is to be called, it should be called before any
      other optional input functions, in order to take effect for any later
      error message.


.. c:function:: int KINSetErrHandlerFn(void * kin_mem, KINErrHandlerFn ehfun, void * eh_data)

   The function :c:func:`KINSetErrHandlerFn` specifies the optional user-defined
   function  to be used in handling error messages.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``ehfun`` -- is the user's CC error handler function (see :numref:`KINSOL.Usage.CC.user_fct_sim.ehFn`).
     * ``eh_data`` -- pointer to user data passed to ``ehfun`` every time it is called.

   **Return value:**
     * ``KIN_SUCCESS`` -- The function ``ehfun`` and data pointer ``eh_data`` have been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      The default internal error handler function directs error messages to the
      file specified by the file pointer ``errfp`` (see :c:func:`KINSetErrFile`
      above).

      Error messages indicating that the KINSOL solver memory is
      ``NULL`` will  always be directed to ``stderr``.


.. c:function:: int KINSetInfoFile(void * kin_mem, FILE * infofp)

   The function :c:func:`KINSetInfoFile` specifies the pointer to the file
   where all informative (non-error) messages should be directed.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``infofp`` -- pointer to output file.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      The default value for ``infofp`` is ``stdout``.


.. c:function:: int KINSetInfoHandlerFn(void * kin_mem, KINInfoHandlerFn ihfun, void * ih_data)

   The function :c:func:`KINSetInfoHandlerFn` specifies the optional
   user-defined function  to be used in handling informative (non-error)
   messages.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``ihfun`` -- is the user's CC information handler function (see :numref:`KINSOL.Usage.CC.user_fct_sim.ihFn`).
     * ``ih_data`` -- pointer to user data passed to ``ihfun`` every time it is called.

   **Return value:**
     * ``KIN_SUCCESS`` -- The function ``ihfun`` and data pointer ``ih_data`` have been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      The default internal information handler function directs informative
      (non-error)  messages to the file specified by the file pointer ``infofp``
      (see  :c:func:`KINSetInfoFile` above).


.. c:function:: int KINSetPrintLevel(void * kin_mem, int printfl)

   The function :c:func:`KINSetPrintLevel` specifies the level of verbosity  of the output.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``printfl`` -- flag indicating the level of verbosity. Must be one of:

       0 -- no information is displayed.

       1 -- for each nonlinear iteration display the following information:

       - the scaled Euclidean :math:`\ell_2` norm of the system function evaluated at the current iterate,
       - the scaled norm of the Newton step (only if using ``KIN_NONE``), and
       - the number of function evaluations performed so far.

       2 -- display level 1 output and the following values for each iteration:

       - :math:`\|F(u)\|_{D_F}` (only for ``KIN_NONE``).
       - :math:`\|F(u)\|_{D_F,\infty}` (for ``KIN_NONE`` and ``KIN_LINESEARCH``).

       3 -- display level 2 output plus

       - additional values used by the global strategy (only if using ``KIN_LINESEARCH``), and
       - statistical information for iterative linear solver modules.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``printfl`` had an illegal value.

   **Notes:**
      The default value for ``printfl`` is :math:`0`.


.. c:function:: int KINSetUserData(void * kin_mem, void * user_data)

   The function :c:func:`KINSetUserData` specifies the pointer to user-defined
   memory  that is to be passed to all user-supplied functions.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``user_data`` -- pointer to the user-defined memory.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      If specified, the pointer to ``user_data`` is passed to all user-supplied
      functions that have it as an argument. Otherwise, a ``NULL`` pointer is
      passed.

   .. warning::
      If ``user_data`` is needed in user linear solver or  preconditioner
      functions, the call to  :c:func:`KINSetUserData` must be made before the
      call to specify the  linear solver module.


.. c:function:: int KINSetNumMaxIters(void * kin_mem, long int mxiter)

   The function :c:func:`KINSetNumMaxIters` specifies the maximum number of
   nonlinear iterations allowed.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``mxiter`` -- maximum number of nonlinear iterations.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The maximum number of iterations was non-positive.

   **Notes:**
      The default value for ``mxiter`` is ``MXITER_DEFAULT`` :math:`=200`.


.. c:function:: int KINSetNoInitSetup(void * kin_mem, booleantype noInitSetup)

   The function :c:func:`KINSetNoInitSetup` specifies whether an initial call
   to the preconditioner or Jacobian setup function should be made or not.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``noInitSetup`` -- flag controlling whether an initial call to the
       preconditioner or Jacobian setup function is made (pass ``SUNFALSE``)
       or not made (pass ``SUNTRUE``).

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      The default value for ``noInitSetup`` is ``SUNFALSE``, meaning that an
      initial call  to the preconditioner or Jacobian setup function will be
      made.  A call to this function is useful when solving a sequence of
      problems, in which  the final preconditioner or Jacobian value from one
      problem is to be used initially  for the next problem.


.. c:function:: int KINSetNoResMon(void * kin_mem, booleantype noNNIResMon)

   The function :c:func:`KINSetNoResMon` specifies whether or not the nonlinear
   residual monitoring scheme is used to control Jacobian updating

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``noNNIResMon`` -- flag controlling whether residual monitoring is used
       (pass ``SUNFALSE``) or not used (pass ``SUNTRUE``).

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      When using a direct solver, the default value for ``noNNIResMon`` is
      ``SUNFALSE``,  meaning that the nonlinear residual will be monitored.

   .. warning::
      Residual monitoring is only available for use with  matrix-based linear
      solver modules.


.. c:function:: int KINSetMaxSetupCalls(void * kin_mem, long int msbset)

   The function :c:func:`KINSetMaxSetupCalls` specifies the maximum number of
   nonlinear iterations that can be performed between calls to the
   preconditioner or Jacobian setup function.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``msbset`` -- maximum number of nonlinear iterations without a call to the preconditioner or Jacobian setup function. Pass 0 to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``msbset`` was negative.

   **Notes:**
      The default value for ``msbset`` is ``MSBSET_DEFAULT=10``.  The
      value of ``msbset`` should be a multiple of ``msbsetsub`` (see
      :c:func:`KINSetMaxSubSetupCalls`).


.. c:function:: int KINSetMaxSubSetupCalls(void * kin_mem, long int msbsetsub)

   The function :c:func:`KINSetMaxSubSetupCalls` specifies the maximum number of
   nonlinear iterations between checks by the residual monitoring algorithm.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``msbsetsub`` -- maximum number of nonlinear iterations without checking the nonlinear residual. Pass 0 to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``msbsetsub`` was negative.

   **Notes:**
      The default value for ``msbsetsub`` is ``MSBSET_SUB_DEFAULT`` :math:`=5`.
      The value of ``msbset`` (see :c:func:`KINSetMaxSetupCalls`) should be a
      multiple  of ``msbsetsub``.

   .. warning::
      Residual monitoring is only available for use with  matrix-based linear
      solver modules.


.. c:function:: int KINSetEtaForm(void * kin_mem, int etachoice)

   The function :c:func:`KINSetEtaForm` specifies the method for computing  the
   value of the :math:`\eta` coefficient used in the calculation of the  linear
   solver convergence tolerance.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``etachoice`` -- flag indicating the method for computing :math:`\eta`. The value must be one of ``KIN_ETACHOICE1`` , ``KIN_ETACHOICE2`` , or ``KIN_ETACONSTANT`` (see Chapter :numref:`KINSOL.Mathematics` for details).

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``etachoice`` had an illegal value.

   **Notes:**
      The default value for ``etachoice`` is ``KIN_ETACHOICE1``.  When using
      either ``KIN_ETACHOICE1`` or ``KIN_ETACHOICE2`` the safeguard

      .. math::
         \eta_n = \max(\eta_n, \eta_{\text{safe}})

      is applied when :math:`\eta_{\text{safe}} > 0.1`. For ``KIN_ETACHOICE1``

      .. math::
         \eta_{\text{safe}} = \eta_{n-1}^{\frac{1+\sqrt{5}}{2}}

      and for ``KIN_ETACHOICE2``

      .. math::
         \eta_{\text{safe}} = \gamma \eta_{n-1}^\alpha

      where :math:`\gamma` and :math:`\alpha` can be set with
      :c:func:`KINSetEtaParams`.

      The following safeguards are always applied when using either
      ``KIN_ETACHOICE1`` or ``KIN_ETACHOICE2`` so that
      :math:`\eta_{\text{min}} \leq \eta_n \leq\eta_{\text{max}}`:

      .. math::
         \begin{aligned}
            \eta_n &= \max(\eta_n, \eta_{\text{min}}) \\
            \eta_n &= \min(\eta_n, \eta_{\text{max}})
         \end{aligned}

      where :math:`\eta_{\text{min}} = 10^{-4}` and :math:`\eta_{\text{max}} = 0.9`.


.. c:function:: int KINSetEtaConstValue(void * kin_mem, realtype eta)

   The function :c:func:`KINSetEtaConstValue` specifies the constant value  for
   :math:`\eta` in the case  ``etachoice = KIN_ETACONSTANT``.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``eta`` -- constant value for :math:`\eta`. Pass :math:`0.0` to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``eta`` had an illegal value

   **Notes:**
      The default value for ``eta`` is :math:`0.1`.  The legal values are
      :math:`0.0 <` ``eta`` :math:`\le 1.0`.


.. c:function:: int KINSetEtaParams(void * kin_mem, realtype egamma, realtype ealpha)

   The function :c:func:`KINSetEtaParams` specifies the parameters
   :math:`\gamma` and  :math:`\alpha` in the formula for :math:`\eta`, in the
   case ``etachoice = KIN_ETACHOICE2``.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``egamma`` -- value of the :math:`\gamma` parameter. Pass :math:`0.0` to indicate the default.
     * ``ealpha`` -- value of the :math:`\alpha` parameter. Pass :math:`0.0` to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional values have been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- One of the arguments ``egamma`` or ``ealpha`` had an illegal value.

   **Notes:**
      The default values for ``egamma`` and ``ealpha`` are :math:`0.9` and
      :math:`2.0`, respectively.  The legal values are :math:`0.0 <` ``egamma``
      :math:`\le 1.0` and  :math:`1.0<` ``ealpha`` :math:`\le 2.0`.


.. c:function:: int KINSetResMonConstValue(void * kin_mem, realtype omegaconst)

   The function :c:func:`KINSetResMonConstValue` specifies the constant value
   for :math:`\omega` when using residual monitoring.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``omegaconst`` -- constant value for :math:`\omega`. Passing :math:`0.0` results in using Eqn. :eq:`KIN_resmon_omega`.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``omegaconst`` had an illegal value

   **Notes:**
      The default value for ``omegaconst`` is :math:`0.9`.  The legal values are
      :math:`0.0 <` ``omegaconst`` :math:`< 1.0`.


.. c:function:: int KINSetResMonParams(void * kin_mem, realtype omegamin, realtype omegamax)

   The function :c:func:`KINSetResMonParams` specifies the parameters
   :math:`\omega_{min}` and  :math:`\omega_{max}` in the formula :eq:`KIN_resmon_omega` for
   :math:`\omega`.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``omegamin`` -- value of the :math:`\omega_{min}` parameter. Pass :math:`0.0` to indicate the default.
     * ``omegamax`` -- value of the :math:`\omega_{max}` parameter. Pass :math:`0.0` to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional values have been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- One of the arguments ``omegamin`` or ``omegamax`` had an illegal value.

   **Notes:**
      The default values for ``omegamin`` and ``omegamax`` are :math:`0.00001`
      and :math:`0.9`,  respectively.  The legal values are :math:`0.0 <`
      ``omegamin`` :math:`<` ``omegamax`` :math:`< 1.0`.


.. c:function:: int KINSetNoMinEps(void * kin_mem, booleantype noMinEps)

   The function :c:func:`KINSetNoMinEps` specifies a flag that controls whether
   or not  the value of :math:`\epsilon`, the scaled linear residual tolerance,
   is  bounded from below.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``noMinEps`` -- flag controlling the bound on :math:`\epsilon`. If ``SUNFALSE`` is passed the value of :math:`\epsilon` is constrained and if ``SUNTRUE`` is passed then :math:`\epsilon` is not constrained.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      The default value for ``noMinEps`` is ``SUNFALSE``, meaning that a
      positive minimum value, equal to :math:`0.01`*``fnormtol``, is applied to
      :math:`\epsilon` (see :c:func:`KINSetFuncNormTol` below).


.. c:function:: int KINSetMaxNewtonStep(void * kin_mem, realtype mxnewtstep)

   The function :c:func:`KINSetMaxNewtonStep` specifies the maximum allowable
   scaled  length of the Newton step.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``mxnewtstep`` -- maximum scaled step length :math:`(\geq 0.0)`.
       Pass :math:`0.0` to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The input value was negative.

   **Notes:**
      The default value of ``mxnewtstep`` is :math:`1000\, \| u_0 \|_{D_u}`,
      where :math:`u_0` is the initial guess.


.. c:function:: int KINSetMaxBetaFails(void * kin_mem, realtype mxnbcf)

   The function :c:func:`KINSetMaxBetaFails` specifies the maximum number of
   :math:`\beta`-condition failures in the linesearch algorithm.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``mxnbcf`` -- maximum number of :math:`\beta` -condition failures.
       Pass :math:`0.0` to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- ``mxnbcf`` was negative.

   **Notes:**
      The default value of ``mxnbcf`` is ``MXNBCF_DEFAULT`` :math:`=10`.


.. c:function:: int KINSetRelErrFunc(void * kin_mem, realtype relfunc)

   The function :c:func:`KINSetRelErrFunc` specifies the relative error in
   computing :math:`F(u)`, which is used in the difference quotient
   approximation to  the Jacobian matrix [see Eq. :eq:`KIN_sigmaDQ_direct` ] or the Jacobian-vector
   product [see Eq. :eq:`KIN_sigmaDQ_iterative` ]. The value stored is
   :math:`\sqrt{\texttt{relfunc}}`.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``relfunc`` -- relative error in :math:`F(u)` (:math:`\texttt{relfunc} \geq 0.0`).
       Pass :math:`0.0` to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The relative error was negative.

   **Notes:**
      The default value for ``relfunc`` is :math:`U` = unit roundoff.


.. c:function:: int KINSetFuncNormTol(void * kin_mem, realtype fnormtol)

   The function :c:func:`KINSetFuncNormTol` specifies the scalar used as a
   stopping  tolerance on the scaled maximum norm of the system function
   :math:`F(u)`.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``fnormtol`` -- tolerance for stopping based on scaled function norm :math:`(\geq 0.0)`.
       Pass :math:`0.0` to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The tolerance was negative.

   **Notes:**
      The default value for ``fnormtol`` is (unit roundoff) :math:`^{1/3}`.


.. c:function:: int KINSetScaledStepTol(void * kin_mem, realtype scsteptol)

   The function :c:func:`KINSetScaledStepTol` specifies the scalar used  as a
   stopping tolerance on the minimum scaled step length.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``scsteptol`` -- tolerance for stopping based on scaled step length :math:`(\geq 0.0)`.
       Pass :math:`0.0` to indicate the default.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The tolerance was non-positive.

   **Notes:**
      The default value for ``scsteptol`` is (unit roundoff) :math:`^{2/3}`.


.. c:function:: int KINSetConstraints(void * kin_mem, N_Vector constraints)

   The function :c:func:`KINSetConstraints` specifies a vector that defines
   inequality constraints for each component of the solution vector :math:`u`.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``constraints`` -- vector of constraint flags. If ``constraints[i]`` is

       - :math:`0.0` then no constraint is imposed on :math:`u_i`.
       - :math:`1.0` then :math:`u_i` will be constrained to be :math:`u_i \ge 0.0`.
       - :math:`-1.0` then :math:`u_i` will be constrained to be :math:`u_i \le 0.0`.
       - :math:`2.0` then :math:`u_i` will be constrained to be :math:`u_i > 0.0`.
       - :math:`-2.0` then :math:`u_i` will be constrained to be :math:`u_i < 0.0`.


   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The constraint vector contains illegal values.

   **Notes:**
      The presence of a non-``NULL`` constraints vector that is not :math:`0.0`
      in  all components will cause constraint checking to be performed. If a
      ``NULL``  vector is supplied, constraint checking will be disabled.  The
      function creates a private copy of the constraints vector. Consequently,
      the user-supplied vector can be freed after the function call, and  the
      constraints can only be changed by calling this function.


.. c:function:: int KINSetSysFunc(void * kin_mem, KINSysFn func)

   The function :c:func:`KINSetSysFunc` specifies the user-provided function
   that evaluates the nonlinear system function :math:`F(u)` or :math:`G(u)`.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``func`` -- user-supplied function that evaluates :math:`F(u)` (or :math:`G(u)` for fixed-point iteration).

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``func`` was ``NULL``.

   **Notes:**
      The nonlinear system function is initially specified through
      :c:func:`KINInit`.  The option of changing the system function is provided
      for a user who wishes  to solve several problems of the same size but with
      different functions.


.. c:function:: int KINSetReturnNewest(void * kin_mem, booleantype ret_newest)

   The function :c:func:`KINSetReturnNewest` specifies if the fixed point
   iteration  should return the newest iteration or the iteration consistent
   with the last  function evaluation.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``ret_newest`` --  ``SUNTRUE`` – return the newest iteration. ``SUNFALSE`` – return the iteration consistent with the last function evaluation.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      The default value of ``ret_newest`` is ``SUNFALSE``.


.. c:function:: int KINSetDamping(void * kin_mem, realtype beta)

   The function :c:func:`KINSetDamping` specifies the value of the damping
   parameter  in the fixed point or Picard iteration.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``beta`` -- the damping parameter value :math:`0 < beta \leq 1.0`.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``beta`` was zero or negative.

   **Notes:**
      This function sets the damping parameter value, which needs to be greater
      than  zero and less than one if damping is to be used. A value :math:`\geq
      1` disables  damping.  The default value of ``beta`` is 1.0, indicating no
      damping.  To set the damping parameter used in Anderson acceleration see
      :c:func:`KINSetDampingAA`.  With the fixed point iteration the difference
      between successive iterations is  used to determine convergence. As such,
      when damping is enabled, the tolerance  used to stop the fixed point
      iteration is scaled by ``beta`` to account for  the effects of damping. If
      ``beta`` is extremely small (close to zero), this  can lead to an
      excessively tight tolerance.


.. c:function:: int KINSetMAA(void * kin_mem, long int maa)

   The function :c:func:`KINSetMAA` specifies the size of the subspace used with
   Anderson acceleration in conjunction with Picard or fixed-point iteration.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``maa`` -- subspace size for various methods. A value of 0 means no acceleration, while a positive value means acceleration will be done.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``maa`` was negative.

   **Notes:**
      This function sets the subspace size, which needs to be :math:`> 0` if
      Anderson  Acceleration is to be used.  It also allocates additional memory
      necessary for Anderson Acceleration.  The default value of ``maa`` is 0,
      indicating no acceleration.  The value of ``maa``  should always be less
      than ``mxiter``.  This function MUST be called before calling
      :c:func:`KINInit`.  If the user calls the function KINSetNumMaxIters, that
      call should be made  before the call to KINSetMAA, as the latter uses the
      value of ``mxiter``.


.. c:function:: int KINSetDampingAA(void * kin_mem, realtype beta)

   The function :c:func:`KINSetDampingAA` specifies the value of the Anderson
   acceleration damping paramter.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``beta`` -- the damping parameter value :math:`0 < beta \leq 1.0`.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``beta`` was zero or negative.

   **Notes:**
      This function sets the damping parameter value, which needs to be greater
      than  zero and less than one if damping is to be used. A value :math:`\geq
      1` disables  damping.  The default value of ``beta`` is 1.0, indicating no
      damping.  When delaying the start of Anderson acceleration with
      :c:func:`KINSetDelayAA`, use  :c:func:`KINSetDamping` to set the damping
      parameter in the fixed point or Picard  iterations before Anderson
      acceleration begins. When using Anderson  acceleration without delay, the
      value provided to :c:func:`KINSetDampingAA` is  applied to all iterations
      and any value provided to :c:func:`KINSetDamping` is  ignored.


.. c:function:: int KINSetDelayAA(void * kin_mem, long int delay)

   The function :c:func:`KINSetDelayAA` specifies the number of iterations to
   delay  the start of Anderson acceleration.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``delay`` -- the number of iterations to delay Anderson acceleration.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``delay`` was less than zero.

   **Notes:**
      The default value of ``delay`` is 0, indicating no delay.


.. c:function:: int KINSetOrthAA(void* kin_mem, int orthaa)

   The function :c:func:`KINSetOrthAA` specifies the orthogonalization routine
   to be used in the QR factorization portion of Anderson acceleration.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``orthaa`` -- the orthogonalization routine parameter. Can be set to any of
        the following

        * ``KIN_ORTH_MGS`` --  Modified Gram Schmidt (default)
        * ``KIN_ORTH_ICWY`` --  Inverse Compact WY Modified Gram Schmidt
        *  ``KIN_ORTH_CGS2`` -- Classical Gram Schmidt with Reorthogonalization
           (CGS2)
        * ``KIN_ORTH_DCGS2`` --  Classical Gram Schmidt with Delayed
          Reorthogonlization

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KIN_ILL_INPUT`` -- The argument ``orthaa`` was not one of the predefined
       orthogonalization routines defined in KINSOL.

   .. note::

      This function *must* be called before calling :c:func:`KINInit`.

      An example of how to use this function can be found in
      ``examples/kinsol/serial/kinAnalytic_fp.c``


.. _KINSOL.Usage.CC.optional_inputs.optin_ls:

Linear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For matrix-based linear solver modules, the KINLS solver interface needs a
function to compute an approximation to the Jacobian matrix :math:`J(u)`. This
function must be of type :c:type:`KINLsJacFn`. The user can supply a Jacobian
function, or if using the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or
:ref:`SUNMATRIX_BAND <SUNMatrix.Band>` modules for :math:`J` can use the default
internal difference quotient approximation that comes with the KINLS solver. To
specify a user-supplied Jacobian function ``jac``, KINLS provides the function
:c:func:`KINSetJacFn`. The KINLS interface passes the pointer ``user_data`` to
the Jacobian function. This allows the user to create an arbitrary structure
with relevant problem data and access it during the execution of the
user-supplied Jacobian function, without using global data in the program. The
pointer ``user_data`` may be specified through :c:func:`KINSetUserData`.

.. c:function:: int KINSetJacFn(void* kin_mem, KINLsJacFn jac)

   The function :c:func:`KINSetJacFn` specifies the Jacobian approximation function to
   be used for a matrix-based solver within the KINLS interface.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``jac`` -- user-defined Jacobian approximation function. See :c:type:`KINLsJacFn` for more
        details.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver interface has not been initialized.

   **Notes:**
      This function must be called after the KINLS linear solver interface has been
      initialized through a call to :c:func:`KINSetLinearSolver`.  By default,
      KINLS uses an internal difference quotient function for the
      :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` and
      :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` modules.  If ``NULL`` is passed to ``jac``,
      this default function is used.  An error will occur if no ``jac`` is supplied when
      using other matrix types.

   .. warning::

      The previous routine :c:func:`KINDlsSetJacFn` is now a wrapper for this routine,
      and may still be used for backward-compatibility.  However, this will be
      deprecated in future releases, so we recommend that users transition to
      the new routine name soon.


When using matrix-free linear solver modules, the KINLS linear solver
interface requires a function to compute an approximation to the product between
the Jacobian matrix :math:`J(u)` and a vector :math:`v`. The user can supply
his/her own Jacobian-times-vector approximation function, or use the internal
difference quotient approximation that comes with the KINLS solver
interface.

A user-defined Jacobian-vector function must be of type :c:type:`KINLsJacTimesVecFn`
and can be specified through a call to :c:func:`KINLsSetJacTimesVecFn` (see
:numref:`KINSOL.Usage.CC.user_fct_sim.jtimesFn` for specification details). The pointer
``user_data`` received through :c:func:`KINSetUserData` (or a pointer to ``NULL`` if
``user_data`` was not specified) is passed to the Jacobian-times-vector function
``jtimes`` each time it is called. This allows the user to create an arbitrary
structure with relevant problem data and access it during the execution of the
user-supplied functions without using global data in the program.

.. c:function:: int KINSetJacTimesVecFn(void * kin_mem, KINLsJacTimesVecFn jtimes)

   The function :c:func:`KINSetJacTimesVecFn` specifies the Jacobian-vector  product function.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``jtimes`` -- user-defined Jacobian-vector product function.

   **Return value:**
     * ``KINLS_SUCCESS`` -- The optional value has been successfully set.
     * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.
     * ``KINLS_SUNLS_FAIL`` -- An error occurred when setting up the system matrix-times-vector routines in the SUNLINSOL object used by the KINLS interface.

   **Notes:**
      The default is to use an internal difference quotient for  ``jtimes``.  If
      ``NULL`` is passed as ``jtimes``, this default  is used.  This function
      must be called after the KINLS linear solver  interface has been
      initialized through a call to  :c:func:`KINSetLinearSolver`.  The function
      type :c:type:`KINLsJacTimesVecFn` is described in
      :numref:`KINSOL.Usage.CC.user_fct_sim.jtimesFn`.  The previous
      routine :c:func:`KINSpilsSetJacTimesVecFn` is now a wrapper for  this
      routine, and may still be used for backward-compatibility.  However, this
      will be deprecated in future releases, so we recommend  that users
      transition to the new routine name soon.


When using the internal difference quotient the user may optionally supply an
alternative system function for use in the Jacobian-vector product approximation
by calling :c:func:`KINSetJacTimesVecSysFn`. The alternative system function
should compute a suitable (and differentiable) approximation of the system
function provided to :c:func:`KINInit`. For example, as done in
:cite:p:`dorr2010numerical` when solving the nonlinear systems that arise in the
implicit integration of ordinary differential equations, the alternative
function may use lagged values when evaluating a nonlinearity to avoid
differencing a potentially non-differentiable factor.

.. c:function:: int KINSetJacTimesVecSysFn(void * kin_mem, KINSysFn jtimesSysFn)

   The function :c:func:`KINSetJacTimesVecSysFn` specifies an alternative system
   function for use in the internal Jacobian-vector product difference quotient
   approximation.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``jtimesSysFn`` -- is the CC function which computes the alternative system function to use in Jacobian-vector product difference quotient approximations. This function has the form ``func(u, fval, user_data)``. (For full details see :numref:`KINSOL.Usage.CC.user_fct_sim.resFn`.)

   **Return value:**
     * ``KINLS_SUCCESS`` -- The optional value has been successfully set.
     * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
     * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.
     * ``KINLS_ILL_INPUT`` -- The internal difference quotient approximation is disabled.

   **Notes:**
      The default is to use the system function provided to :c:func:`KINInit` in
      the  internal difference quotient. If the input system function is
      ``NULL``,  the default is used.  This function must be called after the
      KINLS linear solver interface  has been initialized through a call to
      :c:func:`KINSetLinearSolver`.


When using an iterative linear solver, the user may supply a preconditioning
operator to aid in solution of the system. This operator consists of two
user-supplied functions, ``psetup`` and ``psolve``, that are supplied to KINLS
using the function :c:func:`KINSetPreconditioner`. The ``psetup`` function
supplied to this routine should handle evaluation and preprocessing of any
Jacobian data needed by the user’s preconditioner solve function, ``psolve``.
Both of these functions are fully specified in :numref:`KINSOL.Usage.CC.user_fct_sim`.
The user data pointer received through
:c:func:`KINSetUserData` (or a pointer to ``NULL`` if user data was not
specified) is passed to the ``psetup`` and ``psolve`` functions. This allows the
user to create an arbitrary structure with relevant problem data and access it
during the execution of the user-supplied preconditioner functions without using
global data in the program.

.. c:function:: int KINSetPreconditioner(void * kin_mem, KINLsPrecSetupFn psetup, KINLsPrecSolveFn psolve)

   The function :c:func:`KINSetPreconditioner` specifies the preconditioner
   setup and solve functions.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``psetup`` -- user-defined function to set up the preconditioner. See
        :c:type:`KINLsPrecSetupFn` for more details. Pass ``NULL`` if no setup is
        necessary.
      * ``psolve`` -- user-defined preconditioner solve function. See
        :c:type:`KINLsPrecSolveFn` for more details.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional values have been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.
      * ``KINLS_SUNLS_FAIL`` -- An error occurred when setting up preconditioning
        in the ``SUNLinearSolver`` object used by the KINLS interface.

   **Notes:**
      The default is ``NULL`` for both arguments (i.e., no preconditioning).
      This function must be called after the KINLS linear solver interface has
      been initialized through a call to :c:func:`KINSetLinearSolver`.

   .. warning::

      The previous routine :c:func:`KINSpilsSetPreconditioner` is now a wrapper for
      this routine, and may still be used for backward-compatibility.  However,
      this will be removed in future releases, so we recommend that users
      transition to the new routine name soon.


.. _KINSOL.Usage.CC.optional_output:

Optional output functions
~~~~~~~~~~~~~~~~~~~~~~~~~

KINSOL provides an extensive list of functions that can be used to obtain solver
performance information. :numref:`KINSOL.Usage.CC.optional_output.Table`
lists all optional output functions in KINSOL, which are then described in
detail in the remainder of this section, beginning with those for the main
KINSOL solver and continuing with those for the KINLS linear solver interface.
Where the name of an output from a linear solver module would otherwise conflict
with the name of an optional output from the main solver, a suffix ``LS`` (for
Linear Solver) has been added here (e.g., ``lenrwLS``).

.. _KINSOL.Usage.CC.optional_output.Table:
.. table:: Optional outputs from KINSOL and KINLS

  ======================================================== ==================================
  **Optional output**                                      **Function name**
  ======================================================== ==================================
  **KINSOL main solver**
  Size of KINSOL real and integer workspaces               :c:func:`KINGetWorkSpace`
  Number of function evaluations                           :c:func:`KINGetNumFuncEvals`
  Number of nonlinear iterations                           :c:func:`KINGetNumNonlinSolvIters`
  Number of :math:`\beta`-condition failures               :c:func:`KINGetNumBetaCondFails`
  Number of backtrack operations                           :c:func:`KINGetNumBacktrackOps`
  Scaled norm of :math:`F`                                 :c:func:`KINGetFuncNorm`
  Scaled norm of the step                                  :c:func:`KINGetStepLength`
  **KINLS linear solver interface**
  Size of real and integer workspaces                      :c:func:`KINGetLinWorkSpace`
  No. of Jacobian evaluations                              :c:func:`KINGetNumJacEvals`
  No. of :math:`F` calls for D.Q. Jacobian[-vector] evals. :c:func:`KINGetNumLinFuncEvals`
  No. of linear iterations                                 :c:func:`KINGetNumLinIters`
  No. of linear convergence failures                       :c:func:`KINGetNumLinConvFails`
  No. of preconditioner evaluations                        :c:func:`KINGetNumPrecEvals`
  No. of preconditioner solves                             :c:func:`KINGetNumPrecSolves`
  No. of Jacobian-vector product evaluations               :c:func:`KINGetNumJtimesEvals`
  Last return from a KINLS function                        :c:func:`KINGetLastLinFlag`
  Name of constant associated with a return flag           :c:func:`KINGetLinReturnFlagName`
  ======================================================== ==================================


.. _KINSOL.Usage.CC.optional_output.optout_main:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

KINSOL provides several user-callable functions that can be used to obtain
different quantities that may be of interest to the user, such as solver
workspace requirements and solver performance statistics. These optional output
functions are described next.

.. c:function:: int KINGetWorkSpace(void * kin_mem, long int lenrw, long int leniw)

   The function :c:func:`KINGetWorkSpace` returns the  KINSOL integer and real
   workspace sizes.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``lenrw`` -- the number of ``realtype`` values in the KINSOL workspace.
     * ``leniw`` -- the number of integer values in the KINSOL workspace.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional output values have been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.

   **Notes:**
      KINSOL solver  In terms of the problem size :math:`N`, the actual size of
      the real workspace  is :math:`17 + 5 N` ``realtype`` words. The real workspace
      is increased by  an additional :math:`N` words if constraint checking is
      enabled (see :c:func:`KINSetConstraints`).

      The actual size of the integer
      workspace (without distinction between ``int``  and ``long int``) is
      :math:`22 + 5 N` (increased by :math:`N` if constraint checking is enabled).


.. c:function:: int KINGetNumFuncEvals(void * kin_mem, long int nfevals)

   The function :c:func:`KINGetNumFuncEvals` returns the number of evaluations
   of the system function.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``nfevals`` -- number of calls to the user-supplied function that evaluates :math:`F(u)`.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional output value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.


.. c:function:: int KINGetNumNonlinSolvIters(void * kin_mem, long int nniters)

   The function :c:func:`KINGetNumNonlinSolvIters` returns the number  of
   nonlinear iterations.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``nniters`` -- number of nonlinear iterations.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional output value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.


.. c:function:: int KINGetNumBetaCondFails(void * kin_mem, long int nbcfails)

   The function :c:func:`KINGetNumBetaCondFails` returns the number  of
   :math:`\beta`-condition failures.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``nbcfails`` -- number of :math:`\beta` -condition failures.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional output value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.


.. c:function:: int KINGetNumBacktrackOps(void * kin_mem, long int nbacktr)

   The function :c:func:`KINGetNumBacktrackOps` returns the number of  backtrack
   operations (step length adjustments) performed by the  line search algorithm.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``nbacktr`` -- number of backtrack operations.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional output value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.


.. c:function:: int KINGetFuncNorm(void * kin_mem, realtype fnorm)

   The function :c:func:`KINGetFuncNorm` returns the scaled Euclidean
   :math:`\ell_2` norm of the  nonlinear system function :math:`F(u)` evaluated
   at the current iterate.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``fnorm`` -- current scaled norm of :math:`F(u)`.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional output value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.


.. c:function:: int KINGetStepLength(void * kin_mem, realtype steplength)

   The function :c:func:`KINGetStepLength` returns the scaled Euclidean
   :math:`\ell_2` norm of  the step used during the previous iteration.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``steplength`` -- scaled norm of the Newton step.

   **Return value:**
     * ``KIN_SUCCESS`` -- The optional output value has been successfully set.
     * ``KIN_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.


.. _KINSOL.Usage.CC.optional_output.optout_ls:

KINLS linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the KINLS modules:

.. c:function:: int KINGetLinWorkSpace(void * kin_mem, long int * lenrwLS, long int * leniwLS)

   The function :c:func:`KINGetLinWorkSpace` returns the sizes of the real and
   integer workspaces used by the KINLS linear solver interface.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``lenrwLS`` -- the number of real values in the KINLS workspace.
      * ``leniwLS`` -- the number of integer values in the KINLS workspace.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   **Notes:**
      The workspace requirements reported by this routine correspond only to memory
      allocated within this interface and to memory allocated by the
      ``SUNLinearSolver`` object attached to it.  The template Jacobian
      matrix allocated by the user outside of KINLS is not included in this report.

   .. warning::

      The previous routines :c:func:`KINDlsGetWorkspace` and :c:func:`KINSpilsGetWorkspace`
      are now deprecated.


.. c:function:: int KINGetNumJacEvals(void * kin_mem, long int * njevals)

   The function :c:func:`KINGetNumJacEvals` returns the cumulative number of
   calls to the KINLS Jacobian approximation function.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``njevals`` -- the cumulative number of calls to the Jacobian function total so far.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   .. warning::

      The previous routine :c:func:`KINDlsGetNumJacEvals` is now deprecated,


.. c:function:: int KINGetNumLinFuncEvals(void * kin_mem, long int * nrevalsLS)

   The function :c:func:`KINGetNumLinResEvals` returns the cumulative number of
   calls to the user residual function due to the finite difference Jacobian
   approximation or finite difference Jacobian-vector product approximation.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``nrevalsLS`` -- the cumulative number of calls to the user residual
        function.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   **Notes:**
      The value ``nrevalsLS`` is incremented only if one of the default internal
      difference quotient functions is used.

   .. warning::
      The previous routines :c:func:`KINDlsGetNumRhsEvals` and
      :c:func:`KINSpilsGetNumRhsEvals` are now deprecated.


.. c:function:: int KINGetNumLinIters(void * kin_mem, long int * nliters)

   The function :c:func:`KINGetNumLinIters` returns the cumulative number of
   linear iterations.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``nliters`` -- the current number of linear iterations.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   .. warning::
      The previous routine :c:func:`KINSpilsGetNumLinIters` is now deprecated.


.. c:function:: int KINGetNumLinConvFails(void * kin_mem, long int * nlcfails)

   The function :c:func:`KINGetNumLinConvFails` returns the cumulative number of
   linear convergence failures.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``nlcfails`` -- the current number of linear convergence failures.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   .. warning::
      The previous routine :c:func:`KINSpilsGetNumConvFails` is now deprecated.


.. c:function:: int KINGetNumPrecEvals(void * kin_mem, long int * npevals)

   The function :c:func:`KINGetNumPrecEvals` returns the cumulative number of
   preconditioner evaluations, i.e., the number of calls made to ``psetup``.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``npevals`` -- the cumulative number of calls to ``psetup``.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   .. warning::

      The previous routine :c:func:`KINSpilsGetNumPrecEvals` is now deprecated.


.. c:function:: int KINGetNumPrecSolves(void * kin_mem, long int * npsolves)

   The function :c:func:`KINGetNumPrecSolves` returns the cumulative number of
   calls made to the preconditioner solve function, ``psolve``.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``npsolves`` -- the cumulative number of calls to ``psolve``.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   .. warning::

      The previous routine :c:func:`KINSpilsGetNumPrecSolves` is now deprecated.


.. c:function:: int KINGetNumJtimesEvals(void * kin_mem, long int * njvevals)

   The function :c:func:`KINGetNumJtimesEvals` returns the cumulative number of
   calls made to the Jacobian-vector product function, ``jtimes``.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``njvevals`` -- the cumulative number of calls to ``jtimes``.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   .. warning::

      The previous routine :c:func:`KINSpilsGetNumJtimesEvals` is now deprecated.


.. c:function:: int KINGetLastLinFlag(void * kin_mem, long int * lsflag)

   The function :c:func:`KINGetLastLinFlag` returns the last return value from
   an KINLS routine.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``lsflag`` -- the value of the last return flag from an KINLS function.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer is ``NULL``.
      * ``KINLS_LMEM_NULL`` -- The KINLS linear solver has not been initialized.

   **Notes:**
      If the KINLS setup function failed (i.e., :c:func:`KINSolve` returned
      ``KIN_LSETUP_FAIL``) when using the :ref:`SUNLINSOL_DENSE <SUNLinSol_Dense>`
      or :ref:`SUNLINSOL_BAND <SUNLinSol_Band>` modules, then the value of
      ``lsflag`` is equal to the column index (numbered from one) at which a zero
      diagonal element was encountered during the LU factorization of the (dense or
      banded) Jacobian matrix.

      If the KINLS setup function failed when using another ``SUNLinearSolver``
      object, then ``lsflag`` will be ``SUNLS_PSET_FAIL_UNREC``,
      ``SUNLS_ASET_FAIL_UNREC``, or ``SUNLS_PACKAGE_FAIL_UNREC``.

      If the KINLS solve function failed (:c:func:`KINSolve` returned ``KIN_LSOLVE_FAIL``),
      ``lsflag`` contains the error return flag from the ``SUNLinearSolver``
      object, which will be one of: ``SUNLS_MEM_NULL``, indicating that the
      ``SUNLinearSolver`` memory is ``NULL``; ``SUNLS_ATIMES_FAIL_UNREC``,
      indicating an unrecoverable failure in the :math:`J*v` function;
      ``SUNLS_PSOLVE_FAIL_UNREC``, indicating that the preconditioner solve
      function ``psolve`` failed unrecoverably; ``SUNLS_GS_FAIL``, indicating a
      failure in the Gram-Schmidt procedure (generated only in SPGMR or SPFGMR);
      ``SUNLS_QRSOL_FAIL``, indicating that the matrix :math:`R` was found to be
      singular during the QR solve phase (SPGMR and SPFGMR only); or
      ``SUNLS_PACKAGE_FAIL_UNREC``, indicating an unrecoverable failure in an
      external iterative linear solver package.

   .. warning::

      The previous routines :c:func:`KINDlsGetLastFlag` and :c:func:`KINSpilsGetLastFlag`
      are now deprecated.


.. c:function:: char* KINGetLinReturnFlagName(long int lsflag)

   The function :c:func:`KINGetLinReturnFlagName` returns the name of the KINLS
   constant corresponding to ``lsflag``.

   **Arguments:**
      * ``flag`` -- the flag returned by a call to an KINSOL function

   **Return value:**
      * ``char*`` -- the flag name string or if
        :math:`1 \leq \mathtt{lsflag} \leq N` (LU factorization failed), this
        function returns "NONE".

   .. warning::

      The previous routines :c:func:`KINDlsGetReturnFlagName` and
      :c:func:`KINSpilsGetReturnFlagName` are now deprecated.


.. _KINSOL.Usage.CC.user_fct_sim:

User-supplied functions
-----------------------

The user-supplied functions consist of one function defining the nonlinear system,
(optionally) a function that handles error and warning messages,
(optionally) a function that handles informational messages,
(optionally) one or two functions that provides Jacobian-related information for the linear
solver, and (optionally) one or two functions that define the preconditioner for use in
any of the Krylov iterative algorithms.

.. _KINSOL.Usage.CC.user_fct_sim.resFn:

Problem defining function
~~~~~~~~~~~~~~~~~~~~~~~~~

The user must provide a function of type :c:type:`KINSysFn` defined as follows:

.. c:type:: int (*KINSysFn)(N_Vector u, N_Vector fval, void *user_data)

   This function computes the :math:`F(u)` (or :math:`G(u)` for fixed-point iteration and Anderson acceleration)
   for a given value of the vector :math:`u`.

   **Arguments:**
      * ``u`` -- is the current value of the dependent variable vector, :math:`u`
      * ``fval`` -- is the output vector :math:`F(u)`
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        pointer parameter passed to :c:func:`KINSetUserData`

   **Return value:**
      An :c:type:`KINSysFn` function type should return a value of :math:`0` if
      successful, a positive value if a recoverable error occurred (in which
      case KINSOL will attempt to correct), or a negative value if a nonrecoverable
      error occurred. In the last case, the integrator halts. If a recoverable error
      occurred, the integrator will attempt to correct and retry.

   **Notes:**
      Allocation of memory for ``fval`` is handled within KINSOL.


.. _KINSOL.Usage.CC.user_fct_sim.ehFn:

Error message handler function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an alternative to the default behavior of directing error and warning
messages to the file pointed to by ``errfp`` (see :c:func:`KINSetErrFile`), the
user may provide a function of type :c:type:`KINErrHandlerFn` to process any
such messages.  The function type :c:type:`KINErrHandlerFn` is defined as
follows:

.. c:type:: void (*KINErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *user_data)

   This function processes error and warning messages from KINSOL and its
   sub-modules.

   **Arguments:**
      * ``error_code`` -- is the error code
      * ``module`` -- is the name of the KINSOL module reporting the error
      * ``function`` -- is the name of the function in which the error occurred
      * ``eH_data`` -- is a pointer to user data, the same as the ``eh_data``
        parameter passed to :c:func:`KINSetErrHandlerFn`

   **Return value:**
      This function has no return value.

   **Notes:**
      ``error_code`` is negative for errors and positive (``KIN_WARNING``) for
      warnings. If a function that returns a pointer to memory encounters an error,
      it sets ``error_code`` to 0.

.. _KINSOL.Usage.CC.user_fct_sim.ihFn:

Informational message handler function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an alternative to the default behavior of directing informational (meaning non-error) messages
to the file pointed to by ``infofp`` (see :c:func:`KINSetInfoFile`), the user may
provide a function of type :c:type:`KINInfoHandlerFn` to process any such messages.
The function type :c:type:`KINInfoHandlerFn` is defined as follows:

.. c:type:: void (*KINInfoHandlerFn)(const char *module, const char *function, char *msg, void *ih_data)

   This function processes error and warning messages from KINSOL and its
   sub-modules.

   **Arguments:**
      * ``error_code`` -- is the error code
      * ``module`` -- is the name of the KINSOL module reporting the error
      * ``function`` -- is the name of the function in which the error occurred
      * ``ih_data`` -- is a pointer to user data, the same as the ``ih_data``
        parameter passed to :c:func:`KINSetInfoHandlerFn`

   **Return value:**
      This function has no return value.


.. _KINSOL.Usage.CC.user_fct_sim.jacFn:

Jacobian construction (matrix-based linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a matrix-based linear solver module is used (i.e. a non-``NULL``
``SUNMatrix`` object was supplied to :c:func:`KINSetLinearSolver`), the user may
provide a function of type :c:type:`KINLsJacFn` defined as follows:

.. c:type:: int (*KINLsJacFn)(N_Vector u, N_Vector fu, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2)

   This function computes the Jacobian matrix :math:`J(u)` (or an approximation
   to it).

   **Arguments:**
      * ``u`` -- is the current (unscaled) iterate.
      * ``fu`` -- is the current value of the vector, :math:`F(u)`.
      * ``J`` -- is the output (approximate) Jacobian matrix (of type ``SUNMatrix``),
        :math:`F'(u)`.
      * ``user_data`` - is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`KINSetUserData`.
      * ``tmp1``, ``tmp2``, -- are pointers to memory allocated for
        variables of type ``N_Vector`` which can be used by :c:type:`KINLsJacFn` function
        as temporary storage or work space.

   **Return value:**
      An :c:type:`KINLsJacFn` should return :math:`0` if successful, or a
      non-zero value otherwise.

   **Notes:**
      Information regarding the structure of the specific ``SUNMatrix``
      structure (e.g. number of rows, upper/lower bandwidth, sparsity type) may
      be obtained through using the implementation-specific ``SUNMatrix``
      interface functions (see Chapter :numref:`SUNMatrix` for
      details).

      With direct linear solvers (i.e., linear solvers with type
      ``SUNLINEARSOLVER_DIRECT``), the Jacobian matrix :math:`J(u)` is zeroed
      out prior to calling the user-supplied Jacobian function so only nonzero
      elements need to be loaded into ``J``.

      If the user’s :c:type:`KINLsJacFn` function uses difference quotient
      approximations, it may need to access quantities not in the call list.
      These quantities may include the scale vectors and the unit roundoff. To
      obtain the scale vectors, the user will need to add to ``user_data``
      pointers to ``u_scale`` and/or ``f_scale`` as needed. The unit roundoff
      can be accessed as ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

      **dense:**

      A user-supplied dense Jacobian function must load the ``N`` :math:`\times`
      ``N`` dense matrix ``J`` with an approximation to the Jacobian matrix
      :math:`J(u)` at the point (``u``). The accessor macros ``SM_ELEMENT_D``
      and ``SM_COLUMN_D`` allow the user to read and write dense matrix elements
      without making explicit references to the underlying representation of the
      ``SUNMATRIX_DENSE`` type. ``SM_ELEMENT_D(J, i, j)`` references the (``i``,
      ``j``)-th element of the dense matrix ``J`` (with ``i``, ``j``\ :math:`=
      0\ldots \texttt{N}-1`). This macro is meant for small problems for which
      efficiency of access is not a major concern. Thus, in terms of the indices
      :math:`m` and :math:`n` ranging from :math:`1` to :math:`N`, the Jacobian
      element :math:`J_{m,n}` can be set using the statement ``SM_ELEMENT_D(J,
      m-1, n-1) =`` :math:`J_{m,n}`. Alternatively, ``SM_COLUMN_D(J, j)``
      returns a pointer to the first element of the ``j``-th column of ``J``
      (with ``j``\ :math:`= 0\ldots \texttt{N}-1`), and the elements of the
      ``j``-th column can then be accessed using ordinary array indexing.
      Consequently, :math:`J_{m,n}` can be loaded using the statements ``col_n =
      SM_COLUMN_D(J, n-1);`` ``col_n[m-1] =`` :math:`J_{m,n}`. For large
      problems, it is more efficient to use ``SM_COLUMN_D`` than to use
      ``SM_ELEMENT_D``. Note that both of these macros number rows and columns
      starting from :math:`0`. The ``SUNMATRIX_DENSE`` type and accessor macros
      are documented in :numref:`SUNMatrix.Dense`.

      **banded**:

      A user-supplied banded Jacobian function must load the ``N``
      :math:`\times` ``N`` banded matrix ``J`` with an approximation to the
      Jacobian matrix :math:`J(u)` at the point (``u``). The accessor macros
      ``SM_ELEMENT_B``, ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B`` allow the
      user to read and write banded matrix elements without making specific
      references to the underlying representation of the ``SUNMATRIX_BAND``
      type. ``SM_ELEMENT_B(J, i, j)`` references the (``i``, ``j``)-th element
      of the banded matrix ``J``, counting from :math:`0`. This macro is meant
      for use in small problems for which efficiency of access is not a major
      concern. Thus, in terms of the indices :math:`m` and :math:`n` ranging
      from :math:`1` to :math:`\texttt{N}` with :math:`(m,n)` within the band
      defined by ``mupper`` and ``mlower``, the Jacobian element :math:`J_{m,n}`
      can be loaded using the statement ``SM_ELEMENT_B(J, m-1, n-1) =``
      :math:`J_{m,n}`. The elements within the band are those with ``-mupper``
      :math:`\le` ``m-n`` :math:`\le` ``mlower``. Alternatively,
      ``SM_COLUMN_B(J, j)`` returns a pointer to the diagonal element of the
      ``j``-th column of ``J``, and if we assign this address to ``realtype
      *col_j``, then the ``i``-th element of the ``j``-th column is given by
      ``SM_COLUMN_ELEMENT_B(col_j, i, j)``, counting from :math:`0`. Thus, for
      :math:`(m,n)` within the band, :math:`J_{m,n}` can be loaded by setting
      ``col_n = SM_COLUMN_B(J, n-1);`` and ``SM_COLUMN_ELEMENT_B(col_n, m-1,
      n-1) =`` :math:`J_{m,n}`. The elements of the ``j``-th column can also be
      accessed via ordinary array indexing, but this approach requires knowledge
      of the underlying storage for a band matrix of type ``SUNMATRIX_BAND``.
      The array ``col_n`` can be indexed from :math:`-`\ ``mupper`` to
      ``mlower``. For large problems, it is more efficient to use
      ``SM_COLUMN_B`` and ``SM_COLUMN_ELEMENT_B`` than to use the
      ``SM_ELEMENT_B`` macro. As in the dense case, these macros all number rows
      and columns starting from :math:`0`. The ``SUNMATRIX_BAND`` type and
      accessor macros are documented in :numref:`SUNMatrix.Band`.

      **sparse**:

      A user-supplied sparse Jacobian function must load the ``N``
      :math:`\times` ``N`` compressed-sparse-column or compressed-sparse-row
      matrix ``J`` with an approximation to the Jacobian matrix :math:`J(u)` at
      the point (``u``). Storage for ``J`` already exists on entry to this
      function, although the user should ensure that sufficient space is
      allocated in ``J`` to hold the nonzero values to be set; if the existing
      space is insufficient the user may reallocate the data and index arrays as
      needed. The amount of allocated space in a ``SUNMATRIX_SPARSE`` object may
      be accessed using the macro ``SM_NNZ_S`` or the routine
      ``SUNSparseMatrix_NNZ``. The ``SUNMATRIX_SPARSE`` type and accessor macros
      are documented in :numref:`SUNMatrix.Sparse`.

   .. warning::

      The previous function type :c:func:`KINDlsJacFn` is identical to
      :c:type:`KINLsJacFn`, and may still be used for backward-compatibility.
      However, this will be deprecated in future releases, so we recommend that
      users transition to the new function type name soon.


.. _KINSOL.Usage.CC.user_fct_sim.jtimesFn:

Jacobian-vector product (matrix-free linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a matrix-free linear solver is to be used (i.e., a ``NULL``-valued
``SUNMatrix`` was supplied to :c:func:`KINSetLinearSolver`), the user may
provide a function of type :c:type:`KINLsJacTimesVecFn` in the following form,
to compute matrix-vector products :math:`Jv`. If such a function is not
supplied, the default is a difference quotient approximation to these products.

.. c:type:: int (*KINLsJacTimesVecFn)(N_Vector v, N_Vector Jv, N_Vector u, booleantype* new_u, void* user_data)

   This function computes the product :math:`J v` (or an approximation to it).

   **Arguments:**
      * ``v`` -- is the vector by which the Jacobian must be multplied to the right.
      * ``Jv`` -- is the computed output vector.
      * ``u`` -- is the current value of the dependent variable vector.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`KINSetUserData`.

   **Return value:**
      The value returned by the Jacobian-times-vector function should be 0 if
      successful. If a recoverable failure occurred, the return value should be
      positive.  In this case, KINSOL will attempt to correct by calling the
      preconditioner setup function. If this information is current, KINSOL
      halts.  If the Jacobian-times-vector function encounters an unrecoverable
      error, it should return a negative value, prompting KINSOL to halt.

   **Notes:**
      If a user-defined routine is not given, then an internal ``jtimes`` function,
      using a difference quotient approximation, is used.

      This function must return a value of :math:`J*v` that uses the *current* value
      of :math:`J`, i.e. as evaluated at the current :math:`u`.

      If the user’s :c:type:`KINLsJacTimesVecFn` function uses difference quotient
      approximations, it may need to access quantities not in the call list. These
      might include the scale vectors and the unit roundoff. To obtain the scale
      vectors, the user will need to add to ``user_data`` pointers to ``u_scale``
      and/or ``f_scale`` as needed. The unit roundoff can be accessed as
      ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

   .. warning::

      The previous function type :c:type:`KINSpilsJacTimesVecFn` is identical to
      :c:type:`KINLsJacTimesVecFn`, and may still be used for
      backward-compatibility. However, this will be removed in future
      releases, so we recommend that users transition to the new function type
      name soon.


.. _KINSOL.Usage.CC.user_fct_sim.psolveFn:

Preconditioner solve (iterative linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a user-supplied preconditioner is to be used with a ``SUNLinearSolver``
solver module, then the user must provide a function to solve the linear system
:math:`Pz = r` where :math:`P` is the preconditioner matrix which approximates
(at least crudely) the Jacobian matrix :math:`J = F'(u)`. This function must be
of type :c:type:`KINLsPrecSolveFn`, defined as follows:

.. c:type:: int (*KINLsPrecSolveFn)(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector v, void *user_data)

   This function solves the preconditioning system :math:`Pz = r`.

   **Arguments:**
      * ``u`` -- is the current (unscaled) value of the iterate.
      * ``uscale`` -- is a vector containing diagonal elements of the scaling matrix ``u``
      * ``fval`` -- is the vector :math:`F(u)` evaluated at ``u``
      * ``fscale`` -- is a vector containing diagonal elements of the scaling matrix for ``fval``
      * ``v`` -- on inpuut, ``v`` is set to the right-hand side vector of the linear system, ``r``. On    output, ``v`` must contain the solution ``z`` of the linear system :math:`Pz=r`
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`KINSetUserData`.

   **Return value:**
      The value returned by the preconditioner solve function should be 0 if
      successful, positive for a recoverable error, or negative for an
      unrecoverable error.

   **Notes:**
      If the preconditioner solve function fails recoverably and if the
      preconditioner information (set by the preconditioner setup function) is
      out of date, KINSOL attempts to correct by calling the setup function. If
      the preconditioner data is current, KINSOL halts.


.. _KINSOL.Usage.CC.user_fct_sim.precondFn:

Preconditioner setup (iterative linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the user’s preconditioner requires that any Jacobian-related data be
evaluated or preprocessed, then this needs to be done in a user-supplied
function of type :c:type:`KINLsPrecSetupFn`, defined as follows:

.. c:type:: int (*KINLsPrecSetupFn)(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, void *user_data)

   This function evaluates and/or preprocesses Jacobian-related data needed
   by the preconditioner solve function.

   **Arguments:**
      * ``u`` -- is the current (unscaled) value of the iterate.
      * ``uscale`` -- is a vector containing diagonal elements of the scaling matrix ``u``
      * ``fval`` -- is the vector :math:`F(u)` evaluated at ``u``
      * ``fscale`` -- is a vector containing diagonal elements of the scaling matrix for ``fval``
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`KINSetUserData`.

   **Return value:**
      The value returned by the preconditioner setup function should be 0 if
      successful, positive for a recoverable error (in which case the step will be
      retried), or negative for an unrecoverable error (in which case the
      integration is halted).

   **Notes:**
      The user-supplied preconditioner setup subroutine should compute the right
      preconditioner matrix :math:`P` (stored in the memory block referenced by the
      ``user_data`` pointer) used to form the scaled preconditioned linear system

      .. math:: (D_F J(u) P^{-1} D_u^{-1}) (D_u P x) = - D_F F(u) \, ,

      where :math:`D_u` and :math:`D_F` denote the diagonal scaling matrices whose
      diagonal elements are stored in the vectors ``uscale`` and ``fscale``,
      respectively.

      The preconditioner setup routine will not be called prior to every call made
      to the preconditioner solve function, but will instead be called only as
      often as necessary to achieve convergence of the Newton iteration.

      If the user’s :c:type:`KINLsPrecSetupFn` function uses difference quotient
      approximations, it may need to access quantities not in the call list. These
      might include the scale vectors and the unit roundoff. To obtain the scale
      vectors, the user will need to add to ``user_data`` pointers to ``u_scale``
      and/or ``f_scale`` as needed. The unit roundoff can be accessed as
      ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

      If the preconditioner solve routine requires no preparation, then a
      preconditioner setup function need not be given.


.. _KINSOL.Usage.CC.kin_bbdpre:

A parallel band-block-diagonal preconditioner module
----------------------------------------------------

The efficiency of Krylov iterative methods for the solution of linear systems
can be greatly enhanced through preconditioning. For problems in which the user
cannot define a more effective, problem-specific preconditioner, KINSOL provides
a band-block-diagonal preconditioner module KINBBDPRE, to be used with the
parallel ``N_Vector`` module described in :numref:`NVectors.NVParallel`.

This module provides a preconditioner matrix for KINSOL that is block-diagonal
with banded blocks. The blocking corresponds to the distribution of the
dependent variable vector :math:`u` amongst the processes. Each preconditioner
block is generated from the Jacobian of the local part (associated with the
current process) of a given function :math:`G(u)` approximating :math:`F(u)`
(:math:`G = F` is allowed). The blocks are generated by each process via a
difference quotient scheme, utilizing a specified band structure. This structure
is given by upper and lower half-bandwidths, ``mudq`` and ``mldq``, defined as
the number of non-zero diagonals above and below the main diagonal,
respectively. However, from the resulting approximate Jacobain blocks, only a
matrix of bandwidth ``mukeep`` :math:`+` ``mlkeep`` :math:`+ 1` is retained.

Neither pair of parameters need be the true half-bandwidths of the Jacobian of
the local block of :math:`G`, if smaller values provide a more efficient
preconditioner. Such an efficiency gain may occur if the couplings in the system
outside a certain bandwidth are considerably weaker than those within the band.
Reducing ``mukeep`` and ``mlkeep`` while keeping ``mudq`` and ``mldq`` at their
true values, discards the elements outside the narrower band. Reducing both
pairs has the additional effect of lumping the outer Jacobian elements into the
computed elements within the band, and requires more caution and experimentation
to see whether the lower cost of narrower band matrices offsets the loss of
accuracy in the blocks.

The KINBBDPRE module calls two user-provided functions to construct :math:`P`: a
required function ``Gloc`` (of type :c:type:`KINBBDLocalFn`) which approximates the
nonlinear system function :math:`G(u) \approx F(u)` and which is computed
locally, and an optional function ``Gcomm`` (of type :c:type:`KINBBDCommFn`) which
performs all interprocess communication necessary to evaluate the approximate
function :math:`G`. These are in addition to the user-supplied nonlinear system
function that evaluates :math:`F(u)`. Both functions take as input the same
pointer ``user_data`` as that passed by the user to :c:func:`KINSetUserData` and
passed to the user’s function ``func``, and neither function has a return value.
The user is responsible for providing space (presumably within ``user_data``)
for components of ``u`` that are communicated by ``Gcomm`` from the other
processes, and that are then used by ``Gloc``, which should not do any
communication.

.. c:type:: int (*KINBBDLocalFn)(sunindextype Nlocal, N_Vector u, N_Vector gval, void *user_data)

   This ``Gloc`` function computes :math:`G(u)`, and outputs the resulting vector
   as ``gval``.

   **Arguments:**
      * ``Nlocal`` -- is the local vector length.
      * ``u`` -- is the current value of the iterate.
      * ``gval`` -- is the output vector.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`KINSetUserData`.

   **Return value:**
      An :c:type:`KINBBDLocalFn` function type should return 0 to indicate success,
      or non-zero if an error occured.

   **Notes:**
      This function must assume that all inter-processor communication of data
      needed to calculate ``gval`` has already been done, and this data is
      accessible within ``user_data``.

   The case where :math:`G` is mathematically identical to :math:`F` is allowed.


.. c:type:: int (*KINBBDCommFn)(sunindextype Nlocal, N_Vector u, void *user_data)

   This ``Gcomm`` function performs all inter-processor communications necessary
   for the execution of the ``Gloc`` function above, using the input vectors
   ``u``.

   **Arguments:**
      * ``Nlocal`` -- is the local vector length.
      * ``u`` -- is the current value of the iterate.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`KINSetUserData`.

   **Return value:**
      An :c:type:`KINBBDLocalFn` function type should return 0 to indicate success,
      or non-zero if an error occured.

   **Notes:**
      The ``Gcomm`` function is expected to save communicated data in space defined
      within the structure ``user_data``.

      Each call to the ``Gcomm`` function is preceded by a call to the residual
      function ``func`` with the same ``u`` argument. Thus
      ``Gcomm`` can omit any communications done by ``func`` if relevant to the
      evaluation of ``Gloc``. If all necessary communication was done in ``func``,
      then ``Gcomm = NULL`` can be passed in the call to :c:func:`KINBBDPrecInit`.


Besides the header files required for the integration of the DAE problem (see
:numref:`KINSOL.Usage.CC.header_sim`), to use the KINBBDPRE module, the main program
must include the header file ``kin_bbdpre.h`` which declares the needed function
prototypes.

The following is a summary of the usage of this module and describes the
sequence of calls in the user main program.  Steps that are unchanged from the
user main program presented in :numref:`KINSOL.Usage.CC.skeleton_sim` are not bold.

#. Initialize parallel or multi-threaded environment (*if appropriate*)

#. Create the SUNDIALS context object

#. Set the problem dimensions etc.

#. Create the vector with the initial guess

#. Create matrix object (*if appropriate*)

#. **Create linear solver object** (*if appropriate*)

   When creating the iterative linear solver object, specify the use of right
   preconditioning (``SUN_PREC_RIGHT``) as KINSOL only supports right preconditioning.

#. Create nonlinear solver object (*if appropriate*)

#. Create KINSOL object

#. Initialize KINSOL solver

#. Attach the linear solver (*if appropriate*)

#. **Set linear solver optional inputs** (*if appropriate*)

   Note that the user should not overwrite the preconditioner setup function or
   solve function through calls to :c:func:`KINSetPreconditioner` optional input
   function.

#. **Initialize the KINBBDPRE preconditioner module**

   Call :c:func:`KINBBDPrecInit` to allocate memory and initialize the internal
   preconditioner data. The last two arguments of :c:func:`KINBBDPrecInit` are
   the two user-supplied functions described above.

#. Set optional inputs

#. Solve problem

#. **Get optional outputs**

   Additional optional outputs associated with KINBBDPRE are available by way of
   two routines described below, :c:func:`KINBBDPrecGetWorkSpace` and
   :c:func:`KINBBDPrecGetNumGfnEvals`.

#. Deallocate memory

#. Finalize MPI, if used


The user-callable functions that initialize or re-initialize the KINBBDPRE
preconditioner module are described next.

.. c:function:: int KINBBDPrecInit(void* kin_mem, sunindextype Nlocal, sunindextype mudq, sunindexype mldq, sunindextype mukeep, sunindextype mlkeep, realtype dq_rel_u, KINBBDLocalFn Gloc, KINBBDCommFn Gcomm)

   The function :c:func:`KINBBDPrecInit` initializes and allocates  memory for
   the KINBBDPRE preconditioner.

   **Arguments:**
     * ``kin_mem`` -- pointer to the KINSOL memory block.
     * ``Nlocal`` -- local vector length.
     * ``mudq`` -- upper half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``mldq`` -- lower half-bandwidth to be used in the difference-quotient Jacobian approximation.
     * ``mukeep`` -- upper half-bandwidth of the retained banded approximate Jacobian block.
     * ``mlkeep`` -- lower half-bandwidth of the retained banded approximate Jacobian block.
     * ``dq_rel_u`` -- the relative increment in components of ``u`` used in the difference quotient approximations.
       The default is :math:`\texttt{dq\_rel\_u} = \sqrt{\text{unit roundoff}}` , which can be specified by passing ``dq_rel_u= 0.0``.
     * ``Gloc`` -- the CC function which computes the approximation :math:`G(u) \approx F(u)`.
     * ``Gcomm`` -- the optional CC function which performs all interprocess communication required for the computation of :math:`G(u)`.

   **Return value:**
     * ``KINLS_SUCCESS`` -- The call to :c:func:`KINBBDPrecInit` was successful.
     * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer was ``NULL``.
     * ``KINLS_MEM_FAIL`` -- A memory allocation request has failed.
     * ``KINLS_LMEM_NULL`` -- The KINLS linear solver interface has not been initialized.
     * ``KINLS_ILL_INPUT`` -- The supplied vector implementation was not compatible with the block band preconditioner.

   **Notes:**
     If one of the half-bandwidths ``mudq`` or ``mldq`` to be used in the
     difference-quotient calculation of the approximate Jacobian is negative
     or exceeds the value ``Nlocal-1``, it is replaced with 0 or ``Nlocal-1`` accordingly.

     The half-bandwidths ``mudq`` and ``mldq`` need
     not be the true half-bandwidths of the Jacobian of the local block of :math:`G`,
     when smaller values may provide greater efficiency.

     Also, the half-bandwidths ``mukeep`` and ``mlkeep`` of the retained
     banded approximate Jacobian block may be even smaller, to reduce
     storage and computation costs further.

     For all four half-bandwidths, the values need not be the same for
     every process.



The following two optional output functions are available for use with the
KINBBDPRE module:

.. c:function:: int KINBBDPrecGetWorkSpace(void * kin_mem, long int * lenrwBBDP, long int * leniwBBDP)

   The function :c:func:`KINBBDPrecGetWorkSpace` returns the local sizes of the
   KINBBDPRE real and integer workspaces.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``lenrwBBDP`` -- local number of real values in the KINBBDPRE workspace.
      * ``leniwBBDP`` -- local number of integer values in the KINBBDPRE workspace.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer was ``NULL``.
      * ``KINLS_PMEM_NULL`` -- The KINBBDPRE preconditioner has not been
        initialized.

   **Notes:**
      The workspace requirements reported by this routine correspond only to memory
      allocated within the KINBBDPRE module (the banded matrix approximation,
      banded ``SUNLinearSolver`` object, temporary vectors).  These values
      are local to each process.

      The workspaces referred to here exist in addition
      to those given by the corresponding :c:func:`KINGetLinWorkSpace` function.

.. c:function:: int KINBBDPrecGetNumGfnEvals(void * kin_mem, long int * ngevalsBBDP)

   The function :c:func:`KINBBDPrecGetNumGfnEvals` returns the cumulative number
   of calls to the user ``Gres`` function due to the finite difference
   approximation of the Jacobian blocks used within KINBBDPRE's preconditioner
   setup function.

   **Arguments:**
      * ``kin_mem`` -- pointer to the KINSOL solver object.
      * ``ngevalsBBDP`` -- the cumulative number of calls to the user ``Gres``
        function.

   **Return value:**
      * ``KINLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``KINLS_MEM_NULL`` -- The ``kin_mem`` pointer was ``NULL``.
      * ``KINLS_PMEM_NULL`` -- The KINBBDPRE preconditioner has not been
        initialized.


In addition to the ``ngevalsBBDP`` evaluations of ``Gres``, the costs associated
with KINBBDPRE also includes ``nlinsetups`` LU factorizations, ``nlinsetups``
calls to ``Gcomm``, ``npsolves`` banded backsolve calls, and ``nrevalsLS``
residual function evaluations, where ``nlinsetups`` is an optional KINSOL output
(see :numref:`KINSOL.Usage.CC.optional_output.optout_main`), and ``npsolves`` and
``nrevalsLS`` are linear solver optional outputs (see
:numref:`KINSOL.Usage.CC.optional_output.optout_ls`).

.. _KINSOL.Usage.CC.kinalternative:

Alternative to KINSOL for difficult systems
-------------------------------------------

A nonlinear system :math:`F(u) = 0` may be difficult to solve with KINSOL (or
any other nonlinear system solver) for a variety of reasons. The possible
reasons include high nonlinearity, small region of convergence, and lack of a
good initial guess. For systems with such difficulties, there is an alternative
approach that may be more successful. This is an old idea, but deserves some new
attention.

If the nonlinear system is :math:`F(u) = 0`, consider instead the ODE system
:math:`du/dt = - M^{-1} F(u)`, where :math:`M` is a nonsingular matrix that is
an approximation (even a crude approximation) to the system Jacobian :math:`F_u
= dF/du`. Whatever :math:`M` is, if this ODE is solved until it reaches a steady
state :math:`u^*`, then :math:`u^*` is a zero of the right-hand side of the ODE,
and hence a solution to :math:`F(u) = 0`. There is no issue of having a close
enough initial guess.

A further basis for this choice of ODE is the following: If :math:`M`
approximates :math:`F_u`, then the Jacobian of the ODE system,
:math:`-M^{-1}F_u`, is approximately equal to :math:`-I` where :math:`I` is the
identity matrix. This means that (in a local approximation sense) the solution
modes of the ODE behave like :math:`\exp(-t)`, and that asymptotically the
approach to the steady state goes as :math:`\exp(-t)`. Of course, the closer
:math:`M` is to :math:`F_u`, the better this basis applies.

Using (say) CVODE to solve the above ODE system requires, in addition to the
objective function :math:`F(u)`, the calculation of a suitable matrix :math:`M`
and its inverse, or at least a routine that solves linear systems :math:`Mx =
b`. This is similar to the KINSOL requirement of supplying the system Jacobian
:math:`J` (or solutions to :math:`Jx = b`), but differs in that :math:`M` may be
simpler than :math:`J` and hence easier to deal with. Depending on the nature of
:math:`M`, this may be handled best with a direct solver, or with a
preconditioned Krylov solver. The latter calls for the use of a preconditioner
:math:`P` that may be a crude approximation to :math:`M`, hence even easier to
solve. Note if using ARKODE, the ODE system may be posed in the linearly
implicit from :math:`M du/dt = -F(u)` where :math:`M` is the "mass matrix" for
the system. This use case requires supplying ARKODE with a function to evaluate
:math:`M` or to compute its action on a vector (:math:`Mv = w`) and attaching a
linear solver (direct or iterative) to solve the linear systems :math:`Mx = b`.

The solution of the ODE may be made easier by solving instead the equivalent
DAE, :math:`M du/dt + F(u) = 0`. Applying IDA to this system requires solutions
to linear systems whose matrix is the DAE system Jacobian, :math:`J = F_u +
\alpha M`, where :math:`\alpha` is the scalar coefficient :math:`c_j` supplied
to the user’s Jacobian or preconditioner routine. Selecting a preconditioned
Krylov method requires an approximation to this Jacobian as preconditioner
:math:`P`. Given that :math:`M` approximates :math:`F_u` (possibly crudely), the
appropriate approximation to :math:`J` is :math:`P = M + \alpha M = (1 +
\alpha)M`. Again the user must supply a routine that solves linear systems
:math:`Px = b`, or :math:`M x = b/(1 + \alpha)`. If M is too difficult to solve,
than an approximation :math:`M'` that is easier can be substituted, as long as
it achieves convergence. As always, there is a trade-off between the expense of
solving :math:`M` and the difficulty of achieving convergence in the linear
solver.

For the solution of either the ODE or DAE system above, the chances for
convergence can be improved with a piecewise constant choice for :math:`M`.
Specifically, starting from an initial guess :math:`u_0`, an initial choice for
:math:`M` might be :math:`M_0 = F_u(u_0)`, or some approximation to this
Jacobian. Then one could integrate :math:`M_0 du/dt + F(u) = 0` from :math:`t =
0` to :math:`t = T` for some sizable value :math:`T`, evaluate
:math:`F_u(u(T))`, and take :math:`M_1` to be an approximation to that Jacobian.
Then integrate using :math:`M_1` from :math:`t = T` to :math:`t = 2T`, and
repeat the process until it converges to a steady state.
