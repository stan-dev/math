.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _CVODES.Usage.SIM:

Using CVODES for IVP Solution
=============================

This chapter is concerned with the use of CVODES for the solution of initial
value problems (IVPs). The following sections treat the header files and the
layout of the user’s main program, and provide descriptions of the CVODES
user-callable functions and user-supplied functions.

The sample programs described in the companion document :cite:p:`cvodes_ex` may
also be helpful. Those codes may be used as templates (with the removal of some
lines used in testing) and are included in the CVODES package.

Users with applications written in Fortran should see
:numref:`SUNDIALS.Fortran`, which describes interfacing with CVODES from
Fortran.

The user should be aware that not all ``SUNLinearSolver`` and ``SUNMatrix``
modules are compatible with all ``N_Vector`` implementations. Details on
compatibility are given in the documentation for each ``SUNMatrix`` module
(:numref:`SUNMatrix`) and each ``SUNLinearSolver`` module (:numref:`SUNLinSol`).
For example, ``NVECTOR_PARALLEL`` is not compatible with the dense, banded, or
sparse ``SUNMatrix`` types, or with the corresponding dense, banded, or sparse
``SUNLinearSolver`` modules. Please check :numref:`SUNMatrix` and
:numref:`SUNLinSol` to verify compatibility between these modules. In addition
to that documentation, we note that the CVBANDPRE preconditioning module is only
compatible with the ``NVECTOR_SERIAL``, ``NVECTOR_OPENMP``, and
``NVECTOR_PTHREADS`` vector implementations, and the preconditioner module
CVBBDPRE can only be used with ``NVECTOR_PARALLEL``. It is not recommended to
use a threaded vector module with SuperLU_MT unless it is the ``NVECTOR_OPENMP``
module, and SuperLU_MT is also compiled with OpenMP.

CVODES uses various constants for both input and output. These are defined as
needed in this chapter, but for convenience are also listed separately in
:numref:`CVODES.Constants`.

.. _CVODES.Usage.SIM.file_access:

Access to library and header files
----------------------------------

At this point, it is assumed that the installation of CVODES,
following the procedure described in :numref:`Installation`, has been completed successfully.

Regardless of where the user’s application program resides, its
associated compilation and load commands must make reference to the
appropriate locations for the library and header files required by
CVODES. The relevant library files are

.. code-block::

  <libdir>/libsundials_cvodes.<so|a>
  <libdir>/libsundials_nvec*.<so|a>
  <libdir>/libsundials_sunmat*.<so|a>
  <libdir>/libsundials_sunlinsol*.<so|a>
  <libdir>/libsundials_sunnonlinsol*.<so|a>

where the file extension ``.so`` is typically for shared libraries and
``.a`` for static libraries. The relevant header files are located in the
subdirectories

.. code-block::

  <incdir>/cvodes
  <incdir>/sundials
  <incdir>/nvector
  <incdir>/sunmatrix
  <incdir>/sunlinsol
  <incdir>/sunnonlinsol

The directories ``libdir`` and ``incdir`` are the install library and
include directories, respectively. For a default installation, these are
``<instdir>/lib`` and ``<instdir>/include``, respectively, where ``instdir`` is
the directory where SUNDIALS was installed (:numref:`Installation`).

.. _CVODES.Usage.SIM.data_types:

.. include:: ../../../../shared/Types.rst

.. _CVODES.Usage.SIM.header_sim:

Header files
------------

The calling program must include several header files so that various
macros and data types can be used. The header file that is always
required is:

* ``cvodes/cvodes.h`` the main header file for CVODES, which defines the several types and various constants, and includes function prototypes. This includes the header file for CVLS, ``cvodes/cvodes_ls.h``.

Note that ``cvodes.h`` includes ``sundials_types.h``, which defines the types, ``realtype``, ``sunindextype``, and ``booleantype`` and the constants ``SUNFALSE`` and ``SUNTRUE``.

The calling program must also include an ``N_Vector`` implementation header file, of the form ``nvector/nvector_*.h``. See :numref:`NVectors` for the appropriate name. This file in turn includes the header file ``sundials_nvector.h`` which defines the abstract data type.

If using a non-default nonlinear solver module, or when interacting with a ``SUNNonlinearSolver`` module directly, the calling program must also include a ``SUNNonlinearSolver`` implementation header file, of the form ``sunnonlinsol/sunnonlinsol_*.h`` where is the name of the nonlinear solver module (see :numref:`SUNNonlinSol` for more information).
This file in turn includes the header file which defines the abstract data type.

If using a nonlinear solver that requires the solution of a linear system of the form :eq:`CVODES_Newton` (e.g., the default Newton iteration), then a linear solver module header file will be required.

Other headers may be needed, according to the choice of preconditioner, etc. For example, in the example (see :cite:p:`cvodes_ex`), preconditioning is done with a block-diagonal matrix.
For this, even though the ``SUNLINSOL_SPGMR`` linear solver is used, the header is included for access to the underlying generic dense matrix arithmetic routines.

.. _CVODES.Usage.SIM.skeleton_sim:

A skeleton of the user’s main program
-------------------------------------

The following is a skeleton of the user’s main program (or calling
program) for the integration of an ODE IVP. Most of the steps are
independent of the ``N_Vector``, ``SUNMatrix``, ``SUNLinearSolver``, and
``SUNNonlinearSolver`` implementations used. For the steps that are not, refer
to :numref:`NVectors`, :numref:`SUNMatrix`, :numref:`SUNLinSol`, and
:numref:`SUNNonlinSol` for the specific name of the
function to be called or macro to be referenced.

   #. **Initialize parallel or multi-threaded environment, if appropriate**
      For example, call ``MPI_Init`` to initialize MPI if used, or set the number
      of threads to use within the threaded vector functions if used.

   #. **Create the SUNDIALS context object**
      Call :c:func:`SUNContext_Create` to allocate the ``SUNContext`` object.

   #. **Set problem dimensions etc.**
      This generally includes the problem size ``N``, and may include the local
      vector length ``Nlocal``.

      Note: The variables ``N`` and ``Nlocal`` should be of type ``sunindextype``.

   #. **Set vector of initial values**
      To set the vector of initial values, use the appropriate functions
      defined by the particular ``N_Vector`` implementation.

      For native SUNDIALS vector implementations, use a call of the form ``y0 = N_VMake_***(..., ydata)`` if the array containing the initial values of :math:`y` already exists. Otherwise, create a new vector by making a call of the form ``N_VNew_***(...)``, and then set its elements by accessing the underlying data with a call of the form ``ydata = N_VGetArrayPointer(y0)``.

      For HYPRE and PETSC vector wrappers, first create and
      initialize the underlying vector, and then create an ``N_Vector``
      wrapper with a call of the form ``y0 = N_VMake_***(yvec)``, where ``yvec`` is a HYPRE or PETSC
      vector. Note that calls like ``N_VNew_***(...)`` and ``N_VGetArrayPointer(...)`` are not available
      for these vector wrappers.

      See :numref:`NVectors` for details.

   #. **Create CVODES object**
      Call :c:func:`CVodeCreate` to create the CVODES memory block and to specify the linear
      multistep method. :c:func:`CVodeCreate` returns a pointer to the CVODES memory
      structure.

      See :numref:`CVODES.Usage.SIM.user_callable.cvodemalloc` for details.

   #. **Initialize CVODES solver**
      Call :c:func:`CVodeInit` to provide required problem specifications, allocate internal
      memory for CVODES, and initialize CVODES. :c:func:`CVodeInit` returns a flag, the
      value of which indicates either success or an illegal argument value.

      See :numref:`CVODES.Usage.SIM.user_callable.cvodemalloc` for details.

   #. **Specify integration tolerances**
      Call :c:func:`CVodeSStolerances` or :c:func:`CVodeSVtolerances` to specify either a scalar relative tolerance and scalar absolute tolerance, or a scalar relative tolerance and a vector of absolute tolerances, respectively. Alternatively, call :c:func:`CVodeWFtolerances` to specify a function which sets directly the weights used in evaluating WRMS vector norms.

      See :numref:`CVODES.Usage.SIM.user_callable.cvtolerances` for details.

   #. **Create matrix object**
      If a nonlinear solver requiring a linear solve will be used (e.g.,
      the default Newton iteration) and the linear solver will be a
      matrix-based linear solver, then a template Jacobian matrix must be
      created by calling the appropriate constructor function defined by
      the particular ``SUNMatrix`` implementation.

      For the native SUNDIALS ``SUNMatrix`` implementations, the matrix object may
      be created using a call of the form ``SUN***Matrix(...)`` where ``***`` is
      the name of the matrix (see :numref:`SUNMatrix` for details).

   #. **Create linear solver object**
      If a nonlinear solver requiring a linear solver is chosen (e.g., the
      default Newton iteration), then the desired linear solver object must
      be created by calling the appropriate constructor function defined by
      the particular ``SUNLinearSolver`` implementation.

      For any of the SUNDIALS-supplied ``SUNLinearSolver`` implementations,
      the linear solver object may be created using a call of the form
      ``SUNLinearSolver LS = SUNLinSol_*(...);``
      where ``*`` can be replaced with “Dense”, “SPGMR”, or other options, as
      discussed in :numref:`CVODES.Usage.SIM.user_callable.lin_solv_init` and
      :numref:`SUNLinSol`.

   #. **Set linear solver optional inputs**
      Call functions from the selected linear solver module to change
      optional inputs specific to that linear solver. See the documentation
      for each ``SUNLinearSolver`` module in
      :numref:`SUNLinSol` for details.

   #. **Attach linear solver module**
      If a nonlinear solver requiring a linear solver is chosen (e.g., the
      default Newton iteration), then initialize the CVLS linear solver
      interface by attaching the linear solver object (and matrix object,
      if applicable) with a call ``ier = CVodeSetLinearSolver(cvode_mem, NLS)`` (for details see
      :numref:`CVODES.Usage.SIM.user_callable.lin_solv_init`):

      Alternately, if the CVODES-specific diagonal linear solver module,
      CVDIAG, is desired, initialize the linear solver module and attach
      it to CVODES with the call to :c:func:`CVodeSetLinearSolver`.

   #. **Set optional inputs**
      Call ``CVodeSet***`` functions to change any optional inputs that control the
      behavior of CVODES from their default values. See
      :numref:`CVODES.Usage.SIM.optional_input` for details.

   #. **Create nonlinear solver object** (*optional*)
      If using a non-default nonlinear solver (see :numref:`CVODES.Usage.SIM.nonlin_solv_init`), then create the desired
      nonlinear solver object by calling the appropriate constructor
      function defined by the particular ``SUNNonlinearSolver`` implementation
      (e.g., ``NLS = SUNNonlinSol_***(...);`` where ``***`` is the name of the nonlinear solver (see
      :numref:`SUNNonlinSol` for details).

   #. **Attach nonlinear solver module** (*optional*)
      If using a non-default nonlinear solver, then initialize the
      nonlinear solver interface by attaching the nonlinear solver object
      by calling ``ier = CVodeSetNonlinearSolver`` (see :numref:`CVODES.Usage.SIM.nonlin_solv_init` for details).

   #. **Set nonlinear solver optional inputs** (*optional*)
      Call the appropriate set functions for the selected nonlinear solver
      module to change optional inputs specific to that nonlinear solver.
      These *must* be called after :c:func:`CVodeInit` if using the default nonlinear solver or
      after attaching a new nonlinear solver to CVODES, otherwise the
      optional inputs will be overridden by CVODES defaults. See
      :numref:`SUNNonlinSol` for more information on
      optional inputs.

   #. **Specify rootfinding problem** (*optional*)
      Call :c:func:`CVodeRootInit` to initialize a rootfinding problem to be solved during the
      integration of the ODE system. See :numref:`CVODES.Usage.SIM.cvrootinit`, and
      see :numref:`CVODES.Usage.SIM.optional_input.optin_root` for relevant optional input
      calls.

   #. **Advance solution in time**
      For each point at which output is desired, call ``ier = CVode(cvode_mem,
      tout, yout, tret itask)``. Here ``itask`` specifies the return mode. The
      vector ``yout`` (which can be the same as the vector ``y0`` above) will
      contain :math:`y(t)`. See :c:func:`CVode` for details.

   #. **Get optional outputs**
      Call ``CV*Get*`` functions to obtain optional output. See :numref:`CVODES.Usage.SIM.optional_output` for details.

   #. **Deallocate memory for solution vector**
      Upon completion of the integration, deallocate memory for the vector ``y``
      (or ``yout``) by calling the appropriate destructor function defined by
      the ``N_Vector`` implementation.

   #. **Free solver memory**
      Call :c:func:`CVodeFree` to free the memory allocated by CVODES.

   #. **Free nonlinear solver memory** (*optional*)
      If a non-default nonlinear solver was used, then call :c:func:`SUNNonlinSolFree` to free any
      memory allocated for the ``SUNNonlinearSolver`` object.

   #. **Free linear solver and matrix memory**
      Call :c:func:`SUNLinSolFree` and :c:func:`SUNMatDestroy` to free any memory allocated for the linear solver and
      matrix objects created above.

   #. **Free the SUNContext object**
      Call :c:func:`SUNContext_Free` to free the memory allocated for the ``SUNContext`` object.

   #. **Finalize MPI, if used**
      Call ``MPI_Finalize`` to terminate MPI.


.. _CVODES.Usage.SIM.user_callable:

User-callable functions
-----------------------

This section describes the CVODES functions that are called by the
user to setup and then solve an IVP. Some of these are required.
However, starting with :numref:`CVODES.Usage.SIM.optional_input`, the functions
listed involve optional inputs/outputs or restarting, and those
paragraphs may be skipped for a casual use of CVODES. In any case,
refer to :numref:`CVODES.Usage.SIM.skeleton_sim` for the correct order of these
calls.

On an error, each user-callable function returns a negative value and
sends an error message to the error handler routine, which prints the
message on ``stderr`` by default. However, the user can set a file as error output
or can provide his own error handler function (see :numref:`CVODES.Usage.SIM.optional_input.optin_main`).

.. _CVODES.Usage.SIM.user_callable.cvodemalloc:

CVODES initialization and deallocation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following three functions must be called in the order listed. The
last one is to be called only after the IVP solution is complete, as it
frees the CVODES memory block created and allocated by the first two
calls.

.. c:function:: void* CVodeCreate(int lmm, SUNContext sunctx)

   The function :c:func:`CVodeCreate` instantiates a CVODES solver object and
   specifies the solution method.

   **Arguments:**
      - ``lmm`` -- specifies the linear multistep method and must be one of two possible values: ``CV_ADAMS`` or ``CV_BDF``.
      - ``sunctx`` -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return Value:**
      - If successful, :c:func:`CVodeCreate` returns a pointer to the newly created CVODES memory block.  Otherwise, it returns ``NULL``.

   **Notes:**
      The recommended choices for ``lmm`` are ``CV_ADAMS`` for nonstiff problems
      and ``CV_BDF`` for stiff problems. The default Newton iteration is
      recommended for stiff problems, and the fixed-point solver (previously referred
      to as the functional iteration in this guide) is
      recommended for nonstiff problems. For details on how to attach a
      different nonlinear solver module to CVODES see the description of
      :c:func:`CVodeSetNonlinearSolver`.


.. c:function:: int CVodeInit(void* cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)

   The function ``CVodeInit`` provides required problem and solution
   specifications, allocates internal memory, and initializes CVODES.

   **Arguments:**
      - ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
      - ``f`` -- is the C function which computes the right-hand side function f in the ODE. This function has the form ``f(t, y, ydot, user_data)`` (for full details see :numref:`CVODES.Usage.SIM.user_supplied.rhsFn`).
      - ``t0`` -- is the initial value of t.
      - ``y0`` -- is the initial value of y.

   **Return Value:**
      - ``CV_SUCCESS`` -- The call was successful.
      - ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
      - ``CV_MEM_FAIL`` -- A memory allocation request has failed.
      - ``CV_ILL_INPUT`` -- An input argument to ``CVodeInit`` has an illegal value.

   **Notes:**
      If an error occurred, ``CVodeInit`` also sends an error message to the
      error handler function.

.. c:function:: void CVodeFree(void** cvode_mem);

   The function ``CVodeFree`` frees the memory allocated by
   a previous call to :c:func:`CVodeCreate`.

   **Arguments:**
     - Pointer to the CVODES memory block.

   **Return Value:**
      - The function ``CVodeFree`` has no return value.



.. _CVODES.Usage.SIM.user_callable.cvtolerances:

CVODES tolerance specification functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the following three functions must be called to specify the
integration tolerances (or directly specify the weights used in
evaluating WRMS vector norms). Note that this call must be made after
the call to :c:func:`CVodeInit`.

.. c:function:: int CVodeSStolerances(void* cvode_mem, realtype reltol, realtype abstol)

   The function ``CVodeSStolerances`` specifies scalar relative and absolute  tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``reltol`` -- is the scalar relative error tolerance.
     * ``abstol`` -- is the scalar absolute error tolerance.

   **Return value:**
     * ``CV_SUCCESS`` -- The call was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized.
     * ``CV_NO_MALLOC`` -- The allocation function returned ``NULL``.
     * ``CV_ILL_INPUT`` -- One of the input tolerances was negative.

.. c:function:: int CVodeSVtolerances(void* cvode_mem, realtype reltol, N_Vector abstol)

   The function ``CVodeSVtolerances`` specifies scalar relative tolerance and  vector absolute tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``reltol`` -- is the scalar relative error tolerance.
     * ``abstol`` -- is the vector of absolute error tolerances.

   **Return value:**
     * ``CV_SUCCESS`` -- The call was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized.
     * ``CV_NO_MALLOC`` -- The allocation function returned ``NULL``.
     * ``CV_ILL_INPUT`` -- The relative error tolerance was negative or the absolute tolerance had a negative component.

   **Notes:**
      This choice of tolerances is important when the absolute error tolerance needs to  be different for each component of the state vector y.

.. c:function:: int CVodeWFtolerances(void* cvode_mem, CVEwtFn efun)

   The function ``CVodeWFtolerances`` specifies a user-supplied function ``efun``  that sets the multiplicative error weights W_i for use in the weighted RMS norm, which are normally defined by :eq:`CVODES_errwt`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``efun`` -- is the C function which defines the ``ewt`` vector (see :c:type:`CVEwtFn`).

   **Return value:**
     * ``CV_SUCCESS`` -- The call was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized.
     * ``CV_NO_MALLOC`` -- The allocation function returned ``NULL``.


.. _CVODES.Usage.SIM.user_callable.toladvice:

General advice on choice of tolerances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For many users, the appropriate choices for tolerance values in ``reltol`` and ``abstol`` are a
concern. The following pieces of advice are relevant.

(1) The scalar relative tolerance ``reltol`` is to be set to control relative errors. So
:math:`\texttt{reltol} = 10^{-4}` means that errors are controlled to .01%. We do not recommend
using ``reltol`` larger than :math:`10^{-3}`. On the other hand, ``reltol`` should not be so small
that it is comparable to the unit roundoff of the machine arithmetic (generally
around :math:`10^{-15}`).

(2) The absolute tolerances ``abstol`` (whether scalar or vector) need to be set to control
absolute errors when any components of the solution vector ``y`` may be so small that
pure relative error control is meaningless. For example, if ``y[i]`` starts at some
nonzero value, but in time decays to zero, then pure relative error control on ``y[i]``
makes no sense (and is overly costly) after ``y[i]`` is below some noise level. Then
``abstol`` (if scalar) or ``abstol[i]`` (if a vector) needs to be set to that noise level. If the different
components have different noise levels, then ``abstol`` should be a vector. See the example  ``cvsRoberts_dns``
in the CVODES package, and the discussion of it in the CVODES Examples document
:cite:p:`cvodes_ex`. In that problem, the three components vary betwen 0 and 1,
and have different noise levels; hence the ``abstol`` vector. It is impossible to give any
general advice on ``abstol`` values, because the appropriate noise levels are completely
problem-dependent. The user or modeler hopefully has some idea as to what those
noise levels are.

(3) Finally, it is important to pick all the tolerance values conservatively,
because they control the error committed on each individual time step. The final
(global) errors are some sort of accumulation of those per-step errors. A good
rule of thumb is to reduce the tolerances by a factor of .01 from the actual
desired limits on errors. So if you want .01% accuracy (globally), a good choice
is :math:`\texttt{reltol} = 10^{-6}`. But in any case, it is a good idea to do a few experiments
with the tolerances to see how the computed solution values vary as tolerances
are reduced.

.. _CVODES.Usage.SIM.user_callable.unphysical:

Advice on controlling unphysical negative values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In many applications, some components in the true solution are always positive
or non-negative, though at times very small. In the numerical solution, however,
small negative (hence unphysical) values can then occur. In most cases, these
values are harmless, and simply need to be controlled, not eliminated. The
following pieces of advice are relevant.

(1) The way to control the size of unwanted negative computed values is with
tighter absolute tolerances. Again this requires some knowledge of the noise
level of these components, which may or may not be different for different
components. Some experimentation may be needed.

(2) If output plots or tables are being generated, and it is important to avoid
having negative numbers appear there (for the sake of avoiding a long
explanation of them, if nothing else), then eliminate them, but only in the
context of the output medium. Then the internal values carried by the solver are
unaffected. Remember that a small negative value in ``y`` returned by CVODES, with
magnitude comparable to ``abstol`` or less, is equivalent to zero as far as the computation
is concerned.

(3) The user’s right-hand side routine ``f`` should never change a negative value in
the solution vector ``y`` to a non-negative value, as a "solution" to this problem.
This can cause instability. If the ``f`` routine cannot tolerate a zero or negative
value (e.g. because there is a square root or log of it), then the offending
value should be changed to zero or a tiny positive number in a temporary
variable (not in the input ``y`` vector) for the purposes of computing :math:`f(t,y)`.

(4) Positivity and non-negativity constraints on components can be enforced by
use of the recoverable error return feature in the user-supplied right-hand side
function. However, because this option involves some extra overhead cost, it
should only be exercised if the use of absolute tolerances to control the
computed values is unsuccessful.

.. _CVODES.Usage.SIM.user_callable.lin_solv_init:

Linear solver interface functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As previously explained, if the nonlinear solver requires the solution
of linear systems of the form :eq:`CVODES_Newton` (e.g., the
default Newton iteration), there are two CVODES linear solver
interfaces currently available for this task: CVLS and CVDIAG.

The first corresponds to the main linear solver interface in CVODES,
that supports all valid ``SUNLinearSolver`` modules. Here, matrix-based
``SUNLinearSolver`` modules utilize ``SUNMatrix`` objects to store the
approximate Jacobian matrix :math:`J = \partial{f}/\partial{y}`, the
Newton matrix :math:`M = I-\gamma J`, and factorizations used throughout
the solution process. Conversely, matrix-free ``SUNLinearSolver`` modules
instead use iterative methods to solve the Newton systems of equations,
and only require the *action* of the matrix on a vector, :math:`Mv`.
With most of these methods, preconditioning can be done on the left
only, the right only, on both the left and right, or not at all. The
exceptions to this rule are SPFGMR that supports right
preconditioning only and PCG that performs symmetric
preconditioning. For the specification of a preconditioner, see the
iterative linear solver sections in :numref:`CVODES.Usage.SIM.optional_input`
and :numref:`CVODES.Usage.SIM.user_supplied`.

If preconditioning is done, user-supplied functions define linear
operators corresponding to left and right preconditioner matrices
:math:`P_1` and :math:`P_2` (either of which could be the identity
matrix), such that the product :math:`P_1 P_2` approximates the matrix
:math:`M = I - \gamma J` of :eq:`CVODES_Newtonmat`.

The CVDIAG linear solver interface supports a direct linear solver,
that uses only a diagonal approximation to :math:`J`.

To specify a generic linear solver to CVODES, after the call to :c:func:`CVodeCreate` but
before any calls to :c:func:`CVode`, the user’s program must create the appropriate
:c:type:`SUNLinearSolver` object and call the function :c:func:`CVodeSetLinearSolver`, as documented below. To create the
:c:type:`SUNLinearSolver` object, the user may call one of the SUNDIALS-packaged ``SUNLinearSolver``
module constructor routines via a call of the form ``SUNLinearSolver LS = SUNLinSol_*(...);``

Alternately, a user-supplied :c:type:`SUNLinearSolver` module may be created and used instead. The
use of each of the generic linear solvers involves certain constants,
functions and possibly some macros, that are likely to be needed in the
user code. These are available in the corresponding header file
associated with the specific ``SUNMatrix`` or ``SUNLinearSolver`` module in
question, as described in :numref:`SUNMatrix` and
:numref:`SUNLinSol`.

Once this solver object has been constructed, the user should attach it
to CVODES via a call to :c:func:`CVodeSetLinearSolver`. The first argument passed to this function
is the CVODES memory pointer returned by :c:func:`CVodeCreate`; the second argument is the
desired ``SUNLinearSolver`` object to use for solving linear systems. The
third argument is an optional ``SUNMatrix`` object to accompany
matrix-based ``SUNLinearSolver`` inputs (for matrix-free linear solvers, the
third argument should be ``NULL``). A call to this function initializes the
CVLS linear solver interface, linking it to the main CVODES
integrator, and allows the user to specify additional parameters and
routines pertinent to their choice of linear solver.

To instead specify the CVODES-specific diagonal linear solver
interface, the user’s program must call :c:func:`CVDiag`, as documented below. The first
argument passed to this function is the CVODES memory pointer
returned by :c:func:`CVodeCreate`.

.. c:function:: int CVodeSetLinearSolver(void* cvode_mem, SUNLinearSolver LS, SUNMatrix J)

   The function ``CVodeSetLinearSolver`` attaches a generic ``SUNLinearSolver``  object ``LS`` and corresponding template Jacobian ``SUNMatrix``  object ``J`` (if applicable) to CVODES, initializing the  CVLS linear solver interface.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``LS`` -- ``SUNLinearSolver`` object to use for solving linear systems of the form :eq:`CVODES_Newton`.
     * ``J`` -- ``SUNMatrix`` object for used as a template for the Jacobian (or ``NULL`` if not applicable).

   **Return value:**
     * ``CVLS_SUCCESS`` -- The CVLS initialization was successful.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_ILL_INPUT`` -- The CVLS interface is not compatible with the ``LS`` or ``J`` input objects or is incompatible with the current ``N_Vector`` module.
     * ``CVLS_SUNLS_FAIL`` -- A call to the ``LS`` object failed.
     * ``CVLS_MEM_FAIL`` -- A memory allocation request failed.

   **Notes:**
      If ``LS`` is a matrix-based linear solver, then the template  Jacobian matrix ``J`` will be used in the solve process, so if additional storage is required within the ``SUNMatrix`` object  (e.g. for factorization of a banded matrix), ensure that the input object is allocated with sufficient size (see :numref:`SUNMatrix` for further information).

      When using sparse linear solvers, it is typically much more  efficient to supply ``J`` so that it includes the full sparsity  pattern of the Newton system matrices :math:`M=I-\gamma J`, even if ``J``  itself has zeros in nonzero locations of I.  The reasoning for  this is that M is constructed in-place, on top of the  user-specified values of ``J``, so if the sparsity pattern in  ``J`` is insufficient to store M then it will need to be resized  internally by CVODES.

      The previous routines ``CVDlsSetLinearSolver`` and  ``CVSpilsSetLinearSolver`` are now wrappers for this routine, and may  still be used for backward-compatibility.  However, these will be  deprecated in future releases, so we recommend that users transition  to the new routine name soon.

.. c:function:: int CVDiag(void* cvode_mem)

   The function ``CVDiag`` selects the CVDIAG linear solver.  The user's main program must include the ``cvode_diag.h`` header file.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.

   **Return value:**
     * ``CVDIAG_SUCCESS`` -- The CVDIAG initialization was successful.
     * ``CVDIAG_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVDIAG_ILL_INPUT`` -- The CVDIAG solver is not compatible with the current ``N_Vector`` module.
     * ``CVDIAG_MEM_FAIL`` -- A memory allocation request failed.

   **Notes:**
      The CVDIAG solver is the simplest of all of the available CVODES  linear solvers.  The CVDIAG solver uses an approximate  diagonal Jacobian formed by way of a difference quotient. The user  does *not* have the option of supplying a function to compute an  approximate diagonal Jacobian.


.. _CVODES.Usage.SIM.nonlin_solv_init:

Nonlinear solver interface function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default CVODES uses the ``SUNNonlinearSolver`` implementation of Newton’s
method defined by the :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>` module.
To specify a different nonlinear solver in CVODES, the user’s program must create
a ``SUNNonlinearSolver`` object by calling the appropriate constructor routine.
The user must then attach the ``SUNNonlinearSolver`` object by calling
:c:func:`CVodeSetNonlinearSolver`, as documented below.

When changing the nonlinear solver in CVODES, :c:func:`CVodeSetNonlinearSolver`
must be called after :c:func:`CVodeInit`. If any calls to :c:func:`CVode` have been made,
then CVODES will need to be reinitialized by calling :c:func:`CVodeReInit` to ensure
that the nonlinear solver is initialized correctly before any subsequent calls
to :c:func:`CVode`.

The first argument passed to the routine :c:func:`CVodeSetNonlinearSolver` is
the CVODES memory pointer returned by :c:func:`CVodeCreate` and the second
argument is the ``SUNNonlinearSolver`` object to use for solving the nonlinear
system :eq:`CVODES_Newton` or :eq:`CVODES_nonlinear_fixedpoint`. A call to this function
attaches the nonlinear solver to the main CVODES integrator.

.. c:function:: int CVodeSetNonlinearSolver(void* cvode_mem, SUNNonlinearSolver NLS)

   The function ``CVodeSetNonLinearSolver`` attaches a ``SUNNonlinearSolver``  object (``NLS``) to CVODES.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``NLS`` -- ``SUNNonlinearSolver`` object to use for solving nonlinear systems :eq:`CVODES_nonlinear` or :eq:`CVODES_nonlinear_fixedpoint`.

   **Return value:**
     * ``CV_SUCCESS`` -- The nonlinear solver was successfully attached.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- The ``SUNNonlinearSolver`` object is ``NULL``, does not implement
       the required nonlinear solver operations, is not of the correct type, or the residual
       function, convergence test function, or maximum number of nonlinear iterations could
       not be set.

   **Notes:**
     When forward sensitivity analysis capabilities are enabled and the
     ``CV_STAGGERED`` or ``CV_STAGGERED1`` corrector method is used this
     function sets the nonlinear solver method for correcting state variables
     (see :numref:`CVODES.Usage.FSA.user_callable.nonlin_solv_init` for more details).


.. _CVODES.Usage.SIM.cvrootinit:

Rootfinding initialization function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While solving the IVP, CVODES has the capability to find the roots of
a set of user-defined functions. To activate the root finding algorithm,
call the following function. This is normally called only once, prior to
the first call to :c:func:`CVode`, but if the rootfinding problem is to be changed
during the solution, :c:func:`CVodeRootInit` can also be called prior to a continuation call to :c:func:`CVode`.

.. c:function:: int CVodeRootInit(void* cvode_mem, int nrtfn, CVRootFn g)

   The function ``CVodeRootInit`` specifies that the roots of a set of functions :math:`g_i(t,y)` are to be found while the IVP is being solved.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``nrtfn`` -- is the number of root functions :math:`g_i`.
     * ``g`` -- is the C function which defines the ``nrtfn`` functions :math:`g_i(t,y)`
       whose roots are sought. See :numref:`CVODES.Usage.SIM.user_supplied.rootFn` for details.

   **Return value:**
     * ``CV_SUCCESS`` -- The call was successful.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` argument was ``NULL``.
     * ``CV_MEM_FAIL`` -- A memory allocation failed.
     * ``CV_ILL_INPUT`` -- The function ``g`` is ``NULL``, but ``nrtfn>0``.

   **Notes:**
      If a new IVP is to be solved with a call to ``CVodeReInit``, where the new  IVP has no rootfinding problem but the prior one did, then call  ``CVodeRootInit`` with ``nrtfn=0``.


CVODES solver function
~~~~~~~~~~~~~~~~~~~~~~

This is the central step in the solution process — the call to perform
the integration of the IVP. One of the input arguments (``itask``) specifies one
of two modes as to where CVODES is to return a solution. But these
modes are modified if the user has set a stop time (with :c:func:`CVodeSetStopTime`) or requested
rootfinding.

.. c:function:: int CVode(void* cvode_mem, realtype tout, N_Vector yout, realtype tret, int itask)

   The function ``CVode`` integrates the ODE over an interval in t.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``tout`` -- the next time at which a computed solution is desired.
     * ``yout`` -- the computed solution vector.
     * ``tret`` -- the time reached by the solver (output).
     * ``itask`` --  a flag indicating the job of the solver for the next user step. The ``CV_NORMAL`` option causes the solver to take internal steps until it has reached or just passed the user-specified ``tout`` parameter. The solver then interpolates in order to return an approximate value of :math:`y({tout})`. The ``CV_ONE_STEP`` option tells the solver to take just one internal step and then return the solution at the point reached by that step.

   **Return value:**
     * ``CV_SUCCESS`` -- ``CVode`` succeeded and no roots were found.
     * ``CV_TSTOP_RETURN`` -- ``CVode`` succeeded by reaching the stopping point specified through the optional input function :c:func:`CVodeSetStopTime`.
     * ``CV_ROOT_RETURN`` -- ``CVode`` succeeded and found one or more roots.  In this case, ``tret`` is the location of the root.  If ``nrtfn`` :math:`>1`, call :c:func:`CVodeGetRootInfo` to see which :math:`g_i` were found to have a root.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_NO_MALLOC`` -- The CVODES memory was not allocated by a call to :c:func:`CVodeInit`.
     * ``CV_ILL_INPUT`` -- One of the inputs to ``CVode`` was illegal, or some other input to the solver was illegal or missing. The latter category includes the following situations:

       (a) The tolerances have not been set.

       (b) A component of the error weight vector became zero during internal time-stepping.

       (c) The linear solver initialization function (called by the user after calling :c:func:`CVodeCreate`) failed to set the linear solver-specific ``lsolve`` field in ``cvode_mem``.

       (d) A root of one of the root functions was found both at a point :math:`t` and also very near :math:`t`.

     * ``CV_TOO_CLOSE`` -- The initial time :math:`t_0` and the output time :math:`t_{out}` are too close to each other and the user did not specify an initial step size.
     * ``CV_TOO_MUCH_WORK`` -- The solver took ``mxstep`` internal steps but still could not reach ``tout``. The default value for ``mxstep`` is ``MXSTEP_DEFAULT = 500``.
     * ``CV_TOO_MUCH_ACC`` -- The solver could not satisfy the accuracy demanded by the user for some internal step.
     * ``CV_ERR_FAILURE`` -- Either error test failures occurred too many times (``MXNEF = 7``) during one internal time step, or with :math:`|h| = h_{min}`.
     * ``CV_CONV_FAILURE`` -- Either convergence test failures occurred too many times (``MXNCF = 10``) during one internal time step, or with :math:`|h| = h_{min}`.
     * ``CV_LINIT_FAIL`` -- The linear solver interface's initialization function failed.
     * ``CV_LSETUP_FAIL`` -- The linear solver interface's setup function failed in an unrecoverable manner.
     * ``CV_LSOLVE_FAIL`` -- The linear solver interface's solve function failed in an unrecoverable manner.
     * ``CV_CONSTR_FAIL`` -- The inequality constraints were violated and the solver was unable to recover.
     * ``CV_RHSFUNC_FAIL`` -- The right-hand side function failed in an unrecoverable manner.
     * ``CV_FIRST_RHSFUNC_FAIL`` -- The right-hand side function had a recoverable error at the first call.
     * ``CV_REPTD_RHSFUNC_ERR`` -- Convergence test failures occurred too many times due to repeated recoverable errors in the right-hand side function. This flag will also be returned if the right-hand side function had repeated recoverable errors during the estimation of an initial step size.
     * ``CV_UNREC_RHSFUNC_ERR`` -- The right-hand function had a recoverable error, but no recovery was possible.    This failure mode is rare, as it can occur only if the right-hand side function fails recoverably after an error test failed while at order one.
     * ``CV_RTFUNC_FAIL`` -- The rootfinding function failed.

   **Notes:**
      The vector ``yout`` can occupy the same space as the vector ``y0`` of  initial conditions that was passed to ``CVodeInit``.

      In the ``CV_ONE_STEP`` mode, ``tout`` is used only on the first call,  and only to get the direction and a rough scale of the independent variable.

      If a stop time is enabled (through a call to ``CVodeSetStopTime``), then  ``CVode`` returns the solution at ``tstop``. Once the integrator returns  at a stop time, any future testing for ``tstop`` is disabled (and can be  reenabled only though a new call to ``CVodeSetStopTime``).

      All failure return values are negative and so the test ``flag < 0``  will trap all ``CVode`` failures.

      On any error return in which one or more internal steps were taken by  ``CVode``, the returned values of ``tret`` and ``yout`` correspond to  the farthest point reached in the integration.  On all other error returns,  ``tret`` and ``yout`` are left unchanged from the previous ``CVode``  return.


.. _CVODES.Usage.SIM.optional_input:

Optional input functions
~~~~~~~~~~~~~~~~~~~~~~~~

There are numerous optional input parameters that control the behavior
of the CVODES solver. CVODES provides functions that can be used
to change these optional input parameters from their default values.
:numref:`CVODES.Usage.SIM.optional_input.Table` lists all optional input functions in
CVODES which are then described in detail in the remainder of this
section, begining with those for the main CVODES solver and
continuing with those for the linear solver interfaces. Note that the
diagonal linear solver module has no optional inputs. For the most
casual use of CVODES, the reader can skip to
:numref:`CVODES.Usage.SIM.user_supplied`.

We note that, on an error return, all of the optional input functions
send an error message to the error handler function. All error return
values are negative, so the test will catch all errors. Finally, a call
to a function can be made from the user’s calling program at any time
and, if successful, takes effect immediately.

.. _CVODES.Usage.SIM.optional_input.Table:
.. table:: Optional inputs for CVODES and CVLS

   +-------------------------------+---------------------------------------------+----------------+
   |      **Optional input**       |              **Function name**              |  **Default**   |
   +===============================+=============================================+================+
   | **CVODES main solver**        |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Pointer to an error file      | :c:func:`CVodeSetErrFile`                   | ``stderr``     |
   +-------------------------------+---------------------------------------------+----------------+
   | Error handler function        | :c:func:`CVodeSetErrHandlerFn`              | internal fn.   |
   +-------------------------------+---------------------------------------------+----------------+
   | User data                     | :c:func:`CVodeSetUserData`                  | ``NULL``       |
   +-------------------------------+---------------------------------------------+----------------+
   | Maximum order for BDF method  | :c:func:`CVodeSetMaxOrd`                    | 5              |
   +-------------------------------+---------------------------------------------+----------------+
   | Maximum order for Adams       | :c:func:`CVodeSetMaxOrd`                    | 12             |
   | method                        |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Maximum no. of internal steps | :c:func:`CVodeSetMaxNumSteps`               | 500            |
   | before :math:`t_{out}`        |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Maximum no. of warnings for   | :c:func:`CVodeSetMaxHnilWarns`              | 10             |
   | :math:`t_n+h=t_n`             |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Flag to activate stability    | :c:func:`CVodeSetStabLimDet`                | ``SUNFALSE``   |
   | limit detection               |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Initial step size             | :c:func:`CVodeSetInitStep`                  | estimated      |
   +-------------------------------+---------------------------------------------+----------------+
   | Minimum absolute step size    | :c:func:`CVodeSetMinStep`                   | 0.0            |
   +-------------------------------+---------------------------------------------+----------------+
   | Maximum absolute step size    | :c:func:`CVodeSetMaxStep`                   | :math:`\infty` |
   +-------------------------------+---------------------------------------------+----------------+
   | Value of :math:`t_{stop}`     | :c:func:`CVodeSetStopTime`                  | undefined      |
   +-------------------------------+---------------------------------------------+----------------+
   | Maximum no. of error test     | :c:func:`CVodeSetMaxErrTestFails`           | 7              |
   | failures                      |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Maximum no. of nonlinear      | :c:func:`CVodeSetMaxNonlinIters`            | 3              |
   | iterations                    |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Maximum no. of convergence    | :c:func:`CVodeSetMaxConvFails`              | 10             |
   | failures                      |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Coefficient in the nonlinear  | :c:func:`CVodeSetNonlinConvCoef`            | 0.1            |
   | convergence test              |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Inequality constraints on     | :c:func:`CVodeSetConstraints`               |                |
   | solution                      |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Direction of zero-crossing    | :c:func:`CVodeSetRootDirection`             | both           |
   +-------------------------------+---------------------------------------------+----------------+
   | Disable rootfinding warnings  | :c:func:`CVodeSetNoInactiveRootWarn`        | none           |
   +-------------------------------+---------------------------------------------+----------------+
   | **CVLS linear solver          |                                             |                |
   | interface**                   |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Linear solver setup frequency | :c:func:`CVodeSetLSetupFrequency`           | 20             |
   +-------------------------------+---------------------------------------------+----------------+
   | Jacobian / preconditioner     | :c:func:`CVodeSetJacEvalFrequency`          | 51             |
   | update frequency              |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Jacobian function             | :c:func:`CVodeSetJacFn`                     | DQ             |
   +-------------------------------+---------------------------------------------+----------------+
   | Linear System function        | :c:func:`CVodeSetLinSysFn`                  | internal       |
   +-------------------------------+---------------------------------------------+----------------+
   | Enable or disable linear      | :c:func:`CVodeSetLinearSolutionScaling`     | on             |
   | solution scaling              |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Jacobian-times-vector         | :c:func:`CVodeSetJacTimes`                  | NULL, DQ       |
   | functions                     |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Jacobian-times-vector DQ RHS  | :c:func:`CVodeSetJacTimesRhsFn`             | NULL           |
   | function                      |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Preconditioner functions      | :c:func:`CVodeSetPreconditioner`            | NULL, NULL     |
   +-------------------------------+---------------------------------------------+----------------+
   | Ratio between linear and      | :c:func:`CVodeSetEpsLin`                    | 0.05           |
   | nonlinear tolerances          |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+
   | Newton linear solve tolerance | :c:func:`CVodeSetLSNormFactor`              | vector length  |
   | conversion factor             |                                             |                |
   +-------------------------------+---------------------------------------------+----------------+

.. _CVODES.Usage.SIM.optional_input.optin_main:

Main solver optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The calls listed here can be executed in any order. However, if either
of the functions or is to be called, that call should be first, in order
to take effect for any later error message.

.. c:function:: int CVodeSetErrFile(void* cvode_mem, FILE * errfp)

   The function ``CVodeSetErrFile`` specifies a pointer to the file  where all CVODES messages should be directed when the default  CVODES error handler function is used.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``errfp`` -- pointer to output file.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      The default value for ``errfp`` is ``stderr``.  Passing a value of ``NULL`` disables all future error message output  (except for the case in which the CVODES memory pointer is ``NULL``).  This use of ``CVodeSetErrFile`` is strongly discouraged.

      .. warning::

         If ``CVodeSetErrFile`` is to be called, it should be called before any  other optional input functions, in order to take effect for any later error message.

.. c:function:: int CVodeSetErrHandlerFn(void* cvode_mem, CVErrHandlerFn ehfun, void * eh_data)

   The function ``CVodeSetErrHandlerFn`` specifies the optional user-defined function  to be used in handling error messages.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``ehfun`` -- is the C error handler function.
     * ``eh_data`` -- pointer to user data passed to ``ehfun`` every time it is called.

   **Return value:**
     * ``CV_SUCCESS`` -- The function ``ehfun`` and data pointer ``eh_data`` have been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      Error messages indicating that the CVODES solver memory is ``NULL`` will  always be directed to ``stderr``.

.. c:function:: int CVodeSetUserData(void* cvode_mem, void * user_data)

   The function ``CVodeSetUserData`` specifies the user data block ``user_data``  and attaches it to the main CVODES memory block.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``user_data`` -- pointer to the user data.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      If specified, the pointer to ``user_data`` is passed to all user-supplied functions that have it as an argument. Otherwise, a ``NULL`` pointer is passed.

    .. warning::

      If ``user_data`` is needed in user linear solver or preconditioner functions, the call to ``CVodeSetUserData`` must be made before the call to specify the linear solver.

.. c:function:: int CVodeSetMonitorFn(void* cvode_mem, CVMonitorFn monitorfn)

   The function ``CVodeSetMonitorFn`` specifies a user function,  ``monitorfn``, to be called at some interval of successfully  completed CVODES time steps.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``monitorfn`` -- user-supplied monitor function (``NULL`` by default); a ``NULL`` input will turn off monitoring.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      The frequency with which the monitor function is called can be  set with the function ``CVodeSetMonitorFrequency``.

      .. warning::

         Modifying the solution in this function will result in  undefined behavior. This function is only intended to be used  for monitoring the integrator.  SUNDIALS must be built with the CMake option  ``SUNDIALS_BUILD_WITH_MONITORING``, to utilize this function.  See :numref:`Installation` for more information.

.. c:function:: int CVodeSetMonitorFrequency(void* cvode_mem, long int nst)

   The function ``CVodeSetMonitorFrequency`` specifies the  interval, measured in successfully completed CVODES time-steps,  at which the monitor function should be called.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nst`` -- number of successful steps inbetween calls to the monitor function 0 by default;    a 0 input will turn off monitoring.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized :c:func:`CVodeCreate`.

   **Notes:**
      The monitor function that will be called can be set with  ``CVodeSetMonitorFn``.

      .. warning::

         Modifying the solution in this function will result in undefined behavior. This function is only intended to be used for monitoring the integrator.  SUNDIALS must be built with the CMake option  ``SUNDIALS_BUILD_WITH_MONITORING``, to utilize this function.  See :numref:`Installation` for more information.

.. c:function:: int CVodeSetMaxOrd(void* cvode_mem, int maxord)

   The function ``CVodeSetMaxOrd`` specifies the maximum order of the  linear multistep method.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``maxord`` -- value of the maximum method order.  This must be positive.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- The specified value ``maxord`` is :math:`\leq 0`, or larger than its previous value.

   **Notes:**
      The default value is ``ADAMS_Q_MAX = 12`` for  the Adams-Moulton method and ``BDF_Q_MAX = 5``  for the BDF method.  Since ``maxord`` affects the memory requirements  for the internal CVODES memory block, its value  cannot be increased past its previous value.

      An input value greater than the default will result in the default value.

.. c:function:: int CVodeSetMaxNumSteps(void* cvode_mem, long int mxsteps)

   The function ``CVodeSetMaxNumSteps`` specifies the maximum number  of steps to be taken by the solver in its attempt to reach  the next output time.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``mxsteps`` -- maximum allowed number of steps.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      Passing ``mxsteps`` = 0 results in CVODES using the default value (500).

      Passing ``mxsteps`` < 0 disables the test (not recommended).

.. c:function:: int CVodeSetMaxHnilWarns(void* cvode_mem, int mxhnil)

   The function ``CVodeSetMaxHnilWarns`` specifies the maximum number of  messages issued by the solver warning that :math:`t+h=t` on the next internal step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``mxhnil`` -- maximum number of warning messages :math:`(> 0)`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      The default value is 10.  A negative value for ``mxhnil`` indicates that no warning messages should  be issued.

.. c:function:: int CVodeSetStabLimDet(void* cvode_mem, booleantype stldet)

   The function ``CVodeSetStabLimDet`` indicates if  the BDF stability limit detection algorithm should be used. See :numref:`CVODES.Mathematics.stablimit` for further details.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``stldet`` -- flag controlling stability limit detection (``SUNTRUE`` = on; ``SUNFALSE`` = off).

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- The linear multistep method is not set to ``CV_BDF``.

   **Notes:**
      The default value is ``SUNFALSE``. If ``stldet = SUNTRUE`` when BDF is used  and the method order is greater than or equal to 3, then an internal function, ``CVsldet``,  is called to detect a possible stability limit. If such a limit is detected, then the order is  reduced.

.. c:function:: int CVodeSetInitStep(void* cvode_mem, realtype hin)

   The function ``CVodeSetInitStep`` specifies the initial step size.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``hin`` -- value of the initial step size to be attempted. Pass 0.0 to use the default value.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      By default, CVODES estimates the initial step size to be the solution :math:`h` of the equation :math:`0.5 h^2 \ddot{y} = 1`,  where :math:`\ddot{y}` is an estimated second derivative of the solution at :math:`t_0`.

.. c:function:: int CVodeSetMinStep(void* cvode_mem, realtype hmin)

   The function ``CVodeSetMinStep`` specifies a lower bound on the magnitude  of the step size.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``hmin`` -- minimum absolute value of the step size :math:`(\geq 0.0)`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- Either ``hmin`` is nonpositive or it exceeds the maximum allowable step size.

   **Notes:**
      The default value is 0.0.

.. c:function:: int CVodeSetMaxStep(void* cvode_mem, realtype hmax)

   The function ``CVodeSetMaxStep`` specifies an upper bound on the magnitude  of the step size.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``hmax`` -- maximum absolute value of the step size :math:`( \geq 0.0 )`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- Either ``hmax`` is nonpositive or it is smaller than the minimum allowable step size.

   **Notes:**
      Pass ``hmax`` = 0.0 to obtain the default value :math:`\infty`.

.. c:function:: int CVodeSetStopTime(void* cvode_mem, realtype tstop)

   The function ``CVodeSetStopTime`` specifies the value of the  independent variable :math:`t` past which the solution is not to proceed.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``tstop`` -- value of the independent variable past which the solution should    not proceed.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- The value of ``tstop`` is not beyond the current :math:`t` value, :math:`t_n`.

   **Notes:**
      The default, if this routine is not called, is that no stop time is imposed.

      Once the integrator returns at a stop time, any future testing for ``tstop``  is disabled (and can be reenabled only though a new call to ``CVodeSetStopTime``).

.. c:function:: int CVodeSetMaxErrTestFails(void* cvode_mem, int maxnef)

   The function ``CVodeSetMaxErrTestFails`` specifies the  maximum number of error test failures permitted in attempting one step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``maxnef`` -- maximum number of error test failures allowed on one step :math:`(> 0)`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      The default value is 7.

.. c:function:: int CVodeSetMaxNonlinIters(void* cvode_mem, int maxcor)

   The function ``CVodeSetMaxNonlinIters`` specifies the maximum  number of nonlinear solver iterations permitted per step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``maxcor`` -- maximum number of nonlinear solver iterations allowed per step :math:`(> 0)`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_MEM_FAIL`` -- The ``SUNNonlinearSolver`` module is ``NULL``.

   **Notes:**
      The default value is 3.

.. c:function:: int CVodeSetMaxConvFails(void* cvode_mem, int maxncf)

   The function ``CVodeSetMaxConvFails`` specifies the  maximum number of nonlinear solver convergence failures permitted during  one step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``maxncf`` -- maximum number of allowable nonlinear solver convergence failures per step :math:`(> 0)`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      The default value is 10.

.. c:function:: int CVodeSetNonlinConvCoef(void* cvode_mem, realtype nlscoef)

   The function ``CVodeSetNonlinConvCoef`` specifies the safety factor used in the nonlinear convergence test (see :numref:`CVODES.Mathematics.ivp_sol`).

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nlscoef`` -- coefficient in nonlinear convergence test :math:`(> 0)`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      The default value is 0.1.

.. c:function:: int CVodeSetNlsRhsFn(void* cvode_mem, CVRhsFn f)

   The function ``CVodeSetNlsRhsFn`` specifies an alternative right-hand side function for use in nonlinear system function evaluations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``f`` -- is the alternative C function which computes the right-hand side function :math:`f` in the ODE (for full details see :numref:`CVODES.Usage.SIM.user_supplied.rhsFn`).

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      The default is to use the implicit right-hand side function provided to
      :c:func:`CVodeInit` in nonlinear system function evaluations. If the input
      right-hand side function is ``NULL``, the default is used.

      When using a non-default nonlinear solver, this function must be called after
      :c:func:`CVodeSetNonlinearSolver`.

      When doing forward sensitivity analysis with the simultaneous solver strategy
      and a non-default nonlinear solver, this function must be called after
      :c:func:`CVodeSetNonlinearSolverSensSim`.


.. c:function:: int CVodeSetConstraints(void* cvode_mem, N_Vector constraints)

   The function ``CVodeSetConstraints`` specifies a vector defining  inequality constraints for each component of the solution vector y.

   **Arguments:**
      * ``cvode_mem`` -- pointer to the CVODES memory block.
      * ``constraints`` -- vector of constraint flags. If ``constraints[i]`` is

        * 0.0 then no constraint is imposed on :math:`y_i`.
        * 1.0 then :math:`y_i` will be constrained to be :math:`y_i \ge 0.0`.
        * -1.0  then :math:`y_i` will be constrained to be :math:`y_i \le 0.0`.
        * 2.0 then :math:`y_i` will be constrained to be :math:`y_i > 0.0`.
        * -2.0  then :math:`y_i` will be constrained to be :math:`y_i < 0.0`.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- The constraints vector contains illegal values or the simultaneous corrector option has been selected when doing forward sensitivity analysis.

   **Notes:**
      The presence of a non-``NULL`` constraints vector that is not 0.0 in  all
      components will cause constraint checking to be performed. However, a call
      with 0.0 in all components of ``constraints`` will  result in an illegal
      input return. A ``NULL`` constraints vector will disable  constraint
      checking.

      Constraint checking when doing forward sensitivity analysis with the
      simultaneous corrector option is currently disallowed and will result in an
      illegal input return.


.. _CVODES.Usage.SIM.optional_inputs.optin_ls:

Linear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mathematical explanation of the linear solver methods available to
CVODES is provided in :numref:`CVODES.Mathematics.ivp_sol`. We group the
user-callable routines into four categories: general routines concerning
the overall CVLS linear solver interface, optional inputs for
matrix-based linear solvers, optional inputs for matrix-free linear
solvers, and optional inputs for iterative linear solvers. We note that
the matrix-based and matrix-free groups are mutually exclusive, whereas
the “iterative” tag can apply to either case.

As discussed in :numref:`CVODES.Mathematics.ivp_sol`, CVODES strives to
reuse matrix and preconditioner data for as many solves as possible to
amortize the high costs of matrix construction and factorization. To
that end, CVODES provides user-callable routines to modify this
behavior. Recall that the Newton system matrices are
:math:`M(t,y) = I - \gamma J(t,y)`, where the right-hand side function
has Jacobian matrix :math:`J(t,y) = \dfrac{\partial f(t,y)}{\partial y}`.

The matrix or preconditioner for :math:`M` can only be updated within a
call to the linear solver ‘setup’ routine. In general, the frequency
with which this setup routine is called may be controlled with the ``msbp``
argument to :c:func:`CVodeSetLSetupFrequency`. When this occurs, the validity of :math:`M` for successive
time steps intimately depends on whether the corresponding
:math:`\gamma` and :math:`J` inputs remain valid.

At each call to the linear solver setup routine the decision to update
:math:`M` with a new value of :math:`\gamma`, and to reuse or reevaluate
Jacobian information, depends on several factors including:

-  the success or failure of previous solve attempts,
-  the success or failure of the previous time step attempts,
-  the change in :math:`\gamma` from the value used when constructing
   :math:`M`, and
-  the number of steps since Jacobian information was last evaluated.

The frequency with which to update Jacobian information can be controlled with
the ``msbj`` argument to :c:func:`CVodeSetJacEvalFrequency`. We note that this
is only checked *within* calls to the linear solver setup routine, so values
:math:`<` do not make sense. For linear-solvers with user-supplied
preconditioning the above factors are used to determine whether to recommend
updating the Jacobian information in the preconditioner (i.e., whether to set
``jok`` to ``SUNFALSE`` in calling the user-supplied preconditioner setup
function (see :numref:`CVODES.Usage.SIM.user_supplied.precondFn`). For matrix-based
linear solvers these factors determine whether the matrix
:math:`J(t,y) = \dfrac{\partial f(t,y)}{\partial y}`
should be updated (either with an internal finite
difference approximation or a call to the user-supplied Jacobian function
(see :numref:`CVODES.Usage.SIM.user_supplied.jacFn`); if not then the previous value
is reused and the system matrix :math:`M(t,y) \approx I - \gamma J(t,y)` is recomputed
using the current :math:`\gamma` value.

.. c:function:: int CVodeSetLSetupFrequency(void* cvode_mem, long int msbp)

   The function ``CVodeSetLSetupFrequency`` specifies the frequency of  calls to the linear solver setup function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``msbp`` -- the linear solver setup frequency.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- The frequency ``msbp`` is negative.

   **Notes:**
      Positive values of ``msbp`` specify the linear solver setup frequency. For  example, an input of ``1`` means the setup function will be called every time  step while an input of ``2`` means it will be called called every other time  step. If ``msbp = 0``, the default value of 20 will be used. Otherwise an  error is returned.

.. c:function:: int CVodeSetJacEvalFrequency(void* cvode_mem, long int msbj)

   The function ``CVodeSetJacEvalFrequency`` specifies the frequency for  recomputing the Jacobian or recommending a preconditioner update.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``msbj`` -- the Jacobian re-computation or preconditioner update frequency.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver interface has not been initialized.
     * ``CVLS_ILL_INPUT`` -- The frequency ``msbj`` is negative.

   **Notes:**
      The Jacobian update frequency is only checked within calls to the  linear solver setup routine, as such values of ``msbj`` < ``msbp`` will  result in recomputing the Jacobian every ``msbp`` steps. See :c:func:`CVodeSetLSetupFrequency` for setting the linear solver setup frequency  ``msbp``.  If ``msbj = 0``, the default value of 51 will be used. Otherwise an error is returned. This function must be called after the CVLS linear solver interface has been initialized through a call to :c:func:`CVodeSetLinearSolver`.

When using matrix-based linear solver modules, the CVLS solver interface
needs a function to compute an approximation to the Jacobian matrix :math:`J(t,y)` or
the linear system :math:`M = I - \gamma J`. The function to evaluate :math:`J(t,y)` must
be of type :c:type:`CVLsJacFn`. The user can supply a Jacobian function, or if using
a :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or :ref:`SUNMATRIX_BAND <SUNMatrix.Band>`
matrix :math:`J`, can use the default internal difference quotient
approximation that comes with the CVLS solver. To specify a user-supplied Jacobian function
``jac``, CVLS provides the function :c:func:`CVodeSetJacFn`. The CVLS
interface passes the pointer ``user_data`` to the Jacobian function. This
allows the user to create an arbitrary structure with relevant problem data and
access it during the execution of the user-supplied Jacobian function, without
using global data in the program. The pointer ``user_data`` may be specified
through :c:func:`CVodeSetUserData`.

.. c:function:: int CVodeSetJacFn(void* cvode_mem, CVLsJacFn jac)

   The function ``CVodeSetJacFn`` specifies the Jacobian  approximation function to be used for a matrix-based solver within  the CVLS interface.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``jac`` -- user-defined Jacobian approximation function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver interface has not been initialized.

   **Notes:**
      This function must be called after the CVLS linear solver  interface has been initialized through a call to :c:func:`CVodeSetLinearSolver`.

      By default, CVLS uses an internal difference quotient function for the
      :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` and
      :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` modules.  If ``NULL`` is passed to
      ``jac``,  this default function is used.  An error will occur if no ``jac``
      is supplied when using other matrix types.

      The function type :c:type:`CVLsJacFn` is described in :numref:`CVODES.Usage.SIM.user_supplied.jacFn`.

      The previous routine ``CVDlsSetJacFn`` is now a wrapper for this  routine, and may still be used for backward-compatibility.  However, this will be deprecated in future releases, so we recommend that  users transition to the new routine name soon.


To specify a user-supplied linear system function ``linsys``, CVLS provides
the function :c:func:`CVodeSetLinSysFn`. The CVLS interface passes the pointer
``user_data`` to the linear system function. This allows the user to create an
arbitrary structure with relevant problem data and access it during the
execution of the user-supplied linear system function, without using global data
in the program. The pointer ``user_data`` may be specified through
:c:func:`CVodeSetUserData`.


.. c:function:: int CVodeSetLinSysFn(void* cvode_mem, CVLsLinSysFn linsys)

   The function ``CVodeSetLinSysFn`` specifies the linear system approximation  function to be used for a matrix-based solver within the CVLS interface.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``linsys`` -- user-defined linear system approximation function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver interface has not been initialized.

   **Notes:**
      This function must be called after the CVLS linear solver  interface has been initialized through a call to :c:func:`CVodeSetLinearSolver`.

      By default, CVLS uses an internal linear system function leveraging the  ``SUNMatrix`` API to form the system :math:`M = I - \gamma J` using either an  internal finite difference approximation or user-supplied function to compute the Jacobian. If ``linsys`` is ``NULL``, this default function is used.

      The function type :c:type:`CVLsLinSysFn` is described in :numref:`CVODES.Usage.SIM.user_supplied.jacFn`.

When using a matrix-based linear solver the matrix information will be updated
infrequently to reduce matrix construction and, with direct solvers,
factorization costs. As a result the value of :math:`\gamma` may not be current and,
with BDF methods, a scaling factor is applied to the solution of the linear
system to account for the lagged value of :math:`\gamma`. See
:numref:`SUNLinSol.CVODES.Lagged` for more details. The function
:c:func:`CVodeSetLinearSolutionScaling` can be used to disable this scaling when
necessary, e.g., when providing a custom linear solver that updates the matrix
using the current :math:`\gamma` as part of the solve.

.. c:function:: int CVodeSetLinearSolutionScaling(void* cvode_mem, booleantype onoff)

   The function :c:func:`CVodeSetLinearSolutionScaling` enables or disables scaling  the linear system solution to account for a change in :math:`\gamma` in the linear system. For more details see :numref:`SUNLinSol.CVODES.lagged`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``onoff`` -- flag to enable (``SUNTRUE``) or disable (``SUNFALSE``) scaling.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The flag value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver interface has not been initialized.
     * ``CVLS_ILL_INPUT`` -- The attached linear solver is not matrix-based or the linear multistep method type is not BDF.

   **Notes:**
      This function must be called after the CVLS linear solver  interface has been initialized through a call to  ``CVodeSetLinearSolver``.

      By default scaling is enabled with matrix-based linear solvers when using BDF  methods.

When using matrix-free linear solver modules, the CVLS solver
interface requires a function to compute an approximation to the
product between the Jacobian matrix :math:`J(t,y)` and a vector :math:`v`. The
user can supply a Jacobian-times-vector approximation function or use
the default internal difference quotient function
that comes with the CVLS interface.

A user-defined Jacobian-vector product
function must be of type :c:type:`CVLsJacTimesVecFn` and
can be specified through a call to :c:func:`CVodeSetJacTimes` (see
:numref:`CVODES.Usage.SIM.user_supplied.jtimesFn` for specification details).
The evaluation and processing of any Jacobian-related data needed by
the user's Jacobian-times-vector function may be done in the optional
user-supplied function ``jtsetup`` (see :numref:`CVODES.Usage.SIM.user_supplied.jtsetupFn` for
specification details).
The pointer ``user_data`` received through :c:func:`CVodeSetUserData` (or
a pointer to ``NULL`` if ``user_data`` was not specified)
is passed to the Jacobian-times-vector setup and product functions, ``jtsetup`` and
``jtimes``, each time they are called.  This allows the user to
create an arbitrary structure with relevant problem data and access it
during the execution of the user-supplied functions
without using global data in the program.


.. c:function:: int CVodeSetJacTimes(void* cvode_mem, CVLsJacTimesSetupFn jtsetup, CVLsJacTimesVecFn jtimes)

   The function ``CVodeSetJacTimes`` specifies the Jacobian-vector  setup and product functions.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``jtsetup`` -- user-defined Jacobian-vector setup function.
     * ``jtimes`` -- user-defined Jacobian-vector product function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_SUNLS_FAIL`` -- An error occurred when setting up the system matrix-times-vector routines in the ``SUNLinearSolver`` object used by the CVLS interface.

   **Notes:**
      The default is to use an internal finite difference quotient for  ``jtimes`` and to omit ``jtsetup``.  If ``NULL`` is passed to  ``jtimes``, these defaults are used.  A user may specify  non-``NULL`` ``jtimes`` and ``NULL`` ``jtsetup`` inputs.

      This function must be called after the CVLS linear solver  interface has been initialized through a call to  :c:func:`CVodeSetLinearSolver`.

      The previous routine ``CVSpilsSetJacTimes`` is now a wrapper for this routine, and may still be used for backward-compatibility.  However, this will be deprecated in future releases, so we recommend that users transition to the new routine name soon.


When using the internal difference quotient the user may optionally supply an
alternative right-hand side function for use in the Jacobian-vector product
approximation by calling :c:func:`CVodeSetJacTimesRhsFn`. The alternative right-hand
side function should compute a suitable (and differentiable) approximation to
the right-hand side function provided to :c:func:`CVodeInit`. For example, as done in
:cite:p:`dorr2010numerical`, the alternative function may use lagged values when
evaluating a nonlinearity in the right-hand side to avoid differencing a
potentially non-differentiable factor.

.. c:function:: int CVodeSetJacTimesRhsFn(void* cvode_mem, CVRhsFn jtimesRhsFn)

   The function ``CVodeSetJacTimesRhsFn`` specifies an alternative ODE  right-hand side function for use in the internal Jacobian-vector product  difference quotient approximation.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``jtimesRhsFn`` -- is the C function which computes the alternative ODE right-hand side function to use in Jacobian-vector product difference quotient approximations. This function has the form ``f(t, y, ydot, user_data)`` (for full details see :numref:`CVODES.Usage.SIM.user_supplied.rhsFn`).

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_ILL_INPUT`` -- The internal difference quotient approximation is disabled.

   **Notes:**
      The default is to use the right-hand side function provided to :c:func:`CVodeInit`  in the internal difference quotient. If the input right-hand side function is  ``NULL``, the default is used.

      This function must be called after the CVLS linear solver interface  has been initialized through a call to :c:func:`CVodeSetLinearSolver`.


When using an iterative linear solver, the user may supply a
preconditioning operator to aid in solution of the system.  This
operator consists of two user-supplied functions, ``psetup`` and
``psolve``, that are supplied to CVODES using the function
:c:func:`CVodeSetPreconditioner`.  The ``psetup`` function supplied to
this routine should handle evaluation and preprocessing of any
Jacobian data needed by the user's preconditioner solve function,
``psolve``.  The user data pointer received through
:c:func:`CVodeSetUserData` (or a pointer to ``NULL`` if user data was not
specified) is passed to the ``psetup`` and ``psolve`` functions.
This allows the user to create an arbitrary structure with relevant
problem data and access it during the execution of the user-supplied
preconditioner functions without using global data in the program.

Also, as described in :numref:`CVODES.Mathematics.ivp_sol`, the CVLS interface
requires that iterative linear solvers stop when the norm of the
preconditioned residual satisfies

.. math::

   \|r\| \le \frac{\epsilon_L \epsilon}{10}

where :math:`\epsilon` is the nonlinear solver tolerance, and the default
:math:`\epsilon_L = 0.05`; this value may be modified by the user through
the :c:func:`CVodeSetEpsLin` function.


.. c:function:: int CVodeSetPreconditioner(void* cvode_mem, CVLsPrecSetupFn psetup, CVLsPrecSolveFn psolve)

   The function ``CVodeSetPreconditioner`` specifies the preconditioner  setup and solve functions.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``psetup`` -- user-defined preconditioner setup function. Pass ``NULL`` if no setup is necessary.
     * ``psolve`` -- user-defined preconditioner solve function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional values have been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_SUNLS_FAIL`` -- An error occurred when setting up preconditioning in the    ``SUNLinearSolver`` object used by the CVLS interface.

   **Notes:**
      The default is ``NULL`` for both arguments (i.e., no  preconditioning).

      This function must be called after the CVLS linear solver  interface has been initialized through a call to  :c:func:`CVodeSetLinearSolver`.

      The function type :c:type:`CVLsPrecSolveFn` is described in :numref:`CVODES.Usage.SIM.user_supplied.psolveFn`.

      The function type :c:type:`CVLsPrecSetupFn` is described in :numref:`CVODES.Usage.SIM.user_supplied.precondFn`.

      The previous routine ``CVSpilsSetPreconditioner`` is now a wrapper  for this routine, and may still be used for backward-compatibility.  However, this will be deprecated in future releases, so we recommend  that users transition to the new routine name soon.


.. c:function:: int CVodeSetEpsLin(void* cvode_mem, realtype eplifac)

   The function ``CVodeSetEpsLin`` specifies the factor by  which the Krylov linear solver's convergence test constant is  reduced from the nonlinear solver test constant.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``eplifac`` -- linear convergence safety factor :math:`(\ge 0)`.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.
     * ``CVLS_ILL_INPUT`` -- The factor ``eplifac`` is negative.

   **Notes:**
      The default value is 0.05.

      This function must be called after the CVLS linear solver  interface has been initialized through a call to  :c:func:`CVodeSetLinearSolver`.

      If ``eplifac`` = 0.0 is passed, the default value is used.

      The previous routine ``CVSpilsSetEpsLin`` is now a wrapper for this  routine, and may still be used for backward-compatibility.  However,  this will be deprecated in future releases, so we recommend that  users transition to the new routine name soon.


.. c:function:: int CVodeSetLSNormFactor(void* cvode_mem, realtype nrmfac)

   The function ``CVodeSetLSNormFactor`` specifies the factor to use when  converting from the integrator tolerance (WRMS norm) to the linear solver  tolerance (L2 norm) for Newton linear system solves e.g.,  ``tol_L2 = fac * tol_WRMS``.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nrmfac`` -- the norm conversion factor. If ``nrmfac`` is:

       - :math:`> 0` then the provided value is used.
       - :math:`= 0` then the conversion factor is computed using the vector length, i.e., ``nrmfac = N_VGetLength(y)`` (*default*).
       - :math:`< 0` then the conversion factor is computed using the vector dot product, i.e., ``nrmfac = N_VDotProd(v,v)`` where all the entries of ``v`` are one.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      This function must be called after the CVLS linear solver  interface has been initialized through a call to :c:func:`CVodeSetLinearSolver`.

      Prior to the introduction of ``N_VGetLength`` in SUNDIALS v5.0.0  (CVODES v5.0.0) the value of ``nrmfac`` was computed using the vector  dot product i.e., the ``nrmfac < 0`` case.


.. _CVODES.Usage.SIM.optional_input.optin_root:

Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following functions can be called to set optional inputs to control
the rootfinding algorithm.

.. c:function:: int CVodeSetRootDirection(void* cvode_mem, int * rootdir)

   The function ``CVodeSetRootDirection`` specifies the direction of  zero-crossings to be located and returned.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``rootdir`` -- state array of length ``nrtfn``, the number of root functions :math:`g_i`, as specified in the call to the function :c:func:`CVodeRootInit`.  A value of :math:`0` for ``rootdir[i]`` indicates that crossing in either direction for :math:`g_i` should be reported.  A value of :math:`+1` or :math:`-1` indicates that the solver should report only zero-crossings where :math:`g_i` is increasing or decreasing, respectively.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_ILL_INPUT`` -- rootfinding has not been activated through a call to :c:func:`CVodeRootInit`.

   **Notes:**
      The default behavior is to monitor for both zero-crossing directions.

.. c:function:: int CVodeSetNoInactiveRootWarn(void* cvode_mem)

   The function ``CVodeSetNoInactiveRootWarn`` disables issuing a warning  if some root function appears to be identically zero at the beginning of the integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      CVODES will not report the initial conditions as a possible zero-crossing  (assuming that one or more components :math:`g_i` are zero at the initial time).  However, if it appears that some :math:`g_i` is identically zero at the initial  time (i.e., :math:`g_i` is zero at the initial time and after the first step),  CVODES will issue a warning which can be disabled with this optional input  function.


.. _CVODES.Usage.SIM.optional_dky:

Interpolated output function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An optional function ``CVodeGetDky`` is available to obtain additional output values.
This function should only be called after a successful return from ``CVode`` as it
provides interpolated values either of :math:`y` or of its derivatives
(up to the current order of the integration method) interpolated to any
value of :math:`t` in the last internal step taken by CVODES.

The call to the function has the following form:

.. c:function:: int CVodeGetDky(void* cvode_mem, realtype t, int k, N_Vector dky)

   The function ``CVodeGetDky`` computes the ``k``-th derivative of the function  ``y`` at time ``t``, i.e. :math:`\dfrac{\mathrm d^{k}y}{\mathrm dt^{k}}(t)`, where :math:`t_n - h_u \leq  t \leq t_n`, :math:`t_n` denotes the current internal time reached, and :math:`h_u` is the  last internal step size successfully used by the solver.  The  user may request ``k`` = :math:`0, 1, \ldots, q_u`, where :math:`q_u` is the current order  (optional output ``qlast``).

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``t`` -- the value of the independent variable at which the derivative is to be evaluated.
     * ``k`` -- the derivative order requested.
     * ``dky`` -- vector containing the derivative. This vector must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- ``CVodeGetDky`` succeeded.
     * ``CV_BAD_K`` -- ``k`` is not in the range :math:`0, 1, \ldots, q_u`.
     * ``CV_BAD_T`` -- ``t`` is not in the interval :math:`[t_n - h_u , t_n]`.
     * ``CV_BAD_DKY`` -- The ``dky`` argument was ``NULL``.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      It is only legal to call the function ``CVodeGetDky`` after a  successful return from :c:func:`CVode`. See :c:func:`CVodeGetCurrentTime`, :c:func:`CVodeGetLastOrder`, and :c:func:`CVodeGetLastStep` in the next section for  access to :math:`t_n`, :math:`q_u`, and :math:`h_u`, respectively.


.. _CVODES.Usage.SIM.optional_output:

Optional output functions
~~~~~~~~~~~~~~~~~~~~~~~~~

CVODES provides an extensive set of functions that can be used to
obtain solver performance information. :numref:`CVODES.Usage.SIM.optional_output.Table`
lists all optional output functions in CVODES, which are then described in detail
in the remainder of this section.

Some of the optional outputs, especially the various counters, can be
very useful in determining how successful the CVODES solver is in
doing its job. For example, the counters ``nsteps`` and ``nfevals`` provide a rough measure of
the overall cost of a given run, and can be compared among runs with
differing input options to suggest which set of options is most
efficient. The ratio ``nniters/nsteps`` measures the performance of the nonlinear solver in
solving the nonlinear systems at each time step; typical values for this
range from 1.1 to 1.8. The ratio ``njevals/nniters`` (in the case of a matrix-based linear
solver), and the ratio ``npevals/nniters`` (in the case of an iterative linear solver)
measure the overall degree of nonlinearity in these systems, and also
the quality of the approximate Jacobian or preconditioner being used.
Thus, for example, ``njevals/nniters`` can indicate if a user-supplied Jacobian is
inaccurate, if this ratio is larger than for the case of the
corresponding internal Jacobian. The ratio ``nliters/nniters`` measures the performance of
the Krylov iterative linear solver, and thus (indirectly) the quality of
the preconditioner.

.. _CVODES.Usage.SIM.optional_output.Table:
.. table:: Optional outputs from CVODES, CVLS, and CVDIAG
   :align: center

   +------------------------------------------------+------------------------------------------+
   |              **Optional output**               |            **Function name**             |
   +================================================+==========================================+
   | **CVODES main solver**                         |                                          |
   +------------------------------------------------+------------------------------------------+
   | Size of CVODES real and integer workspaces     | :c:func:`CVodeGetWorkSpace`              |
   +------------------------------------------------+------------------------------------------+
   | Cumulative number of internal steps            | :c:func:`CVodeGetNumSteps`               |
   +------------------------------------------------+------------------------------------------+
   | No. of calls to r.h.s. function                | :c:func:`CVodeGetNumLinSolvSetups`       |
   +------------------------------------------------+------------------------------------------+
   | No. of calls to linear solver setup function   | :c:func:`CVodeGetNumLinSolvSetups`       |
   +------------------------------------------------+------------------------------------------+
   | No. of local error test failures that have     | :c:func:`CVodeGetNumErrTestFails`        |
   | occurred                                       |                                          |
   +------------------------------------------------+------------------------------------------+
   | Order used during the last step                | :c:func:`CVodeGetLastOrder`              |
   +------------------------------------------------+------------------------------------------+
   | Order to be attempted on the next step         | :c:func:`CVodeGetCurrentOrder`           |
   +------------------------------------------------+------------------------------------------+
   | No. of order reductions due to stability limit | :c:func:`CVodeGetNumStabLimOrderReds`    |
   | detection                                      |                                          |
   +------------------------------------------------+------------------------------------------+
   | Actual initial step size used                  | :c:func:`CVodeGetActualInitStep`         |
   +------------------------------------------------+------------------------------------------+
   | Step size used for the last step               | :c:func:`CVodeGetLastStep`               |
   +------------------------------------------------+------------------------------------------+
   | Step size to be attempted on the next step     | :c:func:`CVodeGetCurrentStep`            |
   +------------------------------------------------+------------------------------------------+
   | Current internal time reached by the solver    | :c:func:`CVodeGetCurrentTime`            |
   +------------------------------------------------+------------------------------------------+
   | Suggested factor for tolerance scaling         | :c:func:`CVodeGetTolScaleFactor`         |
   +------------------------------------------------+------------------------------------------+
   | Error weight vector for state variables        | :c:func:`CVodeGetErrWeights`             |
   +------------------------------------------------+------------------------------------------+
   | Estimated local error vector                   | :c:func:`CVodeGetEstLocalErrors`         |
   +------------------------------------------------+------------------------------------------+
   | No. of nonlinear solver iterations             | :c:func:`CVodeGetNumNonlinSolvIters`     |
   +------------------------------------------------+------------------------------------------+
   | No. of nonlinear convergence failures          | :c:func:`CVodeGetNumNonlinSolvConvFails` |
   +------------------------------------------------+------------------------------------------+
   | All CVODES integrator statistics               | :c:func:`CVodeGetIntegratorStats`        |
   +------------------------------------------------+------------------------------------------+
   | CVODES nonlinear solver statistics             | :c:func:`CVodeGetNonlinSolvStats`        |
   +------------------------------------------------+------------------------------------------+
   | Array showing roots found                      | :c:func:`CVodeGetRootInfo`               |
   +------------------------------------------------+------------------------------------------+
   | No. of calls to user root function             | :c:func:`CVodeGetNumGEvals`              |
   +------------------------------------------------+------------------------------------------+
   | Name of constant associated with a return flag | :c:func:`CVodeGetReturnFlagName`         |
   +------------------------------------------------+------------------------------------------+
   | **CVLS linear solver interface**               |                                          |
   +------------------------------------------------+------------------------------------------+
   | Size of real and integer workspaces            | :c:func:`CVodeGetLinWorkSpace`           |
   +------------------------------------------------+------------------------------------------+
   | No. of Jacobian evaluations                    | :c:func:`CVodeGetNumJacEvals`            |
   +------------------------------------------------+------------------------------------------+
   | No. of r.h.s. calls for finite diff.           | :c:func:`CVodeGetNumLinRhsEvals`         |
   | Jacobian[-vector] evals.                       |                                          |
   +------------------------------------------------+------------------------------------------+
   | No. of linear iterations                       | :c:func:`CVodeGetNumLinIters`            |
   +------------------------------------------------+------------------------------------------+
   | No. of linear convergence failures             | :c:func:`CVodeGetNumLinConvFails`        |
   +------------------------------------------------+------------------------------------------+
   | No. of preconditioner evaluations              | :c:func:`CVodeGetNumPrecEvals`           |
   +------------------------------------------------+------------------------------------------+
   | No. of preconditioner solves                   | :c:func:`CVodeGetNumPrecSolves`          |
   +------------------------------------------------+------------------------------------------+
   | No. of Jacobian-vector setup evaluations       | :c:func:`CVodeGetNumJTSetupEvals`        |
   +------------------------------------------------+------------------------------------------+
   | No. of Jacobian-vector product evaluations     | :c:func:`CVodeGetNumJtimesEvals`         |
   +------------------------------------------------+------------------------------------------+
   | Get all linear solver statistics in one        | :c:func:`CVodeGetLinSolvStats`           |
   | function call                                  |                                          |
   +------------------------------------------------+------------------------------------------+
   | Last return from a linear solver function      | :c:func:`CVodeGetLastLinSolvStats`       |
   +------------------------------------------------+------------------------------------------+
   | Name of constant associated with a return flag | :c:func:`CVodeGetLinReturnFlagName`      |
   +------------------------------------------------+------------------------------------------+
   | **CVDIAG linear solver interface**             |                                          |
   +------------------------------------------------+------------------------------------------+
   | Size of CVDIAG real and integer workspaces     | :c:func:`CVDiagGetWorkSpace`             |
   +------------------------------------------------+------------------------------------------+
   | No. of r.h.s. calls for finite diff. Jacobian  | :c:func:`CVDiagGetNumRhsEvals`           |
   | evals.                                         |                                          |
   +------------------------------------------------+------------------------------------------+
   | Last return from a CVDIAG function             | :c:func:`CVDiagGetLastFlag`              |
   +------------------------------------------------+------------------------------------------+
   | Name of constant associated with a return flag | :c:func:`CVDiagGetReturnFlagName`        |
   +------------------------------------------------+------------------------------------------+

.. _CVODES.Usage.SIM.optional_output.optout_main:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CVODES provides several user-callable functions that can be used to
obtain different quantities that may be of interest to the user, such as
solver workspace requirements, solver performance statistics, as well as
additional data from the CVODES memory block (a suggested tolerance
scaling factor, the error weight vector, and the vector of estimated
local errors). Functions are also provided to extract statistics related
to the performance of the CVODES nonlinear solver used. As a
convenience, additional information extraction functions provide the
optional outputs in groups. These optional output functions are
described next.

.. c:function:: int CVodeGetWorkSpace(void* cvode_mem, long int *lenrw, long int *leniw)

   The function ``CVodeGetWorkSpace`` returns the  CVODES real and integer workspace sizes.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``lenrw`` -- the number of ``realtype`` values in the CVODES workspace.
     * ``leniw`` -- the number of integer values in the CVODES workspace.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      In terms of the problem size :math:`N`, the maximum method order :math:`\texttt{maxord}`, and the number :math:`\texttt{nrtfn}` of root functions (see :numref:`CVODES.Usage.SIM.cvrootinit`) the actual size of the real workspace, in ``realtype`` words, is  given by the following:

      * base value: :math:`\texttt{lenrw} = 96 + ( \texttt{maxord} + 5) N_r + 3\texttt{nrtfn}`;

      * using :c:func:`CVodeSVtolerances`: :math:`\texttt{lenrw} = \texttt{lenrw} + N_r`;

      * with constraint checking (see :c:func:`CVodeSetConstraints`):  :math:`\texttt{lenrw} = \texttt{lenrw} + N_r`;

      where :math:`N_r` is the number of real words in one ``N_Vector`` (:math:`\approx N`).

      The size of the integer workspace (without distinction between ``int``  and ``long int`` words) is given by:

      * base value: :math:`\texttt{leniw} = 40 + ( \texttt{maxord} + 5)N_i + \texttt{nrtfn}`;

      * using :c:func:`CVodeSVtolerances`: :math:`\texttt{leniw} = \texttt{leniw} + N_i`;

      * with constraint checking: :math:`\texttt{lenrw} = \texttt{lenrw} + N_i`;

      where :math:`N_i` is the number of integer words in one ``N_Vector``  (= 1 for ``NVECTOR_SERIAL`` and ``2*npes`` for ``NVECTOR_PARALLEL`` and ``npes`` processors).

      For the default value of :math:`\texttt{maxord}`, no rootfinding, no constraints, and  without using :c:func:`CVodeSVtolerances`, these lengths are given roughly by:

      * For the Adams method: :math:`\texttt{lenrw} = 96 + 17N` and :math:`\texttt{leniw} = 57`

      * For the BDF method: :math:`\texttt{lenrw} = 96 + 10N` and :math:`\texttt{leniw} = 50`

      Note that additional memory is allocated if quadratures and/or forward sensitivity
      integration is enabled. See :numref:`CVODES.Usage.purequad.quad_malloc` and :numref:`CVODES.Usage.FSA.user_callable.sensi_malloc`
      for more details.


.. c:function:: int CVodeGetNumSteps(void* cvode_mem, long int *nsteps)

   The function ``CVodeGetNumSteps`` returns the cumulative number of internal  steps taken by the solver (total so far).

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nsteps`` -- number of steps taken by CVODES.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetNumRhsEvals(void* cvode_mem, long int *nfevals)

   The function ``CVodeGetNumRhsEvals`` returns the  number of calls to the user's right-hand side function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfevals`` -- number of calls to the user's ``f`` function.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      The ``nfevals`` value returned by ``CVodeGetNumRhsEvals`` does not  account for calls made to ``f`` by a linear solver or preconditioner  module.



.. c:function:: int CVodeGetNumLinSolvSetups(void* cvode_mem, long int *nlinsetups)

   The function ``CVodeGetNumLinSolvSetups`` returns the  number of calls made to the linear solver's setup function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nlinsetups`` -- number of calls made to the linear solver setup function.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetNumErrTestFails(void* cvode_mem, long int *netfails)

   The function ``CVodeGetNumErrTestFails`` returns the  number of local error test failures that have occurred.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``netfails`` -- number of error test failures.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetLastOrder(void* cvode_mem, int *qlast)

   The function ``CVodeGetLastOrder`` returns the  integration method order used during the last internal step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``qlast`` -- method order used on the last internal step.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetCurrentOrder(void* cvode_mem, int *qcur)

   The function ``CVodeGetCurrentOrder`` returns the  integration method order to be used on the next internal step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``qcur`` -- method order to be used on the next internal step.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetLastStep(void* cvode_mem, realtype *hlast)

   The function ``CVodeGetLastStep`` returns the  integration step size taken on the last internal step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``hlast`` -- step size taken on the last internal step.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetCurrentStep(void* cvode_mem, realtype *hcur)

   The function ``CVodeGetCurrentStep`` returns the  integration step size to be attempted on the next internal step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``hcur`` -- step size to be attempted on the next internal step.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetActualInitStep(void* cvode_mem, realtype *hinused)

   The function ``CVodeGetActualInitStep`` returns the  value of the integration step size used on the first step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``hinused`` -- actual value of initial step size.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      Even if the value of the initial integration step size was specified  by the user through a call to :c:func:`CVodeSetInitStep`, this value might have  been changed by CVODES to ensure that the step size is within the  prescribed bounds (:math:`h_{min} \leq h_0 \leq h_{max}`), or to satisfy the local error test condition.



.. c:function:: int CVodeGetCurrentTime(void* cvode_mem, realtype *tcur)

   The function ``CVodeGetCurrentTime`` returns the  current internal time reached by the solver.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``tcur`` -- current internal time reached.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetNumStabLimOrderReds(void* cvode_mem, long int *nslred)

   The function ``CVodeGetNumStabLimOrderReds`` returns the  number of order reductions dictated by the BDF stability limit sdetection algorithm (see :numref:`CVODES.Mathematics.stablimit`).

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nslred`` -- number of order reductions due to stability limit detection.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      If the stability limit detection algorithm was not initialized  (:c:func:`CVodeSetStabLimDet` was not called), then ``nslred`` = 0.



.. c:function:: int CVodeGetTolScaleFactor(void* cvode_mem, realtype *tolsfac)

   The function ``CVodeGetTolScaleFactor`` returns a  suggested factor by which the user's tolerances  should be scaled when too much accuracy has been  requested for some internal step.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``tolsfac`` -- suggested scaling factor for user-supplied tolerances.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetErrWeights(void* cvode_mem, N_Vector eweight)

   The function ``CVodeGetErrWeights`` returns the solution error weights at  the current time. These are the reciprocals of the :math:`W_i` given by :eq:`CVODES_errwt`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``eweight`` -- solution error weights at the current time.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   .. warning::

      The user must allocate memory for ``eweight``.



.. c:function:: int CVodeGetEstLocalErrors(void* cvode_mem, N_Vector ele)

   The function ``CVodeGetEstLocalErrors`` returns the  vector of estimated local errors.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``ele`` -- estimated local errors.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   .. warning::

      The user must allocate memory for ``ele``.

      The values returned in ``ele`` are valid only if :c:func:`CVode` returned  a non-negative value.

      The ``ele`` vector, togther with the ``eweight`` vector from :c:func:`CVodeGetErrWeights`, can be used to determine how the various  components of the system contributed to the estimated local error  test.  Specifically, that error test uses the RMS norm of a vector  whose components are the products of the components of these two vectors.  Thus, for example, if there were recent error test failures, the components  causing the failures are those with largest values for the products,  denoted loosely as ``eweight[i]*ele[i]``.



.. c:function:: int CVodeGetIntegratorStats(void* cvode_mem, long int *nsteps, long int *nfevals, long int *nlinsetups, long int *netfails, int *qlast, int *qcur, realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur)

   The function ``CVodeGetIntegratorStats`` returns the CVODES integrator statistics as a group.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nsteps`` -- number of steps taken by CVODES.
     * ``nfevals`` -- number of calls to the user's ``f`` function.
     * ``nlinsetups`` -- number of calls made to the linear solver setup function.
     * ``netfails`` -- number of error test failures.
     * ``qlast`` -- method order used on the last internal step.
     * ``qcur`` -- method order to be used on the next internal step.
     * ``hinused`` -- actual value of initial step size.
     * ``hlast`` -- step size taken on the last internal step.
     * ``hcur`` -- step size to be attempted on the next internal step.
     * ``tcur`` -- current internal time reached.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetNumNonlinSolvIters(void* cvode_mem, long int *nniters)

   The function ``CVodeGetNumNonlinSolvIters`` returns the  number of nonlinear iterations performed.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nniters`` -- number of nonlinear iterations performed.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_MEM_FAIL`` -- The ``SUNNonlinearSolver`` module is ``NULL``.



.. c:function:: int CVodeGetNumNonlinSolvConvFails(void* cvode_mem, long int *nncfails)

   The function ``CVodeGetNumNonlinSolvConvFails`` returns the  number of nonlinear convergence failures that have occurred.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nncfails`` -- number of nonlinear convergence failures.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.



.. c:function:: int CVodeGetNonlinSolvStats(void* cvode_mem, long int *nniters, long int  *nncfails)

   The function ``CVodeGetNonlinSolvStats`` returns the  CVODES nonlinear solver statistics as a group.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nniters`` -- number of nonlinear iterations performed.
     * ``nncfails`` -- number of nonlinear convergence failures.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_MEM_FAIL`` -- The ``SUNNonlinearSolver`` module is ``NULL``.



.. c:function:: char* CVodeGetReturnFlagName(int flag)

   The function ``CVodeGetReturnFlagName`` returns the  name of the CVODES constant corresponding to ``flag``.

   **Arguments:**
     * ``flag`` -- return flag from a CVODES function.

   **Return value:**
     * A string containing the name of the corresponding constant



.. _CVODES.Usage.SIM.optional_output.optout_root:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two optional output functions associated with rootfinding.


.. c:function:: int CVodeGetRootInfo(void* cvode_mem, int * rootsfound)

   The function ``CVodeGetRootInfo`` returns an array showing which  functions were found to have a root.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``rootsfound`` -- array of length ``nrtfn`` with the indices of the user functions :math:`g_i` found to have a root.  For :math:`i=0,\ldots,\texttt{nrtfn}-1`, ``rootsfound[i]`` :math:`\ne 0` if :math:`g_i` has a root, and ``rootsfound[i]`` :math:`= 0` if not.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.

   **Notes:**
      Note that, for the components :math:`g_i` for which a root was found, the sign of ``rootsfound[i]`` indicates the direction of  zero-crossing. A value of +1 indicates that :math:`g_i` is increasing,  while a value of -1 indicates a decreasing :math:`g_i`.

   .. warning::

      The user must allocate memory for the vector ``rootsfound``.


.. c:function:: int CVodeGetNumGEvals(void* cvode_mem, long int *ngevals)

   The function ``CVodeGetNumGEvals`` returns the cumulative  number of calls made to the user-supplied root function :math:`g`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``ngevals`` -- number of calls made to the user's function :math:`g` thus far.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.


.. _CVODES.Usage.SIM.optional_output.optout_ls:

CVLS linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the CVLS modules:
workspace requirements, number of calls to the Jacobian routine, number
of calls to the right-hand side routine for finite-difference Jacobian
or Jacobian-vector product approximation, number of linear iterations,
number of linear convergence failures, number of calls to the
preconditioner setup and solve routines, number of calls to the
Jacobian-vector setup and product routines, and last return value from a
linear solver function. Note that, where the name of an output would
otherwise conflict with the name of an optional output from the main
solver, a suffix (for Linear Solver) has been added (e.g. ``lenrwLS``).


.. c:function:: int CVodeGetLinWorkSpace(void* cvode_mem, long int *lenrwLS, long int *leniwLS)

   The function ``CVodeGetLinWorkSpace`` returns the sizes of the real and  integer workspaces used by the CVLS linear solver interface.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``lenrwLS`` -- the number of ``realtype`` values in the CVLS workspace.
     * ``leniwLS`` -- the number of integer values in the CVLS workspace.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output values have been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.

   **Notes:**
      The workspace requirements reported by this routine correspond only  to memory allocated within this interface and to memory allocated by  the ``SUNLinearSolver`` object attached to it.  The template Jacobian  matrix allocated by the user outside of CVLS is not included in  this report.

      The previous routines ``CVDlsGetWorkspace`` and  ``CVSpilsGetWorkspace`` are now wrappers for this routine, and may  still be used for backward-compatibility.  However, these will be  deprecated in future releases, so we recommend that users transition  to the new routine name soon.


.. c:function:: int CVodeGetNumJacEvals(void* cvode_mem, long int *njevals)

   The function ``CVodeGetNumJacEvals`` returns the  number of calls made to the CVLS Jacobian approximation  function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``njevals`` -- the number of calls to the Jacobian function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.

   **Notes:**
      The previous routine ``CVDlsGetNumJacEvals`` is now a wrapper for  this routine, and may still be used for backward-compatibility.  However, this will be deprecated in future releases, so we recommend  that users transition to the new routine name soon.


.. c:function:: int CVodeGetNumLinRhsEvals(void* cvode_mem, long int *nfevalsLS)

   The function ``CVodeGetNumLinRhsEvals`` returns the  number of calls made to the user-supplied right-hand side function  due to the finite difference Jacobian approximation or finite  difference Jacobian-vector product approximation.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfevalsLS`` -- the number of calls made to the user-supplied right-hand side function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.

   **Notes:**
      The value ``nfevalsLS`` is incremented only if one of the default  internal difference quotient functions is used.

      The previous routines ``CVDlsGetNumRhsEvals`` and  ``CVSpilsGetNumRhsEvals`` are now wrappers for this routine, and may  still be used for backward-compatibility.  However, these will be  deprecated in future releases, so we recommend that users transition  to the new routine name soon.


.. c:function:: int CVodeGetNumLinIters(void* cvode_mem, long int *nliters)

   The function ``CVodeGetNumLinIters`` returns the  cumulative number of linear iterations.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nliters`` -- the current number of linear iterations.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.

   **Notes:**
      The previous routine ``CVSpilsGetNumLinIters`` is now a wrapper for  this routine, and may still be used for backward-compatibility.  However, this will be deprecated in future releases, so we recommend  that users transition to the new routine name soon.


.. c:function:: int CVodeGetNumLinConvFails(void* cvode_mem, long int *nlcfails)

   The function ``CVodeGetNumLinConvFails`` returns the  cumulative number of linear convergence failures.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nlcfails`` -- the current number of linear convergence failures.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.

   **Notes:**
      The previous routine ``CVSpilsGetNumConvFails`` is now a wrapper for  this routine, and may still be used for backward-compatibility.  However, this will be deprecated in future releases, so we recommend  that users transition to the new routine name soon.


.. c:function:: int CVodeGetNumPrecEvals(void* cvode_mem, long int *npevals)

   The function ``CVodeGetNumPrecEvals`` returns the  number of preconditioner evaluations, i.e., the number of  calls made to ``psetup`` with ``jok = SUNFALSE``.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``npevals`` -- the current number of calls to ``psetup``.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.

   **Notes:**
      The previous routine ``CVSpilsGetNumPrecEvals`` is now a wrapper for  this routine, and may still be used for backward-compatibility.  However, this will be deprecated in future releases, so we recommend  that users transition to the new routine name soon.


.. c:function:: int CVodeGetNumPrecSolves(void* cvode_mem, long int *npsolves)

   The function ``CVodeGetNumPrecSolves`` returns the  cumulative number of calls made to the preconditioner  solve function, ``psolve``.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``npsolves`` -- the current number of calls to ``psolve``.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.


.. c:function:: int CVodeGetNumJTSetupEvals(void* cvode_mem, long int *njtsetup)

   The function ``CVodeGetNumJTSetupEvals`` returns the  cumulative number of calls made to the Jacobian-vector setup  function ``jtsetup``.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``njtsetup`` -- the current number of calls to ``jtsetup``.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.


.. c:function:: int CVodeGetNumJtimesEvals(void* cvode_mem, long int *njvevals)

   The function ``CVodeGetNumJtimesEvals`` returns the  cumulative number of calls made to the Jacobian-vector function  ``jtimes``.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``njvevals`` -- the current number of calls to ``jtimes``.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.


.. c:function:: int CVodeGetLinSolvStats(void* cvode_mem, long int* njevals, long int* nfevalsLS, long int* nliters, long int* nlcfails, long int* npevals, long int* npsolves, long int* njtsetups, long int* njtimes)

   The function ``CVodeGetLinSolvStats`` returns CVODES linear solver  statistics.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``njevals`` -- the current number of calls to the Jacobian function.
     * ``nfevalsLS`` -- the current number of calls made to the user-supplied right-hand side function by the linear solver.
     * ``nliters`` -- the current number of linear iterations.
     * ``nlcfails`` -- the current number of linear convergence failures.
     * ``npevals`` -- the current number of calls to ``psetup``.
     * ``npsolves`` -- the current number of calls to ``psolve``.
     * ``njtsetup`` -- the current number of calls to ``jtsetup``.
     * ``njtimes`` -- the current number of calls to ``jtimes``.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.


.. c:function:: int CVodeGetLastLinFlag(void* cvode_mem, long int *lsflag)

   The function ``CVodeGetLastLinFlag`` returns the  last return value from a CVLS routine.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``lsflag`` -- the value of the last return flag from a CVLS function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_LMEM_NULL`` -- The CVLS linear solver has not been initialized.

   **Notes:**
      If the CVLS setup function failed (i.e., :c:func:`CVode` returned  ``CV_LSETUP_FAIL``) when using the ``SUNLINSOL_DENSE`` or  ``SUNLINSOL_BAND`` modules, then the value of ``lsflag`` is equal to  the column index (numbered from one) at which a zero diagonal  element was encountered during the LU factorization of the (dense or  banded) Jacobian matrix.

      If the CVLS setup function failed when using another ``SUNLinearSolver``  module, then ``lsflag`` will be ``SUNLS_PSET_FAIL_UNREC``,  ``SUNLS_ASET_FAIL_UNREC``, or  ``SUNLS_PACKAGE_FAIL_UNREC``.

      If the CVLS solve function failed (i.e., :c:func:`CVode` returned  ``CV_LSOLVE_FAIL``), then ``lsflag`` contains the error return  flag from the ``SUNLinearSolver`` object, which will be one of: ``SUNLS_MEM_NULL``, indicating that the ``SUNLinearSolver`` memory is ``NULL``;   ``SUNLS_ATIMES_FAIL_UNREC``, indicating an unrecoverable failure in the  Jv function; ``SUNLS_PSOLVE_FAIL_UNREC``, indicating that the preconditioner solve  function ``psolve`` failed unrecoverably;  ``SUNLS_GS_FAIL``, indicating a failure in the Gram-Schmidt  procedure (SPGMR and SPFGMR only);  ``SUNLS_QRSOL_FAIL``, indicating that the matrix R was found to be  singular during the QR solve phase (SPGMR and SPFGMR only); or  ``SUNLS_PACKAGE_FAIL_UNREC``, indicating an unrecoverable  failure in an external iterative linear solver package.

      The previous routines ``CVDlsGetLastFlag`` and  ``CVSpilsGetLastFlag`` are now wrappers for this routine, and may  still be used for backward-compatibility.  However, these will be  deprecated in future releases, so we recommend that users transition  to the new routine name soon.


.. c:function:: int CVodeGetLinReturnFlagName(long int lsflag)

   The function ``CVodeGetLinReturnFlagName`` returns the name of the CVLS constant corresponding to ``lsflag``.

   **Arguments:**
     * ``lsflag`` -- a return flag from a ``CVLS`` function.

   **Return value:**
     * The return value is a string containing the name of the corresponding constant. If :math:`1 \leq \text{lsflag} \leq N` (LU factorization failed), this routine returns "NONE".

   **Notes:**
      The previous routines ``CVDlsGetReturnFlagName`` and  ``CVSpilsGetReturnFlagName`` are now wrappers for this routine, and may  still be used for backward-compatibility.  However, these will be  deprecated in future releases, so we recommend that users transition  to the new routine name soon.


.. _CVODES.Usage.SIM.optional_output.optout_diag:

Diagonal linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the CVDIAG module:
workspace requirements, number of calls to the right-hand side routine
for finite-difference Jacobian approximation, and last return value from
a CVDIAG function. Note that, where the name of an output would
otherwise conflict with the name of an optional output from the main
solver, a suffix (for Linear Solver) has been added here (e.g. ``lenrwLS``).


.. c:function:: int CVDiagGetWorkSpace(void* cvode_mem, long int *lenrwLS, long int *leniwLS)

   The function ``CVDiagGetWorkSpace`` returns the  CVDIAG real and integer workspace sizes.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``lenrwLS`` -- the number of ``realtype`` values in the CVDIAG workspace.
     * ``leniwLS`` -- the number of integer values in the CVDIAG workspace.

   **Return value:**
     * ``CVDIAG_SUCCESS`` -- The optional output valus have been successfully set.
     * ``CVDIAG_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVDIAG_LMEM_NULL`` -- The CVDIAG linear solver has not been initialized.

   **Notes:**
      In terms of the problem size :math:`N`, the actual size of the real workspace  is roughly :math:`3 N` ``realtype`` words.


.. c:function:: int CVDiagGetNumRhsEvals(void* cvode_mem, long int *nfevalsLS)

   The function ``CVDiagGetNumRhsEvals`` returns the  number of calls made to the user-supplied right-hand side function due to the  finite difference Jacobian approximation.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfevalsLS`` -- the number of calls made to the user-supplied right-hand side function.

   **Return value:**
     * ``CVDIAG_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVDIAG_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVDIAG_LMEM_NULL`` -- The CVDIAG linear solver has not been initialized.

   **Notes:**
      The number of diagonal approximate Jacobians formed is  equal to the number of calls made to the linear solver setup function  (see :c:func:`CVodeGetNumLinSolvSetups`).


.. c:function:: int CVDiagGetLastFlag(void* cvode_mem, long int *lsflag)

   The function ``CVDiagGetLastFlag`` returns the  last return value from a CVDIAG routine.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``lsflag`` -- the value of the last return flag from a CVDIAG function.

   **Return value:**
     * ``CVDIAG_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVDIAG_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CVDIAG_LMEM_NULL`` -- The CVDIAG linear solver has not been initialized.

   **Notes:**
      If the CVDIAG setup function failed (:c:func:`CVode` returned ``CV_LSETUP_FAIL``),  the value of ``lsflag`` is equal to ``CVDIAG_INV_FAIL``, indicating that a  diagonal element with value zero was encountered.  The same value is also returned if the CVDIAG solve function failed  (:c:func:`CVode` returned ``CV_LSOLVE_FAIL``).


.. c:function:: char* CVDiagGetReturnFlagName(long int lsflag)

   The function ``CVDiagGetReturnFlagName`` returns the  name of the CVDIAG constant corresponding to ``lsflag``.

   **Arguments:**
     * ``lsflag`` -- a return flag from a ``CVDIAG`` function.

   **Return value:**
     * A string containing the name of the corresponding constant.


.. _CVODES.Usage.SIM.reinit:

CVODES reinitialization function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function :c:func:`CVodeReInit` reinitializes the main CVODES solver for the solution of
a new problem, where a prior call to :c:func:`CVodeInit` has been made. The new problem must
have the same size as the previous one. :c:func:`CVodeReInit` performs the same input checking
and initializations that does, but does no memory allocation, as it
assumes that the existing internal memory is sufficient for the new
problem. A call to :c:func:`CVodeReInit` deletes the solution history that was stored
internally during the previous integration. Following a successful call
to :c:func:`CVodeReInit`, call :c:func:`CVode` again for the solution of the new problem.

The use of :c:func:`CVodeReInit` requires that the maximum method order, denoted by ``maxord``, be no
larger for the new problem than for the previous problem. This condition
is automatically fulfilled if the multistep method parameter ``lmm`` is
unchanged (or changed from ``CV_ADAMS`` to ``CV_BDF``) and the default value for ``maxord`` is specified.

If there are changes to the linear solver specifications, make the
appropriate calls to either the linear solver objects themselves, or to
the CVLS interface routines, as described in
:numref:`CVODES.Usage.SIM.user_callable.lin_solv_init`. Otherwise, all solver inputs set
previously remain in effect.

One important use of the :c:func:`CVodeReInit` function is in the treating of jump
discontinuities in the RHS function. Except in cases of fairly small
jumps, it is usually more efficient to stop at each point of
discontinuity and restart the integrator with a readjusted ODE model,
using a call to :c:func:`CVodeReInit`. To stop when the location of the discontinuity is
known, simply make that location a value of tout. To stop when the
location of the discontinuity is determined by the solution, use the
rootfinding feature. In either case, it is critical that the RHS
function *not* incorporate the discontinuity, but rather have a smooth
extention over the discontinuity, so that the step across it (and
subsequent rootfinding, if used) can be done efficiently. Then use a
switch within the RHS function (communicated through ``user_data``) that can be
flipped between the stopping of the integration and the restart, so that
the restarted problem uses the new values (which have jumped). Similar
comments apply if there is to be a jump in the dependent variable
vector.


.. c:function:: int CVodeReInit(void* cvode_mem, realtype t0, N_Vector y0)

   The function ``CVodeReInit`` provides required problem specifications  and reinitializes CVODES.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``t0`` -- is the initial value of :math:`t`.
     * ``y0`` -- is the initial value of :math:`y`.

   **Return value:**
     * ``CV_SUCCESS`` -- The call was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory block was not initialized through a previous call to :c:func:`CVodeCreate`.
     * ``CV_NO_MALLOC`` -- Memory space for the CVODES memory block was not allocated through a previous call to :c:func:`CVodeInit`.
     * ``CV_ILL_INPUT`` -- An input argument was an illegal value.

   **Notes:**
      If an error occurred, ``CVodeReInit`` also sends an error message to the  error handler function.


.. _CVODES.Usage.SIM.user_supplied:

User-supplied functions
-----------------------

The user-supplied functions consist of one function defining the ODE,
(optionally) a function that handles error and warning messages,
(optionally) a function that provides the error weight vector,
(optionally) one or two functions that provide Jacobian-related
information for the linear solver, and (optionally) one or two functions
that define the preconditioner for use in any of the Krylov iterative
algorithms.

.. _CVODES.Usage.SIM.user_supplied.rhsFn:

ODE right-hand side
~~~~~~~~~~~~~~~~~~~

The user must provide a function of type defined as follows:

.. c:type:: int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *user_data);

   This function computes the ODE right-hand side for a given value of the independent variable :math:`t` and state vector :math:`y`.

   **Arguments:**
      * ``t`` -- is the current value of the independent variable.
      * ``y`` -- is the current value of the dependent variable vector, :math:`y(t)`.
      * ``ydot`` -- is the output vector :math:`f(t,y)`.
      * ``user_data`` -- is the ``user_data`` pointer passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      A ``CVRhsFn`` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct), or a negative
      value if it failed unrecoverably (in which case the integration is halted and
      ``CV_RHSFUNC_FAIL`` is returned).

   **Notes:**
      Allocation of memory for ``ydot`` is handled within CVODES.

      A recoverable failure error return from the ``CVRhsFn`` is typically used to
      flag a value of the dependent variable :math:`y` that is "illegal" in
      some way (e.g., negative where only a non-negative value is physically
      meaningful). If such a return is made, CVODES will attempt to recover
      (possibly repeating the nonlinear solve, or reducing the step size)
      in order to avoid this recoverable error return.

      For efficiency reasons, the right-hand side function is not evaluated
      at the converged solution of the nonlinear solver. Therefore, in general, a
      recoverable error in that converged value cannot be corrected.  (It may be
      detected when the right-hand side function is called the first time during
      the following integration step, but a successful step cannot be undone.)
      However, if the user program also includes quadrature integration, the
      state variables can be checked for legality in the call to
      :c:type:`CVQuadRhsFn`, which is called at the converged solution of the
      nonlinear system, and therefore CVODES can be flagged to attempt
      to recover from such a situation. Also, if sensitivity analysis is
      performed with one of the staggered methods, the ODE right-hand side
      function is called at the converged solution of the nonlinear system,
      and a recoverable error at that point can be flagged, and CVODES
      will then try to correct it.

      There are two other situations in which recovery is not possible
      even if the right-hand side function returns a recoverable error flag.
      One is when this occurs at the very first call to the ``CVRhsFn``
      (in which case CVODES returns ``CV_FIRST_RHSFUNC_ERR``).
      The other is when a recoverable error is reported by ``CVRhsFn``
      after an error test failure, while the linear multistep method order is
      equal to 1 (in which case CVODES returns ``CV_UNREC_RHSFUNC_ERR``).


.. _CVODES.Usage.SIM.user_supplied.ehFn:

Error message handler function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an alternative to the default behavior of directing error and warning
messages to the file pointed to by ``errfp`` (see :c:func:`CVodeSetErrFile`), the user may provide a
function of type ``CVErrHandlerFn`` to process any such messages. The function type
:c:type:`CVErrHandlerFn` is defined as follows:

.. c:type:: void (*CVErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *eh_data);

   This function processes error and warning message from CVODES and it sub-modules.

   **Arguments:**
      * ``error_code`` is the error code.
      * ``module`` is the name of the CVODES module reporting the error.
      * ``function`` is the name of the function in which the error occurred.
      * ``msg`` is the error message.
      * ``eh_data`` is a pointer to user data, the same as the ``eh_data`` parameter passed to :c:func:`CVodeSetErrHandlerFn`.

   **Return value:**
      * void

   **Notes:**
      ``error_code`` is negative for errors and positive (``CV_WARNING``) for warnings.
      If a function that returns a pointer to memory encounters an error, it sets ``error_code`` to 0.


.. _CVODES.Usage.SIM.user_supplied.monitorfn:

Monitor function
~~~~~~~~~~~~~~~~

A user may provide a function of type ``CVMonitorFn`` to monitor the integrator progress
throughout a simulation. For example, a user may want to check
integrator statistics as a simulation progresses.

.. c:type:: void (*CVMonitorFn)(void* cvode_mem, void* user_data);

   This function is used to monitor the CVODES integrator throughout a simulation.

   **Arguments:**
      * ``cvode_mem`` -- the CVODES memory pointer.
      * ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      Should return 0 if successful, or a negative value if unsuccessful.

   .. warning::

      This function should only be utilized for monitoring the integrator progress (i.e., for debugging).


.. _CVODES.Usage.SIM.user_supplied.ewtsetFn:

Error weight function
~~~~~~~~~~~~~~~~~~~~~

As an alternative to providing the relative and absolute tolerances, the
user may provide a function of type ``CVEwtFn`` to compute a vector containing the
weights in the WRMS norm

.. math::

   \|\ v \|_{\mbox{WRMS}} = \sqrt{\frac1N \sum_{i=1}^N (W_i \cdot v_i)^2}.

These weights will be used in place of those defined by Eq.
:eq:`CVODES_errwt`. The function type is defined as follows:

.. c:type:: int (*CVEwtFn)(N_Vector y, N_Vector ewt, void *user_data);

   This function computes the WRMS error weights for the vector :math:`y`.

   **Arguments:**
      * ``y`` -- the value of the dependent variable vector at which the weight vector is to be computed.
      * ``ewt`` -- the output vector containing the error weights.
      * ``user_data`` a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      Should return 0 if successful, or -1 if unsuccessful.

   **Notes:**
      Allocation of memory for ``ewt`` is handled within CVODES.

   .. warning::

      The error weight vector must have all components positive. It is the user's responsiblity to perform this test and return -1 if it is not satisfied.


.. _CVODES.Usage.SIM.user_supplied.rootFn:

Rootfinding function
~~~~~~~~~~~~~~~~~~~~

If a rootfinding problem is to be solved during the integration of the
ODE system, the user must supply a C function of type ``CVRootFn``, defined as
follows:

.. c:type:: int (*CVRootFn)(realtype t, N_Vector y, realtype *gout, void *user_data);

   This function implements a vector-valued function :math:`g(t,y)` such that the roots of the ``nrtfn`` components :math:`g_i(t,y)` are sought.

   **Arguments:**
      * ``t`` -- the current value of the independent variable.
      * ``y`` -- the current value of the dependent variable vector, :math:`y(t)`.
      * ``gout`` -- the output array of length ``nrtfn`` with components :math:`g_i(t,y)`.
      * ``user_data`` a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      A ``CVRootFn`` should return 0 if successful or a non-zero value if an error occured (in which case the integration is haled and ``CVode`` returns ``CV_RTFUNC_FAIL``.

   **Notes:**
      Allocation of memory for ``gout`` is automatically handled within CVODES.


.. _CVODES.Usage.SIM.user_supplied.jacFn:

Jacobian construction (matrix-based linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a matrix-based linear solver module is used (i.e., a non-``NULL``
``SUNMatrix`` object was supplied to :c:func:`CVodeSetLinearSolver`), the user may optionally provide
a function of type ``CVLsJacFn`` for evaluating the Jacobian of the ODE right-hand
side function (or an approximation of it). ``CVLsJacFn`` is defined as follows:


.. c:type:: int (*CVLsJacFn)(realtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

   This function computes the Jacobian matrix :math:`J = \dfrac{\partial f}{\partial y}` (or an approximation to it).

   **Arguments:**
      * ``t`` -- the current value of the independent variable.
      * ``y`` -- the current value of the dependent variable vector, namely the predicted value of :math:`y(t)`.
      * ``fy`` -- the current value of the vector :math:`f(t,y)`.
      * ``Jac`` -- the output Jacobian matrix.
      * ``user_data`` a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.
      * ``tmp1, tmp2, tmp3`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by a ``CVLsJacFn`` function as temporary storage or work space.

   **Return value:**
      Should return 0 if successful, a positive value if a recoverable error occurred (in which case CVODES will attempt to correct, while CVLS sets ``last_flag`` to ``CVLS_JACFUNC_RECVR``), or a negative value if it failed unrecoverably (in which case the integration is halted, :c:func:`CVode` returns ``CV_LSETUP_FAIL`` and CVLS sets ``last_flag`` to ``CVLS_JACFUNC_UNRECVR``).

   **Notes:**
      Information regarding the structure of the specific ``SUNMatrix``
      structure (e.g. number of rows, upper/lower bandwidth, sparsity type)
      may be obtained through using the implementation-specific ``SUNMatrix``
      interface functions (see :numref:`SUNMatrix` for
      details).

      With direct linear solvers (i.e., linear solvers with type ``SUNLINEARSOLVER_DIRECT``), the
      Jacobian matrix :math:`J(t,y)` is zeroed out prior to calling the
      user-supplied Jacobian function so only nonzero elements need to be
      loaded into ``Jac``.

      With the default nonlinear solver (the native SUNDIALS Newton
      method), each call to the user’s ``CVLsJacFn`` function is preceded by a call to the
      ``CVRhsFn`` user function with the same ``(t,y)`` arguments. Thus, the Jacobian function can
      use any auxiliary data that is computed and saved during the evaluation
      of the ODE right-hand side. In the case of a user-supplied or external
      nonlinear solver, this is also true if the nonlinear system function is
      evaluated prior to calling the :ref:`linear solver setup function <CVODES.Usage.SIM.user_supplied>`.

      If the user’s ``CVLsJacFn`` function uses difference quotient approximations, then it
      may need to access quantities not in the argument list. These include
      the current step size, the error weights, etc. To obtain these, the user
      will need to add a pointer to ``cv_mem`` in ``user_data`` and then use the
      ``CVodeGet*`` functions described in :numref:`CVODES.Usage.SIM.optional_output`.
      The unit roundoff can be accessed as ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

      **Dense**:
      A user-supplied dense Jacobian function must load the :math:`N` by :math:`N`
      dense matrix ``Jac`` with an approximation to the Jacobian matrix :math:`J(t,y)`
      at the point :math:`(t, y)`.  The accessor macros ``SM_ELEMENT_D``
      and ``SM_COLUMN_D`` allow the user to read and write dense matrix
      elements without making explicit references to the underlying
      representation of the SUNMATRIX_DENSE type.
      ``SM_ELEMENT_D(J, i, j)`` references the :math:`(i, j\text{-th})`
      element of the dense matrix ``Jac`` (with :math:`i, j = 0\ldots N-1`).
      This macro is meant for small problems for which efficiency
      of access is not a major concern.  Thus, in terms of the indices :math:`m`
      and :math:`n` ranging from 1 to :math:`N`, the Jacobian element :math:`J_{m,n}` can
      be set using the statement ``SM_ELEMENT_D(J, m-1, n-1) =``
      :math:`J_{m,n}`.  Alternatively, ``SM_COLUMN_D(J, j)`` returns a
      pointer to the first element of the :math:`j`-th column of ``Jac``
      (with :math:`j = 0\ldots N-1`), and the elements of the :math:`j`-th column
      can then be accessed using ordinary array indexing.  Consequently,
      :math:`J(m,n)` can be loaded using the statements
      ``col_n = SM_COLUMN_D(J, n-1); col_n[m-1] =`` :math:`J(m,n)`.
      For large problems, it is more efficient to use ``SM_COLUMN_D``
      than to use ``SM_ELEMENT_D``.  Note that both of these macros
      number rows and columns starting from 0.  The SUNMATRIX_DENSE type
      and accessor macros are documented in :numref:`SUNMatrix.Dense`.

      **Banded**:
      A user-supplied banded Jacobian function must load the :math:`N` by :math:`N` banded matrix
      ``Jac`` with the elements of the Jacobian :math:`J(t,y)` at the point
      :math:`(t,y)`.  The accessor macros ``SM_ELEMENT_B``,
      ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B`` allow the user
      to read and write band matrix elements without making specific
      references to the underlying representation of the SUNMATRIX_BAND
      type.  ``SM_ELEMENT_B(J, i, j)`` references the :math:`(i,j)`,
      element of the band matrix ``Jac``, counting from 0.
      This macro is meant for use in small problems for which efficiency
      of access is not a major concern.  Thus, in terms of the indices :math:`m`
      and :math:`n` ranging from 1 to :math:`N` with :math:`(m,n)` within the band defined
      by ``mupper`` and ``mlower``, the Jacobian element :math:`J(m,n)` can
      be loaded using the statement ``SM_ELEMENT_B(J, m-1, n-1) =``
      :math:`J(m,n)`. The elements within the band are those with ``-mupper``
      :math:`\le m-n \le` ``mlower``. Alternatively,
      ``SM_COLUMN_B(J, j)`` returns a pointer to the diagonal element
      of the :math:`j`-th column of ``Jac``, and if we assign this address
      to ``realtype *col_j``, then the :math:`i`-th element of the
      :math:`j`-th column is given by
      ``SM_COLUMN_ELEMENT_B(col_j, i, j)``, counting from 0.  Thus,
      for :math:`(m,n)` within the band, :math:`J(m,n)` can be loaded by setting
      ``col_n = SM_COLUMN_B(J, n-1); SM_COLUMN_ELEMENT_B(col_n, m-1, n-1) =`` :math:`J(m,n)`. The
      elements of the :math:`j`-th column can also be accessed via ordinary
      array indexing, but this approach requires knowledge of the
      underlying storage for a band matrix of type SUNMATRIX_BAND.
      The array ``col_n`` can be indexed from ``-mupper`` to
      ``mlower``. For large problems, it is more efficient to use
      ``SM_COLUMN_B`` and ``SM_COLUMN_ELEMENT_B`` than to use the
      ``SM_ELEMENT_B`` macro.  As in the dense case, these macros all
      number rows and columns starting from 0.  The SUNMATRIX_BAND type
      and accessor macros are documented in :numref:`SUNMatrix.Band`.

      **Sparse**:
      A user-supplied sparse Jacobian function must load the :math:`N` by :math:`N`
      compressed-sparse-column or compressed-sparse-row matrix ``Jac``
      with an approximation to the Jacobian matrix :math:`J(t,y)` at the point
      :math:`(t,y)`.  Storage for ``Jac`` already exists on entry to
      this function, although the user should ensure that sufficient space
      is allocated in ``Jac`` to hold the nonzero values to be set; if
      the existing space is insufficient the user may reallocate the data
      and index arrays as needed.  The amount of allocated space in a
      SUNMATRIX_SPARSE object may be accessed using the macro
      ``SM_NNZ_S`` or the routine ``SUNSparseMatrix_NNZ``.  The
      SUNMATRIX_SPARSE type and accessor macros are documented in
      :numref:`SUNMatrix.Sparse`.

      The previous function type :c:type:`CVDlsJacFn` is identical to
      :c:type:`CVLsJacFn`, and may still be used for backward-compatibility.
      However, this will be deprecated in future releases, so we recommend
      that users transition to the new function type name soon.


.. _CVODES.Usage.SIM.user_supplied.linsysFn:

Linear system construction (matrix-based linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With matrix-based linear solver modules, as an alternative to optionally
supplying a function for evaluating the Jacobian of the ODE right-hand
side function, the user may optionally supply a function of type ``CVLsLinSysFn``
for evaluating the linear system, :math:`M = I - \gamma J` (or an
approximation of it).  ``CVLsLinSysFn`` is defined as follows:

.. c:type:: int (*CVLsLinSysFn)(realtype t, N_Vector y, N_Vector fy, SUNMatrix M, booleantype jok, booleantype *jcur, realtype gamma, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

   This function computes the linear system matrix :math:`M = I - \gamma J` (or an approximation to it).

   **Arguments:**
      * ``t`` -- the current value of the independent variable.
      * ``y`` -- the current value of the dependent variable vector, namely the predicted value of :math:`y(t)`.
      * ``fy`` -- the current value of the vector :math:`f(t,y)`.
      * ``M`` -- the output linear system matrix.
      * ``jok`` -- an input flag indicating whether the Jacobian-related data needs to be updated. The ``jok`` flag enables reusing of Jacobian data across linear solves however, the user is responsible for storing Jacobian data for reuse. ``jok = SUNFALSE`` means that the Jacobian-related data must be recomputed from scratch. ``jok = SUNTRUE`` means that the Jacobian data, if saved from the previous call to this function, can be reused (with the current value of :math:`\gamma`). A call with ``jok = SUNTRUE`` can only occur after a call with ``jok = SUNFALSE``.
      * ``jcur`` -- a pointer to a flag which should be set to ``SUNTRUE`` if Jacobian data was recomputed, or set to ``SUNFALSE`` if Jacobian data was not recomputed, but saved data was still reused.
      * ``gamma`` -- the scalar :math:`\gamma` appearing in the matrix :math:`M = I - \gamma J`.
      * ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.
      * ``tmp1, tmp2, tmp3`` -- are pointers to memory allocated for variables of type ``N_Vector`` which can be used by a ``CVLsLinSysFn`` function as temporary storage or work space.

   **Return value:**
      Should return 0 if successful, a positive value if a recoverable error occurred (in which case CVODES will attempt to correct, while CVLS sets ``last_flag`` to ``CVLS_JACFUNC_RECVR``), or a negative value if it failed unrecoverably (in which case the integration is halted, :c:func:`CVode` returns ``CV_LSETUP_FAIL`` and CVLS sets ``last_flag`` to ``CVLS_JACFUNC_UNRECVR``).


.. _CVODES.Usage.SIM.user_supplied.jtimesFn:

Jacobian-vector product (matrix-free linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a matrix-free linear solver is to be used (i.e., a ``NULL``-valued
SUNMATRIX was supplied to :c:func:`CVodeSetLinearSolver`, the user may
provide a function of type :c:type:`CVLsJacTimesVecFn` in the following form,
to compute matrix-vector products :math:`Jv`. If such a function is not supplied,
the default is a difference quotient approximation to these products.

.. c:type:: int (*CVLsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp);

   This function computes the product :math:`Jv = \dfrac{\partial f(t,y)}{\partial y} v` (or an approximation to it).

   **Arguments:**
      * ``v`` -- the vector by which the Jacobian must be multiplied.
      * ``Jv`` -- the output vector computed.
      * ``t`` -- the current value of the independent variable.
      * ``y`` -- the current value of the dependent variable vector.
      * ``fy`` -- the current value of the vector :math:`f(t,y)`.
      * ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.
      * ``tmp`` -- a pointer to memory allocated for a variable of type ``N_Vector`` which can be used for work space.

   **Return value:**
      The value returned by the Jacobian-vector product function should be 0 if successful. Any other return value will result in an unrecoverable error of the generic Krylov solver, in which case the integration is halted.

   **Notes:**
      This function must return a value of :math:`Jv` that uses the *current*
      value of :math:`J`, i.e. as evaluated at the current :math:`(t,y)`.

      If the user's :c:type:`CVLsJacTimesVecFn` function uses difference quotient
      approximations, it may need to access quantities not in the argument
      list. These include the current step size, the error weights, etc.
      To obtain these, the user will need to add a pointer to ``cvode_mem``
      to ``user_data`` and then use the ``CVodeGet*`` functions described in
      :numref:`CVODES.Usage.SIM.optional_output.optout_main`. The unit roundoff can be accessed as
      ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

      The previous function type ``CVSpilsJacTimesVecFn`` is identical to
      :c:func:`CVLsJacTimesVecFn`, and may still be used for backward-compatibility.
      However, this will be deprecated in future releases, so we recommend
      that users transition to the new function type name soon.


  .. _CVODES.Usage.SIM.user_supplied.jtsetupFn:

Jacobian-vector product setup (matrix-free linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the user’s Jacobian-times-vector routine requires that any
Jacobian-related data be preprocessed or evaluated, then this needs to
be done in a user-supplied function of type :c:type:`CVLsJacTimesSetupFn`, defined as follows:

.. c:type:: int (*CVLsJacTimesSetupFn)(realtype t, N_Vector y, N_Vector fy, void *user_data);

   This function preprocesses and/or evaluates Jacobian-related data needed by the Jacobian-times-vector routine.

   **Arguments:**
      * ``t`` -- the current value of the independent variable.
      * ``y`` -- the current value of the dependent variable vector.
      * ``fy`` -- the current value of the vector :math:`f(t,y)`.
      * ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      The value returned by the Jacobian-vector setup function
      should be 0 if successful, positive for a recoverable error (in
      which case the step will be retried), or negative for an
      unrecoverable error (in which case the integration is halted).

   **Notes:**
      Each call to the Jacobian-vector setup function is preceded by a call to
      the :c:type:`CVRhsFn` user function with the same :math:`(t,y)` arguments.
      Thus, the setup function can use any auxiliary data that is computed
      and saved during the evaluation of the ODE right-hand side.

      If the user's ``CVLsJacTimesSetupFn`` function uses difference quotient
      approximations, it may need to access quantities not in the argument
      list. These include the current step size, the error weights, etc.
      To obtain these, the user will need to add a pointer to ``cvode_mem``
      to ``user_data`` and then use the ``CVodeGet*`` functions described in
      :numref:`CVODES.Usage.SIM.optional_output.optout_main`. The unit roundoff can be accessed as
      ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

      The previous function type ``CVSpilsJacTimesSetupFn`` is identical
      to :c:type:`CVLsJacTimesSetupFn`, and may still be used for
      backward-compatibility.  However, this will be deprecated in future
      releases, so we recommend that users transition to the new function
      type name soon.


.. _CVODES.Usage.SIM.user_supplied.psolveFn:

Preconditioner solve (iterative linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a user-supplied preconditioner is to be used with a ``SUNLinearSolver``
module, then the user must provide a function to solve the linear
system :math:`Pz = r`, where :math:`P` may be either a left or right
preconditioner matrix. Here :math:`P` should approximate (at least
crudely) the matrix :math:`M = I - \gamma J`, where
:math:`J = \dfrac{\partial f}{\partial y}`. If preconditioning is done on both
sides, the product of the two preconditioner matrices should approximate
:math:`M`. This function must be of type :c:type:`CVLsPrecSolveFn`, defined as follows:

.. c:type:: int (*CVLsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);

   This function solves the preconditioned system :math:`Pz = r`.

   **Arguments:**
      * ``t`` -- the current value of the independent variable.
      * ``y`` -- the current value of the dependent variable vector.
      * ``fy`` -- the current value of the vector :math:`f(t,y)`.
      * ``r`` -- the right-hand side vector of the linear system.
      * ``z`` -- the computed output vector.
      * ``gamma`` -- the scalar :math:`gamma` in the matrix given by :math:`M=I-\gamma J`.
      * ``delta`` -- an input tolerance to be used if an iterative method is employed in the solution. In that case, the residual vector :math:`Res = r - Pz` of the system should be made less than ``delta`` in the weighted :math:`l_2` norm, i.e., :math:`\sqrt{\sum_i (Res_i \cdot ewt_i)^2 } <` ``delta``. To obtain the ``N_Vector``  ``ewt``, call :c:func:`CVodeGetErrWeights`.
      * ``lr`` -- an input flag indicating whether the preconditioner solve function is to use the left preconditioner (``lr = 1``) or the right preconditioner (``lr = 2``).
      * ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      The value returned by the preconditioner solve function is a flag
      indicating whether it was successful.  This value should be 0 if successful,
      positive for a recoverable error (in which case the step will be retried), or
      negative for an unrecoverable error (in which case the integration is halted).

   **Notes:**
      The previous function type ``CVSpilsPrecSolveFn`` is identical to
      :c:type:`CVLsPrecSolveFn`, and may still be used for backward-compatibility.
      However, this will be deprecated in future releases, so we recommend
      that users transition to the new function type name soon.


.. _CVODES.Usage.SIM.user_supplied.precondFn:

Preconditioner setup (iterative linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the user’s preconditioner requires that any Jacobian-related data be
preprocessed or evaluated, then this needs to be done in a user-supplied
function of type , defined as follows:

.. c:type:: int (*CVLsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy, booleantype jok, booleantype *jcurPtr, realtype gamma, void *user_data);

   This function preprocesses and/or evaluates Jacobian-related data needed by the preconditioner.

   **Arguments:**
      * ``t`` -- the current value of the independent variable.
      * ``y`` -- the current value of the dependent variable vector, namely the predicted value of :math:`y(t)`.
      * ``fy`` -- the current value of the vector :math:`f(t,y)`.
      * ``jok`` -- an input flag indicating whether the Jacobian-related data needs to be updated. The ``jok`` argument provides for the reuse of Jacobian data in the preconditioner solve function. ``jok = SUNFALSE`` means that the Jacobian-related data must be recomputed from scratch. ``jok = SUNTRUE`` means that the Jacobian data, if saved from the previous call to this function, can be reused (with the current value of :math:`\gamma`). A call with ``jok = SUNTRUE`` can only occur after a call with ``jok = SUNFALSE``.
      * ``jcur`` -- a pointer to a flag which should be set to ``SUNTRUE`` if Jacobian data was recomputed, or set to ``SUNFALSE`` if Jacobian data was not recomputed, but saved data was still reused.
      * ``gamma`` -- the scalar :math:`\gamma` appearing in the matrix :math:`M = I - \gamma J`.
      * ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      The value returned by the preconditioner setup function is a flag
      indicating whether it was successful.  This value should be 0 if successful,
      positive for a recoverable error (in which case the step will be retried), or
      negative for an unrecoverable error (in which case the integration is halted).

   **Notes:**
      The operations performed by this function might include forming a crude
      approximate Jacobian and performing an LU factorization of the resulting
      approximation to :math:`M=I - \gamma J`.

      With the default nonlinear solver (the native SUNDIALS Newton method), each
      call to the preconditioner setup function is preceded by a call to the
      :c:type:`CVRhsFn` user function with the same :math:`(t,y)` arguments. Thus, the
      preconditioner setup function can use any auxiliary data that is computed and
      saved during the evaluation of the ODE right-hand side. In the case of a
      user-supplied or external nonlinear solver, this is also true if the nonlinear
      system function is evaluated prior to calling the linear solver setup function
      (see :numref:`SUNNonlinSol.API.SUNSuppliedFn` for more information).

      This function is not called in advance of every call to the preconditioner
      solve function, but rather is called only as often as needed to achieve
      convergence in the nonlinear solver.

      If the user's ``CVLsPrecSetupFn`` function uses difference quotient
      approximations, it may need to access quantities not in the call
      list. These include the current step size, the error weights, etc.
      To obtain these, the user will need to add a pointer to ``cvode_mem``
      to ``user_data`` and then use the ``CVodeGet*`` functions described in
      :numref:`CVODES.Usage.SIM.optional_output`. The unit roundoff can be accessed as
      ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

      The previous function type ``CVSpilsPrecSetupFn`` is identical to
      :c:type:`CVLsPrecSetupFn`, and may still be used for backward-compatibility.
      However, this will be deprecated in future releases, so we recommend
      that users transition to the new function type name soon.


.. _CVODES.Usage.purequad:

Integration of pure quadrature equations
========================================

CVODES allows the ODE system to include *pure quadratures*. In this case, it is
more efficient to treat the quadratures separately by excluding them from the
nonlinear solution stage. To do this, begin by excluding the quadrature
variables from the vector ``y`` and excluding the quadrature equations from
within ``res``. Thus a separate vector ``yQ`` of quadrature variables is to
satisfy :math:`(d/dt)`\ ``yQ`` = :math:`f_Q(t,y)`. The following is an overview
of the sequence of calls in a user’s main program in this situation. Steps that
are unchanged from the skeleton program presented in
:numref:`CVODES.Usage.SIM.skeleton_sim` are left unbolded.

  #. Initialize parallel or multi-threaded environment, if appropriate

  #. Create the SUNDIALS context object

  #. Set problem dimensions etc.

  #. Set vector of initial values

  #. Create CVODES object

  #. Initialize CVODES solver

  #. Specify integration tolerances

  #. Create matrix object

  #. Create linear solver object

  #. Set linear solver optional inputs

  #. Attach linear solver module

  #. Set optional inputs

  #. Create nonlinear solver object (*optional*)

  #. Attach nonlinear solver module (*optional*)

  #. Set nonlinear solver optional inputs (*optional*)

  #. **Set vector** ``yQ0`` **of initial values for quadrature variables**
      Typically, the quadrature variables should be initialized to 0.

  #. **Initialize quadrature integration**
      Call :c:func:`CVodeQuadInit` to specify the quadrature equation right-hand
      side function and to allocate internal memory related to quadrature integration.
      See :numref:`CVODES.Usage.purequad.quad_malloc` for details.

  #. **Set optional inputs for quadrature integration**
      Call :c:func:`CVodeSetQuadErrCon` to indicate whether or not quadrature
      variables shoule be used in the step size control mechanism, and to specify
      the integration tolerances for quadrature variables. See
      :numref:`CVODES.Usage.purequad.optional_inputs` for details.

  #. Specify rootfinding problem (*optional*)

  #. Advance solution in time

  #. **Extract quadrature variables**
      Call :c:func:`CVodeGetQuad` to obtain the values of the quadrature
      variables at the current time.

  #. Get optional outputs

  #. **Get quadrature optional outputs**
      Call ``CVodeGetQuad**`` functions to obtain optional output related
      to the integration of quadratures. See :numref:`CVODES.Usage.purequad.optional_output`
      for details.

  #. Deallocate memory for solution vector

  #. Free solver memory

  #. Free nonlinear solver memory (*optional*)

  #. Free linear solver and matrix memory

  #. Free the SUNContext object

  #. Finalize MPI, if used


:c:func:`CVodeQuadInit` can be called and quadrature-related optional inputs
can be set anywhere between the steps creating the CVODES object and advancing
the solution in time.

.. _CVODES.Usage.purequad.quad_malloc:

Quadrature initialization and deallocation functions
----------------------------------------------------

The function :c:func:`CVodeQuadInit` activates integration of quadrature
equations and allocates internal memory related to these calculations. The form
of the call to this function is as follows:

.. c:function:: int CVodeQuadInit(void * cvode_mem, CVQuadRhsFn fQ, N_Vector yQ0)

   The function ``CVodeQuadInit`` provides required problem specifications,  allocates internal memory, and initializes quadrature integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block returned by :c:func:`CVodeCreate`.
     * ``fQ`` -- is the C function which computes :math:`f_Q` , the right-hand side of the quadrature equations.
     * ``yQ0`` -- is the initial value of ``yQ`` typically ``yQ0`` has all zero components.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to ``CVodeQuadInit`` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory was not initialized by a prior call to :c:func:`CVodeCreate`.
     * ``CV_MEM_FAIL`` -- A memory allocation request failed.

   **Notes:**
      If an error occurred, ``CVodeQuadInit`` also sends an error message to the  error handler function.

In terms of the number of quadrature variables :math:`N_q` and maximum method
order ``maxord``, the size of the real workspace is increased as follows:

-  Base value: :math:`\texttt{lenrw} = \texttt{lenrw} + (\texttt{maxord}+5) N_q`
-  If using :c:func:`CVodeSVtolerances` (see :c:func:`CVodeSetQuadErrCon`):  :math:`\texttt{lenrw} = \texttt{lenrw} + N_q`

the size of the integer workspace is increased as follows:

-  Base value: :math:`\texttt{leniw} = \texttt{leniw} + (\texttt{maxord}+5) N_q`
-  If using :c:func:`CVodeSVtolerances`: :math:`\texttt{leniw} = \texttt{leniw} + N_q`

The function :c:func:`CVodeQuadReInit`, useful during the solution of a sequence of
problems of same size, reinitializes the quadrature-related internal memory and
must follow a call to :c:func:`CVodeQuadInit` (and maybe a call to :c:func:`CVodeReInit`).
The number ``Nq`` of quadratures is assumed to be unchanged from the prior call
to :c:func:`CVodeQuadInit`. The call to the :c:func:`CVodeQuadReInit` function has the
following form:

.. c:function:: int CVodeQuadReInit(void * cvode_mem, N_Vector yQ0)

   The function ``CVodeQuadReInit`` provides required problem specifications  and reinitializes the quadrature integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``yQ0`` -- is the initial value of ``yQ``.

   **Return value:**
     * ``CV_SUCCESS`` -- The call to ``CVodeReInit`` was successful.
     * ``CV_MEM_NULL`` -- The CVODES memory was not initialized by a prior call to ``CVodeCreate``.
     * ``CV_NO_QUAD`` -- Memory space for the quadrature integration was not allocated by a prior call to ``CVodeQuadInit``.

   **Notes:**
      If an error occurred, ``CVodeQuadReInit`` also sends an error message to the  error handler function.

.. c:function:: void CVodeQuadFree(void * cvode_mem)

   The function ``CVodeQuadFree`` frees the memory allocated for quadrature integration.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block

  **Return value:**
     * The function has no return value.

  **Notes:**
     In general, ``CVodeQuadFree`` need not be called by the user as it is invoked automatically
     by :c:func:`CVodeFree`.


CVODES solver function
----------------------

Even if quadrature integration was enabled, the call to the main solver function
:c:func:`CVode` is exactly the same as in :numref:`CVODES.Usage.SIM`.
However, in this case the return value ``flag`` can also be one of the
following:

* The quadrature right-hand side function failed in an unrecoverable manner.
* The quadrature right-hand side function failed at the first call.
* Convergence test failures occurred too many times due to repeated recoverable errors in the quadrature right-hand
  side function. This value will also be returned if the quadrature right-hand side function had repeated recoverable
  errors during the estimation of an initial step size (assuming the quadrature variables are included in the error
  tests).
* The quadrature right-hand function had a recoverable error, but no recovery was possible. This failure mode is rare,
  as it can occur only if the quadrature right-hand side function fails recoverably after an error test failed while at
  order one.

.. _CVODES.Usage.purequad.quad_get:

Quadrature extraction functions
-------------------------------

If quadrature integration has been initialized by a call to :c:func:`CVodeQuadInit`,
or reinitialized by a call to :c:func:`CVodeQuadReInit`, then CVODES computes both a
solution and quadratures at time ``t``. However, :c:func:`CVode` will still return
only the solution :math:`y` in ``yout``. Solution quadratures can be obtained
using the following function:

.. c:function:: int CVodeGetQuad(void * cvode_mem, realtype tret, N_Vector yQ)

   The function ``CVodeGetQuad`` returns the quadrature solution vector after a  successful return from ``CVode``.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by ``CVodeInit``.
     * ``tret`` -- the time reached by the solver output.
     * ``yQ`` -- the computed quadrature vector. This vector must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- ``CVodeGetQuad`` was successful.
     * ``CV_MEM_NULL`` -- ``cvode_mem`` was ``NULL``.
     * ``CV_NO_QUAD`` -- Quadrature integration was not initialized.
     * ``CV_BAD_DKY`` -- ``yQ`` is ``NULL``.

   **Notes:**
      In case of an error return, an error message is also sent to the error handler function.


The function :c:func:`CVodeGetQuadDky` computes the ``k``-th derivatives of the
interpolating polynomials for the quadrature variables at time ``t``. This
function is called by :c:func:`CVodeGetQuad` with ``k = 0`` and with the current time
at which :c:func:`CVode` has returned, but may also be called directly by the user.

.. c:function:: int CVodeGetQuadDky(void * cvode_mem, realtype t, int k, N_Vector dkyQ)

   The function ``CVodeGetQuadDky`` returns derivatives of the quadrature solution  vector after a successful return from :c:func:`CVode`.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the memory previously allocated by :c:func:`CVodeInit`.
     * ``t`` -- the time at which quadrature information is requested. The time ``t`` must fall within the interval defined by the last successful step taken by CVODES.
     * ``k`` -- order of the requested derivative. This must be :math:`\leq` ``qlast``.
     * ``dkyQ`` -- the vector containing the derivative. This vector must be allocated by the user.

   **Return value:**
     * ``CV_SUCCESS`` -- ``CVodeGetQuadDky`` succeeded.
     * ``CV_MEM_NULL`` -- The pointer to ``cvode_mem`` was ``NULL``.
     * ``CV_NO_QUAD`` -- Quadrature integration was not initialized.
     * ``CV_BAD_DKY`` -- The vector ``dkyQ`` is ``NULL``.
     * ``CV_BAD_K`` -- ``k`` is not in the range :math:`0, 1, \ldots,` ``qlast``.
     * ``CV_BAD_T`` -- The time ``t`` is not in the allowed range.

   **Notes:**
      In case of an error return, an error message is also sent to the error handler function.


.. _CVODES.Usage.purequad.optional_inputs:

Optional inputs for quadrature integration
------------------------------------------

CVODES provides the following optional input functions to control the integration
of quadrature equations.


.. c:function:: int CVodeSetQuadErrCon(void * cvode_mem, booleantype errconQ)

   The function ``CVodeSetQuadErrCon`` specifies whether or not the  quadrature variables are to be used in the step size control mechanism within CVODES.  If they are, the user must call :c:func:`CVodeQuadSStolerances`  or :c:func:`CVodeQuadSVtolerances` to specify the  integration tolerances for the quadrature variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``errconQ`` -- specifies whether quadrature variables are included ``SUNTRUE`` or not ``SUNFALSE`` in the error control mechanism.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUAD`` -- Quadrature integration has not been initialized.

   **Notes:**
      By default, ``errconQ`` is set to ``SUNFALSE``.

   .. warning::

      It is illegal to call ``CVodeSetQuadErrCon`` before a call  to ``CVodeQuadInit``.


If the quadrature variables are part of the step size control mechanism,
one of the following functions must be called to specify the
integration tolerances for quadrature variables.

.. c:function:: int CVodeQuadSVtolerances(void * cvode_mem, realtype reltolQ, realtype abstolQ)

   The function ``CVodeQuadSStolerances`` specifies scalar relative and absolute  tolerances.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``reltolQ`` --  tolerances is the scalar relative error tolerance.
     * ``abstolQ`` -- is the scalar absolute error tolerance.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional value has been successfully set.
     * ``CV_NO_QUAD`` -- Quadrature integration was not initialized.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_ILL_INPUT`` -- One of the input tolerances was negative.


.. _CVODES.Usage.purequad.optional_output:

Optional outputs for quadrature integration
--------------------------------------------

CVODES provides the following functions that can be used to obtain solver
performance information related to quadrature integration.

.. c:function:: int CVodeGetQuadNumRhsEvals(void * cvode_mem, long int nfQevals)

   The function ``CVodeGetQuadNumRhsEvals`` returns the  number of calls made to the user's quadrature right-hand side function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfQevals`` -- number of calls made to the user's ``fQ`` function.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUAD`` -- Quadrature integration has not been initialized.


.. c:function:: int CVodeGetQuadNumErrTestFails(void * cvode_mem, long int nQetfails)

   The function ``CVodeGetQuadNumErrTestFails`` returns the  number of local error test failures due to quadrature variables.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nQetfails`` -- number of error test failures due to quadrature variables.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUAD`` -- Quadrature integration has not been initialized.


.. c:function:: int CVodeGetQuadErrWeights(void * cvode_mem, N_Vector eQweight)

   The function ``CVodeGetQuadErrWeights`` returns the quadrature error weights  at the current time.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``eQweight`` -- quadrature error weights at the current time.

   **Return value:**
     * ``CV_SUCCESS`` -- The optional output value has been successfully set.
     * ``CV_MEM_NULL`` -- The ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUAD`` -- Quadrature integration has not been initialized.

   **Notes:**
      The user must allocate memory for ``eQweight``.  If quadratures were not included in the error control mechanism (through a  call to ``CVodeSetQuadErrCon`` with ``errconQ = SUNTRUE``),  ``CVodeGetQuadErrWeights`` does not set the ``eQweight`` vector.


.. c:function:: int CVodeGetQuadStats(void * cvode_mem, long int nfQevals, long int nQetfails)

   The function ``CVodeGetQuadStats`` returns the CVODES integrator statistics  as a group.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfQevals`` -- number of calls to the user's ``fQ`` function.
     * ``nQetfails`` -- number of error test failures due to quadrature variables.

   **Return value:**
     * ``CV_SUCCESS`` -- the optional output values have been successfully set.
     * ``CV_MEM_NULL`` -- the ``cvode_mem`` pointer is ``NULL``.
     * ``CV_NO_QUAD`` -- Quadrature integration has not been initialized.


.. _CVODES.Usage.purequad.user_supplied:

User supplied functions for quadrature integration
--------------------------------------------------

For integration of quadrature equations, the user must provide a function
that defines the right-hand side of the quadrature equations (in other words,
the integrand function of the integral that must be evaluated). This function
must be of type :c:type:`CVQuadRhsFn` defined as follows:

.. c:type:: int (*CVQuadRhsFn)(realtype t, N_Vector y, N_Vector yQdot, void *user_data)

   This function computes the quadrature equation right-hand side for a given value
   of the independent variable :math:`t` and state vector :math:`y`.

   **Arguments:**
     * ``t`` -- is the current value of the independent variable.
     * ``y`` -- is the current value of the dependent variable vector, :math:`y(t)`.
     * ``yQdot`` -- is the output vector :math:`f_Q(t,y)`.
     * ``user_data`` -- is the ``user_data`` pointer passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      A ``CVQuadRhsFn`` should return 0 if successful, a positive value if a
      recoverable error occurred (in which case CVODES will attempt to correct),
      or a negative value if it failed unrecoverably (in which case the
      integration is halted and ``CV_QRHSFUNC_FAIL`` is returned).

   **Notes:**
      Allocation of memory for ``yQdot`` is automatically handled within CVODES.

      Both ``y`` and ``yQdot`` are of type ``N_Vector``, but they typically have
      different internal representations. It is the user’s responsibility to
      access the vector data consistently (including the use of the correct
      accessor macros from each ``N_Vector`` implementation). For the sake of
      computational efficiency, the vector functions in the two ``N_Vector``
      implementations provided with CVODES do not perform any consistency checks
      with respect to their ``N_Vector`` arguments.

      There are two situations in which recovery is not possible even if
      ``CVQuadRhsFn`` function returns a recoverable error flag. One is when
      this occurs at the very first call to the ``CVQuadRhsFn`` (in which case
      CVODES returns ``CV_FIRST_QRHSFUNC_ERR``). The other is when a recoverable
      error is reported by ``CVQuadRhsFn`` after an error test failure, while
      the linear multistep method order is equal to 1 (in which case CVODES
      returns ``CV_UNREC_QRHSFUNC_ERR``).


.. _CVODES.Usage.SIM.precond:

Preconditioner modules
----------------------

The efficiency of Krylov iterative methods for the solution of linear
systems can be greatly enhanced through preconditioning. For problems in
which the user cannot define a more effective, problem-specific
preconditioner, CVODES provides a banded preconditioner in the module
CVBANDPRE and a band-block-diagonal preconditioner module
CVBBDPRE.

.. _CVODES.Usage.SIM.precond.cvbandpre:

A serial banded preconditioner module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This preconditioner provides a band matrix preconditioner for use with
iterative ``SUNLinearSolver`` modules through the CVLS linear solver
interface, in a serial setting. It uses difference quotients of the ODE
right-hand side function :math:`f` to generate a band matrix of bandwidth
:math:`m_l + m_u + 1`, where the number of super-diagonals (:math:`m_u`,
the upper half-bandwidth) and sub-diagonals (:math:`m_l`, the lower
half-bandwidth) are specified by the user, and uses this to form a
preconditioner for use with the Krylov linear solver. Although this
matrix is intended to approximate the Jacobian
:math:`\dfrac{\partial f}{\partial y}`, it may be a very crude approximation.
The true Jacobian need not be banded, or its true bandwidth may be
larger than :math:`m_l + m_u + 1`, as long as the banded approximation
generated here is sufficiently accurate to speed convergence as a
preconditioner.

In order to use the CVBANDPRE module, the user need not define any
additional functions. Aside from the header files required for the
integration of the ODE problem (see :numref:`CVODES.Usage.SIM.header_sim`), to use
the CVBANDPRE module, the main program must include the header file
``cvode_bandpre.h`` which declares the needed function prototypes.
The following is a summary of the usage of this module. Steps that are
changed from the skeleton program presented in
:numref:`CVODES.Usage.SIM.skeleton_sim` are shown in bold.

  #. Initialize multi-threaded environment, if appropriate

  #. Create the ``SUNContext`` object.

  #. Set problem dimensions etc.

  #. Set vector of initial values

  #. Create CVODES object

  #. Initialize CVODES solver

  #. Specify integration tolerances

  #. **Create linear solver object**

      When creating the iterative linear solver object, specify the type of
      preconditioning (``SUN_PREC_LEFT`` or ``SUN_PREC_RIGHT``) to use.

  #. Set linear solver optional inputs

  #. Attach linear solver module

  #. **Initialize the CVBANDPRE preconditioner module**

     Specify the upper and lower half-bandwidths (``mu`` and ``ml``, respectively) and call

     .. code-block:: c

        flag = CVBandPrecInit(cvode_mem, N, mu, ml);

     to allocate memory and initialize the internal preconditioner data.

  #. Set optional inputs.

     Note that the user should not overwrite the preconditioner setup function or solve function through calls to the :c:func:`CVodeSetPreconditioner` optional input function.

  #. Create nonlinear solver object

  #. Attach nonlinear solver module

  #. Set nonlinear solver optional inputs

  #. Specify rootfinding problem

  #. Advance solution in time

  #. **Get optional outputs**

     Additional optional outputs associated with CVBANDPRE are available by way of two routines described below, :c:func:`CVBandPrecGetWorkSpace` and :c:func:`CVBandPrecGetNumRhsEvals`.

  #. Deallocate memory for solution vector

  #. Free solver memory

  #. Free nonlinear solver memory

  #. Free linear solver memory

  #. Free the SUNContext object


The CVBANDPRE preconditioner module is initialized and attached by
calling the following function:

.. c:function:: int CVBandPrecInit(void* cvode_mem, sunindextype N, sunindextype mu, sunindextype ml)

   The function ``CVBandPrecInit`` initializes the CVBANDPRE preconditioner  and allocates required (internal) memory for it.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``N`` -- problem dimension.
     * ``mu`` -- upper half-bandwidth of the Jacobian approximation.
     * ``ml`` -- lower half-bandwidth of the Jacobian approximation.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The call to ``CVBandPrecInit`` was successful.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
     * ``CVLS_MEM_FAIL`` -- A memory allocation request has failed.
     * ``CVLS_LMEM_NULL`` -- A CVLS linear solver memory was not attached.
     * ``CVLS_ILL_INPUT`` -- The supplied vector implementation was not compatible with block    band preconditioner.

   **Notes:**
      The banded approximate Jacobian will have nonzero elements only in locations :math:`(i,j)` with :math:`\text{ml} \leq j-i \leq \text{mu}`.


The following two optional output functions are available for use with
the CVBANDPRE module:

.. c:function:: int CVBandPrecGetWorkSpace(void* cvode_mem, long int *lenrwBP, long int *leniwBP)

   The function ``CVBandPrecGetWorkSpace`` returns the sizes of  the CVBANDPRE real and integer workspaces.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``lenrwBP`` -- the number of ``realtype`` values in the CVBANDPRE workspace.
     * ``leniwBP`` -- the number of integer values in the CVBANDPRE workspace.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output values have been successfully set.
     * ``CVLS_PMEM_NULL`` -- The CVBANDPRE preconditioner has not been initialized.

   **Notes:**
      The workspace requirements reported by this routine correspond only  to memory allocated within the CVBANDPRE module (the banded  matrix approximation, banded ``SUNLinearSolver`` object, and temporary vectors).

      The workspaces referred to here exist in addition to those given by the  corresponding function :c:func:`CVodeGetLinWorkSpace`.


.. c:function:: int CVBandPrecGetNumRhsEvals(void* cvode_mem, long int *nfevalsBP)

   The function ``CVBandPrecGetNumRhsEvals`` returns the  number of calls made to the user-supplied right-hand side function for  the finite difference banded Jacobian approximation used within  the preconditioner setup function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``nfevalsBP`` -- the number of calls to the user right-hand side function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_PMEM_NULL`` -- The CVBANDPRE preconditioner has not been initialized.

   **Notes:**
      The counter ``nfevalsBP`` is distinct from the counter ``nfevalsLS`` returned by the corresponding function :c:func:`CVodeGetNumLinRhsEvals` and ``nfevals`` returned by :c:func:`CVodeGetNumRhsEvals`.The total number of right-hand side function evaluations is the sum of all three of these counters.


.. _CVODES.Usage.SIM.precond.cvbbdpre:

A parallel band-block-diagonal preconditioner module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A principal reason for using a parallel ODE solver such as CVODES lies
in the solution of partial differential equations (PDEs).  Moreover,
the use of a Krylov iterative method for the solution of many such
problems is motivated by the nature of the underlying linear system of
equations :eq:`CVODES_Newton` that must be solved at each time step.  The
linear algebraic system is large, sparse, and structured. However, if
a Krylov iterative method is to be effective in this setting, then a
nontrivial preconditioner needs to be used.  Otherwise, the rate of
convergence of the Krylov iterative method is usually unacceptably
slow.  Unfortunately, an effective preconditioner tends to be
problem-specific.

However, we have developed one type of preconditioner that treats a
rather broad class of PDE-based problems.  It has been successfully
used for several realistic, large-scale problems :cite:p:`HiTa:98` and is
included in a software module within the CVODES package. This module
works with the parallel vector module NVECTOR_PARALLEL and is usable with any of
the Krylov iterative linear solvers through the CVLS interface.
It generates a preconditioner that is a block-diagonal matrix with
each block being a band matrix. The blocks need not have the same
number of super- and sub-diagonals and these numbers may vary from
block to block. This Band-Block-Diagonal Preconditioner module is
called CVBBDPRE.

One way to envision these preconditioners is to think of the domain of
the computational PDE problem as being subdivided into :math:`M` non-overlapping
subdomains.  Each of these subdomains is then assigned to one of the
:math:`M` processes to be used to solve the ODE system. The basic idea is
to isolate the preconditioning so that it is local to each process,
and also to use a (possibly cheaper) approximate right-hand side
function. This requires the definition of a new function :math:`g(t,y)`
which approximates the function :math:`f(t, y)` in the definition of the ODE
system :eq:`CVODES_ivp`. However, the user may set :math:`g = f`.  Corresponding
to the domain decomposition, there is a decomposition of the solution
vector :math:`y` into :math:`M` disjoint blocks :math:`y_m`, and a decomposition of :math:`g`
into blocks :math:`g_m`.  The block :math:`g_m` depends both on :math:`y_m` and on
components of blocks :math:`y_{m'}` associated with neighboring subdomains
(so-called ghost-cell data).  Let :math:`\bar{y}_m` denote :math:`y_m` augmented
with those other components on which :math:`g_m` depends.  Then we have

.. math::

   g(t,y) = \begin{bmatrix} g_1(t,\bar{y}_1) & g_2(t,\bar{y}_2) & \cdots & g_M(t,\bar{y}_M) \end{bmatrix}^T

and each of the blocks :math:`g_m(t, \bar{y}_m)` is uncoupled from the others.

The preconditioner associated with this decomposition has the form

.. math::

   P = \begin{bmatrix} P_1 & & & \\ & P_2 & & \\ &  & \ddots & \\ & & & P_M\end{bmatrix}

where

.. math::

   P_m \approx I - \gamma J_m

and :math:`J_m` is a difference quotient approximation to
:math:`\partial g_m/\partial y_m`. This matrix is taken to be banded, with
upper and lower half-bandwidths ``mudq`` and ``mldq`` defined as
the number of non-zero diagonals above and below the main diagonal,
respectively. The difference quotient approximation is computed using
:math:`\texttt{mudq} + \texttt{mldq} + 2` evaluations of :math:`g_m`, but only a matrix
of bandwidth :math:`\texttt{mukeep} + \texttt{mlkeep} + 1` is retained.
Neither pair of parameters need be the true half-bandwidths of the Jacobian of the
local block of :math:`g`, if smaller values provide a more efficient
preconditioner. The solution of the complete linear system

.. math::

   Px = b

reduces to solving each of the equations

.. math::

   P_m x_m = b_m

and this is done by banded LU factorization of :math:`P_m` followed by a banded
backsolve.

Similar block-diagonal preconditioners could be considered with different
treatments of the blocks :math:`P_m`. For example, incomplete LU factorization or
an iterative method could be used instead of banded LU factorization.

The CVBBDPRE module calls two user-provided functions to construct :math:`P`:
a required function ``gloc`` (of type :c:type:`CVLocalFn`) which approximates
the right-hand side function :math:`g(t,y) \approx f(t,y)` and which is computed locally,
and an optional function ``cfn`` (of type :c:type:`CVCommFn`) which performs
all interprocess communication necessary to evaluate the approximate right-hand
side :math:`g`.  These are in addition to the user-supplied right-hand side function
:math:`f`.  Both functions take as input the same pointer ``user_data`` that is passed
by the user to :c:func:`CVodeSetUserData` and that was passed to the user's function :math:`f`.
The user is responsible for providing space (presumably within ``user_data``)
for components of :math:`y` that are communicated between processes by ``cfn``, and
that are then used by ``gloc``, which should not do any communication.

.. c:type:: int (*CVLocalFn)(sunindextype Nlocal, realtype t, N_Vector y, N_Vector glocal, void *user_data);

   This ``gloc`` function computes :math:`g(t,y)`. It loads the vector ``glocal`` as a function of ``t`` and ``y``.

   **Arguments:**
      * ``Nlocal`` -- the local vector length.
      * ``t`` -- the value of the independent variable.
      * ``y`` -- the dependent variable.
      * ``glocal`` -- the output vector.
      * ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      A ``CVLocalFn`` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct), or a negative
      value if it failed unrecoverably (in which case the integration is halted and
      :c:func:`CVode` returns ``CV_LSETUP_FAIL``).

   **Notes:**
      This function must assume that all interprocess communication of data needed to
      calculate ``glocal`` has already been done, and that this data is accessible within
      ``user_data``.

      The case where :math:`g` is mathematically identical to :math:`f` is allowed.


.. c:type:: int (*CVCommFn)(sunindextype Nlocal, realtype t, N_Vector y, void *user_data);

   This ``cfn`` function performs all interprocess communication necessary for the execution of the ``gloc`` function above, using the input vector ``y``.

   **Arguments:**
      * ``Nlocal`` -- the local vector length.
      * ``t`` -- the value of the independent variable.
      * ``y`` -- the dependent variable.
      * ``user_data`` -- a pointer to user data, the same as the ``user_data`` parameter passed to :c:func:`CVodeSetUserData`.

   **Return value:**
      A ``CVCommFn`` should return 0 if successful, a positive value if a recoverable
      error occurred (in which case CVODES will attempt to correct), or a negative
      value if it failed unrecoverably (in which case the integration is halted and
      :c:func:`CVode` returns ``CV_LSETUP_FAIL``).

   **Notes:**
      The ``cfn`` function is expected to save communicated data in space defined
      within the data structure ``user_data``.

      Each call to the ``cfn`` function is preceded by a call to the right-hand side
      function :math:`f` with the same :math:`(t,y)` arguments.  Thus, ``cfn`` can omit
      any communication done by :math:`f` if relevant to the evaluation of ``glocal``.
      If all necessary communication was done in :math:`f`, then ``cfn = NULL``
      can be passed in the call to :c:func:`CVBBDPrecInit` (see below).


Besides the header files required for the integration of the ODE problem
(see :numref:`CVODES.Usage.SIM.header_sim`),  to use the CVBBDPRE module, the main program
must include the header file ``cvode_bbdpre.h`` which declares the needed
function prototypes.

The following is a summary of the usage of this module. Steps that are
changed from the skeleton program presented in
:numref:`CVODES.Usage.SIM.skeleton_sim` are shown in bold.

  #. Initialize MPI environment

  #. Create the ``SUNContext`` object

  #. Set problem dimensions etc.

  #. Set vector of initial values

  #. Create CVODES object

  #. Initialize CVODES solver

  #. Specify integration tolerances

  #. **Create linear solver object**

     When creating the iterative linear solver object, specify the type
     of preconditioning (``SUN_PREC_LEFT`` or ``SUN_PREC_RIGHT``) to use.

  #. Set linear solver optional inputs

  #. Attach linear solver module

  #. **Initialize the CVBBDPRE preconditioner module**

     Specify the upper and lower half-bandwidths ``mudq`` and ``mldq``, and ``mukeep`` and ``mlkeep``, and call

     .. code-block:: c

        flag = CVBBDPrecInit(&cvode_mem, local_N, mudq, mldq,
                             &mukeep, mlkeep, dqrely, gloc, cfn);

     to allocate memory and initialize the internal preconditioner data.
     The last two arguments of :c:func:`CVBBDPrecInit` are the two user-supplied
     functions described above.

  #. Set optional inputs

     Note that the user should not overwrite the preconditioner setup function
     or solve function through calls to the :c:func:`CVodeSetPreconditioner`
     optional input function.

  #. Create nonlinear solver object

  #. Attach nonlinear solver module

  #. Set nonlinear solver optional inputs

  #. Specify rootfinding problem

  #. Advance solution in time

  #. **Get optional outputs**

     Additional optional outputs associated with CVBBDPRE are available by
     way of two routines described below, :c:func:`CVBBDPrecGetWorkSpace`
     and :c:func:`CVBBDPrecGetNumGfnEvals`.

  #. Deallocate memory for solution vector

  #. Free solver memory

  #. Free nonlinear solver memory

  #. Free linear solver memory

  #. Free the ``SUNContext`` object

  #. Finalize MPI



The user-callable functions that initialize or re-initialize
the CVBBDPRE preconditioner module are described next.


.. c:function:: int CVBBDPrecInit(void* cvode_mem, sunindextype local_N, sunindextype mudq, sunindextype mldq, sunindextype mukeep, sunindextype mlkeep, realtype dqrely, CVLocalFn gloc, CVCommFn cfn)

   The function ``CVBBDPrecInit`` initializes and allocates (internal) memory for the CVBBDPRE preconditioner.

   **Arguments:**
      * ``cvode_mem`` -- pointer to the CVODES memory block.
      * ``local_N`` -- local vector length.
      * ``mudq`` -- upper half-bandwidth to be used in the difference quotient Jacobian approximation.
      * ``mldq`` -- lower half-bandwidth to be used in the difference quotient Jacobian approximation.
      * ``mukeep`` -- upper half-bandwidth of the retained banded approximate Jacobian block.
      * ``mlkeep`` -- lower half-bandwidth of the retained banded approximate Jacobian block.
      * ``dqrely`` -- the relative increment in components of :math:`y` used in the difference quotient approximations. The default is :math:`\texttt{dqrely} = \sqrt{\text{unit roundoff}}`, which can be specified by passing ``dqrely = 0.0``.
      * ``gloc`` -- the :c:type:`CVLocalFn` function which computes the approximation :math:`g(t,y) \approx f(t,y)`.
      * ``cfn`` -- the :c:type:`CVCommFn` which performs all interprocess communication required for the computation of :math:`g(t,y)`.

   **Return value:**
      * ``CVLS_SUCCESS`` -- The function was successful
      * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``.
      * ``CVLS_MEM_FAIL`` -- A memory allocation request has failed.
      * ``CVLS_LMEM_NULL`` -- A CVLS linear solver memory was not attached.
      * ``CVLS_ILL_INPUT`` -- The supplied vector implementation was not compatible with block band preconditioner.

   **Notes:**
      If one of the half-bandwidths ``mudq`` or ``mldq`` to be used in the
      difference quotient calculation of the approximate Jacobian is negative or
      exceeds the value ``local_N - 1``, it is replaced by ``0`` or ``local_N - 1`` accordingly.

      The half-bandwidths ``mudq`` and ``mldq`` need not be the true
      half-bandwidths of the Jacobian of the local block of :math:`g`
      when smaller values may provide a greater efficiency.

      Also, the half-bandwidths ``mukeep`` and ``mlkeep`` of the retained
      banded approximate Jacobian block may be even smaller,
      to reduce storage and computational costs further.

      For all four half-bandwidths, the values need not be the
      same on every processor.


The CVBBDPRE module also provides a reinitialization function to allow
solving a sequence of problems of the same size, with the same linear solver
choice, provided there is no change in ``local_N``, ``mukeep``, or ``mlkeep``.
After solving one problem, and after calling :c:func:`CVodeReInit` to
re-initialize CVODES for a subsequent problem, a call to :c:func:`CVBBDPrecReInit`
can be made to change any of the following: the half-bandwidths ``mudq`` and
``mldq`` used in the difference-quotient Jacobian approximations, the relative
increment ``dqrely``, or one of the user-supplied functions ``gloc`` and ``cfn``.
If there is a change in any of the linear solver inputs, an additional call
to the "set" routines provided by the SUNLinearSolver module, and/or
one or more of the corresponding CVLS "set" functions, must
also be made (in the proper order).


.. c:function:: int CVBBDPrecReInit(void* cvode_mem, sunindextype mudq, sunindextype mldq, realtype dqrely)

   The function ``CVBBDPrecReInit`` re-initializes the CVBBDPRE preconditioner.

   **Arguments:**
      * ``cvode_mem`` -- pointer to the CVODES memory block.
      * ``mudq`` -- upper half-bandwidth to be used in the difference quotient Jacobian approximation.
      * ``mldq`` -- lower half-bandwidth to be used in the difference quotient Jacobian approximation.
      * ``dqrely`` -- the relative increment in components of

   **Return value:**
      * ``CVLS_SUCCESS`` -- The function was successful
      * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer is ``NULL``. ``cvode_mem`` pointer was ``NULL``.
      * ``CVLS_LMEM_NULL`` -- A CVLS linear solver memory was not attached.
      * ``CVLS_PMEM_NULL`` -- The function :c:func:`CVBBDPrecInit` was not previously called

   **Notes:**
      If one of the half-bandwidths ``mudq`` or ``mldq`` is negative or  exceeds the value ``local_N-1``, it is replaced by ``0`` or  ``local_N-1`` accordingly.


The following two optional output functions are available for use with
the CVBBDPRE module:


.. c:function:: int CVBBDPrecGetWorkSpace(void* cvode_mem, long int *lenrwBBDP, long int *leniwBBDP)

   The function ``CVBBDPrecGetWorkSpace`` returns the local  CVBBDPRE real and integer workspace sizes.

   **Arguments:**
      * ``cvode_mem`` -- pointer to the CVODES memory block.
      * ``lenrwBBDP`` -- local number of ``realtype`` values in the CVBBDPRE workspace.
      * ``leniwBBDP`` -- local number of integer values in the CVBBDPRE workspace.

   **Return value:**
      * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
      * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer was ``NULL``.
      * ``CVLS_PMEM_NULL`` -- The CVBBDPRE preconditioner has not been initialized.

   **Notes:**
      The workspace requirements reported by this routine correspond only  to memory allocated within the CVBBDPRE module (the banded  matrix approximation, banded ``SUNLinearSolver`` object, temporary vectors).  These values are local to each process.  The workspaces referred to here exist in addition to those given by the  corresponding function ``CVodeGetLinWorkSpace``.

.. c:function:: int CVBBDPrecGetNumGfnEvals(void* cvode_mem, long int *ngevalsBBDP)

   The function ``CVBBDPrecGetNumGfnEvals`` returns the  number of calls made to the user-supplied ``gloc`` function due to the  finite difference approximation of the Jacobian blocks used within  the preconditioner setup function.

   **Arguments:**
     * ``cvode_mem`` -- pointer to the CVODES memory block.
     * ``ngevalsBBDP`` -- the number of calls made to the user-supplied ``gloc`` function due to the
       finite difference approximation of the Jacobian blocks used within
       the preconditioner setup function.

   **Return value:**
     * ``CVLS_SUCCESS`` -- The optional output value has been successfully set.
     * ``CVLS_MEM_NULL`` --  The ``cvode_mem`` pointer was ``NULL``.
     * ``CVLS_PMEM_NULL`` -- The CVBBDPRE preconditioner has not been initialized.



In addition to the ``ngevalsBBDP`` ``gloc`` evaluations,
the costs associated with CVBBDPRE also include ``nlinsetups`` LU
factorizations, ``nlinsetups`` calls to ``cfn``, ``npsolves`` banded
backsolve calls, and ``nfevalsLS`` right-hand side function evaluations,
where ``nlinsetups`` is an optional CVODES output and ``npsolves`` and
``nfevalsLS`` are linear solver optional outputs (see :numref:`CVODES.Usage.SIM.optional_output`).
