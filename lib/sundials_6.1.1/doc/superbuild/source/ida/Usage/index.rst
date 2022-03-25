.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _IDA.Usage.CC:

**************************
Using IDA for IVP Solution
**************************

This chapter is concerned with the use of IDA for the integration of DAEs.

The following sections treat the header files and the layout of the user’s main
program, and provide descriptions of the IDA user-callable functions and
user-supplied functions. The sample programs described in the companion document
:cite:p:`ida_ex` may also be helpful. Those codes may be used as templates (with
the removal of some lines used in testing) and are included in the IDA package.

IDA uses various constants for both input and output. These are defined as
needed in this chapter, but for convenience are also listed separately in
:numref:`IDA.Constants`.

The user should be aware that not all ``SUNLinearSolver`` and ``SUNMatrix``
objects are compatible with all ``N_Vector`` implementations. Details on
compatibility are given in the documentation for each ``SUNMatrix`` (Chapter
:numref:`SUNMatrix`) and ``SUNLinearSolver`` (Chapter :numref:`SUNLinSol`)
implementation. For example, ``NVECTOR_PARALLEL`` is not compatible with the
dense, banded, or sparse ``SUNMatrix`` types, or with the corresponding dense,
banded, or sparse ``SUNLinearSolver`` objects. Please check Chapters
:numref:`SUNMatrix` and :numref:`SUNLinSol` to verify compatibility between
these objects. In addition to that documentation, we note that the IDABBDPRE
preconditioner can only be used with ``NVECTOR_PARALLEL``. It is not recommended
to use a threaded vector object with SuperLU_MT unless it is the
``NVECTOR_OPENMP`` module, and SuperLU_MT is also compiled with OpenMP.

.. _IDA.Usage.CC.file_access:

Access to library and header files
----------------------------------

At this point, it is assumed that the installation of IDA, following the
procedure described in :numref:`Installation`, has been completed successfully.

Regardless of where the user’s application program resides, its associated
compilation and load commands must make reference to the appropriate locations
for the library and header files required by IDA. The relevant library files are

.. code-block::

  <libdir>/libsundials_ida.<so|a>
  <libdir>/libsundials_nvec*.<so|a>
  <libdir>/libsundials_sunmat*.<so|a>
  <libdir>/libsundials_sunlinsol*.<so|a>
  <libdir>/libsundials_sunnonlinsol*.<so|a>

where the file extension ``.so`` is typically for shared libraries and ``.a``
for static libraries. The relevant header files are located in the
subdirectories

.. code-block::

  <incdir>/ida
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


.. _IDA.Usage.CC.header_sim:

Header files
------------

The calling program must include several header files so that various macros and
data types can be used. The header file that is always required is:

* ``ida/ida.h`` the main header file for IDA, which defines the types and
  various constants, and includes function prototypes. This includes the
  header file for IDALS, ``ida/ida_ls.h``.

Note that ``ida.h`` includes ``sundials_types.h``, which defines the types,
``realtype``, ``sunindextype``, and ``booleantype`` and the constants
``SUNFALSE`` and ``SUNTRUE``.

The calling program must also include an ``N_Vector`` implementation
header file, of the form ``nvector/nvector_*.h`` (see Chapter :numref:`NVectors`
for more information). This file in turn includes the header file
``sundials_nvector.h`` which defines the abstract vector data type.

If using a non-default nonlinear solver object, or when interacting with a
``SUNNonlinearSolver`` object directly, the calling program must also include a
``SUNNonlinearSolver`` implementation header file, of the form
``sunnonlinsol/sunnonlinsol_*.h`` where ``*`` is the name of the nonlinear
solver (see Chapter :numref:`SUNNonlinSol` for more information). This file in
turn includes the header file ``sundials_nonlinearsolver.h`` which defines the
abstract nonlinear linear solver data type.

If using a nonlinear solver that requires the solution of a linear system of the
form :eq:`IDA_DAE_nls` (e.g., the default Newton iteration), the calling program
must also include a ``SUNLinearSolver`` implementation header file, of the from
``sunlinsol/sunlinsol_*.h`` where ``*`` is the name of the linear solver
(see Chapter :numref:`SUNLinSol` for more information).  This file in
turn includes the header file ``sundials_linearsolver.h`` which defines the
abstract linear solver data type.

If the linear solver is matrix-based, the linear solver header will also include
a header file of the from ``sunmatrix/sunmatrix_*.h`` where ``*`` is the name of
the matrix implementation compatible with the linear solver. The matrix header
file provides access to the relevant matrix functions/macros and in turn
includes the header file ``sundials_matrix.h`` which defines the abstract matrix
data type.

Other headers may be needed, according to the choice of preconditioner, etc. For
example, in the example ``idaFoodWeb_kry_p`` (see :cite:p:`ida_ex`),
preconditioning is done with a block-diagonal matrix. For this, even though the
``SUNLINSOL_SPGMR`` linear solver is used, the header
``sundials/sundials_dense.h`` is included for access to the underlying generic
dense matrix arithmetic routines.

.. _IDA.Usage.CC.skeleton_sim:

A skeleton of the user’s main program
-------------------------------------

The following is a skeleton of the user’s main program (or calling program) for
the integration of a DAE IVP. Most of the steps are independent of the
``N_Vector``, ``SUNMatrix``, ``SUNLinearSolver``, and
``SUNNonlinearSolver`` implementations used. For the steps that are not,
refer to Chapters :numref:`NVectors`, :numref:`SUNMatrix`, :numref:`SUNLinSol`,
and :numref:`SUNNonlinSol` for the specific name of the function to be called or
macro to be referenced.

#. **Initialize parallel or multi-threaded environment** (*if appropriate*)

   For example, call ``MPI_Init`` to initialize MPI if used.

#. **Create the SUNDIALS context object**

   Call :c:func:`SUNContext_Create` to allocate the ``SUNContext`` object.

#. **Create the vector of initial values**

   Construct an ``N_Vector`` of initial values using the appropriate functions
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

   Set the vector ``yp0`` of initial conditions for :math:`\dot{y}` similarly.

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

#. **Create nonlinear solver object** (*if appropriate*)

   If using a non-default nonlinear solver, then the desired nonlinear solver
   object must be created by calling the appropriate constructor defined by the
   particular ``SUNNonlinearSolver`` implementation.

   For any of the native SUNDIALS ``SUNNonLinearSolver`` implementations, the
   nonlinear solver object may be created using a call of the form
   ``SUNNonlinearSolver NLS = SUNNonlinSol_***(...);`` where ``***`` is the name
   of the nonlinear solver (see :numref:`SUNNonlinSol` for details).

#. **Create IDA object**

   Call :c:func:`IDACreate` to create the IDA solver object.

#. **Initialize IDA solver**

   Call :c:func:`IDAInit` to provide the initial condition vectors created
   above, set the DAE residual function, and initialize IDA.

#. **Specify integration tolerances**

   Call one of the following functions to set the integration tolerances:

   * :c:func:`IDASStolerances` to specify scalar relative and absolute
     tolerances.

   * :c:func:`IDASVtolerances` to specify a scalar relative tolerance and
     a vector of absolute tolerances.

   * :c:func:`IDAWFtolerances` to specify a function which sets directly the
     weights used in evaluating WRMS vector norms.

   See :numref:`IDA.Usage.CC.callable_fct_sim.toladvice` for general advice on
   selecting tolerances and :numref:`IDA.Usage.CC.callable_fct_sim.unphysical` for
   advice on controlling unphysical values.

#. **Attach the linear solver** (*if appropriate*)

   If a linear solver was created above, initialize the IDALS linear solver
   interface by attaching the linear solver object (and matrix object,
   if applicable) with :c:func:`IDASetLinearSolver`.

#. **Set linear solver optional inputs** (*if appropriate*)

   See :numref:`IDA.Usage.CC.optional_input_ls.Table` for IDALS optional inputs
   and Chapter :numref:`SUNLinSol` for linear solver specific optional inputs.

#. **Attach nonlinear solver module** (*if appropriate*)

   If a nonlinear solver was created above, initialize the IDANLS nonlinear
   solver interface by attaching the nonlinear solver object with
   :c:func:`IDASetNonlinearSolver`.

#. **Set nonlinear solver optional inputs** (*if appropriate*)

   See :numref:`IDA.Usage.CC.optional_input_nls.Table` for IDANLS optional inputs
   and Chapter :numref:`SUNNonlinSol` for nonlinear solver specific optional
   inputs. Note, solver specific optional inputs *must* be called after
   :c:func:`IDASetNonlinearSolver`, otherwise the optional inputs will be
   overridden by IDA defaults.

#. **Specify rootfinding problem** (*optional*)

   Call :c:func:`IDARootInit` to initialize a rootfinding problem to be solved
   during the integration of the ODE system. See
   :numref:`IDA.Usage.CC.optional_input_root.Table` for relevant optional input
   calls.

#. **Set optional inputs**

   Call ``IDASet***`` functions to change any optional inputs that control the
   behavior of IDA from their default values. See
   :numref:`IDA.Usage.CC.optional_input` for details.

#. **Correct initial values** (*optional*)

   Call :c:func:`IDACalcIC` to correct the initial values ``y0`` and ``yp0``
   passed to :c:func:`IDAInit`. See :numref:`IDA.Usage.CC.optional_input_ic.Table`
   for relevant optional input calls.

#. **Advance solution in time**

   For each point at which output is desired, call ``ier = IDASolve(ida_mem,
   tout,  &tret, yret, ypret, itask)``. Here ``itask`` specifies the return
   mode. The vector ``yret`` (which can be the same as the vector ``y0`` above)
   will contain :math:`y(t)`, while the vector ``ypret`` (which can be the same
   as the vector ``yp0`` above) will contain :math:`\dot{y}(t)`.

   See :c:func:`IDASolve` for details.

#. **Get optional outputs**

   Call ``IDAGet***`` functions to obtain optional output. See
   :numref:`IDA.Usage.CC.optional_output` for details.

#. **Deallocate memory**

   Upon completion of the integration call the following, as necessary, to free
   any objects or memory allocated above:

   * Call :c:func:`N_VDestroy` to free vector objects.
   * Call :c:func:`SUNMatDestroy` to free matrix objects.
   * Call :c:func:`SUNLinSolFree` to free linear solvers objects.
   * Call :c:func:`SUNNonlinSolFree` to free nonlinear solvers objects.
   * Call :c:func:`IDAFree` to free the memory allocated by IDA.

#. **Free the SUNContext object**

   Call :c:func:`SUNContext_Free` to free the memory allocated for the ``SUNContext`` object.

#. **Finalize MPI, if used**

   Call ``MPI_Finalize`` to terminate MPI.


.. _IDA.Usage.CC.callable_fct_sim:

User-callable functions
-----------------------

This section describes the IDA functions that are called by the user to setup
and then solve an IVP. Some of these are required.  However, starting with
:numref:`IDA.Usage.CC.optional_input`, the functions listed involve optional
inputs/outputs or restarting, and those paragraphs may be skipped for a casual
use of IDA. In any case, refer to :numref:`IDA.Usage.CC.skeleton_sim` for the
correct order of these calls.

On an error, each user-callable function returns a negative value and sends an
error message to the error handler routine, which prints the message on
``stderr`` by default. However, the user can set a file as error output or can
provide his own error handler function (see
:numref:`IDA.Usage.CC.optional_input.optin_main`).

.. _IDA.Usage.CC.callable_fct_sim.idamalloc:

IDA initialization and deallocation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. c:function:: void* IDACreate(SUNContext sunctx)

   The function ``IDACreate`` instantiates an IDA solver object.

   **Arguments:**
      - ``sunctx`` -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      * ``void*`` pointer the IDA solver object.

.. c:function:: int IDAInit(void* ida_mem, IDAResFn res, realtype t0, N_Vector y0, N_Vector yp0)

   The function ``IDAInit`` provides required problem and solution
   specifications, allocates internal memory, and initializes IDA.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``res`` -- is the function which computes the residual function
        :math:`F(t, y, \dot{y})` for the DAE. For full details see
        :c:type:`IDAResFn`.
      * ``t0`` -- is the initial value of :math:`t`.
      * ``y0`` -- is the initial value of :math:`y`.
      * ``yp0`` -- is the initial value of :math:`\dot{y}`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The call was successful.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
      * ``IDA_MEM_FAIL`` -- A memory allocation request has failed.
      * ``IDA_ILL_INPUT`` -- An input argument to :c:func:`IDAInit` has an illegal
        value.

   **Notes:**
      If an error occurred, :c:func:`IDAInit` also sends an error message to the
      error handler function.

.. c:function:: void IDAFree(void** ida_mem)

   The function ``IDAFree`` frees the pointer allocated by a previous call to
   :c:func:`IDACreate`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.

   **Return value:**
      * ``void``


.. _IDA.Usage.CC.callable_fct_sim.idatolerances:

IDA tolerance specification functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the following three functions must be called to specify the integration
tolerances (or directly specify the weights used in evaluating WRMS vector
norms). Note that this call must be made after the call to :c:func:`IDAInit`.

.. c:function:: int IDASStolerances(void* ida_mem, realtype reltol, realtype abstol)

   The function ``IDASStolerances`` specifies scalar relative and absolute
   tolerances.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``reltol`` -- is the scalar relative error tolerance.
      * ``abstol`` -- is the scalar absolute error tolerance.

   **Return value:**
      * ``IDA_SUCCESS`` -- The call was successful
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
      * ``IDA_NO_MALLOC`` -- The allocation function :c:func:`IDAInit` has not been
        called.
      * ``IDA_ILL_INPUT`` -- One of the input tolerances was negative.

.. c:function:: int IDASVtolerances(void* ida_mem, realtype reltol, N_Vector abstol)

   The function ``IDASVtolerances`` specifies scalar relative tolerance and
   vector absolute tolerances.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``reltol`` -- is the scalar relative error tolerance.
      * ``abstol`` -- is the vector of absolute error tolerances.

   **Return value:**
      * ``IDA_SUCCESS`` -- The call was successful
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
      * ``IDA_NO_MALLOC`` -- The allocation function :c:func:`IDAInit` has not been
        called.
      * ``IDA_ILL_INPUT`` -- The relative error tolerance was negative or the
        absolute tolerance vector had a negative component.

   **Notes:**
      This choice of tolerances is important when the absolute error tolerance
      needs to be different for each component of the state vector :math:`y`.

.. c:function:: int IDAWFtolerances(void* ida_mem, IDAEwtFn efun)

   The function ``IDAWFtolerances`` specifies a user-supplied function ``efun``
   that sets the multiplicative error weights :math:`W_i` for use in the
   weighted RMS norm, which are normally defined by :eq:`IDA_errwt`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
        :c:func:`IDACreate`
      * ``efun`` -- is the function which defines the ``ewt`` vector. For full
        details see :c:type:`IDAEwtFn`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The call was successful
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
      * ``IDA_NO_MALLOC`` -- The allocation function :c:func:`IDAInit` has not been
        called.


.. _IDA.Usage.CC.callable_fct_sim.toladvice:

General advice on choice of tolerances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For many users, the appropriate choices for tolerance values in ``reltol`` and
``abstol`` are a concern. The following pieces of advice are relevant.

#. The scalar relative tolerance ``reltol`` is to be set to control relative
   errors. So ``reltol`` of :math:`10^{-4}` means that errors are controlled to
   .01%. We do not recommend using ``reltol`` larger than :math:`10^{-3}`. On
   the other hand, ``reltol`` should not be so small that it is comparable to
   the unit roundoff of the machine arithmetic (generally around
   :math:`10^{-15}`).

#. The absolute tolerances ``abstol`` (whether scalar or vector) need to be set
   to control absolute errors when any components of the solution vector ``y``
   may be so small that pure relative error control is meaningless. For example,
   if ``y[i]`` starts at some nonzero value, but in time decays to zero, then
   pure relative error control on ``y[i]`` makes no sense (and is overly costly)
   after ``y[i]`` is below some noise level. Then ``abstol`` (if a scalar) or
   ``abstol[i]`` (if a vector) needs to be set to that noise level. If the
   different components have different noise levels, then ``abstol`` should be a
   vector. See the example ``idaRoberts_dns`` in the IDA package, and the
   discussion of it in the IDA Examples document :cite:p:`ida_ex`. In that
   problem, the three components vary betwen 0 and 1, and have different noise
   levels; hence the ``abstol`` vector. It is impossible to give any general
   advice on ``abstol`` values, because the appropriate noise levels are
   completely problem-dependent. The user or modeler hopefully has some idea as
   to what those noise levels are.

#. Finally, it is important to pick all the tolerance values conservatively,
   because they control the error committed on each individual time step. The
   final (global) errors are some sort of accumulation of those per-step errors.
   A good rule of thumb is to reduce the tolerances by a factor of .01 from the
   actual desired limits on errors. So if you want .01% accuracy (globally), a
   good choice is to is a ``reltol`` of :math:`10^{-6}`. But in any case, it is
   a good idea to do a few experiments with the tolerances to see how the
   computed solution values vary as tolerances are reduced.

.. _IDA.Usage.CC.callable_fct_sim.unphysical:

Advice on controlling unphysical negative values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In many applications, some components in the true solution are always positive
or non-negative, though at times very small. In the numerical solution, however,
small negative (hence unphysical) values can then occur. In most cases, these
values are harmless, and simply need to be controlled, not eliminated. The
following pieces of advice are relevant.

#. The way to control the size of unwanted negative computed values is with
   tighter absolute tolerances. Again this requires some knowledge of the noise
   level of these components, which may or may not be different for different
   components. Some experimentation may be needed.

#. If output plots or tables are being generated, and it is important to avoid
   having negative numbers appear there (for the sake of avoiding a long
   explanation of them, if nothing else), then eliminate them, but only in the
   context of the output medium. Then the internal values carried by the solver
   are unaffected. Remember that a small negative value in ``yret`` returned by
   IDA, with magnitude comparable to ``abstol`` or less, is equivalent to zero
   as far as the computation is concerned.

#. The user’s residual function ``res`` should never change a negative value in
   the solution vector ``yy`` to a non-negative value, as a "solution" to this
   problem. This can cause instability. If the ``res`` routine cannot tolerate a
   zero or negative value (e.g., because there is a square root or log of it),
   then the offending value should be changed to zero or a tiny positive number
   in a temporary variable (not in the input ``yy`` vector) for the purposes of
   computing :math:`F(t,y,\dot{y})`.

#. IDA provides the option of enforcing positivity or non-negativity on
   components. Also, such constraints can be enforced by use of the recoverable
   error return feature in the user-supplied residual function. However, because
   these options involve some extra overhead cost, they should only be exercised
   if the use of absolute tolerances to control the computed values is
   unsuccessful.

.. _IDA.Usage.CC.callable_fct_sim.lin_solv_init:

Linear solver interface functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As previously explained, if the nonlinear solver requires the solution of linear
systems of the form :eq:`IDA_DAE_Newtoncorr`, e.g., the default Newton solver, then
the solution of these linear systems is handled with the IDALS linear solver
interface. This interface supports all valid ``SUNLinearSolver`` objects.
Here, a matrix-based ``SUNLinearSolver`` utilizes ``SUNMatrix``
objects to store the Jacobian matrix :math:`J = \dfrac{\partial{F}}{\partial{y}} + \alpha
\dfrac{\partial{F}}{\partial{\dot{y}}}` and factorizations used throughout the solution
process. Conversely, matrix-free ``SUNLinearSolver`` object instead use
iterative methods to solve the linear systems of equations, and only require the
*action* of the Jacobian on a vector, :math:`Jv`.

With most iterative linear solvers, preconditioning can be done on the left
only, on the right only, on both the left and the right, or not at all. The
exceptions to this rule are SPFGMR that supports right preconditioning only and
PCG that performs symmetric preconditioning. However, in IDA only left
preconditioning is supported. For the specification of a preconditioner, see the
iterative linear solver sections in :numref:`IDA.Usage.CC.optional_input` and
:numref:`IDA.Usage.CC.user_fct_sim`. A preconditioner matrix :math:`P` must
approximate the Jacobian :math:`J`, at least crudely.

To attach a generic linear solver to IDA, after the call to :c:func:`IDACreate`
but before any calls to :c:func:`IDASolve`, the user’s program must create the
appropriate ``SUNLinearSolver`` object and call the function
:c:func:`IDASetLinearSolver`. To create the ``SUNLinearSolver`` object,
the user may call one of the SUNDIALS-packaged ``SUNLinearSolver``
constructors via a call of the form

.. code-block:: c

   SUNLinearSolver LS = SUNLinSol_*(...);

Alternately, a user-supplied ``SUNLinearSolver`` object may be created and
used instead. The use of each of the generic linear solvers involves certain
constants, functions and possibly some macros, that are likely to be needed in
the user code. These are available in the corresponding header file associated
with the specific ``SUNMatrix`` or ``SUNLinearSolver`` object in
question, as described in Chapters :numref:`SUNMatrix` and :numref:`SUNLinSol`.

Once this solver object has been constructed, the user should attach it to IDA
via a call to :c:func:`IDASetLinearSolver`. The first argument passed to this
function is the IDA memory pointer returned by :c:func:`IDACreate`; the second
argument is the desired ``SUNLinearSolver`` object to use for solving
systems. The third argument is an optional ``SUNMatrix`` object to
accompany matrix-based ``SUNLinearSolver`` inputs (for matrix-free linear
solvers, the third argument should be ``NULL``). A call to this function
initializes the IDALS linear solver interface, linking it to the main IDA
integrator, and allows the user to specify additional parameters and routines
pertinent to their choice of linear solver.


.. c:function:: int IDASetLinearSolver(void* ida_mem, SUNLinearSolver LS, SUNMatrix J)

   The function ``IDASetLinearSolver`` attaches a ``SUNLinearSolver``
   object ``LS`` and corresponding template Jacobian ``SUNMatrix`` object
   ``J`` (if applicable) to IDA, initializing the IDALS linear solver interface.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``LS`` -- ``SUNLinearSolver`` object to use for solving linear
        systems of the form :eq:`IDA_DAE_Newtoncorr`.
      * ``J`` -- ``SUNMatrix`` object for used as a template for the Jacobian
        or ``NULL`` if not applicable.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The IDALS initialization was successful.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_ILL_INPUT`` -- The IDALS interface is not compatible with the
        ``LS`` or ``J`` input objects or is incompatible with the
        ``N_Vector`` object passed to :c:func:`IDAInit`.
      * ``IDALS_SUNLS_FAIL`` -- A call to the ``LS`` object failed.
      * ``IDALS_MEM_FAIL`` -- A memory allocation request failed.

   **Notes:**
      If ``LS`` is a matrix-based linear solver, then the template Jacobian matrix
      ``J`` will be used in the solve process, so if additional storage is required
      within the ``SUNMatrix`` object (e.g., for factorization of a banded
      matrix), ensure that the input object is allocated with sufficient size (see
      the documentation of the particular ``SUNMatrix`` in Chapter
      :numref:`SUNMatrix` for further information).

   .. warning::

      The previous routines ``IDADlsSetLinearSolver`` and
      ``IDASpilsSetLinearSolver`` are now wrappers for this routine, and may
      still be used for backward-compatibility.  However, these will be
      deprecated in future releases, so we recommend that users transition to
      the new routine name soon.


.. _IDA.Usage.CC.nonlin_solv_init:

Nonlinear solver interface function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default IDA uses the ``SUNNonlinearSolver`` implementation of Newton’s method
(see :numref:`SUNNonlinSol.Newton`). To attach a different nonlinear solver in
IDA, the user’s program must create a ``SUNNonlinearSolver`` object by calling
the appropriate constructor routine. The user must then attach the
``SUNNonlinearSolver`` object to IDA by calling :c:func:`IDASetNonlinearSolver`.

When changing the nonlinear solver in IDA, :c:func:`IDASetNonlinearSolver` must
be called after :c:func:`IDAInit`. If any calls to :c:func:`IDASolve` have been
made, then IDA will need to be reinitialized by calling :c:func:`IDAReInit` to
ensure that the nonlinear solver is initialized correctly before any subsequent
calls to :c:func:`IDASolve`.

The first argument passed to :c:func:`IDASetNonlinearSolver` is the IDA memory
pointer returned by :c:func:`IDACreate` and the second argument is the
``SUNNonlinearSolver`` object to use for solving the nonlinear system
:eq:`IDA_DAE_nls`. A call to this function attaches the nonlinear solver to the main
IDA integrator. We note that at present, the ``SUNNonlinearSolver`` object
*must be of type* ``SUNNONLINEARSOLVER_ROOTFIND``.

.. c:function:: int IDASetNonlinearSolver(void* ida_mem, SUNNonlinearSolver NLS)

   The function ``IDASetNonLinearSolver`` attaches a ``SUNNonlinearSolver``  object (``NLS``) to IDA.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``NLS`` -- ``SUNNonlinearSolver`` object to use for solving nonlinear systems.

   **Return value:**
      * ``IDA_SUCCESS`` -- The nonlinear solver was successfully attached.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- The ``SUNNonlinearSolver`` object is ``NULL`` , does
        not implement the required nonlinear solver operations, is not of the
        correct type, or the residual function, convergence test function, or
        maximum number of nonlinear iterations could not be set.


.. _IDA.Usage.CC.initialcondition:

Initial condition calculation function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:c:func:`IDACalcIC` calculates corrected initial conditions for the DAE system
for certain index-one problems including a class of systems of semi-implicit
form (see :numref:`IDA.Mathematics.ivp_sol` and :cite:p:`BHP:98`). It uses a Newton
iteration combined with a linesearch algorithm. Calling :c:func:`IDACalcIC` is
optional. It is only necessary when the initial conditions do not satisfy the
given system. Thus if ``y0`` and ``yp0`` are known to satisfy
:math:`F(t_0, y_0, \dot{y}_0) = 0`, then a call to :c:func:`IDACalcIC` is
generally *not* necessary.

A call to the function :c:func:`IDACalcIC` must be preceded by successful calls
to :c:func:`IDACreate` and :c:func:`IDAInit` (or :c:func:`IDAReInit`), and by a
successful call to the linear system solver specification function. The call to
:c:func:`IDACalcIC` should precede the call(s) to :c:func:`IDASolve` for the
given problem.

.. c:function:: int IDACalcIC(void* ida_mem, int icopt, realtype tout1)

   The function ``IDACalcIC`` corrects the initial values ``y0`` and ``yp0`` at
   time ``t0``.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``icopt`` -- is one of the following two options for the initial condition
        calculation.

        * ``IDA_YA_YDP_INIT`` directs :c:func:`IDACalcIC` to compute the algebraic
          components of :math:`y` and differential components of :math:`\dot{y}`,
          given the differential components of :math:`y`. This option requires that
          the ``N_Vector id`` was set through :c:func:`IDASetId`, specifying the
          differential and algebraic components.
        * ``IDA_Y_INIT`` directs :c:func:`IDACalcIC` to compute all components of
          :math:`y`, given :math:`\dot{y}`. In this case, ``id`` is not required.

      * ``tout1`` -- is the first value of :math:`t` at which a solution will be
        requested (from :c:func:`IDASolve`). This value is needed here only to
        determine the direction of integration and rough scale in the independent
        variable :math:`t`.

   **Return value:**
      * ``IDA_SUCCESS`` -- :c:func:`IDACalcIC` succeeded.
      * ``IDA_MEM_NULL`` -- The argument ``ida_mem`` was ``NULL``.
      * ``IDA_NO_MALLOC`` -- The allocation function :c:func:`IDAInit` has not been
        called.
      * ``IDA_ILL_INPUT`` -- One of the input arguments was illegal.
      * ``IDA_LSETUP_FAIL`` -- The linear solver's setup function failed in an
        unrecoverable manner.
      * ``IDA_LINIT_FAIL`` -- The linear solver's initialization function failed.
      * ``IDA_LSOLVE_FAIL`` -- The linear solver's solve function failed in an
        unrecoverable manner.
      * ``IDA_BAD_EWT`` -- Some component of the error weight vector is zero
        (illegal), either for the input value of ``y0`` or a corrected value.
      * ``IDA_FIRST_RES_FAIL`` -- The user's residual function returned a
        recoverable error flag on the first call, but :c:func:`IDACalcIC` was
        unable to recover.
      * ``IDA_RES_FAIL`` -- The user's residual function returned a nonrecoverable
        error flag.
      * ``IDA_NO_RECOVERY`` -- The user's residual function, or the linear solver's
        setup or solve function had a recoverable error, but :c:func:`IDACalcIC`
        was unable to recover.
      * ``IDA_CONSTR_FAIL`` -- :c:func:`IDACalcIC` was unable to find a solution
        satisfying the inequality constraints.
      * ``IDA_LINESEARCH_FAIL`` -- The linesearch algorithm failed to find a
        solution with a step larger than ``steptol`` in weighted RMS norm, and
        within the allowed number of backtracks.
      * ``IDA_CONV_FAIL`` -- :c:func:`IDACalcIC` failed to get convergence of the
        Newton iterations.

   **Notes:**
      :c:func:`IDACalcIC` will correct the values of :math:`y(t_0)` and
      :math:`\dot{y}(t_0)` which were specified in the previous call to
      :c:func:`IDAInit` or :c:func:`IDAReInit`. To obtain the corrected values,
      call :c:func:`IDAGetConsistentIC`.


.. _IDA.Usage.CC.idarootinit:

Rootfinding initialization function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While solving the IVP, IDA has the capability to find the roots of a set of
user-defined functions. To activate the root finding algorithm, call the
following function. This is normally called only once, prior to the first call
to :c:func:`IDASolve`, but if the rootfinding problem is to be changed during
the solution, :c:func:`IDARootInit` can also be called prior to a continuation
call to :c:func:`IDASolve`.

.. c:function:: int IDARootInit(void* ida_mem, int nrtfn, IDARootFn g)

   The function ``IDARootInit`` specifies that the roots of a set of functions
   :math:`g_i(t,y)` are to be found while the IVP is being solved.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object
      * ``nrtfn`` -- is the number of root functions
      * ``g`` -- is the function which defines the ``nrtfn`` functions
        :math:`g_i(t,y,\dot{y})` whose roots are sought. See :c:type:`IDARootFn`
        for more details.

   **Return value:**
      * ``IDA_SUCCESS`` -- The call was successful
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
      * ``IDA_MEM_FAIL`` -- A memory allocation failed.
      * ``IDA_ILL_INPUT`` -- The function ``g`` is ``NULL``, but ``nrtfn > 0``.

   **Notes:**
      If a new IVP is to be solved with a call to :c:func:`IDAReInit`, where the
      new IVP has no rootfinding problem but the prior one did, then call
      :c:func:`IDARootInit` with ``nrtfn = 0``.


.. _IDA.Usage.CC.ida:

IDA solver function
~~~~~~~~~~~~~~~~~~~~~

This is the central step in the solution process, the call to perform the
integration of the DAE. The input arguments (``itask``) specifies one of two
modes as to where IDA is to return a solution. These modes are modified if
the user has set a stop time (with :c:func:`IDASetStopTime`) or requested
rootfinding (with :c:func:`IDARootInit`).


.. c:function:: int IDASolve(void * ida_mem, realtype tout, realtype tret, N_Vector yret, N_Vector ypret, int itask)

   The function ``IDASolve`` integrates the DAE over an interval in t.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``tout`` -- the next time at which a computed solution is desired.
      * ``tret`` -- the time reached by the solver output.
      * ``yret`` -- the computed solution vector y.
      * ``ypret`` -- the computed solution vector :math:`\dot{y}`.
      * ``itask`` -- a flag indicating the job of the solver for the next user step

        * ``IDA_NORMAL`` -- the solver will take internal steps until it has
          reached or just passed the user specified ``tout`` parameter. The solver
          then interpolates in order to return approximate values of
          :math:`y(t_{out})` and :math:`\dot{y}(t_{out})`.
        * ``IDA_ONE_STEP`` -- the solver will just take one internal step and
          return the solution at the point reached by that step.

   **Return value:**
      * ``IDA_SUCCESS`` -- The call was successful.
      * ``IDA_TSTOP_RETURN`` -- :c:func:`IDASolve` succeeded by reaching the stop
        point specified through the optional input function
        :c:func:`IDASetStopTime`.
      * ``IDA_ROOT_RETURN`` -- :c:func:`IDASolve` succeeded and found one or more
        roots. In this case, ``tret`` is the location of the root. If ``nrtfn`` >1,
        call :c:func:`IDAGetRootInfo` to see which :math:`g_i` were found to have a
        root.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
      * ``IDA_ILL_INPUT`` -- One of the inputs to :c:func:`IDASolve` was illegal,
        or some other input to the solver was either illegal or missing. The latter
        category includes the following situations:

        * The tolerances have not been set.
        * A component of the error weight vector became zero during internal
          time-stepping.
        * The linear solver initialization function called by the user after
          calling :c:func:`IDACreate` failed to set the linear solver-specific
          ``lsolve`` field in ``ida_mem``.
        * A root of one of the root functions was found both at a point :math:`t`
          and also very near :math:`t`.

        In any case, the user should see the printed error message for details.

      * ``IDA_TOO_MUCH_WORK`` -- The solver took ``mxstep`` internal steps but
        could not reach ``tout``. The default value for ``mxstep`` is
        ``MXSTEP_DEFAULT = 500``.
      * ``IDA_TOO_MUCH_ACC`` -- The solver could not satisfy the accuracy demanded
        by the user for some internal step.
      * ``IDA_ERR_FAIL`` -- Error test failures occurred too many times (``MXNEF =
        10``) during one internal time step or occurred with
        :math:`|h| = h_{\text{min}}`.
      * ``IDA_CONV_FAIL`` -- Convergence test failures occurred too many times
        (``MXNCF = 10``) during one internal time step or occurred with
        :math:`|h| = h_{\text{min}}`.
      * ``IDA_LINIT_FAIL`` -- The linear solver's initialization function failed.
      * ``IDA_LSETUP_FAIL`` -- The linear solver's setup function failed in an
        unrecoverable manner.
      * ``IDA_LSOLVE_FAIL`` -- The linear solver's solve function failed in an
        unrecoverable manner.
      * ``IDA_CONSTR_FAIL`` -- The inequality constraints were violated and the
        solver was unable to recover.
      * ``IDA_REP_RES_ERR`` -- The user's residual function repeatedly returned a
        recoverable error flag, but the solver was unable to recover.
      * ``IDA_RES_FAIL`` -- The user's residual function returned a nonrecoverable
        error flag.
      * ``IDA_RTFUNC_FAIL`` -- The rootfinding function failed.

   **Notes:**
      The vectors ``yret`` and ``ypret`` can occupy the same space as the initial
      condition vectors ``y0`` and ``yp0``, respectively, that were passed to
      :c:func:`IDAInit`.

      In the ``IDA_ONE_STEP`` mode, ``tout`` is used on the first call only, and
      only to get the direction and rough scale of the independent variable.

      If a stop time is enabled (through a call to :c:func:`IDASetStopTime`), then
      :c:func:`IDASolve` returns the solution at ``tstop``. Once the integrator
      returns at a stop time, any future testing for ``tstop`` is disabled (and
      can be reenabled only though a new call to :c:func:`IDASetStopTime`).

      All failure return values are negative and therefore a test ``flag < 0`` will
      trap all :c:func:`IDASolve` failures.

      On any error return in which one or more internal steps were taken by
      :c:func:`IDASolve`, the returned values of ``tret``, ``yret``, and ``ypret``
      correspond to the farthest point reached in the integration.  On all other
      error returns, these values are left unchanged from the previous
      :c:func:`IDASolve` return.


.. _IDA.Usage.CC.optional_input:

Optional input functions
~~~~~~~~~~~~~~~~~~~~~~~~

There are numerous optional input parameters that control the behavior of the
IDA solver. IDA provides functions that can be used to change these optional
input parameters from their default values. The main inputs are divided in the
following categories:

* :numref:`IDA.Usage.CC.optional_input_main.Table` list the main IDA optional
  inputs,

* :numref:`IDA.Usage.CC.optional_input_ls.Table` lists the IDALS linear solver
  interface optional inputs,

* :numref:`IDA.Usage.CC.optional_input_nls.Table` lists the IDANLS nonlinear solver
  interface optional inputs,

* :numref:`IDA.Usage.CC.optional_input_ic.Table` lists the initial condition
  calculation optional inputs, and

* :numref:`IDA.Usage.CC.optional_input_root.Table` lists the rootfinding optional
  inputs.

These optional inputs are described in detail in the remainder of this section.
For the most casual use of IDA, the reader can skip to
:numref:`IDA.Usage.CC.user_fct_sim`.

We note that, on an error return, all of the optional input functions also send
an error message to the error handler function. All error return values are
negative, so the test ``flag < 0`` will catch all errors.

The optional input calls can, unless otherwise noted, be executed in any order.
However, if the user’s program calls either :c:func:`IDASetErrFile` or
:c:func:`IDASetErrHandlerFn`, then that call should appear first, in order to
take effect for any later error message. Finally, a call to an ``IDASet***``
function can, unless otherwise noted, be made at any time from the user’s
calling program and, if successful, takes effect immediately.


.. _IDA.Usage.CC.optional_input.optin_main:

Main solver optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _IDA.Usage.CC.optional_input_main.Table:

.. table:: Optional inputs for IDA

   +--------------------------------------------------------------------+---------------------------------+----------------+
   | **Optional input**                                                 | **Function name**               | **Default**    |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Pointer to an error file                                           | :c:func:`IDASetErrFile`         | ``stderr``     |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Error handler function                                             | :c:func:`IDASetErrHandlerFn`    | internal fn.   |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | User data                                                          | :c:func:`IDASetUserData`        | NULL           |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Maximum order for BDF method                                       | :c:func:`IDASetMaxOrd`          | 5              |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Maximum no. of internal steps before :math:`t_{{\scriptsize out}}` | :c:func:`IDASetMaxNumSteps`     | 500            |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Initial step size                                                  | :c:func:`IDASetInitStep`        | estimated      |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Maximum absolute step size                                         | :c:func:`IDASetMaxStep`         | :math:`\infty` |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Value of :math:`t_{stop}`                                          | :c:func:`IDASetStopTime`        | :math:`\infty` |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Maximum no. of error test failures                                 | :c:func:`IDASetMaxErrTestFails` | 10             |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Suppress alg. vars. from error test                                | :c:func:`IDASetSuppressAlg`     | ``SUNFALSE``   |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Variable types (differential/algebraic)                            | :c:func:`IDASetId`              | NULL           |
   +--------------------------------------------------------------------+---------------------------------+----------------+
   | Inequality constraints on solution                                 | :c:func:`IDASetConstraints`     | NULL           |
   +--------------------------------------------------------------------+---------------------------------+----------------+


.. c:function:: int IDASetErrFile(void * ida_mem, FILE * errfp)

   The function ``IDASetErrFile`` specifies the file pointer where all IDA
   messages should be directed when using the default IDA error handler
   function.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``errfp`` -- pointer to output file.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      The default value for ``errfp`` is ``stderr``.  Passing a value ``NULL``
      disables all future error message output (except for the case in which the
      IDA memory pointer is ``NULL``).  This use of :c:func:`IDASetErrFile` is
      strongly discouraged.

   .. warning::

      If :c:func:`IDASetErrFile` is to be called, it should be called before any
      other optional input functions, in order to take effect for any later
      error message.

.. c:function:: int IDASetErrHandlerFn(void * ida_mem, IDAErrHandlerFn ehfun, void * eh_data)

   The function ``IDASetErrHandlerFn`` specifies the optional user-defined
   function to be used in handling error messages.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``ehfun`` -- is the user's error handler function. See
        :c:type:`IDAErrHandlerFn` for more details.
      * ``eh_data`` -- pointer to user data passed to ``ehfun`` every time it is
        called.

   **Return value:**
      * ``IDA_SUCCESS`` -- The function ``ehfun`` and data pointer ``eh_data`` have
        been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      Error messages indicating that the IDA solver memory is ``NULL`` will always
      be directed to ``stderr``.

.. c:function:: int IDASetUserData(void * ida_mem, void * user_data)

   The function ``IDASetUserData`` attaches a user-defined data pointer to the
   main IDA solver object.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``user_data`` -- pointer to the user data.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      If specified, the pointer to ``user_data`` is passed to all user-supplied
      functions that have it as an argument. Otherwise, a ``NULL`` pointer is
      passed.

   .. warning::

      If ``user_data`` is needed in user linear solver or preconditioner
      functions, the call to :c:func:`IDASetUserData` must be made before the
      call to specify the linear solver.

.. c:function:: int IDASetMaxOrd(void * ida_mem, int maxord)

   The function ``IDASetMaxOrd`` specifies the maximum order of the linear
   multistep method.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``maxord`` -- value of the maximum method order. This must be positive.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- The input value ``maxord`` is :math:`\leq` 0 , or
        larger than the max order value when :c:func:`IDAInit` was called.

   **Notes:**
      The default value is 5. If the input value exceeds 5, the value 5 will be
      used. If called before :c:func:`IDAInit`, ``maxord`` limits the memory
      requirements for the internal IDA memory block and its value cannot be
      increased past the value set when :c:func:`IDAInit` was called.

.. c:function:: int IDASetMaxNumSteps(void * ida_mem, long int mxsteps)

   The function ``IDASetMaxNumSteps`` specifies the maximum number of steps to
   be taken by the solver in its attempt to reach the next output time.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``mxsteps`` -- maximum allowed number of steps.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      Passing ``mxsteps`` = 0 results in IDA using the default value (500).
      Passing ``mxsteps`` < 0 disables the test (not recommended).

.. c:function:: int IDASetInitStep(void * ida_mem, realtype hin)

   The function ``IDASetInitStep`` specifies the initial step size.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``hin`` -- value of the initial step size to be attempted. Pass 0.0 to have
        IDA use the default value.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      By default, IDA estimates the initial step as the solution of
      :math:`\|h \dot{y} \|_{{\scriptsize WRMS}} = 1/2`, with an added restriction
      that :math:`|h| \leq .001|t_{\text{out}} - t_0|`.

.. c:function:: int IDASetMaxStep(void * ida_mem, realtype hmax)

   The function ``IDASetMaxStep`` specifies the maximum absolute value of the
   step size.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``hmax`` -- maximum absolute value of the step size.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- Either ``hmax`` is not positive or it is smaller than
        the minimum allowable step.

   **Notes:**
      Pass ``hmax = 0`` to obtain the default value :math:`\infty`.

.. c:function:: int IDASetStopTime(void * ida_mem, realtype tstop)

   The function ``IDASetStopTime`` specifies the value of the independent
   variable :math:`t` past which the solution is not to proceed.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``tstop`` -- value of the independent variable past which the solution
        should not proceed.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- The value of ``tstop`` is not beyond the current
        :math:`t` value, :math:`t_n`.

   **Notes:**
      The default, if this routine is not called, is that no stop time is imposed.
      Once the integrator returns at a stop time, any future testing for ``tstop``
      is disabled (and can be reenabled only though a new call to
      :c:func:`IDASetStopTime`).

.. c:function:: int IDASetMaxErrTestFails(void * ida_mem, int maxnef)

   The function ``IDASetMaxErrTestFails`` specifies the maximum number of error
   test failures in attempting one step.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``maxnef`` -- maximum number of error test failures allowed on one step
        (>0).

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      The default value is 10.

.. c:function:: int IDASetSuppressAlg(void * ida_mem, booleantype suppressalg)

   The function ``IDASetSuppressAlg`` indicates whether or not to suppress
   algebraic variables in the local error test.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``suppressalg`` -- indicates whether to suppress (``SUNTRUE``) or include
        (``SUNFALSE``) the algebraic variables in the local error test.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      The default value is ``SUNFALSE``.  If ``suppressalg = SUNTRUE`` is selected,
      then the ``id`` vector must be set (through :c:func:`IDASetId`) to specify
      the algebraic components.  In general, the use of this option (with
      ``suppressalg = SUNTRUE``) is *discouraged* when solving DAE systems of index
      1, whereas it is generally *encouraged* for systems of index 2 or more. See
      pp. 146-147 of :cite:p:`BCP:96` for more on this issue.

.. c:function:: int IDASetId(void * ida_mem, N_Vector id)

   The function ``IDASetId`` specifies algebraic/differential components in the
   :math:`y` vector.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``id`` -- a vector of values identifying the components of :math:`y` as
        differential or algebraic variables. A value of 1.0 indicates a
        differential variable, while 0.0 indicates an algebraic variable.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      The vector ``id`` is required if the algebraic variables are to be suppressed
      from the local error test (see :c:func:`IDASetSuppressAlg`) or if
      :c:func:`IDACalcIC` is to be called with ``icopt`` = ``IDA_YA_YDP_INIT``.

.. c:function:: int IDASetConstraints(void * ida_mem, N_Vector constraints)

   The function ``IDASetConstraints`` specifies a vector defining inequality
   constraints for each component of the solution vector :math:`y`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``constraints`` -- vector of constraint flags.

        * If ``constraints[i] = 0``,  no constraint is imposed on :math:`y_i`.
        * If ``constraints[i] = 1``,  :math:`y_i` will be constrained to be :math:`y_i \ge 0.0`.
        * If ``constraints[i] = -1``, :math:`y_i` will be constrained to be :math:`y_i \le 0.0`.
        * If ``constraints[i] = 2``,  :math:`y_i` will be constrained to be :math:`y_i > 0.0`.
        * If ``constraints[i] = -2``, :math:`y_i` will be constrained to be :math:`y_i < 0.0`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- The constraints vector contains illegal values.

   **Notes:**
      The presence of a non-``NULL`` constraints vector that is not 0.0 in all
      components will cause constraint checking to be performed.  However, a call
      with 0.0 in all components of constraints vector will result in an illegal
      input return. A ``NULL`` input will disable constraint checking.


.. _IDA.Usage.CC.optional_inputs.optin_ls:

Linear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _IDA.Usage.CC.optional_input_ls.Table:

.. table:: Optional inputs for the IDALS linear solver interface

   +-------------------------------------------------+---------------------------------------+---------------+
   | **Optional input**                              | **Function name**                     | **Default**   |
   +-------------------------------------------------+---------------------------------------+---------------+
   | Jacobian function                               | :c:func:`IDASetJacFn`                 | DQ            |
   +-------------------------------------------------+---------------------------------------+---------------+
   | Enable or disable linear solution scaling       | :c:func:`IDASetLinearSolutionScaling` | on            |
   +-------------------------------------------------+---------------------------------------+---------------+
   | Jacobian-times-vector function                  | :c:func:`IDASetJacTimes`              | NULL, DQ      |
   +-------------------------------------------------+---------------------------------------+---------------+
   | Preconditioner functions                        | :c:func:`IDASetPreconditioner`        | NULL, NULL    |
   +-------------------------------------------------+---------------------------------------+---------------+
   | Ratio between linear and nonlinear tolerances   | :c:func:`IDASetEpsLin`                | 0.05          |
   +-------------------------------------------------+---------------------------------------+---------------+
   | Increment factor used in DQ :math:`Jv` approx.  | :c:func:`IDASetIncrementFactor`       | 1.0           |
   +-------------------------------------------------+---------------------------------------+---------------+
   | Jacobian-times-vector DQ Res function           | :c:func:`IDASetJacTimesResFn`         | NULL          |
   +-------------------------------------------------+---------------------------------------+---------------+
   | Newton linear solve tolerance conversion factor | :c:func:`IDASetLSNormFactor`          | vector length |
   +-------------------------------------------------+---------------------------------------+---------------+

The mathematical explanation of the linear solver methods available to IDA is
provided in :numref:`IDA.Mathematics.ivp_sol`. We group the user-callable routines
into four categories: general routines concerning the overall IDALS linear
solver interface, optional inputs for matrix-based linear solvers, optional
inputs for matrix-free linear solvers, and optional inputs for iterative linear
solvers. We note that the matrix-based and matrix-free groups are mutually
exclusive, whereas the “iterative” tag can apply to either case.

When using matrix-based linear solver modules, the IDALS solver interface needs
a function to compute an approximation to the Jacobian matrix
:math:`J(t,y,\dot{y})`. This function must be of type :c:type:`IDALsJacFn`. The
user can supply a Jacobian function or, if using a
:ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` matrix
:math:`J`, can use the default internal difference quotient approximation that
comes with the IDALS interface. To specify a user-supplied Jacobian function
``jac``, IDALS provides the function :c:func:`IDASetJacFn`. The IDALS interface
passes the pointer ``user_data`` to the Jacobian function. This allows the user
to create an arbitrary structure with relevant problem data and access it during
the execution of the user-supplied Jacobian function, without using global data
in the program. The pointer ``user_data`` may be specified through
:c:func:`IDASetUserData`.

.. c:function:: int IDASetJacFn(void * ida_mem, IDALsJacFn jac)

   The function ``IDASetJacFn`` specifies the Jacobian approximation function to
   be used for a matrix-based solver within the IDALS interface.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``jac`` -- user-defined Jacobian approximation function. See
        :c:type:`IDALsJacFn` for more details.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver interface has not been
        initialized.

   **Notes:**
      This function must be called after the IDALS linear solver interface has been
      initialized through a call to :c:func:`IDASetLinearSolver`.  By default,
      IDALS uses an internal difference quotient function for the
      :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` and
      :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` modules.  If ``NULL`` is passed to
      ``jac``, this default function is used.
      An error will occur if no ``jac`` is supplied when using other matrix types.

   .. warning::

      The previous routine ``IDADlsSetJacFn`` is now a wrapper for this routine,
      and may still be used for backward-compatibility.  However, this will be
      deprecated in future releases, so we recommend that users transition to
      the new routine name soon.


When using a matrix-based linear solver the matrix information will be updated
infrequently to reduce matrix construction and, with direct solvers,
factorization costs. As a result the value of :math:`\alpha` may not be current
and a scaling factor is applied to the solution of the linear system to account
for the lagged value of :math:`\alpha`. See :numref:`SUNLinSol.IDA.Lagged` for
more details. The function :c:func:`IDASetLinearSolutionScaling` can be used to
disable this scaling when necessary, e.g., when providing a custom linear solver
that updates the matrix using the current :math:`\alpha` as part of the solve.

.. c:function:: int IDASetLinearSolutionScaling(void * ida_mem, booleantype onoff)

   The function ``IDASetLinearSolutionScaling`` enables or disables scaling the
   linear system solution to account for a change in :math:`\alpha` in the
   linear system. For more details see :numref:`SUNLinSol.IDA.Lagged`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``onoff`` -- flag to enable (``SUNTRUE``) or disable (``SUNFALSE``) scaling

   **Return value:**
      * ``IDALS_SUCCESS`` -- The flag value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver interface has not been
        initialized.
      * ``IDALS_ILL_INPUT`` -- The attached linear solver is not matrix-based.

   **Notes:**
      This function must be called after the IDALS linear solver interface has been
      initialized through a call to :c:func:`IDASetLinearSolver`.  By default
      scaling is enabled with matrix-based linear solvers.


When using matrix-free linear solver modules, the IDALS solver interface
requires a function to compute an approximation to the product between the
Jacobian matrix :math:`J(t,y,\dot{y})` and a vector :math:`v`. The user can
supply a Jacobian-times-vector approximation function, or use the default
internal difference quotient function that comes with the IDALS solver
interface.

A user-defined Jacobian-vector product function must be of type
:c:type:`IDALsJacTimesVecFn` and can be specified through a call to
:c:func:`IDASetJacTimes`. The evaluation and processing of any Jacobian-related
data needed by the user’s Jacobian-vector product function may be done in the
optional user-supplied function ``jtsetup`` (see
:numref:`IDA.Usage.CC.user_fct_sim.jtsetupFn` for specification details). The
pointer ``user_data`` received through :c:func:`IDASetUserData` (or a pointer to
``NULL`` if ``user_data`` was not specified) is passed to the Jacobian-vector
product setup and product functions, ``jtsetup`` and ``jtimes``, each time they
are called.  This allows the user to create an arbitrary structure with relevant
problem data and access it during the execution of the user-supplied functions
without using global data in the program.

.. c:function:: int IDASetJacTimes(void * ida_mem, IDALsJacTimesSetupFn jsetup, IDALsJacTimesVecFn jtimes)

   The function ``IDASetJacTimes`` specifies the Jacobian-vector product setup
   and product functions.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``jtsetup`` -- user-defined function to set up the Jacobian-vector
        product. See :c:type:`IDALsJacTimesSetupFn` for more details. Pass ``NULL``
        if no setup is necessary.
      * ``jtimes`` -- user-defined Jacobian-vector product function. See
        :c:type:`IDALsJacTimesVecFn` for more details.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
      * ``IDALS_SUNLS_FAIL`` -- An error occurred when setting up the system
        matrix-times-vector routines in the ``SUNLinearSolver`` object used by the
        IDALS interface.

   **Notes:**
      The default is to use an internal finite difference quotient for ``jtimes``
      and to omit ``jtsetup``.  If ``NULL`` is passed to ``jtimes``, these defaults
      are used.  A user may specify non-``NULL`` ``jtimes`` and ``NULL``
      ``jtsetup`` inputs.  This function must be called after the IDALS linear
      solver interface has been initialized through a call to
      :c:func:`IDASetLinearSolver`.

   .. warning::

      The previous routine ``IDASpilsSetJacTimes`` is now a wrapper for this
      routine, and may still be used for backward-compatibility.  However, this
      will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.


When using the default difference-quotient approximation to the Jacobian-vector
product, the user may specify the factor to use in setting increments for the
finite-difference approximation, via a call to :c:func:`IDASetIncrementFactor`.

.. c:function:: int IDASetIncrementFactor(void * ida_mem, realtype dqincfac)

   The function ``IDASetIncrementFactor`` specifies the increment factor to be
   used in the difference-quotient approximation to the product :math:`Jv`.
   Specifically, :math:`Jv` is approximated via the formula

   .. math::
      Jv = \frac{1}\sigma\left[F(t,\tilde{y},\tilde{\dot{y}}) - F(t,y,\dot{y})\right],

   where :math:`\tilde{y} = y + \sigma v`,
   :math:`\tilde{\dot{y}} = \dot{y} + c_j  \sigma v`, :math:`c_j` is a BDF
   parameter proportional to the step size,
   :math:`\sigma = \mathtt{dqincfac} \sqrt{N}`, and :math:`N` is the number of
   equations in the DAE system.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``dqincfac`` -- user-specified increment factor positive.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
      * ``IDALS_ILL_INPUT`` -- The specified value of ``dqincfac`` is
        :math:`\le 0`.

   **Notes:**
      The default value is 1.0.  This function must be called after the IDALS
      linear solver interface has been initialized through a call to
      :c:func:`IDASetLinearSolver`.

   .. warning::

      The previous routine :c:func:`IDASpilsSetIncrementFactor` is now a wrapper
      for this routine, and may still be used for backward-compatibility.
      However, this will be deprecated in future releases, so we recommend that
      users transition to the new routine name soon.


Additionally, when using the internal difference quotient, the user may also
optionally supply an alternative residual function for use in the
Jacobian-vector product approximation by calling
:c:func:`IDASetJacTimesResFn`. The alternative residual function should compute
a suitable (and differentiable) approximation to the residual function provided
to :c:func:`IDAInit`. For example, as done in :cite:p:`dorr2010numerical` for an
ODE in explicit form, the alternative function may use lagged values when
evaluating a nonlinearity to avoid differencing a potentially non-differentiable
factor.

.. c:function:: int IDASetJacTimesResFn(void * ida_mem, IDAResFn jtimesResFn)

   The function ``IDASetJacTimesResFn`` specifies an alternative DAE residual
   function for use in the internal Jacobian-vector product difference quotient
   approximation.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``jtimesResFn`` -- is the function which computes the alternative DAE
        residual function to use in Jacobian-vector product difference quotient
        approximations. See :c:type:`IDAResFn` for more details.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
      * ``IDALS_ILL_INPUT`` -- The internal difference quotient approximation is
        disabled.

   **Notes:**
      The default is to use the residual function provided to :c:func:`IDAInit` in
      the internal difference quotient. If the input resudual function is ``NULL``,
      the default is used.  This function must be called after the IDALS linear
      solver interface has been initialized through a call to
      :c:func:`IDASetLinearSolver`.


When using an iterative linear solver, the user may supply a preconditioning
operator to aid in solution of the system. This operator consists of two
user-supplied functions, ``psetup`` and ``psolve``, that are supplied to IDA
using the function :c:func:`IDASetPreconditioner`. The ``psetup`` function
supplied to this routine should handle evaluation and preprocessing of any
Jacobian data needed by the user’s preconditioner solve function,
``psolve``. Both of these functions are fully specified in
:numref:`IDA.Usage.CC.user_fct_sim.psolveFn` and
:numref:`IDA.Usage.CC.user_fct_sim.precondFn`).  The user data pointer received
through :c:func:`IDASetUserData` (or ``NULL`` if a user data pointer was not
specified) is passed to the ``psetup`` and ``psolve`` functions. This allows the
user to create an arbitrary structure with relevant problem data and access it
during the execution of the user-supplied preconditioner functions without using
global data in the program.

.. c:function:: int IDASetPreconditioner(void * ida_mem, IDALsPrecSetupFn psetup, IDALsPrecSolveFn psolve)

   The function ``IDASetPreconditioner`` specifies the preconditioner setup and
   solve functions.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``psetup`` -- user-defined function to set up the preconditioner. See
        :c:type:`IDALsPrecSetupFn` for more details. Pass ``NULL`` if no setup is
        necessary.
      * ``psolve`` -- user-defined preconditioner solve function. See
        :c:type:`IDALsPrecSolveFn` for more details.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional values have been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
      * ``IDALS_SUNLS_FAIL`` -- An error occurred when setting up preconditioning
        in the ``SUNLinearSolver`` object used by the IDALS interface.

   **Notes:**
      The default is ``NULL`` for both arguments (i.e., no preconditioning).  This
      function must be called after the IDALS linear solver interface has been
      initialized through a call to :c:func:`IDASetLinearSolver`.

   .. warning::

      The previous routine ``IDASpilsSetPreconditioner`` is now a wrapper for
      this routine, and may still be used for backward-compatibility.  However,
      this will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.


Also, as described in :numref:`IDA.Mathematics.ivp_sol`, the IDALS interface
requires that iterative linear solvers stop when the norm of the preconditioned
residual satisfies

.. math::
   \|r\| \le \frac{\epsilon_L \epsilon}{10}

where :math:`\epsilon` is the nonlinear solver tolerance, and the default
:math:`\epsilon_L = 0.05`; this value may be modified by the user through the
:c:func:`IDASetEpsLin` function.

.. c:function:: int IDASetEpsLin(void * ida_mem, realtype eplifac)

   The function ``IDASetEpsLin`` specifies the factor by which the Krylov linear
   solver's convergence test constant is reduced from the nonlinear iteration
   test constant.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``eplifac`` -- linear convergence safety factor :math:`\geq 0.0`.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.
      * ``IDALS_ILL_INPUT`` -- The factor ``eplifac`` is negative.

   **Notes:**
      The default value is :math:`0.05`.  This function must be called after the
      IDALS linear solver interface has been initialized through a call to
      :c:func:`IDASetLinearSolver`.  If ``eplifac`` :math:`= 0.0` is passed, the
      default value is used.

   .. warning::

      The previous routine ``IDASpilsSetEpsLin`` is now a wrapper for this
      routine, and may still be used for backward-compatibility.  However, this
      will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.

.. c:function:: int IDASetLSNormFactor(void * ida_mem, realtype nrmfac)

   The function ``IDASetLSNormFactor`` specifies the factor to use when
   converting from the integrator tolerance (WRMS norm) to the linear solver
   tolerance (L2 norm) for Newton linear system solves e.g.,
   ``tol_L2 = fac * tol_WRMS``.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nrmfac`` -- the norm conversion factor.

        * If ``nrmfac > 0``, the provided value is used.
        * If ``nrmfac = 0`` then the conversion factor is computed using the
          vector length i.e., ``nrmfac = N_VGetLength(y)`` (*default*).
        * If ``nrmfac < 0`` then the conversion factor is computed using the
          vector dot product ``nrmfac = N_VDotProd(v,v)`` where all the entries of
          ``v`` are one.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      This function must be called after the IDALS linear solver interface has been
      initialized through a call to :c:func:`IDASetLinearSolver`.  Prior to the
      introduction of :c:func:`N_VGetLength` in SUNDIALS v5.0.0 (IDA v5.0.0) the
      value of ``nrmfac`` was computed using :c:func:`N_VDotProd` i.e., the
      ``nrmfac < 0`` case.


.. _IDA.Usage.CC.optional_inputs.optin_nls:

Nonlinear solver interface optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _IDA.Usage.CC.optional_input_nls.Table:

.. table:: Optional inputs for the IDANLS nonlinear solver interface

   +----------------------------------------------------+--------------------------------+-------------+
   | **Optional input**                                 | **Function name**              | **Default** |
   +----------------------------------------------------+--------------------------------+-------------+
   | Maximum no. of nonlinear iterations                | :c:func:`IDASetMaxNonlinIters` | 4           |
   +----------------------------------------------------+--------------------------------+-------------+
   | Maximum no. of convergence failures                | :c:func:`IDASetMaxConvFails`   | 10          |
   +----------------------------------------------------+--------------------------------+-------------+
   | Coeff. in the nonlinear convergence test           | :c:func:`IDASetNonlinConvCoef` | 0.33        |
   +----------------------------------------------------+--------------------------------+-------------+
   | Residual function for nonlinear system evaluations | :c:func:`IDASetNlsResFn`       | ``NULL``    |
   +----------------------------------------------------+--------------------------------+-------------+

The following functions can be called to set optional inputs controlling the
nonlinear solver.

.. c:function:: int IDASetMaxNonlinIters(void * ida_mem, int maxcor)

   The function ``IDASetMaxNonlinIters`` specifies the maximum number of
   nonlinear solver iterations in one solve attempt.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``maxcor`` -- maximum number of nonlinear solver iterations allowed in one
        solve attempt (>0).

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_MEM_FAIL`` -- The ``SUNNonlinearSolver`` object is ``NULL``.

   **Notes:**
      The default value is 4.


.. c:function:: int IDASetMaxConvFails(void * ida_mem, int maxncf)

   The function ``IDASetMaxConvFails`` specifies the maximum number of nonlinear
   solver convergence failures in one step.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``maxncf`` -- maximum number of allowable nonlinear solver convergence
        failures in one step (>0).

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      The default value is 10.


.. c:function:: int IDASetNonlinConvCoef(void * ida_mem, realtype nlscoef)

   The function ``IDASetNonlinConvCoef`` specifies the safety factor in the
   nonlinear convergence test; see :eq:`IDA_DAE_nls_test`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nlscoef`` -- coefficient in nonlinear convergence test (>0.0).

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- The value of ``nlscoef`` is :math:`\leq 0.0`.

   **Notes:**
      The default value is 0.33.


.. c:function:: int IDASetNlsResFn(void * ida_mem, IDAResFn res)

   The function ``IDASetNlsResFn`` specifies an alternative residual function
   for use in nonlinear system function evaluations.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``res`` -- the alternative function which computes the DAE residual
        function :math:`F(t, y, \dot{y})`. See :c:type:`IDAResFn` for more details.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional function has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      The default is to use the residual function provided to :c:func:`IDAInit`
      in nonlinear system function evaluations. If the input residual function
      is ``NULL``, the default is used.

      When using a non-default nonlinear solver, this function must be called
      after :c:func:`IDASetNonlinearSolver`.


.. _IDA.Usage.CC.optional_input.optin_iccalc:

Initial condition calculation optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _IDA.Usage.CC.optional_input_ic.Table:

.. table:: Optional inputs for IDA initial condition calculation

   +---------------------------------------------+----------------------------------+------------------------+
   | **Optional input**                          | **Function name**                | **Default**            |
   +---------------------------------------------+----------------------------------+------------------------+
   | Coeff. in the nonlinear convergence test    | :c:func:`IDASetNonlinConvCoefIC` | 0.0033                 |
   +---------------------------------------------+----------------------------------+------------------------+
   | Maximum no. of steps                        | :c:func:`IDASetMaxNumStepsIC`    | 5                      |
   +---------------------------------------------+----------------------------------+------------------------+
   | Maximum no. of Jacobian/precond. evals.     | :c:func:`IDASetMaxNumJacsIC`     | 4                      |
   +---------------------------------------------+----------------------------------+------------------------+
   | Maximum no. of Newton iterations            | :c:func:`IDASetMaxNumItersIC`    | 10                     |
   +---------------------------------------------+----------------------------------+------------------------+
   | Max. linesearch backtracks per Newton iter. | :c:func:`IDASetMaxBacksIC`       | 100                    |
   +---------------------------------------------+----------------------------------+------------------------+
   | Turn off linesearch                         | :c:func:`IDASetLineSearchOffIC`  | ``SUNFALSE``           |
   +---------------------------------------------+----------------------------------+------------------------+
   | Lower bound on Newton step                  | :c:func:`IDASetStepToleranceIC`  | uround\ :math:`^{2/3}` |
   +---------------------------------------------+----------------------------------+------------------------+

The following functions can be called just prior to calling :c:func:`IDACalcIC`
to set optional inputs controlling the initial condition calculation.

.. c:function:: int IDASetNonlinConvCoefIC(void * ida_mem, realtype epiccon)

   The function ``IDASetNonlinConvCoefIC`` specifies the positive constant in
   the Newton iteration convergence test within the initial condition
   calculation.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``epiccon`` -- coefficient in the Newton convergence test :math:`(>0)`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- The ``epiccon`` factor is :math:`\leq 0.0`.

   **Notes:**
      The default value is :math:`0.01 \cdot 0.33`.  This test uses a weighted RMS
      norm (with weights defined by the tolerances).  For new initial value vectors
      :math:`y` and :math:`\dot{y}` to be accepted, the norm of
      :math:`J^{-1}F(t_0, y, \dot{y})` must be :math:`\leq \mathtt{epiccon}`, where
      :math:`J` is the system Jacobian.

.. c:function:: int IDASetMaxNumStepsIC(void * ida_mem, int maxnh)

   The function ``IDASetMaxNumStepsIC`` specifies the maximum number of steps
   allowed when ``icopt = IDA_YA_YDP_INIT`` in :c:func:`IDACalcIC`, where
   :math:`h` appears in the system Jacobian, :math:`J = \dfrac{\partial F}{\partial
   y} + \left(\dfrac1h\right)\dfrac{\partial F}{\partial \dot{y}}`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``maxnh`` -- maximum allowed number of values for :math:`h`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- ``maxnh`` is non-positive.

   **Notes:**
      The default value is :math:`5`.

.. c:function:: int IDASetMaxNumJacsIC(void * ida_mem, int maxnj)

   The function ``IDASetMaxNumJacsIC`` specifies the maximum number of the
   approximate Jacobian or preconditioner evaluations allowed when the Newton
   iteration appears to be slowly converging.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``maxnj`` -- maximum allowed number of Jacobian or preconditioner
        evaluations.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- ``maxnj`` is non-positive.

   **Notes:**
      The default value is :math:`4`.

.. c:function:: int IDASetMaxNumItersIC(void * ida_mem, int maxnit)

   The function ``IDASetMaxNumItersIC`` specifies the maximum number of Newton
   iterations allowed in any one attempt to solve the initial conditions
   calculation problem.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``maxnit`` -- maximum number of Newton iterations.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- ``maxnit`` is non-positive.

   **Notes:**
      The default value is :math:`10`.

.. c:function:: int IDASetMaxBacksIC(void * ida_mem, int maxbacks)

   The function ``IDASetMaxBacksIC`` specifies the maximum number of linesearch
   backtracks allowed in any Newton iteration, when solving the initial
   conditions calculation problem.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``maxbacks`` -- maximum number of linesearch backtracks per Newton step.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- ``maxbacks`` is non-positive.

   **Notes:**
      The default value is :math:`100`.

.. c:function:: int IDASetLineSearchOffIC(void * ida_mem, booleantype lsoff)

   The function ``IDASetLineSearchOffIC`` specifies whether to turn on or off
   the linesearch algorithm.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``lsoff`` -- a flag to turn off (``SUNTRUE``) or keep (``SUNFALSE``) the
        linesearch algorithm.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**

   The default value is ``SUNFALSE``.

.. c:function:: int IDASetStepToleranceIC(void * ida_mem, int steptol)

   The function ``IDASetStepToleranceIC`` specifies a positive lower bound on
   the Newton step.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``steptol`` -- Minimum allowed WRMS-norm of the Newton step :math:`(> 0.0)`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- The ``steptol`` tolerance is :math:`\leq 0.0`.

   **Notes:**
      The default value is :math:`(\text{unit roundoff})^{2/3}`.



.. _IDA.Usage.CC.optional_input.optin_root:

Rootfinding optional input functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _IDA.Usage.CC.optional_input_root.Table:

.. table:: Optional inputs for IDA rootfinding

   +------------------------------+------------------------------------+-------------+
   | **Optional input**           | **Function name**                  | **Default** |
   +------------------------------+------------------------------------+-------------+
   | Direction of zero-crossing   | :c:func:`IDASetRootDirection`      | both        |
   +------------------------------+------------------------------------+-------------+
   | Disable rootfinding warnings | :c:func:`IDASetNoInactiveRootWarn` | none        |
   +------------------------------+------------------------------------+-------------+

The following functions can be called to set optional inputs to control the
rootfinding algorithm.


.. c:function:: int IDASetRootDirection(void * ida_mem, int * rootdir)

   The function ``IDASetRootDirection`` specifies the direction of
   zero-crossings to be located and returned to the user.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``rootdir`` -- state array of length ``nrtfn`` , the number of root
        functions :math:`g_i` , as specified in the call to the function
        :c:func:`IDARootInit`.

        * A value of :math:`0` for ``rootdir[i]`` indicates that crossing in either
          direction should be reported for :math:`g_i`.

        * A value of :math:`+1` or :math:`-1` for ``rootdir[i]`` indicates that the
          solver should report only zero-crossings where :math:`g_i` is increasing
          or decreasing, respectively.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_ILL_INPUT`` -- rootfinding has not been activated through a call to
        :c:func:`IDARootInit`.

   **Notes:**
      The default behavior is to locate both zero-crossing directions.

.. c:function:: int IDASetNoInactiveRootWarn(void * ida_mem)

   The function ``IDASetNoInactiveRootWarn`` disables issuing a warning if some
   root function appears to be identically zero at the beginning of the
   integration.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      IDA will not report the initial conditions as a possible zero-crossing
      (assuming that one or more components :math:`g_i` are zero at the initial
      time).  However, if it appears that some :math:`g_i` is identically zero at
      the initial time (i.e., :math:`g_i` is zero at the initial time and after the
      first step), IDA will issue a warning which can be disabled with this
      optional input function.


.. _IDA.Usage.CC.optional_dky:

Interpolated output function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An optional function :c:func:`IDAGetDky` is available to obtain additional
output values. This function must be called after a successful return from
:c:func:`IDASolve` and provides interpolated values of :math:`y` or its
derivatives of order up to the last internal order used for any value of
:math:`t` in the last internal step taken by IDA.

.. c:function:: int IDAGetDky(void * ida_mem, realtype t, int k, N_Vector dky)

   The function ``IDAGetDky`` computes the interpolated values of the
   :math:`k^{th}` derivative of :math:`y` for any value of :math:`t` in the last
   internal step taken by IDA.  The value of :math:`k` must be non-negative and
   smaller than the last internal order used. A value of :math:`0` for :math:`k`
   means that the :math:`y` is interpolated.  The value of :math:`t` must
   satisfy :math:`t_n - h_u \le t \le t_n`, where :math:`t_n` denotes the
   current internal time reached, and :math:`h_u` is the last internal step size
   used successfully.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``t`` -- time at which to interpolate.
      * ``k`` -- integer specifying the order of the derivative of :math:`y`
        wanted.
      * ``dky`` -- vector containing the interpolated :math:`k^{th}` derivative of
        :math:`y(t)`.

   **Return value:**
      * ``IDA_SUCCESS`` -- ``IDAGetDky`` succeeded.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` argument was ``NULL``.
      * ``IDA_BAD_T`` -- ``t`` is not in the interval :math:`[t_n - h_u , t_n]`.
      * ``IDA_BAD_K`` -- ``k`` is not one of
        :math:`{0, 1, \ldots, k_{\text{last}}}`.
      * ``IDA_BAD_DKY`` -- ``dky`` is ``NULL``.

   **Notes:**
      It is only legal to call the function :c:func:`IDAGetDky` after a successful
      return from :c:func:`IDASolve`. Functions :c:func:`IDAGetCurrentTime`,
      :c:func:`IDAGetLastStep` and :c:func:`IDAGetLastOrder` can be used to access
      :math:`t_n`, :math:`h_u`, and :math:`k_{\text{last}}`.



.. _IDA.Usage.CC.optional_output:

Optional output functions
~~~~~~~~~~~~~~~~~~~~~~~~~

IDA provides an extensive list of functions that can be used to obtain solver
performance information. :numref:`IDA.Usage.CC.optional_output.Table` lists
all optional output functions in IDA, which are then described in detail in the
remainder of this section.

Some of the optional outputs, especially the various counters, can be very
useful in determining how successful the IDA solver is in doing its job. For
example, the counters ``nsteps`` and ``nrevals`` provide a rough measure of the
overall cost of a given run, and can be compared among runs with differing input
options to suggest which set of options is most efficient. The ratio
``nniters/nsteps`` measures the performance of the nonlinear solver in solving
the nonlinear systems at each time step; typical values for this range from 1.1
to 1.8. The ratio ``njevals/nniters`` (in the case of a matrix-based linear
solver), and the ratio ``npevals/nniters`` (in the case of an iterative linear
solver) measure the overall degree of nonlinearity in these systems, and also
the quality of the approximate Jacobian or preconditioner being used. Thus, for
example, ``njevals/nniters`` can indicate if a user-supplied Jacobian is
inaccurate, if this ratio is larger than for the case of the corresponding
internal Jacobian. The ratio ``nliters/nniters`` measures the performance of the
Krylov iterative linear solver, and thus (indirectly) the quality of the
preconditioner.

.. _IDA.Usage.CC.optional_output.Table:
.. table:: Optional outputs for IDA, IDALS, and IDANLS
  :align: center

  +--------------------------------------------------------------------+---------------------------------------+
  | **Optional output**                                                | **Function name**                     |
  +====================================================================+=======================================+
  | Size of IDA real and integer workspace                             | :c:func:`IDAGetWorkSpace`             |
  +--------------------------------------------------------------------+---------------------------------------+
  | Cumulative number of internal steps                                | :c:func:`IDAGetNumSteps`              |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of calls to residual function                                  | :c:func:`IDAGetNumResEvals`           |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of calls to linear solver setup function                       | :c:func:`IDAGetNumLinSolvSetups`      |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of local error test failures that have occurred                | :c:func:`IDAGetNumErrTestFails`       |
  +--------------------------------------------------------------------+---------------------------------------+
  | Order used during the last step                                    | :c:func:`IDAGetLastOrder`             |
  +--------------------------------------------------------------------+---------------------------------------+
  | Order to be attempted on the next step                             | :c:func:`IDAGetCurrentOrder`          |
  +--------------------------------------------------------------------+---------------------------------------+
  | Actual initial step size used                                      | :c:func:`IDAGetActualInitStep`        |
  +--------------------------------------------------------------------+---------------------------------------+
  | Step size used for the last step                                   | :c:func:`IDAGetLastStep`              |
  +--------------------------------------------------------------------+---------------------------------------+
  | Step size to be attempted on the next step                         | :c:func:`IDAGetCurrentStep`           |
  +--------------------------------------------------------------------+---------------------------------------+
  | Current internal time reached by the solver                        | :c:func:`IDAGetCurrentTime`           |
  +--------------------------------------------------------------------+---------------------------------------+
  | Suggested factor for tolerance scaling                             | :c:func:`IDAGetTolScaleFactor`        |
  +--------------------------------------------------------------------+---------------------------------------+
  | Error weight vector for state variables                            | :c:func:`IDAGetErrWeights`            |
  +--------------------------------------------------------------------+---------------------------------------+
  | Estimated local errors                                             | :c:func:`IDAGetEstLocalErrors`        |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of nonlinear solver iterations                                 | :c:func:`IDAGetNumNonlinSolvIters`    |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of nonlinear convergence failures                              | :c:func:`IDAGetNumNonlinSolvConvFails`|
  +--------------------------------------------------------------------+---------------------------------------+
  | Array showing roots found                                          | :c:func:`IDAGetRootInfo`              |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of calls to user root function                                 | :c:func:`IDAGetNumGEvals`             |
  +--------------------------------------------------------------------+---------------------------------------+
  | Name of constant associated with a return flag                     | :c:func:`IDAGetReturnFlagName`        |
  +--------------------------------------------------------------------+---------------------------------------+
  | Number of backtrack operations                                     | :c:func:`IDAGetNumBacktrackOps`       |
  +--------------------------------------------------------------------+---------------------------------------+
  | Corrected initial conditions                                       | :c:func:`IDAGetConsistentIC`          |
  +--------------------------------------------------------------------+---------------------------------------+
  | Size of real and integer workspace                                 | :c:func:`IDAGetLinWorkSpace`          |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of Jacobian evaluations                                        | :c:func:`IDAGetNumJacEvals`           |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of residual calls for finite diff. Jacobian-vector evals.      | :c:func:`IDAGetNumLinResEvals`        |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of linear iterations                                           | :c:func:`IDAGetNumLinIters`           |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of linear convergence failures                                 | :c:func:`IDAGetNumLinConvFails`       |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of preconditioner evaluations                                  | :c:func:`IDAGetNumPrecEvals`          |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of preconditioner solves                                       | :c:func:`IDAGetNumPrecSolves`         |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of Jacobian-vector setup evaluations                           | :c:func:`IDAGetNumJTSetupEvals`       |
  +--------------------------------------------------------------------+---------------------------------------+
  | No. of Jacobian-vector product evaluations                         | :c:func:`IDAGetNumJtimesEvals`        |
  +--------------------------------------------------------------------+---------------------------------------+
  | Last return from a linear solver function                          | :c:func:`IDAGetLastLinFlag`           |
  +--------------------------------------------------------------------+---------------------------------------+
  | Name of constant associated with a return flag                     | :c:func:`IDAGetLinReturnFlagName`     |
  +--------------------------------------------------------------------+---------------------------------------+


.. _IDA.Usage.CC.optional_output.optout_main:

Main solver optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IDA provides several user-callable functions that can be used to obtain
different quantities that may be of interest to the user, such as solver
workspace requirements, solver performance statistics, as well as additional
data from the IDA solver object (a suggested tolerance scaling factor, the error
weight vector, and the vector of estimated local errors). Also provided are
functions to extract statistics related to the performance of the
nonlinear solver being used. As a convenience, additional extraction functions
provide the optional outputs in groups. These optional output functions are
described next.

.. c:function:: int IDAGetWorkSpace(void * ida_mem, long int lenrw, long int leniw)

   The function ``IDAGetWorkSpace`` returns the IDA real and integer workspace
   sizes.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``lenrw`` -- number of real values in the IDA workspace.
      * ``leniw`` -- number of integer values in the IDA workspace.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      In terms of the problem size :math:`N`, the maximum method order
      ``maxord``, and the number of root functions ``nrtfn`` (see
      :numref:`IDA.Usage.CC.idarootinit`), the actual size of the real workspace, in
      :c:type:`realtype` words, is given by the following:

      * base value:
        :math:`\mathtt{lenrw} = 55 + (m + 6) * N_r + 3 * \mathtt{nrtfn}`;
      * with :c:func:`IDASVtolerances`:
        :math:`\mathtt{lenrw} = \mathtt{lenrw} + N_r`;
      * with constraint checking (see :c:func:`IDASetConstraints`):
        :math:`\mathtt{lenrw} = \mathtt{lenrw} + N_r`;
      * with ``id`` specified (see :c:func:`IDASetId`):
        :math:`\mathtt{lenrw} = \mathtt{lenrw} + N_r`;

      where :math:`m = \max(\mathtt{maxord}, 3)`, and :math:`N_r` is the number
      of real words in one ``N_Vector`` :math:`(\approx N)`.

      The size of the integer workspace (without distinction between ``int`` and
      ``long int`` words) is  given by:

      * base value: :math:`\mathtt{leniw} = 38 + (m + 6) * N_i + \mathtt{nrtfn}`;
      * with :c:func:`IDASVtolerances`:
        :math:`\mathtt{leniw} = \mathtt{leniw} + N_i`;
      * with constraint checking: :math:`\mathtt{lenrw} = \mathtt{lenrw} + N_i`;
      * with ``id`` specified (see :c:func:`IDASetId`):
        :math:`\mathtt{lenrw} = \mathtt{lenrw} + N_i`;

      where :math:`N_i` is the number of integer words in one ``N_Vector`` (= 1
      for the serial ``N_Vector`` and ``2 * npes`` for the parallel ``N_Vector``
      on ``npes`` processors). For the default value of ``maxord``, with no
      rootfinding, no ``id``, no constraints, and with no call to
      :c:func:`IDASVtolerances`, these lengths are given roughly by
      :math:`\mathtt{lenrw} = 55 + 11 * N` and :math:`\mathtt{leniw} = 49`.

.. c:function:: int IDAGetNumSteps(void * ida_mem, long int * nsteps)

   The function ``IDAGetNumSteps`` returns the cumulative number of internal
   steps taken by the solver (total so far).

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nsteps`` -- number of steps taken by IDA.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetNumResEvals(void * ida_mem, long int * nrevals)

   The function ``IDAGetNumResEvals`` returns the number of calls to the user's
   residual evaluation function.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nrevals`` -- number of calls to the user's ``res`` function.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      The ``nrevals`` value returned by :c:func:`IDAGetNumResEvals` does not
      account for calls made to ``res`` from a linear solver or preconditioner
      module.

.. c:function:: int IDAGetNumLinSolvSetups(void * ida_mem, long int * nlinsetups)

   The function ``IDAGetNumLinSolvSetups`` returns the cumulative number of
   calls made to the linear solver's setup function (total so far).

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nlinsetups`` -- number of calls made to the linear solver setup function.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetNumErrTestFails(void * ida_mem, long int * netfails)

   The function ``IDAGetNumErrTestFails`` returns the cumulative number of local
   error test failures that have occurred (total so far).

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``netfails`` -- number of error test failures.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetLastOrder(void * ida_mem, int * klast)

   The function ``IDAGetLastOrder`` returns the integration method order used
   during the last internal step.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``klast`` -- method order used on the last internal step.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetCurrentOrder(void * ida_mem, int * kcur)

   The function ``IDAGetCurrentOrder`` returns the integration method order to
   be used on the next internal step.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``kcur`` -- method order to be used on the next internal step.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetLastStep(void * ida_mem, realtype * hlast)

   The function ``IDAGetLastStep`` returns the integration step size taken on
   the last internal step (if from :c:func:`IDASolve`), or the last value of the
   artificial step size :math:`h` (if from :c:func:`IDACalcIC`).

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``hlast`` -- step size taken on the last internal step by IDA, or last
        artificial step size used in :c:func:`IDACalcIC` , whichever was called
        last.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetCurrentStep(void * ida_mem, realtype * hcur)

   The function ``IDAGetCurrentStep`` returns the integration step size to be
   attempted on the next internal step.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``hcur`` -- step size to be attempted on the next internal step.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetActualInitStep(void * ida_mem, realtype * hinused)

   The function ``IDAGetActualInitStep`` returns the value of the integration
   step size used on the first step.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``hinused`` -- actual value of initial step size.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**

   Even if the value of the initial integration step size was specified by the
   user through a call to :c:func:`IDASetInitStep`, this value might have been
   changed by IDA to ensure that the step size is within the prescribed bounds
   :math:`(h_{min} \le h_0 \le h_{max})`, or to meet the local error test.

.. c:function:: int IDAGetCurrentTime(void * ida_mem, realtype * tcur)

   The function ``IDAGetCurrentTime`` returns the current internal time reached
   by the solver.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``tcur`` -- current internal time reached.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetTolScaleFactor(void * ida_mem, realtype * tolsfac)

   The function ``IDAGetTolScaleFactor`` returns a suggested factor by which the
   user's tolerances should be scaled when too much accuracy has been requested
   for some internal step.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``tolsfac`` -- suggested scaling factor for user tolerances.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetErrWeights(void * ida_mem, N_Vector eweight)

   The function ``IDAGetErrWeights`` returns the solution error weights at the
   current time. These are the :math:`W_i` given by :eq:`IDA_errwt` (or by the
   user's :c:type:`IDAEwtFn`).

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``eweight`` -- solution error weights at the current time.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   .. warning::

      The user must allocate space for ``eweight``.

.. c:function:: int IDAGetEstLocalErrors(void * ida_mem, N_Vector ele)

   The function ``IDAGetEstLocalErrors`` returns the estimated local errors.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``ele`` -- estimated local errors at the current time.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   .. warning::

      The user must allocate space for ``ele``. The values returned in ``ele``
      are only valid if :c:func:`IDASolve` returned a non-negative value.

   .. note::

      The ``ele`` vector, togther with the ``eweight`` vector from
      :c:func:`IDAGetErrWeights`, can be used to determine how the various
      components of the system contributed to the estimated local error test.
      Specifically, that error test uses the RMS norm of a vector whose
      components are the products of the components of these two vectors.  Thus,
      for example, if there were recent error test failures, the components
      causing the failures are those with largest values for the products,
      denoted loosely as ``eweight[i]*ele[i]``.

.. c:function:: int IDAGetIntegratorStats(void *ida_mem, long int *nsteps, long int *nrevals, long int *nlinsetups, long int *netfails, int *qlast, int *qcur, realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur)

    The function ``IDAGetIntegratorStats`` returns the IDA integrator stats in
    one function call.

    **Arguments:**
       * ``ida_mem`` -- pointer to the IDA solver object.
       * ``nsteps`` -- cumulative number of steps taken by IDA.
       * ``nrevals`` -- cumulative number of calls to the user's ``res`` functions.
       * ``nlinsetups`` -- cumulative number of calls made to the linear solver
         setup function.
       * ``netfails`` -- cumulative number of error test failures.
       * ``klast`` -- method order used on the last internal step.
       * ``kcur`` -- method order to be used on the next internal step.
       * ``hinused`` -- actual value of initial step size.
       * ``hlast`` -- step sized taken on the last internal step.
       * ``hcur`` -- step size to be attempted on the next internal step.
       * ``tcur`` -- current internal time reached.

    **Return value:**
       * ``IDA_SUCCESS`` -- The optional output values have been successfully set.
       * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetNumNonlinSolvIters(void * ida_mem, long int * nniters)

   The function ``IDAGetNumNonlinSolvIters`` returns the cumulative number of
   nonlinear iterations performed.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nniters`` -- number of nonlinear iterations performed.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_MEM_FAIL`` -- The ``SUNNonlinearSolver`` object is ``NULL``.

.. c:function:: int IDAGetNumNonlinSolvConvFails(void * ida_mem, long int * nncfails)

   The function ``IDAGetNumNonlinSolvConvFails`` returns the cumulative number
   of nonlinear convergence failures that have occurred.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nncfails`` -- number of nonlinear convergence failures.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetNonlinSolvStats(void * ida_mem, long int * nniters, long int * nncfails)

   The function ``IDAGetNonlinSolvStats`` returns the IDA nonlinear solver
   statistics as a group.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nniters`` -- cumulative number of nonlinear iterations performed.
      * ``nncfails`` -- cumulative number of nonlinear convergence failures.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDA_MEM_FAIL`` -- The ``SUNNonlinearSolver`` object is ``NULL``.

.. c:function:: char* IDAGetReturnFlagName(long int flag)

   The function ``IDAGetReturnFlagName`` returns the name of the IDA constant
   corresponding to ``flag``.

   **Arguments:**
      * ``flag`` -- the flag returned by a call to an IDA function

   **Return value:**
      * ``char*`` -- the flag name string


.. _IDA.Usage.CC.optional_output.optout_iccalc:

Initial condition calculation optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. c:function:: int IDAGetNumBacktrackOps(void * ida_mem, long int * nbacktr)

   The function ``IDAGetNumBacktrackOps`` returns the number of backtrack
   operations done in the linesearch algorithm in :c:func:`IDACalcIC`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nbacktr`` -- the cumulative number of backtrack operations.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

.. c:function:: int IDAGetConsistentIC(void * ida_mem, N_Vector yy0_mod, N_Vector yp0_mod)

   The function ``IDAGetConsistentIC`` returns the corrected initial conditions
   calculated by :c:func:`IDACalcIC`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``yy0_mod`` -- consistent solution vector.
      * ``yp0_mod`` -- consistent derivative vector.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_ILL_INPUT`` -- The function was not called before the first call to
        :c:func:`IDASolve`.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      If the consistent solution vector or consistent derivative vector is not
      desired, pass ``NULL`` for the corresponding argument.

   .. warning::

      The user must allocate space for ``yy0_mod`` and ``yp0_mod`` (if not
      ``NULL``).



.. _IDA.Usage.CC.optional_output.optout_root:

Rootfinding optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are two optional output functions associated with rootfinding.

.. c:function:: int IDAGetRootInfo(void * ida_mem, int * rootsfound)

   The function ``IDAGetRootInfo`` returns an array showing which functions were
   found to have a root.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``rootsfound`` -- array of length ``nrtfn`` with the indices of the user
        functions :math:`g_i` found to have a root. For
        :math:`\mathtt{i} = 0, \ldots, \mathtt{nrtfn} -1`,
        :math:`\mathtt{rootsfound[i]} \ne 0` if :math:`g_i` has a root, and
        :math:`= 0` if not.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output values have been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.

   **Notes:**
      Note that, for the components :math:`g_i` for which a root was found, the
      sign of ``rootsfound[i]`` indicates the direction of zero-crossing. A value
      of :math:`+1` indicates that :math:`g_i` is increasing, while a value of
      :math:`-1` indicates a decreasing :math:`g_i`.

   .. warning::

      The user must allocate memory for the vector ``rootsfound``.

.. c:function:: int IDAGetNumGEvals(void * ida_mem, long int * ngevals)

   The function ``IDAGetNumGEvals`` returns the cumulative number of calls to
   the user root function :math:`g`.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``ngevals`` -- number of calls to the user's function :math:`g` so far.

   **Return value:**
      * ``IDA_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDA_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.


.. _IDA.Usage.CC.optional_output.optout_ls:

IDALS linear solver interface optional output functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following optional outputs are available from the IDALS modules:

.. c:function:: int IDAGetLinWorkSpace(void * ida_mem, long int * lenrwLS, long int * leniwLS)

   The function ``IDAGetLinWorkSpace`` returns the sizes of the real and integer
   workspaces used by the IDALS linear solver interface.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``lenrwLS`` -- the number of real values in the IDALS workspace.
      * ``leniwLS`` -- the number of integer values in the IDALS workspace.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   **Notes:**
      The workspace requirements reported by this routine correspond only to memory
      allocated within this interface and to memory allocated by the
      ``SUNLinearSolver`` object attached to it.  The template Jacobian
      matrix allocated by the user outside of IDALS is not included in this report.

   .. warning::

      The previous routines ``IDADlsGetWorkspace`` and ``IDASpilsGetWorkspace``
      are now wrappers for this routine, and may still be used for
      backward-compatibility.  However, these will be deprecated in future
      releases, so we recommend that users transition to the new routine name
      soon.

.. c:function:: int IDAGetNumJacEvals(void * ida_mem, long int * njevals)

   The function ``IDAGetNumJacEvals`` returns the cumulative number of calls to
   the IDALS Jacobian approximation function.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``njevals`` -- the cumulative number of calls to the Jacobian function total so far.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   .. warning::

      The previous routine ``IDADlsGetNumJacEvals`` is now a wrapper for this
      routine, and may still be used for backward-compatibility.  However, this
      will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.

.. c:function:: int IDAGetNumLinResEvals(void * ida_mem, long int * nrevalsLS)

   The function ``IDAGetNumLinResEvals`` returns the cumulative number of calls
   to the user residual function due to the finite difference Jacobian
   approximation or finite difference Jacobian-vector product approximation.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nrevalsLS`` -- the cumulative number of calls to the user residual
        function.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   **Notes:**
      The value ``nrevalsLS`` is incremented only if one of the default internal
      difference quotient functions is used.

   .. warning::

      The previous routines ``IDADlsGetNumRhsEvals`` and
      ``IDASpilsGetNumRhsEvals`` are now wrappers for this routine, and may
      still be used for backward-compatibility.  However, these will be
      deprecated in future releases, so we recommend that users transition to
      the new routine name soon.

.. c:function:: int IDAGetNumLinIters(void * ida_mem, long int * nliters)

   The function ``IDAGetNumLinIters`` returns the cumulative number of linear
   iterations.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nliters`` -- the current number of linear iterations.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   .. warning::

      The previous routine ``IDASpilsGetNumLinIters`` is now a wrapper for this
      routine, and may still be used for backward-compatibility.  However, this
      will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.

.. c:function:: int IDAGetNumLinConvFails(void * ida_mem, long int * nlcfails)

   The function ``IDAGetNumLinConvFails`` returns the cumulative number of
   linear convergence failures.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``nlcfails`` -- the current number of linear convergence failures.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   .. warning::

      The previous routine ``IDASpilsGetNumConvFails`` is now a wrapper for this
      routine, and may still be used for backward-compatibility.  However, this
      will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.

.. c:function:: int IDAGetNumPrecEvals(void * ida_mem, long int * npevals)

   The function ``IDAGetNumPrecEvals`` returns the cumulative number of
   preconditioner evaluations, i.e., the number of calls made to ``psetup``.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``npevals`` -- the cumulative number of calls to ``psetup``.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   .. warning::

      The previous routine ``IDASpilsGetNumPrecEvals`` is now a wrapper for this
      routine, and may still be used for backward-compatibility.  However, this
      will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.

.. c:function:: int IDAGetNumPrecSolves(void * ida_mem, long int * npsolves)

   The function ``IDAGetNumPrecSolves`` returns the cumulative number of calls
   made to the preconditioner solve function, ``psolve``.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``npsolves`` -- the cumulative number of calls to ``psolve``.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   .. warning::

      The previous routine ``IDASpilsGetNumPrecSolves`` is now a wrapper for
      this routine, and may still be used for backward-compatibility.  However,
      this will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.

.. c:function:: int IDAGetNumJTSetupEvals(void * ida_mem, long int * njtsetup)

   The function ``IDAGetNumJTSetupEvals`` returns the cumulative number of calls
   made to the Jacobian-vector product setup function ``jtsetup``.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``njtsetup`` -- the current number of calls to ``jtsetup``.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   .. warning::

      The previous routine ``IDASpilsGetNumJTSetupEvals`` is now a wrapper for
      this routine, and may still be used for backward-compatibility.  However,
      this will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.

.. c:function:: int IDAGetNumJtimesEvals(void * ida_mem, long int * njvevals)

   The function ``IDAGetNumJtimesEvals`` returns the cumulative number of calls
   made to the Jacobian-vector product function, ``jtimes``.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``njvevals`` -- the cumulative number of calls to ``jtimes``.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   .. warning::

      The previous routine ``IDASpilsGetNumJtimesEvals`` is now a wrapper for
      this routine, and may still be used for backward-compatibility.  However,
      this will be deprecated in future releases, so we recommend that users
      transition to the new routine name soon.

.. c:function:: int IDAGetLastLinFlag(void * ida_mem, long int * lsflag)

   The function ``IDAGetLastLinFlag`` returns the last return value from an
   IDALS routine.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``lsflag`` -- the value of the last return flag from an IDALS function.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer is ``NULL``.
      * ``IDALS_LMEM_NULL`` -- The IDALS linear solver has not been initialized.

   **Notes:**
      If the IDALS setup function failed (i.e., :c:func:`IDASolve` returned
      ``IDA_LSETUP_FAIL``) when using the :ref:`SUNLINSOL_DENSE <SUNLinSol_Dense>`
      or :ref:`SUNLINSOL_BAND <SUNLinSol_Band>` modules, then the value of
      ``lsflag`` is equal to the column index (numbered from one) at which a zero
      diagonal element was encountered during the LU factorization of the (dense or
      banded) Jacobian matrix.  If the IDALS setup function failed when using
      another ``SUNLinearSolver`` object, then ``lsflag`` will be
      ``SUNLS_PSET_FAIL_UNREC``, ``SUNLS_ASET_FAIL_UNREC``, or
      ``SUNLS_PACKAGE_FAIL_UNREC``.  If the IDALS solve function failed
      (:c:func:`IDASolve` returned ``IDA_LSOLVE_FAIL``), ``lsflag`` contains the
      error return flag from the ``SUNLinearSolver`` object, which will be one of:
      ``SUNLS_MEM_NULL``, indicating that the ``SUNLinearSolver`` memory is
      ``NULL``; ``SUNLS_ATIMES_FAIL_UNREC``, indicating an unrecoverable failure in
      the :math:`J*v` function; ``SUNLS_PSOLVE_FAIL_UNREC``, indicating that the
      preconditioner solve function ``psolve`` failed unrecoverably;
      ``SUNLS_GS_FAIL``, indicating a failure in the Gram-Schmidt procedure
      (generated only in SPGMR or SPFGMR); ``SUNLS_QRSOL_FAIL``, indicating that
      the matrix :math:`R` was found to be singular during the QR solve phase
      (SPGMR and SPFGMR only); or ``SUNLS_PACKAGE_FAIL_UNREC``, indicating an
      unrecoverable failure in an external iterative linear solver package.

   .. warning::

      The previous routines ``IDADlsGetLastFlag`` and ``IDASpilsGetLastFlag``
      are now wrappers for this routine, and may still be used for
      backward-compatibility.  However, these will be deprecated in future
      releases, so we recommend that users transition to the new routine name
      soon.

.. c:function:: char* IDAGetLinReturnFlagName(long int lsflag)

   The function ``IDAGetLinReturnFlagName`` returns the name of the IDALS
   constant corresponding to ``lsflag``.

   **Arguments:**
      * ``flag`` -- the flag returned by a call to an IDA function

   **Return value:**
      * ``char*`` -- the flag name string or if
        :math:`1 \leq \mathtt{lsflag} \leq N` (LU factorization failed), this
        function returns "NONE".

   .. warning::

      The previous routines ``IDADlsGetReturnFlagName`` and
      ``IDASpilsGetReturnFlagName`` are now wrappers for this routine, and may
      still be used for backward-compatibility.  However, these will be
      deprecated in future releases, so we recommend that users transition to
      the new routine name soon.


.. _IDA.Usage.CC.reinit:

IDA reinitialization function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The function :c:func:`IDAReInit` reinitializes the main IDA solver for the
solution of a new problem, where a prior call to :c:func:`IDAInit` has been
made. The new problem must have the same size as the previous
one. :c:func:`IDAReInit` performs the same input checking and initializations
that :c:func:`IDAInit` does, but does no memory allocation, as it assumes that
the existing internal memory is sufficient for the new problem. A call to
:c:func:`IDAReInit` deletes the solution history that was stored internally
during the previous integration. Following a successful call to
:c:func:`IDAReInit`, call :c:func:`IDASolve` again for the solution of the new
problem.

The use of :c:func:`IDAReInit` requires that the maximum method order,
``maxord``, is no larger for the new problem than for the problem specified in
the last call to :c:func:`IDAInit`. In addition, the same ``N_Vector`` module
set for the previous problem will be reused for the new problem.

If there are changes to the linear solver specifications, make the appropriate
calls to either the linear solver objects themselves, or to the IDALS interface
routines, as described in :numref:`IDA.Usage.CC.callable_fct_sim.lin_solv_init`.

If there are changes to any optional inputs, make the appropriate ``IDASet***``
calls, as described in :numref:`IDA.Usage.CC.optional_input.optin_main`. Otherwise,
all solver inputs set previously remain in effect.

One important use of the :c:func:`IDAReInit` function is in the treating of jump
discontinuities in the residual function. Except in cases of fairly small jumps,
it is usually more efficient to stop at each point of discontinuity and restart
the integrator with a readjusted DAE model, using a call to :c:func:`IDAReInit`.
To stop when the location of the discontinuity is known, simply make that
location a value of :math:`t_{\text{out}}`. To stop when the location of the
discontinuity is determined by the solution, use the rootfinding feature. In
either case, it is critical that the residual function *not* incorporate the
discontinuity, but rather have a smooth extention over the discontinuity, so
that the step across it (and subsequent rootfinding, if used) can be done
efficiently. Then use a switch within the residual function (communicated
through ``user_data``) that can be flipped between the stopping of the
integration and the restart, so that the restarted problem uses the new values
(which have jumped).  Similar comments apply if there is to be a jump in the
dependent variable vector.

.. c:function:: int IDAReInit(void * ida_mem, realtype t0, N_Vector y0, N_Vector yp0)

   The function ``IDAReInit`` provides required problem specifications and
   reinitializes IDA.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``t0`` -- is the initial value of :math:`t`.
      * ``y0`` -- is the initial value of :math:`y`.
      * ``yp0`` -- is the initial value of :math:`\dot{y}`.

   **Return value:**
      * ``IDA_SUCCESS`` -- The call to was successful.
      * ``IDA_MEM_NULL`` -- The IDA solver object was not initialized through a
        previous call to :c:func:`IDACreate`.
      * ``IDA_NO_MALLOC`` -- Memory space for the IDA solver object was not
        allocated through a previous call to :c:func:`IDAInit`.
      * ``IDA_ILL_INPUT`` -- An input argument to :c:func:`IDAReInit` has an
        illegal value.

   **Notes:**
      If an error occurred, :c:func:`IDAReInit` also sends an error message to the
      error handler function.


.. _IDA.Usage.CC.user_fct_sim:

User-supplied functions
-----------------------

The user-supplied functions consist of one function defining the DAE residual,
(optionally) a function that handles error and warning messages, (optionally) a
function that provides the error weight vector, (optionally) one or two
functions that provide Jacobian-related information for the linear solver, and
(optionally) one or two functions that define the preconditioner for use in any
of the Krylov iteration algorithms.

.. _IDA.Usage.CC.user_fct_sim.resFn:

DAE residual function
~~~~~~~~~~~~~~~~~~~~~

The user must provide a function of type :c:type:`IDAResFn` defined as follows:

.. c:type:: int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)

   This function computes the problem residual for given values of the
   independent variable :math:`t`, state vector :math:`y`, and derivative
   :math:`\dot{y}`.

   **Arguments:**
      * ``tt`` -- is the current value of the independent variable.
      * ``yy`` -- is the current value of the dependent variable vector,
        :math:`y(t)`.
      * ``yp`` -- is the current value of :math:`\dot{y}(t)`.
      * ``rr`` -- is the output residual vector :math:`F(t,y,\dot{y})`.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        pointer parameter passed to :c:func:`IDASetUserData`.

   **Return value:**
      An :c:type:`IDAResFn` function type should return a value of :math:`0` if
      successful, a positive value if a recoverable error occurred (e.g., ``yy``
      has an illegal value), or a negative value if a nonrecoverable error
      occurred. In the last case, the integrator halts. If a recoverable error
      occurred, the integrator will attempt to correct and retry.

   **Notes:**
      A recoverable failure error return from the :c:type:`IDAResFn` is typically
      used to flag a value of the dependent variable :math:`y` that is "illegal" in
      some way (e.g., negative where only a non-negative value is physically
      meaningful).  If such a return is made, IDA will attempt to recover (possibly
      repeating the nonlinear solve, or reducing the step size) in order to avoid
      this recoverable error return.

      For efficiency reasons, the DAE residual function is not evaluated at the
      converged solution of the nonlinear solver. Therefore, in general, a
      recoverable error in that converged value cannot be corrected.  (It may be
      detected when the residual function is called the first time during the
      following integration step, but a successful step cannot be undone.)


.. _IDA.Usage.CC.user_fct_sim.ehFn:

Error message handler function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an alternative to the default behavior of directing error and warning
messages to the file pointed to by ``errfp`` (see :c:func:`IDASetErrFile`), the
user may provide a function of type :c:type:`IDAErrHandlerFn` to process any
such messages.  The function type :c:type:`IDAErrHandlerFn` is defined as
follows:

.. c:type:: void (*IDAErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *user_data)

   This function processes error and warning messages from IDA and its
   sub-modules.

   **Arguments:**
      * ``error_code`` -- is the error code.
      * ``module`` -- is the name of the IDA module reporting the error.
      * ``function`` -- is the name of the function in which the error occurred.
      * ``eH_data`` -- is a pointer to user data, the same as the ``eh_data``
        parameter passed to :c:func:`IDASetErrHandlerFn`.

   **Return value:**
      This function has no return value.

   **Notes:**
      ``error_code`` is negative for errors and positive (``IDA_WARNING``) for
      warnings. If a function that returns a pointer to memory encounters an error,
      it sets ``error_code`` to 0.


.. _IDA.Usage.CC.user_fct_sim.ewtsetFn:

Error weight function
~~~~~~~~~~~~~~~~~~~~~

.. c:type:: int (*IDAEwtFn)(N_Vector y, N_Vector ewt, void *user_data)

   This function computes the WRMS error weights for the vector :math:`y`.

   **Arguments:**
      * ``y`` -- is the value of the dependent variable vector at which the weight
        vector is to be computed.
      * ``ewt`` -- is the output vector containing the error weights.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.

   **Return value:**
      * ``0`` -- if it the error weights were successfully set.
      * ``-1`` -- if any error occured.

   **Notes:**
      Allocation of memory for ``ewt`` is handled within IDA.

   .. warning::

      The error weight vector must have all components positive. It is the
      user's responsiblity to perform this test and return -1 if it is not
      satisfied.


.. _IDA.Usage.CC.user_fct_sim.rootFn:

Rootfinding function
~~~~~~~~~~~~~~~~~~~~

If a rootfinding problem is to be solved during the integration of the DAE
system, the user must supply a function of type :c:type:`IDARootFn`, defined
as follows:

.. c:type:: int (*IDARootFn)(realtype t, N_Vector y, N_Vector yp, realtype *gout, void *user_data)

   This function computes a vector-valued function :math:`g(t,y,\dot{y})` such
   that the roots of the ``nrtfn`` components :math:`g_i(t,y,\dot{y})` are to be
   found during the integration.

   **Arguments:**
      * ``t`` -- is the current value of the independent variable.
      * ``y`` -- is the current value of the dependent variable vector,
        :math:`y(t)`.
      * ``yp`` -- is the current value of :math:`\dot{y}(t)`, the
        :math:`t-\text{derivative}` of :math:`y`.
      * ``gout`` -- is the output array, of length ``nrtfn``, with components
        :math:`g_i(t,y,\dot{y})`.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.

   **Return value:**
      ``0`` if successful or non-zero if an error occured (in which case the
      integration is halted and :c:func:`IDASolve` returs ``IDA_RTFUNC_FAIL``).

   **Notes:**
      Allocation of memory for ``gout`` is handled within IDA.


.. _IDA.Usage.CC.user_fct_sim.jacFn:

Jacobian construction (matrix-based linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a matrix-based linear solver module is used (i.e. a non-``NULL``
``SUNMatrix`` object was supplied to :c:func:`IDASetLinearSolver`), the
user may provide a function of type :c:type:`IDALsJacFn` defined as follows:

.. c:type:: int (*IDALsJacFn)(realtype t, realtype c_j, N_Vector y, N_Vector yp, N_Vector r, SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)

   This function computes the Jacobian matrix :math:`J` of the DAE system (or an
   approximation to it), defined by :eq:`IDA_DAE_Jacobian`.

   **Arguments:**
      * ``tt`` -- is the current value of the independent variable :math:`t`.
      * ``cj`` -- is the scalar in the system Jacobian, proportional to the inverse
        of the step size (:math:`\alpha` in :eq:`IDA_DAE_Jacobian`).
      * ``yy`` -- is the current value of the dependent variable vector,
        :math:`y(t)`.
      * ``yp`` -- is the current value of :math:`\dot{y}(t)`.
      * ``rr`` -- is the current value of the residual vector
        :math:`F(t,y,\dot{y})`.
      * ``Jac`` -- is the output (approximate) Jacobian matrix (of type
        ``SUNMatrix``),
        :math:`J = \dfrac{\partial{F}}{\partial{y}} + cj ~ \dfrac{\partial{F}}{\partial{\dot{y}}}`.
      * ``user_data`` - is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.
      * ``tmp1``, ``tmp2``, and ``tmp3`` -- are pointers to memory allocated for
        variables of type ``N_Vector`` which can be used by ``IDALsJacFn`` function
        as temporary storage or work space.

   **Return value:**
      An :c:type:`IDALsJacFn` should return :math:`0` if successful, a positive
      value if a recoverable error occurred, or a negative value if a
      nonrecoverable error occurred.

      In the case of a recoverable eror return, the integrator will attempt to
      recover by reducing the stepsize, and hence changing :math:`\alpha` in
      :eq:`IDA_DAE_Jacobian`.

   **Notes:**
      Information regarding the structure of the specific ``SUNMatrix``
      structure (e.g., number of rows, upper/lower bandwidth, sparsity type) may be
      obtained through using the implementation-specific ``SUNMatrix``
      interface functions (see Chapter :numref:`SUNMatrix` for details).

      With direct linear solvers (i.e., linear solvers with type
      ``SUNLINEARSOLVER_DIRECT``), the Jacobian matrix :math:`J(t,y,\dot{y})` is
      zeroed out prior to calling the user-supplied Jacobian function so only
      nonzero elements need to be loaded into ``Jac``.

      With the default nonlinear solver (the native SUNDIALS Newton method), each
      call to the user’s ``IDALsJacFn`` function is preceded by a call to the
      ``IDAResFn`` user function with the same (``tt``, ``yy``, ``yp``)
      arguments. Thus the Jacobian function can use any auxiliary data that is
      computed and saved during the evaluation of the DAE residual.  In the case of
      a user-supplied or external nonlinear solver, this is also true if the
      residual function is evaluated prior to calling the linear solver setup
      function (see :numref:`SUNNonlinSol.API.SUNSuppliedFn` for more information).

      If the user’s :c:type:`IDALsJacFn` function uses difference quotient
      approximations, it may need to access quantities not in the call list. These
      quantities may include the current stepsize, the error weights, etc. To
      obtain these, the user will need to add a pointer to ``ida_mem`` to
      ``user_data`` and then use the ``IDAGet*`` functions described in
      :numref:`IDA.Usage.CC.optional_output.optout_main`. The unit roundoff can be
      accessed as ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

      **dense:**

      A user-supplied dense Jacobian function must load the ``Neq`` :math:`\times`
      ``Neq`` dense matrix ``Jac`` with an approximation to the Jacobian matrix
      :math:`J(t,y,\dot{y})` at the point (``tt``, ``yy``, ``yp``). The accessor
      macros ``SM_ELEMENT_D`` and ``SM_COLUMN_D`` allow the user to read and write
      dense matrix elements without making explicit references to the underlying
      representation of the ``SUNMATRIX_DENSE`` type. ``SM_ELEMENT_D(J, i, j)``
      references the (``i``, ``j``)-th element of the dense matrix ``Jac`` (with
      ``i``, ``j``\ :math:`= 0\ldots \texttt{N}-1`). This macro is meant for small
      problems for which efficiency of access is not a major concern. Thus, in
      terms of the indices :math:`m` and :math:`n` ranging from :math:`1` to
      :math:`N`, the Jacobian element :math:`J_{m,n}` can be set using the
      statement ``SM_ELEMENT_D(J, m-1, n-1) =`` :math:`J_{m,n}`. Alternatively,
      ``SM_COLUMN_D(J, j)`` returns a pointer to the first element of the ``j``-th
      column of ``Jac`` (with ``j``\ :math:`= 0\ldots \texttt{N}-1`), and the
      elements of the ``j``-th column can then be accessed using ordinary array
      indexing. Consequently, :math:`J_{m,n}` can be loaded using the statements
      ``col_n = SM_COLUMN_D(J, n-1);`` ``col_n[m-1] =`` :math:`J_{m,n}`.  For large
      problems, it is more efficient to use ``SM_COLUMN_D`` than to use
      ``SM_ELEMENT_D``. Note that both of these macros number rows and columns
      starting from :math:`0`. The ``SUNMATRIX_DENSE`` type and accessor macros are
      documented in :numref:`SUNMatrix.Dense`.

      **banded**:

      A user-supplied banded Jacobian function must load the ``Neq`` :math:`\times`
      ``Neq`` banded matrix ``Jac`` with an approximation to the Jacobian matrix
      :math:`J(t,y,\dot{y})` at the point (``tt``, ``yy``, ``yp``). The accessor
      macros ``SM_ELEMENT_B``, ``SM_COLUMN_B``, and ``SM_COLUMN_ELEMENT_B`` allow
      the user to read and write banded matrix elements without making specific
      references to the underlying representation of the ``SUNMATRIX_BAND`` type.
      ``SM_ELEMENT_B(J, i, j)`` references the (``i``, ``j``)-th element of the
      banded matrix ``Jac``, counting from :math:`0`. This macro is meant for use
      in small problems for which efficiency of access is not a major
      concern. Thus, in terms of the indices :math:`m` and :math:`n` ranging from
      :math:`1` to :math:`\texttt{N}` with :math:`(m,n)` within the band defined by
      ``mupper`` and ``mlower``, the Jacobian element :math:`J_{m,n}` can be loaded
      using the statement ``SM_ELEMENT_B(J, m-1, n-1) =`` :math:`J_{m,n}`. The
      elements within the band are those with ``-mupper`` :math:`\le` ``m-n``
      :math:`\le` ``mlower``. Alternatively, ``SM_COLUMN_B(J, j)`` returns a
      pointer to the diagonal element of the ``j``-th column of ``Jac``, and if we
      assign this address to ``realtype *col_j``, then the ``i``-th element of the
      ``j``-th column is given by ``SM_COLUMN_ELEMENT_B(col_j, i, j)``, counting
      from :math:`0`. Thus, for :math:`(m,n)` within the band, :math:`J_{m,n}` can
      be loaded by setting ``col_n = SM_COLUMN_B(J, n-1);`` and
      ``SM_COLUMN_ELEMENT_B(col_n, m-1, n-1) =`` :math:`J_{m,n}`. The elements of
      the ``j``-th column can also be accessed via ordinary array indexing, but
      this approach requires knowledge of the underlying storage for a band matrix
      of type ``SUNMATRIX_BAND``. The array ``col_n`` can be indexed from
      :math:`-`\ ``mupper`` to ``mlower``. For large problems, it is more efficient
      to use ``SM_COLUMN_B`` and ``SM_COLUMN_ELEMENT_B`` than to use the
      ``SM_ELEMENT_B`` macro. As in the dense case, these macros all number rows
      and columns starting from :math:`0`. The ``SUNMATRIX_BAND`` type and accessor
      macros are documented in :numref:`SUNMatrix.Band`.

      **sparse**:

      A user-supplied sparse Jacobian function must load the ``Neq`` :math:`\times`
      ``Neq`` compressed-sparse-column or compressed-sparse-row matrix ``Jac`` with
      an approximation to the Jacobian matrix :math:`J(t,y,\dot{y})` at the point
      (``tt``, ``yy``, ``yp``). Storage for ``Jac`` already exists on entry to this
      function, although the user should ensure that sufficient space is allocated
      in ``Jac`` to hold the nonzero values to be set; if the existing space is
      insufficient the user may reallocate the data and index arrays as needed. The
      amount of allocated space in a ``SUNMATRIX_SPARSE`` object may be accessed
      using the macro ``SM_NNZ_S`` or the routine ``SUNSparseMatrix_NNZ``. The
      ``SUNMATRIX_SPARSE`` type and accessor macros are documented in
      :numref:`SUNMatrix.Sparse`.

   .. warning::

      The previous function type ``IDADlsJacFn`` is identical to ``IDALsJacFn``,
      and may still be used for backward-compatibility. However, this will be
      deprecated in future releases, so we recommend that users transition to
      the new function type name soon.


.. _IDA.Usage.CC.user_fct_sim.jtimesFn:

Jacobian-vector product (matrix-free linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a matrix-free linear solver is to be used (i.e., a ``NULL``-valued
``SUNMatrix`` was supplied to :c:func:`IDASetLinearSolver`), the user may
provide a function of type :c:type:`IDALsJacTimesVecFn` in the following form, to
compute matrix-vector products :math:`Jv`. If such a function is not supplied,
the default is a difference quotient approximation to these products.

.. c:type:: int (*IDALsJacTimesVecFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector v, N_Vector Jv, realtype cj, void *user_data, N_Vector tmp1, N_Vector tmp2)

   This function computes the product :math:`Jv` of the DAE system Jacobian
   :math:`J` (or an approximation to it) and a given vector ``v``, where
   :math:`J` is defined by :eq:`IDA_DAE_Jacobian`.

   **Arguments:**
      * ``tt`` -- is the current value of the independent variable.
      * ``yy`` -- is the current value of the dependent variable vector,
        :math:`y(t)`.
      * ``yp`` -- is the current value of :math:`\dot{y}(t)`.
      * ``rr`` -- is the current value of the residual vector
        :math:`F(t,y,\dot{y})`.
      * ``v`` -- is the vector by which the Jacobian must be multiplied to the
        right.
      * ``Jv`` -- is the computed output vector.
      * ``cj`` -- is the scalar in the system Jacobian, proportional to the inverse
        of the step size (:math:`\alpha` in :eq:`IDA_DAE_Jacobian`).
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.
      * ``tmp1`` and ``tmp2`` -- are pointers to memory allocated for variables of
        type ``N_Vector`` which can be used by :c:type:`IDALsJacTimesVecFn` as
        temporary storage or work space.

   **Return value:**
      The value returned by the Jacobian-times-vector function should be 0 if
      successful.  A nonzero value indicates that a nonrecoverable error occurred.

   **Notes:**
      This function must return a value of :math:`Jv` that uses an approximation
      to the **current** value of :math:`J`, i.e. as evaluated at the current
      :math:`(t,y,\dot{y})`.

      If the user’s ``IDALsJacTimesVecFn`` function uses difference quotient
      approximations, it may need to access quantities not in the call list. These
      include the current stepsize, the error weights, etc. To obtain these, the
      user will need to add a pointer to ``ida_mem`` to ``user_data`` and then use
      the ``IDAGet*`` functions described in
      :numref:`IDA.Usage.CC.optional_output.optout_main`. The unit roundoff can be
      accessed as ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

   .. warning::

      The previous function type ``IDASpilsJacTimesVecFn`` is identical to
      ``IDALsJacTimesVecFn``, and may still be used for
      backward-compatibility. However, this will be deprecated in future
      releases, so we recommend that users transition to the new function type
      name soon.


.. _IDA.Usage.CC.user_fct_sim.jtsetupFn:

Jacobian-vector product setup (matrix-free linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the user's Jacobian-vector product function requires that any
Jacobian-related data be preprocessed or evaluated, then this needs to be done
in a user-supplied function of type :c:type:`IDALsJacTimesSetupFn`, defined as
follows:

.. c:type:: int (*IDALsJacTimesSetupFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, ealtype cj, void *user_data);

   This function setups any data needed by :math:`Jv` product function (see
   :c:type:`IDALsJacTimesVecFn`).

   **Arguments:**
      * ``tt`` -- is the current value of the independent variable.
      * ``yy`` -- is the current value of the dependent variable vector,
        :math:`y(t)`.
      * ``yp`` -- is the current value of :math:`\dot{y}(t)`.
      * ``rr`` -- is the current value of the residual vector
        :math:`F(t,y,\dot{y})`.
      * ``cj`` -- is the scalar in the system Jacobian, proportional to the inverse
        of the step size (:math:`\alpha` in :eq:`IDA_DAE_Jacobian`).
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.

   **Return value:**
      The value returned by the Jacobian-vector setup function should be 0 if
      successful, positive for a recoverable error (in which case the step will be
      retried), or negative for an unrecoverable error (in which case the
      integration is halted).

   **Notes:**
      Each call to the Jacobian-vector product setup function is preceded by a call
      to the :c:type:`IDAResFn` user function with the same :math:`(t,y,\dot{y})`
      arguments. Thus, the setup function can use any auxiliary data that is
      computed and saved during the evaluation of the DAE residual.

      If the user’s :c:type:`IDALsJacTimesVecFn` function uses difference quotient
      approximations, it may need to access quantities not in the call list. These
      include the current stepsize, the error weights, etc. To obtain these, the
      user will need to add a pointer to ``ida_mem`` to ``user_data`` and then use
      the ``IDAGet*`` functions described in
      :numref:`IDA.Usage.CC.optional_output.optout_main`. The unit roundoff can be
      accessed as ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.

   .. warning::

      The previous function type ``IDASpilsJacTimesSetupFn`` is identical to
      ``IDALsJacTimesSetupFn``, and may still be used for
      backward-compatibility. However, this will be deprecated in future
      releases, so we recommend that users transition to the new function type
      name soon.



.. _IDA.Usage.CC.user_fct_sim.psolveFn:

Preconditioner solve (iterative linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a user-supplied preconditioner is to be used with a ``SUNLinearSolver``
solver module, then the user must provide a function to solve the linear system
:math:`Pz = r` where :math:`P` is a left preconditioner matrix which
approximates (at least crudely) the Jacobian matrix :math:`J =
\dfrac{\partial{F}}{\partial{y}} + cj ~ \dfrac{\partial{F}}{\partial{\dot{y}}}`. This function
must be of type :c:type:`IDALsPrecSolveFn`, defined as follows:


.. c:type:: int (*IDALsPrecSolveFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, realtype cj, realtype delta, void *user_data)

   This function solves the preconditioning system :math:`Pz = r`.

   **Arguments:**
      * ``tt`` -- is the current value of the independent variable.
      * ``yy`` -- is the current value of the dependent variable vector,
        :math:`y(t)`.
      * ``yp`` -- is the current value of :math:`\dot{y}(t)`.
      * ``rr`` -- is the current value of the residual vector
        :math:`F(t,y,\dot{y})`.
      * ``rvec`` -- is the right-hand side vector :math:`r` of the linear system to
        be solved.
      * ``zvec`` -- is the computed output vector.
      * ``cj`` -- is the scalar in the system Jacobian, proportional to the inverse
        of the step size (:math:`\alpha` in :eq:`IDA_DAE_Jacobian`).
      * ``delta`` -- is an input tolerance to be used if an iterative method is
        employed in the solution.  In that case, the residual vector :math:`Res =
        r - P z` of the system should be made less than ``delta`` in weighted
        :math:`l_2` norm, i.e., :math:`\sqrt{\displaystyle\sum_i (Res_i \cdot ewt_i)^2 } <`
        ``delta``. To obtain the ``N_Vector`` ``ewt``, call
        :c:func:`IDAGetErrWeights`.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.

   **Return value:**
      The value returned by the preconditioner solve function should be 0 if
      successful, positive for a recoverable error (in which case the step will be
      retried), or negative for an unrecoverable error (in which case the
      integration is halted).


.. _IDA.Usage.CC.user_fct_sim.precondFn:

Preconditioner setup (iterative linear solvers)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the user’s preconditioner requires that any Jacobian-related data be
evaluated or preprocessed, then this needs to be done in a user-supplied
function of type :c:type:`IDALsPrecSetupFn`, defined as follows:

.. c:type:: int (*IDALsPrecSetupFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, realtype cj, void *user_data)

   This function solves the preconditioning system :math:`Pz = r`.

   **Arguments:**
      * ``tt`` -- is the current value of the independent variable.
      * ``yy`` -- is the current value of the dependent variable vector,
        :math:`y(t)`.
      * ``yp`` -- is the current value of :math:`\dot{y}(t)`.
      * ``rr`` -- is the current value of the residual vector
        :math:`F(t,y,\dot{y})`.
      * ``cj`` -- is the scalar in the system Jacobian, proportional to the inverse
        of the step size (:math:`\alpha` in :eq:`IDA_DAE_Jacobian`).
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.

   **Return value:**
      The value returned by the preconditioner setup function should be 0 if
      successful, positive for a recoverable error (in which case the step will be
      retried), or negative for an unrecoverable error (in which case the
      integration is halted).

   **Notes:**
      With the default nonlinear solver (the native SUNDIALS Newton method), each
      call to the preconditioner setup function is preceded by a call to the
      :c:type:`IDAResFn` user function with the same :math:`(t,y,\dot{y})`
      arguments. Thus the preconditioner setup function can use any auxiliary data
      that is computed and saved during the evaluation of the DAE residual. In the
      case of a user-supplied or external nonlinear solver, this is also true if
      the residual function is evaluated prior to calling the linear solver setup
      function (see :numref:`SUNNonlinSol.API.SUNSuppliedFn` for more information).

      This function is not called in advance of every call to the preconditioner
      solve function, but rather is called only as often as needed to achieve
      convergence in the nonlinear solver.

      If the user’s :c:type:`IDALsPrecSetupFn` function uses difference quotient
      approximations, it may need to access quantities not in the call list. These
      include the current stepsize, the error weights, etc. To obtain these, the
      user will need to add a pointer to ``ida_mem`` to ``user_data`` and then use
      the ``IDAGet*`` functions described in
      :numref:`IDA.Usage.CC.optional_output.optout_main`. The unit roundoff can be
      accessed as ``UNIT_ROUNDOFF`` defined in ``sundials_types.h``.


.. _IDA.Usage.CC.precond:

Preconditioner modules
----------------------

A principal reason for using a parallel DAE solver such as IDA lies in the
solution of partial differential equations (PDEs). Moreover, the use of a Krylov
iterative method for the solution of many such problems is motivated by the
nature of the underlying linear system of equations :eq:`IDA_DAE_Newtoncorr` that
must be solved at each time step. The linear algebraic system is large, sparse,
and structured. However, if a Krylov iterative method is to be effective in this
setting, then a nontrivial preconditioner needs to be used. Otherwise, the rate
of convergence of the Krylov iterative method is usually unacceptably
slow. Unfortunately, an effective preconditioner tends to be problem-specific.

However, we have developed one type of preconditioner that treats a rather broad
class of PDE-based problems. It has been successfully used for several
realistic, large-scale problems :cite:p:`HiTa:98` and is included in a software
module within the IDA package. This module works with the parallel vector module
:ref:`NVECTOR_PARALLEL <NVectors.NVParallel>` and generates a preconditioner
that is a block-diagonal matrix with each block being a band matrix. The blocks
need not have the same number of super- and sub-diagonals, and these numbers may
vary from block to block. This Band-Block-Diagonal Preconditioner module is
called IDABBDPRE.

One way to envision these preconditioners is to think of the domain of the
computational PDE problem as being subdivided into :math:`M` non-overlapping
sub-domains. Each of these sub-domains is then assigned to one of the :math:`M`
processors to be used to solve the DAE system. The basic idea is to isolate the
preconditioning so that it is local to each processor, and also to use a
(possibly cheaper) approximate residual function. This requires the definition
of a new function :math:`G(t,y,\dot{y})` which approximates the function
:math:`F(t, y, \dot{y})` in the definition of the DAE system :eq:`IDA_DAE`. However,
the user may set :math:`G = F`. Corresponding to the domain decomposition, there
is a decomposition of the solution vectors :math:`y` and :math:`\dot{y}` into
:math:`M` disjoint blocks :math:`y_m` and :math:`\dot{y}_m`, and a decomposition
of :math:`G` into blocks :math:`G_m`. The block :math:`G_m` depends on
:math:`y_m` and :math:`\dot{y}_m`, and also on components of :math:`y_{m'}` and
:math:`\dot{y}_{m'}` associated with neighboring sub-domains (so-called
ghost-cell data). Let :math:`\bar{y}_m` and :math:`\bar{\dot{y}}_m` denote
:math:`y_m` and :math:`\dot{y}_m` (respectively) augmented with those other
components on which :math:`G_m` depends. Then we have

.. math::

   G(t,y,\dot{y}) = [G_1(t,\bar{y}_1,\bar{\dot{y}}_1), G_2(t,\bar{y}_2,\bar{\dot{y}}_2),
                  \ldots, G_M(t,\bar{y}_M,\bar{\dot{y}}_M)]^T ~,

and each of the blocks :math:`G_m(t,\bar{y}_m,\bar{\dot{y}}_m)` is uncoupled
from the others.

The preconditioner associated with this decomposition has the form

.. math:: P= \begin{bmatrix} P_1 & & &\\ & P_2 & &\\ & & \ddots & \\ & & & P_M\end{bmatrix}

where

.. math::

   P_m \approx \frac{\partial G_m}{\partial y_m}
     + \alpha \frac{\partial G_m}{\partial \dot{y}_m}

This matrix is taken to be banded, with upper and lower half-bandwidths ``mudq``
and ``mldq`` defined as the number of non-zero diagonals above and below the
main diagonal, respectively. The difference quotient approximation is computed
using ``mudq`` :math:`+` ``mldq`` :math:`+ 2` evaluations of :math:`G_m`, but
only a matrix of bandwidth ``mukeep`` :math:`+` ``mlkeep`` :math:`+ 1` is
retained.

Neither pair of parameters need be the true half-bandwidths of the Jacobians of
the local block of :math:`G`, if smaller values provide a more efficient
preconditioner. Such an efficiency gain may occur if the couplings in the DAE
system outside a certain bandwidth are considerably weaker than those within the
band. Reducing ``mukeep`` and ``mlkeep`` while keeping ``mudq`` and ``mldq`` at
their true values, discards the elements outside the narrower band. Reducing
both pairs has the additional effect of lumping the outer Jacobian elements into
the computed elements within the band, and requires more caution and
experimentation.

The solution of the complete linear system

.. math:: Px = b

reduces to solving each of the equations

.. math:: P_m x_m = b_m

and this is done by banded LU factorization of :math:`P_m` followed by a banded
backsolve.

Similar block-diagonal preconditioners could be considered with different
treatment of the blocks :math:`P_m`. For example, incomplete LU factorization or
an iterative method could be used instead of banded LU factorization.


.. _IDA.Usage.CC.precond.idabbdpre:

A parallel band-block-diagonal preconditioner module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The IDABBDPRE module calls two user-provided functions to construct :math:`P`: a
required function ``Gres`` (of type :c:type:`IDABBDLocalFn`) which approximates
the residual function :math:`G(t,y,\dot{y}) \approx F(t,y,\dot{y})` and which is
computed locally, and an optional function ``Gcomm`` (of type
:c:type:`IDABBDCommFn`) which performs all inter-process communication necessary
to evaluate the approximate residual :math:`G`. These are in addition to the
user-supplied residual function ``res``. Both functions take as input the same
pointer ``user_data`` as passed by the user to :c:func:`IDASetUserData` and
passed to the user’s function ``res``. The user is responsible for providing
space (presumably within ``user_data``) for components of ``yy`` and ``yp`` that
are communicated by ``Gcomm`` from the other processors, and that are then used
by ``Gres``, which should not do any communication.

.. c:type:: int (*IDABBDLocalFn)(sunindextype Nlocal, realtype tt, N_Vector yy, N_Vector yp, N_Vector gval, void *user_data)

   This ``Gres`` function computes :math:`G(t,y,\dot{y})`. It loads the vector
   ``gval`` as a function of ``tt``, ``yy``, and ``yp``.

   **Arguments:**
      * ``Nlocal`` -- is the local vector length.
      * ``tt`` -- is the value of the independent variable.
      * ``yy`` -- is the dependent variable.
      * ``yp`` -- is the derivative of the dependent variable.
      * ``gval`` -- is the output vector.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.

   **Return value:**
      An :c:type:`IDABBDLocalFn` function type should return 0 to indicate success,
      1 for a recoverable error, or -1 for a non-recoverable error.

   **Notes:**
      This function must assume that all inter-processor communication of data
      needed to calculate ``gval`` has already been done, and this data is
      accessible within ``user_data``.

      The case where :math:`G` is mathematically identical to :math:`F` is allowed.

.. c:type:: int (*IDABBDCommFn)(sunindextype Nlocal, realtype tt, N_Vector yy, N_Vector yp, void *user_data)

   This ``Gcomm`` function performs all inter-processor communications necessary
   for the execution of the ``Gres`` function above, using the input vectors
   ``yy`` and ``yp``.

   **Arguments:**
      * ``Nlocal`` -- is the local vector length.
      * ``tt`` -- is the value of the independent variable.
      * ``yy`` -- is the dependent variable.
      * ``yp`` -- is the derivative of the dependent variable.
      * ``gval`` -- is the output vector.
      * ``user_data`` -- is a pointer to user data, the same as the ``user_data``
        parameter passed to :c:func:`IDASetUserData`.

   **Return value:**
      An :c:type:`IDABBDCommFn` function type should return 0 to indicate success,
      1 for a recoverable error, or -1 for a non-recoverable error.

   **Notes:**
      The ``Gcomm`` function is expected to save communicated data in space defined
      within the structure ``user_data``.

      Each call to the ``Gcomm`` function is preceded by a call to the residual
      function ``res`` with the same :math:`(t,y,\dot{y})` arguments. Thus
      ``Gcomm`` can omit any communications done by ``res`` if relevant to the
      evaluation of ``Gres``. If all necessary communication was done in ``res``,
      then ``Gcomm = NULL`` can be passed in the call to :c:func:`IDABBDPrecInit`.

Besides the header files required for the integration of the DAE problem (see
:numref:`IDA.Usage.CC.header_sim`), to use the IDABBDPRE module, the main program
must include the header file ``ida_bbdpre.h`` which declares the needed function
prototypes.

The following is a summary of the usage of this module and describes the
sequence of calls in the user main program.  Steps that are unchanged from the
user main program presented in :numref:`IDA.Usage.CC.skeleton_sim` are not bold.

#. Initialize parallel or multi-threaded environment (*if appropriate*)

#. Create the vector of initial values

#. Create matrix object (*if appropriate*)

#. **Create linear solver object** (*if appropriate*)

   When creating the iterative linear solver object, specify the use of left
   preconditioning (``SUN_PREC_LEFT``) as IDA only supports left preconditioning.

#. Create nonlinear solver object (*if appropriate*)

#. Create IDA object

#. Initialize IDA solver

#. Specify integration tolerances

#. Attach the linear solver (*if appropriate*)

#. **Set linear solver optional inputs** (*if appropriate*)

   Note that the user should not overwrite the preconditioner setup function or
   solve function through calls to :c:func:`IDASetPreconditioner` optional input
   function.

#. **Initialize the IDABBDPRE preconditioner module**

   Call :c:func:`IDABBDPrecInit` to allocate memory and initialize the internal
   preconditioner data. The last two arguments of :c:func:`IDABBDPrecInit` are
   the two user-supplied functions described above.

#. Attach nonlinear solver module (*if appropriate*)

#. Set nonlinear solver optional inputs (*if appropriate*)

#. Specify rootfinding problem (*optional*)

#. Set optional inputs

#. Advance solution in time

#. **Get optional outputs**

   Additional optional outputs associated with IDABBDPRE are available by way of
   two routines described below, :c:func:`IDABBDPrecGetWorkSpace` and
   :c:func:`IDABBDPrecGetNumGfnEvals`.

#. Deallocate memory

#. Finalize MPI, if used


The user-callable functions that initialize or re-initialize the IDABBDPRE
preconditioner module are described next.

.. c:function:: int IDABBDPrecInit(void *ida_mem, sunindextype Nlocal, sunindextype mudq, sunindextype mldq, sunindextype mukeep, sunindextype mlkeep, realtype dq_rel_yy, IDABBDLocalFn Gres, IDABBDCommFn Gcomm);

   The function ``IDABBDPrecInit`` initializes and allocates (internal) memory
   for the IDABBDPRE preconditioner.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``Nlocal`` -- local vector dimension.
      * ``mudq`` -- upper half-bandwidth to be used in the difference-quotient
        Jacobian approximation.
      * ``mldq`` -- lower half-bandwidth to be used in the difference-quotient
        Jacobian approximation.
      * ``mukeep`` -- upper half-bandwidth of the retained banded approximate
        Jacobian block.
      * ``mlkeep`` -- lower half-bandwidth of the retained banded approximate
        Jacobian block.
      * ``dq_rel_yy`` -- the relative increment in components of ``y`` used in the
        difference quotient approximations. The default is
        :math:`\mathtt{dq\_rel\_yy} = \sqrt{\text{unit roundoff}}` , which can be
        specified by passing :math:`\mathtt{dq\_rel\_yy} = 0.0`.
      * ``Gres`` -- the function which computes the local residual approximation
        :math:`G(t,y,\dot{y})`.
      * ``Gcomm`` -- the optional function which performs all inter-process
        communication required for the computation of :math:`G(t,y,\dot{y})`.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The call was successful.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer was ``NULL``.
      * ``IDALS_MEM_FAIL`` -- A memory allocation request has failed.
      * ``IDALS_LMEM_NULL`` -- An IDALS linear solver memory was not attached.
      * ``IDALS_ILL_INPUT`` -- The supplied vector implementation was not
        compatible with the block band preconditioner.

   **Notes:**
      If one of the half-bandwidths ``mudq`` or ``mldq`` to be used in the
      difference-quotient calculation of the approximate Jacobian is negative or
      exceeds the value ``Nlocal-1``, it is replaced by 0 or
      ``Nlocal-1`` accordingly.  The half-bandwidths ``mudq`` and ``mldq``
      need not be the true half-bandwidths of the Jacobian of the local block of
      :math:`G`, when smaller values may provide a greater efficiency.  Also, the
      half-bandwidths ``mukeep`` and ``mlkeep`` of the retained banded approximate
      Jacobian block may be even smaller, to reduce storage and computation costs
      further.  For all four half-bandwidths, the values need not be the same on
      every processor.


The IDABBDPRE module also provides a reinitialization function to allow for a
sequence of problems of the same size, with the same linear solver choice,
provided there is no change in ``local_N``, ``mukeep``, or ``mlkeep``. After
solving one problem, and after calling :c:func:`IDAReInit` to re-initialize IDA
for a subsequent problem, a call to :c:func:`IDABBDPrecReInit` can be made to
change any of the following: the half-bandwidths ``mudq`` and ``mldq`` used in
the difference-quotient Jacobian approximations, the relative increment
``dq_rel_yy``, or one of the user-supplied functions ``Gres`` and ``Gcomm``. If
there is a change in any of the linear solver inputs, an additional call to the
“Set”routines provided by the ``SUNLinearSolver`` object, and/or one or more of
the corresponding ``IDASet***`` functions, must also be made (in the proper
order).

.. c:function:: int IDABBDPrecReInit(void * ida_mem, sunindextype mudq, sunindextype mldq, realtype dq_rel_yy)

   The function ``IDABBDPrecReInit`` reinitializes the IDABBDPRE preconditioner.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``mudq`` -- upper half-bandwidth to be used in the difference-quotient
        Jacobian approximation.
      * ``Mldq`` -- lower half-bandwidth to be used in the difference-quotient
        Jacobian approximation.
      * ``dq_rel_yy`` -- the relative increment in components of ``y`` used in the
        difference quotient approximations. The default is
        :math:`\mathtt{dq\_rel\_yy} = \sqrt{\text{unit roundoff}}` , which can be
        specified by passing :math:`\mathtt{dq\_rel\_yy} = 0.0`.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The call was successful.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer was ``NULL``.
      * ``IDALS_LMEM_NULL`` -- An IDALS linear solver memory was not attached.
      * ``IDALS_PMEM_NULL`` -- The function :c:func:`IDABBDPrecInit` was not
        previously called.

   **Notes:**
      If one of the half-bandwidths ``mudq`` or ``mldq`` is negative or exceeds the
      value ``Nlocal - 1``, it is replaced by 0 or ``Nlocal - 1``, accordingly.


The following two optional output functions are available for use with the
IDABBDPRE module:

.. c:function:: int IDABBDPrecGetWorkSpace(void * ida_mem, long int * lenrwBBDP, long int * leniwBBDP)

   The function ``IDABBDPrecGetWorkSpace`` returns the local sizes of the
   IDABBDPRE real and integer workspaces.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``lenrwBBDP`` -- local number of real values in the IDABBDPRE workspace.
      * ``leniwBBDP`` -- local number of integer values in the IDABBDPRE workspace.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer was ``NULL``.
      * ``IDALS_PMEM_NULL`` -- The IDABBDPRE preconditioner has not been
        initialized.

   **Notes:**
      The workspace requirements reported by this routine correspond only to memory
      allocated within the IDABBDPRE module (the banded matrix approximation,
      banded ``SUNLinearSolver`` object, temporary vectors).  These values
      are local to each process.  The workspaces referred to here exist in addition
      to those given by the corresponding function :c:func:`IDAGetLinWorkSpace`.

.. c:function:: int IDABBDPrecGetNumGfnEvals(void * ida_mem, long int * ngevalsBBDP)

   The function ``IDABBDPrecGetNumGfnEvals`` returns the cumulative number of
   calls to the user ``Gres`` function due to the finite difference
   approximation of the Jacobian blocks used within IDABBDPRE's preconditioner
   setup function.

   **Arguments:**
      * ``ida_mem`` -- pointer to the IDA solver object.
      * ``ngevalsBBDP`` -- the cumulative number of calls to the user ``Gres``
        function.

   **Return value:**
      * ``IDALS_SUCCESS`` -- The optional output value has been successfully set.
      * ``IDALS_MEM_NULL`` -- The ``ida_mem`` pointer was ``NULL``.
      * ``IDALS_PMEM_NULL`` -- The IDABBDPRE preconditioner has not been
        initialized.


In addition to the ``ngevalsBBDP`` evaluations of ``Gres``, the costs associated
with IDABBDPRE also includes ``nlinsetups`` LU factorizations, ``nlinsetups``
calls to ``Gcomm``, ``npsolves`` banded backsolve calls, and ``nrevalsLS``
residual function evaluations, where ``nlinsetups`` is an optional IDA output
(see :numref:`IDA.Usage.CC.optional_output.optout_main`), and ``npsolves`` and
``nrevalsLS`` are linear solver optional outputs (see
:numref:`IDA.Usage.CC.optional_output.optout_ls`).
